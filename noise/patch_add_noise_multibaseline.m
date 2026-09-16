function noisy_patches = patch_add_noise_multibaseline(clean_patches, coh_maps, varargin)
%PATCH_ADD_NOISE_MULTIBASELINE  Add InSAR phase noise to a multi-baseline phase stack.
%
%   noisy_patches = PATCH_ADD_NOISE_MULTIBASELINE(clean_patches, coh_maps)
%   noisy_patches = PATCH_ADD_NOISE_MULTIBASELINE(clean_patches, coh_maps, 'param', value, ...)
%
%   Inputs:
%     clean_patches  {1 x K} cell array of [M x N] clean wrapped phase images
%                    (one per baseline), or a single [M x N] array that is
%                    broadcast to all K baselines.
%                    Values are expected in [-pi, pi].
%     coh_maps       {1 x K} cell array of [M x N] coherence maps in (0, 1),
%                    one per baseline.  Also accepts a single [M x N] array or
%                    scalar (broadcast to all baselines), or a {1 x K} cell
%                    produced by gen_coherence_map_multibaseline.
%
%   Name-Value Parameters:
%     'base_seed'  Integer base seed.  Baseline k receives seed (base_seed + k),
%                  guaranteeing INDEPENDENT noise realisations across baselines.
%                  Use [] to skip RNG seeding entirely (default 42).
%     'n_phi'      Number of PDF sample points on [-pi, pi]  (default 1000).
%     'n_bins'     Number of coherence quantisation bins      (default 50).
%
%   Output:
%     noisy_patches  {1 x K} cell array of [M x N] noisy wrapped phases
%                    in [-pi, pi], one per baseline.
%
%   Noise model (unchanged from patch_add_noise.m / patch_add_noise_spatialcoh.m):
%
%     p(phi; gamma) = (1 - gamma^2) / (2*pi)
%                   * 1 / (1 - gamma^2*cos^2(phi))
%                   * [1 + gamma*cos(phi)*acos(-gamma*cos(phi))
%                         / sqrt(1 - gamma^2*cos^2(phi))]
%
%   Noise is sampled per-pixel via inverse-CDF, identical to
%   patch_add_noise_spatialcoh.m, so results are fully backward-compatible
%   when K == 1.
%
%   Key design decisions:
%     1. DIFFERENT SEEDS PER BASELINE: each baseline k uses RNG seed
%        (base_seed + k), making noise realisations statistically independent.
%        This models the fact that, in a real SAR stack, acquisitions on
%        different dates / passes produce independent thermal noise.
%     2. SAME PHASE SIGNAL: the same clean_phase (derived from DEM elevation)
%        can be passed for all baselines — only the dem2phase_ratio differs
%        per baseline and should be applied BEFORE calling this function.
%     3. SPATIAL COH BROADCAST: coh_maps can be the output of
%        gen_coherence_map_multibaseline (same spatial pattern, different mean
%        per baseline) or entirely independent per-baseline maps.
%
%   Things to consider in your workflow (not handled here):
%     * Scale each clean phase by the appropriate dem2phase_ratio before
%       calling this function:
%           ratio_k = 4*pi/lambda_k * B_k / R_k / sind(theta_k);
%           clean_k = angle(exp(1i * dem_phase * ratio_k));
%     * Use gen_coherence_map_multibaseline to produce coh_maps so that the
%       spatial decorrelation pattern is physically consistent across baselines.
%     * After adding noise, call calc_coherence to re-estimate coherence from
%       noisy vs. clean phase and save it as the network input feature.
%
%   Example:
%     baselines  = [100, 300, 600];      % perpendicular baselines [m]
%     dem_phase  = randn(256);           % placeholder unwrapped phase [rad]
%
%     % Build per-baseline coherence maps (same spatial pattern, lower mean
%     % for longer baselines).
%     coh_maps = gen_coherence_map_multibaseline(256, 256, baselines, ...
%                    'max_coh', 0.9, 'decay_alpha', 0.5, 'seed', 1, ...
%                    'mode', 'region', 'n_regions', 6);
%
%     % Compute per-baseline InSAR parameters (ALOS-2 example).
%     lambda    = 0.236;  slant_r = 868142;  theta = 38.77;
%     clean_phs = cell(1, 3);
%     for k = 1:3
%         ratio        = 4*pi/lambda * baselines(k) / slant_r / sind(theta);
%         clean_phs{k} = angle(exp(1i * dem_phase * ratio));
%     end
%
%     % Add independent noise to each baseline.
%     noisy = patch_add_noise_multibaseline(clean_phs, coh_maps, 'base_seed', 42);
%
%   See also: gen_coherence_map_multibaseline, patch_add_noise_spatialcoh,
%             calc_coherence, gen_dataset_from_dem

    p = inputParser;
    addRequired(p,  'clean_patches', @(x) iscell(x) || isnumeric(x));
    addRequired(p,  'coh_maps',      @(x) iscell(x) || isnumeric(x));
    addParameter(p, 'base_seed', 42,   @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'n_phi',     1000, @(x) isnumeric(x) && isscalar(x) && x > 100);
    addParameter(p, 'n_bins',    50,   @(x) isnumeric(x) && isscalar(x) && x >= 1);
    parse(p, clean_patches, coh_maps, varargin{:});
    opt = p.Results;

    % ── Normalise inputs to cell arrays ───────────────────────────────────
    if ~iscell(clean_patches)
        % Single array: broadcast to all K baselines determined by coh_maps size.
        if iscell(coh_maps)
            K = numel(coh_maps);
        else
            K = 1;
        end
        clean_patches = repmat({clean_patches}, 1, K);
    end

    K = numel(clean_patches);

    if ~iscell(coh_maps)
        % Scalar or single array: broadcast to all K baselines.
        coh_maps = repmat({coh_maps}, 1, K);
    end

    assert(numel(coh_maps) == K, ...
        ['patch_add_noise_multibaseline: numel(coh_maps) (%d) must equal ' ...
         'numel(clean_patches) (%d).'], numel(coh_maps), K);

    % ── Validate individual arrays ─────────────────────────────────────────
    ref_size = size(clean_patches{1});
    for k = 1:K
        assert(isequal(size(clean_patches{k}), ref_size), ...
            'patch_add_noise_multibaseline: all clean_patches must have the same size.');

        % Allow scalar coherence (uniform map); otherwise must match patch size.
        cm = coh_maps{k};
        if ~isscalar(cm)
            assert(isequal(size(cm), ref_size), ...
                ['patch_add_noise_multibaseline: coh_maps{%d} size must match ' ...
                 'clean_patches size (%dx%d).'], k, ref_size(1), ref_size(2));
        end
    end

    % ── Add noise independently for each baseline ─────────────────────────
    noisy_patches = cell(1, K);

    for k = 1:K
        % Unique RNG seed per baseline: base_seed + k.
        % This ensures noise is drawn from an independent random state for
        % every baseline, even when the same clean phase and coherence map
        % are used.  Using [] skips seeding (non-reproducible but valid).
        if ~isempty(opt.base_seed)
            seed_k = opt.base_seed + k;
        else
            seed_k = [];
        end

        noisy_patches{k} = patch_add_noise_spatialcoh( ...
            clean_patches{k}, coh_maps{k}, ...
            'seed',   seed_k,     ...
            'n_phi',  opt.n_phi,  ...
            'n_bins', opt.n_bins);
    end
end

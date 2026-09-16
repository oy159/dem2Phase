function coh_maps = gen_coherence_map_multibaseline(rows, cols, baselines, varargin)
%GEN_COHERENCE_MAP_MULTIBASELINE  Generate coherence maps for a multi-baseline InSAR stack.
%
%   coh_maps = GEN_COHERENCE_MAP_MULTIBASELINE(rows, cols, baselines)
%   coh_maps = GEN_COHERENCE_MAP_MULTIBASELINE(rows, cols, baselines, 'param', value, ...)
%
%   Inputs:
%     rows, cols   output map size (pixels)
%     baselines    [1 x K] perpendicular baseline lengths [m], all > 0
%
%   Name-Value Parameters:
%     'min_coh'        absolute minimum coherence (default 0.01)
%     'max_coh'        coherence upper bound at the reference (shortest) baseline
%                      (default 0.95)
%     'spatial_scale'  spatial correlation / boundary-smoothing length [px] (default 30)
%     'seed'           RNG seed for the shared spatial base field; [] = no reset
%                      (default [])
%     'mode'           'smooth' or 'region' — same modes as gen_coherence_map
%                      (default 'region')
%     'n_regions'      number of Voronoi regions for 'region' mode (default 4)
%     'baseline_ref'   reference baseline [m] that anchors max_coh; defaults to
%                      min(baselines), i.e. the shortest baseline is the cleanest.
%     'decay_alpha'    exponential decay rate (default 0.5).
%                      The upper-coherence bound for baseline B is:
%                        max_k = min_coh + (max_coh - min_coh)
%                                * exp(-decay_alpha * (B - B_ref) / B_ref)
%                      alpha = 0  → all baselines get identical max_coh
%                      alpha = 1  → doubling the baseline halves the dynamic range
%     'thermal_floor'  minimum coherence added to every pixel to model
%                      thermal / receiver-noise decorrelation (default 0.0).
%                      Useful for simulating high-noise systems; set to
%                      ~1/(2*N_looks) for a realistic thermal floor.
%
%   Output:
%     coh_maps   {1 x K} cell array of [rows x cols] double arrays in
%                [min_coh, max_k], one map per baseline.
%
%   Design rationale:
%     A single random spatial base field (normalised to [0, 1]) is generated
%     ONCE and then linearly rescaled to [min_coh, max_k] for each baseline k.
%     This ensures that:
%       (a) the SPATIAL PATTERN is identical across baselines — pixels that
%           correspond to decorrelated land cover (water, vegetation) remain
%           relatively decorrelated for all baselines, matching real InSAR.
%       (b) the MEAN coherence decreases monotonically with |baseline|,
%           consistent with geometric / volumetric decorrelation theory.
%
%   Example:
%     baselines = [100, 200, 400];   % three perpendicular baselines [m]
%     coh = gen_coherence_map_multibaseline(256, 256, baselines, ...
%               'max_coh', 0.9, 'decay_alpha', 0.6, 'seed', 7, ...
%               'mode', 'region', 'n_regions', 6);
%     % coh{1} has the highest mean; coh{3} the lowest
%     figure;
%     for k = 1:3
%         subplot(1,3,k); imagesc(coh{k}); caxis([0 1]); colorbar;
%         title(sprintf('B=%d m, mean=%.2f', baselines(k), mean(coh{k}(:))));
%     end
%
%   See also: gen_coherence_map, patch_add_noise_multibaseline, calc_coherence

    p = inputParser;
    addRequired(p,  'rows',          @(x) isnumeric(x) && isscalar(x) && x > 0);
    addRequired(p,  'cols',          @(x) isnumeric(x) && isscalar(x) && x > 0);
    addRequired(p,  'baselines',     @(x) isnumeric(x) && isvector(x) && numel(x) >= 1 && all(x > 0));
    addParameter(p, 'min_coh',       0.01,     @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'max_coh',       0.95,     @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'spatial_scale', 30,       @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'seed',          [],       @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'mode',          'region', @(x) ismember(x, {'smooth', 'region'}));
    addParameter(p, 'n_regions',     4,        @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'baseline_ref',  [],       @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    addParameter(p, 'decay_alpha',   0.5,      @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'thermal_floor', 0.0,      @(x) isnumeric(x) && isscalar(x) && x >= 0);
    parse(p, rows, cols, baselines, varargin{:});
    opt = p.Results;

    assert(opt.min_coh > 0,          'gen_coherence_map_multibaseline: min_coh must be > 0.');
    assert(opt.max_coh < 1,          'gen_coherence_map_multibaseline: max_coh must be < 1.');
    assert(opt.min_coh < opt.max_coh,'gen_coherence_map_multibaseline: min_coh must be < max_coh.');
    assert(opt.thermal_floor < opt.min_coh || opt.thermal_floor == 0, ...
        'gen_coherence_map_multibaseline: thermal_floor must be 0 or < min_coh.');

    baselines = baselines(:)';   % enforce row vector
    K = numel(baselines);

    if isempty(opt.baseline_ref)
        opt.baseline_ref = min(baselines);
    end

    % Seed the RNG before generating the shared spatial base field.
    if ~isempty(opt.seed)
        rng(opt.seed);
    end

    % ── 1. Shared normalised spatial base field [0, 1] ────────────────────
    base_field = generate_base_field(rows, cols, opt.mode, opt.spatial_scale, opt.n_regions);

    % ── 2. Per-baseline coherence maps ────────────────────────────────────
    coh_maps = cell(1, K);
    B_ref    = opt.baseline_ref;

    for k = 1:K
        B = baselines(k);

        % Exponential decay of the dynamic range with baseline length.
        % At B == B_ref the decay factor is 1 (full [min,max] range).
        % Longer baselines shrink the upper bound toward min_coh.
        decay = exp(-opt.decay_alpha * (B - B_ref) / B_ref);
        decay = max(decay, 0);

        max_k = opt.min_coh + (opt.max_coh - opt.min_coh) * decay;
        max_k = min(max_k, 0.9999);            % hard upper cap
        max_k = max(max_k, opt.min_coh + 1e-4); % ensure non-zero range

        % Warn when the resulting max coherence is very low (very long baseline)
        if max_k < 0.1
            warning('gen_coherence_map_multibaseline:lowCoherence', ...
                ['Baseline %g m yields max_coh ≈ %.3f. ' ...
                 'Consider increasing baseline_ref or reducing decay_alpha.'], ...
                B, max_k);
        end

        coh_k = opt.min_coh + base_field * (max_k - opt.min_coh);

        % Apply thermal noise floor: coherence is at least thermal_floor
        if opt.thermal_floor > 0
            coh_k = max(coh_k, opt.thermal_floor);
        end

        % Final clamp to (0, 1) — strict open interval required by noise PDF
        coh_k = min(max(coh_k, 1e-4), 0.9999);

        coh_maps{k} = coh_k;
    end
end

% ── Private helpers ────────────────────────────────────────────────────────

function field = generate_base_field(rows, cols, mode, sigma, n_regions)
%GENERATE_BASE_FIELD  Return a [rows x cols] normalised field in [0, 1].
    switch mode
        case 'smooth'
            % Gaussian-filtered white noise: smooth gradual coherence variations.
            noise  = randn(rows, cols);
            kernel = gaussian_kernel_2d(sigma);
            field  = conv2(noise, kernel, 'same');

        case 'region'
            % Voronoi piecewise-constant regions with softened boundaries,
            % mimicking distinct land-cover types (bare soil, vegetation, water).
            centers_r  = 1 + rand(n_regions, 1) * (rows - 1);
            centers_c  = 1 + rand(n_regions, 1) * (cols - 1);
            coh_values = rand(n_regions, 1);   % uniform levels per region

            [C, R] = meshgrid(1:cols, 1:rows);

            dist2 = zeros(rows, cols, n_regions);
            for k = 1:n_regions
                dist2(:,:,k) = (R - centers_r(k)).^2 + (C - centers_c(k)).^2;
            end
            [~, idx] = min(dist2, [], 3);
            field = reshape(coh_values(idx(:)), rows, cols);

            % Lightly blur region boundaries (same sigma scaling as gen_coherence_map)
            blur_sigma = max(1, sigma / 5);
            kernel     = gaussian_kernel_2d(blur_sigma);
            field      = conv2(field, kernel, 'same');
    end

    % Normalise to [0, 1] so per-baseline rescaling is unambiguous.
    f_min = min(field(:));
    f_max = max(field(:));
    if f_max > f_min
        field = (field - f_min) / (f_max - f_min);
    else
        field = 0.5 * ones(rows, cols);
    end
end

function h = gaussian_kernel_2d(sigma)
%GAUSSIAN_KERNEL_2D  Normalised 2-D Gaussian kernel (no Image Processing Toolbox).
    r  = ceil(3 * sigma);
    x  = -r:r;
    h1 = exp(-x.^2 / (2 * sigma^2));
    h  = h1' * h1;
    h  = h / sum(h(:));
end

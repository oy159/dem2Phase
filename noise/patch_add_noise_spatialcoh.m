function noisy_patch = patch_add_noise_spatialcoh(clean_patch, coh_map, varargin)
%PATCH_ADD_NOISE_SPATIALCOH  Add InSAR phase noise with spatially-varying coherence.
%
%   noisy_patch = PATCH_ADD_NOISE_SPATIALCOH(clean_patch, coh_map)
%   noisy_patch = PATCH_ADD_NOISE_SPATIALCOH(clean_patch, coh_map, 'param', value, ...)
%
%   Inputs:
%     clean_patch  [M x N] clean wrapped phase in [-pi, pi]
%     coh_map      [M x N] coherence map in (0, 1), or a scalar coherence
%
%   Name-Value Parameters:
%     'n_phi'   number of PDF sample points on [-pi, pi]  (default 1000)
%     'n_bins'  coherence quantisation bins               (default 50)
%     'seed'    random seed; [] = no reset                (default [])
%
%   Output:
%     noisy_patch  [M x N] noisy wrapped phase in [-pi, pi]
%
%   Per-pixel phase noise PDF (identical to patch_add_noise.m, unchanged):
%
%     p(phi; gamma) = (1 - gamma^2) / (2*pi)
%                   * 1 / (1 - gamma^2 * cos^2(phi))
%                   * [1 + gamma*cos(phi)*acos(-gamma*cos(phi))
%                       / sqrt(1 - gamma^2*cos^2(phi))]
%
%   Algorithm:
%     1. Quantise coh_map into n_bins coherence levels.
%     2. For each bin, compute PDF/CDF using the formula above, then draw
%        noise samples via inverse-CDF (the same method as patch_add_noise.m).
%     3. Map sampled noise back to each pixel's spatial location.
%     4. Add noise to clean_patch and wrap to [-pi, pi].
%
%   When coh_map is a scalar this function is equivalent to patch_add_noise.m.
%
%   Example:
%     coh_map     = gen_coherence_map(256, 256, 'seed', 1);
%     clean_phase = angle(exp(1i * randn(256)));   % dummy clean phase
%     noisy_phase = patch_add_noise_spatialcoh(clean_phase, coh_map, 'seed', 42);
%
%   See also: gen_coherence_map, patch_add_noise, calc_coherence

    p = inputParser;
    addRequired(p,  'clean_patch', @(x) isnumeric(x) && ismatrix(x));
    addRequired(p,  'coh_map',     @(x) isnumeric(x));
    addParameter(p, 'n_phi',  1000, @(x) isnumeric(x) && isscalar(x) && x > 100);
    addParameter(p, 'n_bins', 50,   @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'seed',   [],   @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    parse(p, clean_patch, coh_map, varargin{:});
    opt = p.Results;

    rng_cleanup = scoped_rng(opt.seed); %#ok<NASGU>

    [M, N] = size(clean_patch);

    % Accept scalar coherence (uniform case, equivalent to patch_add_noise.m)
    if isscalar(coh_map)
        coh_map = coh_map * ones(M, N);
    end

    assert(isequal(size(coh_map), [M, N]), ...
        'patch_add_noise_spatialcoh: coh_map must have the same size as clean_patch.');
    assert(all(coh_map(:) > 0 & coh_map(:) < 1), ...
        'patch_add_noise_spatialcoh: coherence values must all be in (0, 1).');

    % ------------------------------------------------------------------
    % Quantise coherence map into n_bins discrete levels
    % ------------------------------------------------------------------
    coh_lo   = min(coh_map(:));
    coh_hi   = max(coh_map(:));
    n_bins   = opt.n_bins;

    if coh_lo == coh_hi
        % Uniform coherence — single bin
        bin_centers = coh_lo;
        n_bins      = 1;
    else
        bin_edges   = linspace(coh_lo, coh_hi, n_bins + 1);
        bin_centers = 0.5 * (bin_edges(1:end-1) + bin_edges(2:end));
    end

    % Map each pixel to a bin index (1 … n_bins)
    range = coh_hi - coh_lo + eps;
    bin_idx = min(n_bins, max(1, ...
        floor((coh_map - coh_lo) / range * n_bins) + 1));

    % ------------------------------------------------------------------
    % phi axis (shared across all bins)
    % ------------------------------------------------------------------
    phi_range = 2 * pi;
    phi_axis  = linspace(-phi_range/2, phi_range/2, opt.n_phi);

    % ------------------------------------------------------------------
    % Sample noise bin-by-bin using the SAME PDF as patch_add_noise.m
    % ------------------------------------------------------------------
    phase_noise = zeros(M, N);

    for b = 1:n_bins
        mask     = (bin_idx == b);
        n_pixels = sum(mask(:));
        if n_pixels == 0
            continue;
        end

        gamma = bin_centers(b);

        % ---- PDF (formula preserved from patch_add_noise.m) ----------
        cos_phi  = cos(phi_axis);
        denom    = 1 - gamma^2 * cos_phi.^2;
        denom    = max(denom, eps);          % guard numerical underflow

        phase_noise_pdf = ((1 - gamma^2) ./ (2 * pi)) ...
            .* (1 ./ denom) ...
            .* (1 + (gamma .* cos_phi .* acos(-gamma .* cos_phi)) ...
                ./ sqrt(denom));

        phase_noise_pdf = max(phase_noise_pdf, 0);   % remove rounding negatives

        % ---- CDF (same normalisation as patch_add_noise.m) -----------
        phase_noise_cdf = cumsum(phase_noise_pdf);
        phase_noise_cdf = phase_noise_cdf / sum(phase_noise_pdf);

        % ---- ICDF sampling (same method as patch_add_noise.m) --------
        quantiles  = rand(1, n_pixels);
        noise_vals = interp1(phase_noise_cdf, phi_axis, quantiles, 'linear', 0);
        noise_vals(isnan(noise_vals)) = 0;

        phase_noise(mask) = noise_vals;
    end

    % ------------------------------------------------------------------
    % Add noise and wrap to [-pi, pi]
    % ------------------------------------------------------------------
    noisy_patch = angle(exp(1i * (clean_patch + phase_noise)));
end

function cleanup = scoped_rng(seed)
%SCOPED_RNG Seed locally without changing the caller's random stream.
    cleanup = [];
    if ~isempty(seed)
        previous_state = rng;
        cleanup = onCleanup(@() rng(previous_state));
        rng(seed, 'twister');
    end
end

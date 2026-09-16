function coh_map = gen_coherence_map(rows, cols, varargin)
%GEN_COHERENCE_MAP  Generate a spatially-varying InSAR coherence map.
%
%   coh_map = GEN_COHERENCE_MAP(rows, cols)
%   coh_map = GEN_COHERENCE_MAP(rows, cols, 'param', value, ...)
%
%   Inputs:
%     rows, cols      output map size (pixels)
%
%   Name-Value Parameters:
%     'min_coh'        minimum coherence value            (default 0.3)
%     'max_coh'        maximum coherence value            (default 0.9)
%     'spatial_scale'  spatial correlation length (px)    (default 30)
%     'seed'           random seed; [] for random         (default [])
%     'mode'           'smooth' or 'region'               (default 'smooth')
%     'n_regions'      number of Voronoi regions          (default 4)
%
%   Output:
%     coh_map  [rows x cols] double array in [min_coh, max_coh]
%
%   Modes:
%     'smooth'  Spatially-correlated random field via Gaussian-filtered
%               white noise. Controls realistic gradual transitions.
%     'region'  Piecewise-constant Voronoi regions (e.g. bare soil vs.
%               vegetation) with lightly smoothed boundaries.
%
%   Example:
%     % Smooth coherence map
%     coh = gen_coherence_map(256, 256, 'spatial_scale', 50, 'seed', 1);
%     figure; imagesc(coh); colorbar; caxis([0 1]);
%     title('Spatially-varying coherence (smooth)');
%
%     % Region-based coherence map
%     coh = gen_coherence_map(256, 256, 'mode', 'region', 'n_regions', 5, 'seed', 1);
%
%   See also: patch_add_noise_spatialcoh, calc_coherence

    p = inputParser;
    addRequired(p,  'rows',           @(x) isnumeric(x) && isscalar(x) && x > 0);
    addRequired(p,  'cols',           @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'min_coh',        0.3,      @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'max_coh',        0.9,      @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'spatial_scale',  30,       @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'seed',           [],       @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'mode',           'smooth', @(x) ismember(x, {'smooth', 'region'}));
    addParameter(p, 'n_regions',      4,        @(x) isnumeric(x) && isscalar(x) && x >= 1);
    parse(p, rows, cols, varargin{:});
    opt = p.Results;

    assert(opt.min_coh > 0,          'gen_coherence_map: min_coh must be > 0.');
    assert(opt.max_coh < 1,          'gen_coherence_map: max_coh must be < 1.');
    assert(opt.min_coh < opt.max_coh,'gen_coherence_map: min_coh must be less than max_coh.');

    if ~isempty(opt.seed)
        rng(opt.seed);
    end

    switch opt.mode
        case 'smooth'
            coh_map = smooth_coh_map(rows, cols, opt.min_coh, opt.max_coh, ...
                                     opt.spatial_scale);
        case 'region'
            coh_map = region_coh_map(rows, cols, opt.min_coh, opt.max_coh, ...
                                     opt.n_regions, opt.spatial_scale);
    end
end

% =========================================================================
function coh_map = smooth_coh_map(rows, cols, min_coh, max_coh, sigma)
% Gaussian-filtered white noise -> spatially correlated coherence field.
    noise        = randn(rows, cols);
    kernel       = gaussian_kernel_2d(sigma);
    smooth_field = conv2(noise, kernel, 'same');

    % Normalize to [0, 1]
    f_min = min(smooth_field(:));
    f_max = max(smooth_field(:));
    if f_max > f_min
        smooth_field = (smooth_field - f_min) / (f_max - f_min);
    else
        smooth_field = 0.5 * ones(rows, cols);
    end

    coh_map = min_coh + smooth_field * (max_coh - min_coh);
end

% =========================================================================
function coh_map = region_coh_map(rows, cols, min_coh, max_coh, n_regions, sigma)
% Voronoi region-based coherence map (e.g. land-cover types).
    centers_r  = 1 + rand(n_regions, 1) * (rows - 1);
    centers_c  = 1 + rand(n_regions, 1) * (cols - 1);
    coh_values = min_coh + rand(n_regions, 1) * (max_coh - min_coh);

    [C, R] = meshgrid(1:cols, 1:rows);   % R(i,j)=row index, C(i,j)=col index

    % Nearest-center assignment (Voronoi)
    dist2 = zeros(rows, cols, n_regions);
    for k = 1:n_regions
        dist2(:,:,k) = (R - centers_r(k)).^2 + (C - centers_c(k)).^2;
    end
    [~, idx]  = min(dist2, [], 3);
    coh_map   = reshape(coh_values(idx(:)), rows, cols);

    % Lightly smooth region boundaries
    blur_sigma = max(1, sigma / 5);
    kernel     = gaussian_kernel_2d(blur_sigma);
    coh_map    = conv2(coh_map, kernel, 'same');
    coh_map    = min(max(coh_map, min_coh), max_coh);
end

% =========================================================================
function h = gaussian_kernel_2d(sigma)
% 2-D normalised Gaussian kernel (no Image Processing Toolbox required).
    r  = ceil(3 * sigma);
    x  = -r:r;
    h1 = exp(-x.^2 / (2 * sigma^2));
    h  = h1' * h1;
    h  = h / sum(h(:));
end

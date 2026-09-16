function [coh_maps, terrain] = gen_coherence_map_multibaseline_terrainaware(dem_patch, baselines, varargin)
%GEN_COHERENCE_MAP_MULTIBASELINE_TERRAINAWARE Terrain-conditioned coherence.
%
% The deterministic component is derived from DEM slope, local incidence,
% curvature, roughness, layover and shadow. A shared smooth residual field
% represents scene-level uncertainty, while an independent residual field is
% generated for every interferometric edge. This avoids both extremes of
% identical and fully independent coherence maps across baselines.

    p = inputParser;
    addRequired(p, 'dem_patch', @(x) isnumeric(x) && ismatrix(x) && ~isempty(x));
    addRequired(p, 'baselines', @(x) isnumeric(x) && isvector(x) && all(x > 0));
    addParameter(p, 'pixel_spacing_m', 10, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'incidence_angle_deg', 45, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 90);
    addParameter(p, 'look_azimuth_deg', 90, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'min_coh', 0.01, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'max_coh', 0.95, @(x) isnumeric(x) && isscalar(x) && x < 1);
    addParameter(p, 'baseline_ref', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    addParameter(p, 'baseline_decay_model', 'exponential', ...
        @(x) ischar(x) || isstring(x));
    addParameter(p, 'critical_baseline_m', [], ...
        @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    addParameter(p, 'decay_alpha', 0.35, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'terrain_normalization_model', 'patch_percentile', ...
        @(x) ischar(x) || isstring(x));
    addParameter(p, 'slope_scale_deg', 45, ...
        @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'roughness_scale_m', 50, ...
        @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'curvature_scale_per_m', 0.01, ...
        @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'slope_weight', 0.9, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'roughness_weight', 0.8, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'curvature_weight', 0.35, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'shared_residual_std', 0.12, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'edge_residual_std', 0.08, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'residual_scale_px', 24, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'seed', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    parse(p, dem_patch, baselines, varargin{:});
    opt = p.Results;

    assert(opt.min_coh < opt.max_coh, ...
        'dem2phase:InvalidCoherenceRange', 'min_coh must be smaller than max_coh.');
    rng_cleanup = scoped_rng(opt.seed); %#ok<NASGU>

    dem = double(dem_patch);
    invalid = ~isfinite(dem);
    if any(invalid(:))
        valid_values = dem(~invalid);
        assert(~isempty(valid_values), 'dem2phase:AllInvalidDem', ...
            'dem_patch contains no finite elevation samples.');
        dem(invalid) = median(valid_values);
    end

    spacing = opt.pixel_spacing_m;
    % MATLAB returns the column/x derivative first and row/y derivative second.
    [dz_dx, dz_dy] = gradient(dem, spacing, spacing);
    slope_rad = atan(hypot(dz_dx, dz_dy));
    slope_deg = rad2deg(slope_rad);
    aspect_rad = atan2(dz_dy, dz_dx);

    look_azimuth_rad = deg2rad(opt.look_azimuth_deg);
    theta_rad = deg2rad(opt.incidence_angle_deg);
    directional_gradient = dz_dx*cos(look_azimuth_rad) + dz_dy*sin(look_azimuth_rad);
    range_slope_rad = atan(directional_gradient);

    normal_x = -dz_dx ./ sqrt(1 + dz_dx.^2 + dz_dy.^2);
    normal_y = -dz_dy ./ sqrt(1 + dz_dx.^2 + dz_dy.^2);
    normal_z = 1 ./ sqrt(1 + dz_dx.^2 + dz_dy.^2);
    look_x = sin(theta_rad)*cos(look_azimuth_rad);
    look_y = sin(theta_rad)*sin(look_azimuth_rad);
    look_z = cos(theta_rad);
    cos_local_incidence = normal_x*look_x + normal_y*look_y + normal_z*look_z;
    cos_local_incidence = min(max(cos_local_incidence, -1), 1);
    local_incidence_deg = rad2deg(acos(cos_local_incidence));

    local_mean = box_mean(dem, 7);
    local_mean_sq = box_mean(dem.^2, 7);
    roughness_m = sqrt(max(local_mean_sq - local_mean.^2, 0));
    curvature = del2(dem, spacing, spacing);
    tpi_m = dem - box_mean(dem, 15);

    normalization_model = char(opt.terrain_normalization_model);
    if strcmp(normalization_model, 'absolute_scales')
        roughness_unit = min(max(roughness_m/opt.roughness_scale_m, 0), 1);
        curvature_unit = min(max(abs(curvature)/opt.curvature_scale_per_m, 0), 1);
    else
        assert(strcmp(normalization_model, 'patch_percentile'), ...
            'dem2phase:UnsupportedTerrainNormalizationModel', ...
            'Terrain normalization must be absolute_scales or patch_percentile.');
        roughness_unit = robust_unit_interval(roughness_m);
        curvature_unit = robust_unit_interval(abs(curvature));
    end
    slope_term = (slope_deg / opt.slope_scale_deg).^2;
    terrain_quality = exp(-opt.slope_weight*slope_term ...
        - opt.roughness_weight*roughness_unit ...
        - opt.curvature_weight*curvature_unit);

    orientation_quality = sqrt(max(cos_local_incidence, 0));
    terrain_quality = terrain_quality .* orientation_quality;

    % First-order side-looking geometry masks. These are diagnostic proxies,
    % not a substitute for ray tracing through the full DEM.
    layover_mask = range_slope_rad > theta_rad;
    shadow_mask = range_slope_rad < -(pi/2 - theta_rad) | cos_local_incidence <= 0;
    invalid_geometry_mask = layover_mask | shadow_mask | invalid;
    terrain_quality(invalid_geometry_mask) = 0;
    terrain_quality = min(max(terrain_quality, 0), 1);

    baselines = double(baselines(:)');
    if isempty(opt.baseline_ref)
        baseline_ref = min(baselines);
    else
        baseline_ref = opt.baseline_ref;
    end

    shared_residual = smooth_standard_field(size(dem), opt.residual_scale_px);
    coh_maps = cell(1, numel(baselines));
    for k = 1:numel(baselines)
        decay_model = char(opt.baseline_decay_model);
        if strcmp(decay_model, 'critical_baseline')
            assert(~isempty(opt.critical_baseline_m), ...
                'dem2phase:MissingCriticalBaseline', ...
                'critical_baseline_m is required for critical_baseline decay.');
            baseline_decay = max(0, 1 - abs(baselines(k))/opt.critical_baseline_m);
        else
            assert(strcmp(decay_model, 'exponential'), ...
                'dem2phase:UnsupportedBaselineDecayModel', ...
                'baseline_decay_model must be critical_baseline or exponential.');
            baseline_decay = exp(-opt.decay_alpha * ...
                max(baselines(k) - baseline_ref, 0) / baseline_ref);
        end
        edge_residual = smooth_standard_field(size(dem), opt.residual_scale_px/2);
        residual_multiplier = exp(opt.shared_residual_std*shared_residual ...
            + opt.edge_residual_std*edge_residual);
        quality_k = min(max(terrain_quality .* residual_multiplier, 0), 1);
        max_k = opt.min_coh + (opt.max_coh - opt.min_coh)*baseline_decay;
        coh_k = opt.min_coh + quality_k*(max_k - opt.min_coh);
        coh_k(invalid_geometry_mask) = opt.min_coh;
        coh_maps{k} = min(max(coh_k, 1e-4), 0.9999);
    end

    terrain.slope_deg = slope_deg;
    terrain.aspect_rad = aspect_rad;
    terrain.range_slope_deg = rad2deg(range_slope_rad);
    terrain.local_incidence_deg = local_incidence_deg;
    terrain.roughness_m = roughness_m;
    terrain.curvature = curvature;
    terrain.tpi_m = tpi_m;
    terrain.terrain_quality = terrain_quality;
    terrain.layover_mask = layover_mask;
    terrain.shadow_mask = shadow_mask;
    terrain.invalid_geometry_mask = invalid_geometry_mask;
    terrain.pixel_spacing_m = spacing;
    terrain.look_azimuth_deg = opt.look_azimuth_deg;
    terrain.baseline_decay_model = char(opt.baseline_decay_model);
    terrain.critical_baseline_m = opt.critical_baseline_m;
    terrain.normalization_model = normalization_model;
    terrain.slope_scale_deg = opt.slope_scale_deg;
    terrain.roughness_scale_m = opt.roughness_scale_m;
    terrain.curvature_scale_per_m = opt.curvature_scale_per_m;
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

function result = box_mean(data, window_size)
    kernel = ones(window_size, window_size) / (window_size^2);
    result = conv2(data, kernel, 'same');
end

function unit = robust_unit_interval(data)
    values = sort(data(isfinite(data)));
    if isempty(values)
        unit = zeros(size(data));
        return
    end
    lo = values(max(1, round(0.02*numel(values))));
    hi = values(max(1, round(0.98*numel(values))));
    if hi <= lo
        unit = zeros(size(data));
    else
        unit = min(max((data - lo) / (hi - lo), 0), 1);
    end
end

function field = smooth_standard_field(array_size, sigma)
    sigma = max(sigma, 0.5);
    radius = ceil(3*sigma);
    x = -radius:radius;
    h1 = exp(-(x.^2)/(2*sigma^2));
    kernel = h1'*h1;
    kernel = kernel/sum(kernel(:));
    field = conv2(randn(array_size), kernel, 'same');
    field = field - mean(field(:));
    scale = std(field(:));
    if scale > eps
        field = field/scale;
    end
end

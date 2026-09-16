%% gen_dataset_from_dem_v13.m
%
% V13 Multi-Baseline InSAR Dataset Generator
%
% Generates multi-baseline InSAR training datasets from DEM / DSM GeoTIFF files.
% Based on the v13 design plan: stores only essential data, computes real/imag
% and wrap_count on-the-fly in the Python loader.
%
% Input/output paths are relative to dem2phase/ and configured in the JSON
% profile. GeoTIFF names must contain *_dem* or *_DSM*.
%
% Functions required (all in dem2phase/):
%   gen_coherence_map_multibaseline_terrainaware.m
%   patch_add_noise_multibaseline.m
%   calc_coherence.m
%
% Default output: one patch_groups/<patch_name>.mat file containing [K,H,W]
% stacks, edge geometry and a valid_edge_mask. Legacy per-edge output can be
% enabled in the JSON profile for compatibility.
%
% -------------------------------------------------------------------------
clc; close all;

%% ══════════════════════════════════════════════════════════════════════════
%  USER CONFIGURATION  –  edit these values before running the script
%  ══════════════════════════════════════════════════════════════════════════

% --- Multi-UAV InSAR physics ---------------------------------------------
% JSON profiles define frequency, geometry, measurement mode, UAV count,
% interferometric graph and the perpendicular baseline assigned to each edge.
sim_config_filename = 'uav_p_500m_monostatic.json';

% Coherence estimation window (for calc_coherence after adding noise)
coh_est_win = 7;

%% ══════════════════════════════════════════════════════════════════════════
%  END USER CONFIGURATION
%  ══════════════════════════════════════════════════════════════════════════

%% ── 0.  Paths ─────────────────────────────────────────────────────────────
script_dir    = fileparts(mfilename('fullpath'));   % dem2phase/scripts/
dem2phase_dir = fileparts(script_dir);             % dem2phase/
addpath(genpath(dem2phase_dir));  % adds dem2phase/ and all subdirs

sim_config_path = fullfile(dem2phase_dir, 'configs', sim_config_filename);
sim_cfg = load_sim_config(sim_config_path);

% Optional paired Phase-4 sweep override. Environment variables allow the
% sweep runner to reuse identical crops and non-Phase-4 random draws.
sweep_scenario = getenv('DEM2PHASE_PHASE4_SCENARIO');
sweep_level = getenv('DEM2PHASE_PHASE4_LEVEL');
if ~isempty(sweep_scenario)
    assert(~isempty(sweep_level), 'dem2phase:MissingPhase4SweepLevel', ...
        'DEM2PHASE_PHASE4_LEVEL is required when a sweep scenario is set.');
    sim_cfg = apply_phase4_sweep_override(sim_cfg, sweep_scenario, sweep_level);
    sweep_patch_count = str2double(getenv('DEM2PHASE_PHASE4_PATCHES'));
    if isfinite(sweep_patch_count) && sweep_patch_count >= 1
        sim_cfg.dataset.generation.patches_per_dem = floor(sweep_patch_count);
    end
    sim_cfg = validate_sim_config(sim_cfg);
end

generation_cfg = sim_cfg.dataset.generation;
patches_per_dem = double(generation_cfg.patches_per_dem);
patch_size = double(generation_cfg.patch_size);
n_precomp_scales = double(generation_cfg.n_precomp_scales);
interp_scale_range = double(generation_cfg.interp_scale_range(:)');
rng_seed = double(generation_cfg.rng_seed);
max_wrap_count = double(generation_cfg.max_wrap_count);
void_frac_thresh = double(generation_cfg.void_fraction_threshold);

physics_lambda = sim_cfg.physics.wavelength_m;
physics_slant_range = double(sim_cfg.geometry.slant_range_m);
physics_incidence_angle = double(sim_cfg.geometry.incidence_angle_deg);
physics_baseline_lengths = double(sim_cfg.interferometry.baseline_perp_m(:)');
phase_path_multiplicity = double(sim_cfg.interferometry.phase_path_multiplicity);
dem2phase_ratios = double(sim_cfg.physics.dem2phase_ratios_rad_per_m(:)');
ambiguity_heights_m = double(sim_cfg.physics.ambiguity_heights_m(:)');
K = numel(physics_baseline_lengths);
[short_baseline_length, short_baseline_idx] = min(physics_baseline_lengths);

coh_min = double(sim_cfg.coherence.min);
coh_max = double(sim_cfg.coherence.max);
coh_decay_alpha = double(sim_cfg.coherence.decay_alpha);
coh_ref_baseline = short_baseline_length;
if strcmp(char(sim_cfg.coherence.baseline_decay_model), 'critical_baseline')
    baseline_coherence_factor = max(0, 1 - abs(physics_baseline_lengths) / ...
        sim_cfg.physics.critical_baseline_m);
else
    baseline_coherence_factor = exp(-coh_decay_alpha * ...
        (physics_baseline_lengths-coh_ref_baseline)/coh_ref_baseline);
end
max_coh_gen = baseline_coherence_factor*(coh_max-coh_min) + coh_min;

dem_dir = resolve_project_path(dem2phase_dir, char(sim_cfg.dataset.dem_directory));
output_root = resolve_project_path(dem2phase_dir, char(sim_cfg.dataset.output_directory));
source_dem_pixel_spacing_m = double(sim_cfg.dataset.source_dem_pixel_spacing_m);
look_azimuth_deg = double(sim_cfg.geometry.look_azimuth_deg);
terrain_sampling = sim_cfg.dataset.terrain_sampling;
storage_cfg = sim_cfg.dataset.storage;
write_legacy_per_edge = logical(storage_cfg.write_legacy_per_edge);
edge_sampling = sim_cfg.interferometry.edge_sampling;
noise_cfg = sim_cfg.noise;
phase4_cfg = sim_cfg.phase4_errors;
failure_cfg = sim_cfg.phase4_failures;
landcover_cfg = sim_cfg.dataset.landcover;
landcover_sampling = landcover_cfg.sampling;
landcover_enabled = logical(landcover_cfg.enabled);

%% ── 1.  Discover DEM / DSM files ─────────────────────────────────────────
all_tif   = dir(fullfile(dem_dir, '**', '*.tif'));
is_dem    = contains({all_tif.name}, '_dem', 'IgnoreCase', true);
is_dsm    = contains({all_tif.name}, '_DSM', 'IgnoreCase', true);
dem_files = all_tif(is_dem | is_dsm);
split_manifest_path = resolve_project_path(dem2phase_dir, ...
    char(sim_cfg.dataset.split_manifest));
split_assignments = load_dem_split_manifest(split_manifest_path, ...
    string({dem_files.name}));

fprintf('\n========================================\n');
fprintf(' V13 Multi-Baseline Dataset Generator\n');
fprintf('========================================\n');
fprintf(' DEM/DSM directory : %s\n', dem_dir);
fprintf(' Files found       : %d\n', numel(dem_files));
fprintf(' Baselines (K=%d)  : %s m\n', K, mat2str(physics_baseline_lengths));
fprintf(' UAVs / edges      : %d / %d\n', sim_cfg.interferometry.num_uavs, K);
fprintf(' Measurement mode  : %s\n', sim_cfg.interferometry.measurement_mode);
fprintf(' Coherence model   : %s\n', sim_cfg.coherence.model);
fprintf(' Terrain sampling  : %s\n', terrain_sampling.mode);
fprintf(' Storage mode      : %s (legacy per-edge: %d)\n', ...
    storage_cfg.mode, write_legacy_per_edge);
fprintf(' Noise model       : %s\n', noise_cfg.model);
fprintf(' Phase-4 errors    : %s (enabled: %d)\n', ...
    phase4_cfg.profile, logical(phase4_cfg.enabled));
fprintf(' Compound failures: %s (enabled: %d)\n', ...
    failure_cfg.profile, logical(failure_cfg.enabled));
fprintf(' Land cover        : %s (enabled: %d)\n', ...
    landcover_cfg.source, landcover_enabled);
fprintf(' Split manifest    : %s\n', split_manifest_path);
fprintf(' Short baseline    : baseline %d (%.3f m)\n', short_baseline_idx, short_baseline_length);
fprintf(' Output directory  : %s\n', output_root);
fprintf('========================================\n');

if numel(dem_files) == 0
    fprintf('[WARN] No *_dem*.tif or *_DSM*.tif files found in:\n  %s\n', dem_dir);
    fprintf('       Place GeoTIFF files there and re-run.\n');
    return
end

for k = 1:numel(dem_files)
    meta = parse_dem_metadata(dem_files(k).name);
    fprintf('  [%02d] %-45s  sensor=%-8s  tile=%-12s  strip=%-4s  type=%s\n', ...
        k, dem_files(k).name, meta.sensor, meta.tile, meta.strip, meta.file_type);
end
fprintf('\n');

%% ── 2.  Report validated per-edge physics ────────────────────────────────
fprintf('Physics parameters:\n');
fprintf('  Center frequency : %.3f MHz\n', sim_cfg.radar.center_frequency_hz / 1e6);
fprintf('  Lambda           : %.6f m\n', physics_lambda);
fprintf('  Altitude         : %.1f m\n', sim_cfg.geometry.altitude_m);
fprintf('  Slant range      : %.1f m\n', physics_slant_range);
fprintf('  Range resolution : %.3f m\n', sim_cfg.physics.range_resolution_m);
fprintf('  Critical baseline: %.3f m\n', sim_cfg.physics.critical_baseline_m);
fprintf('  Incidence angle  : %.2f deg\n', physics_incidence_angle);
fprintf('  Path multiplicity: %d (phase coefficient = %d*2pi/lambda)\n', ...
    phase_path_multiplicity, phase_path_multiplicity);
fprintf('Per-baseline dem2phase ratios:\n');
for k = 1:K
    fprintf('  Edge %d (%6.3f m): ratio = %.6f rad/m, h_amb = %.3f m\n', ...
                k, physics_baseline_lengths(k), dem2phase_ratios(k), ambiguity_heights_m(k));
    fprintf('    Coherence map upper bound <= %.4f\n', max_coh_gen(k));
end
fprintf('\n');

%% ── 3.  Create output sub-directories ────────────────────────────────────
subdirs = {'patch_groups', fullfile('patch_groups', 'train'), ...
    fullfile('patch_groups', 'test'), ...
    fullfile('patch_groups', 'validation')};
if write_legacy_per_edge
    subdirs = [subdirs, {'wrappedphase_withnoise', ...
        'wrappedphase_withoutnoise', 'unwrapped_phase', ...
        'coherence_estimated', 'coherence_true', 'terrain_features'}];
end
for s = 1:numel(subdirs)
    d = fullfile(output_root, subdirs{s});
    if ~exist(d, 'dir'); mkdir(d); end
end
copyfile(sim_config_path, fullfile(output_root, 'simulation_config.json'));
save(fullfile(output_root, 'simulation_config_derived.mat'), 'sim_cfg');
effective_config_fid = fopen(fullfile(output_root, ...
    'simulation_config_effective.json'), 'w');
assert(effective_config_fid >= 0, 'dem2phase:ConfigWriteFailed', ...
    'Could not create simulation_config_effective.json.');
effective_config_cleanup = onCleanup(@() fclose(effective_config_fid));
fprintf(effective_config_fid, '%s', jsonencode(sim_cfg, 'PrettyPrint', true));
clear effective_config_cleanup

%% ── 4.  Open log files (summary + detail) ────────────────────────────────
log_ts = datetime('now', 'Format','yyyy-MM-dd-HH-mm-ss');

% Summary log: one row per source DEM/DSM file
log_path = fullfile(output_root, sprintf('generation_log_%s.csv', log_ts));
log_fid  = fopen(log_path, 'w');
fprintf(log_fid, ['fi,source_file,file_type,sensor,tile,strip,' ...
    'dem_rows,dem_cols,dem_size_bytes,' ...
    'num_available,num_file,' ...
    'patches_saved,attempts_total,' ...
    'void_skipped,wrap_skipped,terrain_skipped,' ...
    'short_baseline_id,short_baseline_length_m,short_wrap_mean,short_wrap_max,' ...
    'num_void_frac_mean,elapsed_s,rng_seed,timestamp\n']);

% Detail log: one row per saved patch
detail_path = fullfile(output_root, sprintf('generation_detail_%s.csv', log_ts));
detail_fid  = fopen(detail_path, 'w');
baseline_wrap_cols = arrayfun(@(k) sprintf('max_wrap_count_baseline%d', k), 1:K, 'UniformOutput', false);
detail_wrap_cols_csv = strjoin(baseline_wrap_cols, ',');
fprintf(detail_fid, ['patch_global_id,patch_name,source_file,file_type,sensor,' ...
    'tile,strip,dem_rows,dem_cols,dem_size_bytes,' ...
    'interp_scale,scaled_rows,scaled_cols,' ...
    'crop_r0,crop_c0,crop_r1,crop_c1,' ...
    'max_wrap_count_all_baselines,short_baseline_id,short_baseline_length_m,' ...
    'short_baseline_max_wrap_count,short_baseline_mean_wrap_count,' ...
    detail_wrap_cols_csv, ',active_edge_count,valid_edge_mask,mean_node_snr_db,' ...
    'coherence_seed,slc_seed,phase4_error_seed,edge_mask_seed,failure_seed,' ...
    'accept_probability,dropout_present,sync_anomaly_present,' ...
    'low_coherence_fraction,compound_failure_present,' ...
    'landcover_sampling_multiplier,priority_landcover_code,' ...
    'terrain_class,ground_pixel_spacing_m,' ...
    'mean_slope_deg,p95_slope_deg,layover_fraction,shadow_fraction,mean_terrain_quality,' ...
    'num_available,num_void_frac,num_mean,num_min,' ...
    'patch_size,timestamp\n']);

% File-detail log: one row per saved baseline sample (with concrete filenames)
file_detail_path = fullfile(output_root, sprintf('generation_file_detail_%s.csv', log_ts));
file_detail_fid  = fopen(file_detail_path, 'w');
fprintf(file_detail_fid, ['patch_global_id,patch_name,baseline_idx,baseline_length_m,' ...
    'source_file,valid_edge,max_wrap_count,mean_wrap_count,patch_group_file,' ...
    'wrappedphase_withnoise_file,wrappedphase_withoutnoise_file,unwrapped_phase_file,' ...
    'coherence_estimated_file,coherence_true_file,' ...
    'timestamp\n']);

fprintf('Summary log : %s\n', log_path);
fprintf('Detail  log : %s\n\n', detail_path);
fprintf('File-detail log : %s\n\n', file_detail_path);

%% ── 5.  Main loop: DEM/DSM → patches ────────────────────────────────────
total_saved     = 0;
total_skip      = 0;
total_void_skip = 0;
total_terrain_skip = 0;
total_short_wrap_max_sum = 0;
total_short_wrap_count = 0;
total_short_wrap_max = 0;

% Set random seed and record the actual seed used for reproducibility.
rng(rng_seed);
rng_state      = rng;
rng_seed_used  = rng_state.Seed;
component_seed_counter = 0;
rng_seed_path  = fullfile(output_root, sprintf('rng_seed_%s.mat', log_ts));
save(rng_seed_path, 'rng_state', 'rng_seed_used');
fprintf('RNG seed     : %d  (saved to %s)\n\n', rng_seed_used, rng_seed_path);

tic
for fi = 1:numel(dem_files)
    dem_path  = fullfile(dem_files(fi).folder, dem_files(fi).name);
    [~, dem_stem, ~] = fileparts(dem_files(fi).name);
    meta = parse_dem_metadata(dem_files(fi).name);
    split_idx = find(strcmpi(split_assignments.dem_file, ...
        dem_files(fi).name));
    assert(isscalar(split_idx), 'dem2phase:InvalidDemSplitLookup', ...
        'Expected one split assignment for %s.', dem_files(fi).name);
    dataset_split = char(split_assignments.split(split_idx));
    geographic_region = char(split_assignments.region(split_idx));
    geographic_tile = char(split_assignments.geographic_tile(split_idx));
    source_quality_role = char(split_assignments.quality_role(split_idx));
    source_dataset = char(split_assignments.source(split_idx));
    patches_target = double(split_assignments.patches_per_dem(split_idx));

    fprintf('──────────────────────────────────────\n');
    fprintf('[%02d/%02d] %s\n', fi, numel(dem_files), dem_files(fi).name);
    fprintf('         sensor=%s  tile=%s  strip=%s  type=%s\n', ...
        meta.sensor, meta.tile, meta.strip, meta.file_type);
    fprintf('         split=%s  region=%s  geotile=%s  quality=%s  target=%d patches\n', ...
        dataset_split, geographic_region, geographic_tile, ...
        source_quality_role, patches_target);

    % ── 5a. Read file ─────────────────────────────────────────────────────
    dem_raw = double(imread(dem_path));
    [dem_rows, dem_cols] = size(dem_raw);
    [landcover_full, file_landcover_enabled, landcover_path] = ...
        resolve_landcover_for_dem(dem2phase_dir, dem_stem, ...
        landcover_cfg, landcover_enabled);
    if file_landcover_enabled
        assert(isequal(size(landcover_full), [dem_rows, dem_cols]), ...
            'dem2phase:LandcoverDemShapeMismatch', ...
            'Aligned land-cover raster must match DEM %s.', dem_files(fi).name);
        fprintf('         Landcover: %s\n', landcover_path);
    elseif landcover_enabled
        fprintf('         Landcover: neutral fallback (aligned map unavailable)\n');
    end
    fprintf('         Size    : %d x %d px  (%s)\n', dem_rows, dem_cols, ...
        format_bytes(dem_files(fi).bytes));

    % ── 5a'. Load NUM file (ASTGTM stacking-count, optional) ──────────────
    num_map      = [];
    num_file     = '';
    num_avail    = 0;
    num_path_try = get_num_path(dem_path);
    if ~isempty(num_path_try)
        num_map   = uint8(imread(num_path_try));
        num_file  = num_path_try;
        num_avail = 1;
        fprintf('         NUM     : %s\n', num_path_try);
    end

    % ── 5b. Random patch extraction ───────────────────────────────────────
    patch_count         = 0;
    attempt             = 0;
    max_attempts        = patches_target * double(terrain_sampling.max_attempts_factor);
    void_skipped        = 0;
    wrap_skipped        = 0;
    terrain_skipped     = 0;
    num_void_frac_sum   = 0;
    num_void_frac_count = 0;
    short_wrap_max_sum  = 0;
    short_wrap_count    = 0;
    short_wrap_max_file = 0;
    wrap_max_sum_per_baseline = zeros(1, K);
    file_tic            = tic;

    % ── 5b'. Pre-compute scaled DEM images
    precomp_sc   = linspace(interp_scale_range(1), interp_scale_range(2), n_precomp_scales);
    precomp_imgs = cell(n_precomp_scales, 1);
    [xg_base, yg_base] = meshgrid(1:dem_cols, 1:dem_rows);
    for si = 1:n_precomp_scales
        sc = precomp_sc(si);
        rc = round(sc * dem_rows);
        cc = round(sc * dem_cols);
        if rc < patch_size + 1 || cc < patch_size + 1
            precomp_imgs{si} = dem_raw;
            precomp_sc(si)   = 1.0;
        else
            xlin = linspace(1, dem_cols, cc);
            ylin = linspace(1, dem_rows, rc);
            [xxg, yyg] = meshgrid(xlin, ylin);
            precomp_imgs{si} = interp2(xg_base, yg_base, dem_raw, xxg, yyg, 'bicubic');
        end
    end
    fprintf('         Scaled  : %d pre-computed images at scales [%s]\n', ...
        n_precomp_scales, num2str(precomp_sc, '%.2f '));

    while patch_count < patches_target && attempt < max_attempts
        attempt = attempt + 1;

        % Pick one of the pre-computed scaled images at random
        si          = randi(n_precomp_scales);
        scale       = precomp_sc(si);
        dem_interp  = precomp_imgs{si};
        [rows_scaled, cols_scaled] = size(dem_interp);

        % Random crop
        [mr, mc] = size(dem_interp);
        r0 = randi([1, max(1, mr - patch_size)]);
        c0 = randi([1, max(1, mc - patch_size)]);
        r1 = r0 + patch_size - 1;
        c1 = c0 + patch_size - 1;

        if r1 > mr || c1 > mc
            continue
        end

        patch_dem = dem_interp(r0:r1, c0:c1);

        if size(patch_dem,1) ~= patch_size || size(patch_dem,2) ~= patch_size
            continue
        end

        % Normalise DC offset
        patch_dem = patch_dem - min(patch_dem(:));

        if file_landcover_enabled
            source_rows = round(1 + ((r0:r1)-1) * ...
                (dem_rows-1)/max(rows_scaled-1, 1));
            source_cols = round(1 + ((c0:c1)-1) * ...
                (dem_cols-1)/max(cols_scaled-1, 1));
            patch_landcover = landcover_full(source_rows, source_cols);
        else
            patch_landcover = zeros(size(patch_dem), 'uint8');
        end

        % ── 5c. Generate per-baseline clean wrapped phases ────────────────
        % Each baseline k has its own unwrapped phase:
        %   unw_k = patch_dem * ratio_k
        % and wrapped phase:
        %   wrap_k = angle(exp(1i * unw_k))

        clean_phases_unw = cell(1, K);
        clean_phases_wrap = cell(1, K);
        max_wrap_all = 0;
        max_wrap_per_baseline = zeros(1, K);
        mean_wrap_per_baseline = zeros(1, K);

        for k = 1:K
            unw_k = patch_dem * dem2phase_ratios(k);
            wrap_k = angle(exp(1i * unw_k));
            clean_phases_unw{k} = unw_k;
            clean_phases_wrap{k} = wrap_k;

            % Compute wrap count for this baseline
            wrap_count_k = round((unw_k - wrap_k) / (2*pi));
            max_wrap_k = max(wrap_count_k(:));
            mean_wrap_k = mean(wrap_count_k(:));
            max_wrap_per_baseline(k) = max_wrap_k;
            mean_wrap_per_baseline(k) = mean_wrap_k;
            if max_wrap_k > max_wrap_all
                max_wrap_all = max_wrap_k;
            end
        end

        short_wrap_max_patch = max_wrap_per_baseline(short_baseline_idx);
        short_wrap_mean_patch = mean_wrap_per_baseline(short_baseline_idx);

        % Wrap count filter: reject if ANY baseline exceeds max_wrap_count
        if max_wrap_all > max_wrap_count
            wrap_skipped = wrap_skipped + 1;
            total_skip   = total_skip + 1;
            continue
        end

        % NUM quality filter (only when _num.tif was found alongside the DEM)
        num_void_frac = NaN;
        num_mean_val  = NaN;
        num_min_val   = NaN;
        if ~isempty(num_map)
            scale_r = (dem_rows - 1) / max(rows_scaled - 1, 1);
            scale_c = (dem_cols - 1) / max(cols_scaled - 1, 1);
            nr0 = max(1,        round(1 + (r0 - 1) * scale_r));
            nr1 = min(dem_rows, round(1 + (r1 - 1) * scale_r));
            nc0 = max(1,        round(1 + (c0 - 1) * scale_c));
            nc1 = min(dem_cols, round(1 + (c1 - 1) * scale_c));
            if nr1 >= nr0 && nc1 >= nc0
                num_crop      = double(num_map(nr0:nr1, nc0:nc1));
                num_void_frac = mean(num_crop(:) == 0);
                num_mean_val  = mean(num_crop(:));
                num_min_val   = min(num_crop(:));
                if num_void_frac > void_frac_thresh
                    void_skipped        = void_skipped + 1;
                    total_skip          = total_skip + 1;
                    total_void_skip     = total_void_skip + 1;
                    continue
                end
                num_void_frac_sum   = num_void_frac_sum   + num_void_frac;
                num_void_frac_count = num_void_frac_count + 1;
            end
        end

        % ── 5d. Terrain-conditioned coherence and sampling ────────────────
        % Bicubic up-scaling reduces the effective ground spacing.
        ground_pixel_spacing_m = source_dem_pixel_spacing_m / scale;
        component_seed_counter = component_seed_counter + 1;
        coherence_seed = derive_component_seed(rng_seed_used, component_seed_counter);
        [coh_maps_input, terrain_features] = ...
            gen_coherence_map_multibaseline_terrainaware( ...
                patch_dem, physics_baseline_lengths, ...
                'pixel_spacing_m', ground_pixel_spacing_m, ...
                'incidence_angle_deg', physics_incidence_angle, ...
                'look_azimuth_deg', look_azimuth_deg, ...
                'min_coh', coh_min, 'max_coh', coh_max, ...
                'baseline_ref', coh_ref_baseline, ...
                'baseline_decay_model', char(sim_cfg.coherence.baseline_decay_model), ...
                'critical_baseline_m', double(sim_cfg.physics.critical_baseline_m), ...
                'decay_alpha', coh_decay_alpha, ...
                'terrain_normalization_model', ...
                    char(sim_cfg.coherence.terrain_normalization_model), ...
                'slope_scale_deg', double(sim_cfg.coherence.slope_scale_deg), ...
                'roughness_scale_m', double(sim_cfg.coherence.roughness_scale_m), ...
                'curvature_scale_per_m', ...
                    double(sim_cfg.coherence.curvature_scale_per_m), ...
                'slope_weight', double(sim_cfg.coherence.slope_weight), ...
                'roughness_weight', double(sim_cfg.coherence.roughness_weight), ...
                'curvature_weight', double(sim_cfg.coherence.curvature_weight), ...
                'shared_residual_std', double(sim_cfg.coherence.shared_residual_std), ...
                'edge_residual_std', double(sim_cfg.coherence.edge_residual_std), ...
                'residual_scale_px', double(sim_cfg.coherence.residual_scale_px), ...
                'seed', coherence_seed);

        coh_maps_terrain = coh_maps_input;
        if file_landcover_enabled
            [coh_maps_input, landcover_info] = apply_landcover_coherence( ...
                coh_maps_input, patch_landcover, landcover_cfg, coh_min);
        else
            landcover_info.codes = patch_landcover;
            landcover_info.factor_map = ones(size(patch_dem));
            landcover_info.dominant_code = 0;
            landcover_info.water_fraction = 0;
            if landcover_enabled
                landcover_info.source = 'aligned-map-missing-neutral';
                landcover_info.profile = char(landcover_cfg.profile);
            else
                landcover_info.source = 'disabled';
                landcover_info.profile = 'none';
            end
            landcover_info.factors_are_calibrated = false;
        end

        terrain_summary = summarize_terrain(terrain_features);
        accept_probability = terrain_acceptance_probability( ...
            terrain_summary.class_name, terrain_sampling);
        [landcover_multiplier, priority_landcover_code] = ...
            landcover_sampling_multiplier(patch_landcover, landcover_sampling);
        accept_probability = min(1, accept_probability * landcover_multiplier);
        if rand() > accept_probability
            terrain_skipped = terrain_skipped + 1;
            total_terrain_skip = total_terrain_skip + 1;
            total_skip = total_skip + 1;
            continue
        end

        patch_count = patch_count + 1;
        total_saved = total_saved + 1;
        global_id   = total_saved;
        short_wrap_max_sum = short_wrap_max_sum + short_wrap_max_patch;
        short_wrap_count = short_wrap_count + 1;
        short_wrap_max_file = max(short_wrap_max_file, short_wrap_max_patch);
        total_short_wrap_max_sum = total_short_wrap_max_sum + short_wrap_max_patch;
        total_short_wrap_count = total_short_wrap_count + 1;
        total_short_wrap_max = max(total_short_wrap_max, short_wrap_max_patch);
        wrap_max_sum_per_baseline = wrap_max_sum_per_baseline + max_wrap_per_baseline;

        patch_name = sprintf('%s_patch_%05d', dem_stem, global_id);

        % ── 5e. Generate distributed observation noise ───────────────────
        component_seed_counter = component_seed_counter + 1;
        slc_seed = derive_component_seed(rng_seed_used, component_seed_counter);
        component_seed_counter = component_seed_counter + 1;
        phase4_error_seed = derive_component_seed(rng_seed_used, component_seed_counter);
        component_seed_counter = component_seed_counter + 1;
        edge_mask_seed = derive_component_seed(rng_seed_used, component_seed_counter);
        component_seed_counter = component_seed_counter + 1;
        failure_seed = derive_component_seed(rng_seed_used, component_seed_counter);
        if strcmp(char(noise_cfg.model), 'complex_slc_nodes')
            snr_range = double(noise_cfg.node_snr_db_range(:)');
            node_snr_db = snr_range(1) + diff(snr_range) * ...
                rand(1, double(sim_cfg.interferometry.num_uavs));
            [noisy_phases, node_slc, coh_maps_true, noise_info] = ...
                simulate_distributed_slc_star(clean_phases_wrap, ...
                    coh_maps_input, sim_cfg.interferometry.edge_index, ...
                    double(sim_cfg.interferometry.num_uavs), ...
                    'node_snr_db', node_snr_db, 'seed', slc_seed);
            noisy_phases_node_noise_only = noisy_phases;
            if logical(phase4_cfg.enabled)
                [node_slc, noisy_phases, edge_coreg_valid, phase4_info] = ...
                    apply_distributed_uav_errors(node_slc, ...
                        sim_cfg.interferometry.edge_index, phase4_cfg, ...
                        'wavelength_m', physics_lambda, ...
                        'phase_path_multiplicity', phase_path_multiplicity, ...
                        'altitude_m', double(sim_cfg.geometry.altitude_m), ...
                        'ground_pixel_spacing_m', ground_pixel_spacing_m, ...
                        'seed', phase4_error_seed);
            else
                edge_coreg_valid = true(K, patch_size, patch_size);
                phase4_info.profile = 'disabled';
                phase4_info.seed = phase4_error_seed;
            end
            multilook_phases = multilook_node_interferograms(node_slc, ...
                sim_cfg.interferometry.edge_index, double(noise_cfg.multilook_window));
        else
            noisy_phases = patch_add_noise_multibaseline( ...
                clean_phases_wrap, coh_maps_input, 'base_seed', slc_seed);
            node_slc = [];
            coh_maps_true = coh_maps_input;
            node_snr_db = NaN;
            noise_info.model = 'legacy_phase_pdf';
            noise_info.node_snr_db = node_snr_db;
            noisy_phases_node_noise_only = noisy_phases;
            edge_coreg_valid = true(K, patch_size, patch_size);
            phase4_info.profile = 'unsupported_legacy_noise_path';
            phase4_info.seed = phase4_error_seed;
            multilook_phases = phase_only_multilook(noisy_phases, ...
                double(noise_cfg.multilook_window));
        end

        % ── 5f. Estimate deployable coherence from observed node SLCs ────
        coh_maps_est = cell(1, K);
        for k = 1:K
            if strcmp(char(noise_cfg.model), 'complex_slc_nodes')
                master_slc = squeeze(node_slc( ...
                    sim_cfg.interferometry.edge_index(1,k),:,:));
                secondary_slc = squeeze(node_slc( ...
                    sim_cfg.interferometry.edge_index(2,k),:,:));
                [coh_est_k, ~] = calc_coherence(master_slc, secondary_slc, ...
                    'method', 'slc', 'win_size', coh_est_win);
                coherence_observation_method = 'observed_node_slc_pair';
            else
                [coh_est_k, ~] = calc_coherence( ...
                    noisy_phases{k}, clean_phases_wrap{k}, ...
                    'method', 'phase', 'win_size', coh_est_win);
                coherence_observation_method = 'legacy_clean_phase_residual';
            end
            coh_maps_est{k} = coh_est_k;
        end

        % ── 5g. Select observed edges and save one [K,H,W] patch group ───
        initial_edge_mask = sample_valid_edge_mask( ...
            K, short_baseline_idx, edge_sampling, edge_mask_seed);
        if strcmp(char(noise_cfg.model), 'complex_slc_nodes')
            [valid_edge_mask, failure_labels] = apply_compound_failure_mask( ...
                initial_edge_mask, sim_cfg.interferometry.edge_index, ...
                physics_baseline_lengths, coh_maps_true, phase4_info, ...
                failure_cfg, failure_seed);
        else
            valid_edge_mask = initial_edge_mask;
            failure_labels.profile = 'unsupported_legacy_noise_path';
            failure_labels.enabled = false;
            failure_labels.dropout_present = false;
            failure_labels.sync_anomaly_present = false;
            failure_labels.low_coherence_fraction = NaN;
            failure_labels.compound_failure_present = false;
        end
        valid_edge_mask_text = sprintf('%d', valid_edge_mask);
        group_filename = sprintf('%s.mat', patch_name);
        group_relative_path = fullfile('patch_groups', dataset_split, ...
            group_filename);
        group_metadata.patch_global_id = global_id;
        group_metadata.patch_name = patch_name;
        group_metadata.source_file = dem_files(fi).name;
        group_metadata.dataset_split = dataset_split;
        group_metadata.split_group = geographic_tile;
        group_metadata.geographic_tile = geographic_tile;
        group_metadata.geographic_region = geographic_region;
        group_metadata.source_dataset = source_dataset;
        group_metadata.source_quality_role = source_quality_role;
        group_metadata.num_available = logical(num_avail);
        group_metadata.num_void_fraction = num_void_frac;
        group_metadata.num_mean = num_mean_val;
        group_metadata.num_min = num_min_val;
        group_metadata.file_type = meta.file_type;
        group_metadata.sensor = meta.sensor;
        group_metadata.tile = meta.tile;
        group_metadata.strip = meta.strip;
        group_metadata.interp_scale = scale;
        group_metadata.crop_bounds_rc = [r0, c0, r1, c1];
        group_metadata.ground_pixel_spacing_m = ground_pixel_spacing_m;
        group_metadata.terrain_class = terrain_summary.class_name;
        group_metadata.noise = noise_info;
        group_metadata.landcover_source = landcover_info.source;
        group_metadata.landcover_profile = landcover_info.profile;
        group_metadata.dominant_landcover_code = landcover_info.dominant_code;
        group_metadata.water_fraction = landcover_info.water_fraction;
        group_metadata.rng_batch_seed = rng_seed_used;
        group_metadata.rng_component_seeds = struct( ...
            'coherence', coherence_seed, 'slc', slc_seed, ...
            'phase4_error', phase4_error_seed, 'edge_mask', edge_mask_seed, ...
            'failure', failure_seed);
        group_metadata.phase4_error_profile = phase4_info.profile;
        group_metadata.failure_labels = failure_labels;
        group_metadata.coherence_observation_method = ...
            coherence_observation_method;
        group_metadata.accept_probability = accept_probability;
        group_metadata.landcover_sampling_multiplier = landcover_multiplier;
        group_metadata.priority_landcover_code = priority_landcover_code;
        group_metadata.timestamp = char(datetime('now', ...
            'Format','yyyy-MM-dd HH:mm:ss.SSS'));
        patch_group = build_patch_group(clean_phases_wrap, noisy_phases, ...
            clean_phases_unw, coh_maps_est, coh_maps_true, valid_edge_mask, ...
            sim_cfg, terrain_features, group_metadata);
        patch_group.coherence_scene = stack_edge_cells( ...
            coh_maps_input, valid_edge_mask, char(storage_cfg.numeric_type));
        patch_group.coherence_terrain_only = stack_edge_cells( ...
            coh_maps_terrain, valid_edge_mask, char(storage_cfg.numeric_type));
        patch_group.coherence_observed = patch_group.coherence_estimated;
        patch_group.landcover_codes = uint8(patch_landcover);
        patch_group.landcover_factor = cast(landcover_info.factor_map, ...
            char(storage_cfg.numeric_type));
        patch_group.wrappedphase_multilook = stack_edge_cells( ...
            multilook_phases, valid_edge_mask, char(storage_cfg.numeric_type));
        patch_group.wrappedphase_node_noise_only = stack_edge_cells( ...
            noisy_phases_node_noise_only, valid_edge_mask, char(storage_cfg.numeric_type));
        patch_group.coregistration_valid_mask = uint8(edge_coreg_valid);
        patch_group.coregistration_valid_mask(~valid_edge_mask,:,:) = 0;
        if logical(phase4_cfg.enabled) && strcmp(char(noise_cfg.model), 'complex_slc_nodes')
            patch_group.phase4_errors = phase4_info;
        end
        patch_group.failure_labels = failure_labels;
        if logical(noise_cfg.save_node_slc) && ~isempty(node_slc)
            patch_group.node_slc = cast(node_slc, char(storage_cfg.numeric_type));
        end
        save_patch_group(fullfile(output_root, group_relative_path), patch_group);

        % Optional legacy per-edge data and one provenance row per edge.
        for k = 1:K
            basename = sprintf('%s_baseline%d', patch_name, k);
            fn_wrapped_withnoise = '';
            fn_wrapped_withoutnoise = '';
            fn_unwrapped_phase = '';
            fn_coherence_estimated = '';
            fn_coherence_true = '';
            if write_legacy_per_edge
                fn_wrapped_withnoise = sprintf('%s_wrappedphase_withnoise.mat', basename);
                fn_wrapped_withoutnoise = sprintf('%s_wrappedphase_withoutnoise.mat', basename);
                fn_unwrapped_phase = sprintf('%s_unwrapped_phase.mat', basename);
                fn_coherence_estimated = sprintf('%s_coherence_estimated.mat', basename);
                fn_coherence_true = sprintf('%s_coherence_true.mat', basename);
                save_mat(output_root, 'wrappedphase_withnoise', basename, ...
                    'wrappedphase_withnoise', noisy_phases{k});
                save_mat(output_root, 'wrappedphase_withoutnoise', basename, ...
                    'wrappedphase_withoutnoise', clean_phases_wrap{k});
                save_mat(output_root, 'unwrapped_phase', basename, ...
                    'unwrapped_phase', clean_phases_unw{k});
                save_mat(output_root, 'coherence_estimated', basename, ...
                    'coherence_estimated', coh_maps_est{k});
                save_mat(output_root, 'coherence_true', basename, ...
                    'coherence_true', coh_maps_true{k});
            end

            % File-detail log entry (one row per saved baseline sample)
            fprintf(file_detail_fid, ...
                '%d,%s,%d,%.6f,%s,%d,%.0f,%.4f,%s,%s,%s,%s,%s,%s,%s\n', ...
                global_id, patch_name, k, physics_baseline_lengths(k), dem_files(fi).name, ...
                valid_edge_mask(k), max_wrap_per_baseline(k), mean_wrap_per_baseline(k), ...
                group_relative_path, ...
                fn_wrapped_withnoise, fn_wrapped_withoutnoise, fn_unwrapped_phase, ...
                fn_coherence_estimated, fn_coherence_true, ...
                datetime('now', 'Format','yyyy-MM-dd HH:mm:ss.SSS'));
        end

        if write_legacy_per_edge
            save_mat(output_root, 'terrain_features', patch_name, ...
                'terrain_features', terrain_features);
        end

        % ── 5h. Write detail log entry ────────────────────────────────────
        fprintf(detail_fid, '%d,%s,%s,%s,%s,%s,%s,%d,%d,%d,%.4f,%d,%d,%d,%d,%d,%d,%d,%d,%.6f,%d,%.4f,', ...
            global_id, patch_name, dem_files(fi).name, meta.file_type, meta.sensor, ...
            meta.tile, meta.strip, ...
            dem_rows, dem_cols, dem_files(fi).bytes, ...
            scale, rows_scaled, cols_scaled, ...
            r0, c0, r1, c1, ...
            max_wrap_all, short_baseline_idx, short_baseline_length, ...
            short_wrap_max_patch, short_wrap_mean_patch);
        for kb = 1:K
            fprintf(detail_fid, '%.0f,', max_wrap_per_baseline(kb));
        end
        fprintf(detail_fid, ['%d,%s,%.4f,%u,%u,%u,%u,%u,%.4f,%d,%d,%.6f,%d,' ...
            '%.4f,%d,%s,' ...
            '%.4f,%.4f,%.4f,%.6f,%.6f,%.4f,%d,%.4f,%.4f,%.4f,%d,%s\n'], ...
            nnz(valid_edge_mask), valid_edge_mask_text, mean(node_snr_db), ...
            coherence_seed, slc_seed, phase4_error_seed, edge_mask_seed, failure_seed, ...
            accept_probability, failure_labels.dropout_present, ...
            failure_labels.sync_anomaly_present, ...
            failure_labels.low_coherence_fraction, ...
            failure_labels.compound_failure_present, ...
            landcover_multiplier, priority_landcover_code, ...
            terrain_summary.class_name, ground_pixel_spacing_m, ...
            terrain_summary.mean_slope_deg, terrain_summary.p95_slope_deg, ...
            terrain_summary.layover_fraction, terrain_summary.shadow_fraction, ...
            terrain_summary.mean_terrain_quality, ...
            num_avail, num_void_frac, num_mean_val, num_min_val, ...
            patch_size, datetime('now', 'Format','yyyy-MM-dd HH:mm:ss.SSS'));

        if mod(patch_count, 50) == 0
            fprintf('         Saved %d / %d  (skipped: %d)\n', ...
                patch_count, patches_target, total_skip);
        end
    end % while

    if short_wrap_count > 0
        short_wrap_mean_file = short_wrap_max_sum / short_wrap_count;
    else
        short_wrap_mean_file = NaN;
    end
    fprintf(['         Done : %d patches  (skipped: %d wrap, %d void, ' ...
        '%d terrain-weighted, attempts: %d)\n'], ...
        patch_count, wrap_skipped, void_skipped, terrain_skipped, attempt);
    fprintf('         Short baseline (B%d=%.3f m) wrap max/mean: %.0f / %.2f\n', ...
        short_baseline_idx, short_baseline_length, short_wrap_max_file, short_wrap_mean_file);
    for k = 1:K
        if patch_count > 0
            mean_wrap_k_file = wrap_max_sum_per_baseline(k) / patch_count;
        else
            mean_wrap_k_file = NaN;
        end
        fprintf('         Baseline %d (%.3f m) max-wrap mean: %.2f\n', ...
            k, physics_baseline_lengths(k), mean_wrap_k_file);
    end

    % Write summary log row (one per source file)
    file_elapsed = toc(file_tic);
    if num_void_frac_count > 0
        num_void_frac_mean = num_void_frac_sum / num_void_frac_count;
    else
        num_void_frac_mean = NaN;
    end
    fprintf(log_fid, '%d,%s,%s,%s,%s,%s,%d,%d,%d,%d,%s,%d,%d,%d,%d,%d,%d,%.6f,%.4f,%.0f,%.4f,%.1f,%d,%s\n', ...
        fi, dem_files(fi).name, meta.file_type, meta.sensor, meta.tile, meta.strip, ...
        dem_rows, dem_cols, dem_files(fi).bytes, ...
        num_avail, num_file, ...
        patch_count, attempt, void_skipped, wrap_skipped, terrain_skipped, ...
        short_baseline_idx, short_baseline_length, short_wrap_mean_file, short_wrap_max_file, ...
        num_void_frac_mean, file_elapsed, rng_seed_used, ...
        datetime('now', 'Format','yyyy-MM-dd HH:mm:ss.SSS'));
end % for fi

fclose(log_fid);
fclose(detail_fid);
fclose(file_detail_fid);
rng_final_state = rng;
save(rng_seed_path, 'rng_final_state', 'component_seed_counter', '-append');

elapsed = toc;
fprintf('\n========================================\n');
fprintf(' Total patches saved : %d\n', total_saved);
fprintf(' Total baselines/patch: %d\n', K);
fprintf(' Total samples        : %d\n', total_saved * K);
fprintf(' Skipped (all)       : %d\n', total_skip);
fprintf('   wrap (>%d)        : %d\n', max_wrap_count, ...
    total_skip - total_void_skip - total_terrain_skip);
fprintf('   void-fill         : %d\n', total_void_skip);
fprintf('   terrain-weighted  : %d\n', total_terrain_skip);
if total_short_wrap_count > 0
    total_short_wrap_mean = total_short_wrap_max_sum / total_short_wrap_count;
else
    total_short_wrap_mean = NaN;
end
fprintf(' Short baseline      : B%d=%.3f m\n', short_baseline_idx, short_baseline_length);
fprintf('   short-wrap max    : %.0f\n', total_short_wrap_max);
fprintf('   short-wrap mean   : %.2f\n', total_short_wrap_mean);
fprintf(' Elapsed time        : %.1f s\n', elapsed);
fprintf(' RNG seed used       : %d\n', rng_seed_used);
fprintf(' Output directory    : %s\n', output_root);
fprintf(' Summary log         : %s\n', log_path);
fprintf(' Detail  log         : %s\n', detail_path);
fprintf(' File detail log     : %s\n', file_detail_path);
fprintf(' RNG seed file       : %s\n', rng_seed_path);
fprintf('========================================\n');

%% ── Helper functions ──────────────────────────────────────────────────────

function path_out = resolve_project_path(project_root, configured_path)
%RESOLVE_PROJECT_PATH Resolve absolute paths or paths relative to dem2phase/.
    is_windows_absolute = ~isempty(regexp(configured_path, ...
        '^[A-Za-z]:[\\/]', 'once')) || startsWith(configured_path, '\\');
    if is_windows_absolute || startsWith(configured_path, '/')
        path_out = configured_path;
    else
        path_out = fullfile(project_root, configured_path);
    end
end

% -------------------------------------------------------------------------
function [codes, enabled_for_file, aligned_path] = resolve_landcover_for_dem( ...
        project_root, dem_stem, landcover, globally_enabled)
%RESOLVE_LANDCOVER_FOR_DEM Never reuse one tile's land cover for another DEM.
    codes = [];
    enabled_for_file = false;
    aligned_path = '';
    if ~globally_enabled
        return
    end
    aligned_dir = resolve_project_path(project_root, ...
        char(landcover.aligned_directory));
    aligned_path = fullfile(aligned_dir, ...
        [dem_stem char(landcover.aligned_suffix)]);
    if exist(aligned_path, 'file')
        data = load(aligned_path, 'landcover_codes');
        assert(isfield(data, 'landcover_codes'), ...
            'dem2phase:MissingLandcoverCodes', ...
            'Aligned land-cover file lacks landcover_codes: %s', aligned_path);
        codes = data.landcover_codes;
        enabled_for_file = true;
        return
    end
    if strcmp(char(landcover.missing_policy), 'error')
        error('dem2phase:MissingLandcoverForDem', ...
            'No aligned land-cover map exists for DEM %s.', dem_stem);
    end
    aligned_path = '';
end

% -------------------------------------------------------------------------
function summary = summarize_terrain(terrain)
%SUMMARIZE_TERRAIN Assign a physical terrain stratum and useful diagnostics.
    slope_values = terrain.slope_deg(isfinite(terrain.slope_deg));
    slope_values = sort(slope_values(:));
    if isempty(slope_values)
        mean_slope = NaN;
        p95_slope = NaN;
    else
        mean_slope = mean(slope_values);
        p95_idx = max(1, min(numel(slope_values), ceil(0.95*numel(slope_values))));
        p95_slope = slope_values(p95_idx);
    end
    layover_fraction = mean(terrain.layover_mask(:));
    shadow_fraction = mean(terrain.shadow_mask(:));
    invalid_fraction = mean(terrain.invalid_geometry_mask(:));

    if invalid_fraction >= 0.01
        class_name = 'geometry_hazard';
    elseif p95_slope < 5
        class_name = 'flat';
    elseif p95_slope < 20
        class_name = 'rolling';
    else
        class_name = 'steep';
    end

    summary.class_name = class_name;
    summary.mean_slope_deg = mean_slope;
    summary.p95_slope_deg = p95_slope;
    summary.layover_fraction = layover_fraction;
    summary.shadow_fraction = shadow_fraction;
    summary.mean_terrain_quality = mean(terrain.terrain_quality(:));
end

% -------------------------------------------------------------------------
function probability = terrain_acceptance_probability(class_name, sampling)
%TERRAIN_ACCEPTANCE_PROBABILITY Weighted random sampling within real strata.
    mode = char(sampling.mode);
    if strcmp(mode, 'uniform')
        probability = 1;
        return
    end
    assert(strcmp(mode, 'weighted'), ...
        'dem2phase:UnsupportedTerrainSamplingMode', ...
        'terrain_sampling.mode must be uniform or weighted.');
    probabilities = sampling.acceptance_probability;
    assert(isfield(probabilities, class_name), ...
        'dem2phase:MissingTerrainAcceptanceProbability', ...
        'No acceptance probability is configured for terrain class %s.', class_name);
    probability = double(probabilities.(class_name));
    assert(isfinite(probability) && probability > 0 && probability <= 1, ...
        'dem2phase:InvalidTerrainAcceptanceProbability', ...
        'Terrain acceptance probabilities must lie in (0, 1].');
end

function seed = derive_component_seed(batch_seed, counter)
%DERIVE_COMPONENT_SEED Create deterministic, distinct local RNG substreams.
    modulus = 2147483646;
    seed = mod(double(batch_seed) + 104729*double(counter), modulus) + 1;
end

function save_mat(root, subdir, patch_name, patch_type, data)
%SAVE_MAT  Save data to a .mat file named <patch_name>_<patch_type>.mat
    filename = sprintf('%s_%s.mat', patch_name, patch_type);
    save(fullfile(root, subdir, filename), 'data');
end

% -------------------------------------------------------------------------
function meta = parse_dem_metadata(filename)
%PARSE_DEM_METADATA  Extract sensor, tile, strip and file-type from a filename.
    meta.sensor    = 'unknown';
    meta.tile      = 'N/A';
    meta.strip     = 'N/A';
    meta.file_type = 'DEM';

    [~, stem, ~] = fileparts(filename);
    stem_up = upper(stem);

    % ── File type ──────────────────────────────────────────────────────────
    if contains(stem_up, '_DSM')
        meta.file_type = 'DSM';
    else
        meta.file_type = 'DEM';
    end

    % ── Sensor ─────────────────────────────────────────────────────────────
    sensors = {'COPERNICUS', 'ASTGTM', 'SRTM', 'ALPSMLC', 'ALOS2', ...
        'ALOS', 'COPDEM', 'TDX', 'TANDEM'};
    for s = 1:numel(sensors)
        if contains(stem_up, upper(sensors{s}))
            meta.sensor = sensors{s};
            break
        end
    end

    tok = regexp(stem, 'system_(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        meta.sensor = sprintf('system_%s', tok{1});
    end

    % ── Tile (lat/lon grid cell) ───────────────────────────────────────────
    tok = regexp(stem, '[NnSs]\d{2,3}[EeWw]\d{3}', 'match', 'once');
    if isempty(tok)
        copernicus_tok = regexp(stem, ...
            '([NnSs]\d{2})_00_([EeWw]\d{3})_00', 'tokens', 'once');
        if ~isempty(copernicus_tok)
            tok = [copernicus_tok{1} copernicus_tok{2}];
        end
    end
    if ~isempty(tok)
        meta.tile = upper(tok);
    end

    % ── Strip / track index ───────────────────────────────────────────────
    tok = regexp(stem, '(?:DSM|dsm|track|Track|_T)_?(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        meta.strip = tok{1};
    else
        stem_no_tile = regexprep(stem, '[NnSs]\d+[EeWw]\d+', '');
        toks = regexp(stem_no_tile, '(?<![A-Za-z\d])(\d+)(?![A-Za-z\d])', 'tokens');
        if ~isempty(toks)
            meta.strip = toks{end}{1};
        end
    end
end

% -------------------------------------------------------------------------
function str = format_bytes(n)
%FORMAT_BYTES  Human-readable file size string.
    if n >= 1e6
        str = sprintf('%.1f MB', n/1e6);
    elseif n >= 1e3
        str = sprintf('%.1f KB', n/1e3);
    else
        str = sprintf('%d B', n);
    end
end

% -------------------------------------------------------------------------
function num_path = get_num_path(dem_path)
%GET_NUM_PATH  Return path of the companion _num.tif, or '' if absent.
    [d, stem, ext] = fileparts(dem_path);
    num_stem = regexprep(stem, '_[Dd][Ee][Mm]$', '_num');
    if strcmp(num_stem, stem)
        num_stem = [stem '_num'];
    end
    candidate = fullfile(d, [num_stem ext]);
    if exist(candidate, 'file')
        num_path = candidate;
    else
        num_path = '';
    end
end

% -------------------------------------------------------------------------
function stack = stack_edge_cells(edge_cells, mask, numeric_type)
%STACK_EDGE_CELLS Convert a cell stack to [K,H,W] and zero inactive edges.
    stack = cast(permute(cat(3, edge_cells{:}), [3, 1, 2]), numeric_type);
    stack(~logical(mask), :, :) = 0;
end

% -------------------------------------------------------------------------
function multilook = phase_only_multilook(phases, window_size)
%PHASE_ONLY_MULTILOOK Compatibility approximation for legacy phase samples.
    kernel = ones(window_size)/(window_size^2);
    multilook = cellfun(@(phase) angle(conv2(exp(1i*phase), ...
        kernel, 'same')), phases, 'UniformOutput', false);
end

% Run analyze_wrap_count_distribution.m separately after pointing its
% dataset_root at output_root. Automatic execution here could analyze the
% legacy training_dataset_v13 instead of the configured UAV dataset.

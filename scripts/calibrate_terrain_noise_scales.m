%% calibrate_terrain_noise_scales.m
% Estimate dataset-level absolute terrain scales from configured DEM files.
clc;

script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
addpath(genpath(project_root));
cfg_path = fullfile(project_root, 'configs', 'uav_p_500m_monostatic.json');
cfg = load_sim_config(cfg_path);

dem_dir = resolve_path(project_root, char(cfg.dataset.dem_directory));
files = [dir(fullfile(dem_dir, '*_dem*.tif')); ...
         dir(fullfile(dem_dir, '*_DSM*.tif'))];
assert(~isempty(files), 'dem2phase:NoCalibrationDem', ...
    'No configured DEM/DSM GeoTIFF files were found.');

rng(20260910);
patch_size = 256;
num_patches_per_file = 96;
scale_range = [2, 4];
sample_stride = 8;
slope_samples = [];
roughness_samples = [];
curvature_samples = [];
patch_records = struct([]);
record_idx = 0;

for file_idx = 1:numel(files)
    dem = double(imread(fullfile(files(file_idx).folder, files(file_idx).name)));
    [rows, cols] = size(dem);
    for patch_idx = 1:num_patches_per_file
        scale = scale_range(1) + diff(scale_range)*rand();
        source_size = max(16, ceil(patch_size/scale));
        if rows < source_size || cols < source_size; continue; end
        r0 = randi([1, rows-source_size+1]);
        c0 = randi([1, cols-source_size+1]);
        source = dem(r0:r0+source_size-1, c0:c0+source_size-1);
        if any(~isfinite(source(:))); continue; end
        patch = resize_bicubic(source, patch_size);
        spacing_m = double(cfg.dataset.source_dem_pixel_spacing_m) * ...
            source_size/patch_size;
        [~, terrain] = gen_coherence_map_multibaseline_terrainaware( ...
            patch, cfg.interferometry.baseline_perp_m(1), ...
            'pixel_spacing_m', spacing_m, ...
            'incidence_angle_deg', double(cfg.geometry.incidence_angle_deg), ...
            'look_azimuth_deg', double(cfg.geometry.look_azimuth_deg), ...
            'shared_residual_std', 0, 'edge_residual_std', 0, ...
            'residual_scale_px', 0.5, 'seed', 1);

        sample_idx = 1:sample_stride:numel(patch);
        slope_samples = [slope_samples; terrain.slope_deg(sample_idx)']; %#ok<AGROW>
        roughness_samples = [roughness_samples; terrain.roughness_m(sample_idx)']; %#ok<AGROW>
        curvature_samples = [curvature_samples; abs(terrain.curvature(sample_idx))']; %#ok<AGROW>
        record_idx = record_idx + 1;
        patch_records(record_idx).source_file = files(file_idx).name;
        patch_records(record_idx).scale = scale;
        patch_records(record_idx).pixel_spacing_m = spacing_m;
        patch_records(record_idx).mean_slope_deg = mean(terrain.slope_deg(:));
        patch_records(record_idx).mean_roughness_m = mean(terrain.roughness_m(:));
        patch_records(record_idx).mean_abs_curvature_per_m = ...
            mean(abs(terrain.curvature(:)));
    end
end

assert(~isempty(slope_samples), 'dem2phase:NoCalibrationSamples', ...
    'No valid terrain samples were collected.');
percentiles = [50, 75, 90, 95, 98, 99];
stats.percentiles = percentiles;
stats.slope_deg = empirical_percentiles(slope_samples, percentiles);
stats.roughness_m = empirical_percentiles(roughness_samples, percentiles);
stats.abs_curvature_per_m = empirical_percentiles(curvature_samples, percentiles);
stats.recommended.slope_scale_deg = stats.slope_deg(4);
stats.recommended.roughness_scale_m = stats.roughness_m(4);
stats.recommended.curvature_scale_per_m = stats.abs_curvature_per_m(4);
stats.num_patches = record_idx;
stats.num_pixel_samples = numel(slope_samples);
stats.seed = 20260910;

output_dir = fullfile(project_root, 'outputs');
if ~exist(output_dir, 'dir'); mkdir(output_dir); end
mat_path = fullfile(output_dir, 'terrain_noise_scale_calibration.mat');
json_path = fullfile(output_dir, 'terrain_noise_scale_calibration.json');
save(mat_path, 'stats', 'patch_records', '-v7');
fid = fopen(json_path, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', jsonencode(stats, PrettyPrint=true));
clear cleanup

fprintf('Collected %d patches and %d sampled pixels.\n', ...
    stats.num_patches, stats.num_pixel_samples);
fprintf('P95 slope       : %.6f deg\n', stats.recommended.slope_scale_deg);
fprintf('P95 roughness   : %.6f m\n', stats.recommended.roughness_scale_m);
fprintf('P95 |curvature| : %.9f 1/m\n', ...
    stats.recommended.curvature_scale_per_m);
fprintf('JSON=%s\nMAT=%s\n', json_path, mat_path);

function values = empirical_percentiles(data, percentages)
    sorted = sort(data(isfinite(data)));
    indices = max(1, min(numel(sorted), ceil(percentages/100*numel(sorted))));
    values = sorted(indices);
end

function resized = resize_bicubic(source, target_size)
    [rows, cols] = size(source);
    [x0, y0] = meshgrid(1:cols, 1:rows);
    [xq, yq] = meshgrid(linspace(1, cols, target_size), ...
        linspace(1, rows, target_size));
    resized = interp2(x0, y0, source, xq, yq, 'bicubic');
end

function path_out = resolve_path(project_root, configured_path)
    is_absolute = ~isempty(regexp(configured_path, '^[A-Za-z]:[\\/]', 'once')) ...
        || startsWith(configured_path, '\\') || startsWith(configured_path, '/');
    if is_absolute
        path_out = configured_path;
    else
        path_out = fullfile(project_root, configured_path);
    end
end

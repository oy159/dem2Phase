%% generate_preview_interferograms.m
% Generate one deterministic, terrain-rich preview using the nominal UAV profile.
clc; close all;

script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
addpath(genpath(project_root));
cfg = load_sim_config(fullfile(project_root, 'configs', ...
    'uav_p_500m_monostatic.json'));

dem_dir = resolve_preview_path(project_root, char(cfg.dataset.dem_directory));
files = dir(fullfile(dem_dir, '*_DSM*.tif'));
if isempty(files)
    files = dir(fullfile(dem_dir, '*_dem*.tif'));
end
assert(~isempty(files), 'dem2phase:NoPreviewDem', ...
    'No matching DEM/DSM GeoTIFF was found in %s.', dem_dir);

dem = double(imread(fullfile(files(1).folder, files(1).name)));
source_window = 96;
preview_size = 256;
[source_patch, crop_bounds] = select_representative_patch(dem, source_window, 0.65);
patch_dem = resize_bicubic(source_patch, preview_size);
patch_dem = patch_dem - min(patch_dem(:));
pixel_spacing_m = double(cfg.dataset.source_dem_pixel_spacing_m) * ...
    source_window / preview_size;

baselines = double(cfg.interferometry.baseline_perp_m(:)');
ratios = double(cfg.physics.dem2phase_ratios_rad_per_m(:)');
[coherence_terrain, terrain] = gen_coherence_map_multibaseline_terrainaware( ...
    patch_dem, baselines, ...
    'pixel_spacing_m', pixel_spacing_m, ...
    'incidence_angle_deg', double(cfg.geometry.incidence_angle_deg), ...
    'look_azimuth_deg', double(cfg.geometry.look_azimuth_deg), ...
    'min_coh', double(cfg.coherence.min), ...
    'max_coh', double(cfg.coherence.max), ...
    'baseline_ref', min(baselines), ...
    'baseline_decay_model', char(cfg.coherence.baseline_decay_model), ...
    'critical_baseline_m', double(cfg.physics.critical_baseline_m), ...
    'decay_alpha', double(cfg.coherence.decay_alpha), ...
    'terrain_normalization_model', ...
        char(cfg.coherence.terrain_normalization_model), ...
    'slope_scale_deg', double(cfg.coherence.slope_scale_deg), ...
    'roughness_scale_m', double(cfg.coherence.roughness_scale_m), ...
    'curvature_scale_per_m', double(cfg.coherence.curvature_scale_per_m), ...
    'slope_weight', double(cfg.coherence.slope_weight), ...
    'roughness_weight', double(cfg.coherence.roughness_weight), ...
    'curvature_weight', double(cfg.coherence.curvature_weight), ...
    'shared_residual_std', double(cfg.coherence.shared_residual_std), ...
    'edge_residual_std', double(cfg.coherence.edge_residual_std), ...
    'residual_scale_px', double(cfg.coherence.residual_scale_px), ...
    'seed', 20260910);

landcover_path = fullfile(project_root, char(cfg.dataset.landcover.aligned_file));
landcover_data = load(landcover_path, 'landcover_codes');
source_landcover = landcover_data.landcover_codes( ...
    crop_bounds(1):crop_bounds(2), crop_bounds(3):crop_bounds(4));
patch_landcover = resize_nearest(source_landcover, preview_size);
[coherence_scene, landcover_info] = apply_landcover_coherence( ...
    coherence_terrain, patch_landcover, cfg.dataset.landcover, ...
    double(cfg.coherence.min));

clean_unwrapped = arrayfun(@(k) patch_dem*ratios(k), ...
    1:numel(baselines), 'UniformOutput', false);
clean_wrapped = cellfun(@(x) angle(exp(1i*x)), clean_unwrapped, ...
    'UniformOutput', false);
node_snr_db = ones(1, double(cfg.interferometry.num_uavs))*22;
[noisy_wrapped, node_slc, coherence_true, noise_info] = ...
    simulate_distributed_slc_star(clean_wrapped, coherence_scene, ...
        cfg.interferometry.edge_index, double(cfg.interferometry.num_uavs), ...
        'node_snr_db', node_snr_db, 'seed', 20260911);
phase_error = cellfun(@(noisy, clean) angle(exp(1i*(noisy-clean))), ...
    noisy_wrapped, clean_wrapped, 'UniformOutput', false);
multilook_wrapped = multilook_node_interferograms(node_slc, ...
    cfg.interferometry.edge_index, double(cfg.noise.multilook_window));
multilook_error = cell(1, numel(baselines));
for k = 1:numel(baselines)
    multilook_error{k} = angle(exp(1i*(multilook_wrapped{k}-clean_wrapped{k})));
end

output_dir = fullfile(project_root, 'outputs');
if ~exist(output_dir, 'dir'); mkdir(output_dir); end
png_path = fullfile(output_dir, 'interferogram_preview_p500_500m.png');
mat_path = fullfile(output_dir, 'interferogram_preview_p500_500m.mat');

fig = figure('Color', 'w', 'Position', [50, 50, 2100, 1500]);
layout = tiledlayout(fig, 4, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
title(layout, sprintf(['Distributed UAV InSAR preview | 500 MHz, H=500 m, ' ...
    'incidence=45 deg, node SNR=22 dB']), 'FontWeight', 'bold');

nexttile; imagesc(patch_dem); axis image off; colorbar;
title(sprintf('DEM relief (m)\nsource crop [%d:%d, %d:%d]', crop_bounds));
colormap(gca, turbo);
for k = 1:numel(baselines)
    nexttile; imagesc(clean_wrapped{k}, [-pi, pi]); axis image off;
    title(sprintf('Clean phase | B_{perp}=%.2f m\nh_{amb}=%.1f m', ...
        baselines(k), cfg.physics.ambiguity_heights_m(k)));
    colormap(gca, hsv); colorbar;
end

nexttile; imagesc(worldcover_display_indices(patch_landcover), [0.5, 11.5]);
axis image off; title('ESA WorldCover 2021 classes');
colormap(gca, worldcover_colormap());
cb = colorbar; cb.Ticks = 1:11; cb.TickLabels = ...
    {'10','20','30','40','50','60','70','80','90','95','100'};
for k = 1:numel(baselines)
    nexttile; imagesc(noisy_wrapped{k}, [-pi, pi]); axis image off;
    title(sprintf('Noisy interferogram | edge %d-%d', ...
        cfg.interferometry.edge_index(1,k), cfg.interferometry.edge_index(2,k)));
    colormap(gca, hsv); colorbar;
end

nexttile; imagesc(terrain.slope_deg, [0, min(60, max(terrain.slope_deg(:)))]);
axis image off; title('Terrain slope (deg)'); colormap(gca, turbo); colorbar;
for k = 1:numel(baselines)
    nexttile; imagesc(coherence_true{k}, [0, 1]); axis image off;
    title(sprintf('True coherence | mean %.2f', mean(coherence_true{k}(:))));
    colormap(gca, parula); colorbar;
end

nexttile; imagesc(terrain.invalid_geometry_mask, [0, 1]); axis image off;
title('Layover / shadow proxy'); colormap(gca, gray); colorbar;
for k = 1:numel(baselines)
    nexttile; imagesc(multilook_wrapped{k}, [-pi, pi]); axis image off;
    rms_circular = sqrt(mean(multilook_error{k}(:).^2));
    title(sprintf('%dx%d complex multilook | error RMS %.2f rad', ...
        cfg.noise.multilook_window, cfg.noise.multilook_window, rms_circular));
    colormap(gca, hsv); colorbar;
end

exportgraphics(fig, png_path, 'Resolution', 150);
save(mat_path, 'patch_dem', 'clean_unwrapped', 'clean_wrapped', ...
    'noisy_wrapped', 'phase_error', 'multilook_wrapped', 'multilook_error', ...
    'coherence_terrain', 'coherence_scene', 'coherence_true', ...
    'terrain', 'patch_landcover', 'landcover_info', 'node_slc', ...
    'noise_info', 'baselines', 'crop_bounds', ...
    'pixel_spacing_m', 'cfg', '-v7');
fprintf('PNG=%s\nMAT=%s\n', png_path, mat_path);

function path_out = resolve_preview_path(project_root, configured_path)
    is_absolute = ~isempty(regexp(configured_path, '^[A-Za-z]:[\\/]', 'once')) ...
        || startsWith(configured_path, '\\') || startsWith(configured_path, '/');
    if is_absolute
        path_out = configured_path;
    else
        path_out = fullfile(project_root, configured_path);
    end
end

function [best_patch, bounds] = select_representative_patch(dem, window_size, quantile)
    [rows, cols] = size(dem);
    assert(rows >= window_size && cols >= window_size, ...
        'dem2phase:PreviewDemTooSmall', 'DEM is smaller than the preview window.');
    row_starts = unique(round(linspace(1, rows-window_size+1, 13)));
    col_starts = unique(round(linspace(1, cols-window_size+1, 13)));
    candidates = {};
    candidate_bounds = zeros(0, 4);
    scores = zeros(0, 1);
    for r0 = row_starts
        for c0 = col_starts
            candidate = dem(r0:r0+window_size-1, c0:c0+window_size-1);
            if any(~isfinite(candidate(:))); continue; end
            relief = max(candidate(:)) - min(candidate(:));
            roughness = std(candidate(:));
            score = relief + 2*roughness;
            candidates{end+1,1} = candidate; %#ok<AGROW>
            candidate_bounds(end+1,:) = ...
                [r0, r0+window_size-1, c0, c0+window_size-1]; %#ok<AGROW>
            scores(end+1,1) = score; %#ok<AGROW>
        end
    end
    assert(~isempty(candidates), 'dem2phase:NoValidPreviewPatch', ...
        'No finite preview patch could be selected.');
    [~, order] = sort(scores);
    selected_rank = max(1, min(numel(order), round(quantile*numel(order))));
    selected_idx = order(selected_rank);
    best_patch = candidates{selected_idx};
    bounds = candidate_bounds(selected_idx,:);
end

function resized = resize_bicubic(source, target_size)
    [rows, cols] = size(source);
    [x0, y0] = meshgrid(1:cols, 1:rows);
    [xq, yq] = meshgrid(linspace(1, cols, target_size), ...
        linspace(1, rows, target_size));
    resized = interp2(x0, y0, source, xq, yq, 'bicubic');
end

function resized = resize_nearest(source, target_size)
    row_idx = round(linspace(1, size(source,1), target_size));
    col_idx = round(linspace(1, size(source,2), target_size));
    resized = source(row_idx, col_idx);
end

function indices = worldcover_display_indices(codes)
    class_codes = [10,20,30,40,50,60,70,80,90,95,100];
    indices = zeros(size(codes));
    for idx = 1:numel(class_codes)
        indices(codes == class_codes(idx)) = idx;
    end
end

function colors = worldcover_colormap()
    colors = [0.0 0.39 0.0; 1.0 0.73 0.0; 1.0 1.0 0.30; ...
        0.94 0.59 0.0; 0.77 0.0 0.0; 0.71 0.71 0.71; ...
        1.0 1.0 1.0; 0.0 0.39 0.78; 0.0 0.59 0.59; ...
        0.0 0.39 0.39; 0.98 0.78 0.98];
end

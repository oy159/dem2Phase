%% generate_phase4_preview.m
% Compare node-noise-only and Phase-4-corrupted interferograms.
clc; close all;
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
dataset_root = fullfile(project_root, 'data', 'pilot_dataset_uav_p_500m_phase4');
files = dir(fullfile(dataset_root, 'patch_groups', '*.mat'));
assert(~isempty(files), 'dem2phase:NoPhase4Pilot', ...
    'Generate the Phase-4 pilot before creating its preview.');

data = load(fullfile(files(1).folder, files(1).name), ...
    'wrappedphase_node_noise_only', 'wrappedphase_withnoise', ...
    'coregistration_valid_mask', 'baseline_perp_m', 'ambiguity_height_m');
num_edges = size(data.wrappedphase_withnoise, 1);
figure('Color', 'w', 'Position', [50, 50, 1350, 300*num_edges]);
layout = tiledlayout(num_edges, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for edge_idx = 1:num_edges
    before = squeeze(data.wrappedphase_node_noise_only(edge_idx,:,:));
    after = squeeze(data.wrappedphase_withnoise(edge_idx,:,:));
    difference = angle(exp(1i*(double(after)-double(before))));
    valid = logical(squeeze(data.coregistration_valid_mask(edge_idx,:,:)));
    difference(~valid) = NaN;

    nexttile; imagesc(before, [-pi, pi]); axis image off; colormap(gca, hsv);
    title(sprintf('B_{\\perp}=%.2f m, h_{amb}=%.1f m | node noise', ...
        data.baseline_perp_m(edge_idx), data.ambiguity_height_m(edge_idx)));
    nexttile; imagesc(after, [-pi, pi]); axis image off; colormap(gca, hsv);
    title('with sync + trajectory + coreg');
    nexttile; imagesc(difference, [-pi, pi]); axis image off; colormap(gca, turbo);
    colorbar; title('wrapped Phase-4 difference');
end
title(layout, 'Distributed UAV InSAR Phase-4 error preview');
output_path = fullfile(project_root, 'outputs', 'phase4_error_preview.png');
exportgraphics(gcf, output_path, 'Resolution', 160);
fprintf('Saved %s\n', output_path);

%% generate_compound_failure_preview.m
% Visualize one positive four-factor compound-failure patch.
clc; close all;
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
root = fullfile(project_root, 'data', 'phase4_sweep', ...
    'compound_failure_failure', 'patch_groups');
files = dir(fullfile(root, '*.mat'));
selected = [];
for idx = 1:numel(files)
    candidate = load(fullfile(files(idx).folder, files(idx).name), ...
        'failure_labels');
    if candidate.failure_labels.compound_failure_present
        selected = idx;
        break
    end
end
assert(~isempty(selected), 'dem2phase:NoPositiveCompoundFailure', ...
    'No positive compound-failure patch is available.');
data = load(fullfile(files(selected).folder, files(selected).name));

figure('Color', 'w', 'Position', [40, 40, 1500, 950]);
tiledlayout(2, 4, 'TileSpacing', 'loose', 'Padding', 'loose');
for edge_idx = 1:4
    nexttile;
    phase = squeeze(data.wrappedphase_withnoise(edge_idx,:,:));
    if data.valid_edge_mask(edge_idx)
        state = 'ACTIVE';
        imagesc(phase, [-pi, pi]); axis image off; colormap(gca, hsv);
    else
        state = 'UAV DROPOUT';
        axis image off;
        text(0.5, 0.5, 'NO OBSERVATION', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold', ...
            'Color', [0.65 0.1 0.1]);
    end
    title(sprintf('B_{\\perp}=%.2f m — %s', ...
        data.baseline_perp_m(edge_idx), state));
end
nexttile;
coherence = squeeze(data.coherence_true(end,:,:));
imagesc(coherence, [0,1]); axis image off; colormap(gca, turbo); colorbar;
title(sprintf('Longest-edge coherence | low=%.1f%%', ...
    100*data.failure_labels.low_coherence_fraction));
nexttile;
imagesc(coherence < data.failure_labels.low_coherence_threshold);
axis image off; colormap(gca, gray); title('Low-coherence mask');
nexttile([1 2]);
plot(data.phase4_errors.sync_phase_error_rad', 'LineWidth', 1.2);
xlabel('Azimuth row'); ylabel('Synchronization phase error (rad)'); grid on;
legend(compose('UAV %d', 1:size(data.phase4_errors.sync_phase_error_rad,1)), ...
    'Location', 'best'); title('Node synchronization errors and jumps');
output_path = fullfile(project_root, 'outputs', 'compound_failure_preview.png');
exportgraphics(gcf, output_path, 'Resolution', 160);
fprintf('Saved %s\n', output_path);

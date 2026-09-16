%% analyze_wrap_count_distribution.m
% Analyze wrap-count distribution in training_dataset_v13 by baseline idx.
clear; clc; close all;

%% 1) 配置
script_dir    = fileparts(mfilename('fullpath'));   % dem2phase/scripts
dem2phase_dir = fileparts(script_dir);              % dem2phase
repo_root     = fileparts(dem2phase_dir);           % repo root

dataset_root = fullfile(repo_root, 'data', 'training_dataset_v13');
dataset_dir  = fullfile(dataset_root, 'unwrapped_phase');

% baseline 过滤：
%   []      -> 全部 baseline
%   [1 3]   -> 只分析 baseline1 和 baseline3
baseline_idx_filter = [];

% 直方图上限（若真实值更大会自动扩展到真实最大值）
max_possible_wrap = 50;

files = dir(fullfile(dataset_dir, '*.mat'));
if isempty(files)
    error('No .mat files found under: %s', dataset_dir);
end

fprintf('Found %d unwrapped phase files.\n', numel(files));

%% 2) 从文件名提取 baseline idx
num_files = numel(files);
file_baseline_idx = nan(num_files, 1);

for i = 1:num_files
    tok = regexp(files(i).name, '_baseline(\d+)_unwrapped_phase\.mat$', 'tokens', 'once');
    if ~isempty(tok)
        file_baseline_idx(i) = str2double(tok{1});
    end
end

if all(isnan(file_baseline_idx))
    error(['Failed to parse baseline idx from filenames. Expected format: ' ...
           '<patch_name>_baseline<k>_unwrapped_phase.mat']);
end

all_baselines = unique(file_baseline_idx(~isnan(file_baseline_idx)))';
if isempty(baseline_idx_filter)
    selected_baselines = all_baselines;
else
    selected_baselines = intersect(all_baselines, baseline_idx_filter, 'stable');
end

if isempty(selected_baselines)
    error('No overlap between baseline_idx_filter and dataset baselines. Available: %s', mat2str(all_baselines));
end

fprintf('Available baseline idx: %s\n', mat2str(all_baselines));
fprintf('Selected baseline idx: %s\n', mat2str(selected_baselines));

%% 3) 初始化统计容器
edges = 0:1:max_possible_wrap;
n_sel = numel(selected_baselines);

max_wrap_list_all = [];
pixel_wrap_hist_all = zeros(1, numel(edges)-1);

max_wrap_by_baseline = cell(1, n_sel);
mean_wrap_by_baseline = cell(1, n_sel);
pixel_hist_by_baseline = zeros(n_sel, numel(edges)-1);

detail_file_names = {};
detail_baseline_idx = [];
detail_max_wrap = [];
detail_mean_wrap = [];

%% 4) 遍历文件
tic;
processed = 0;
for i = 1:num_files
    bidx = file_baseline_idx(i);
    [is_sel, bpos] = ismember(bidx, selected_baselines);
    if ~is_sel
        continue
    end

    filepath = fullfile(dataset_dir, files(i).name);
    loaded = load(filepath);
    unw_phase = loaded.data;

    wrap_map = round((unw_phase - angle(exp(1j * unw_phase))) / (2*pi));
    max_wrap = max(wrap_map(:));
    mean_wrap = mean(wrap_map(:));

    max_wrap_list_all(end+1,1) = max_wrap; %#ok<SAGROW>
    max_wrap_by_baseline{bpos}(end+1,1) = max_wrap; %#ok<SAGROW>
    mean_wrap_by_baseline{bpos}(end+1,1) = mean_wrap; %#ok<SAGROW>

    if max_wrap >= max_possible_wrap
        extra_edges = (max_possible_wrap+1):max_wrap;
        if ~isempty(extra_edges)
            edges = [edges extra_edges]; %#ok<AGROW>
            pixel_wrap_hist_all(end+1:end+numel(extra_edges)) = 0; %#ok<AGROW>
            pixel_hist_by_baseline(:, end+1:end+numel(extra_edges)) = 0;
            max_possible_wrap = max_wrap;
        end
    end

    counts = histcounts(wrap_map(:), edges);
    pixel_wrap_hist_all = pixel_wrap_hist_all + counts;
    pixel_hist_by_baseline(bpos, :) = pixel_hist_by_baseline(bpos, :) + counts;

    detail_file_names{end+1,1} = files(i).name; %#ok<SAGROW>
    detail_baseline_idx(end+1,1) = bidx; %#ok<SAGROW>
    detail_max_wrap(end+1,1) = max_wrap; %#ok<SAGROW>
    detail_mean_wrap(end+1,1) = mean_wrap; %#ok<SAGROW>

    processed = processed + 1;
    if mod(processed, 500) == 0
        fprintf('Processed %d selected files...\n', processed);
    end
end
elapsed_time = toc;

if isempty(max_wrap_list_all)
    error('筛选后没有可分析文件。');
end

%% 5) 输出统计结果
fprintf('\nAnalysis done. Elapsed time: %.2f s\n', elapsed_time);
fprintf('----------------------------------------\n');
fprintf('Analyzed files: %d\n', numel(max_wrap_list_all));
fprintf('Global max wrap count: %d cycles\n', max(max_wrap_list_all));
fprintf('Global mean(max-wrap per file): %.3f cycles\n', mean(max_wrap_list_all));
fprintf('----------------------------------------\n');

for bi = 1:n_sel
    bidx = selected_baselines(bi);
    v = max_wrap_by_baseline{bi};
    m = mean_wrap_by_baseline{bi};
    if isempty(v)
        fprintf('baseline %d: no files\n', bidx);
        continue
    end
    fprintf('baseline %d -> files=%d, max-wrap(max)=%.0f, max-wrap(mean)=%.3f, max-wrap(p95)=%.3f, mean-wrap(mean)=%.3f\n', ...
        bidx, numel(v), max(v), mean(v), prctile(v,95), mean(m));
end

%% 6) 输出逐文件结果（便于人工判定基线合理性）
detail_tbl = table(detail_file_names, detail_baseline_idx, detail_max_wrap, detail_mean_wrap, ...
    'VariableNames', {'file_name', 'baseline_idx', 'max_wrap_count', 'mean_wrap_count'});
detail_out = fullfile(dataset_root, sprintf('analyze_wrapcount_detail_%s.csv', datestr(now, 'yyyymmdd_HHMMSS')));
writetable(detail_tbl, detail_out);
fprintf('\nPer-file stats saved: %s\n', detail_out);

%% 7) 可视化
figure('Position', [80, 80, 1300, 850], 'Name', 'Wrap Count Analysis by Baseline');

subplot(2, 2, 1);
histogram(max_wrap_list_all, 'BinMethod', 'integers', 'FaceColor', '#0072BD', 'EdgeColor', 'w');
title('All Selected Files: Max Wrap Count Distribution');
xlabel('Max Wrap Count (cycles)');
ylabel('Number of Files');
grid on;

subplot(2, 2, 2);
bar(edges(1:end-1), pixel_wrap_hist_all, 'FaceColor', '#D95319', 'EdgeColor', 'none');
title('All Selected Files: Pixel-wise Wrap Count Distribution');
xlabel('Wrap Count (cycles)');
ylabel('Number of Pixels (Log Scale)');
set(gca, 'YScale', 'log');
grid on;

subplot(2, 2, 3);
group_data = [];
group_label = [];
for bi = 1:n_sel
    v = max_wrap_by_baseline{bi};
    group_data = [group_data; v]; %#ok<AGROW>
    group_label = [group_label; repmat(selected_baselines(bi), numel(v), 1)]; %#ok<AGROW>
end
boxplot(group_data, group_label);
title('Max Wrap Count per File by Baseline');
xlabel('Baseline idx');
ylabel('Max Wrap Count (cycles)');
grid on;

subplot(2, 2, 4);
mean_max_wrap = zeros(1, n_sel);
for bi = 1:n_sel
    v = max_wrap_by_baseline{bi};
    mean_max_wrap(bi) = mean(v);
end
bar(selected_baselines, mean_max_wrap, 0.6, 'FaceColor', '#77AC30');
title('Mean of Max Wrap Count by Baseline');
xlabel('Baseline idx');
ylabel('Mean Max Wrap Count (cycles)');
grid on;

function manifest = build_failure_balanced_manifest(dataset_roots, output_csv)
%BUILD_FAILURE_BALANCED_MANIFEST Build mask-aware inverse-frequency sampling weights.
    script_dir = fileparts(mfilename('fullpath'));
    project_root = fileparts(script_dir);
    if nargin < 1 || isempty(dataset_roots)
        dataset_roots = {fullfile(project_root, 'data', ...
            'pilot_dataset_uav_p_500m_phase4'), ...
            fullfile(project_root, 'data', 'phase4_sweep')};
    end
    if ischar(dataset_roots) || isstring(dataset_roots)
        dataset_roots = cellstr(dataset_roots);
    end
    if nargin < 2 || isempty(output_csv)
        output_csv = fullfile(project_root, 'data', ...
            'failure_balanced_manifest.csv');
    end

    files = struct([]);
    for root_idx = 1:numel(dataset_roots)
        found = dir(fullfile(dataset_roots{root_idx}, '**', '*.mat'));
        in_patch_groups = contains({found.folder}, ...
            [filesep 'patch_groups'], 'IgnoreCase', true);
        found = found(in_patch_groups);
        files = [files; found]; %#ok<AGROW>
    end
    assert(~isempty(files), 'dem2phase:NoTrainingPatches', ...
        'No patch-group MAT files were found.');
    num_samples = numel(files);
    path = strings(num_samples,1);
    scenario = strings(num_samples,1);
    terrain = strings(num_samples,1);
    split_group = strings(num_samples,1);
    split = strings(num_samples,1);
    region = strings(num_samples,1);
    geographic_tile = strings(num_samples,1);
    source_dataset = strings(num_samples,1);
    source_quality_role = strings(num_samples,1);
    num_available = false(num_samples,1);
    num_mean = nan(num_samples,1);
    dropout = false(num_samples,1);
    sync_anomaly = false(num_samples,1);
    low_coherence = false(num_samples,1);
    compound = false(num_samples,1);
    active_edges = zeros(num_samples,1);
    landcover_code = zeros(num_samples,1);
    failure_factor_count = zeros(num_samples,1);

    for idx = 1:num_samples
        path(idx) = fullfile(files(idx).folder, files(idx).name);
        data = load(path(idx), 'failure_labels', 'valid_edge_mask', 'metadata');
        assert(isfield(data, 'failure_labels'), ...
            'dem2phase:MissingFailureLabels', ...
            'Sample %s uses an obsolete schema.', path(idx));
        labels = data.failure_labels;
        dropout(idx) = logical(labels.dropout_present);
        sync_anomaly(idx) = logical(labels.sync_anomaly_present);
        low_coherence(idx) = logical(labels.low_coherence_present);
        compound(idx) = logical(labels.compound_failure_present);
        failure_factor_count(idx) = double(labels.failure_factor_count);
        active_edges(idx) = nnz(data.valid_edge_mask);
        terrain(idx) = string(data.metadata.terrain_class);
        landcover_code(idx) = double(data.metadata.dominant_landcover_code);
        if isfield(data.metadata, 'split_group')
            split_group(idx) = string(data.metadata.split_group);
        else
            split_group(idx) = string(data.metadata.source_file);
        end
        if isfield(data.metadata, 'dataset_split')
            split(idx) = string(data.metadata.dataset_split);
        else
            split(idx) = "unassigned";
        end
        if isfield(data.metadata, 'geographic_region')
            region(idx) = string(data.metadata.geographic_region);
        else
            region(idx) = "unknown";
        end
        if isfield(data.metadata, 'geographic_tile')
            geographic_tile(idx) = string(data.metadata.geographic_tile);
        else
            geographic_tile(idx) = split_group(idx);
        end
        if isfield(data.metadata, 'source_dataset')
            source_dataset(idx) = string(data.metadata.source_dataset);
            source_quality_role(idx) = string(data.metadata.source_quality_role);
            num_available(idx) = logical(data.metadata.num_available);
            num_mean(idx) = double(data.metadata.num_mean);
        else
            source_dataset(idx) = "unknown";
            source_quality_role(idx) = "unknown";
        end
        scenario(idx) = string(data.metadata.phase4_error_profile);
    end

    stratum = compose('D%d_S%d_L%d_K%d_T%s', dropout, sync_anomaly, ...
        low_coherence, active_edges, terrain);
    sample_weight = ones(num_samples,1);
    train_mask = split == "train";
    if any(train_mask)
        sample_weight(train_mask) = compute_balanced_sample_weights( ...
            stratum(train_mask), 10);
    end
    manifest = table(path, scenario, split_group, geographic_tile, split, ...
        region, source_dataset, source_quality_role, num_available, num_mean, terrain, ...
        landcover_code, active_edges, dropout, sync_anomaly, low_coherence, ...
        compound, failure_factor_count, stratum, sample_weight);
    [group_id, group_name] = findgroups(split_group);
    split_count = splitapply(@(value) numel(unique(value)), split, group_id);
    leaking = group_name(split_count > 1);
    assert(isempty(leaking), 'dem2phase:SplitGroupLeakage', ...
        'Source DEMs assigned to more than one split: %s', ...
        strjoin(leaking, ', '));
    writetable(manifest, output_csv);
    count_csv = replace(output_csv, '.csv', '_strata.csv');
    [groups, split_names, stratum_names] = findgroups(split, stratum);
    sample_count = splitapply(@numel, split, groups);
    counts = table(split_names, stratum_names, sample_count, ...
        'VariableNames', {'split', 'stratum', 'sample_count'});
    writetable(counts, count_csv);
    fprintf('Wrote %d samples across %d strata to %s\n', ...
        height(manifest), height(counts), output_csv);
end

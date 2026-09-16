function [mask, labels] = apply_compound_failure_mask(initial_mask, edge_index, ...
    baselines_m, coherence_maps, phase4_info, cfg, seed)
%APPLY_COMPOUND_FAILURE_MASK Apply explicit UAV/edge failures and labels.
    rng_cleanup = scoped_rng(seed); %#ok<NASGU>
    mask = logical(initial_mask(:)');
    num_edges = numel(mask);
    num_nodes = max(edge_index(:));
    labels.profile = char(cfg.profile);
    labels.enabled = logical(cfg.enabled);
    labels.initial_edge_mask = uint8(mask);
    labels.edge_dropout_mask = zeros(1, num_edges, 'uint8');
    labels.uav_available_mask = ones(1, num_nodes, 'uint8');
    labels.low_coherence_threshold = double(cfg.low_coherence_threshold);
    labels.low_coherence_min_fraction = double(cfg.low_coherence_min_fraction);
    labels.sync_jump_min_abs_rad = double(cfg.sync_jump_min_abs_rad);

    [~, shortest_idx] = min(baselines_m);
    [~, longest_idx] = max(baselines_m);
    if logical(cfg.enabled)
        if logical(cfg.always_include_shortest)
            mask(shortest_idx) = true;
        end
        if logical(cfg.always_include_longest)
            mask(longest_idx) = true;
        end
        protected = unique([shortest_idx, longest_idx]);
        droppable = setdiff(1:num_edges, protected);
        if logical(cfg.force_secondary_uav_dropout) && ~isempty(droppable)
            drop_edge = droppable(randi(numel(droppable)));
            dropped_uav = edge_index(2,drop_edge);
            incident = any(edge_index == dropped_uav, 1);
            mask(incident) = false;
            labels.edge_dropout_mask(incident) = 1;
            labels.uav_available_mask(dropped_uav) = 0;
        end
    end

    low_count = 0;
    observed_count = 0;
    for edge_idx = find(mask)
        values = double(coherence_maps{edge_idx});
        low_count = low_count + nnz(values < cfg.low_coherence_threshold);
        observed_count = observed_count + numel(values);
    end
    labels.low_coherence_fraction = low_count/max(observed_count, 1);
    labels.low_coherence_present = labels.low_coherence_fraction >= ...
        cfg.low_coherence_min_fraction;
    jump_values = arrayfun(@(node) abs(node.sync_jump_rad), ...
        phase4_info.node_parameters);
    labels.sync_jump_node_mask = uint8(jump_values >= cfg.sync_jump_min_abs_rad);
    labels.sync_anomaly_present = any(labels.sync_jump_node_mask);
    labels.dropout_present = any(labels.edge_dropout_mask);
    labels.shortest_baseline_active = mask(shortest_idx);
    labels.longest_baseline_active = mask(longest_idx);
    labels.compound_failure_present = labels.dropout_present && ...
        labels.sync_anomaly_present && labels.low_coherence_present && ...
        labels.longest_baseline_active;
    labels.failure_factor_count = double(labels.dropout_present) + ...
        double(labels.sync_anomaly_present) + ...
        double(labels.low_coherence_present) + ...
        double(labels.longest_baseline_active);
    labels.final_edge_mask = uint8(mask);
end

function cleanup = scoped_rng(seed)
    cleanup = [];
    if ~isempty(seed)
        previous_state = rng;
        cleanup = onCleanup(@() rng(previous_state));
        rng(seed, 'twister');
    end
end

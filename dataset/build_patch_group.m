function group = build_patch_group(clean_wrap, noisy_wrap, clean_unwrapped, ...
    coherence_estimated, coherence_true, valid_edge_mask, cfg, ...
    terrain_features, metadata)
%BUILD_PATCH_GROUP Pack one multi-edge observation into [K,H,W] arrays.
%
% Inactive edges are zero-filled and must always be interpreted together
% with valid_edge_mask. Geometry remains stored for all configured edges.

    stacks = {clean_wrap, noisy_wrap, clean_unwrapped, ...
        coherence_estimated, coherence_true};
    num_edges = numel(clean_wrap);
    assert(num_edges > 0 && all(cellfun(@numel, stacks) == num_edges), ...
        'dem2phase:PatchGroupEdgeCountMismatch', ...
        'Every multi-edge field must contain the same non-zero edge count.');
    mask = logical(valid_edge_mask(:)');
    assert(numel(mask) == num_edges && any(mask), ...
        'dem2phase:InvalidValidEdgeMask', ...
        'valid_edge_mask must match the edge count and retain at least one edge.');

    target_size = size(clean_wrap{1});
    for field_idx = 1:numel(stacks)
        for edge_idx = 1:num_edges
            assert(isequal(size(stacks{field_idx}{edge_idx}), target_size), ...
                'dem2phase:PatchGroupShapeMismatch', ...
                'All edge rasters must share the same spatial dimensions.');
        end
    end

    numeric_type = char(cfg.dataset.storage.numeric_type);
    group.wrappedphase_withoutnoise = stack_and_mask(clean_wrap, mask, numeric_type);
    group.wrappedphase_withnoise = stack_and_mask(noisy_wrap, mask, numeric_type);
    group.unwrapped_phase = stack_and_mask(clean_unwrapped, mask, numeric_type);
    group.coherence_estimated = stack_and_mask(coherence_estimated, mask, numeric_type);
    group.coherence_true = stack_and_mask(coherence_true, mask, numeric_type);
    group.valid_edge_mask = uint8(mask);
    group.edge_index = uint16(cfg.interferometry.edge_index);
    group.baseline_perp_m = cast(cfg.interferometry.baseline_perp_m(:)', numeric_type);
    group.ambiguity_height_m = cast(cfg.physics.ambiguity_heights_m(:)', numeric_type);
    group.dem2phase_ratio_rad_per_m = cast( ...
        cfg.physics.dem2phase_ratios_rad_per_m(:)', numeric_type);
    group.phase_path_multiplicity = uint8(cfg.interferometry.phase_path_multiplicity);
    group.terrain_features = cast_terrain_features(terrain_features, numeric_type);
    group.metadata = metadata;
    group.schema_version = uint16(1);
end

function result = stack_and_mask(values, mask, numeric_type)
    result = permute(cat(3, values{:}), [3, 1, 2]);
    result = cast(result, numeric_type);
    result(~mask, :, :) = 0;
end

function terrain_out = cast_terrain_features(terrain_in, numeric_type)
    terrain_out = terrain_in;
    names = fieldnames(terrain_in);
    for idx = 1:numel(names)
        value = terrain_in.(names{idx});
        if isnumeric(value) && ~isscalar(value)
            terrain_out.(names{idx}) = cast(value, numeric_type);
        end
    end
end

function mask = sample_valid_edge_mask(num_edges, shortest_edge_idx, options, seed)
%SAMPLE_VALID_EDGE_MASK Draw a reproducible variable-K observation mask.

    arguments
        num_edges (1, 1) double {mustBeInteger, mustBePositive}
        shortest_edge_idx (1, 1) double {mustBeInteger, mustBePositive}
        options struct
        seed = []
    end
    rng_cleanup = scoped_rng(seed); %#ok<NASGU>
    assert(shortest_edge_idx <= num_edges, ...
        'dem2phase:InvalidShortestEdge', ...
        'shortest_edge_idx must lie within the edge list.');

    if ~logical(options.enabled)
        mask = true(1, num_edges);
        return
    end

    min_edges = double(options.min_active_edges);
    max_edges = double(options.max_active_edges);
    active_count = randi([min_edges, max_edges]);
    mask = false(1, num_edges);

    if logical(options.always_include_shortest)
        mask(shortest_edge_idx) = true;
    end
    remaining_count = active_count - nnz(mask);
    candidates = find(~mask);
    if remaining_count > 0
        selected = candidates(randperm(numel(candidates), remaining_count));
        mask(selected) = true;
    end
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

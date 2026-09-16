function graph = build_interferometric_graph(num_uavs, edge_mode, reference_uav)
%BUILD_INTERFEROMETRIC_GRAPH Construct the ordered interferometric edge list.

    arguments
        num_uavs (1, 1) double {mustBeInteger, mustBeGreaterThanOrEqual(num_uavs, 2)}
        edge_mode (1, :) char
        reference_uav (1, 1) double {mustBeInteger, mustBePositive}
    end

    assert(reference_uav <= num_uavs, ...
        'dem2phase:InvalidReferenceUav', ...
        'reference_uav must not exceed num_uavs.');

    switch edge_mode
        case 'star'
            secondary = setdiff(1:num_uavs, reference_uav, 'stable');
            edge_index = [repmat(reference_uav, 1, numel(secondary)); secondary];
        case 'chain'
            edge_index = [1:(num_uavs - 1); 2:num_uavs];
        case 'complete'
            pairs = nchoosek(1:num_uavs, 2);
            edge_index = pairs';
        otherwise
            error('dem2phase:UnsupportedEdgeMode', ...
                'Unsupported edge_mode=%s. Use star, chain, or complete.', edge_mode);
    end

    graph.edge_mode = edge_mode;
    graph.reference_uav = reference_uav;
    graph.edge_index = edge_index;
end


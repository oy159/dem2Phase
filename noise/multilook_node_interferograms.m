function multilook_phases = multilook_node_interferograms(node_slc, edge_index, window_size)
%MULTILOOK_NODE_INTERFEROGRAMS Complex boxcar multilook for every graph edge.

    assert(ndims(node_slc) == 3, 'dem2phase:InvalidNodeSlcShape', ...
        'node_slc must have shape [N,H,W].');
    assert(size(edge_index,1) == 2, 'dem2phase:InvalidEdgeIndexShape', ...
        'edge_index must have shape [2,K].');
    assert(window_size >= 1 && window_size == fix(window_size) && ...
        mod(window_size,2) == 1, 'dem2phase:InvalidMultilookWindow', ...
        'window_size must be a positive odd integer.');

    kernel = ones(window_size)/(window_size^2);
    num_edges = size(edge_index,2);
    multilook_phases = cell(1, num_edges);
    for edge_idx = 1:num_edges
        master = squeeze(node_slc(edge_index(1,edge_idx),:,:));
        secondary = squeeze(node_slc(edge_index(2,edge_idx),:,:));
        interferogram = master .* conj(secondary);
        multilook_phases{edge_idx} = angle(conv2(interferogram, kernel, 'same'));
    end
end

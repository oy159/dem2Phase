function [noisy_phases, node_slc, effective_coherence, info] = ...
    simulate_distributed_slc_star(clean_phases, scene_coherence, edge_index, ...
    num_uavs, varargin)
%SIMULATE_DISTRIBUTED_SLC_STAR Shared-node complex SLC noise for a star graph.
%
% The reference UAV provides one common complex speckle and receiver-noise
% realization to every interferometric edge. Each secondary UAV combines a
% phase-shifted correlated component, an edge decorrelation component and its
% own receiver noise. Thus edges sharing the master are statistically coupled.

    p = inputParser;
    addRequired(p, 'clean_phases', @(x) iscell(x) && ~isempty(x));
    addRequired(p, 'scene_coherence', @(x) iscell(x) && ~isempty(x));
    addRequired(p, 'edge_index', @(x) isnumeric(x) && size(x,1) == 2);
    addRequired(p, 'num_uavs', @(x) isnumeric(x) && isscalar(x) && x >= 2);
    addParameter(p, 'node_snr_db', 25, @(x) isnumeric(x) && isvector(x));
    addParameter(p, 'seed', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    parse(p, clean_phases, scene_coherence, edge_index, num_uavs, varargin{:});
    opt = p.Results;
    rng_cleanup = scoped_rng(opt.seed); %#ok<NASGU>

    num_edges = numel(clean_phases);
    assert(numel(scene_coherence) == num_edges && size(edge_index,2) == num_edges, ...
        'dem2phase:SlcEdgeCountMismatch', ...
        'Phase, coherence and edge_index counts must match.');
    reference_uav = edge_index(1,1);
    assert(all(edge_index(1,:) == reference_uav) && ...
        numel(unique(edge_index(2,:))) == num_edges, ...
        'dem2phase:SlcGraphNotStar', ...
        'edge_index must describe one directed star with unique secondary UAVs.');

    snr_db = double(opt.node_snr_db(:)');
    if isscalar(snr_db)
        snr_db = repmat(snr_db, 1, num_uavs);
    end
    assert(numel(snr_db) == num_uavs && all(isfinite(snr_db)), ...
        'dem2phase:NodeSnrCountMismatch', ...
        'node_snr_db must be scalar or contain one value per UAV.');
    snr_linear = 10.^(snr_db/10);

    image_size = size(clean_phases{1});
    common_speckle = complex_normal(image_size);
    node_slc_cells = cell(1, num_uavs);
    node_slc_cells{reference_uav} = common_speckle + ...
        complex_normal(image_size)/sqrt(snr_linear(reference_uav));
    effective_coherence = cell(1, num_edges);

    for edge_idx = 1:num_edges
        secondary_uav = edge_index(2,edge_idx);
        phase = double(clean_phases{edge_idx});
        gamma_scene = min(max(double(scene_coherence{edge_idx}), 0), 0.9999);
        assert(isequal(size(phase), image_size) && ...
            isequal(size(gamma_scene), image_size), ...
            'dem2phase:SlcShapeMismatch', ...
            'All phase and coherence rasters must share one size.');

        independent_scatter = complex_normal(image_size);
        secondary_clean = gamma_scene .* common_speckle .* exp(-1i*phase) + ...
            sqrt(max(1 - gamma_scene.^2, 0)) .* independent_scatter;
        node_slc_cells{secondary_uav} = secondary_clean + ...
            complex_normal(image_size)/sqrt(snr_linear(secondary_uav));

        thermal_factor = 1/sqrt((1 + 1/snr_linear(reference_uav)) * ...
            (1 + 1/snr_linear(secondary_uav)));
        effective_coherence{edge_idx} = gamma_scene * thermal_factor;
    end

    assert(all(~cellfun(@isempty, node_slc_cells)), ...
        'dem2phase:UnobservedUavNode', ...
        'Every UAV must participate in the configured star graph.');
    node_slc = permute(cat(3, node_slc_cells{:}), [3, 1, 2]);
    noisy_phases = cell(1, num_edges);
    for edge_idx = 1:num_edges
        master = edge_index(1,edge_idx);
        secondary = edge_index(2,edge_idx);
        noisy_phases{edge_idx} = angle(node_slc_cells{master} .* ...
            conj(node_slc_cells{secondary}));
    end

    info.node_snr_db = snr_db;
    info.reference_uav = reference_uav;
    info.model = 'complex_slc_nodes';
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

function values = complex_normal(array_size)
    values = (randn(array_size) + 1i*randn(array_size))/sqrt(2);
end

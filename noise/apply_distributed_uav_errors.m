function [node_slc_out, edge_phases, edge_valid, info] = ...
    apply_distributed_uav_errors(node_slc, edge_index, cfg, varargin)
%APPLY_DISTRIBUTED_UAV_ERRORS Propagate node-level Phase-4 errors to edges.
%
% Synchronization and LOS trajectory errors are phase fields. Attitude and
% coregistration errors create slowly varying subpixel warps. Every error is
% generated per UAV node, so edges sharing a UAV inherit correlated errors.

    p = inputParser;
    addRequired(p, 'node_slc', @(x) isnumeric(x) && ndims(x) == 3);
    addRequired(p, 'edge_index', @(x) isnumeric(x) && size(x,1) == 2);
    addRequired(p, 'cfg', @isstruct);
    addParameter(p, 'wavelength_m', 1, @(x) isscalar(x) && x > 0);
    addParameter(p, 'phase_path_multiplicity', 2, @(x) isscalar(x) && x > 0);
    addParameter(p, 'altitude_m', 500, @(x) isscalar(x) && x > 0);
    addParameter(p, 'ground_pixel_spacing_m', 10, @(x) isscalar(x) && x > 0);
    addParameter(p, 'seed', [], @(x) isempty(x) || (isscalar(x) && isnumeric(x)));
    parse(p, node_slc, edge_index, cfg, varargin{:});
    opt = p.Results;
    rng_cleanup = scoped_rng(opt.seed); %#ok<NASGU>

    [num_nodes, rows, cols] = size(node_slc);
    assert(all(edge_index(:) >= 1 & edge_index(:) <= num_nodes), ...
        'dem2phase:InvalidErrorEdgeIndex', 'edge_index references an absent UAV.');
    sync_cfg = cfg.synchronization;
    trajectory_cfg = cfg.trajectory;
    attitude_cfg = cfg.attitude;
    coreg_cfg = cfg.coregistration;

    row_axis = linspace(-0.5, 0.5, rows)';
    row_unit = linspace(0, 1, rows)';
    [x_grid, y_grid] = meshgrid(1:cols, 1:rows);
    sync_phase = zeros(num_nodes, rows, 'single');
    trajectory_phase = zeros(num_nodes, rows, 'single');
    los_error_m = zeros(num_nodes, rows, 'single');
    range_displacement_px = zeros(num_nodes, rows, 'single');
    azimuth_displacement_px = zeros(num_nodes, rows, 'single');
    node_valid = false(num_nodes, rows, cols);
    node_slc_out = zeros(size(node_slc), 'like', node_slc);
    node_parameters = repmat(struct(), 1, num_nodes);

    for node_idx = 1:num_nodes
        bias = sync_cfg.phase_bias_std_rad * randn;
        drift = sync_cfg.linear_drift_std_rad * randn;
        random_walk = cumsum(sync_cfg.random_walk_std_rad_per_row * randn(rows,1));
        random_walk = random_walk - mean(random_walk);
        jump = zeros(rows,1);
        jump_row = 0;
        jump_amplitude = 0;
        if rand < sync_cfg.jump_probability_per_node
            jump_row = randi([2, rows]);
            jump_amplitude = sync_cfg.jump_std_rad * randn;
            jump(jump_row:end) = jump_amplitude;
        end
        sync_vector = bias + drift*row_axis + random_walk + jump;

        los_bias = trajectory_cfg.los_bias_std_m * randn;
        los_drift = trajectory_cfg.los_drift_std_m * randn;
        vibration_amplitude = trajectory_cfg.vibration_std_m * randn;
        cycles = trajectory_cfg.vibration_cycles_range(1) + ...
            diff(trajectory_cfg.vibration_cycles_range)*rand;
        vibration_phase = 2*pi*rand;
        los_vector = los_bias + los_drift*row_axis + vibration_amplitude * ...
            sin(2*pi*cycles*row_unit + vibration_phase);
        trajectory_vector = opt.phase_path_multiplicity*2*pi/opt.wavelength_m * ...
            los_vector;

        roll_deg = attitude_cfg.roll_std_deg * randn;
        pitch_deg = attitude_cfg.pitch_std_deg * randn;
        range_shift = coreg_cfg.range_shift_std_px*randn + ...
            opt.altitude_m*tand(roll_deg)/opt.ground_pixel_spacing_m;
        azimuth_shift = coreg_cfg.azimuth_shift_std_px*randn + ...
            opt.altitude_m*tand(pitch_deg)/opt.ground_pixel_spacing_m;
        range_drift = coreg_cfg.linear_drift_std_px*randn;
        azimuth_drift = coreg_cfg.linear_drift_std_px*randn;
        range_map = range_shift + range_drift*row_axis;
        azimuth_map = azimuth_shift + azimuth_drift*row_axis;
        range_map = repmat(range_map, 1, cols);
        azimuth_map = repmat(azimuth_map, 1, cols);

        total_phase = repmat(sync_vector + trajectory_vector, 1, cols);
        phased = squeeze(node_slc(node_idx,:,:)) .* exp(1i*total_phase);
        query_x = x_grid - range_map;
        query_y = y_grid - azimuth_map;
        warped = interp2(x_grid, y_grid, real(phased), query_x, query_y, ...
            'linear', 0) + 1i*interp2(x_grid, y_grid, imag(phased), ...
            query_x, query_y, 'linear', 0);
        valid = query_x >= 1 & query_x <= cols & query_y >= 1 & query_y <= rows;

        node_slc_out(node_idx,:,:) = warped;
        sync_phase(node_idx,:) = single(sync_vector);
        trajectory_phase(node_idx,:) = single(trajectory_vector);
        los_error_m(node_idx,:) = single(los_vector);
        range_displacement_px(node_idx,:) = single(range_map(:,1));
        azimuth_displacement_px(node_idx,:) = single(azimuth_map(:,1));
        node_valid(node_idx,:,:) = valid;
        node_parameters(node_idx).sync_bias_rad = bias;
        node_parameters(node_idx).sync_drift_rad = drift;
        node_parameters(node_idx).sync_jump_row = jump_row;
        node_parameters(node_idx).sync_jump_rad = jump_amplitude;
        node_parameters(node_idx).los_bias_m = los_bias;
        node_parameters(node_idx).los_drift_m = los_drift;
        node_parameters(node_idx).vibration_amplitude_m = vibration_amplitude;
        node_parameters(node_idx).vibration_cycles = cycles;
        node_parameters(node_idx).roll_deg = roll_deg;
        node_parameters(node_idx).pitch_deg = pitch_deg;
        node_parameters(node_idx).range_shift_px = range_shift;
        node_parameters(node_idx).azimuth_shift_px = azimuth_shift;
    end

    num_edges = size(edge_index,2);
    edge_phases = cell(1, num_edges);
    edge_valid = false(num_edges, rows, cols);
    for edge_idx = 1:num_edges
        master = edge_index(1,edge_idx);
        secondary = edge_index(2,edge_idx);
        edge_phases{edge_idx} = angle(squeeze(node_slc_out(master,:,:)) .* ...
            conj(squeeze(node_slc_out(secondary,:,:))));
        edge_valid(edge_idx,:,:) = node_valid(master,:,:) & node_valid(secondary,:,:);
    end

    info.profile = char(cfg.profile);
    info.seed = opt.seed;
    info.sync_phase_error_rad = sync_phase;
    info.trajectory_phase_error_rad = trajectory_phase;
    info.los_range_error_m = los_error_m;
    info.range_displacement_px = range_displacement_px;
    info.azimuth_displacement_px = azimuth_displacement_px;
    info.node_parameters = node_parameters;
    info.profile_axis = 'azimuth_row';
end

function cleanup = scoped_rng(seed)
    cleanup = [];
    if ~isempty(seed)
        previous_state = rng;
        cleanup = onCleanup(@() rng(previous_state));
        rng(seed, 'twister');
    end
end

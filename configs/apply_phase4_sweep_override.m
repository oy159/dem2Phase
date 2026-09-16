function cfg = apply_phase4_sweep_override(cfg, scenario, level)
%APPLY_PHASE4_SWEEP_OVERRIDE Build paired Phase-4 ablation configurations.
    scenario = char(lower(string(scenario)));
    level = char(lower(string(level)));
    level_names = {'nominal', 'challenging', 'failure'};
    level_multipliers = [1, 3, 10];
    level_idx = find(strcmp(level_names, level), 1);
    assert(~isempty(level_idx), 'dem2phase:InvalidPhase4SweepLevel', ...
        'Phase-4 level must be nominal, challenging, or failure.');
    multiplier = level_multipliers(level_idx);

    base = cfg.phase4_errors;
    sync = zero_synchronization(base.synchronization);
    trajectory = zero_trajectory(base.trajectory);
    attitude = zero_attitude(base.attitude);
    coregistration = zero_coregistration(base.coregistration);
    switch scenario
        case 'sync_only'
            sync = scale_synchronization(base.synchronization, multiplier);
        case 'trajectory_only'
            trajectory = scale_numeric_fields(base.trajectory, multiplier, ...
                {'vibration_cycles_range'});
        case 'coreg_only'
            attitude = scale_numeric_fields(base.attitude, multiplier, {});
            coregistration = scale_numeric_fields(base.coregistration, multiplier, {});
        case 'combined'
            sync = scale_synchronization(base.synchronization, multiplier);
            trajectory = scale_numeric_fields(base.trajectory, multiplier, ...
                {'vibration_cycles_range'});
            attitude = scale_numeric_fields(base.attitude, multiplier, {});
            coregistration = scale_numeric_fields(base.coregistration, multiplier, {});
        case 'compound_failure'
            sync = scale_synchronization(base.synchronization, multiplier);
            sync.jump_probability_per_node = 1;
            trajectory = scale_numeric_fields(base.trajectory, multiplier, ...
                {'vibration_cycles_range'});
            attitude = scale_numeric_fields(base.attitude, multiplier, {});
            coregistration = scale_numeric_fields(base.coregistration, multiplier, {});
            cfg.phase4_failures.enabled = true;
        otherwise
            error('dem2phase:InvalidPhase4SweepScenario', ...
                ['Phase-4 scenario must be sync_only, trajectory_only, ' ...
                 'coreg_only, combined, or compound_failure.']);
    end

    cfg.phase4_errors.enabled = true;
    cfg.phase4_errors.profile = sprintf('pband_uav_%s_%s_v1_unvalidated', ...
        scenario, level);
    cfg.phase4_errors.synchronization = sync;
    cfg.phase4_errors.trajectory = trajectory;
    cfg.phase4_errors.attitude = attitude;
    cfg.phase4_errors.coregistration = coregistration;
    cfg.dataset.output_directory = fullfile('data', 'phase4_sweep', ...
        sprintf('%s_%s', scenario, level));
end

function value = scale_synchronization(value, multiplier)
    value = scale_numeric_fields(value, multiplier, ...
        {'jump_probability_per_node'});
    value.jump_probability_per_node = min(1, ...
        value.jump_probability_per_node*multiplier);
end

function value = scale_numeric_fields(value, multiplier, excluded)
    names = fieldnames(value);
    for idx = 1:numel(names)
        name = names{idx};
        if isnumeric(value.(name)) && ~ismember(name, excluded)
            value.(name) = value.(name)*multiplier;
        end
    end
end

function value = zero_synchronization(value)
    value.phase_bias_std_rad = 0;
    value.linear_drift_std_rad = 0;
    value.random_walk_std_rad_per_row = 0;
    value.jump_probability_per_node = 0;
    value.jump_std_rad = 0;
end

function value = zero_trajectory(value)
    value.los_bias_std_m = 0;
    value.los_drift_std_m = 0;
    value.vibration_std_m = 0;
end

function value = zero_attitude(value)
    value.roll_std_deg = 0;
    value.pitch_std_deg = 0;
end

function value = zero_coregistration(value)
    value.range_shift_std_px = 0;
    value.azimuth_shift_std_px = 0;
    value.linear_drift_std_px = 0;
end

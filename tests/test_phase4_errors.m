function tests = test_phase4_errors
%TEST_PHASE4_ERRORS Tests node-level synchronization, motion and coreg errors.
    tests = functiontests(localfunctions);
end

function testZeroErrorsPreserveNodeSlc(testCase)
    cfg = zero_error_config();
    node_slc = complex(randn(3,24,20), randn(3,24,20));
    edges = [1 1; 2 3];
    [output, ~, valid, info] = apply_distributed_uav_errors( ...
        node_slc, edges, cfg, 'seed', 11);

    verifyEqual(testCase, output, node_slc, 'AbsTol', 1e-12);
    verifyTrue(testCase, all(valid(:)));
    verifyEqual(testCase, info.sync_phase_error_rad, zeros(3,24,'single'));
end

function testSeedReplaysAndRestoresCallerStream(testCase)
    cfg = nominal_error_config();
    node_slc = ones(3,32,28);
    edges = [1 1; 2 3];
    rng(719);
    expected_first = rand;
    expected_second = rand;

    rng(719);
    actual_first = rand;
    [out_a, phase_a, valid_a] = apply_distributed_uav_errors( ...
        node_slc, edges, cfg, 'seed', 91);
    actual_second = rand;
    [out_b, phase_b, valid_b] = apply_distributed_uav_errors( ...
        node_slc, edges, cfg, 'seed', 91);

    verifyEqual(testCase, actual_first, expected_first);
    verifyEqual(testCase, actual_second, expected_second);
    verifyEqual(testCase, out_a, out_b);
    verifyEqual(testCase, phase_a, phase_b);
    verifyEqual(testCase, valid_a, valid_b);
end

function testNodePhasePropagatesAsEdgeDifference(testCase)
    cfg = zero_error_config();
    cfg.synchronization.phase_bias_std_rad = 0.2;
    node_slc = ones(3,16,12);
    edges = [1 1; 2 3];
    [~, phases, ~, info] = apply_distributed_uav_errors( ...
        node_slc, edges, cfg, 'seed', 33);

    expected_row = double(info.sync_phase_error_rad(1,:)) - ...
        double(info.sync_phase_error_rad(2,:));
    expected = repmat(expected_row(:), 1, size(node_slc,3));
    residual = angle(exp(1i*(phases{1}-expected)));
    verifyLessThan(testCase, max(abs(residual(:))), 1e-6);
end

function testSubpixelWarpProducesValidityMask(testCase)
    cfg = zero_error_config();
    cfg.coregistration.range_shift_std_px = 2;
    node_slc = ones(3,24,20);
    edges = [1 1; 2 3];
    [~, ~, valid] = apply_distributed_uav_errors( ...
        node_slc, edges, cfg, 'seed', 8);

    verifyLessThan(testCase, mean(valid(:)), 1);
    verifyGreaterThan(testCase, mean(valid(:)), 0.5);
end

function testSweepAblationsEnableOnlyRequestedComponents(testCase)
    here = fileparts(mfilename('fullpath'));
    full_cfg = load_sim_config(fullfile(here, '..', 'configs', ...
        'uav_p_500m_monostatic.json'));
    sync_cfg = apply_phase4_sweep_override(full_cfg, 'sync_only', 'challenging');
    coreg_cfg = apply_phase4_sweep_override(full_cfg, 'coreg_only', 'failure');

    verifyEqual(testCase, sync_cfg.phase4_errors.synchronization.phase_bias_std_rad, ...
        3*full_cfg.phase4_errors.synchronization.phase_bias_std_rad);
    verifyEqual(testCase, sync_cfg.phase4_errors.trajectory.los_bias_std_m, 0);
    verifyEqual(testCase, sync_cfg.phase4_errors.coregistration.range_shift_std_px, 0);
    verifyEqual(testCase, coreg_cfg.phase4_errors.synchronization.phase_bias_std_rad, 0);
    verifyEqual(testCase, coreg_cfg.phase4_errors.attitude.roll_std_deg, ...
        10*full_cfg.phase4_errors.attitude.roll_std_deg);
    verifyEqual(testCase, coreg_cfg.dataset.output_directory, ...
        fullfile('data', 'phase4_sweep', 'coreg_only_failure'));
end

function testSweepRejectsUnknownLevel(testCase)
    here = fileparts(mfilename('fullpath'));
    full_cfg = load_sim_config(fullfile(here, '..', 'configs', ...
        'uav_p_500m_monostatic.json'));
    verifyError(testCase, @() apply_phase4_sweep_override( ...
        full_cfg, 'combined', 'extreme'), 'dem2phase:InvalidPhase4SweepLevel');
end

function testCompoundFailurePreservesProtectedBaselines(testCase)
    full_cfg = full_config();
    cfg = full_cfg.phase4_failures;
    cfg.enabled = true;
    initial = logical([1 1 0 0]);
    edges = [1 1 1 1; 2 3 4 5];
    baselines = [0.75 1.2 1.8 3.0];
    coherence = repmat({ones(8)*0.1}, 1, 4);
    phase4_info.node_parameters = repmat(struct('sync_jump_rad', 0.2), 1, 5);
    [mask, labels] = apply_compound_failure_mask(initial, edges, baselines, ...
        coherence, phase4_info, cfg, 17);

    verifyTrue(testCase, mask(1));
    verifyTrue(testCase, mask(4));
    verifyEqual(testCase, nnz(labels.edge_dropout_mask), 1);
    verifyEqual(testCase, nnz(labels.uav_available_mask == 0), 1);
    verifyTrue(testCase, labels.compound_failure_present);
    verifyEqual(testCase, labels.failure_factor_count, 4);
end

function testHighCoherencePreventsCompoundLabel(testCase)
    full_cfg = full_config();
    cfg = full_cfg.phase4_failures;
    cfg.enabled = true;
    edges = [1 1 1 1; 2 3 4 5];
    coherence = repmat({ones(6)*0.9}, 1, 4);
    phase4_info.node_parameters = repmat(struct('sync_jump_rad', 0.2), 1, 5);
    [~, labels] = apply_compound_failure_mask(true(1,4), edges, ...
        [0.75 1.2 1.8 3.0], coherence, phase4_info, cfg, 5);

    verifyFalse(testCase, labels.low_coherence_present);
    verifyFalse(testCase, labels.compound_failure_present);
end

function cfg = nominal_error_config()
    full_cfg = full_config();
    cfg = full_cfg.phase4_errors;
end

function cfg = full_config()
    here = fileparts(mfilename('fullpath'));
    cfg = load_sim_config(fullfile(here, '..', 'configs', ...
        'uav_p_500m_monostatic.json'));
end

function cfg = zero_error_config()
    cfg = nominal_error_config();
    cfg.synchronization.phase_bias_std_rad = 0;
    cfg.synchronization.linear_drift_std_rad = 0;
    cfg.synchronization.random_walk_std_rad_per_row = 0;
    cfg.synchronization.jump_probability_per_node = 0;
    cfg.synchronization.jump_std_rad = 0;
    cfg.trajectory.los_bias_std_m = 0;
    cfg.trajectory.los_drift_std_m = 0;
    cfg.trajectory.vibration_std_m = 0;
    cfg.attitude.roll_std_deg = 0;
    cfg.attitude.pitch_std_deg = 0;
    cfg.coregistration.range_shift_std_px = 0;
    cfg.coregistration.azimuth_shift_std_px = 0;
    cfg.coregistration.linear_drift_std_px = 0;
end

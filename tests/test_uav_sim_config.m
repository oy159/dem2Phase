function tests = test_uav_sim_config
%TEST_UAV_SIM_CONFIG Unit tests for configurable multi-UAV InSAR physics.
    tests = functiontests(localfunctions);
end

function testNominalPBandProfile(testCase)
    here = fileparts(mfilename('fullpath'));
    config_path = fullfile(here, '..', 'configs', 'uav_p_500m_monostatic.json');
    cfg = load_sim_config(config_path);

    verifyEqual(testCase, cfg.physics.wavelength_m, 299792458/500e6, ...
        'RelTol', 1e-12);
    verifyEqual(testCase, cfg.physics.range_resolution_m, 299792458/(2*200e6), ...
        'RelTol', 1e-12);
    verifyGreaterThan(testCase, cfg.physics.critical_baseline_m, 280);
    verifyLessThan(testCase, cfg.physics.critical_baseline_m, 285);
    verifyEqual(testCase, char(cfg.coherence.terrain_normalization_model), ...
        'absolute_scales');
    verifyEqual(testCase, cfg.coherence.roughness_scale_m, 71.823859, ...
        'AbsTol', 1e-9);
    verifyEqual(testCase, cfg.interferometry.num_edges, 4);
    verifyEqual(testCase, cfg.interferometry.edge_index, [1 1 1 1; 2 3 4 5]);
    verifyEqual(testCase, cfg.physics.ambiguity_heights_m, ...
        [199.8616387, 124.9135242, 83.2756828, 49.9654097], ...
        'RelTol', 1e-8);
end

function testSingleTransmitterUsesSinglePathContribution(testCase)
    here = fileparts(mfilename('fullpath'));
    config_path = fullfile(here, '..', 'configs', 'uav_p_500m_monostatic.json');
    cfg = jsondecode(fileread(config_path));
    cfg.interferometry.measurement_mode = 'single_tx_multireceiver';
    cfg.interferometry.phase_path_multiplicity = 1;
    cfg.interferometry.target_ambiguity_height_m = ...
        2 * cfg.interferometry.target_ambiguity_height_m;
    cfg = validate_sim_config(cfg);

    verifyEqual(testCase, cfg.physics.ambiguity_heights_m, ...
        [399.7232773, 249.8270483, 166.5513655, 99.9308193], ...
        'RelTol', 1e-8);
end

function testRejectsWrongPathMultiplicity(testCase)
    here = fileparts(mfilename('fullpath'));
    config_path = fullfile(here, '..', 'configs', 'uav_p_500m_monostatic.json');
    cfg = jsondecode(fileread(config_path));
    cfg.interferometry.phase_path_multiplicity = 1;

    verifyError(testCase, @() validate_sim_config(cfg), ...
        'dem2phase:PathMultiplicityMismatch');
end

function tests = test_distributed_slc_noise
%TEST_DISTRIBUTED_SLC_NOISE Tests shared-node complex SLC simulation.
    tests = functiontests(localfunctions);
end

function testEffectiveCoherenceIncludesNodeSnr(testCase)
    phases = repmat({zeros(64)}, 1, 2);
    scene = repmat({ones(64)*0.8}, 1, 2);
    edges = [1 1; 2 3];
    snr_db = [10, 20, 30];
    [noisy, slc, effective, info] = simulate_distributed_slc_star( ...
        phases, scene, edges, 3, 'node_snr_db', snr_db, 'seed', 5);

    expected1 = 0.8/sqrt((1 + 1/10) * (1 + 1/100));
    expected2 = 0.8/sqrt((1 + 1/10) * (1 + 1/1000));
    verifyEqual(testCase, effective{1}, ones(64)*expected1, 'AbsTol', 1e-12);
    verifyEqual(testCase, effective{2}, ones(64)*expected2, 'AbsTol', 1e-12);
    verifySize(testCase, slc, [3, 64, 64]);
    verifyNumElements(testCase, noisy, 2);
    verifyEqual(testCase, info.node_snr_db, snr_db);
end

function testSharedMasterCreatesCrossEdgeErrorCorrelation(testCase)
    n = 192;
    phases = repmat({zeros(n)}, 1, 2);
    scene = repmat({ones(n)*0.9999}, 1, 2);
    edges = [1 1; 2 3];
    noisy = simulate_distributed_slc_star(phases, scene, edges, 3, ...
        'node_snr_db', [10, 80, 80], 'seed', 17);
    rho = corrcoef(sin(noisy{1}(:)), sin(noisy{2}(:)));

    verifyGreaterThan(testCase, rho(1,2), 0.5);
end

function testSeedIsReproducible(testCase)
    phases = {ones(32)*0.4, ones(32)*-0.7};
    scene = {ones(32)*0.75, ones(32)*0.6};
    edges = [1 1; 2 3];
    [phase_a, slc_a] = simulate_distributed_slc_star( ...
        phases, scene, edges, 3, 'node_snr_db', 25, 'seed', 29);
    [phase_b, slc_b] = simulate_distributed_slc_star( ...
        phases, scene, edges, 3, 'node_snr_db', 25, 'seed', 29);

    verifyEqual(testCase, phase_a, phase_b);
    verifyEqual(testCase, slc_a, slc_b);
end

function testLocalSeedDoesNotResetCallerStream(testCase)
    phases = {zeros(16), zeros(16)};
    scene = {ones(16)*0.8, ones(16)*0.7};
    edges = [1 1; 2 3];
    rng(321);
    expected_first = rand;
    expected_second = rand;

    rng(321);
    actual_first = rand;
    simulate_distributed_slc_star(phases, scene, edges, 3, ...
        'node_snr_db', 25, 'seed', 77);
    actual_second = rand;

    verifyEqual(testCase, actual_first, expected_first);
    verifyEqual(testCase, actual_second, expected_second);
end

function testComplexMultilookReducesPhaseError(testCase)
    n = 160;
    clean = repmat({ones(n)*0.7}, 1, 2);
    scene = repmat({ones(n)*0.65}, 1, 2);
    edges = [1 1; 2 3];
    [raw, slc] = simulate_distributed_slc_star(clean, scene, edges, 3, ...
        'node_snr_db', 22, 'seed', 43);
    looked = multilook_node_interferograms(slc, edges, 7);
    raw_error = angle(exp(1i*(raw{1}-clean{1})));
    looked_error = angle(exp(1i*(looked{1}-clean{1})));

    verifyLessThan(testCase, sqrt(mean(looked_error(:).^2)), ...
        sqrt(mean(raw_error(:).^2)));
end

function tests = test_terrainaware_coherence
%TEST_TERRAINAWARE_COHERENCE Tests for DEM-conditioned coherence maps.
    tests = functiontests(localfunctions);
end

function testFlatTerrainHasHigherCoherenceThanRoughTerrain(testCase)
    n = 96;
    flat = zeros(n);
    [x, y] = meshgrid(1:n, 1:n);
    rough = 0.15*x + 8*sin(x/3).*sin(y/4);
    baselines = [0.75, 1.2, 1.8, 3.0];

    flat_maps = gen_coherence_map_multibaseline_terrainaware(flat, baselines, ...
        'pixel_spacing_m', 1, 'seed', 7, ...
        'shared_residual_std', 0, 'edge_residual_std', 0);
    rough_maps = gen_coherence_map_multibaseline_terrainaware(rough, baselines, ...
        'pixel_spacing_m', 1, 'seed', 7, ...
        'shared_residual_std', 0, 'edge_residual_std', 0);

    verifyGreaterThan(testCase, mean(flat_maps{1}(:)), mean(rough_maps{1}(:)));
end

function testLongerBaselinesHaveLowerMeanCoherence(testCase)
    [x, y] = meshgrid(1:80, 1:80);
    dem = 0.05*x + 0.08*y;
    baselines = [0.75, 1.2, 1.8, 3.0];
    maps = gen_coherence_map_multibaseline_terrainaware(dem, baselines, ...
        'pixel_spacing_m', 2, 'seed', 11);
    means = cellfun(@(v) mean(v(:)), maps);

    verifyTrue(testCase, all(diff(means) < 0));
    verifyGreaterThanOrEqual(testCase, min(cellfun(@(v) min(v(:)), maps)), 0);
    verifyLessThanOrEqual(testCase, max(cellfun(@(v) max(v(:)), maps)), 1);
end

function testMeterScaleBaselinesHaveMildPhysicalDecorrelation(testCase)
    dem = zeros(64);
    baselines = [0.75, 1.2, 1.8, 3.0];
    maps = gen_coherence_map_multibaseline_terrainaware(dem, baselines, ...
        'pixel_spacing_m', 10, 'seed', 12, ...
        'baseline_decay_model', 'critical_baseline', ...
        'critical_baseline_m', 282.65, ...
        'shared_residual_std', 0, 'edge_residual_std', 0);
    means = cellfun(@(value) mean(value(:)), maps);

    verifyTrue(testCase, all(diff(means) < 0));
    verifyLessThan(testCase, means(1)-means(end), 0.02);
end

function testAbsoluteScaleIsIndependentOfDistantPatchContent(testCase)
    [x, ~] = meshgrid(1:128, 1:128);
    dem_a = x;
    dem_b = dem_a;
    dem_b(1:50,:) = dem_b(1:50,:) + 30*sin(x(1:50,:));
    common = {'pixel_spacing_m', 1, 'seed', 4, ...
        'terrain_normalization_model', 'absolute_scales', ...
        'slope_scale_deg', 45, 'roughness_scale_m', 10, ...
        'curvature_scale_per_m', 0.1, ...
        'shared_residual_std', 0, 'edge_residual_std', 0};
    [~, terrain_a] = gen_coherence_map_multibaseline_terrainaware( ...
        dem_a, 1, common{:});
    [~, terrain_b] = gen_coherence_map_multibaseline_terrainaware( ...
        dem_b, 1, common{:});

    verifyEqual(testCase, terrain_a.terrain_quality(100,64), ...
        terrain_b.terrain_quality(100,64), 'AbsTol', 1e-12);
end

function testEdgeResidualsAreCorrelatedButNotIdentical(testCase)
    [x, y] = meshgrid(1:72, 1:72);
    dem = 0.03*x + 0.02*y + 0.4*sin(x/9);
    maps = gen_coherence_map_multibaseline_terrainaware(dem, [1, 1], ...
        'pixel_spacing_m', 1, 'seed', 21, 'decay_alpha', 0, ...
        'shared_residual_std', 0.18, 'edge_residual_std', 0.05);
    rho = corrcoef(maps{1}(:), maps{2}(:));

    verifyGreaterThan(testCase, rho(1,2), 0.5);
    verifyLessThan(testCase, rho(1,2), 0.99999);
end

function testLocalSeedDoesNotResetCallerStream(testCase)
    rng(456);
    expected_first = rand;
    expected_second = rand;

    rng(456);
    actual_first = rand;
    gen_coherence_map_multibaseline_terrainaware(zeros(24), [1, 2], ...
        'seed', 88);
    actual_second = rand;

    verifyEqual(testCase, actual_first, expected_first);
    verifyEqual(testCase, actual_second, expected_second);
end

function testSteepFacingSlopeCreatesLayoverProxy(testCase)
    [~, y] = meshgrid(1:64, 1:64);
    dem = 2*y;
    [~, terrain] = gen_coherence_map_multibaseline_terrainaware(dem, 1, ...
        'pixel_spacing_m', 1, 'incidence_angle_deg', 45, ...
        'look_azimuth_deg', 90, 'seed', 3);

    verifyGreaterThan(testCase, mean(terrain.layover_mask(:)), 0.9);
end

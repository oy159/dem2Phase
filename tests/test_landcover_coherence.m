function tests = test_landcover_coherence
%TEST_LANDCOVER_COHERENCE Tests categorical surface coherence factors.
    tests = functiontests(localfunctions);
end

function testWorldCoverThreeDegreeTileNaming(testCase)
    verifyEqual(testCase, worldcover_tile_token(27, 85), 'N27E084');
    verifyEqual(testCase, worldcover_tile_token(28, 83), 'N27E081');
    verifyEqual(testCase, worldcover_tile_token(32, 78), 'N30E078');
    verifyEqual(testCase, worldcover_tile_token(46, 9), 'N45E009');
    verifyEqual(testCase, worldcover_tile_token(-33, -70), 'S33W072');
end

function testAllSplitDemsHaveAlignedWorldCover(testCase)
    here = fileparts(mfilename('fullpath'));
    project_root = fullfile(here, '..');
    manifest_path = fullfile(project_root, 'data', 'landcover', ...
        'worldcover_alignment_manifest.csv');
    records = readtable(manifest_path, 'TextType', 'string', 'Delimiter', ',');
    verifyEqual(testCase, height(records), 19);
    verifyEqual(testCase, numel(unique(records.dem_file)), 19);
    verifyEqual(testCase, nnz(records.split == "train"), 12);
    verifyEqual(testCase, nnz(records.split == "test"), 5);
    verifyEqual(testCase, nnz(records.split == "validation"), 2);
    for idx = 1:height(records)
        aligned_path = fullfile(project_root, records.aligned_file(idx));
        verifyEqual(testCase, exist(aligned_path, 'file'), 2);
        data = load(aligned_path, 'landcover_codes', ...
            'alignment_method', 'wc_tile');
        verifySize(testCase, data.landcover_codes, ...
            [records.rows(idx), records.cols(idx)]);
        verifyEqual(testCase, data.alignment_method, ...
            'nearest class at DEM cell center');
        verifyEqual(testCase, string(data.wc_tile), ...
            records.worldcover_tile(idx));
    end
end

function testPrioritySamplingUsesLargestEligibleMultiplier(testCase)
    cfg = nominal_config();
    codes = uint8(40*ones(10));
    codes(1:2,:) = 50;
    codes(1,1:2) = 80;
    [multiplier, code] = landcover_sampling_multiplier( ...
        codes, cfg.dataset.landcover.sampling);

    verifyEqual(testCase, multiplier, 3.0);
    verifyEqual(testCase, code, 80);
end

function testSamplingDoesNotInventAbsentClass(testCase)
    cfg = nominal_config();
    codes = uint8(10*ones(20));
    [multiplier, code] = landcover_sampling_multiplier( ...
        codes, cfg.dataset.landcover.sampling);

    verifyEqual(testCase, multiplier, 1.0);
    verifyEqual(testCase, code, 0);
end

function testWaterAndForestFactorsAreApplied(testCase)
    cfg = nominal_config();
    codes = uint8([60 10; 80 255]);
    input = {ones(2)*0.8};
    [output, info] = apply_landcover_coherence( ...
        input, codes, cfg.dataset.landcover, 0.01);

    verifyEqual(testCase, output{1}(1,1), 0.8, 'AbsTol', 1e-12);
    verifyEqual(testCase, output{1}(1,2), 0.8*0.85, 'AbsTol', 1e-12);
    verifyEqual(testCase, output{1}(2,1), 0.8*0.08, 'AbsTol', 1e-12);
    verifyEqual(testCase, output{1}(2,2), 0.8, 'AbsTol', 1e-12);
    verifyEqual(testCase, info.water_fraction, 0.25);
    verifyFalse(testCase, info.factors_are_calibrated);
end

function testMinimumCoherenceFloor(testCase)
    cfg = nominal_config();
    output = apply_landcover_coherence({ones(3)*0.02}, ...
        uint8(ones(3)*80), cfg.dataset.landcover, 0.01);
    verifyEqual(testCase, output{1}, ones(3)*0.01, 'AbsTol', 1e-12);
end

function testAlignedWorldCoverHasPlausibleClasses(testCase)
    cfg = nominal_config();
    here = fileparts(mfilename('fullpath'));
    project_root = fileparts(here);
    data = load(fullfile(project_root, ...
        char(cfg.dataset.landcover.aligned_file)), ...
        'landcover_codes', 'source_name', 'source_license');
    present = unique(data.landcover_codes(:))';

    verifySize(testCase, data.landcover_codes, [3600, 3600]);
    verifyTrue(testCase, all(ismember(present, ...
        uint8(cfg.dataset.landcover.class_codes))));
    verifyTrue(testCase, ismember(uint8(10), present));
    verifyTrue(testCase, ismember(uint8(80), present));
    verifyEqual(testCase, data.source_name, 'ESA WorldCover 10m 2021 v200');
    verifyEqual(testCase, data.source_license, 'CC BY 4.0');
end

function cfg = nominal_config()
    here = fileparts(mfilename('fullpath'));
    cfg = load_sim_config(fullfile(here, '..', 'configs', ...
        'uav_p_500m_monostatic.json'));
end

function tests = test_patch_group
%TEST_PATCH_GROUP Tests grouped multi-edge storage and observation masks.
    tests = functiontests(localfunctions);
end

function testBuildsKHWAndZeroFillsInactiveEdges(testCase)
    cfg = nominal_config();
    k = numel(cfg.interferometry.baseline_perp_m);
    fields = arrayfun(@(idx) ones(8, 7)*idx, 1:k, 'UniformOutput', false);
    mask = logical([1, 0, 1, 0]);
    terrain.slope_deg = ones(8, 7);
    terrain.layover_mask = false(8, 7);
    metadata.patch_name = 'unit_patch';

    group = build_patch_group(fields, fields, fields, fields, fields, ...
        mask, cfg, terrain, metadata);

    verifySize(testCase, group.wrappedphase_withnoise, [k, 8, 7]);
    verifyClass(testCase, group.wrappedphase_withnoise, 'single');
    verifyEqual(testCase, group.valid_edge_mask, uint8(mask));
    verifyEqual(testCase, squeeze(group.wrappedphase_withnoise(2,:,:)), ...
        zeros(8, 7, 'single'));
    verifyEqual(testCase, squeeze(group.wrappedphase_withnoise(3,:,:)), ...
        ones(8, 7, 'single')*3);
    verifyEqual(testCase, group.edge_index, uint16([1 1 1 1; 2 3 4 5]));
end

function testVariableEdgeMaskKeepsShortestEdge(testCase)
    cfg = nominal_config();
    rng(14);
    observed_counts = zeros(1, 30);
    for idx = 1:numel(observed_counts)
        mask = sample_valid_edge_mask(4, 1, cfg.interferometry.edge_sampling);
        verifyTrue(testCase, mask(1));
        observed_counts(idx) = nnz(mask);
    end
    verifyGreaterThanOrEqual(testCase, min(observed_counts), 2);
    verifyLessThanOrEqual(testCase, max(observed_counts), 4);
    verifyGreaterThan(testCase, numel(unique(observed_counts)), 1);
end

function testEdgeMaskSeedReplaysWithoutFreezingBatchStream(testCase)
    cfg = nominal_config();
    rng(123);
    before = rng;
    mask_a = sample_valid_edge_mask(4, 1, cfg.interferometry.edge_sampling, 99);
    after_local_call = rng;
    mask_b = sample_valid_edge_mask(4, 1, cfg.interferometry.edge_sampling, 99);

    verifyEqual(testCase, mask_a, mask_b);
    verifyEqual(testCase, after_local_call, before);
    batch_draw_a = rand(1, 8);
    batch_draw_b = rand(1, 8);
    verifyNotEqual(testCase, batch_draw_a, batch_draw_b);
end

function testRoundTripStandardMatFile(testCase)
    cfg = nominal_config();
    fields = repmat({rand(5, 6)}, 1, 4);
    terrain.slope_deg = rand(5, 6);
    metadata.patch_name = 'roundtrip';
    group = build_patch_group(fields, fields, fields, fields, fields, ...
        true(1, 4), cfg, terrain, metadata);
    temp_file = [tempname, '.mat'];
    cleanup = onCleanup(@() delete_if_present(temp_file));
    save_patch_group(temp_file, group);
    loaded = load(temp_file);

    verifyEqual(testCase, loaded.valid_edge_mask, uint8(ones(1, 4)));
    verifySize(testCase, loaded.coherence_true, [4, 5, 6]);
    verifyEqual(testCase, loaded.metadata.patch_name, 'roundtrip');
    clear cleanup
end

function cfg = nominal_config()
    here = fileparts(mfilename('fullpath'));
    cfg = load_sim_config(fullfile(here, '..', 'configs', ...
        'uav_p_500m_monostatic.json'));
end

function delete_if_present(filename)
    if exist(filename, 'file')
        delete(filename);
    end
end

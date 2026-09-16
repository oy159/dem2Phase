function tests = test_training_loader
%TEST_TRAINING_LOADER Tests strict masks and balanced sampling weights.
    tests = functiontests(localfunctions);
end

function testInverseFrequencyWeightsEqualizeStrata(testCase)
    strata = [repmat("common", 9, 1); "rare"];
    [weights, counts] = compute_balanced_sample_weights(strata, 20);
    verifyEqual(testCase, height(counts), 2);
    verifyGreaterThan(testCase, weights(end), weights(1));
    verifyEqual(testCase, sum(weights(strata == "common")), ...
        sum(weights(strata == "rare")), 'AbsTol', 1e-12);
    verifyEqual(testCase, mean(weights), 1, 'AbsTol', 1e-12);
end

function testStrictLoaderRejectsLeakyCoherence(testCase)
    group = minimal_group();
    group.metadata.coherence_observation_method = 'legacy_clean_phase_residual';
    filename = [tempname, '.mat'];
    cleanup = onCleanup(@() delete_if_present(filename));
    save(filename, '-struct', 'group');
    verifyError(testCase, @() load_patch_group_for_training(filename), ...
        'dem2phase:LeakyCoherenceInput');
    clear cleanup
end

function testStrictLoaderBuildsCombinedPixelMask(testCase)
    group = minimal_group();
    filename = [tempname, '.mat'];
    cleanup = onCleanup(@() delete_if_present(filename));
    save(filename, '-struct', 'group');
    sample = load_patch_group_for_training(filename);
    verifySize(testCase, sample.valid_pixel_mask, [2,4,3]);
    verifyFalse(testCase, any(sample.valid_pixel_mask(2,:,:), 'all'));
    verifyEqual(testCase, nnz(sample.valid_pixel_mask), 11);
    clear cleanup
end

function group = minimal_group()
    group.wrappedphase_withnoise = zeros(2,4,3,'single');
    group.wrappedphase_withnoise(1,:,:) = 0.5;
    group.coherence_observed = zeros(2,4,3,'single');
    group.coherence_observed(1,:,:) = 0.8;
    group.valid_edge_mask = uint8([1,0]);
    group.coregistration_valid_mask = zeros(2,4,3,'uint8');
    group.coregistration_valid_mask(1,:,:) = 1;
    group.coregistration_valid_mask(1,1,1) = 0;
    group.edge_index = uint16([1,1;2,3]);
    group.baseline_perp_m = single([0.75,1.2]);
    group.ambiguity_height_m = single([200,125]);
    group.failure_labels.compound_failure_present = false;
    group.metadata.coherence_observation_method = 'observed_node_slc_pair';
    group.metadata.dataset_split = 'train';
    group.metadata.split_group = 'N00E000';
    group.metadata.geographic_tile = 'N00E000';
    group.metadata.geographic_region = 'test_region';
    group.metadata.source_dataset = 'test_source';
    group.metadata.source_quality_role = 'primary';
    group.unwrapped_phase = zeros(2,4,3,'single');
    group.coherence_true = zeros(2,4,3,'single');
end

function delete_if_present(filename)
    if exist(filename, 'file')
        delete(filename);
    end
end

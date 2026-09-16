function tests = test_dem_split_manifest
%TEST_DEM_SPLIT_MANIFEST Whole-DEM split validation tests.
    tests = functiontests(localfunctions);
end

function testProjectManifestHasDisjointCompleteSplits(testCase)
    here = fileparts(mfilename('fullpath'));
    path = fullfile(here, '..', 'configs', 'dem_split_manifest.csv');
    table = load_dem_split_manifest(path);
    verifyEqual(testCase, height(table), 19);
    verifyEqual(testCase, nnz(table.split == "train"), 12);
    verifyEqual(testCase, nnz(table.split == "test"), 5);
    verifyEqual(testCase, nnz(table.split == "validation"), 2);
    verifyEqual(testCase, sum(table.patches_per_dem(table.split == "validation")), 16);
    verifyEqual(testCase, sum(table.patches_per_dem), 234);
    verifyEqual(testCase, numel(unique(lower(table.dem_file))), height(table));
end

function testRejectsOneDemInTwoSplits(testCase)
    path = [tempname '.csv'];
    cleanup = onCleanup(@() delete_if_present(path));
    fid = fopen(path, 'w');
    fprintf(fid, ['dem_file,geographic_tile,split,region,quality_role,' ...
        'patches_per_dem,latitude_deg,longitude_deg,source\n' ...
        'a.tif,N00E000,train,a,primary,1,0,0,x\n' ...
        'a.tif,N00E000,test,a,primary,1,0,0,x\n']);
    fclose(fid);
    verifyError(testCase, @() load_dem_split_manifest(path), ...
        'dem2phase:DuplicateDemSplit');
end

function testRejectsCrossSensorGeographicLeakage(testCase)
    path = [tempname '.csv'];
    cleanup = onCleanup(@() delete_if_present(path));
    fid = fopen(path, 'w');
    fprintf(fid, ['dem_file,geographic_tile,split,region,quality_role,' ...
        'patches_per_dem,latitude_deg,longitude_deg,source\n' ...
        'alos.tif,N27E088,train,a,secondary,1,27,88,alos\n' ...
        'cop.tif,N27E088,test,a,primary,1,27,88,cop\n']);
    fclose(fid);
    verifyError(testCase, @() load_dem_split_manifest(path), ...
        'dem2phase:GeographicSplitLeakage');
end

function delete_if_present(path)
    if exist(path, 'file'); delete(path); end
end

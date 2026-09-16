%% prepare_worldcover_for_dem.m
% Align the correct ESA WorldCover categorical tile to every configured DEM.
% The split manifest is the source of truth. Categorical codes are sampled at
% DEM cell centres with nearest-neighbour indexing; they are never interpolated.
clc;

script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
addpath(genpath(project_root));
cfg = load_sim_config(fullfile(project_root, 'configs', ...
    'uav_p_500m_monostatic.json'));
split_path = resolve_path(project_root, char(cfg.dataset.split_manifest));
assignments = load_dem_split_manifest(split_path);
dem_root = resolve_path(project_root, char(cfg.dataset.dem_directory));
dem_files = dir(fullfile(dem_root, '**', '*.tif'));
raw_root = fullfile(project_root, 'data', 'landcover', 'raw');
aligned_root = resolve_path(project_root, ...
    char(cfg.dataset.landcover.aligned_directory));
if ~exist(aligned_root, 'dir'); mkdir(aligned_root); end

records = table('Size', [height(assignments), 11], ...
    'VariableTypes', {'string','string','string','string','string','string', ...
    'double','double','double','double','string'}, ...
    'VariableNames', {'dem_file','split','region','worldcover_tile', ...
    'raw_file','aligned_file','rows','cols','unknown_fraction', ...
    'boundary_clamped_fraction','codes_present'});

for idx = 1:height(assignments)
    dem_name = char(assignments.dem_file(idx));
    matches = find(strcmpi({dem_files.name}, dem_name));
    assert(isscalar(matches), 'dem2phase:DemSplitLookupFailed', ...
        'Expected exactly one DEM named %s.', dem_name);
    dem_path = fullfile(dem_files(matches).folder, dem_files(matches).name);
    wc_tile = worldcover_tile_token(assignments.latitude_deg(idx), ...
        assignments.longitude_deg(idx));
    wc_name = sprintf('ESA_WorldCover_10m_2021_v200_%s_Map.tif', wc_tile);
    worldcover_path = fullfile(raw_root, wc_name);
    assert(exist(worldcover_path, 'file') == 2, ...
        'dem2phase:MissingLandcoverInput', ...
        'Missing WorldCover tile %s for DEM %s.', wc_name, dem_name);

    dem_info = georasterinfo(dem_path);
    worldcover_info = georasterinfo(worldcover_path);
    dem_ref = dem_info.RasterReference;
    worldcover_ref = worldcover_info.RasterReference;

    dem_cols = 1:dem_ref.RasterSize(2);
    dem_rows = 1:dem_ref.RasterSize(1);
    [~, target_lon] = intrinsicToGeographic(dem_ref, dem_cols, ...
        ones(size(dem_cols)));
    [target_lat, ~] = intrinsicToGeographic(dem_ref, ...
        ones(size(dem_rows)), dem_rows);
    assert(strcmpi(worldcover_ref.RowsStartFrom, 'west') && ...
        strcmpi(worldcover_ref.ColumnsStartFrom, 'north'), ...
        'dem2phase:UnsupportedWorldCoverOrientation', ...
        'WorldCover reference must increase eastward and southward.');
    % Compute the separable geographic-cell transform explicitly. MATLAB's
    % geographicToIntrinsic may wrap an exact western-boundary longitude by
    % 360 degrees when it appears in an ASTER point-grid vector.
    source_x = 0.5 + (target_lon-worldcover_ref.LongitudeLimits(1)) / ...
        worldcover_ref.CellExtentInLongitude;
    source_y = 0.5 + (worldcover_ref.LatitudeLimits(2)-target_lat) / ...
        worldcover_ref.CellExtentInLatitude;
    % Point-registered 3601x3601 ASTER tiles include a shared boundary row
    % or column. Permit at most two source pixels outside the half-open
    % WorldCover tile and clamp only that shared edge.
    boundary_tolerance_px = 2;
    assert(all(source_x >= 0.5-boundary_tolerance_px & ...
        source_x <= worldcover_ref.RasterSize(2)+0.5+boundary_tolerance_px) && ...
        all(source_y >= 0.5-boundary_tolerance_px & ...
        source_y <= worldcover_ref.RasterSize(1)+0.5+boundary_tolerance_px), ...
        'dem2phase:WorldCoverDoesNotContainDem', ...
        'WorldCover tile %s does not fully contain DEM %s.', wc_tile, dem_name);
    boundary_clamped_fraction = (nnz(source_x < 0.5 | ...
        source_x > worldcover_ref.RasterSize(2)+0.5) + ...
        nnz(source_y < 0.5 | source_y > worldcover_ref.RasterSize(1)+0.5)) / ...
        (numel(source_x) + numel(source_y));
    source_cols = min(max(round(source_x), 1), worldcover_ref.RasterSize(2));
    source_rows = min(max(round(source_y), 1), worldcover_ref.RasterSize(1));

    row_bounds = [min(source_rows), max(source_rows)];
    col_bounds = [min(source_cols), max(source_cols)];
    subset = imread(worldcover_path, 'PixelRegion', ...
        {[row_bounds(1), row_bounds(2)], [col_bounds(1), col_bounds(2)]});
    subset = uint8(subset);
    landcover_codes = subset(source_rows-row_bounds(1)+1, ...
        source_cols-col_bounds(1)+1);
    assert(isequal(size(landcover_codes), dem_ref.RasterSize), ...
        'dem2phase:LandcoverAlignmentShapeMismatch', ...
        'Aligned WorldCover raster must match DEM %s.', dem_name);

    source_name = 'ESA WorldCover 10m 2021 v200';
    source_url = ['https://esa-worldcover.s3.eu-central-1.amazonaws.com/' ...
        'v200/2021/map/' wc_name];
    source_license = 'CC BY 4.0';
    attribution = ['© ESA WorldCover project 2021 / Contains modified ' ...
        'Copernicus Sentinel data (2021) processed by ESA WorldCover consortium'];
    alignment_method = 'nearest class at DEM cell center';
    [~, dem_stem] = fileparts(dem_name);
    aligned_name = [dem_stem char(cfg.dataset.landcover.aligned_suffix)];
    output_path = fullfile(aligned_root, aligned_name);
    save(output_path, 'landcover_codes', 'dem_ref', 'source_name', ...
        'source_url', 'source_license', 'attribution', ...
        'alignment_method', 'wc_tile', '-v7');

    codes = unique(landcover_codes(:))';
    records.dem_file(idx) = string(dem_name);
    records.split(idx) = assignments.split(idx);
    records.region(idx) = assignments.region(idx);
    records.worldcover_tile(idx) = string(wc_tile);
    records.raw_file(idx) = string(fullfile('data', 'landcover', 'raw', wc_name));
    records.aligned_file(idx) = string(fullfile('data', 'landcover', aligned_name));
    records.rows(idx) = size(landcover_codes, 1);
    records.cols(idx) = size(landcover_codes, 2);
    records.unknown_fraction(idx) = mean(landcover_codes(:) == 0);
    records.boundary_clamped_fraction(idx) = boundary_clamped_fraction;
    records.codes_present(idx) = string(mat2str(codes));
    fprintf('[%d/%d] %s <- %s, codes=%s\n', idx, height(assignments), ...
        dem_name, wc_tile, mat2str(codes));
end

manifest_path = fullfile(aligned_root, 'worldcover_alignment_manifest.csv');
writetable(records, manifest_path);
fprintf('Aligned %d DEMs. Manifest: %s\n', height(records), manifest_path);

function output = resolve_path(root, configured)
    if ~isempty(regexp(configured, '^[A-Za-z]:[\\/]', 'once')) || ...
            startsWith(configured, '/') || startsWith(configured, '\\')
        output = configured;
    else
        output = fullfile(root, configured);
    end
end

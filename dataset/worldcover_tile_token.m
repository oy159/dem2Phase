function token = worldcover_tile_token(latitude_deg, longitude_deg)
%WORLDCOVER_TILE_TOKEN Return the 3-degree tile containing a coordinate.
% The token identifies the lower-left corner using WorldCover naming.
    validateattributes(latitude_deg, {'numeric'}, ...
        {'scalar','finite','>=',-90,'<=',90});
    validateattributes(longitude_deg, {'numeric'}, ...
        {'scalar','finite','>=',-180,'<=',180});
    lat0 = 3*floor(double(latitude_deg)/3);
    lon0 = 3*floor(double(longitude_deg)/3);
    if lat0 >= 0; lat_h = 'N'; else; lat_h = 'S'; end
    if lon0 >= 0; lon_h = 'E'; else; lon_h = 'W'; end
    token = sprintf('%s%02d%s%03d', lat_h, abs(lat0), lon_h, abs(lon0));
end

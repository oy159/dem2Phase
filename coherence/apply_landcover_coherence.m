function [adjusted_maps, info] = apply_landcover_coherence( ...
    coherence_maps, landcover_codes, landcover_cfg, min_coherence)
%APPLY_LANDCOVER_COHERENCE Apply categorical surface factors to coherence.
%
% Factors are explicitly configuration-driven and are not universal P-band
% constants. Unknown codes retain unknown_factor and provenance is preserved.

    assert(iscell(coherence_maps) && ~isempty(coherence_maps), ...
        'dem2phase:InvalidCoherenceStack', ...
        'coherence_maps must be a non-empty cell array.');
    codes = double(landcover_cfg.class_codes(:)');
    factors = double(landcover_cfg.coherence_factors(:)');
    factor_map = ones(size(landcover_codes))*double(landcover_cfg.unknown_factor);
    for idx = 1:numel(codes)
        factor_map(landcover_codes == codes(idx)) = factors(idx);
    end

    adjusted_maps = cell(size(coherence_maps));
    for edge_idx = 1:numel(coherence_maps)
        assert(isequal(size(coherence_maps{edge_idx}), size(landcover_codes)), ...
            'dem2phase:LandcoverCoherenceShapeMismatch', ...
            'landcover_codes must match every coherence raster.');
        adjusted_maps{edge_idx} = max(min_coherence, ...
            coherence_maps{edge_idx} .* factor_map);
    end

    present_codes = unique(landcover_codes(:))';
    fractions = arrayfun(@(code) mean(landcover_codes(:) == code), present_codes);
    [~, dominant_idx] = max(fractions);
    info.codes = uint8(landcover_codes);
    info.factor_map = factor_map;
    info.present_codes = present_codes;
    info.fractions = fractions;
    info.dominant_code = present_codes(dominant_idx);
    info.water_fraction = mean(landcover_codes(:) == 80);
    info.source = char(landcover_cfg.source);
    info.profile = char(landcover_cfg.profile);
    info.factors_are_calibrated = false;
end

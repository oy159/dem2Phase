function [multiplier, priority_code] = landcover_sampling_multiplier(codes, sampling)
%LANDCOVER_SAMPLING_MULTIPLIER Boost acceptance of configured surface types.
%
% The largest eligible multiplier is used. This changes only the probability
% of accepting a real crop; it never paints or randomly assigns a class.

    multiplier = 1;
    priority_code = 0;
    if ~logical(sampling.enabled)
        return
    end

    priority_codes = double(sampling.priority_codes(:)');
    minimum_fraction = double(sampling.minimum_fraction(:)');
    multipliers = double(sampling.acceptance_multiplier(:)');
    assert(numel(priority_codes) == numel(minimum_fraction) && ...
        numel(priority_codes) == numel(multipliers), ...
        'dem2phase:InvalidLandcoverSampling', ...
        'Land-cover sampling arrays must have equal lengths.');

    fractions = arrayfun(@(code) mean(codes(:) == code), priority_codes);
    eligible = fractions >= minimum_fraction;
    if any(eligible)
        eligible_indices = find(eligible);
        [multiplier, local_idx] = max(multipliers(eligible));
        priority_code = priority_codes(eligible_indices(local_idx));
    end
end

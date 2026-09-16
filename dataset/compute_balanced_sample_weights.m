function [weights, counts] = compute_balanced_sample_weights(strata, max_weight_ratio)
%COMPUTE_BALANCED_SAMPLE_WEIGHTS Inverse-frequency weights with clipping.
    if nargin < 2
        max_weight_ratio = 10;
    end
    strata = string(strata(:));
    assert(~isempty(strata) && max_weight_ratio >= 1, ...
        'dem2phase:InvalidBalancedSamplingInput', ...
        'Strata must be nonempty and max_weight_ratio must be >= 1.');
    [unique_strata, ~, group_idx] = unique(strata, 'stable');
    group_counts = accumarray(group_idx, 1);
    raw = numel(strata) ./ (numel(unique_strata) * group_counts(group_idx));
    raw = min(max(raw, 1/max_weight_ratio), max_weight_ratio);
    weights = raw/mean(raw);
    counts = table(unique_strata, group_counts, ...
        'VariableNames', {'stratum', 'sample_count'});
end

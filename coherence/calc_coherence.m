function [coh_map_est, coh_mean] = calc_coherence(input1, input2, varargin)
%CALC_COHERENCE  Estimate coherence for verifying InSAR simulation.
%
%   [coh_map_est, coh_mean] = CALC_COHERENCE(slc1, slc2)
%   [coh_map_est, coh_mean] = CALC_COHERENCE(slc1, slc2, 'win_size', 7)
%   [coh_map_est, coh_mean] = CALC_COHERENCE(noisy_phase, clean_phase, ...
%                                             'method', 'phase')
%
%   Inputs:
%     input1, input2   Complex SLC pair  (method = 'slc',   default)
%                   OR real noisy/clean phase images (method = 'phase')
%
%   Name-Value Parameters:
%     'method'    'slc'   — coherence from complex SLC pair   (default)
%                 'phase' — coherence from phase residual variance
%     'win_size'  estimation window size (pixels)             (default 5)
%
%   Outputs:
%     coh_map_est  [M x N] estimated coherence map in [0, 1]
%     coh_mean     scalar mean coherence over the whole image
%
%   -----------------------------------------------------------------
%   Method 'slc' (standard InSAR coherence estimator):
%
%     gamma(x,y) = |<s1 * conj(s2)>_w| / sqrt(<|s1|^2>_w * <|s2|^2>_w)
%
%     where <.>_w denotes spatial averaging over a win_size x win_size window.
%
%   Method 'phase' (verification without SLC data):
%
%     1. Compute the wrapped phase residual:
%           r(x,y) = angle(exp(i * (noisy_phase - clean_phase)))
%     2. Estimate the local variance of r over the estimation window.
%     3. Invert the variance through a pre-computed look-up table derived
%        from the same PDF as patch_add_noise.m:
%
%           sigma^2(gamma) = integral_{-pi}^{pi} phi^2 * p(phi; gamma) dphi
%
%     The inversion gives an estimate of the local coherence gamma.
%     Since this uses a single noise realisation the estimate has a
%     standard deviation of roughly sigma / sqrt(win_size^2 - 1).
%   -----------------------------------------------------------------
%
%   Typical verification workflow:
%     coh_in      = gen_coherence_map(256, 256, 'seed', 1);
%     noisy       = patch_add_noise_spatialcoh(clean_phase, coh_in, 'seed', 42);
%     [coh_est, ~] = calc_coherence(noisy, clean_phase, 'method', 'phase', ...
%                                   'win_size', 11);
%     figure; subplot(1,2,1); imagesc(coh_in);  title('Input coherence');
%             subplot(1,2,2); imagesc(coh_est); title('Estimated coherence');
%             colormap jet; colorbar;
%
%   See also: gen_coherence_map, patch_add_noise_spatialcoh

    p = inputParser;
    addRequired(p,  'input1',    @isnumeric);
    addRequired(p,  'input2',    @isnumeric);
    addParameter(p, 'method',   'slc', @(x) ismember(x, {'slc', 'phase'}));
    addParameter(p, 'win_size',  5,    @(x) isnumeric(x) && isscalar(x) && x >= 1);
    parse(p, input1, input2, varargin{:});
    opt = p.Results;

    assert(isequal(size(opt.input1), size(opt.input2)), ...
        'calc_coherence: input1 and input2 must have the same size.');

    switch opt.method
        case 'slc'
            coh_map_est = coherence_from_slc(opt.input1, opt.input2, opt.win_size);
        case 'phase'
            coh_map_est = coherence_from_phase(opt.input1, opt.input2, opt.win_size);
    end

    coh_mean = mean(coh_map_est(:));
end

% =========================================================================
function coh_map = coherence_from_slc(slc1, slc2, win_size)
%   Sliding-window coherence estimation from a complex SLC pair.
    w      = ones(win_size) / win_size^2;
    cross  = conv2(slc1 .* conj(slc2),  w, 'same');
    power1 = conv2(abs(slc1).^2,         w, 'same');
    power2 = conv2(abs(slc2).^2,         w, 'same');
    denom  = sqrt(power1 .* power2);
    coh_map = abs(cross) ./ (denom + eps);
    coh_map = min(max(real(coh_map), 0), 1);
end

% =========================================================================
function coh_map = coherence_from_phase(noisy_phase, clean_phase, win_size)
%   Estimate coherence from the variance of the wrapped phase noise residual.
    residual = angle(exp(1i * (noisy_phase - clean_phase)));

    % Local mean and second moment of the residual
    w    = ones(win_size) / win_size^2;
    mu   = conv2(residual,      w, 'same');
    mu2  = conv2(residual.^2,   w, 'same');
    var_r = max(mu2 - mu.^2, 0);

    % Pre-compute variance vs coherence LUT (uses same PDF as patch_add_noise.m)
    [var_lut, coh_lut] = phase_variance_lut();

    % Invert local variance to coherence estimate
    coh_flat  = interp1(var_lut, coh_lut, var_r(:)', 'linear', 'extrap');
    coh_map   = reshape(min(max(coh_flat, 0.01), 0.99), size(noisy_phase));
end

% =========================================================================
function [var_lut, coh_lut] = phase_variance_lut()
%   Pre-compute a monotone variance-coherence look-up table.
%   Uses the SAME PDF formula as patch_add_noise.m (unchanged).

    coh_lut  = linspace(0.01, 0.99, 200);
    var_lut  = zeros(1, 200);

    phi_range = 2 * pi;
    phi_axis  = linspace(-phi_range/2, phi_range/2, 1000);
    d_phi     = phi_axis(2) - phi_axis(1);

    for k = 1:numel(coh_lut)
        gamma   = coh_lut(k);

        % ---- PDF (preserved from patch_add_noise.m) ------------------
        cos_phi  = cos(phi_axis);
        denom    = max(1 - gamma^2 * cos_phi.^2, eps);
        pdf_k    = ((1 - gamma^2) / (2 * pi)) ...
                   .* (1 ./ denom) ...
                   .* (1 + (gamma .* cos_phi .* acos(-gamma .* cos_phi)) ...
                       ./ sqrt(denom));
        pdf_k    = max(pdf_k, 0);
        pdf_k    = pdf_k / (sum(pdf_k) * d_phi);   % normalise to unit integral

        var_lut(k) = sum(phi_axis.^2 .* pdf_k) * d_phi;
    end

    % Ensure the LUT is strictly monotone (higher coh -> lower variance).
    % Sort by ascending variance so interp1 works correctly.
    [var_lut, sort_idx] = sort(var_lut, 'ascend');
    coh_lut = coh_lut(sort_idx);
end

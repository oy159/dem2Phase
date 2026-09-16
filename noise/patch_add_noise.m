function current_patch_withnoise = patch_add_noise(current_patch , coherence)

%% 等效视数为1
% 计算pdf
phi_range =  2 * pi; % 相位范围
phi_axis = linspace(-phi_range/2 , phi_range/2 , 1000);
phase_noise_pdf = ((1-coherence^2) ./ (2 * pi)) .* (1 ./ (1-coherence^2 * cos(phi_axis).^2)) ...
                    .* (1 + (coherence .* cos(phi_axis) .* acos(-coherence .* cos(phi_axis))) ./ sqrt(1 - coherence^2 .* cos(phi_axis).^2));
% figure;
% plot(phi_axis , phase_noise_pdf);
% title('相位噪声概率分布函数（等效视数为1）');
% grid on;axis tight;

% 计算cdf
phase_noise_cdf = cumsum(phase_noise_pdf); % calculate the cumulative sum of the sorted pdf values
phase_noise_cdf = phase_noise_cdf / sum(phase_noise_pdf); % normalize the cdf
% figure;
% plot(phi_axis , phase_noise_cdf);
% title('相位噪声概率累积分布函数（等效视数为1）');
% grid on;axis tight;

% 计算icdf
load random_generate_state.mat
rng(s);
quantiles = rand(1 , size(current_patch , 1) * size(current_patch , 2));
phase_noise_icdf = interp1(phase_noise_cdf, phi_axis, quantiles);
phase_noise_icdf(find(isnan(phase_noise_icdf) == 1)) = 0;
if sum(isnan(phase_noise_icdf)) > 0
    fprintf('nan_error\n');
end
current_patch_withnoise = reshape(phase_noise_icdf , size(current_patch , 1) , size(current_patch , 2)) + current_patch;
current_patch_withnoise = angle(exp(1i * (current_patch_withnoise)));
% figure;
% imagesc(current_patch_withnoise); axis tight equal; colormap jet;
end






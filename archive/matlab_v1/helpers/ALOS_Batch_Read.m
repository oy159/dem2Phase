%%
clc;
clear all;
close all;

base_path = pwd;
sub_dir = dir([base_path , '\dem\']);
folder_name = {sub_dir.name};
lambda =[0.236,0.031, 0.005 ];
baseline = [230, 237, 16];
slant_range = [868142, 640708, 845521];
theta = [38.77, 41.07, 39.67];
ratio_dem2phase = 4 * pi ./ lambda .* baseline ./ slant_range ./ sind(theta); % 高程到相位转换比率
save_path = 'E:\trans_unet\data\training_dataset_temp';

patch_size = 256; % 划分的patch大小
stride = 650; % 滑动步长
% 训练集用的数字是210  验证集 1000
% stride:550 对应5292个样本 stride:1650 对应972个样本

tic
m = 3600; % 单张dem尺寸
n = 3600; 
row_patch_number = floor((m - patch_size) / stride) + 1;
col_patch_number = floor((n - patch_size) / stride) + 1;
number_patches = row_patch_number * col_patch_number * length(0.4 : 0.05 : 0.95) * 3 * 10; % 总切片数量
coherence = 0.4; % 相位相干性
maximum_wrap_count_number = 30;
wrap_count_statistics = zeros(3, 3, row_patch_number, col_patch_number, 3);
temp = 0;
% 对dem下面各个文件夹中的原始DEM数据进行裁剪
temp = [];
% type_options = {'wrapped_phase'}; % 不同类型patch                 
type_options = {'wrapped_phase', 'wrapped_phase_real', ...
    'wrapped_phase_imag', 'wrapped_phase_withnoise', ...
    'unwrapped_phase', 'wrap_count', 'wrapped_phase_gradient',...
    'wrapped_phase_gradient_vertical', 'wrapped_phase_gradient_horizontal', ...      
    'wrapped_phase_gradient_horizontal', 'wrapped_phase_gradient_vertical'}; % 不同类型patch     
for ll = 1 : length(type_options)                    
    current_patch_type = type_options{ll};
    row_patch_number = floor((m - patch_size) / stride) + 1;
    col_patch_number = floor((n - patch_size) / stride) + 1;
    row_index = 1 : stride : (row_patch_number * stride);
    col_index = 1 : stride : (col_patch_number * stride);
    save_flag = 0;
    count = 0;
    for kk = 1 : 3 % 遍历系统参数
        for ii = 3 : length(folder_name) % 遍历DSM        
            for jj = 0.4 : 0.05 : 0.95  % 遍历相干性
                coherence = jj;
                current_file_name = [folder_name{ii}];
                dem_name = sprintf('system_%d_DSM_%d__coherence_%.2f' , kk, ii-2, jj)
                dem_path = [base_path, '\dem\' , current_file_name];           
                current_dem = double(imread(dem_path)) * ratio_dem2phase(kk);
                [x_grid, y_grid] = meshgrid(1 : size(current_dem, 2) , 1 : size(current_dem, 1));
                [x_interp_grid , y_interp_grid] = meshgrid(linspace(1 , size(current_dem , 2) , round(3.5 * size(current_dem , 2))) , linspace(1 , size(current_dem , 1) , round(3.5 * size(current_dem , 1))));
                current_dem = interp2(x_grid , y_grid , current_dem , x_interp_grid , y_interp_grid);
                for oo = 1 : length(row_index)
                    for pp = 1 : length(col_index)
                        count = count + 1;
                        current_patch_dem = current_dem(row_index(oo) : (row_index(oo) + patch_size - 1) , col_index(pp) : (col_index(pp) + patch_size - 1));
                        current_patch_dem = current_patch_dem - min(current_patch_dem(:)); 
                        current_patch_wrapcount = round((current_patch_dem - angle(exp(1i * current_patch_dem)))/2/pi); 
                        if  max(current_patch_wrapcount(:)) <= 20                      
                            temp = [temp , max(current_patch_wrapcount(:))];
                            current_patch_wrapped_phase = angle(exp(1i * current_patch_dem));      
                            current_patch_wrapped_vertical_gradient = [diff(current_patch_wrapcount , 1 , 1) ; zeros(1 , 256)];
                            current_patch_wrapped_vertical_gradient = abs(current_patch_wrapped_vertical_gradient);
                            current_patch_wrapped_vertical_gradient_map = zeros(size(current_patch_wrapped_vertical_gradient));
                            current_patch_wrapped_vertical_gradient_map(find(current_patch_wrapped_vertical_gradient >= 1)) = 1;
                            current_patch_wrapped_vertical_gradient_map(find(current_patch_wrapped_vertical_gradient < 1)) = 0;    

                            current_patch_wrapped_horizontal_gradient = [diff(current_patch_wrapped_vertical_gradient_map , 1 , 2) , zeros(256 , 1)];
                            current_patch_wrapped_horizontal_gradient = abs(current_patch_wrapped_horizontal_gradient);
                            current_patch_wrapped_horizontal_gradient_map = zeros(size(current_patch_wrapped_horizontal_gradient));
                            current_patch_wrapped_horizontal_gradient_map(find(current_patch_wrapped_horizontal_gradient >= 1)) = 1;
                            current_patch_wrapped_horizontal_gradient_map(find(current_patch_wrapped_horizontal_gradient < 1)) = 0;            
                            current_patch_wrapped_gradient_map = current_patch_wrapped_vertical_gradient_map + current_patch_wrapped_horizontal_gradient_map;
                            current_patch_wrapped_gradient_map(find(current_patch_wrapped_gradient_map >= 2)) = 0;    
                                switch current_patch_type
                                    case 'wrapped_phase'
                                        current_wrapped_phase_patch = angle(exp(1i * current_patch_dem));
                                        current_save_path = [save_path , '\wrappedphase_withoutnoise\'];
                                        get_patch(current_wrapped_phase_patch , current_save_path , dem_name , coherence , count , current_patch_type);
                                    case 'wrapped_phase_real'
                                        current_wrapped_phase_real_patch = real(exp(1i * current_patch_dem));
                                        current_save_path = [save_path , '\wrappedphase_withoutnoise_real\'];
                                        get_patch(current_wrapped_phase_real_patch , current_save_path , dem_name , coherence , count , current_patch_type);                    
                                    case 'wrapped_phase_imag'
                                        current_wrapped_phase_imag_patch = imag(exp(1i * current_patch_dem));
                                        current_save_path = [save_path , '\wrappedphase_withoutnoise_imag\'];
                                        get_patch(current_wrapped_phase_imag_patch , current_save_path , dem_name , coherence , count , current_patch_type);                                          
                                    case 'wrapped_phase_withnoise'
                                        current_wrapped_phase_patch = angle(exp(1i * current_patch_dem));
                                        current_save_path = {[save_path , '\wrappedphase_withnoise\']; ...
                                                             [save_path , '\wrappedphase_withnoise_real\']; ...
                                                             [save_path , '\wrappedphase_withnoise_imag\']};
                                        get_patch(current_wrapped_phase_patch , current_save_path , dem_name , coherence , count , current_patch_type);                
                                    case 'unwrapped_phase'
                                        current_save_path = [save_path , '\unwrapped_phase'];
                                        get_patch(current_patch_dem , current_save_path , dem_name , coherence , count , current_patch_type);                                         
                                    case 'wrap_count'
                                        current_save_path = [save_path , '\wrapped_number'];                                               
                                        get_patch(current_patch_wrapcount , current_save_path , dem_name , coherence , count , current_patch_type);                                        
                                    case 'wrapped_phase_gradient'
                                        current_save_path = [save_path , '\unwrappedphase_gradient'];    
                                        get_patch(current_patch_wrapped_gradient_map , current_save_path , dem_name , coherence , count , current_patch_type);                                          
                                    case 'wrapped_phase_gradient_vertical'
                                        current_save_path = [save_path , '\unwrappedphase_vertical_gradient'];
                                        get_patch(current_patch_wrapped_vertical_gradient_map , current_save_path , dem_name , coherence , count , current_patch_type);                                         
                                    case 'wrapped_phase_gradient_horizontal'                        
                                        current_save_path = [save_path , '\unwrappedphase_horizontal_gradient'];      
                                        get_patch(current_patch_wrapped_horizontal_gradient_map , current_save_path , dem_name , coherence , count , current_patch_type);                                          
                                end 
                        end
                    end
                end
            end
        end
    end
end
toc

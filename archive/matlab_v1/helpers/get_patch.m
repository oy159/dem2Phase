% input_dem-输入的完整DEM
% current_save_path-当前裁剪patch所属文件夹
% dem_name-当前DEM名字
% patch_size-裁剪图像大小
% stride-裁剪时滑动步长
% coherence-相位噪声相干性
% current_patch_type-当前裁剪patch类型
function get_patch(input, save_path, input_name, coherence, count, patch_type)
    current_patch = input;
    switch patch_type
       case 'wrapped_phase_real'
           wrappedphase_withoutnoise_real_patch = current_patch;
           save([save_path, '\' , input_name , '_wrapped_phase_real_', num2str(count) ,'.mat'] , 'wrappedphase_withoutnoise_real_patch');
       case 'wrapped_phase'
           wrappedphase_withoutnoise_patch = current_patch;
           save([save_path, '\' , input_name , '_wrapped_phase_', num2str(count) ,'.mat'] , 'wrappedphase_withoutnoise_patch');
       case 'wrapped_phase_imag'
           wrappedphase_withoutnoise_imag_patch = current_patch;     
           save([save_path, '\' , input_name , '_wrapped_phase_imag_', num2str(count),'.mat'] , 'wrappedphase_withoutnoise_imag_patch');                   
       case 'wrapped_phase_withnoise'    
            wrappedphasewithnoise_patch = patch_add_noise(current_patch , coherence); % 添加噪声
            save([save_path{1}, '\' , input_name , '_wrapped_phase_withnoise_', num2str(count),'.mat'] , 'wrappedphasewithnoise_patch');                
            wrappedphasewithnoise_real_patch = real(exp(1i * wrappedphasewithnoise_patch)); % 含噪声的缠绕相位图的实部
            save([save_path{2}, '\' , input_name , '_wrapped_phase_withnoise_real_', num2str(count),'.mat'] , 'wrappedphasewithnoise_real_patch');                
            wrappedphasewithnoise_imag_patch = imag(exp(1i * wrappedphasewithnoise_patch)); % 含噪声的缠绕相位图的虚部
            save([save_path{3}, '\' , input_name , '_wrapped_phase_withnoise_imag_', num2str(count),'.mat'] , 'wrappedphasewithnoise_imag_patch');                                                                  
       case 'unwrapped_phase'     
           unwrappedphase_patch = current_patch;               
           save([save_path, '\' , input_name , '_unwrapped_phase_', num2str(count),'.mat'] , 'unwrappedphase_patch');    
       case 'wrap_count'
           wrapped_number_patch = current_patch;
           save([save_path, '\' , input_name , '_wrap_count_', num2str(count), '.mat'] , 'wrapped_number_patch');                 
       case 'wrapped_phase_gradient'      
           wrappedphase_gradient_patch = current_patch;
           save([save_path, '\' , input_name , '_phase_gradient_', num2str(count),'.mat'] , 'wrappedphase_gradient_patch');          
       case 'wrapped_phase_gradient_vertical'
           wrappedphase_vertical_gradient__patch = current_patch;
           save([save_path, '\' , input_name , '_phase_gradient_vertical_', num2str(count),'.mat'] , 'wrappedphase_vertical_gradient__patch');   
       case 'wrapped_phase_gradient_horizontal'
           wrappedphase_horizontal_gradient_patch = current_patch;
           save([save_path, '\' , input_name , '_phase_gradient_horizontal_', num2str(count),'.mat'] , 'wrappedphase_horizontal_gradient_patch');                
    end
end


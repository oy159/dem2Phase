
base_path = 'E:\DEM\ALOS2\RipeData\N035E100_N040E105_PATCH\unwrappedphase';
addpath('E:\Code\TraditionUnwrapper\相位解缠\图割相位解缠');

sub_dir = dir(base_path);
file_name = {sub_dir.name};

for i =3:length(file_name) i

    fn_path = [base_path,'\',file_name{i}];
    load(fn_path);
    absolutephase = unwrappedphase;
    M = 256;
% %%%%%%%%    随机倍数插值处理
    multipliers = [3,4,5];
    randindex_mm = randi([1 length(multipliers)],1);
    multiplier = multipliers(randindex_mm);
    [m,n] = size(absolutephase);
    [x,y]=meshgrid(1:m,1:n);
    xlin = linspace(1,m,multiplier*m);
    ylin = linspace(1,n,multiplier*n);
    [xx,yy]=meshgrid(xlin,ylin);
    absolutephase = interp2(x,y,absolutephase,xx,yy,'bicubic');
% %%%%%%%%    随机裁剪
    [mm,nn] = size(absolutephase);
    mm = mm-M-1;
    nn = nn-M-1;
    randindex_xx = randi([1 mm],1);
    randindex_yy = randi([1 nn],1);
    absolutephase = absolutephase(randindex_xx:randindex_xx+M-1,randindex_yy:randindex_yy+M-1);

% %%%%%%%%    随机生成噪声相位    
    cohlevel = [0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9];
    randindex = randi([1 length(cohlevel)],1);
    co = cohlevel(randindex);
    [x1,x2] = insarpair_v2(ones(M), co*ones(M), absolutephase, 0);	
    
    slc1 = x1;
    slc2 = x2;

    wrappedphasewithnoise  = angle((x1.*conj(x2)));	
    wrappedphase = (exp(1i*absolutephase));
    
    realwithnoise = cos(wrappedphasewithnoise);
    imagwithnoise = sin(wrappedphasewithnoise);
    
    realnoiseless = real(wrappedphase);
    imagnoiseless = imag(wrappedphase);
    
    co100str = ['C',num2str(co*100),'_'];

    rwn_name = replace(file_name{i},'up',[co100str,'rwn']); % 'up' is old,'wpn' is new
    save_path_2 = ['D:\Datasets\Denoising\ALL\realwithnoise\',rwn_name];
    save (save_path_2,'realwithnoise');
    
    iwn_name = replace(file_name{i},'up',[co100str,'iwn']); % 'up' is old,'wpn' is new
    save_path_2 = ['D:\Datasets\Denoising\ALL\imagwithnoise\',iwn_name];
    save (save_path_2,'imagwithnoise');
    
    rnl_name = replace(file_name{i},'up',[co100str,'rnl']); % 'up' is old,'wpn' is new
    save_path_2 = ['D:\Datasets\Denoising\ALL\realnoiseless\',rnl_name];
    save (save_path_2,'realnoiseless');
    
    inl_name = replace(file_name{i},'up',[co100str,'inl']); % 'up' is old,'wpn' is new
    save_path_2 = ['D:\Datasets\Denoising\ALL\imagnoiseless\',inl_name];
    save (save_path_2,'imagnoiseless');
    
    slc1_name = replace(file_name{i},'up',[co100str,'slc1']); % 'up' is old,'wpn' is new
    save_path_2 = ['D:\Datasets\Denoising\ALL\slc1\',slc1_name];
    save (save_path_2,'slc1');
    
    slc2_name = replace(file_name{i},'up',[co100str,'slc2']); % 'up' is old,'wpn' is new
    save_path_2 = ['D:\Datasets\Denoising\ALL\slc2\',slc2_name];
    save (save_path_2,'slc2');

end



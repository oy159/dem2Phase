function save_patch_group(filename, group)
%SAVE_PATCH_GROUP Save flat MAT variables for MATLAB and scipy compatibility.
    parent_dir = fileparts(filename);
    if ~exist(parent_dir, 'dir')
        mkdir(parent_dir);
    end
    save(filename, '-struct', 'group', '-v7');
end

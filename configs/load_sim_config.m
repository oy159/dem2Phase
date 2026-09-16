function cfg = load_sim_config(config_path)
%LOAD_SIM_CONFIG Load and validate a dem2phase simulation JSON file.

    arguments
        config_path (1, :) char
    end

    if ~isfile(config_path)
        error('dem2phase:ConfigNotFound', ...
            'Simulation config does not exist: %s', config_path);
    end

    cfg = jsondecode(fileread(config_path));
    cfg = validate_sim_config(cfg);
    cfg.config_path = config_path;
end


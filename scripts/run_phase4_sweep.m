function run_phase4_sweep()
%RUN_PHASE4_SWEEP Generate a paired 4-scenario x 3-level Phase-4 pilot.
    script_dir = fileparts(mfilename('fullpath'));
    project_root = fileparts(script_dir);
    manifest = jsondecode(fileread(fullfile(project_root, 'configs', ...
        'phase4_sweep_matrix.json')));
    scenarios = cellstr(manifest.scenarios);
    levels = fieldnames(manifest.levels);
    generator_path = fullfile(script_dir, 'gen_dataset_from_dem_v13.m');

    old_scenario = getenv('DEM2PHASE_PHASE4_SCENARIO');
    old_level = getenv('DEM2PHASE_PHASE4_LEVEL');
    old_patches = getenv('DEM2PHASE_PHASE4_PATCHES');
    cleanup = onCleanup(@() restore_environment( ...
        old_scenario, old_level, old_patches));
    setenv('DEM2PHASE_PHASE4_PATCHES', num2str(manifest.patches_per_dem));

    for sweep_scenario_idx = 1:numel(scenarios)
        for sweep_level_idx = 1:numel(levels)
            setenv('DEM2PHASE_PHASE4_SCENARIO', scenarios{sweep_scenario_idx});
            setenv('DEM2PHASE_PHASE4_LEVEL', levels{sweep_level_idx});
            fprintf('\nPHASE-4 SWEEP: %s / %s\n', ...
                scenarios{sweep_scenario_idx}, levels{sweep_level_idx});
            run(generator_path);
        end
    end
    clear cleanup
end

function restore_environment(scenario, level, patches)
    setenv('DEM2PHASE_PHASE4_SCENARIO', scenario);
    setenv('DEM2PHASE_PHASE4_LEVEL', level);
    setenv('DEM2PHASE_PHASE4_PATCHES', patches);
end

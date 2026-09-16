%% gen_dataset_from_dem.m
%
% Demo: Generate an InSAR training dataset from DEM / DSM GeoTIFF files.
%
% Data layout expected:
%   <repo_root>/data/dem/          <- input files: *_dem.tif and *_DSM*.tif
%   <repo_root>/data/patches/      <- output root (sub-folders created automatically)
%
% Functions required (all in dem2phase/):
%   gen_coherence_map.m
%   patch_add_noise_spatialcoh.m
%   calc_coherence.m
%
% Outputs (in addition to phase patches):
%   data/patches/coherence_estimated/     <- calc_coherence result after noise
%   data/patches/generation_log_<ts>.csv    <- per-file summary log
%   data/patches/generation_detail_<ts>.csv <- per-patch detail log (with num stats)
%
% See docs/gen_dataset_from_dem_usage.md for the full usage guide.
%
% -------------------------------------------------------------------------
clc; clear; close all;

%% ══════════════════════════════════════════════════════════════════════════
%  USER CONFIGURATION  –  edit these values before running the script
%  ══════════════════════════════════════════════════════════════════════════

% --- Output volume --------------------------------------------------------
patches_per_dem     = 200;    % random patches to draw from each DEM/DSM file
patch_size          = 256;    % output patch dimensions [px]

% --- Scale pre-computation (key speed / RAM trade-off) --------------------
% interp2 is run ONCE per scale level per DEM (not once per attempt).
% n_precomp_scales evenly-spaced scale levels are pre-computed for each DEM,
% then each patch attempt picks one at random.
%
% Peak RAM per DEM (3600×3600 input, interp_scale_range = [2 4]):
%   n=1  →  ~0.9 GB   (single midpoint scale, fastest start-up)
%   n=3  →  ~3.0 GB   (three levels: min / mid / max)
%   n=5  →  ~5.0 GB   (five levels, most scale diversity)
n_precomp_scales    = 3;      % number of discrete scale levels per DEM

interp_scale_range  = [2, 4]; % [min, max] bicubic up-scaling factor

% --- Reproducibility ------------------------------------------------------
% Set rng_seed to a positive integer for fully reproducible output, or to
% 'shuffle' to seed from the system clock (different result every run).
% The actual seed used is always saved to rng_seed_<ts>.mat and printed.
rng_seed            = 42;     % positive integer, or 'shuffle'

% --- Patch quality filters ------------------------------------------------
max_wrap_count      = 40;     % discard patches with more wrap cycles than this
void_frac_thresh    = 0.05;   % skip if fraction of void pixels (num=0) exceeds this

% --- Output file types ----------------------------------------------------
% Disk usage: ~1.20 MB / patch with all 6 types + coherence_estimated.
% Comment out types you do not need to reduce storage.
%
%   All 6 types + coherence  →  ~1.20 MB / patch  →  ~600 MB / DEM  →  n × 0.6 GB total
%   4 essential types only   →  ~0.40 MB / patch  →  ~200 MB / DEM  →  n × 0.2 GB total
type_options = { ...
    'wrapped_phase', ...
    'wrapped_phase_real', ...
    'wrapped_phase_imag', ...
    'wrapped_phase_withnoise', ...   % also saves _real and _imag sub-files
    'unwrapped_phase', ...
    'wrap_count' };
% (coherence_estimated is always saved regardless)

%% ══════════════════════════════════════════════════════════════════════════
%  END USER CONFIGURATION
%  ══════════════════════════════════════════════════════════════════════════

%% ── 0.  Paths ─────────────────────────────────────────────────────────────
script_dir    = fileparts(mfilename('fullpath'));   % dem2phase/scripts/
dem2phase_dir = fileparts(script_dir);             % dem2phase/
repo_root     = fileparts(dem2phase_dir);          % repo root

dem_dir     = fullfile(repo_root, 'data', 'dem');
output_root = fullfile(repo_root, 'data', 'patches');

addpath(genpath(dem2phase_dir));  % adds dem2phase/ and all subdirs (noise/, coherence/, …)

%% ── 1.  Discover DEM / DSM files ─────────────────────────────────────────
% Accept any .tif whose name contains "_dem" OR "_DSM" (case-insensitive).
all_tif   = dir(fullfile(dem_dir, '*.tif'));
is_dem    = contains({all_tif.name}, '_dem', 'IgnoreCase', true);
is_dsm    = contains({all_tif.name}, '_DSM', 'IgnoreCase', true);
dem_files = all_tif(is_dem | is_dsm);

fprintf('\n========================================\n');
fprintf(' DEM/DSM directory : %s\n', dem_dir);
fprintf(' Files found       : %d\n', numel(dem_files));
fprintf('========================================\n');

if numel(dem_files) == 0
    fprintf('[WARN] No *_dem*.tif or *_DSM*.tif files found in:\n  %s\n', dem_dir);
    fprintf('       Place GeoTIFF files there and re-run.\n');
    return
end

for k = 1:numel(dem_files)
    meta = parse_dem_metadata(dem_files(k).name);
    fprintf('  [%02d] %-45s  sensor=%-8s  tile=%-12s  strip=%-4s  type=%s\n', ...
        k, dem_files(k).name, meta.sensor, meta.tile, meta.strip, meta.file_type);
end
fprintf('\n');

%% ── 2.  Parameters ────────────────────────────────────────────────────────
% (patch_size, patches_per_dem, interp_scale_range, n_precomp_scales,
%  max_wrap_count, void_frac_thresh, and type_options are set at the top.)

% InSAR system parameters – one row per system (multi-baseline diversity).
% Systems: [1] ALOS-2 L-band  [2] TerraSAR-X X-band  [3] Ka-band
% dem2phase_ratio is drawn per patch by randomly selecting one system.
insar_lambda      = [0.236,   0.031,   0.005  ];  % wavelength  [m]
insar_baseline    = [230,     237,     16     ];   % perpendicular baseline [m]
insar_slant_range = [868142,  640708,  845521 ];   % slant range [m]
insar_theta_deg   = [38.77,   41.07,   39.67  ];   % incidence angle [deg]

% Coherence map generation settings (region mode, full dynamic range)
coh_min            = 0.01;   % near-zero coherence (e.g. water surfaces)
coh_max            = 0.95;
coh_spatial_scale  = 30;     % boundary smoothing length [px]
coh_n_regions_mean = 7;      % Gaussian mean for random n_regions (range 1–20)
coh_n_regions_std  = 4;      % Gaussian std  for random n_regions

% Coherence estimation window (for calc_coherence after adding noise)
coh_est_win = 7;

%% ── 3.  Create output sub-directories ────────────────────────────────────
subdir_map = struct( ...
    'wrapped_phase',                'wrappedphase_withoutnoise', ...
    'wrapped_phase_real',           'wrappedphase_withoutnoise_real', ...
    'wrapped_phase_imag',           'wrappedphase_withoutnoise_imag', ...
    'wrapped_phase_withnoise',      'wrappedphase_withnoise', ...
    'wrapped_phase_withnoise_real', 'wrappedphase_withnoise_real', ...
    'wrapped_phase_withnoise_imag', 'wrappedphase_withnoise_imag', ...
    'unwrapped_phase',              'unwrapped_phase', ...
    'wrap_count',                   'wrapped_number', ...
    'coherence_estimated',          'coherence_estimated' );

all_subdirs = fieldnames(subdir_map);
for s = 1:numel(all_subdirs)
    d = fullfile(output_root, subdir_map.(all_subdirs{s}));
    if ~exist(d, 'dir'); mkdir(d); end
end

%% ── 4.  Open log files (summary + detail) ────────────────────────────────
log_ts = datestr(now, 'yyyymmdd_HHMMSS');

% Summary log: one row per source DEM/DSM file
log_path = fullfile(output_root, sprintf('generation_log_%s.csv', log_ts));
log_fid  = fopen(log_path, 'w');
fprintf(log_fid, ['fi,source_file,file_type,sensor,tile,strip,' ...
    'dem_rows,dem_cols,dem_size_bytes,' ...
    'num_available,num_file,' ...
    'patches_saved,attempts_total,' ...
    'void_skipped,wrap_skipped,' ...
    'num_void_frac_mean,elapsed_s,rng_seed,timestamp\n']);

% Detail log: one row per saved patch
detail_path = fullfile(output_root, sprintf('generation_detail_%s.csv', log_ts));
detail_fid  = fopen(detail_path, 'w');
fprintf(detail_fid, ['patch_global_id,patch_name,source_file,file_type,sensor,' ...
    'tile,strip,dem_rows,dem_cols,dem_size_bytes,' ...
    'interp_scale,scaled_rows,scaled_cols,' ...
    'crop_r0,crop_c0,crop_r1,crop_c1,' ...
    'max_wrap_count,coh_input_mean,coh_input_min,coh_input_max,' ...
    'coh_estimated_mean,coh_estimated_min,coh_estimated_max,' ...
    'dem2phase_ratio,lambda_m,baseline_m,slant_range_m,theta_deg,' ...
    'coh_n_regions,insar_sys_idx,' ...
    'num_available,num_void_frac,num_mean,num_min,' ...
    'patch_size,timestamp\n']);

fprintf('Summary log : %s\n', log_path);
fprintf('Detail  log : %s\n\n', detail_path);

%% ── 5.  Main loop: DEM/DSM → patches ────────────────────────────────────
total_saved     = 0;
total_skip      = 0;
total_void_skip = 0;

% Set random seed and record the actual seed used for reproducibility.
rng(rng_seed);
rng_state      = rng;          % capture current state (includes actual seed)
rng_seed_used  = rng_state.Seed;
rng_seed_path  = fullfile(output_root, sprintf('rng_seed_%s.mat', log_ts));
save(rng_seed_path, 'rng_state', 'rng_seed_used');
fprintf('RNG seed     : %d  (saved to %s)\n\n', rng_seed_used, rng_seed_path);

tic
for fi = 1:numel(dem_files)
    dem_path  = fullfile(dem_dir, dem_files(fi).name);
    [~, dem_stem, ~] = fileparts(dem_files(fi).name);
    meta = parse_dem_metadata(dem_files(fi).name);

    fprintf('──────────────────────────────────────\n');
    fprintf('[%02d/%02d] %s\n', fi, numel(dem_files), dem_files(fi).name);
    fprintf('         sensor=%s  tile=%s  strip=%s  type=%s\n', ...
        meta.sensor, meta.tile, meta.strip, meta.file_type);

    % ── 5a. Read file ─────────────────────────────────────────────────────
    dem_raw = double(imread(dem_path));
    [dem_rows, dem_cols] = size(dem_raw);
    fprintf('         Size    : %d x %d px  (%s)\n', dem_rows, dem_cols, ...
        format_bytes(dem_files(fi).bytes));

    % ── 5a'. Load NUM file (ASTGTM stacking-count, optional) ──────────────
    num_map      = [];
    num_file     = '';
    num_avail    = 0;
    num_path_try = get_num_path(dem_path);
    if ~isempty(num_path_try)
        num_map   = uint8(imread(num_path_try));
        num_file  = num_path_try;
        num_avail = 1;
        fprintf('         NUM     : %s\n', num_path_try);
    end

    % ── 5b. Random patch extraction ───────────────────────────────────────
    patch_count         = 0;
    attempt             = 0;
    max_attempts        = patches_per_dem * 5;
    void_skipped        = 0;
    wrap_skipped        = 0;
    num_void_frac_sum   = 0;
    num_void_frac_count = 0;
    file_tic            = tic;

    % ── 5b'. Pre-compute scaled DEM images (raw elevation; ratio applied per patch)
    % interp2 is called once per scale level here, not once per attempt.
    % Each patch attempt picks one of these images at random, giving scale
    % diversity while reducing interp2 calls from ~max_attempts down to
    % n_precomp_scales (typically 3).
    precomp_sc   = linspace(interp_scale_range(1), interp_scale_range(2), n_precomp_scales);
    precomp_imgs = cell(n_precomp_scales, 1);
    [xg_base, yg_base] = meshgrid(1:dem_cols, 1:dem_rows);
    for si = 1:n_precomp_scales
        sc = precomp_sc(si);
        rc = round(sc * dem_rows);
        cc = round(sc * dem_cols);
        if rc < patch_size + 1 || cc < patch_size + 1
            precomp_imgs{si} = dem_raw;
            precomp_sc(si)   = 1.0;
        else
            xlin = linspace(1, dem_cols, cc);
            ylin = linspace(1, dem_rows, rc);
            [xxg, yyg] = meshgrid(xlin, ylin);
            precomp_imgs{si} = interp2(xg_base, yg_base, dem_raw, xxg, yyg, 'bicubic');
        end
    end
    fprintf('         Scaled  : %d pre-computed images at scales [%s]\n', ...
        n_precomp_scales, num2str(precomp_sc, '%.2f '));

    while patch_count < patches_per_dem && attempt < max_attempts
        attempt = attempt + 1;

        % Pick one of the pre-computed scaled images at random
        si          = randi(n_precomp_scales);
        scale       = precomp_sc(si);
        dem_interp  = precomp_imgs{si};
        [rows_scaled, cols_scaled] = size(dem_interp);

        % Random crop
        [mr, mc] = size(dem_interp);
        r0 = randi([1, max(1, mr - patch_size)]);
        c0 = randi([1, max(1, mc - patch_size)]);
        r1 = r0 + patch_size - 1;
        c1 = c0 + patch_size - 1;

        if r1 > mr || c1 > mc
            continue
        end

        patch_dem = dem_interp(r0:r1, c0:c1);

        if size(patch_dem,1) ~= patch_size || size(patch_dem,2) ~= patch_size
            continue
        end

        % Normalise DC offset
        patch_dem = patch_dem - min(patch_dem(:));

        % ── Per-patch InSAR system (multi-baseline diversity) ─────────────
        sys_idx         = randi(numel(insar_lambda));
        lambda          = insar_lambda(sys_idx);
        baseline        = insar_baseline(sys_idx);
        slant_range     = insar_slant_range(sys_idx);
        theta_deg       = insar_theta_deg(sys_idx);
        dem2phase_ratio = 4*pi / lambda * baseline / slant_range / sind(theta_deg);
        patch_dem       = patch_dem * dem2phase_ratio;

        % Wrap count filter
        wrap_count_map = round((patch_dem - angle(exp(1i * patch_dem))) / (2*pi));
        if max(wrap_count_map(:)) > max_wrap_count
            wrap_skipped = wrap_skipped + 1;
            total_skip   = total_skip + 1;
            continue
        end

        % NUM quality filter (only when _num.tif was found alongside the DEM)
        num_void_frac = NaN;
        num_mean_val  = NaN;
        num_min_val   = NaN;
        if ~isempty(num_map)
            % Map crop coordinates back to original DEM space using the same
            % linspace grid that interp2 used, so the mapping is exact.
            scale_r = (dem_rows - 1) / max(rows_scaled - 1, 1);
            scale_c = (dem_cols - 1) / max(cols_scaled - 1, 1);
            nr0 = max(1,        round(1 + (r0 - 1) * scale_r));
            nr1 = min(dem_rows, round(1 + (r1 - 1) * scale_r));
            nc0 = max(1,        round(1 + (c0 - 1) * scale_c));
            nc1 = min(dem_cols, round(1 + (c1 - 1) * scale_c));
            if nr1 >= nr0 && nc1 >= nc0
                num_crop      = double(num_map(nr0:nr1, nc0:nc1));
                num_void_frac = mean(num_crop(:) == 0);
                num_mean_val  = mean(num_crop(:));
                num_min_val   = min(num_crop(:));
                if num_void_frac > void_frac_thresh
                    void_skipped        = void_skipped + 1;
                    total_skip          = total_skip + 1;
                    total_void_skip     = total_void_skip + 1;
                    continue
                end
                num_void_frac_sum   = num_void_frac_sum   + num_void_frac;
                num_void_frac_count = num_void_frac_count + 1;
            end
        end

        % ── 5c. Phase representations ─────────────────────────────────────
        patch_count = patch_count + 1;
        total_saved = total_saved + 1;
        global_id   = total_saved;

        patch_name = sprintf('%s_patch_%05d', dem_stem, global_id);

        wrapped_phase   = angle(exp(1i * patch_dem));
        wrapped_phase_r = real(exp(1i * patch_dem));
        wrapped_phase_i = imag(exp(1i * patch_dem));

        % Input coherence map – region mode with random number of Voronoi regions
        % n_reg sampled from truncated Gaussian in [1, 20]
        n_reg = max(1, min(20, round(coh_n_regions_mean + coh_n_regions_std * randn())));
        coh_map_input = gen_coherence_map(patch_size, patch_size, ...
            'min_coh', coh_min, 'max_coh', coh_max, ...
            'spatial_scale', coh_spatial_scale, ...
            'mode', 'region', 'n_regions', n_reg);

        % Add spatially-varying noise
        noisy_phase   = patch_add_noise_spatialcoh(wrapped_phase, coh_map_input);
        noisy_phase_r = real(exp(1i * noisy_phase));
        noisy_phase_i = imag(exp(1i * noisy_phase));

        % Estimated coherence map (re-computed from noisy vs. clean phase)
        [coh_map_est, coh_est_mean] = calc_coherence(noisy_phase, wrapped_phase, ...
            'method', 'phase', 'win_size', coh_est_win);

        % ── 5d. Save patch types ──────────────────────────────────────────
        for ti = 1:numel(type_options)
            ptype = type_options{ti};
            switch ptype
                case 'wrapped_phase'
                    save_mat(output_root, subdir_map.wrapped_phase, patch_name, ...
                        ptype, wrapped_phase);
                case 'wrapped_phase_real'
                    save_mat(output_root, subdir_map.wrapped_phase_real, patch_name, ...
                        ptype, wrapped_phase_r);
                case 'wrapped_phase_imag'
                    save_mat(output_root, subdir_map.wrapped_phase_imag, patch_name, ...
                        ptype, wrapped_phase_i);
                case 'wrapped_phase_withnoise'
                    save_mat(output_root, subdir_map.wrapped_phase_withnoise, patch_name, ...
                        ptype, noisy_phase);
                    save_mat(output_root, subdir_map.wrapped_phase_withnoise_real, patch_name, ...
                        'wrapped_phase_withnoise_real', noisy_phase_r);
                    save_mat(output_root, subdir_map.wrapped_phase_withnoise_imag, patch_name, ...
                        'wrapped_phase_withnoise_imag', noisy_phase_i);
                case 'unwrapped_phase'
                    save_mat(output_root, subdir_map.unwrapped_phase, patch_name, ...
                        ptype, patch_dem);
                case 'wrap_count'
                    save_mat(output_root, subdir_map.wrap_count, patch_name, ...
                        ptype, wrap_count_map);
            end
        end

        % ── 5e. Save coherence map ───────────────────────────────────────
        save_mat(output_root, subdir_map.coherence_estimated, patch_name, ...
            'coherence_estimated', coh_map_est);

        % ── 5f. Write detail log entry ────────────────────────────────────
        fprintf(detail_fid, '%d,%s,%s,%s,%s,%s,%s,%d,%d,%d,%.4f,%d,%d,%d,%d,%d,%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.6f,%.4f,%.1f,%.1f,%.4f,%d,%d,%d,%.4f,%.4f,%.4f,%d,%s\n', ...
            global_id, patch_name, dem_files(fi).name, meta.file_type, meta.sensor, ...
            meta.tile, meta.strip, ...
            dem_rows, dem_cols, dem_files(fi).bytes, ...
            scale, rows_scaled, cols_scaled, ...
            r0, c0, r1, c1, ...
            max(wrap_count_map(:)), ...
            mean(coh_map_input(:)), min(coh_map_input(:)), max(coh_map_input(:)), ...
            coh_est_mean, min(coh_map_est(:)), max(coh_map_est(:)), ...
            dem2phase_ratio, lambda, baseline, slant_range, theta_deg, ...
            n_reg, sys_idx, ...
            num_avail, num_void_frac, num_mean_val, num_min_val, ...
            patch_size, datestr(now, 'yyyy-mm-dd HH:MM:SS'));

        if mod(patch_count, 50) == 0
            fprintf('         Saved %d / %d  (skipped: %d)\n', ...
                patch_count, patches_per_dem, total_skip);
        end
    end % while

    fprintf('         Done : %d patches  (skipped: %d wrap, %d void, attempts: %d)\n', ...
        patch_count, wrap_skipped, void_skipped, attempt);

    % Write summary log row (one per source file)
    file_elapsed = toc(file_tic);
    if num_void_frac_count > 0
        num_void_frac_mean = num_void_frac_sum / num_void_frac_count;
    else
        num_void_frac_mean = NaN;
    end
    fprintf(log_fid, '%d,%s,%s,%s,%s,%s,%d,%d,%d,%d,%s,%d,%d,%d,%d,%.4f,%.1f,%d,%s\n', ...
        fi, dem_files(fi).name, meta.file_type, meta.sensor, meta.tile, meta.strip, ...
        dem_rows, dem_cols, dem_files(fi).bytes, ...
        num_avail, num_file, ...
        patch_count, attempt, void_skipped, wrap_skipped, ...
        num_void_frac_mean, file_elapsed, rng_seed_used, ...
        datestr(now, 'yyyy-mm-dd HH:MM:SS'));
end % for fi

fclose(log_fid);
fclose(detail_fid);

elapsed = toc;
fprintf('\n========================================\n');
fprintf(' Total patches saved : %d\n', total_saved);
fprintf(' Skipped (all)       : %d\n', total_skip);
fprintf('   wrap (>%d)        : %d\n', max_wrap_count, total_skip - total_void_skip);
fprintf('   void-fill         : %d\n', total_void_skip);
fprintf(' Elapsed time        : %.1f s\n', elapsed);
fprintf(' RNG seed used       : %d\n', rng_seed_used);
fprintf(' Output directory    : %s\n', output_root);
fprintf(' Summary log         : %s\n', log_path);
fprintf(' Detail  log         : %s\n', detail_path);
fprintf(' RNG seed file       : %s\n', rng_seed_path);
fprintf('========================================\n');

%% ── Helper functions ──────────────────────────────────────────────────────

function save_mat(root, subdir, patch_name, patch_type, data)
%SAVE_MAT  Save data to a .mat file named <patch_name>_<patch_type>.mat
    filename = sprintf('%s_%s.mat', patch_name, patch_type);
    save(fullfile(root, subdir, filename), 'data');
end

% -------------------------------------------------------------------------
function meta = parse_dem_metadata(filename)
%PARSE_DEM_METADATA  Extract sensor, tile, strip and file-type from a filename.
%
%   Supported naming patterns (case-insensitive):
%     ASTGTM_<TILE>_dem.tif           -> sensor=ASTGTM, tile=<TILE>, strip=N/A
%     ASTGTMV003_<TILE>_dem.tif       -> sensor=ASTGTM (version tag ignored)
%     ALPSMLC30_<TILE>_DSM.tif        -> sensor=ALPSMLC (ALOS AW3D30)
%     SRTM_<TILE>_dem.tif             -> sensor=SRTM
%     ALOS_<TILE>_dem.tif             -> sensor=ALOS
%     N<lat>E<lon>_dem.tif            -> sensor=unknown, tile=N<lat>E<lon>
%     *_DSM_<IDX>*.tif                -> file_type=DSM, strip=<IDX>
%     system_<S>_DSM_<IDX>*.tif       -> sensor=system_<S>, strip=<IDX>
%     <anything>                       -> graceful fallback

    meta.sensor    = 'unknown';
    meta.tile      = 'N/A';
    meta.strip     = 'N/A';
    meta.file_type = 'DEM';

    [~, stem, ~] = fileparts(filename);
    stem_up = upper(stem);

    % ── File type ──────────────────────────────────────────────────────────
    if contains(stem_up, '_DSM')
        meta.file_type = 'DSM';
    else
        meta.file_type = 'DEM';
    end

    % ── Sensor ─────────────────────────────────────────────────────────────
    % ALPSMLC is the ALOS World 3D (AW3D30) product; listed before ALOS so
    % the more-specific name matches first.
    sensors = {'ASTGTM', 'SRTM', 'ALPSMLC', 'ALOS2', 'ALOS', 'COPDEM', 'TDX', 'TANDEM'};
    for s = 1:numel(sensors)
        if contains(stem_up, upper(sensors{s}))
            meta.sensor = sensors{s};
            break
        end
    end

    % system_<N>_DSM pattern from original ALOS_Batch_Read.m
    tok = regexp(stem, 'system_(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        meta.sensor = sprintf('system_%s', tok{1});
    end

    % ── Tile (lat/lon grid cell) ───────────────────────────────────────────
    % Matches 2- or 3-digit latitude (e.g. N30, N027) and 3-digit longitude
    % (e.g. E103, E087), covering both SRTM/ASTER and ALOS AW3D30 styles.
    tok = regexp(stem, '[NnSs]\d{2,3}[EeWw]\d{3}', 'match', 'once');
    if ~isempty(tok)
        meta.tile = upper(tok);
    end

    % ── Strip / track index ───────────────────────────────────────────────
    % DSM_<idx> or Track<idx> or _T<idx>_ or a trailing number.
    tok = regexp(stem, '(?:DSM|dsm|track|Track|_T)_?(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        meta.strip = tok{1};
    else
        % Fallback: last standalone number in the stem, excluding any digits
        % that belong to the lat/lon tile token (e.g. N29E102) to avoid
        % extracting coordinate sub-numbers like "02" from "E102".
        stem_no_tile = regexprep(stem, '[NnSs]\d+[EeWw]\d+', '');
        toks = regexp(stem_no_tile, '(?<![A-Za-z\d])(\d+)(?![A-Za-z\d])', 'tokens');
        if ~isempty(toks)
            meta.strip = toks{end}{1};
        end
    end
end

% -------------------------------------------------------------------------
function str = format_bytes(n)
%FORMAT_BYTES  Human-readable file size string.
    if n >= 1e6
        str = sprintf('%.1f MB', n/1e6);
    elseif n >= 1e3
        str = sprintf('%.1f KB', n/1e3);
    else
        str = sprintf('%d B', n);
    end
end

% -------------------------------------------------------------------------
function num_path = get_num_path(dem_path)
%GET_NUM_PATH  Return path of the companion _num.tif, or '' if absent.
%
%   ASTGTM tiles follow the convention:
%       <prefix>_dem.tif  ->  <prefix>_num.tif
%   For files without a _dem suffix the function tries appending _num
%   before the extension as a fallback.  Returns '' when no file is found,
%   so callers can safely skip NUM-based filtering.
    [d, stem, ext] = fileparts(dem_path);
    % Replace trailing _dem (case-insensitive) with _num
    num_stem = regexprep(stem, '_[Dd][Ee][Mm]$', '_num');
    if strcmp(num_stem, stem)
        % No _dem suffix found – try appending _num as fallback
        num_stem = [stem '_num'];
    end
    candidate = fullfile(d, [num_stem ext]);
    if exist(candidate, 'file')
        num_path = candidate;
    else
        num_path = '';
    end
end

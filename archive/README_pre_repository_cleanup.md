# dem2phase

MATLAB toolbox for generating synthetic InSAR training datasets from DEM / DSM GeoTIFF files.

The active generator also has a Linux/cloud Python implementation under
[`python/`](python/README.md). It preserves the grouped MAT schema and current
P-band multi-UAV physics, supports deterministic multi-process generation and
resume, and includes MATLAB/Python statistical comparison tools.

The current mainline is a focused, patch-level simulator for distributed UAV
InSAR. Its nominal profile uses 500 MHz P band, 500 m altitude, a five-UAV star
graph and four configurable effective perpendicular baselines
`[0.75, 1.20, 1.80, 3.00] m`. It includes terrain/land-cover-aware coherence,
shared-node complex SLC noise, multilooking, synchronization/trajectory/coregistration
errors, edge dropout and mask-aware failure labels. It is not a raw-echo SAR
simulator.

The complete Chinese technical reference, including equations, parameter tables,
MAT schema and known limitations, is in
[`docs/current_code_technical_reference_zh.md`](docs/current_code_technical_reference_zh.md).

---

## Directory layout

```
dem2phase/
├── scripts/        Runnable user scripts
│   ├── gen_dataset_from_dem.m        Main dataset generation script
│   ├── validate_dataset.m            QA / validation script
│   └── visualize_patches.m           Patch visualisation script
│
├── coherence/      Coherence map generation and estimation
│   ├── gen_coherence_map.m           Single-baseline spatially-varying coherence map
│   ├── gen_coherence_map_multibaseline.m  Multi-baseline coherence stack
│   └── calc_coherence.m              Coherence estimation from SLC pair or phase
│
├── noise/          Phase noise injection functions
│   ├── patch_add_noise_spatialcoh.m  Noise with spatially-varying coherence (main)
│   ├── patch_add_noise_multibaseline.m  Noise for a multi-baseline phase stack
│   └── patch_add_noise.m             Legacy uniform-coherence noise (scalar γ)
│
└── legacy/         Old scripts kept for reference (hardcoded paths, not maintained)
    ├── ALOS_Batch_Read.m
    ├── addnoise.m
    └── get_patch.m
```

---

## Quick start

### 1 — Place DEM / DSM files

The default UAV profile scans the project directory itself, so the bundled
`ALPSMLC30_N027E085_DSM.tif` is immediately usable. Change
`dataset.dem_directory` in `configs/uav_p_500m_monostatic.json` to point at a
different DEM collection. Relative paths are resolved from `dem2phase/`.

```
dem2phase/
```

### 2 — Generate a dataset

Open MATLAB, navigate to any directory, then run:

```matlab
run('C:/Users/oycr/Documents/research/dem2phase/scripts/gen_dataset_from_dem_v13.m')
```

The script loads `configs/uav_p_500m_monostatic.json`, discovers matching
GeoTIFFs in the configured directory, draws
terrain-weighted patches, adds terrain- and baseline-dependent phase noise, and
saves the dataset below `dem2phase/data/` by default.

Key dataset parameters are configured in the JSON profile:

| Parameter | Default | Description |
|---|---|---|
| `dataset.generation.patches_per_dem` | 32 | Accepted patches per DEM file |
| `dataset.generation.patch_size` | 256 | Output patch size in pixels |
| `dataset.generation.max_wrap_count` | 40 | Reject patches with more wrap cycles |
| `dataset.generation.void_fraction_threshold` | 0.05 | Reject patches with more than 5% void pixels |
| `dataset.generation.rng_seed` | 42 | Reproducible batch seed |
| `dataset.generation.n_precomp_scales` | 3 | Pre-computed bicubic scale levels per DEM |
| `dataset.generation.interp_scale_range` | `[2, 4]` | Min/max upscaling factor |

The fixed batch seed reproduces the dataset, while candidate- and module-specific
counter-derived seeds keep random realizations different within the same batch.

### Fixed whole-DEM splits

`configs/dem_split_manifest.csv` assigns every geographic 1-degree tile to
exactly one of `train`, `test`, or `validation`. Random cropping is allowed
within a DEM, but all scales, overlapping crops, baselines, noise realizations,
and cross-sensor versions of the same tile inherit the same split. The
generator writes directly to
`patch_groups/train`, `patch_groups/test`, and `patch_groups/validation`.

The current corpus contains 19 disjoint geographic tiles and 234 patches:
12 training tiles (162 patches), five Himalayan test tiles (56 patches), and
two small out-of-region validation tiles from the Alps and Andes (16 patches).
Copernicus GLO-30 is the primary elevation source, ALOS AW3D30 is secondary
geographic augmentation, and four ASTER GDEM V3 tiles with NUM companions are
low-weight robustness data. The manifest loader rejects both repeated files
and cross-sensor geographic leakage. Balancing weights use training data only.

### 3 — Validate the output

```matlab
run('<repo_root>/dem2phase/scripts/validate_dataset.m')
```

Checks file counts, value ranges, real/imag consistency, and prints a QA report.

### 4 — Visualise patches

```matlab
% Edit idx_list inside the script first, then:
run('<repo_root>/dem2phase/scripts/visualize_patches.m')
```

Opens one figure per requested index showing all available data types as subplots.

---

## Data layout

The current V13 generator writes one grouped MAT per patch. The main training
observations are `wrappedphase_withnoise`, `coherence_observed`,
`valid_edge_mask` and `coregistration_valid_mask`; labels include
`unwrapped_phase`, `coherence_true`, graph geometry and structured failure
annotations. Use the strict loader so inactive edges and invalid resampling pixels
cannot enter a loss:

```matlab
sample = load_patch_group_for_training('path/to/patch.mat');
```

`coherence_observed` is estimated only from the two disturbed node SLCs. The
loader deliberately rejects obsolete samples whose coherence used clean phase.
`build_failure_balanced_manifest` creates inverse-frequency weights over dropout,
synchronization anomaly, low-coherence state, active edge count and terrain type.

The directory description below documents the legacy per-field format.

```
<repo_root>/data/
├── dem/                               Input GeoTIFF files (*_dem.tif, *_DSM*.tif)
└── patches/
    ├── wrappedphase_withoutnoise/     Clean wrapped phase     [-π, π]
    ├── wrappedphase_withoutnoise_real/  cos(phase)            [-1, 1]
    ├── wrappedphase_withoutnoise_imag/  sin(phase)            [-1, 1]
    ├── wrappedphase_withnoise/        Noisy wrapped phase     [-π, π]
    ├── wrappedphase_withnoise_real/   cos(noisy phase)        [-1, 1]
    ├── wrappedphase_withnoise_imag/   sin(noisy phase)        [-1, 1]
    ├── unwrapped_phase/               DEM-derived phase (unwrapped)
    ├── wrapped_number/                Integer wrap count (≥ 0)
    ├── coherence_estimated/           Estimated coherence     [0, 1]
    ├── generation_log_<ts>.csv        Per-file summary log
    └── generation_detail_<ts>.csv     Per-patch detail log
```

Each `.mat` file stores a single `[patch_size × patch_size]` `double` array in a variable named `data`.

For `scripts/gen_dataset_from_dem_v13.m`, logs now include explicit multi-baseline wrap statistics:
- `generation_log_*.csv` includes short-baseline summary fields:
  `short_baseline_id`, `short_baseline_length_m`, `short_wrap_mean`, `short_wrap_max`.
- `generation_detail_*.csv` includes per-patch multi-baseline fields:
  `max_wrap_count_all_baselines`, short-baseline wrap stats, and
  `max_wrap_count_baseline1..K` (one column per baseline).
- `generation_file_detail_*.csv` includes one row per saved baseline sample:
  file-level names (`*_wrappedphase_withnoise.mat`, `*_wrappedphase_withoutnoise.mat`,
  `*_unwrapped_phase.mat`, `*_coherence_estimated.mat`) with
  `baseline_idx`, `baseline_length_m`, `max_wrap_count`, `mean_wrap_count`.

`scripts/analyze_wrap_count_distribution.m` supports baseline-index split analysis via
`baseline_idx_filter` (`[]` means all baselines), and exports
`analyze_wrapcount_detail_*.csv` with per-file baseline wrap statistics.

---

## Function reference

### `scripts/gen_dataset_from_dem.m`

Main runnable script. Reads DEM GeoTIFFs, pre-computes bicubic-upscaled versions at several scale levels, draws random `patch_size × patch_size` crops, filters by wrap count and void fraction, adds spatially-varying InSAR phase noise via `patch_add_noise_spatialcoh`, estimates coherence via `calc_coherence`, and saves all patch types to disk.

### `scripts/validate_dataset.m`

QA script. Checks directory completeness, patch count balance across types, per-patch value ranges, real/imag consistency (`cos²+sin²=1`), and log file statistics.

### `scripts/visualize_patches.m`

Visualisation script. Sorts patches by their numeric ID suffix, then opens one `imagesc` figure per requested index with up to 9 subplots (one per data type).

---

### `coherence/gen_coherence_map.m`

Generate a spatially-varying coherence map for a single baseline.

```matlab
coh_map = gen_coherence_map(rows, cols)
coh_map = gen_coherence_map(rows, cols, 'min_coh', 0.3, 'max_coh', 0.9, ...
              'spatial_scale', 30, 'mode', 'smooth', 'seed', 1)
```

| Mode | Description |
|---|---|
| `'smooth'` | Gaussian-filtered white noise → gradual coherence transitions |
| `'region'` | Voronoi regions with lightly smoothed boundaries (land-cover types) |

### `coherence/gen_coherence_map_multibaseline.m`

Generate a `{1×K}` cell array of coherence maps for a multi-baseline InSAR stack. All maps share the same spatial pattern; mean coherence decreases exponentially with baseline length.

```matlab
coh_maps = gen_coherence_map_multibaseline(256, 256, [100, 300, 600], ...
               'max_coh', 0.9, 'decay_alpha', 0.5, 'seed', 1)
```

The upper-coherence bound for baseline `B` is:

```
max_k = min_coh + (max_coh - min_coh) * exp(-decay_alpha * (B - B_ref) / B_ref)
```

### `coherence/calc_coherence.m`

Estimate coherence from a complex SLC pair (`method = 'slc'`, standard sliding-window estimator) or from a noisy/clean phase pair (`method = 'phase'`, variance inversion through a pre-computed LUT).

```matlab
[coh_map, coh_mean] = calc_coherence(slc1, slc2)
[coh_map, coh_mean] = calc_coherence(noisy_phase, clean_phase, 'method', 'phase', 'win_size', 11)
```

---

### `noise/patch_add_noise_spatialcoh.m`

Add InSAR phase noise with a **spatially-varying** coherence map. The per-pixel noise PDF is the standard InSAR phase noise distribution:

```
p(φ; γ) = (1−γ²)/(2π) · 1/(1−γ²cos²φ) · [1 + γcosφ·acos(−γcosφ)/√(1−γ²cos²φ)]
```

Coherence is quantised into `n_bins` levels; noise is sampled via inverse-CDF for each bin.

```matlab
noisy_phase = patch_add_noise_spatialcoh(clean_phase, coh_map)
noisy_phase = patch_add_noise_spatialcoh(clean_phase, coh_map, 'seed', 42, 'n_bins', 50)
```

Passing a scalar `coh_map` is equivalent to the legacy `patch_add_noise.m`.

### `noise/patch_add_noise_multibaseline.m`

Apply `patch_add_noise_spatialcoh` across a baseline stack with **independent per-baseline RNG seeds** (`base_seed + k`), modelling statistically independent thermal noise between acquisitions.

```matlab
noisy = patch_add_noise_multibaseline(clean_phs, coh_maps)
noisy = patch_add_noise_multibaseline(clean_phs, coh_maps, 'base_seed', 42)
```

`clean_phs` and `coh_maps` are each a `{1×K}` cell array or a single array broadcast to all `K` baselines.

### `noise/patch_add_noise.m`  *(legacy)*

Original noise injection function using a scalar coherence value and a saved RNG state file (`random_generate_state.mat`). Superseded by `patch_add_noise_spatialcoh`.

---

## Typical single-baseline workflow

```matlab
addpath(genpath('dem2phase'));

% 1. Load a DEM patch and convert to unwrapped phase
dem_patch = double(imread('tile.tif'));
ratio     = 4*pi/0.056 * 150 / 800000 / sind(40);   % example parameters
unw_phase = dem_patch * ratio;

% 2. Clean wrapped phase
clean_phase = angle(exp(1i * unw_phase));

% 3. Generate spatially-varying coherence map
coh_map = gen_coherence_map(256, 256, 'min_coh', 0.4, 'max_coh', 0.9, 'seed', 1);

% 4. Add phase noise
noisy_phase = patch_add_noise_spatialcoh(clean_phase, coh_map, 'seed', 42);

% 5. Estimate coherence from noisy / clean pair
[coh_est, ~] = calc_coherence(noisy_phase, clean_phase, 'method', 'phase', 'win_size', 11);
```

## Typical multi-baseline workflow

```matlab
baselines = [100, 300, 600];

% Per-baseline coherence maps (shared spatial pattern, decreasing mean)
coh_maps = gen_coherence_map_multibaseline(256, 256, baselines, ...
               'max_coh', 0.9, 'decay_alpha', 0.5, 'seed', 1);

% Per-baseline clean wrapped phases
lambda = 0.056; slant_r = 800000; theta = 40;
clean_phs = cell(1, numel(baselines));
for k = 1:numel(baselines)
    ratio         = 4*pi/lambda * baselines(k) / slant_r / sind(theta);
    clean_phs{k}  = angle(exp(1i * dem_patch * ratio));
end

% Independent noise per baseline
noisy = patch_add_noise_multibaseline(clean_phs, coh_maps, 'base_seed', 42);
```

## Configurable multi-UAV profile

`scripts/gen_dataset_from_dem_v13.m` now loads its radar and interferometric
geometry from `configs/uav_p_500m_monostatic.json`. The default profile uses:

- 500 MHz center frequency (`lambda ~= 0.60 m`)
- 500 m altitude, 45 degree incidence, and 707.1 m center slant range
- five UAVs and four star-graph interferometric edges
- perpendicular baselines `[0.75, 1.20, 1.80, 3.00] m`
- monostatic image pairs (`phase_path_multiplicity = 2`)
- ambiguity heights approximately `[200, 125, 83.3, 50] m`

Change the JSON profile rather than editing physics constants in the generator.
The validator checks the UAV/edge count, baseline range, measurement mode, path
multiplicity, and target ambiguity heights. For a single-transmitter,
multi-receiver system, use `measurement_mode = single_tx_multireceiver` and
`phase_path_multiplicity = 1`; the ambiguity heights are twice those of the
monostatic-pair profile for otherwise identical geometry.

The exact configuration and its derived physics are copied into every generated
dataset as `simulation_config.json` and `simulation_config_derived.mat`.

### Terrain-aware coherence and sampling

The V13 UAV generator no longer assigns random Voronoi regions as its primary
coherence model. `coherence/gen_coherence_map_multibaseline_terrainaware.m`
derives a deterministic terrain-quality field from DEM slope, aspect, local
incidence, roughness, curvature, TPI, and first-order layover/shadow proxies.
It then adds a smooth residual shared by all baselines plus an independent
residual for each interferometric edge. Longer baselines receive stronger
decorrelation.

Patch selection is random *within real terrain strata*, not random over an
invented block map. The default weighted policy accepts flat, rolling, steep,
and geometry-hazard patches with probabilities `0.35`, `0.70`, `1.0`, and
`1.0`, respectively. This avoids a flat-area-dominated corpus while retaining
geographic variety. The policy is configurable under
`dataset.terrain_sampling`; set its mode to `uniform` to disable weighting.
It is combined with `dataset.landcover.sampling`: a crop containing enough of
a configured real class (currently water, built-up, wetland, or cropland)
receives an acceptance multiplier. The largest eligible multiplier is used
and the final probability is capped at one. This over-samples rare real
surfaces without inventing a categorical block map.

Generated datasets now also contain:

- `patch_groups/<patch_name>.mat`: one flat MAT file per geographic patch;
- `[K,H,W]` phase and coherence arrays, including simulator
  `coherence_true` and observed `coherence_observed`;
- `edge_index`, `baseline_perp_m`, `ambiguity_height_m` and
  `valid_edge_mask` for graph-aware training;
- embedded terrain features: slope, aspect, local incidence, roughness,
  curvature, TPI, terrain quality and geometry masks;
- terrain class and diagnostics in `generation_detail_*.csv`.

Inactive edges are zero-filled and must be consumed together with
`valid_edge_mask`. The default configuration randomly retains 2--4 edges per
patch and always keeps the shortest baseline. Set
`interferometry.edge_sampling.enabled=false` for a fixed K. Legacy per-edge
folders can be restored with `dataset.storage.write_legacy_per_edge=true`, but
are disabled by default to avoid duplicating the dataset.

### Reproducible random streams

`dataset.generation.rng_seed` controls crop locations, scales, SNR draws, and
accept/reject decisions for the complete batch. Each candidate or saved patch
also receives deterministic, distinct component seeds for terrain-coherence
residuals, node SLC noise, and the valid-edge mask. Seeded component functions
restore the caller RNG state when they return, so replaying one component does
not freeze or restart random values elsewhere in the same batch.

The batch seed, final RNG state, and component counter are saved in
`rng_seed_*.mat`; every newly generated patch stores component seeds for
coherence, SLC, Phase-4 errors, edge sampling, and explicit failures in its
metadata and `generation_detail_*.csv`. Per-DEM patch counts are controlled by
`configs/dem_split_manifest.csv`; `dataset.generation.patches_per_dem` is the
fallback when a split-specific override is not used.

### Distributed complex-SLC noise

The default `noise.model=complex_slc_nodes` builds one complex SLC for every
UAV node. In the star graph, all interferometric edges reuse the same master
UAV speckle and receiver-noise realization, while each secondary UAV has its
own decorrelated scattering term and receiver noise. This creates cross-edge
error correlation naturally; the old independent phase-PDF path remains
available as `legacy_phase_pdf`.

Node SNR is sampled from the configured `node_snr_db_range` (15--30 dB by
default). `coherence_scene` records terrain/baseline decorrelation before
receiver noise, while `coherence_true` includes the thermal SNR factor. When
`noise.save_node_slc=true`, each patch group also stores complex
`node_slc[N,H,W]` for downstream interferogram reconstruction and diagnostics.
The pilot profile sets this to `false` because the derived phase/coherence
products are sufficient for routine QA and the node cube substantially
increases storage.

### Phase-4 synchronization, trajectory, and coregistration errors

`noise/apply_distributed_uav_errors.m` applies errors to UAV nodes before
forming interferometric edges. The current configurable MVP includes phase
bias, linear drift, random walk and sparse jumps; LOS position bias, drift and
vibration; roll/pitch-induced displacement; and range/azimuth subpixel shift
with slow drift. Edges sharing a UAV therefore inherit correlated errors.

Each patch retains `wrappedphase_node_noise_only` for ablation, the final
`wrappedphase_withnoise`, an edge-level `coregistration_valid_mask`, compact
`[N,H]` node error profiles, scalar node parameters, and a dedicated error
seed. The nominal profile is explicitly called
`pband_uav_nominal_v1_unvalidated`: its values are starting assumptions for
stress testing, not calibrated flight-system specifications.

Run `run_phase4_sweep` to generate the paired matrix declared in
`configs/phase4_sweep_matrix.json`. It contains `sync_only`,
`trajectory_only`, `coreg_only`, and `combined` scenarios at `nominal`,
`challenging` (3x), and `failure` (10x) levels. All cells deliberately reuse
seed 42, so matching patch IDs have identical terrain, SNR, coherence residual,
node noise, and edge mask. `analyze_phase4_sweep` verifies this pairing and
writes `data/phase4_sweep/phase4_sweep_summary.csv`.

The additional `compound_failure/failure` profile preserves the shortest and
longest baselines, drops one non-protected secondary UAV, forces synchronization
jump candidates, and labels low-coherence area from the real scene coherence.
It does not manufacture low coherence. See `docs/phase4_label_schema.md` for
field shapes, units, masks, and target definitions; run
`analyze_compound_failure_pilot` for QA.
The current complex-SLC implementation intentionally accepts star graphs only;
chain and complete graphs require a globally consistent node/edge phase model
and currently fall back to the legacy mode.

The nominal profile uses a physical critical-baseline approximation instead
of the earlier short-baseline-normalized exponential decay. With 500 MHz,
200 MHz bandwidth, 707.1 m slant range, and 45 degree incidence, the derived
slant-range resolution is about 0.749 m and the critical perpendicular
baseline is about 283 m. Therefore 0.75--3 m baselines receive only mild direct
geometric decorrelation. The `exponential` model remains selectable for stress
tests, but should not be treated as the nominal physics.

`noise.multilook_window=7` additionally saves
`wrappedphase_multilook[K,H,W]`, computed by complex averaging of node-derived
interferograms. Dense fringes on longer baselines can still decorrelate under a
fixed spatial window even when the baseline is far below the critical value;
this is a resolution/phase-gradient effect rather than baseline decorrelation.

`dataset.source_dem_pixel_spacing_m` must match the source DEM. The generator
automatically divides this spacing by its bicubic up-scaling factor when
computing terrain derivatives.

Terrain normalization now defaults to dataset-level `absolute_scales`, so the
same physical roughness receives the same penalty in every patch. The current
calibration sampled 96 scale/location combinations (786,432 pixels) from the
bundled ALPSMLC30 tile. Its P95 values were 20.58 degrees slope, 71.82 m local
roughness, and 0.006326 1/m absolute curvature. The nominal model keeps the
physical 45-degree slope reference and uses the latter two P95 values as fixed
roughness/curvature scales. The former `patch_percentile` normalization remains
available only as a fallback.

Run `scripts/calibrate_terrain_noise_scales.m` whenever the training DEM
collection or ground resolution changes. Calibration must use training tiles
only; using validation/test tiles would leak their terrain distribution into
the generator configuration. The calibration JSON/MAT records the seed,
sample count, percentiles, and recommended scales.

### Real land-cover conditioning

All 19 configured DEMs are aligned with the official ESA WorldCover 10 m
2021 v200 maps. Eight 3-degree categorical EPSG:4326 COGs cover the 19 DEM
tiles; `scripts/prepare_worldcover_for_dem.m` selects the correct source tile
from each DEM coordinate and samples each DEM cell center with nearest-neighbor
assignment. It also handles ASTER's point-registered 3601-by-3601 shared tile
boundaries without interpolating class codes. The configuration uses
`missing_policy=error`, so a newly added DEM cannot silently fall back to
another tile or a neutral surface factor.

`coherence_terrain_only` stores coherence before land-cover adjustment.
`coherence_scene` stores terrain coherence multiplied by the configured class
factor, while `landcover_codes` and `landcover_factor` are embedded in every
patch group. This supports terrain-only versus terrain-plus-surface ablation.

The default factors are initial P-band single-pass hypotheses, not values
measured by WorldCover and not universal constants. The profile is explicitly
named `initial_pband_singlepass_v1_unvalidated`. Water receives a strong
reduction, bare ground no additional reduction, and vegetation/built classes
intermediate factors. These values require real-data fitting or sensitivity
sweeps before publication.

---

## Requirements

- MATLAB R2019b or later  
- No additional toolboxes required (Image Processing Toolbox is **not** needed)
- GeoTIFF reading uses the built-in `imread` (MATLAB reads single-band GeoTIFFs natively from R2019b)

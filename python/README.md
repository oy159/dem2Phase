# dem2phase Python cloud generator

This package ports the active MATLAB generation path to Python:
`terrain_aware + complex_slc_nodes + grouped_mat + Phase-4`.

## Reproducibility contract

The Python generator is array-reproducible for a fixed package version,
configuration, input rasters, and batch seed. Random streams are derived with
BLAKE2b and PCG64DXSM from the geographic tile, DEM filename, candidate index,
and component name. Results therefore do not depend on process scheduling or
the `--workers` value.

MATLAB and Python do not emit identical random arrays from the same integer
seed. Cross-language acceptance is based on schema, deterministic physics
tolerances, and dataset-level statistics. `compare-matlab` uses an absolute
mean-coherence tolerance of 0.03 and a P05/P50/P95 tolerance of 0.05.

## Local installation

```bash
cd python
python -m pip install -e '.[test]'
python -m pytest
```

Generate and validate the current 19-DEM dataset from the repository root:

```bash
python -m dem2phase generate \
  --config configs/uav_p_500m_monostatic.json \
  --seed 42 --workers 8 \
  --output data/pilot_dataset_uav_p_500m_phase4_python

python -m dem2phase validate \
  --dataset data/pilot_dataset_uav_p_500m_phase4_python

python -m dem2phase compare-matlab \
  --matlab-dataset data/pilot_dataset_uav_p_500m_phase4 \
  --python-dataset data/pilot_dataset_uav_p_500m_phase4_python
```

Generate only the predefined test split while preserving the full manifest's
global patch IDs:

```bash
python -m dem2phase generate \
  --config configs/uav_p_500m_monostatic.json \
  --split test --seed 42 --workers 4 \
  --output data/local_test_python
```

Render a patch overview, a terrain/land-cover panel, and a JSON summary:

```bash
python -m dem2phase visualize \
  --dataset data/local_test_python --split test --index 1 \
  --colormap jet \
  --output data/local_test_python/visualizations/test_patch_001.png
```

`jet` is the default colormap. Pass another Matplotlib colormap name through
`--colormap` when a cyclic or perceptually uniform view is needed.

Select by filename text or render particular one-based edge slots:

```bash
python -m dem2phase visualize \
  --dataset data/local_test_python --split test --match N32E078 \
  --edge 1 --edge 4 --output visualizations/N32E078.png

python -m dem2phase visualize \
  --patch data/local_test_python/patch_groups/test/example.mat \
  --output visualizations/example.png
```

## Cloud 10x profile and local recovery

`configs/uav_p_500m_cloud_10x.json` inherits the production physics from the
local profile and multiplies every manifest patch count by 10. It generates
2,340 grouped MAT samples: 1,620 train, 560 test, and 160 validation samples.
The expected output size is approximately 29--31 GiB. The local profile stays
at 234 samples. Cloud candidates whose source-DEM footprint overlaps an
already accepted patch by more than 25% are rejected; the measured overlap and
`spatial_overlap` rejection reason are retained in `generation_plan.csv`.

```bash
python -m dem2phase generate \
  --config configs/uav_p_500m_cloud_10x.json \
  --seed 42 --workers 8 \
  --output /data/dem2phase_cloud_10x
```

Every run writes `input_inventory.csv` with the byte size and SHA-256 of the
split manifest, selected DEMs, ASTER NUM rasters, and aligned WorldCover MAT
files. `provenance.json` also records the effective config hash, generator
source fingerprint, Git commit, RNG algorithm, Python version, and the pinned
NumPy/SciPy/Pandas/Rasterio versions. `recovery_recipe.json` contains the exact
seed and command templates.

Cloud MAT files do not have to be downloaded to recover a deterministic
subset. With the same repository revision and locally retained input rasters,
regenerate one DEM shard directly on the local machine:

```bash
python -m dem2phase generate \
  --config configs/uav_p_500m_cloud_10x.json \
  --seed 42 \
  --dem-file Copernicus_DSM_COG_10_N27_00_E084_00_DEM.tif \
  --output data/recovered_N27E084
```

`--dem-file` preserves the DEM's position and global offsets from the full
manifest. Patch names, global IDs, crop coordinates, component seeds, and all
numeric arrays therefore match the corresponding full cloud run. Repeat the
option to recover several DEM shards. A different generator fingerprint or an
input SHA-256 mismatch means exact recovery is not guaranteed.

Prepare WorldCover into a separate validation directory and compare it with
the MATLAB-aligned rasters without overwriting either set:

```bash
python -m dem2phase prepare-landcover \
  --config configs/uav_p_500m_monostatic.json \
  --raw-root data/landcover/raw \
  --output data/landcover_python_validation

python python/tools/compare_landcover.py \
  data/landcover data/landcover_python_validation
```

The comparison is pixelwise. On the current 19 DEMs the Python cell-centre
implementation differs on 0.0264% of 246,268,804 pixels (maximum 0.061% for
one DEM), at floating-point half-pixel class boundaries; all four ASTER
rasters are identical.

Rebuild and verify the deterministic MATLAB/Python golden fixture:

```matlab
export_cross_language_golden( ...
    'python/tests/fixtures/matlab_cross_language_golden.mat')
```

```bash
python python/tools/compare_cross_language_golden.py \
  --fixture python/tests/fixtures/matlab_cross_language_golden.mat \
  --config configs/uav_p_500m_monostatic.json \
  --report python/tests/fixtures/golden_comparison_report.json
```

This checks bicubic crop, phase physics, terrain features, categorical nearest
sampling, convolution and linear registration. The current fixture passes;
its maximum DEM, phase and deterministic-coherence differences are 1.11 mm,
1.39e-4 rad and 2.22e-4 respectively. Exact categorical and mask fields match.

The output directory must be empty. Use `--resume` only for the same config
hash, batch seed, package version, and manifest. Existing completed MAT files
are retained and missing files are regenerated atomically. If the plan is
already complete and every accepted MAT exists, resume returns before loading
any DEM (0.762 seconds in the current 234-sample local test).

For a one-DEM cloud smoke test, keep the production configuration unchanged:

```bash
python -m dem2phase generate --config configs/uav_p_500m_monostatic.json \
  --dem-file ALPSMLC30_N027E085_DSM.tif --patches-per-dem 1 \
  --seed 42 --workers 1 --output data/cloud_smoke
```

## Docker

Build once:

```bash
docker build -t dem2phase-cloud:0.1.0 python
```

Run from the repository root on Linux:

```bash
docker run --rm \
  --user "$(id -u):$(id -g)" \
  -v "$PWD:/workspace" -w /workspace \
  dem2phase-cloud:0.1.0 generate \
  --config configs/uav_p_500m_monostatic.json \
  --seed 42 --workers 8 \
  --output data/pilot_dataset_uav_p_500m_phase4_python
```

Only Python 3.11 plus the locked packages are required in the container. The
input DEMs, ASTER NUM files, aligned WorldCover files, configuration, and split
manifest are mounted as volumes; MATLAB is not required. Generated patch
groups stay in the mounted output volume or cloud object-storage mount.
Passing the host UID/GID prevents bind-mounted results from being owned by
root.

## Commands

- `generate`: deterministic multi-process grouped-MAT generation.
- `prepare-landcover`: nearest-neighbour alignment from raw WorldCover COGs.
- `validate`: schema, inactive-edge, split-leakage, plan, and seed validation.
- `compare-matlab`: pooled MATLAB/Python coherence and edge-count comparison;
  per-source numbers are diagnostic because some ASTER strata have only three
  samples.
- `tools/compare_cross_language_golden.py`: strict deterministic cross-language
  numerical regression test.

Every run records `run_state.json`, `provenance.json`, `input_inventory.csv`,
`recovery_recipe.json`, the effective and derived configurations,
`rng_manifest.mat`, `generation_plan.csv`, three generation logs, and
failure-balanced training manifests.

## Verified container smoke test

The image was built and tested on Ubuntu 24.04 with Docker 29.4.0. A real
`ALPSMLC30_N027E085_DSM.tif` input generated one four-baseline 256-by-256 MAT
sample, passed container validation, and passed a same-seed `--resume` run.
Running with the host UID/GID produced host-owned bind-mount files. The image
size is approximately 167 MB. Windows and Linux produced
the same semantic configuration hash and package source fingerprint; a
complete one-DEM container resume took 0.95 seconds after the terrain-boundary
fix.

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

Only the input DEMs, ASTER NUM files, aligned WorldCover files, configuration,
and split manifest need to be uploaded. Generated patch groups stay in the
mounted output volume or cloud object-storage mount. Passing the host UID/GID
prevents bind-mounted results from being owned by root.

## Commands

- `generate`: deterministic multi-process grouped-MAT generation.
- `prepare-landcover`: nearest-neighbour alignment from raw WorldCover COGs.
- `validate`: schema, inactive-edge, split-leakage, plan, and seed validation.
- `compare-matlab`: pooled MATLAB/Python coherence and edge-count comparison;
  per-source numbers are diagnostic because some ASTER strata have only three
  samples.
- `tools/compare_cross_language_golden.py`: strict deterministic cross-language
  numerical regression test.

Every run records `run_state.json`, `provenance.json`, the effective and
derived configurations, `rng_manifest.mat`, `generation_plan.csv`, three
generation logs, and failure-balanced training manifests.

## Verified container smoke test

The image was built and tested on Ubuntu 24.04 with Docker 29.4.0. A real
`ALPSMLC30_N027E085_DSM.tif` input generated one four-baseline 256-by-256 MAT
sample, passed container validation, and passed a same-seed `--resume` run.
Running with the host UID/GID produced host-owned bind-mount files. The image
size is approximately 167 MB. Windows and Linux produced
the same semantic configuration hash and package source fingerprint; a
complete one-DEM container resume took 0.95 seconds after the terrain-boundary
fix.

from __future__ import annotations

import json
import copy
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import rasterio
from rasterio.transform import from_origin
from scipy.io import loadmat, savemat

from dem2phase.config import load_config
from dem2phase.core import sample_edge_mask, terrain_coherence, wrap_phase
from dem2phase.generator import config_hash, generate_dataset
from dem2phase.io import load_split_manifest
from dem2phase.landcover import worldcover_token
from dem2phase.rng import derive_seed, generator
from dem2phase.validation import validate_dataset
from dem2phase.visualization import visualize_patch
import dem2phase.generator as generator_module


def test_hierarchical_seed_is_stable_and_distinct():
    assert derive_seed(42, "N27E085", 1, "slc") == derive_seed(42, "N27E085", 1, "slc")
    assert len({derive_seed(42, "N27E085", 1, name)
                for name in ("crop", "slc", "coherence", "failure")}) == 4
    assert np.array_equal(generator(42, "x").normal(size=10),
                          generator(42, "x").normal(size=10))


def test_wrap_and_edge_mask_replay():
    phase = np.linspace(-20, 20, 100).reshape(10, 10)
    assert np.all((wrap_phase(phase) >= -np.pi) & (wrap_phase(phase) <= np.pi))
    options = {"enabled": True, "min_active_edges": 2, "max_active_edges": 4,
               "always_include_shortest": True}
    left = sample_edge_mask(4, 0, options, generator(12, "mask"))
    right = sample_edge_mask(4, 0, options, generator(12, "mask"))
    assert np.array_equal(left, right) and left[0] and 2 <= left.sum() <= 4


def test_worldcover_three_degree_tokens():
    assert worldcover_token(29.5, 102.5) == "N27E102"
    assert worldcover_token(-32.5, -70.0) == "S33W072"


def test_manifest_rejects_cross_sensor_leakage(tmp_path: Path):
    frame = pd.DataFrame([
        _manifest_row("a_DEM.tif", "N27E085", "train"),
        _manifest_row("b_DEM.tif", "N27E085", "test"),
    ])
    path = tmp_path / "split.csv"
    frame.to_csv(path, index=False)
    with pytest.raises(ValueError, match="leakage"):
        load_split_manifest(path)


def test_one_and_two_workers_are_identical(tmp_path: Path):
    config_path = _fixture_project(tmp_path)
    cfg = load_config(config_path)
    one, two = tmp_path / "out_one", tmp_path / "out_two"
    generate_dataset(cfg, one, 42, workers=1)
    generate_dataset(cfg, two, 42, workers=2)
    report_one, report_two = validate_dataset(one), validate_dataset(two)
    assert report_one["samples"] == report_two["samples"] == 4
    assert report_one["split_counts"] == {"train": 2, "test": 2}
    assert report_one["array_hashes"] == report_two["array_hashes"]
    plan_one = pd.read_csv(one / "generation_plan.csv")
    plan_two = pd.read_csv(two / "generation_plan.csv")
    pd.testing.assert_frame_equal(plan_one, plan_two)
    original_hashes = report_one["array_hashes"]
    missing = next((one / "patch_groups").glob("*/*.mat"))
    missing.unlink()
    generate_dataset(cfg, one, 42, workers=2, resume=True)
    assert validate_dataset(one)["array_hashes"] == original_hashes


def test_single_dem_smoke_override(tmp_path: Path):
    cfg = load_config(_fixture_project(tmp_path))
    output = tmp_path / "smoke"
    generate_dataset(cfg, output, 42, dem_files=["fixture_N00E000_DEM.tif"],
                     patches_per_dem=1)
    report = validate_dataset(output)
    assert report["samples"] == 1
    assert report["split_counts"] == {"train": 1}
    sample_path = next((output / "patch_groups" / "train").glob("*.mat"))
    sample = loadmat(sample_path, squeeze_me=True, struct_as_record=False)
    node_parameters = sample["phase4_errors"].node_parameters
    assert isinstance(node_parameters, np.ndarray)
    assert node_parameters.size == 3
    assert all(hasattr(node, "sync_jump_rad") for node in node_parameters.ravel())


def test_selected_dem_shard_preserves_full_run_identity(tmp_path: Path):
    cfg = load_config(_fixture_project(tmp_path))
    full, shard = tmp_path / "full", tmp_path / "shard"
    generate_dataset(cfg, full, 42, workers=1)
    selected = "fixture_N00E001_DEM.tif"
    generate_dataset(cfg, shard, 42, workers=1, dem_files=[selected])

    full_files = sorted((full / "patch_groups" / "test").glob("*.mat"))
    shard_files = sorted((shard / "patch_groups" / "test").glob("*.mat"))
    assert [path.name for path in shard_files] == [path.name for path in full_files]
    assert [path.name for path in shard_files] == [
        "fixture_N00E001_DEM_patch_00003.mat",
        "fixture_N00E001_DEM_patch_00004.mat",
    ]
    array_fields = ("wrappedphase_withoutnoise", "wrappedphase_withnoise",
                    "unwrapped_phase", "coherence_observed", "coherence_true",
                    "valid_edge_mask", "landcover_codes")
    for full_path, shard_path in zip(full_files, shard_files):
        full_sample = loadmat(full_path, squeeze_me=True, struct_as_record=False)
        shard_sample = loadmat(shard_path, squeeze_me=True, struct_as_record=False)
        for field in array_fields:
            assert np.array_equal(full_sample[field], shard_sample[field])
        assert (full_sample["metadata"].patch_global_id ==
                shard_sample["metadata"].patch_global_id)
    assert (shard / "input_inventory.csv").is_file()
    assert (shard / "recovery_recipe.json").is_file()


def test_split_selection_preserves_global_identity_and_visualizes(tmp_path: Path):
    cfg = load_config(_fixture_project(tmp_path))
    output = tmp_path / "test_only"
    generate_dataset(cfg, output, 42, workers=1, splits=["test"])
    report = validate_dataset(output)
    assert report["samples"] == 2
    files = sorted((output / "patch_groups" / "test").glob("*.mat"))
    assert [path.name for path in files] == [
        "fixture_N00E001_DEM_patch_00003.mat",
        "fixture_N00E001_DEM_patch_00004.mat",
    ]
    rendered = visualize_patch(None, output, split="test", index=1, dpi=40)
    assert Path(rendered["overview_png"]).is_file()
    assert Path(rendered["terrain_png"]).is_file()
    assert Path(rendered["summary_json"]).is_file()


def test_config_extends_and_patch_multiplier(tmp_path: Path):
    base_path = _fixture_project(tmp_path)
    child = tmp_path / "configs" / "cloud.json"
    child.write_text(json.dumps({
        "extends": base_path.name,
        "name": "fixture_cloud",
        "dataset": {"generation": {"patch_multiplier": 10}},
    }), encoding="utf-8")
    cfg = load_config(child)
    assert cfg["name"] == "fixture_cloud"
    assert cfg["dataset"]["generation"]["patch_size"] == 16
    assert cfg["dataset"]["generation"]["patch_multiplier"] == 10


def test_source_overlap_fraction():
    overlap = generator_module._source_overlap_fraction
    assert overlap((0, 0, 10, 10), (20, 20, 30, 30)) == 0
    assert overlap((0, 0, 10, 10), (5, 0, 15, 10)) == pytest.approx(0.5)
    assert overlap((0, 0, 10, 10), (2, 2, 8, 8)) == pytest.approx(1.0)


def test_complete_resume_skips_dem_computation(tmp_path: Path, monkeypatch):
    cfg = load_config(_fixture_project(tmp_path))
    output = tmp_path / "complete"
    generate_dataset(cfg, output, 42, workers=1)

    def unexpected_generation(_task):
        raise AssertionError("a complete resume must not recompute a DEM")

    monkeypatch.setattr(generator_module, "_generate_dem", unexpected_generation)
    result = generate_dataset(cfg, output, 42, workers=1, resume=True)
    assert result["resumed_complete"] is True


def test_config_hash_ignores_recomputed_float_roundoff(tmp_path: Path):
    cfg = load_config(_fixture_project(tmp_path))
    perturbed = copy.deepcopy(cfg)
    perturbed["physics"]["wavelength_m"] = np.nextafter(
        perturbed["physics"]["wavelength_m"], np.inf
    )
    perturbed["interferometry"]["edge_index"] = [[99], [100]]
    assert config_hash(cfg) == config_hash(perturbed)

    changed = copy.deepcopy(cfg)
    changed["radar"]["center_frequency_hz"] += 1
    assert config_hash(cfg) != config_hash(changed)


def test_matlab_cross_language_golden():
    root = Path(__file__).resolve().parents[1]
    subprocess.run([
        sys.executable, str(root / "tools" / "compare_cross_language_golden.py"),
        "--fixture", str(root / "tests" / "fixtures" / "matlab_cross_language_golden.mat"),
        "--config", str(root.parent / "configs" / "uav_p_500m_monostatic.json"),
    ], check=True, capture_output=True, text=True)


def _manifest_row(name: str, tile: str, split: str) -> dict:
    return {"dem_file": name, "geographic_tile": tile, "split": split,
            "region": "fixture", "quality_role": "primary_copernicus",
            "patches_per_dem": 2, "latitude_deg": 0, "longitude_deg": 0,
            "source": "fixture", "expected_bytes": 0, "sha256": "", "download_url": ""}


def _fixture_project(root: Path) -> Path:
    (root / "configs").mkdir()
    (root / "data" / "dem").mkdir(parents=True)
    (root / "data" / "landcover").mkdir(parents=True)
    records = []
    for idx, (tile, split) in enumerate((("N00E000", "train"), ("N00E001", "test"))):
        name = f"fixture_{tile}_DEM.tif"
        dem = (np.add.outer(np.arange(32), np.arange(32)) + idx * 10).astype(np.float32)
        path = root / "data" / "dem" / name
        with rasterio.open(path, "w", driver="GTiff", height=32, width=32, count=1,
                           dtype="float32", crs="EPSG:4326",
                           transform=from_origin(idx, 1, 1 / 32, 1 / 32)) as dst:
            dst.write(dem, 1)
        savemat(root / "data" / "landcover" / f"{path.stem}_worldcover2021_aligned.mat",
                {"landcover_codes": np.full((32, 32), 60, np.uint8)})
        records.append(_manifest_row(name, tile, split))
    pd.DataFrame(records).to_csv(root / "configs" / "dem_split_manifest.csv", index=False)
    cfg = {
        "schema_version": 1,
        "dataset": {"dem_directory": "data/dem", "split_manifest": "configs/dem_split_manifest.csv",
            "require_split_assignment": True, "output_directory": "data/out",
            "source_dem_pixel_spacing_m": 30,
            "generation": {"patches_per_dem": 2, "patch_size": 16, "n_precomp_scales": 2,
                "interp_scale_range": [1, 2], "rng_seed": 42, "max_wrap_count": 100,
                "void_fraction_threshold": 0.05},
            "landcover": {"enabled": True, "aligned_file": "", "aligned_directory": "data/landcover",
                "aligned_suffix": "_worldcover2021_aligned.mat", "missing_policy": "error",
                "source": "fixture", "profile": "fixture", "unknown_factor": 1,
                "class_codes": [60, 80], "coherence_factors": [1, 0.1],
                "sampling": {"enabled": False, "priority_codes": [80],
                    "minimum_fraction": [0.01], "acceptance_multiplier": [2]}},
            "storage": {"mode": "grouped_mat", "write_legacy_per_edge": False,
                        "numeric_type": "single"},
            "terrain_sampling": {"mode": "uniform", "max_attempts_factor": 4,
                "acceptance_probability": {"flat": 1, "rolling": 1, "steep": 1,
                                           "geometry_hazard": 1}}},
        "radar": {"center_frequency_hz": 500e6, "bandwidth_hz": 200e6},
        "geometry": {"altitude_m": 500, "incidence_angle_deg": 45,
                     "look_azimuth_deg": 90, "slant_range_m": 707.1067811865476},
        "coherence": {"model": "terrain_aware", "baseline_decay_model": "critical_baseline",
            "terrain_normalization_model": "absolute_scales", "min": 0.01, "max": 0.95,
            "decay_alpha": 0.35, "slope_weight": 0.9, "roughness_weight": 0.8,
            "curvature_weight": 0.35, "slope_scale_deg": 45, "roughness_scale_m": 70,
            "curvature_scale_per_m": 0.006, "shared_residual_std": 0.12,
            "edge_residual_std": 0.08, "residual_scale_px": 4},
        "noise": {"model": "complex_slc_nodes", "node_snr_db_range": [15, 30],
                  "multilook_window": 3, "save_node_slc": False},
        "phase4_errors": {"enabled": True, "profile": "fixture",
            "synchronization": {"phase_bias_std_rad": 0.02, "linear_drift_std_rad": 0.03,
                "random_walk_std_rad_per_row": 0.002, "jump_probability_per_node": 0.1,
                "jump_std_rad": 0.15},
            "trajectory": {"los_bias_std_m": 0.01, "los_drift_std_m": 0.005,
                "vibration_std_m": 0.002, "vibration_cycles_range": [1, 5]},
            "attitude": {"roll_std_deg": 0.02, "pitch_std_deg": 0.02},
            "coregistration": {"range_shift_std_px": 0.05, "azimuth_shift_std_px": 0.05,
                               "linear_drift_std_px": 0.02}},
        "phase4_failures": {"enabled": False, "profile": "fixture",
            "force_secondary_uav_dropout": True, "always_include_shortest": True,
            "always_include_longest": True, "low_coherence_threshold": 0.25,
            "low_coherence_min_fraction": 0.05, "sync_jump_min_abs_rad": 0.1},
        "interferometry": {"measurement_mode": "monostatic_pair", "phase_path_multiplicity": 2,
            "num_uavs": 3, "edge_mode": "star", "reference_uav": 1,
            "baseline_perp_m": [0.75, 1.2],
            "edge_sampling": {"enabled": True, "min_active_edges": 1,
                              "max_active_edges": 2, "always_include_shortest": True}}
    }
    path = root / "configs" / "fixture.json"
    path.write_text(json.dumps(cfg), encoding="utf-8")
    return path

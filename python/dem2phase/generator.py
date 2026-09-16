"""End-to-end deterministic dataset generation and manifests."""

from __future__ import annotations

import concurrent.futures
import copy
import hashlib
import json
import math
import os
import platform
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from . import __version__
from .config import resolve_path
from .core import (apply_landcover, apply_phase4, coherence_from_slc, failure_mask,
                   multilook_phases, sample_edge_mask, simulate_node_slc, stack,
                   summarize_terrain, terrain_coherence, wrap_phase)
from .io import (atomic_savemat, cubic_crop, discover_dems, load_aligned_landcover,
                 load_split_manifest, nearest_source_indices, read_raster,
                 spline_coefficients)
from .rng import RNG_ALGORITHM, derive_seed, generator


PLAN_COLUMNS = [
    "manifest_index", "candidate_index", "accepted_index", "accepted", "reject_reason",
    "patch_global_id", "patch_name", "patch_group_file", "source_file", "geographic_tile",
    "split", "region", "source_dataset", "source_quality_role", "interp_scale",
    "scaled_rows", "scaled_cols", "crop_r0", "crop_c0", "crop_r1", "crop_c1",
    "candidate_seed", "coherence_seed", "slc_seed", "phase4_error_seed",
    "edge_mask_seed", "failure_seed", "accept_seed", "node_snr_db", "valid_edge_mask",
    "terrain_class", "dominant_landcover_code", "landcover_sampling_multiplier",
    "priority_landcover_code", "num_available", "num_void_fraction",
    "num_mean", "num_min", "max_wrap_count", "active_edge_count", "dropout_present",
    "sync_anomaly_present", "low_coherence_present", "compound_failure_present",
    "failure_factor_count", "accept_probability", "ground_pixel_spacing_m",
]


def _config_for_hash(cfg: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in cfg.items() if key not in {"config_path", "project_root"}}


def config_hash(cfg: dict[str, Any]) -> str:
    # Hash user-controlled semantics only. Derived trigonometric values can
    # differ by an ulp across libm/Python builds and must not make an otherwise
    # identical Windows/Linux job impossible to resume.
    value = copy.deepcopy(_config_for_hash(cfg))
    value.pop("physics", None)
    interferometry = value.get("interferometry", {})
    interferometry.pop("edge_index", None)
    interferometry.pop("num_edges", None)
    raw = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def source_fingerprint() -> str:
    digest = hashlib.sha256()
    package_root = Path(__file__).resolve().parent
    for path in sorted(package_root.glob("*.py")):
        digest.update(path.name.encode("utf-8"))
        digest.update(path.read_bytes())
    return digest.hexdigest()


def _git_commit(root: Path) -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root,
                                       stderr=subprocess.DEVNULL, text=True).strip()
    except (OSError, subprocess.SubprocessError):
        return "unknown"


def _state(cfg: dict[str, Any], seed: int, manifest: pd.DataFrame) -> dict[str, Any]:
    return {"schema_version": 1, "generator": "dem2phase-python",
            "generator_version": __version__, "config_sha256": config_hash(cfg),
            "generator_fingerprint": source_fingerprint(),
            "batch_seed": int(seed), "rng_algorithm": RNG_ALGORITHM,
            "manifest_rows": int(len(manifest)),
            "dem_files": list(manifest["dem_file"].astype(str)),
            "target_patches": int(manifest["patches_per_dem"].sum())}


def _prepare_output(output: Path, state: dict[str, Any], resume: bool) -> None:
    state_path = output / "run_state.json"
    if output.exists() and any(output.iterdir()):
        if not resume:
            raise FileExistsError(f"Output is not empty: {output}; use --resume")
        if not state_path.exists():
            raise ValueError("Cannot resume an output without run_state.json")
        saved = json.loads(state_path.read_text(encoding="utf-8"))
        for key in ("schema_version", "generator_version", "config_sha256", "batch_seed",
                    "generator_fingerprint", "rng_algorithm", "manifest_rows", "dem_files",
                    "target_patches"):
            if saved.get(key) != state.get(key):
                raise ValueError(f"Resume state mismatch for {key}")
    else:
        output.mkdir(parents=True, exist_ok=True)
        for split in ("train", "test", "validation"):
            (output / "patch_groups" / split).mkdir(parents=True, exist_ok=True)
        state_path.write_text(json.dumps(state, indent=2), encoding="utf-8")


def _complete_plan(output: Path, target_patches: int) -> pd.DataFrame | None:
    """Return an existing complete plan when every accepted MAT is present."""
    plan_path = output / "generation_plan.csv"
    if not plan_path.is_file():
        return None
    plan = pd.read_csv(plan_path)
    if "accepted" not in plan or "patch_group_file" not in plan:
        return None
    accepted = plan[plan["accepted"].astype(str).str.lower() == "true"]
    if len(accepted) != target_patches:
        return None
    if not all((output / Path(relative)).is_file()
               for relative in accepted["patch_group_file"].astype(str)):
        return None
    return plan


def generate_dataset(cfg: dict[str, Any], output: Path, seed: int,
                     workers: int = 1, resume: bool = False,
                     dem_files: list[str] | None = None,
                     patches_per_dem: int | None = None) -> dict[str, Any]:
    manifest_path = resolve_path(cfg, cfg["dataset"]["split_manifest"])
    dem_root = resolve_path(cfg, cfg["dataset"]["dem_directory"])
    manifest = load_split_manifest(manifest_path)
    if dem_files:
        unknown = sorted(set(dem_files) - set(manifest["dem_file"]))
        if unknown:
            raise ValueError(f"Requested DEMs are absent from the split manifest: {unknown}")
        order = {name: index for index, name in enumerate(dem_files)}
        manifest = manifest[manifest["dem_file"].isin(dem_files)].copy()
        manifest["_selection_order"] = manifest["dem_file"].map(order)
        manifest = manifest.sort_values("_selection_order").drop(columns="_selection_order")
    if patches_per_dem is not None:
        if patches_per_dem < 1:
            raise ValueError("patches_per_dem override must be positive")
        manifest["patches_per_dem"] = int(patches_per_dem)
    dem_paths = discover_dems(dem_root, manifest, allow_extras=bool(dem_files))
    state = _state(cfg, seed, manifest)
    _prepare_output(output, state, resume)
    if resume and _complete_plan(output, state["target_patches"]) is not None:
        return {"patches": state["target_patches"], "dems": int(len(manifest)),
                "output": str(output), "config_sha256": state["config_sha256"],
                "resumed_complete": True}
    effective = _config_for_hash(cfg)
    (output / "simulation_config_effective.json").write_text(
        json.dumps(effective, indent=2), encoding="utf-8")
    atomic_savemat(output / "simulation_config_derived.mat", {"sim_cfg": effective})
    provenance = {**state, "python": platform.python_version(), "platform": platform.platform(),
                  "numpy": np.__version__, "git_commit": _git_commit(Path(cfg["project_root"]))}
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    atomic_savemat(output / "rng_manifest.mat", {
        "rng_seed_used": np.uint64(seed), "rng_algorithm": RNG_ALGORITHM,
        "config_sha256": state["config_sha256"]})

    offsets = np.concatenate(([0], np.cumsum(manifest["patches_per_dem"].to_numpy(dtype=int))))
    tasks = [(cfg, row.to_dict(), str(path), str(output), int(seed), index,
              int(offsets[index]), resume)
             for index, (path, (_, row)) in enumerate(zip(dem_paths, manifest.iterrows()))]
    results: list[dict[str, Any]] = []
    if workers <= 1:
        results = [_generate_dem(task) for task in tasks]
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(_generate_dem, task) for task in tasks]
            results = [future.result() for future in futures]
    results.sort(key=lambda item: item["manifest_index"])
    plan_rows = [row for result in results for row in result["plan"]]
    plan = pd.DataFrame(plan_rows, columns=PLAN_COLUMNS).sort_values(
        ["manifest_index", "candidate_index"])
    plan.to_csv(output / "generation_plan.csv", index=False, float_format="%.12g")
    accepted = plan[plan["accepted"] == True].copy()  # noqa: E712
    if len(accepted) != state["target_patches"]:
        raise RuntimeError(f"Generated {len(accepted)} of {state['target_patches']} target patches")
    _write_logs(output, results, accepted, seed)
    _write_balanced_manifest(output, accepted)
    return {"patches": int(len(accepted)), "dems": int(len(manifest)),
            "output": str(output), "config_sha256": state["config_sha256"]}


def _generate_dem(task: tuple[Any, ...]) -> dict[str, Any]:
    cfg, row, dem_path_text, output_text, batch_seed, manifest_index, offset, resume = task
    dem_path, output = Path(dem_path_text), Path(output_text)
    dem, _ = read_raster(dem_path)
    rows, cols = dem.shape
    coefficients = spline_coefficients(dem)
    patch_size = int(cfg["dataset"]["generation"]["patch_size"])
    n_scales = int(cfg["dataset"]["generation"]["n_precomp_scales"])
    scale_min, scale_max = cfg["dataset"]["generation"]["interp_scale_range"]
    scales = np.linspace(float(scale_min), float(scale_max), n_scales)
    target = int(row["patches_per_dem"])
    max_attempts = int(math.ceil(target * float(cfg["dataset"]["terrain_sampling"]["max_attempts_factor"])))
    numeric_type = np.float32 if cfg["dataset"]["storage"]["numeric_type"] == "single" else np.float64
    baselines = np.asarray(cfg["interferometry"]["baseline_perp_m"], dtype=np.float64)
    ratios = np.asarray(cfg["physics"]["dem2phase_ratios_rad_per_m"], dtype=np.float64)
    edge_index = np.asarray(cfg["interferometry"]["edge_index"], dtype=np.uint16)
    shortest = int(np.argmin(baselines))
    stem = dem_path.stem
    lc_cfg = cfg["dataset"]["landcover"]
    lc_path = resolve_path(cfg, lc_cfg["aligned_directory"]) / f"{stem}{lc_cfg['aligned_suffix']}"
    if lc_cfg["enabled"] and not lc_path.exists() and lc_cfg["missing_policy"] == "error":
        raise FileNotFoundError(f"Missing aligned land cover: {lc_path}")
    landcover = load_aligned_landcover(lc_path) if lc_cfg["enabled"] and lc_path.exists() else np.zeros(dem.shape, np.uint8)
    if landcover.shape != dem.shape:
        raise ValueError(f"Land-cover shape does not match {dem_path.name}")
    num_path = _num_path(dem_path)
    num_map = read_raster(num_path, np.uint8)[0] if num_path else None
    plan: list[dict[str, Any]] = []
    accepted_count = 0
    for candidate in range(1, max_attempts + 1):
        base = (row["geographic_tile"], dem_path.name, candidate)
        candidate_seed = derive_seed(batch_seed, *base, "candidate")
        seeds = {name: derive_seed(batch_seed, *base, name) for name in
                 ("coherence", "slc", "phase4_error", "edge_mask", "failure", "accept")}
        crop_rng = generator(batch_seed, *base, "crop")
        scale = float(scales[int(crop_rng.integers(0, n_scales))])
        scaled_rows, scaled_cols = int(round(scale * rows)), int(round(scale * cols))
        r0 = int(crop_rng.integers(0, max(1, scaled_rows - patch_size)))
        c0 = int(crop_rng.integers(0, max(1, scaled_cols - patch_size)))
        bounds = (r0, c0, r0 + patch_size - 1, c0 + patch_size - 1)
        base_plan = {"manifest_index": manifest_index + 1, "candidate_index": candidate,
                     "accepted_index": 0, "accepted": False, "reject_reason": "",
                     "patch_global_id": 0, "patch_name": "", "patch_group_file": "",
                     "source_file": dem_path.name, "geographic_tile": row["geographic_tile"],
                     "split": row["split"], "region": row["region"],
                     "source_dataset": row["source"], "source_quality_role": row["quality_role"],
                     "interp_scale": scale, "scaled_rows": scaled_rows, "scaled_cols": scaled_cols,
                     "crop_r0": r0 + 1, "crop_c0": c0 + 1,
                     "crop_r1": bounds[2] + 1, "crop_c1": bounds[3] + 1,
                     "candidate_seed": candidate_seed, **{f"{k}_seed": v for k, v in seeds.items()},
                     "node_snr_db": "", "valid_edge_mask": "", "terrain_class": "",
                     "dominant_landcover_code": 0, "landcover_sampling_multiplier": 1.0,
                     "priority_landcover_code": 0, "num_available": num_map is not None,
                     "num_void_fraction": np.nan, "num_mean": np.nan, "num_min": np.nan,
                     "max_wrap_count": np.nan, "active_edge_count": 0,
                     "dropout_present": False, "sync_anomaly_present": False,
                     "low_coherence_present": False, "compound_failure_present": False,
                     "failure_factor_count": 0, "accept_probability": np.nan,
                     "ground_pixel_spacing_m": float(cfg["dataset"]["source_dem_pixel_spacing_m"]) / scale}
        patch_dem = cubic_crop(coefficients, (scaled_rows, scaled_cols), bounds)
        patch_dem -= np.nanmin(patch_dem)
        clean_unwrapped = [patch_dem * ratio for ratio in ratios]
        clean_wrapped = [wrap_phase(item) for item in clean_unwrapped]
        wraps = [np.rint((u - w) / (2 * math.pi)) for u, w in zip(clean_unwrapped, clean_wrapped)]
        max_wrap = int(max(np.max(item) for item in wraps))
        base_plan["max_wrap_count"] = max_wrap
        if max_wrap > int(cfg["dataset"]["generation"]["max_wrap_count"]):
            base_plan["reject_reason"] = "wrap"
            plan.append(base_plan)
            continue
        source_rows = nearest_source_indices(r0, bounds[2], rows, scaled_rows)
        source_cols = nearest_source_indices(c0, bounds[3], cols, scaled_cols)
        patch_lc = landcover[np.ix_(source_rows, source_cols)]
        if num_map is not None:
            num_crop = num_map[np.ix_(source_rows, source_cols)].astype(np.float64)
            base_plan["num_void_fraction"] = float(np.mean(num_crop == 0))
            base_plan["num_mean"] = float(num_crop.mean())
            base_plan["num_min"] = float(num_crop.min())
            if base_plan["num_void_fraction"] > float(cfg["dataset"]["generation"]["void_fraction_threshold"]):
                base_plan["reject_reason"] = "void"
                plan.append(base_plan)
                continue
        coherence_terrain, terrain = terrain_coherence(
            patch_dem, baselines, cfg, base_plan["ground_pixel_spacing_m"],
            np.random.Generator(np.random.PCG64DXSM(seeds["coherence"])))
        coherence_scene, landcover_info = apply_landcover(coherence_terrain, patch_lc, cfg)
        terrain_summary = summarize_terrain(terrain)
        terrain_cfg = cfg["dataset"]["terrain_sampling"]
        terrain_probability = 1.0 if terrain_cfg["mode"] == "uniform" else float(
            terrain_cfg["acceptance_probability"][terrain_summary["class_name"]])
        lc_multiplier, priority_code = _landcover_multiplier(patch_lc, lc_cfg["sampling"])
        accept_probability = min(1.0, terrain_probability * lc_multiplier)
        base_plan.update({"terrain_class": terrain_summary["class_name"],
                          "dominant_landcover_code": landcover_info["dominant_code"],
                          "landcover_sampling_multiplier": lc_multiplier,
                          "priority_landcover_code": priority_code,
                          "accept_probability": accept_probability})
        if np.random.Generator(np.random.PCG64DXSM(seeds["accept"])).random() > accept_probability:
            base_plan["reject_reason"] = "terrain_weighted"
            plan.append(base_plan)
            continue
        accepted_count += 1
        global_id = offset + accepted_count
        patch_name = f"{stem}_patch_{global_id:05d}"
        relative = Path("patch_groups") / row["split"] / f"{patch_name}.mat"
        snr_rng = generator(batch_seed, *base, "node_snr")
        snr_min, snr_max = cfg["noise"]["node_snr_db_range"]
        snr = float(snr_min) + (float(snr_max) - float(snr_min)) * snr_rng.random(
            int(cfg["interferometry"]["num_uavs"]))
        noisy_only, nodes, coherence_true = simulate_node_slc(
            clean_wrapped, coherence_scene, edge_index, int(cfg["interferometry"]["num_uavs"]),
            snr, np.random.Generator(np.random.PCG64DXSM(seeds["slc"])))
        if cfg["phase4_errors"]["enabled"]:
            nodes, noisy, coreg_valid, phase4 = apply_phase4(
                nodes, edge_index, cfg, base_plan["ground_pixel_spacing_m"],
                np.random.Generator(np.random.PCG64DXSM(seeds["phase4_error"])))
        else:
            noisy, coreg_valid = noisy_only, np.ones((len(baselines), patch_size, patch_size), bool)
            phase4 = {"profile": "disabled", "node_parameters": [
                {"sync_jump_rad": 0.0} for _ in range(int(cfg["interferometry"]["num_uavs"]))]}
        phase4["seed"] = seeds["phase4_error"]
        coherence_observed = [coherence_from_slc(nodes[int(edge_index[0, k]) - 1],
                                                 nodes[int(edge_index[1, k]) - 1], 7)
                              for k in range(len(baselines))]
        multilook = multilook_phases(nodes, edge_index, int(cfg["noise"]["multilook_window"]))
        mask = sample_edge_mask(len(baselines), shortest, cfg["interferometry"]["edge_sampling"],
                                np.random.Generator(np.random.PCG64DXSM(seeds["edge_mask"])))
        mask, labels = failure_mask(mask, edge_index, baselines, coherence_true, phase4,
                                    cfg["phase4_failures"],
                                    np.random.Generator(np.random.PCG64DXSM(seeds["failure"])))
        metadata = _metadata(row, dem_path, global_id, patch_name, scale, bounds,
                             base_plan["ground_pixel_spacing_m"], terrain_summary, landcover_info,
                             batch_seed, seeds, phase4, labels, snr, base_plan)
        payload = _build_payload(cfg, clean_wrapped, noisy, clean_unwrapped,
                                 coherence_observed, coherence_true, coherence_scene,
                                 coherence_terrain, noisy_only, multilook, mask, edge_index,
                                 terrain, metadata, patch_lc, landcover_info, coreg_valid,
                                 phase4, labels, numeric_type)
        destination = output / relative
        if not (resume and destination.exists()):
            atomic_savemat(destination, payload)
        base_plan.update({"accepted_index": accepted_count, "accepted": True,
                          "patch_global_id": global_id, "patch_name": patch_name,
                          "patch_group_file": relative.as_posix(), "reject_reason": "",
                          "node_snr_db": ";".join(f"{item:.12g}" for item in snr),
                          "valid_edge_mask": "".join("1" if item else "0" for item in mask),
                          "active_edge_count": int(mask.sum()),
                          "dropout_present": labels["dropout_present"],
                          "sync_anomaly_present": labels["sync_anomaly_present"],
                          "low_coherence_present": labels["low_coherence_present"],
                          "compound_failure_present": labels["compound_failure_present"],
                          "failure_factor_count": labels["failure_factor_count"]})
        plan.append(base_plan)
        if accepted_count == target:
            break
    return {"manifest_index": manifest_index, "source_file": dem_path.name,
            "target": target, "accepted": accepted_count, "attempts": len(plan), "plan": plan}


def _build_payload(cfg: dict[str, Any], clean: list[np.ndarray], noisy: list[np.ndarray],
                   unwrapped: list[np.ndarray], observed: list[np.ndarray], true: list[np.ndarray],
                   scene: list[np.ndarray], terrain_only: list[np.ndarray],
                   noisy_only: list[np.ndarray], multilook: list[np.ndarray], mask: np.ndarray,
                   edge_index: np.ndarray, terrain: dict[str, Any], metadata: dict[str, Any],
                   landcover: np.ndarray, lc_info: dict[str, Any], coreg_valid: np.ndarray,
                   phase4: dict[str, Any], labels: dict[str, Any], dtype: np.dtype) -> dict[str, Any]:
    terrain_out = {key: (np.asarray(value, dtype=dtype) if isinstance(value, np.ndarray) and
                         np.issubdtype(value.dtype, np.number) else value)
                   for key, value in terrain.items()}
    coreg_valid = coreg_valid.astype(np.uint8)
    coreg_valid[~mask] = 0
    return {"wrappedphase_withoutnoise": stack(clean, mask, dtype),
            "wrappedphase_withnoise": stack(noisy, mask, dtype),
            "unwrapped_phase": stack(unwrapped, mask, dtype),
            "coherence_estimated": stack(observed, mask, dtype),
            "coherence_observed": stack(observed, mask, dtype),
            "coherence_true": stack(true, mask, dtype),
            "coherence_scene": stack(scene, mask, dtype),
            "coherence_terrain_only": stack(terrain_only, mask, dtype),
            "wrappedphase_multilook": stack(multilook, mask, dtype),
            "wrappedphase_node_noise_only": stack(noisy_only, mask, dtype),
            "valid_edge_mask": mask.astype(np.uint8), "edge_index": edge_index,
            "baseline_perp_m": np.asarray(cfg["interferometry"]["baseline_perp_m"], dtype=dtype),
            "ambiguity_height_m": np.asarray(cfg["physics"]["ambiguity_heights_m"], dtype=dtype),
            "dem2phase_ratio_rad_per_m": np.asarray(cfg["physics"]["dem2phase_ratios_rad_per_m"], dtype=dtype),
            "phase_path_multiplicity": np.uint8(cfg["interferometry"]["phase_path_multiplicity"]),
            "terrain_features": terrain_out, "metadata": metadata,
            "schema_version": np.uint16(1), "landcover_codes": landcover.astype(np.uint8),
            "landcover_factor": np.asarray(lc_info["factor_map"], dtype=dtype),
            "coregistration_valid_mask": coreg_valid, "phase4_errors": phase4,
            "failure_labels": labels}


def _metadata(row: dict[str, Any], path: Path, global_id: int, patch_name: str,
              scale: float, bounds: tuple[int, int, int, int], spacing: float,
              terrain: dict[str, Any], lc: dict[str, Any], batch_seed: int,
              seeds: dict[str, int], phase4: dict[str, Any], labels: dict[str, Any],
              snr: np.ndarray, plan: dict[str, Any]) -> dict[str, Any]:
    sensor = "COPERNICUS" if path.name.startswith("Copernicus") else (
        "ASTGTM" if path.name.startswith("ASTGTM") else "ALPSMLC")
    return {"patch_global_id": global_id, "patch_name": patch_name, "source_file": path.name,
            "dataset_split": row["split"], "split_group": row["geographic_tile"],
            "geographic_tile": row["geographic_tile"], "geographic_region": row["region"],
            "source_dataset": row["source"], "source_quality_role": row["quality_role"],
            "num_available": bool(plan["num_available"]),
            "num_void_fraction": plan["num_void_fraction"], "num_mean": plan["num_mean"],
            "num_min": plan["num_min"], "file_type": "DEM" if "_dem" in path.name.lower() else "DSM",
            "sensor": sensor, "tile": row["geographic_tile"], "strip": "N/A",
            "interp_scale": scale, "crop_bounds_rc": np.asarray(bounds, dtype=np.int64) + 1,
            "ground_pixel_spacing_m": spacing, "terrain_class": terrain["class_name"],
            "noise": {"model": "complex_slc_nodes", "node_snr_db": snr},
            "landcover_source": lc["source"], "landcover_profile": lc["profile"],
            "dominant_landcover_code": lc["dominant_code"], "water_fraction": lc["water_fraction"],
            "rng_batch_seed": batch_seed, "rng_component_seeds": seeds,
            "phase4_error_profile": phase4["profile"], "failure_labels": labels,
            "coherence_observation_method": "observed_node_slc_pair",
            "accept_probability": plan["accept_probability"],
            "landcover_sampling_multiplier": plan["landcover_sampling_multiplier"],
            "priority_landcover_code": plan["priority_landcover_code"],
            "generator": "dem2phase-python", "generator_version": __version__}


def _landcover_multiplier(codes: np.ndarray, sampling: dict[str, Any]) -> tuple[float, int]:
    if not sampling["enabled"]:
        return 1.0, 0
    best_multiplier, best_code = 1.0, 0
    for code, fraction, multiplier in zip(sampling["priority_codes"], sampling["minimum_fraction"],
                                           sampling["acceptance_multiplier"]):
        value = float(multiplier)
        if np.mean(codes == int(code)) >= float(fraction) and value > best_multiplier:
            best_multiplier, best_code = value, int(code)
    return best_multiplier, best_code


def _num_path(path: Path) -> Path | None:
    stem = re.sub(r"_dem$", "_num", path.stem, flags=re.IGNORECASE)
    candidate = path.with_name(stem + path.suffix)
    return candidate if candidate.exists() and candidate != path else None


def _write_logs(output: Path, results: list[dict[str, Any]], accepted: pd.DataFrame,
                seed: int) -> None:
    summary = pd.DataFrame([{key: result[key] for key in
                             ("manifest_index", "source_file", "target", "accepted", "attempts")}
                            for result in results])
    summary["rng_seed"] = seed
    summary.to_csv(output / "generation_log.csv", index=False)
    accepted.to_csv(output / "generation_detail.csv", index=False, float_format="%.12g")
    file_rows = []
    for _, row in accepted.iterrows():
        for baseline_idx in range(1, 1 + len(str(row["valid_edge_mask"]))):
            file_rows.append({"patch_global_id": row["patch_global_id"],
                              "patch_name": row["patch_name"], "baseline_idx": baseline_idx,
                              "valid_edge": int(str(row["valid_edge_mask"])[baseline_idx - 1]),
                              "patch_group_file": row["patch_group_file"]})
    pd.DataFrame(file_rows).to_csv(output / "generation_file_detail.csv", index=False)


def _write_balanced_manifest(output: Path, accepted: pd.DataFrame) -> None:
    frame = accepted.copy()
    frame["path"] = frame["patch_group_file"].map(lambda value: str((output / value).resolve()))
    frame["scenario"] = "pband_uav_nominal_v1_unvalidated"
    frame["split_group"] = frame["geographic_tile"]
    frame["terrain"] = frame["terrain_class"]
    frame["landcover_code"] = frame["dominant_landcover_code"]
    frame["active_edges"] = frame["active_edge_count"]
    frame["dropout"] = frame["dropout_present"].astype(int)
    frame["sync_anomaly"] = frame["sync_anomaly_present"].astype(int)
    frame["low_coherence"] = frame["low_coherence_present"].astype(int)
    frame["compound"] = frame["compound_failure_present"].astype(int)
    frame["stratum"] = frame.apply(lambda r: f"D{r.dropout}_S{r.sync_anomaly}_L{r.low_coherence}_K{r.active_edges}_T{r.terrain}", axis=1)
    frame["sample_weight"] = 1.0
    train = frame["split"] == "train"
    if train.any():
        counts = frame.loc[train, "stratum"].value_counts()
        raw = len(counts.index) and len(frame.loc[train]) / (len(counts.index) * frame.loc[train, "stratum"].map(counts))
        raw = raw.clip(0.1, 10.0)
        frame.loc[train, "sample_weight"] = raw / raw.mean()
    columns = ["path", "scenario", "split_group", "geographic_tile", "split", "region",
               "source_dataset", "source_quality_role", "num_available", "num_mean", "terrain",
               "landcover_code", "active_edges", "dropout", "sync_anomaly", "low_coherence",
               "compound", "failure_factor_count", "stratum", "sample_weight"]
    frame[columns].to_csv(output / "failure_balanced_manifest.csv", index=False)
    frame.groupby(["split", "stratum"]).size().rename("sample_count").reset_index().to_csv(
        output / "failure_balanced_manifest_strata.csv", index=False)

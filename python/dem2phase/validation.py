"""Dataset schema, reproducibility hashing, and MATLAB/Python comparisons."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.io import loadmat


REQUIRED_FIELDS = {
    "wrappedphase_withoutnoise", "wrappedphase_withnoise", "unwrapped_phase",
    "coherence_observed", "coherence_true", "valid_edge_mask",
    "coregistration_valid_mask", "edge_index", "baseline_perp_m",
    "ambiguity_height_m", "failure_labels", "metadata", "schema_version",
    "landcover_codes", "landcover_factor",
}
ARRAY_FIELDS = sorted(REQUIRED_FIELDS - {"failure_labels", "metadata"})


def _field_names(path: Path) -> set[str]:
    return {name for name in loadmat(path, variable_names=None).keys() if not name.startswith("__")}


def validate_dataset(root: Path) -> dict[str, Any]:
    files = sorted((root / "patch_groups").glob("*/*.mat"))
    if not files:
        raise ValueError("No grouped MAT samples found")
    split_counts: dict[str, int] = {}
    groups: dict[str, set[str]] = {}
    hashes: dict[str, str] = {}
    for path in files:
        data = loadmat(path, squeeze_me=True, struct_as_record=False)
        missing = REQUIRED_FIELDS - set(data)
        if missing:
            raise ValueError(f"{path} lacks {sorted(missing)}")
        phase = np.asarray(data["wrappedphase_withnoise"])
        coherence = np.asarray(data["coherence_observed"])
        mask = np.asarray(data["valid_edge_mask"]).astype(bool).ravel()
        if phase.ndim != 3 or phase.shape != coherence.shape or phase.shape[0] != mask.size:
            raise ValueError(f"Invalid [K,H,W] schema in {path}")
        if not mask.any() or np.any(phase[~mask] != 0) or np.any(coherence[~mask] != 0):
            raise ValueError(f"Invalid inactive-edge zero filling in {path}")
        metadata = data["metadata"]
        split = str(metadata.dataset_split).strip()
        group = str(metadata.split_group).strip()
        split_counts[split] = split_counts.get(split, 0) + 1
        groups.setdefault(group, set()).add(split)
        digest = hashlib.sha256()
        for field in ARRAY_FIELDS:
            value = np.ascontiguousarray(np.asarray(data[field]))
            digest.update(field.encode())
            digest.update(str(value.dtype).encode())
            digest.update(np.asarray(value.shape, np.int64).tobytes())
            digest.update(value.tobytes())
        hashes[path.relative_to(root).as_posix()] = digest.hexdigest()
    leaking = {group: list(splits) for group, splits in groups.items() if len(splits) > 1}
    if leaking:
        raise ValueError(f"Geographic split leakage: {leaking}")
    plan_path = root / "generation_plan.csv"
    if not plan_path.exists():
        raise ValueError("generation_plan.csv is missing")
    plan = pd.read_csv(plan_path)
    accepted = plan[plan["accepted"] == True]  # noqa: E712
    if len(accepted) != len(files) or accepted["patch_group_file"].duplicated().any():
        raise ValueError("Generation plan does not match saved patch groups")
    seed_columns = ["coherence_seed", "slc_seed", "phase4_error_seed",
                    "edge_mask_seed", "failure_seed", "accept_seed"]
    seed_values = accepted[seed_columns].astype("uint64").to_numpy().ravel()
    if np.unique(seed_values).size != seed_values.size:
        raise ValueError("Accepted samples contain duplicate component seeds")
    report = {"samples": len(files), "split_counts": split_counts,
              "geographic_groups": len(groups), "leaking_groups": 0,
              "unique_component_seeds": int(seed_values.size),
              "array_hashes": hashes}
    (root / "validation_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def compare_datasets(matlab_root: Path, python_root: Path) -> dict[str, Any]:
    matlab_files = sorted((matlab_root / "patch_groups").glob("*/*.mat"))
    python_files = sorted((python_root / "patch_groups").glob("*/*.mat"))
    matlab_by_source = _group_by_source(matlab_files)
    python_by_source = _group_by_source(python_files)
    shared = sorted(set(matlab_by_source) & set(python_by_source))
    if not shared:
        raise ValueError("Datasets have no shared source DEMs")
    rows = []
    for source in shared:
        left = _distribution(matlab_by_source[source])
        right = _distribution(python_by_source[source])
        rows.append({"source_file": source, "matlab_samples": len(matlab_by_source[source]),
                     "python_samples": len(python_by_source[source]),
                     "mean_coherence_abs_diff": abs(left["mean_coherence"] - right["mean_coherence"]),
                     "p05_coherence_abs_diff": abs(left["p05"] - right["p05"]),
                     "p50_coherence_abs_diff": abs(left["p50"] - right["p50"]),
                     "p95_coherence_abs_diff": abs(left["p95"] - right["p95"]),
                     "mean_active_edges_abs_diff": abs(left["active_edges"] - right["active_edges"])})
    frame = pd.DataFrame(rows)
    matlab_all = _distribution(matlab_files)
    python_all = _distribution(python_files)
    pooled = {
        "mean_coherence_abs_diff": abs(matlab_all["mean_coherence"] - python_all["mean_coherence"]),
        "p05_coherence_abs_diff": abs(matlab_all["p05"] - python_all["p05"]),
        "p50_coherence_abs_diff": abs(matlab_all["p50"] - python_all["p50"]),
        "p95_coherence_abs_diff": abs(matlab_all["p95"] - python_all["p95"]),
        "mean_active_edges_abs_diff": abs(matlab_all["active_edges"] - python_all["active_edges"]),
    }
    pooled_quantile_max = max(pooled[key] for key in
                              ("p05_coherence_abs_diff", "p50_coherence_abs_diff",
                               "p95_coherence_abs_diff"))
    report = {"sources_compared": len(frame),
              "pooled": pooled,
              "passes_nominal_thresholds": bool(pooled["mean_coherence_abs_diff"] <= 0.03 and
                                                  pooled_quantile_max <= 0.05),
              "per_source_diagnostics": rows,
              "note": "Per-source values are diagnostics only; small ASTER strata are not acceptance gates."}
    (python_root / "matlab_python_comparison.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    return report


def _group_by_source(paths: list[Path]) -> dict[str, list[Path]]:
    result: dict[str, list[Path]] = {}
    for path in paths:
        metadata = loadmat(path, variable_names=["metadata"], squeeze_me=True,
                           struct_as_record=False)["metadata"]
        result.setdefault(str(metadata.source_file).strip(), []).append(path)
    return result


def _distribution(paths: list[Path]) -> dict[str, float]:
    values, active = [], []
    for path in paths:
        data = loadmat(path, variable_names=["coherence_observed", "valid_edge_mask"])
        coherence = np.asarray(data["coherence_observed"])
        mask = np.asarray(data["valid_edge_mask"]).astype(bool).ravel()
        values.append(coherence[mask].ravel())
        active.append(mask.sum())
    joined = np.concatenate(values)
    p05, p50, p95 = np.quantile(joined, [0.05, 0.5, 0.95])
    return {"mean_coherence": float(joined.mean()), "p05": float(p05),
            "p50": float(p50), "p95": float(p95), "active_edges": float(np.mean(active))}

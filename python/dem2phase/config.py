"""Configuration validation and derived InSAR physics."""

from __future__ import annotations

import copy
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

C0 = 299_792_458.0


def _require(mapping: dict[str, Any], names: list[str], context: str) -> None:
    missing = [name for name in names if name not in mapping]
    if missing:
        raise ValueError(f"Missing {context} fields: {', '.join(missing)}")


def load_config(path: str | Path) -> dict[str, Any]:
    path = Path(path).resolve()
    with path.open("r", encoding="utf-8") as handle:
        cfg = json.load(handle)
    cfg = validate_and_derive(cfg)
    cfg["config_path"] = str(path)
    cfg["project_root"] = str(path.parent.parent)
    return cfg


def validate_and_derive(value: dict[str, Any]) -> dict[str, Any]:
    cfg = copy.deepcopy(value)
    _require(cfg, ["dataset", "radar", "geometry", "coherence", "noise",
                   "phase4_errors", "phase4_failures", "interferometry"], "config")
    _require(cfg["dataset"], ["dem_directory", "split_manifest", "output_directory",
                              "source_dem_pixel_spacing_m", "generation", "landcover",
                              "storage", "terrain_sampling"], "dataset")
    _require(cfg["radar"], ["center_frequency_hz", "bandwidth_hz"], "radar")
    _require(cfg["geometry"], ["altitude_m", "incidence_angle_deg", "slant_range_m"],
             "geometry")
    intr = cfg["interferometry"]
    _require(intr, ["measurement_mode", "phase_path_multiplicity", "num_uavs",
                    "edge_mode", "reference_uav", "baseline_perp_m"], "interferometry")

    frequency = float(cfg["radar"]["center_frequency_hz"])
    bandwidth = float(cfg["radar"]["bandwidth_hz"])
    altitude = float(cfg["geometry"]["altitude_m"])
    theta_deg = float(cfg["geometry"]["incidence_angle_deg"])
    slant_range = float(cfg["geometry"]["slant_range_m"])
    baselines = np.asarray(intr["baseline_perp_m"], dtype=np.float64)
    num_uavs = int(intr["num_uavs"])
    multiplicity = int(intr["phase_path_multiplicity"])
    if not (frequency > 0 and 0 < bandwidth < frequency):
        raise ValueError("Radar frequency/bandwidth is invalid")
    if not (altitude > 0 and slant_range >= altitude and 0 < theta_deg < 90):
        raise ValueError("Geometry is invalid")
    if np.any(~np.isfinite(baselines)) or np.any(baselines <= 0):
        raise ValueError("All perpendicular baselines must be positive")
    if intr["edge_mode"] != "star":
        raise NotImplementedError("Python v1 supports only a star graph")
    if cfg["noise"]["model"] != "complex_slc_nodes":
        raise NotImplementedError("Python v1 supports only complex_slc_nodes")
    if cfg["coherence"]["model"] != "terrain_aware":
        raise NotImplementedError("Python v1 supports only terrain_aware coherence")
    if cfg["dataset"]["storage"]["mode"] != "grouped_mat":
        raise NotImplementedError("Python v1 supports only grouped_mat storage")
    expected_multiplicity = 1 if intr["measurement_mode"] == "single_tx_multireceiver" else 2
    if multiplicity != expected_multiplicity:
        raise ValueError("phase_path_multiplicity does not match measurement_mode")
    if num_uavs - 1 != baselines.size:
        raise ValueError("A star graph requires num_uavs - 1 baselines")
    reference = int(intr["reference_uav"])
    secondaries = [node for node in range(1, num_uavs + 1) if node != reference]
    edge_index = np.asarray([[reference] * len(secondaries), secondaries], dtype=np.uint16)

    wavelength = C0 / frequency
    range_resolution = C0 / (2.0 * bandwidth)
    theta = math.radians(theta_deg)
    critical_baseline = wavelength * slant_range * math.tan(theta) / (2.0 * range_resolution)
    ratios = multiplicity * 2.0 * math.pi / wavelength * baselines / (
        slant_range * math.sin(theta)
    )
    ambiguity = 2.0 * math.pi / ratios
    if "target_ambiguity_height_m" in intr:
        target = np.asarray(intr["target_ambiguity_height_m"], dtype=np.float64)
        if target.shape != ambiguity.shape or np.any(np.abs(target - ambiguity) / target > 0.01):
            raise ValueError("Configured target ambiguity heights differ by more than 1%")
    cfg["physics"] = {
        "speed_of_light_mps": C0,
        "wavelength_m": wavelength,
        "range_resolution_m": range_resolution,
        "critical_baseline_m": critical_baseline,
        "dem2phase_ratios_rad_per_m": ratios.tolist(),
        "ambiguity_heights_m": ambiguity.tolist(),
    }
    intr["baseline_perp_m"] = baselines.tolist()
    intr["num_edges"] = int(baselines.size)
    intr["edge_index"] = edge_index.tolist()
    return cfg


def resolve_path(cfg: dict[str, Any], configured: str | Path) -> Path:
    path = Path(configured)
    if path.is_absolute():
        return path
    return Path(cfg["project_root"]) / path

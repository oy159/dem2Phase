"""Numerical core matching the active MATLAB dem2phase generation path."""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy import ndimage, signal


def wrap_phase(phase: np.ndarray) -> np.ndarray:
    return np.angle(np.exp(1j * phase))


def box_mean(data: np.ndarray, window: int) -> np.ndarray:
    kernel = np.ones((window, window), dtype=np.float64) / float(window * window)
    return signal.convolve2d(data, kernel, mode="same", boundary="fill", fillvalue=0)


def matlab_del2(data: np.ndarray, spacing: float) -> np.ndarray:
    """MATLAB-compatible 2-D ``del2`` including extrapolated boundaries."""
    data = np.asarray(data, dtype=np.float64)
    if min(data.shape) < 3:
        raise ValueError("del2 requires at least three samples per dimension")
    d2x = np.empty_like(data)
    d2x[:, 1:-1] = data[:, 2:] - 2 * data[:, 1:-1] + data[:, :-2]
    d2x[:, 0] = d2x[:, 1]
    d2x[:, -1] = d2x[:, -2]
    d2y = np.empty_like(data)
    d2y[1:-1, :] = data[2:, :] - 2 * data[1:-1, :] + data[:-2, :]
    d2y[0, :] = d2y[1, :]
    d2y[-1, :] = d2y[-2, :]
    return (d2x + d2y) / (4 * float(spacing) ** 2)


def _smooth_standard_field(shape: tuple[int, int], sigma: float,
                           rng: np.random.Generator) -> np.ndarray:
    sigma = max(float(sigma), 0.5)
    radius = int(math.ceil(3 * sigma))
    x = np.arange(-radius, radius + 1, dtype=np.float64)
    h = np.exp(-(x * x) / (2 * sigma * sigma))
    kernel = np.outer(h, h)
    kernel /= kernel.sum()
    field = signal.convolve2d(rng.standard_normal(shape), kernel, mode="same",
                              boundary="fill", fillvalue=0)
    field -= field.mean()
    scale = field.std(ddof=1)
    return field / scale if scale > np.finfo(np.float64).eps else field


def terrain_coherence(dem_patch: np.ndarray, baselines: np.ndarray,
                      cfg: dict[str, Any], pixel_spacing_m: float,
                      rng: np.random.Generator) -> tuple[list[np.ndarray], dict[str, Any]]:
    dem = np.asarray(dem_patch, dtype=np.float64).copy()
    invalid = ~np.isfinite(dem)
    if invalid.all():
        raise ValueError("DEM patch contains no finite values")
    if invalid.any():
        dem[invalid] = np.median(dem[~invalid])
    spacing = float(pixel_spacing_m)
    dz_dy, dz_dx = np.gradient(dem, spacing, spacing)
    slope_rad = np.arctan(np.hypot(dz_dx, dz_dy))
    slope_deg = np.degrees(slope_rad)
    aspect = np.arctan2(dz_dy, dz_dx)
    look_az = math.radians(float(cfg["geometry"].get("look_azimuth_deg", 90)))
    theta = math.radians(float(cfg["geometry"]["incidence_angle_deg"]))
    directional = dz_dx * math.cos(look_az) + dz_dy * math.sin(look_az)
    range_slope = np.arctan(directional)
    norm = np.sqrt(1 + dz_dx * dz_dx + dz_dy * dz_dy)
    normal_x, normal_y, normal_z = -dz_dx / norm, -dz_dy / norm, 1 / norm
    look_x = math.sin(theta) * math.cos(look_az)
    look_y = math.sin(theta) * math.sin(look_az)
    look_z = math.cos(theta)
    cos_local = np.clip(normal_x * look_x + normal_y * look_y + normal_z * look_z, -1, 1)
    local_incidence = np.degrees(np.arccos(cos_local))
    local_mean = box_mean(dem, 7)
    roughness = np.sqrt(np.maximum(box_mean(dem * dem, 7) - local_mean * local_mean, 0))
    curvature = matlab_del2(dem, spacing)
    tpi = dem - box_mean(dem, 15)
    coh_cfg = cfg["coherence"]
    if coh_cfg["terrain_normalization_model"] == "absolute_scales":
        rough_unit = np.clip(roughness / float(coh_cfg["roughness_scale_m"]), 0, 1)
        curv_unit = np.clip(np.abs(curvature) / float(coh_cfg["curvature_scale_per_m"]), 0, 1)
    else:
        rough_unit = _robust_unit(roughness)
        curv_unit = _robust_unit(np.abs(curvature))
    slope_term = (slope_deg / float(coh_cfg["slope_scale_deg"])) ** 2
    terrain_quality = np.exp(-float(coh_cfg["slope_weight"]) * slope_term
                             - float(coh_cfg["roughness_weight"]) * rough_unit
                             - float(coh_cfg["curvature_weight"]) * curv_unit)
    terrain_quality *= np.sqrt(np.maximum(cos_local, 0))
    layover = range_slope > theta
    shadow = (range_slope < -(math.pi / 2 - theta)) | (cos_local <= 0)
    invalid_geometry = layover | shadow | invalid
    terrain_quality[invalid_geometry] = 0
    terrain_quality = np.clip(terrain_quality, 0, 1)
    shared = _smooth_standard_field(dem.shape, float(coh_cfg["residual_scale_px"]), rng)
    maps: list[np.ndarray] = []
    for baseline in np.asarray(baselines, dtype=np.float64):
        if coh_cfg["baseline_decay_model"] == "critical_baseline":
            decay = max(0.0, 1.0 - abs(baseline) / float(cfg["physics"]["critical_baseline_m"]))
        else:
            ref = float(np.min(baselines))
            decay = math.exp(-float(coh_cfg["decay_alpha"]) * max(baseline - ref, 0) / ref)
        edge = _smooth_standard_field(dem.shape, float(coh_cfg["residual_scale_px"]) / 2, rng)
        residual = np.exp(float(coh_cfg["shared_residual_std"]) * shared
                          + float(coh_cfg["edge_residual_std"]) * edge)
        quality = np.clip(terrain_quality * residual, 0, 1)
        minimum, maximum = float(coh_cfg["min"]), float(coh_cfg["max"])
        upper = minimum + (maximum - minimum) * decay
        coherence = minimum + quality * (upper - minimum)
        coherence[invalid_geometry] = minimum
        maps.append(np.clip(coherence, 1e-4, 0.9999))
    terrain = {
        "slope_deg": slope_deg, "aspect_rad": aspect,
        "range_slope_deg": np.degrees(range_slope),
        "local_incidence_deg": local_incidence, "roughness_m": roughness,
        "curvature": curvature, "tpi_m": tpi, "terrain_quality": terrain_quality,
        "layover_mask": layover, "shadow_mask": shadow,
        "invalid_geometry_mask": invalid_geometry,
        "pixel_spacing_m": spacing,
        "look_azimuth_deg": float(cfg["geometry"].get("look_azimuth_deg", 90)),
        "baseline_decay_model": coh_cfg["baseline_decay_model"],
        "critical_baseline_m": float(cfg["physics"]["critical_baseline_m"]),
        "normalization_model": coh_cfg["terrain_normalization_model"],
        "slope_scale_deg": float(coh_cfg["slope_scale_deg"]),
        "roughness_scale_m": float(coh_cfg["roughness_scale_m"]),
        "curvature_scale_per_m": float(coh_cfg["curvature_scale_per_m"]),
    }
    return maps, terrain


def _robust_unit(data: np.ndarray) -> np.ndarray:
    values = np.sort(data[np.isfinite(data)])
    if values.size == 0:
        return np.zeros_like(data)
    lo = values[max(0, round(0.02 * values.size) - 1)]
    hi = values[max(0, round(0.98 * values.size) - 1)]
    return np.zeros_like(data) if hi <= lo else np.clip((data - lo) / (hi - lo), 0, 1)


def summarize_terrain(terrain: dict[str, Any]) -> dict[str, Any]:
    slopes = np.sort(np.asarray(terrain["slope_deg"])[np.isfinite(terrain["slope_deg"])])
    p95 = float(slopes[max(0, math.ceil(0.95 * slopes.size) - 1)]) if slopes.size else math.nan
    invalid_fraction = float(np.mean(terrain["invalid_geometry_mask"]))
    if invalid_fraction >= 0.01:
        name = "geometry_hazard"
    elif p95 < 5:
        name = "flat"
    elif p95 < 20:
        name = "rolling"
    else:
        name = "steep"
    return {"class_name": name, "mean_slope_deg": float(np.mean(slopes)),
            "p95_slope_deg": p95, "layover_fraction": float(np.mean(terrain["layover_mask"])),
            "shadow_fraction": float(np.mean(terrain["shadow_mask"])),
            "mean_terrain_quality": float(np.mean(terrain["terrain_quality"]))}


def apply_landcover(coherence: list[np.ndarray], codes: np.ndarray,
                    cfg: dict[str, Any]) -> tuple[list[np.ndarray], dict[str, Any]]:
    lc = cfg["dataset"]["landcover"]
    factor = np.full(codes.shape, float(lc["unknown_factor"]), dtype=np.float64)
    for code, value in zip(lc["class_codes"], lc["coherence_factors"]):
        factor[codes == int(code)] = float(value)
    minimum = float(cfg["coherence"]["min"])
    adjusted = [np.maximum(minimum, item * factor) for item in coherence]
    present, counts = np.unique(codes, return_counts=True)
    dominant = int(present[np.argmax(counts)]) if present.size else 0
    return adjusted, {"codes": codes.astype(np.uint8), "factor_map": factor,
                      "present_codes": present, "fractions": counts / codes.size,
                      "dominant_code": dominant, "water_fraction": float(np.mean(codes == 80)),
                      "source": lc["source"], "profile": lc["profile"],
                      "factors_are_calibrated": False}


def simulate_node_slc(clean: list[np.ndarray], coherence: list[np.ndarray],
                      edge_index: np.ndarray, num_uavs: int, snr_db: np.ndarray,
                      rng: np.random.Generator) -> tuple[list[np.ndarray], np.ndarray, list[np.ndarray]]:
    shape = clean[0].shape
    complex_normal = lambda: (rng.standard_normal(shape) + 1j * rng.standard_normal(shape)) / math.sqrt(2)
    common = complex_normal()
    nodes = [None] * num_uavs
    reference = int(edge_index[0, 0]) - 1
    snr_linear = np.power(10.0, np.asarray(snr_db) / 10.0)
    nodes[reference] = common + complex_normal() / math.sqrt(snr_linear[reference])
    effective: list[np.ndarray] = []
    for idx, phase in enumerate(clean):
        secondary = int(edge_index[1, idx]) - 1
        gamma = np.clip(coherence[idx], 0, 0.9999)
        secondary_clean = gamma * common * np.exp(-1j * phase) + np.sqrt(
            np.maximum(1 - gamma * gamma, 0)) * complex_normal()
        nodes[secondary] = secondary_clean + complex_normal() / math.sqrt(snr_linear[secondary])
        thermal = 1 / math.sqrt((1 + 1 / snr_linear[reference]) *
                                (1 + 1 / snr_linear[secondary]))
        effective.append(gamma * thermal)
    node_array = np.stack(nodes, axis=0)
    phases = edge_phases(node_array, edge_index)
    return phases, node_array, effective


def edge_phases(nodes: np.ndarray, edge_index: np.ndarray) -> list[np.ndarray]:
    return [np.angle(nodes[int(edge_index[0, k]) - 1] *
                     np.conj(nodes[int(edge_index[1, k]) - 1]))
            for k in range(edge_index.shape[1])]


def apply_phase4(nodes: np.ndarray, edge_index: np.ndarray, cfg: dict[str, Any],
                 spacing: float, rng: np.random.Generator) -> tuple[np.ndarray, list[np.ndarray], np.ndarray, dict[str, Any]]:
    errors = cfg["phase4_errors"]
    n, rows, cols = nodes.shape
    row_axis = np.linspace(-0.5, 0.5, rows)[:, None]
    row_unit = np.linspace(0, 1, rows)[:, None]
    sync_phase = np.zeros((n, rows), np.float32)
    trajectory_phase = np.zeros((n, rows), np.float32)
    los_error = np.zeros((n, rows), np.float32)
    range_disp = np.zeros((n, rows), np.float32)
    azimuth_disp = np.zeros((n, rows), np.float32)
    valid = np.zeros((n, rows, cols), bool)
    output = np.zeros_like(nodes)
    params: list[dict[str, Any]] = []
    yy, xx = np.meshgrid(np.arange(rows), np.arange(cols), indexing="ij")
    sync, trajectory = errors["synchronization"], errors["trajectory"]
    attitude, coreg = errors["attitude"], errors["coregistration"]
    phase_scale = int(cfg["interferometry"]["phase_path_multiplicity"]) * 2 * math.pi / float(cfg["physics"]["wavelength_m"])
    for node in range(n):
        bias = float(sync["phase_bias_std_rad"]) * rng.standard_normal()
        drift = float(sync["linear_drift_std_rad"]) * rng.standard_normal()
        walk = np.cumsum(float(sync["random_walk_std_rad_per_row"]) * rng.standard_normal(rows))
        walk -= walk.mean()
        jump = np.zeros(rows)
        jump_row, jump_amp = 0, 0.0
        if rng.random() < float(sync["jump_probability_per_node"]):
            jump_row = int(rng.integers(2, rows + 1))
            jump_amp = float(sync["jump_std_rad"]) * rng.standard_normal()
            jump[jump_row - 1:] = jump_amp
        sync_vector = bias + drift * row_axis[:, 0] + walk + jump
        los_bias = float(trajectory["los_bias_std_m"]) * rng.standard_normal()
        los_drift = float(trajectory["los_drift_std_m"]) * rng.standard_normal()
        vibration_amp = float(trajectory["vibration_std_m"]) * rng.standard_normal()
        cycle_range = trajectory["vibration_cycles_range"]
        cycles = float(cycle_range[0]) + (float(cycle_range[1]) - float(cycle_range[0])) * rng.random()
        vibration_phase = 2 * math.pi * rng.random()
        los = los_bias + los_drift * row_axis[:, 0] + vibration_amp * np.sin(
            2 * math.pi * cycles * row_unit[:, 0] + vibration_phase)
        trajectory_vector = phase_scale * los
        roll = float(attitude["roll_std_deg"]) * rng.standard_normal()
        pitch = float(attitude["pitch_std_deg"]) * rng.standard_normal()
        range_shift = float(coreg["range_shift_std_px"]) * rng.standard_normal() + float(
            cfg["geometry"]["altitude_m"]) * math.tan(math.radians(roll)) / spacing
        azimuth_shift = float(coreg["azimuth_shift_std_px"]) * rng.standard_normal() + float(
            cfg["geometry"]["altitude_m"]) * math.tan(math.radians(pitch)) / spacing
        range_drift = float(coreg["linear_drift_std_px"]) * rng.standard_normal()
        azimuth_drift = float(coreg["linear_drift_std_px"]) * rng.standard_normal()
        range_map = range_shift + range_drift * row_axis[:, 0]
        azimuth_map = azimuth_shift + azimuth_drift * row_axis[:, 0]
        phased = nodes[node] * np.exp(1j * (sync_vector + trajectory_vector)[:, None])
        query_x = xx - range_map[:, None]
        query_y = yy - azimuth_map[:, None]
        coords = np.asarray([query_y, query_x])
        output[node] = ndimage.map_coordinates(phased.real, coords, order=1, mode="constant",
                                                cval=0, prefilter=False) + 1j * ndimage.map_coordinates(
            phased.imag, coords, order=1, mode="constant", cval=0, prefilter=False)
        valid[node] = (query_x >= 0) & (query_x <= cols - 1) & (query_y >= 0) & (query_y <= rows - 1)
        sync_phase[node], trajectory_phase[node], los_error[node] = sync_vector, trajectory_vector, los
        range_disp[node], azimuth_disp[node] = range_map, azimuth_map
        params.append({"sync_bias_rad": bias, "sync_drift_rad": drift,
                       "sync_jump_row": jump_row, "sync_jump_rad": jump_amp,
                       "los_bias_m": los_bias, "los_drift_m": los_drift,
                       "vibration_amplitude_m": vibration_amp, "vibration_cycles": cycles,
                       "roll_deg": roll, "pitch_deg": pitch,
                       "range_shift_px": range_shift, "azimuth_shift_px": azimuth_shift})
    edge_valid = np.stack([valid[int(edge_index[0, k]) - 1] & valid[int(edge_index[1, k]) - 1]
                           for k in range(edge_index.shape[1])])
    info = {"profile": errors["profile"], "sync_phase_error_rad": sync_phase,
            "trajectory_phase_error_rad": trajectory_phase, "los_range_error_m": los_error,
            "range_displacement_px": range_disp, "azimuth_displacement_px": azimuth_disp,
            "node_parameters": params, "profile_axis": "azimuth_row"}
    return output, edge_phases(output, edge_index), edge_valid, info


def coherence_from_slc(a: np.ndarray, b: np.ndarray, window: int) -> np.ndarray:
    kernel = np.ones((window, window), np.float64) / (window * window)
    cross = signal.convolve2d(a * np.conj(b), kernel, mode="same", boundary="fill")
    p1 = signal.convolve2d(np.abs(a) ** 2, kernel, mode="same", boundary="fill")
    p2 = signal.convolve2d(np.abs(b) ** 2, kernel, mode="same", boundary="fill")
    return np.clip(np.real(np.abs(cross) / (np.sqrt(p1 * p2) + np.finfo(float).eps)), 0, 1)


def multilook_phases(nodes: np.ndarray, edge_index: np.ndarray, window: int) -> list[np.ndarray]:
    kernel = np.ones((window, window), np.float64) / (window * window)
    result = []
    for k in range(edge_index.shape[1]):
        a, b = nodes[int(edge_index[0, k]) - 1], nodes[int(edge_index[1, k]) - 1]
        result.append(np.angle(signal.convolve2d(a * np.conj(b), kernel, mode="same", boundary="fill")))
    return result


def sample_edge_mask(k: int, shortest: int, options: dict[str, Any],
                     rng: np.random.Generator) -> np.ndarray:
    if not options["enabled"]:
        return np.ones(k, bool)
    count = int(rng.integers(int(options["min_active_edges"]), int(options["max_active_edges"]) + 1))
    mask = np.zeros(k, bool)
    if options["always_include_shortest"]:
        mask[shortest] = True
    candidates = np.flatnonzero(~mask)
    needed = count - int(mask.sum())
    if needed:
        mask[rng.choice(candidates, size=needed, replace=False)] = True
    return mask


def failure_mask(mask: np.ndarray, edge_index: np.ndarray, baselines: np.ndarray,
                 coherence: list[np.ndarray], phase4: dict[str, Any], cfg: dict[str, Any],
                 rng: np.random.Generator) -> tuple[np.ndarray, dict[str, Any]]:
    mask = mask.copy()
    initial = mask.astype(np.uint8)
    shortest, longest = int(np.argmin(baselines)), int(np.argmax(baselines))
    dropout = np.zeros(mask.size, np.uint8)
    uav_available = np.ones(int(edge_index.max()), np.uint8)
    if cfg["enabled"]:
        if cfg["always_include_shortest"]: mask[shortest] = True
        if cfg["always_include_longest"]: mask[longest] = True
        droppable = np.setdiff1d(np.arange(mask.size), [shortest, longest])
        if cfg["force_secondary_uav_dropout"] and droppable.size:
            dropped_edge = int(rng.choice(droppable))
            dropped_uav = int(edge_index[1, dropped_edge])
            incident = np.any(edge_index == dropped_uav, axis=0)
            mask[incident], dropout[incident], uav_available[dropped_uav - 1] = False, 1, 0
    observed = np.concatenate([coherence[k].ravel() for k in np.flatnonzero(mask)])
    low_fraction = float(np.mean(observed < float(cfg["low_coherence_threshold"]))) if observed.size else 0.0
    jump_values = np.asarray([abs(node["sync_jump_rad"]) for node in phase4["node_parameters"]])
    jump_mask = (jump_values >= float(cfg["sync_jump_min_abs_rad"])).astype(np.uint8)
    labels = {"profile": cfg["profile"], "enabled": bool(cfg["enabled"]),
              "initial_edge_mask": initial, "edge_dropout_mask": dropout,
              "uav_available_mask": uav_available,
              "low_coherence_threshold": float(cfg["low_coherence_threshold"]),
              "low_coherence_min_fraction": float(cfg["low_coherence_min_fraction"]),
              "sync_jump_min_abs_rad": float(cfg["sync_jump_min_abs_rad"]),
              "low_coherence_fraction": low_fraction,
              "low_coherence_present": low_fraction >= float(cfg["low_coherence_min_fraction"]),
              "sync_jump_node_mask": jump_mask, "sync_anomaly_present": bool(jump_mask.any()),
              "dropout_present": bool(dropout.any()), "shortest_baseline_active": bool(mask[shortest]),
              "longest_baseline_active": bool(mask[longest])}
    labels["compound_failure_present"] = bool(labels["dropout_present"] and
        labels["sync_anomaly_present"] and labels["low_coherence_present"] and
        labels["longest_baseline_active"])
    labels["failure_factor_count"] = float(labels["dropout_present"] +
        labels["sync_anomaly_present"] + labels["low_coherence_present"] +
        labels["longest_baseline_active"])
    labels["final_edge_mask"] = mask.astype(np.uint8)
    return mask, labels


def stack(items: list[np.ndarray], mask: np.ndarray, dtype: np.dtype) -> np.ndarray:
    result = np.stack(items).astype(dtype, copy=False)
    result[~mask] = 0
    return result

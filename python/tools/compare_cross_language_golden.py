"""Recompute and compare the deterministic MATLAB/Python golden fixture."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import ndimage, signal
from scipy.io import loadmat

from dem2phase.config import load_config
from dem2phase.core import terrain_coherence, wrap_phase
from dem2phase.io import cubic_crop, nearest_source_indices, spline_coefficients


def _error(actual: np.ndarray, expected: np.ndarray, margin: int = 0) -> dict[str, float]:
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    if margin:
        selection = (slice(margin, -margin), slice(margin, -margin))
        actual, expected = actual[selection], expected[selection]
    difference = np.abs(actual - expected)
    scale = np.maximum(np.abs(expected), np.finfo(np.float64).eps)
    return {"max_abs": float(np.max(difference)),
            "mean_abs": float(np.mean(difference)),
            "max_rel": float(np.max(difference / scale))}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fixture", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--diagnostic", action="store_true")
    args = parser.parse_args()
    ref = loadmat(args.fixture, squeeze_me=True, struct_as_record=False)
    cfg = load_config(args.config)

    source = np.asarray(ref["source_dem"], dtype=np.float64)
    scaled_shape = (int(ref["scaled_rows"]), int(ref["scaled_cols"]))
    bounds_one = np.asarray(ref["crop_bounds_1based"], dtype=int).ravel()
    bounds = tuple(int(value - 1) for value in bounds_one)
    patch = cubic_crop(spline_coefficients(source), scaled_shape, bounds)
    patch -= patch.min()
    baselines = np.asarray(cfg["interferometry"]["baseline_perp_m"], dtype=np.float64)
    ratios = np.asarray(cfg["physics"]["dem2phase_ratios_rad_per_m"], dtype=np.float64)
    unwrapped = np.stack([patch * ratio for ratio in ratios])
    wrapped = np.stack([wrap_phase(item) for item in unwrapped])

    source_rows = nearest_source_indices(bounds[0], bounds[2], source.shape[0], scaled_shape[0])
    source_cols = nearest_source_indices(bounds[1], bounds[3], source.shape[1], scaled_shape[1])
    landcover = np.asarray(ref["source_landcover"], dtype=np.uint8)[np.ix_(source_rows, source_cols)]
    deterministic_cfg = json.loads(json.dumps(cfg))
    deterministic_cfg["coherence"]["shared_residual_std"] = 0.0
    deterministic_cfg["coherence"]["edge_residual_std"] = 0.0
    coherence, terrain = terrain_coherence(
        patch, baselines, deterministic_cfg, float(ref["pixel_spacing_m"]),
        np.random.Generator(np.random.PCG64DXSM(1)),
    )

    convolution = signal.convolve2d(
        np.asarray(ref["convolution_input"]), np.asarray(ref["convolution_kernel"]),
        mode="same", boundary="fill", fillvalue=0,
    )
    base_complex = np.asarray(ref["base_complex"])
    rows, cols = base_complex.shape
    yy, xx = np.meshgrid(np.arange(rows), np.arange(cols), indexing="ij")
    range_map = np.asarray(ref["range_displacement_px"]).reshape(rows, 1)
    azimuth_map = np.asarray(ref["azimuth_displacement_px"]).reshape(rows, 1)
    query_x, query_y = xx - range_map, yy - azimuth_map
    coords = np.asarray([query_y, query_x])
    registered = ndimage.map_coordinates(base_complex.real, coords, order=1, mode="constant",
                                          cval=0, prefilter=False) + 1j * ndimage.map_coordinates(
        base_complex.imag, coords, order=1, mode="constant", cval=0, prefilter=False)
    registration_mask = ((query_x >= 0) & (query_x <= cols - 1) &
                         (query_y >= 0) & (query_y <= rows - 1))

    fields = {
        "dem_patch": _error(patch, ref["dem_patch"]),
        "unwrapped_phase": _error(unwrapped, ref["unwrapped_phase"]),
        "wrapped_phase": _error(wrapped, ref["wrapped_phase"]),
        "convolution_same": _error(convolution, ref["convolution_same"]),
        "registered_complex": _error(registered, ref["registered_complex"]),
        "coherence_deterministic": _error(np.stack(coherence), ref["coherence_deterministic"]),
    }
    for name in ("slope_deg", "aspect_rad", "range_slope_deg", "local_incidence_deg",
                 "roughness_m", "curvature", "tpi_m", "terrain_quality"):
        fields[f"terrain.{name}"] = _error(terrain[name], getattr(ref["terrain"], name))

    exact = {
        "source_row_index": np.array_equal(source_rows + 1, np.asarray(ref["source_row_index"]).ravel()),
        "source_col_index": np.array_equal(source_cols + 1, np.asarray(ref["source_col_index"]).ravel()),
        "landcover_patch": np.array_equal(landcover, np.asarray(ref["landcover_patch"])),
        "layover_mask": np.array_equal(terrain["layover_mask"], ref["terrain"].layover_mask),
        "shadow_mask": np.array_equal(terrain["shadow_mask"], ref["terrain"].shadow_mask),
        "registration_valid_mask": np.array_equal(registration_mask, ref["registration_valid_mask"]),
    }
    physics_names = {
        "wavelength_m": "wavelength_m",
        "critical_baseline_m": "critical_baseline_m",
        "range_resolution_m": "range_resolution_m",
        "ambiguity_heights_m": "ambiguity_heights_m",
        "dem2phase_ratios_rad_per_m": "phase_ratios_rad_per_m",
    }
    physics = {name: _error(np.asarray(cfg["physics"][name]), np.asarray(ref[fixture_name]))
               for name, fixture_name in physics_names.items()}
    report = {"fixture_schema_version": int(ref["fixture_schema_version"]),
              "fields": fields, "exact_fields": exact, "physics": physics}
    tolerances = {
        "dem_patch": 1.5e-3,
        "unwrapped_phase": 2.0e-4,
        "wrapped_phase": 2.0e-4,
        "convolution_same": 1.0e-12,
        "registered_complex": 1.0e-12,
        "coherence_deterministic": 3.0e-4,
        "terrain.slope_deg": 3.0e-3,
        "terrain.aspect_rad": 3.0e-3,
        "terrain.range_slope_deg": 3.5e-3,
        "terrain.local_incidence_deg": 3.5e-3,
        "terrain.roughness_m": 3.0e-4,
        "terrain.curvature": 7.0e-6,
        "terrain.tpi_m": 1.2e-3,
        "terrain.terrain_quality": 3.0e-4,
    }
    numeric_pass = all(fields[name]["max_abs"] <= tolerance
                       for name, tolerance in tolerances.items())
    physics_pass = max(item["max_abs"] for item in physics.values()) <= 1e-10
    report["atol"] = tolerances
    report["passes"] = bool(all(exact.values()) and numeric_pass and physics_pass)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    if not args.diagnostic and not report["passes"]:
        raise SystemExit("Golden comparison failed")


if __name__ == "__main__":
    main()

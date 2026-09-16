"""GeoTIFF/MATLAB I/O helpers and manifest validation."""

from __future__ import annotations

import os
import re
import tempfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import rasterio
from scipy import ndimage
from scipy.io import loadmat, savemat


TILE_PATTERN = re.compile(r"^[NS]\d{2,3}[EW]\d{3}$")


def load_split_manifest(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    required = {"dem_file", "geographic_tile", "split", "region", "quality_role",
                "patches_per_dem", "source"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Split manifest lacks: {sorted(missing)}")
    if frame["dem_file"].duplicated().any():
        raise ValueError("A DEM filename appears more than once")
    if not frame["geographic_tile"].map(lambda item: bool(TILE_PATTERN.match(str(item)))).all():
        raise ValueError("Invalid geographic_tile token")
    leaking = frame.groupby("geographic_tile")["split"].nunique()
    if (leaking > 1).any():
        raise ValueError(f"Geographic split leakage: {list(leaking[leaking > 1].index)}")
    if not frame["split"].isin(["train", "test", "validation"]).all():
        raise ValueError("split must be train, test, or validation")
    return frame


def discover_dems(root: Path, manifest: pd.DataFrame, allow_extras: bool = False) -> list[Path]:
    lookup = {path.name: path for path in root.rglob("*.tif")
              if re.search(r"_(dem|dsm)(?:\.|$)", path.name, re.IGNORECASE)}
    expected = list(manifest["dem_file"])
    missing = [name for name in expected if name not in lookup]
    extras = sorted(set(lookup) - set(expected))
    if missing or (extras and not allow_extras):
        raise ValueError(f"DEM/manifest mismatch; missing={missing}, extras={extras}")
    return [lookup[name] for name in expected]


def read_raster(path: Path, dtype: np.dtype = np.float64) -> tuple[np.ndarray, dict[str, Any]]:
    with rasterio.open(path) as dataset:
        data = dataset.read(1).astype(dtype)
        profile = dataset.profile.copy()
    return data, profile


def load_aligned_landcover(path: Path) -> np.ndarray:
    value = loadmat(path, squeeze_me=True)
    for name in ("landcover_codes", "codes", "data"):
        if name in value:
            return np.asarray(value[name], dtype=np.uint8)
    candidates = [item for key, item in value.items() if not key.startswith("__") and
                  isinstance(item, np.ndarray) and item.ndim == 2]
    if len(candidates) != 1:
        raise ValueError(f"Cannot identify land-cover codes in {path}")
    return np.asarray(candidates[0], dtype=np.uint8)


def cubic_crop(coefficients: np.ndarray, scaled_shape: tuple[int, int],
               bounds: tuple[int, int, int, int]) -> np.ndarray:
    """Evaluate a MATLAB-style endpoint grid crop from cached spline coefficients."""
    r0, c0, r1, c1 = bounds
    rows, cols = coefficients.shape
    sy, sx = scaled_shape
    y = np.linspace(0, rows - 1, sy, dtype=np.float64)[r0:r1 + 1]
    x = np.linspace(0, cols - 1, sx, dtype=np.float64)[c0:c1 + 1]
    yy, xx = np.meshgrid(y, x, indexing="ij")
    return ndimage.map_coordinates(coefficients, [yy, xx], order=3, mode="nearest",
                                    prefilter=False)


def spline_coefficients(dem: np.ndarray) -> np.ndarray:
    return ndimage.spline_filter(dem, order=3, output=np.float64, mode="nearest")


def nearest_source_indices(start: int, stop: int, source_size: int,
                           scaled_size: int) -> np.ndarray:
    values = np.arange(start, stop + 1, dtype=np.float64)
    return np.clip(np.rint(values * (source_size - 1) / max(scaled_size - 1, 1)).astype(int),
                   0, source_size - 1)


def matlab_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): matlab_safe(item) for key, item in value.items()}
    if isinstance(value, list) and value and isinstance(value[0], dict):
        field_names = list(dict.fromkeys(
            field for item in value for field in item.keys()
        ))
        struct_array = np.empty((1, len(value)),
                                dtype=[(str(field), object) for field in field_names])
        for index, item in enumerate(value):
            for field in field_names:
                struct_array[0, index][field] = matlab_safe(item.get(field, ""))
        return struct_array
    if value is None:
        return ""
    return value


def atomic_savemat(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    os.close(handle)
    try:
        savemat(temporary, matlab_safe(payload), format="5", do_compression=False,
                long_field_names=True, oned_as="row", appendmat=False)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)

"""Prepare nearest-neighbour WorldCover rasters for configured DEMs."""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import Window

from .config import resolve_path
from .io import atomic_savemat, discover_dems, load_split_manifest


def worldcover_token(latitude: float, longitude: float) -> str:
    south = 3 * math.floor(latitude / 3)
    west = 3 * math.floor(longitude / 3)
    return f"{'N' if south >= 0 else 'S'}{abs(south):02d}{'E' if west >= 0 else 'W'}{abs(west):03d}"


def align_at_dem_centres(dem_ds: rasterio.DatasetReader,
                         source_ds: rasterio.DatasetReader) -> tuple[np.ndarray, float]:
    """Match MATLAB's explicit nearest sampling at every DEM cell centre."""
    dt = dem_ds.transform
    st = source_ds.transform
    if dem_ds.crs != source_ds.crs or str(dem_ds.crs) != "EPSG:4326":
        raise ValueError("DEM and WorldCover must both use EPSG:4326")
    if any(abs(value) > 1e-15 for value in (dt.b, dt.d, st.b, st.d)):
        raise ValueError("Rotated rasters are not supported")
    if not (dt.a > 0 and dt.e < 0 and st.a > 0 and st.e < 0):
        raise ValueError("Expected eastward columns and southward rows")

    columns = np.arange(dem_ds.width, dtype=np.float64)
    rows = np.arange(dem_ds.height, dtype=np.float64)
    if dem_ds.tags().get("AREA_OR_POINT", "Area").lower() == "point":
        # Preserve MATLAB GeographicPostingsReference's operation order:
        # first posting + integer * sample spacing. This matters at exact
        # half-pixel ties in the categorical source grid.
        target_lon = (dt.c + 0.5 * dt.a) + columns * dt.a
        target_lat = (dt.f + 0.5 * dt.e) + rows * dt.e
    else:
        target_lon = dt.c + (columns + 0.5) * dt.a
        target_lat = dt.f + (rows + 0.5) * dt.e
    source_x = 0.5 + (target_lon - st.c) / st.a
    source_y = 0.5 + (st.f - target_lat) / (-st.e)
    tolerance = 2.0
    if (np.any(source_x < 0.5 - tolerance) or
            np.any(source_x > source_ds.width + 0.5 + tolerance) or
            np.any(source_y < 0.5 - tolerance) or
            np.any(source_y > source_ds.height + 0.5 + tolerance)):
        raise ValueError("WorldCover tile does not fully contain the DEM")
    clamped = (np.count_nonzero((source_x < 0.5) |
                                (source_x > source_ds.width + 0.5)) +
               np.count_nonzero((source_y < 0.5) |
                                (source_y > source_ds.height + 0.5)))
    clamped_fraction = clamped / (source_x.size + source_y.size)

    # MATLAB round is half-away-from-zero. Intrinsic coordinates here are
    # positive, so floor(x + 0.5) is its exact integer equivalent.
    source_cols = np.clip(np.floor(source_x + 0.5).astype(np.int64),
                          1, source_ds.width)
    source_rows = np.clip(np.floor(source_y + 0.5).astype(np.int64),
                          1, source_ds.height)
    row_min, row_max = int(source_rows.min()), int(source_rows.max())
    col_min, col_max = int(source_cols.min()), int(source_cols.max())
    subset = source_ds.read(
        1,
        window=Window(col_min - 1, row_min - 1,
                      col_max - col_min + 1, row_max - row_min + 1),
    )
    aligned = subset[np.ix_(source_rows - row_min, source_cols - col_min)]
    return np.asarray(aligned, dtype=np.uint8), float(clamped_fraction)


def prepare_landcover(cfg: dict[str, Any], raw_root: Path | None = None,
                      aligned_root: Path | None = None) -> pd.DataFrame:
    manifest = load_split_manifest(resolve_path(cfg, cfg["dataset"]["split_manifest"]))
    dem_paths = discover_dems(resolve_path(cfg, cfg["dataset"]["dem_directory"]), manifest)
    aligned_root = aligned_root or resolve_path(
        cfg, cfg["dataset"]["landcover"]["aligned_directory"])
    aligned_root.mkdir(parents=True, exist_ok=True)
    raw_root = raw_root or aligned_root / "raw"
    raw_lookup: dict[str, Path] = {}
    for path in raw_root.glob("*.tif"):
        match = re.search(r"_([NS]\d{2}[EW]\d{3})_Map", path.name)
        if match:
            raw_lookup[match.group(1)] = path
    records = []
    for dem_path, (_, row) in zip(dem_paths, manifest.iterrows()):
        token = worldcover_token(float(row["latitude_deg"]), float(row["longitude_deg"]))
        if token not in raw_lookup:
            raise FileNotFoundError(f"Missing raw WorldCover tile {token}")
        with rasterio.open(dem_path) as dem_ds, rasterio.open(raw_lookup[token]) as src:
            destination, clamped_fraction = align_at_dem_centres(dem_ds, src)
        output = aligned_root / f"{dem_path.stem}{cfg['dataset']['landcover']['aligned_suffix']}"
        atomic_savemat(output, {"landcover_codes": destination,
                                "source": cfg["dataset"]["landcover"]["source"],
                                "worldcover_tile": token, "dem_file": dem_path.name})
        records.append({"dem_file": dem_path.name, "split": row["split"],
                        "region": row["region"], "worldcover_tile": token,
                        "raw_file": str(raw_lookup[token]), "aligned_file": str(output),
                        "rows": destination.shape[0], "cols": destination.shape[1],
                        "unknown_fraction": float(np.mean(destination == 0)),
                        "boundary_clamped_fraction": clamped_fraction,
                        "codes_present": ";".join(map(str, np.unique(destination)))})
    frame = pd.DataFrame(records)
    frame.to_csv(aligned_root / "worldcover_alignment_manifest_python.csv", index=False)
    return frame

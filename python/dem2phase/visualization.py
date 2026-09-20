"""Static diagnostic figures for grouped dem2phase MAT patches."""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any

import numpy as np
from scipy.io import loadmat


_PATCH_ID = re.compile(r"_patch_(\d+)\.mat$", re.IGNORECASE)


def _patch_sort_key(path: Path) -> tuple[int, str]:
    match = _PATCH_ID.search(path.name)
    return (int(match.group(1)) if match else 2**63 - 1, path.name.lower())


def _resolve_patch(patch: Path | None, dataset: Path | None, split: str,
                   index: int, match: str | None) -> tuple[Path, Path]:
    if patch is not None:
        path = patch.resolve()
        root = path.parents[2] if len(path.parents) >= 3 and path.parent.parent.name == "patch_groups" else path.parent
        if not path.is_file():
            raise FileNotFoundError(path)
        return path, root
    assert dataset is not None
    root = dataset.resolve()
    files = sorted((root / "patch_groups" / split).glob("*.mat"), key=_patch_sort_key)
    if match:
        files = [item for item in files if match.lower() in item.name.lower()]
    if not files:
        raise FileNotFoundError(f"No matching MAT patches in {root / 'patch_groups' / split}")
    if index < 1 or index > len(files):
        raise IndexError(f"Patch index {index} is outside 1..{len(files)}")
    return files[index - 1], root


def _stack(payload: dict[str, Any], key: str) -> np.ndarray:
    if key not in payload:
        raise KeyError(f"Required MAT field {key!r} is missing")
    value = np.asarray(payload[key]).squeeze()
    return value[None, ...] if value.ndim == 2 else value


def _field(struct: Any, name: str, default: Any = None) -> Any:
    return getattr(struct, name, default) if struct is not None else default


def _json_scalar(value: Any) -> Any:
    array = np.asarray(value).squeeze()
    if array.size != 1:
        return array.tolist()
    item = array.item()
    return item.decode() if isinstance(item, bytes) else item


def _finite_limits(array: np.ndarray, lower: float = 1, upper: float = 99) -> tuple[float, float]:
    values = array[np.isfinite(array)]
    if values.size == 0:
        return 0.0, 1.0
    lo, hi = np.percentile(values, [lower, upper])
    if not hi > lo:
        hi = lo + 1.0
    return float(lo), float(hi)


def _circular_error(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return np.angle(np.exp(1j * (left - right)))


def visualize_patch(patch: Path | None, dataset: Path | None, split: str = "test",
                    index: int = 1, match: str | None = None,
                    edges: list[int] | None = None, output: Path | None = None,
                    dpi: int = 150, show: bool = False,
                    colormap: str = "jet") -> dict[str, Any]:
    path, root = _resolve_patch(patch, dataset, split, index, match)
    payload = loadmat(path, squeeze_me=True, struct_as_record=False)
    noisy = _stack(payload, "wrappedphase_withnoise").astype(np.float64)
    clean = _stack(payload, "wrappedphase_withoutnoise").astype(np.float64)
    unwrapped = _stack(payload, "unwrapped_phase").astype(np.float64)
    coherence = _stack(payload, "coherence_observed").astype(np.float64)
    edge_mask = np.asarray(payload["valid_edge_mask"]).astype(bool).reshape(-1)
    pixel_mask = _stack(payload, "coregistration_valid_mask").astype(bool)
    baselines = np.asarray(payload["baseline_perp_m"]).reshape(-1)
    ambiguity = np.asarray(payload["ambiguity_height_m"]).reshape(-1)
    if noisy.shape[0] != edge_mask.size:
        raise ValueError(f"Edge dimension {noisy.shape[0]} does not match mask {edge_mask.size}")
    selected = ([item - 1 for item in edges] if edges else np.flatnonzero(edge_mask).tolist())
    if not selected:
        raise ValueError("The patch has no selected active edges")
    invalid = [item + 1 for item in selected if item < 0 or item >= edge_mask.size]
    if invalid:
        raise IndexError(f"Edge slots outside 1..{edge_mask.size}: {invalid}")

    if output is None:
        output = root / "visualizations" / f"{path.stem}_overview.png"
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    terrain_output = output.with_name(f"{output.stem}_terrain{output.suffix}")
    summary_output = output.with_name(f"{output.stem}_summary.json")

    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    columns = len(selected)
    fig, axes = plt.subplots(5, columns, figsize=(4.2 * columns, 17), squeeze=False,
                             constrained_layout=True)
    metadata = payload.get("metadata")
    patch_name = str(_json_scalar(_field(metadata, "patch_name", path.stem)))
    fig.suptitle(f"{patch_name} | active edges {np.flatnonzero(edge_mask) + 1}", fontsize=14)
    edge_summaries = []
    for column, edge in enumerate(selected):
        valid = pixel_mask[edge] & edge_mask[edge]
        residual = _circular_error(noisy[edge], clean[edge])
        valid_values = valid & np.isfinite(residual)
        phase_lo, phase_hi = _finite_limits(unwrapped[edge][valid] if valid.any() else unwrapped[edge])
        panels = (
            (noisy[edge], "Noisy wrapped phase", colormap, -math.pi, math.pi),
            (clean[edge], "Clean wrapped phase", colormap, -math.pi, math.pi),
            (unwrapped[edge], "Unwrapped phase GT", colormap, phase_lo, phase_hi),
            (coherence[edge], "Observed coherence", colormap, 0.0, 1.0),
            (residual, "Circular noisy-clean error", colormap, -math.pi, math.pi),
        )
        for row, (data, title, cmap, vmin, vmax) in enumerate(panels):
            image = axes[row, column].imshow(
                data, cmap=cmap, vmin=vmin, vmax=vmax,
                interpolation="nearest", resample=False,
            )
            axes[row, column].set_title(f"E{edge + 1} | {title}")
            axes[row, column].set_xticks([])
            axes[row, column].set_yticks([])
            fig.colorbar(image, ax=axes[row, column], fraction=0.046, pad=0.02)
        edge_summaries.append({
            "edge_slot": edge + 1,
            "active": bool(edge_mask[edge]),
            "baseline_perp_m": float(baselines[edge]),
            "ambiguity_height_m": float(ambiguity[edge]),
            "valid_pixel_fraction": float(valid.mean()),
            "mean_observed_coherence": float(np.mean(coherence[edge][valid])) if valid.any() else None,
            "circular_noise_mae_rad": float(np.mean(np.abs(residual[valid_values]))) if valid_values.any() else None,
            "unwrapped_phase_min_rad": float(np.min(unwrapped[edge][valid])) if valid.any() else None,
            "unwrapped_phase_max_rad": float(np.max(unwrapped[edge][valid])) if valid.any() else None,
        })
    fig.savefig(output, dpi=dpi)

    terrain = payload.get("terrain_features")
    terrain_panels = [
        (np.asarray(payload["landcover_codes"]), "WorldCover code", colormap, None, None),
        (np.asarray(payload["landcover_factor"]), "Land-cover coherence factor", colormap, 0, 1),
        (np.asarray(_field(terrain, "slope_deg")), "Slope (deg)", colormap, None, None),
        (np.asarray(_field(terrain, "local_incidence_deg")), "Local incidence (deg)", colormap, None, None),
        (np.asarray(_field(terrain, "roughness_m")), "Roughness (m)", colormap, None, None),
        (np.asarray(_field(terrain, "tpi_m")), "TPI (m)", colormap, None, None),
        (np.asarray(_field(terrain, "terrain_quality")), "Terrain quality", colormap, 0, 1),
        (np.asarray(_field(terrain, "layover_mask")), "Layover mask", "gray_r", 0, 1),
        (np.asarray(_field(terrain, "shadow_mask")), "Shadow mask", "gray_r", 0, 1),
    ]
    terrain_fig, terrain_axes = plt.subplots(3, 3, figsize=(13, 12), constrained_layout=True)
    terrain_fig.suptitle(f"{patch_name} | terrain and land cover", fontsize=14)
    for axis, (data, title, cmap, vmin, vmax) in zip(terrain_axes.ravel(), terrain_panels):
        if vmin is None:
            vmin, vmax = _finite_limits(data)
        image = axis.imshow(
            data, cmap=cmap, vmin=vmin, vmax=vmax,
            interpolation="nearest", resample=False,
        )
        axis.set_title(title)
        axis.set_xticks([])
        axis.set_yticks([])
        terrain_fig.colorbar(image, ax=axis, fraction=0.046, pad=0.02)
    terrain_fig.savefig(terrain_output, dpi=dpi)

    summary = {
        "patch_file": str(path),
        "patch_name": patch_name,
        "dataset_split": str(_json_scalar(_field(metadata, "dataset_split", split))),
        "geographic_tile": str(_json_scalar(_field(metadata, "geographic_tile", "unknown"))),
        "source_file": str(_json_scalar(_field(metadata, "source_file", "unknown"))),
        "active_edge_slots": (np.flatnonzero(edge_mask) + 1).tolist(),
        "rendered_edge_slots": [item + 1 for item in selected],
        "colormap": colormap,
        "edges": edge_summaries,
        "overview_png": str(output),
        "terrain_png": str(terrain_output),
    }
    summary_output.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    summary["summary_json"] = str(summary_output)
    if show:
        plt.show()
    plt.close(fig)
    plt.close(terrain_fig)
    return summary


def browse_patches(dataset: Path, split: str = "test", index: int = 1,
                   match: str | None = None, colormap: str = "jet") -> dict[str, Any]:
    """Open a no-output A/D keyboard browser ordered by global patch ID."""
    root = dataset.resolve()
    files = sorted((root / "patch_groups" / split).glob("*.mat"), key=_patch_sort_key)
    if match:
        files = [item for item in files if match.lower() in item.name.lower()]
    if not files:
        raise FileNotFoundError(f"No matching MAT patches in {root / 'patch_groups' / split}")
    if index < 1 or index > len(files):
        raise IndexError(f"Initial index {index} is outside 1..{len(files)}")

    import matplotlib.pyplot as plt

    state = {"index": index - 1}
    # ``compressed`` specifically removes excess whitespace around fixed-aspect
    # image axes, which is important for the 5-row diagnostic matrix.
    fig = plt.figure(figsize=(15, 13.5), layout="compressed")

    def render() -> None:
        path = files[state["index"]]
        payload = loadmat(path, squeeze_me=True, struct_as_record=False)
        noisy = _stack(payload, "wrappedphase_withnoise").astype(np.float64)
        clean = _stack(payload, "wrappedphase_withoutnoise").astype(np.float64)
        unwrapped = _stack(payload, "unwrapped_phase").astype(np.float64)
        coherence = _stack(payload, "coherence_observed").astype(np.float64)
        edge_mask = np.asarray(payload["valid_edge_mask"]).astype(bool).reshape(-1)
        active = np.flatnonzero(edge_mask).tolist()
        terrain = payload.get("terrain_features")
        metadata = payload.get("metadata")
        patch_name = str(_json_scalar(_field(metadata, "patch_name", path.stem)))
        tile = str(_json_scalar(_field(metadata, "geographic_tile", "unknown")))
        baselines = np.asarray(payload["baseline_perp_m"]).reshape(-1)
        ambiguity = np.asarray(payload["ambiguity_height_m"]).reshape(-1)

        fig.clear()
        columns = max(1, len(active)) + 1
        fig.set_size_inches(max(9.0, 3.0 * columns), 13.5, forward=True)
        axes = fig.subplots(5, columns, squeeze=False)
        fig.suptitle(
            f"[{state['index'] + 1}/{len(files)}] {patch_name} | tile {tile} | "
            "A: previous   D: next   Q/Esc: quit",
            fontsize=13,
        )
        row_names = ("Noisy wrapped", "Clean wrapped", "Unwrapped GT",
                     "Observed coherence", "Circular noise error")
        for column, edge in enumerate(active):
            residual = _circular_error(noisy[edge], clean[edge])
            lo, hi = _finite_limits(unwrapped[edge])
            panels = (
                (noisy[edge], -math.pi, math.pi),
                (clean[edge], -math.pi, math.pi),
                (unwrapped[edge], lo, hi),
                (coherence[edge], 0.0, 1.0),
                (residual, -math.pi, math.pi),
            )
            for row, (data, vmin, vmax) in enumerate(panels):
                axes[row, column].imshow(
                    data, cmap=colormap, vmin=vmin, vmax=vmax,
                    interpolation="nearest", resample=False,
                )
                axes[row, column].set_title(
                    f"E{edge + 1} B={float(baselines[edge]):g}m "
                    f"ha={float(ambiguity[edge]):.1f}m\n{row_names[row]}", fontsize=8.5)
                axes[row, column].set_xticks([])
                axes[row, column].set_yticks([])

        terrain_panels = (
            (np.asarray(payload["landcover_codes"]), "WorldCover"),
            (np.asarray(_field(terrain, "slope_deg")), "Slope (deg)"),
            (np.asarray(_field(terrain, "local_incidence_deg")), "Local incidence"),
            (np.asarray(_field(terrain, "terrain_quality")), "Terrain quality"),
            (np.asarray(_field(terrain, "layover_mask")) +
             2 * np.asarray(_field(terrain, "shadow_mask")), "Layover=1 / shadow=2"),
        )
        terrain_column = columns - 1
        for row, (data, title) in enumerate(terrain_panels):
            cmap = "gray_r" if row == 4 else colormap
            lo, hi = (0.0, 2.0) if row == 4 else _finite_limits(data)
            axes[row, terrain_column].imshow(
                data, cmap=cmap, vmin=lo, vmax=hi,
                interpolation="nearest", resample=False,
            )
            axes[row, terrain_column].set_title(title, fontsize=9)
            axes[row, terrain_column].set_xticks([])
            axes[row, terrain_column].set_yticks([])
        fig.canvas.manager.set_window_title(f"dem2phase patch browser - {state['index'] + 1}/{len(files)}")
        fig.canvas.draw_idle()

    def on_key(event: Any) -> None:
        key = (event.key or "").lower()
        if key in {"d", "right"}:
            state["index"] = (state["index"] + 1) % len(files)
            render()
        elif key in {"a", "left"}:
            state["index"] = (state["index"] - 1) % len(files)
            render()
        elif key in {"q", "escape"}:
            plt.close(fig)

    fig.canvas.mpl_connect("key_press_event", on_key)
    render()
    plt.show()
    plt.close(fig)
    return {"dataset": str(root), "split": split, "patches": len(files),
            "initial_index": index, "colormap": colormap, "files_written": 0}

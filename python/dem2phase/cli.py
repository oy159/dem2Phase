"""Command-line interface for cloud generation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .config import load_config, resolve_path
from .generator import generate_dataset
from .landcover import prepare_landcover
from .validation import compare_datasets, validate_dataset
from .visualization import visualize_patch


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="dem2phase")
    commands = parser.add_subparsers(dest="command", required=True)
    generate = commands.add_parser("generate", help="Generate grouped MAT samples")
    generate.add_argument("--config", required=True, type=Path)
    generate.add_argument("--seed", type=int)
    generate.add_argument("--workers", type=int, default=1)
    generate.add_argument("--output", type=Path)
    generate.add_argument("--resume", action="store_true")
    generate.add_argument("--dem-file", action="append",
                          help="Generate only this manifest DEM; repeat for multiple DEMs")
    generate.add_argument("--patches-per-dem", type=int,
                          help="Override manifest patch count for smoke tests")
    generate.add_argument("--split", action="append",
                          choices=("train", "test", "validation"),
                          help="Generate only this predefined split; repeat to select several")
    landcover = commands.add_parser("prepare-landcover", help="Align raw WorldCover to DEMs")
    landcover.add_argument("--config", required=True, type=Path)
    landcover.add_argument("--raw-root", type=Path)
    landcover.add_argument("--output", type=Path,
                           help="Write aligned MAT files to a separate directory")
    validate = commands.add_parser("validate", help="Validate a generated dataset")
    validate.add_argument("--dataset", required=True, type=Path)
    compare = commands.add_parser("compare-matlab", help="Compare MATLAB/Python distributions")
    compare.add_argument("--matlab-dataset", required=True, type=Path)
    compare.add_argument("--python-dataset", required=True, type=Path)
    visualize = commands.add_parser("visualize", help="Render one grouped MAT patch")
    source = visualize.add_mutually_exclusive_group(required=True)
    source.add_argument("--patch", type=Path, help="Path to one grouped MAT file")
    source.add_argument("--dataset", type=Path, help="Generated dataset root")
    visualize.add_argument("--split", choices=("train", "test", "validation"), default="test")
    visualize.add_argument("--index", type=int, default=1,
                           help="One-based patch index after filename sorting")
    visualize.add_argument("--match", help="Select the first patch whose filename contains this text")
    visualize.add_argument("--edge", type=int, action="append",
                           help="One-based edge slot to render; repeat for several")
    visualize.add_argument("--output", type=Path, help="Output PNG path")
    visualize.add_argument("--dpi", type=int, default=150)
    visualize.add_argument("--colormap", default="jet",
                           help="Matplotlib colormap for continuous and categorical maps")
    visualize.add_argument("--show", action="store_true", help="Open an interactive matplotlib window")
    return parser


def main(argv: list[str] | None = None) -> None:
    args = _parser().parse_args(argv)
    if args.command == "generate":
        cfg = load_config(args.config)
        seed = int(args.seed if args.seed is not None else cfg["dataset"]["generation"]["rng_seed"])
        output = args.output or resolve_path(cfg, cfg["dataset"]["output_directory"] + "_python")
        result = generate_dataset(cfg, output.resolve(), seed, max(1, args.workers), args.resume,
                                  args.dem_file, args.patches_per_dem, args.split)
    elif args.command == "prepare-landcover":
        cfg = load_config(args.config)
        frame = prepare_landcover(cfg, args.raw_root, args.output)
        result = {"aligned": len(frame)}
    elif args.command == "validate":
        result = validate_dataset(args.dataset.resolve())
    elif args.command == "compare-matlab":
        result = compare_datasets(args.matlab_dataset.resolve(), args.python_dataset.resolve())
    else:
        result = visualize_patch(args.patch, args.dataset, args.split, args.index, args.match,
                                 args.edge, args.output, args.dpi, args.show, args.colormap)
    print(json.dumps(result, indent=2))

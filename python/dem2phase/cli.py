"""Command-line interface for cloud generation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .config import load_config, resolve_path
from .generator import generate_dataset
from .landcover import prepare_landcover
from .validation import compare_datasets, validate_dataset


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
    return parser


def main(argv: list[str] | None = None) -> None:
    args = _parser().parse_args(argv)
    if args.command == "generate":
        cfg = load_config(args.config)
        seed = int(args.seed if args.seed is not None else cfg["dataset"]["generation"]["rng_seed"])
        output = args.output or resolve_path(cfg, cfg["dataset"]["output_directory"] + "_python")
        result = generate_dataset(cfg, output.resolve(), seed, max(1, args.workers), args.resume,
                                  args.dem_file, args.patches_per_dem)
    elif args.command == "prepare-landcover":
        cfg = load_config(args.config)
        frame = prepare_landcover(cfg, args.raw_root, args.output)
        result = {"aligned": len(frame)}
    elif args.command == "validate":
        result = validate_dataset(args.dataset.resolve())
    else:
        result = compare_datasets(args.matlab_dataset.resolve(), args.python_dataset.resolve())
    print(json.dumps(result, indent=2))

"""Compare two directories of aligned WorldCover MAT rasters pixel by pixel."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from scipy.io import loadmat


def _codes(path: Path) -> np.ndarray:
    return np.asarray(loadmat(path)["landcover_codes"], dtype=np.uint8)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    args = parser.parse_args()

    files = sorted(args.candidate.glob("*_aligned.mat"))
    if not files:
        raise SystemExit("No aligned MAT files found in candidate directory")

    total_pixels = 0
    total_mismatches = 0
    worst_name = ""
    worst_fraction = -1.0
    for candidate in files:
        reference = args.reference / candidate.name
        if not reference.is_file():
            raise FileNotFoundError(reference)
        expected = _codes(reference)
        actual = _codes(candidate)
        if expected.shape != actual.shape:
            raise ValueError(
                f"Shape mismatch for {candidate.name}: {expected.shape} != {actual.shape}"
            )
        mismatches = int(np.count_nonzero(expected != actual))
        fraction = mismatches / expected.size
        total_pixels += expected.size
        total_mismatches += mismatches
        if fraction > worst_fraction:
            worst_name = candidate.name
            worst_fraction = fraction
        print(f"{candidate.name}: {mismatches}/{expected.size} ({fraction:.8%})")

    print(f"TOTAL: {total_mismatches}/{total_pixels} "
          f"({total_mismatches / total_pixels:.8%})")
    print(f"WORST: {worst_name} ({worst_fraction:.8%})")


if __name__ == "__main__":
    main()

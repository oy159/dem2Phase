"""Stable hierarchical random streams independent of worker scheduling."""

from __future__ import annotations

import hashlib
from typing import Any

import numpy as np


RNG_ALGORITHM = "numpy-pcg64dxsm-blake2b-v1"


def derive_seed(batch_seed: int, *parts: Any) -> int:
    """Return a stable unsigned 64-bit seed for a logical component."""
    payload = "\x1f".join([str(int(batch_seed)), *(str(p) for p in parts)])
    digest = hashlib.blake2b(
        payload.encode("utf-8"), digest_size=8, person=b"dem2phase"
    ).digest()
    return int.from_bytes(digest, "little", signed=False)


def generator(batch_seed: int, *parts: Any) -> np.random.Generator:
    return np.random.Generator(np.random.PCG64DXSM(derive_seed(batch_seed, *parts)))

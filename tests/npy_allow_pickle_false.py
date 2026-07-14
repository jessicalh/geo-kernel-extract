#!/usr/bin/env python3
"""Fail unless every named producer NPY is readable without pickle."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def main() -> int:
    if len(sys.argv) < 2:
        raise SystemExit("usage: npy_allow_pickle_false.py ARRAY.npy [...]")
    for argument in sys.argv[1:]:
        path = Path(argument)
        array = np.load(path, allow_pickle=False)
        if array.dtype.hasobject:
            raise AssertionError(f"object dtype is forbidden: {path}")
        if not isinstance(array, np.ndarray):
            raise AssertionError(f"not an ndarray: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

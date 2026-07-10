"""Forcing tests for the smoke validator's SDK trust boundary."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import types

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
VALIDATOR_PATH = REPO_ROOT / "tests" / "validate_smoke.py"
VALIDATOR_SPEC = importlib.util.spec_from_file_location(
    "nmr_extract_validate_smoke", VALIDATOR_PATH)
assert VALIDATOR_SPEC is not None
assert VALIDATOR_SPEC.loader is not None
validate_smoke = importlib.util.module_from_spec(VALIDATOR_SPEC)
VALIDATOR_SPEC.loader.exec_module(validate_smoke)


class _RingContributionsProbe:
    # The second row is the real non-degenerate near-axis geometry from the
    # committed with-DFT smoke baseline.  rho is below 0.1 A but far above the
    # producer's 1e-10 A azimuth degeneracy guard.
    rho = np.array([0.0, 0.07470494696997479], dtype=np.float64)
    cos_phi = np.array([1.0, 0.577915491669793], dtype=np.float64)
    sin_phi = np.array([0.0, 0.8160966146774914], dtype=np.float64)
    n_pairs = 2


class _ProteinProbe:
    n_atoms = 1
    protein_id = "ring-azimuth-probe"
    ring_contributions = _RingContributionsProbe()
    mopac = None


def test_sdk_validator_uses_producer_ring_axis_threshold(
        tmp_path, monkeypatch):
    sdk = types.ModuleType("nmr_extract")
    sdk.load = lambda _path: _ProteinProbe()
    monkeypatch.setitem(sys.modules, "nmr_extract", sdk)

    assert validate_smoke.validate_sdk_load(tmp_path, "probe") == 0


def test_sdk_import_failure_is_a_validation_failure(
        tmp_path, monkeypatch, capsys):
    # `None` makes Python raise ModuleNotFoundError even if another test
    # imported the editable package earlier in this pytest process.
    monkeypatch.setitem(sys.modules, "nmr_extract", None)

    assert validate_smoke.validate_sdk_load(tmp_path, "probe") == 1
    output = capsys.readouterr().out
    assert "FAIL: SDK validation unavailable" in output

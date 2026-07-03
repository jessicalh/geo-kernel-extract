"""REST coverage for the shipped experimental shielding ML runtime."""

from __future__ import annotations

import os
from pathlib import Path

import pytest


MODEL_ENV = "H5READER_EXPERIMENTAL_SHIELDING_ML_MODEL"
MANIFEST_ENV = "H5READER_EXPERIMENTAL_SHIELDING_ML_MANIFEST"
HELPER_ENV = "H5READER_EXPERIMENTAL_SHIELDING_ML_HELPER"


def _dev_runtime_present() -> bool:
    return all(Path(os.environ.get(name, "")).is_file()
               for name in (MODEL_ENV, MANIFEST_ENV, HELPER_ENV))


@pytest.mark.skipif(
    not _dev_runtime_present(),
    reason="experimental shielding ML dev runtime env vars are not present",
)
def test_experimental_shielding_ml_runtime_manifest_is_reported(rest):
    state = rest.client.get("/ui/state").json()
    ml = state["experimentalShieldingMl"]

    assert ml["available"] is True
    assert ml["runtime"] == "development"

    manifest = ml["manifest"]
    assert manifest["name"] == "Experimental Shielding ML"
    assert manifest["bundleVersion"] == "0.1"
    assert manifest["bundleDate"] == "2026-07-03"
    assert manifest["target"] == "total ORCA shielding tensor as sigma_iso plus traceless T2"
    assert manifest["training"]["fold"] == "full720_static90_traj_eval"
    assert manifest["training"]["labelVocabPolicy"] == "train_full720_only_with_UNK"

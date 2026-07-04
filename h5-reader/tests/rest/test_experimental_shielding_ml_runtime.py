"""REST coverage for the shipped experimental shielding ML runtime."""

from __future__ import annotations

import os
from pathlib import Path

import pytest


MODEL_ENV = "H5READER_EXPERIMENTAL_SHIELDING_ML_MODEL"
MANIFEST_ENV = "H5READER_EXPERIMENTAL_SHIELDING_ML_MANIFEST"
HELPER_ENV = "H5READER_EXPERIMENTAL_SHIELDING_ML_HELPER"


def _dev_runtime_present() -> bool:
    model = Path(os.environ.get(MODEL_ENV, ""))
    manifest = Path(os.environ.get(MANIFEST_ENV, ""))
    helper = Path(os.environ.get(HELPER_ENV, ""))
    if not all(path.is_file() for path in (model, manifest, helper)):
        return False
    required_near_helper = (
        "c10.dll",
        "torch.dll",
        "torch_cpu.dll",
        "torch_global_deps.dll",
        "libiomp5md.dll",
        "libiompstubs5md.dll",
        "uv.dll",
    )
    return model.with_name("model_no_mopac_no_tripeptide.ts").is_file() and all(
        (helper.parent / name).is_file() for name in required_near_helper
    )


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
    assert manifest["bundleVersion"] == "0.2"
    assert manifest["bundleDate"] == "2026-07-04"
    assert manifest["target"] == "total ORCA shielding tensor as sigma_iso plus traceless T2"
    assert manifest["training"]["fold"] == "full720_static90_traj_eval"
    assert manifest["training"]["labelVocabPolicy"] == "train_full720_only_with_UNK"

    models = {model["id"]: model for model in manifest["models"]}
    assert set(models) == {"full", "no_mopac_no_tripeptide"}
    assert models["full"]["modelFile"] == "model.ts"
    assert models["full"]["inputPreset"] == "full"
    assert models["no_mopac_no_tripeptide"]["modelFile"] == "model_no_mopac_no_tripeptide.ts"
    assert models["no_mopac_no_tripeptide"]["inputPreset"] == "no_mopac_no_tripeptide"

    assert ml["inputProfile"]["loaded"] is True
    assert ml["inputProfile"]["mopacPresent"] is True
    assert ml["selectedModel"]["id"] == "full"
    assert ml["selectedModel"]["modelFile"] == "model.ts"
    assert ml["selectedModel"]["reason"] == "mopac_features_available"

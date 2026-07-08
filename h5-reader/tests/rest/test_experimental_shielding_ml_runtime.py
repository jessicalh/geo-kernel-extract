"""REST coverage for the shipped experimental shielding ML runtime."""

from __future__ import annotations

import os
from pathlib import Path

import pytest


MODEL_ENV = "H5READER_EXPERIMENTAL_SHIELDING_ML_MODEL"
MANIFEST_ENV = "H5READER_EXPERIMENTAL_SHIELDING_ML_MANIFEST"
HELPER_ENV = "H5READER_EXPERIMENTAL_SHIELDING_ML_HELPER"
EXPECT_STALE_DEV_FALLBACK_ENV = "H5READER_EXPECT_STALE_DEV_ML_FALLBACK"


REQUIRED_RUNTIME_FILES = (
    "model.ts",
    "model_no_mopac_no_tripeptide.ts",
    "manifest.json",
    "infer.exe",
    "c10.dll",
    "torch.dll",
    "torch_cpu.dll",
    "torch_global_deps.dll",
    "libiomp5md.dll",
    "libiompstubs5md.dll",
    "uv.dll",
)


def _dev_runtime_present() -> bool:
    model = Path(os.environ.get(MODEL_ENV, ""))
    manifest = Path(os.environ.get(MANIFEST_ENV, ""))
    helper = Path(os.environ.get(HELPER_ENV, ""))
    if not all(path.is_file() for path in (model, manifest, helper)):
        return False
    return model.with_name("model_no_mopac_no_tripeptide.ts").is_file() and all(
        (helper.parent / name).is_file()
        for name in REQUIRED_RUNTIME_FILES
        if name not in {"model.ts", "model_no_mopac_no_tripeptide.ts", "manifest.json", "infer.exe"}
    )


def _installed_runtime_present() -> bool:
    binary = Path(os.environ.get("H5READER_BINARY", ""))
    if not binary.is_file():
        return False
    runtime = binary.parent / "ml" / "experimental_shielding_ml"
    return all((runtime / name).is_file() for name in REQUIRED_RUNTIME_FILES)


def _runtime_present() -> bool:
    return _dev_runtime_present() or _installed_runtime_present()


def _expect_stale_dev_fallback() -> bool:
    value = os.environ.get(EXPECT_STALE_DEV_FALLBACK_ENV, "")
    return value.lower() not in {"", "0", "false", "no"}


@pytest.mark.skipif(
    not _expect_stale_dev_fallback() and not _runtime_present(),
    reason="experimental shielding ML runtime is not present",
)
def test_experimental_shielding_ml_runtime_manifest_is_reported(rest):
    if _expect_stale_dev_fallback() and not _installed_runtime_present():
        pytest.fail("expected installed experimental shielding ML runtime next to H5READER_BINARY")

    state = rest.client.get("/ui/state").json()
    ml = state["experimentalShieldingMl"]

    assert ml["available"] is True
    assert ml["runtime"] in {"development", "installed"}
    if _expect_stale_dev_fallback():
        assert ml["runtime"] == "installed"
        assert "model.ts" in ml["developmentMissing"]

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

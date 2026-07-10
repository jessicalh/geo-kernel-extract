"""REST coverage for the shipped experimental shielding ML runtime."""

from __future__ import annotations

import os
import math
import time
from pathlib import Path

import pytest


MODEL_ENV = "H5READER_EXPERIMENTAL_SHIELDING_ML_MODEL"
MANIFEST_ENV = "H5READER_EXPERIMENTAL_SHIELDING_ML_MANIFEST"
HELPER_ENV = "H5READER_EXPERIMENTAL_SHIELDING_ML_HELPER"
EXPECT_STALE_DEV_FALLBACK_ENV = "H5READER_EXPECT_STALE_DEV_ML_FALLBACK"
EXPECT_CPU_FALLBACK_ENV = "H5READER_EXPECT_ML_CPU_FALLBACK"


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

REQUIRED_ROCM_FILES = (
    "infer.exe",
    "c10.dll",
    "c10_hip.dll",
    "torch_cpu.dll",
    "torch_hip.dll",
    "caffe2_nvrtc.dll",
    "aotriton_v2.dll",
    "amdhip64_7.dll",
    "amd_comgr0702.dll",
    "hiprtc0702.dll",
    "MIOpen.dll",
    "rocblas.dll",
    "libhipblaslt.dll",
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
    rocm = runtime / "rocm"
    return (
        all((runtime / name).is_file() for name in REQUIRED_RUNTIME_FILES)
        and all((rocm / name).is_file() for name in REQUIRED_ROCM_FILES)
        and (rocm / "rocblas" / "library" / "TensileLibrary_lazy_gfx1151.dat").is_file()
        and (rocm / "hipblaslt" / "library").is_dir()
    )


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
    assert ml["inferenceReady"] is True
    assert ml["devicePreference"] == "auto"
    assert ml["configuredDevice"] == "rocm"
    assert ml["cpuFallbackConfigured"] is True
    assert ml["runtime"] in {"development", "installed"}
    if _expect_stale_dev_fallback():
        assert ml["runtime"] == "installed"
        assert "model.ts" in ml["developmentMissing"]

    manifest = ml["manifest"]
    assert manifest["name"] == "Experimental Shielding ML"
    assert manifest["bundleVersion"] == "0.2"
    assert manifest["bundleDate"] == "2026-07-04"
    assert manifest["inferenceSchemaVersion"] == 1
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


@pytest.mark.skipif(not _runtime_present(), reason="experimental shielding ML runtime is not present")
def test_experimental_shielding_ml_produces_a_dashboard_sample(rest):
    response = rest.client.post(
        "/dashboard/metric",
        json={
            "descriptor_id": "ml:experimental_shielding_iso",
            "anchor": {"atom": 16},
            "modes": ["strip.scalar"],
        },
    )
    assert response.status_code == 200, response.text
    signal_id = response.json()["id"]
    try:
        deadline = time.monotonic() + 60.0
        value = None
        while time.monotonic() < deadline:
            tracks = rest.client.get("/dashboard/display").json()["strip_tracks"]
            track = next(
                (
                    item
                    for item in tracks
                    if item["signal_id"] == signal_id
                    and item["descriptor_id"] == "ml:experimental_shielding_iso"
                ),
                None,
            )
            if track and track["valid"] and track["valid"][0] == 1:
                value = track["values"][0]
                break
            time.sleep(0.05)

        assert value is not None, "Experimental Shielding ML did not produce frame 0"
        assert math.isfinite(value)
        assert abs(value) < 10_000.0

        state = rest.client.get("/ui/state").json()
        if state["protein"].upper() == "1P9J":
            # Independent eager-Python graph/model oracle for frame 0, atom 16.
            # This pins coordinates, feature order, categorical IDs, edges, and
            # radial basis across the C++ bridge while allowing CPU/ROCm noise.
            assert abs(value - 35.75434494018555) < 0.01

        ml = state["experimentalShieldingMl"]
        assert ml["inferenceRunning"] is False
        assert ml["activeDevice"] in {"cpu", "rocm"}
        expect_cpu_fallback = os.environ.get(EXPECT_CPU_FALLBACK_ENV, "").lower() not in {
            "",
            "0",
            "false",
            "no",
        }
        if expect_cpu_fallback:
            assert ml["configuredDevice"] == "rocm"
            assert ml["activeDevice"] == "cpu"
            assert ml["usingCpuFallback"] is True
        elif ml["configuredDevice"] == "rocm" and ml["activeDevice"] == "cpu":
            assert ml["usingCpuFallback"] is True
        else:
            assert ml["usingCpuFallback"] is False
    finally:
        rest.client.post("/dashboard/metric/remove", json={"id": signal_id})

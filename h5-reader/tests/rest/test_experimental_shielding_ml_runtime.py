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
DEVICE_ENV = "H5READER_EXPERIMENTAL_SHIELDING_ML_DEVICE"
EXPECT_STALE_DEV_FALLBACK_ENV = "H5READER_EXPECT_STALE_DEV_ML_FALLBACK"
EXPECT_CPU_FALLBACK_ENV = "H5READER_EXPECT_ML_CPU_FALLBACK"
EXPECT_ACTIVE_DEVICE_ENV = "H5READER_EXPECT_ML_ACTIVE_DEVICE"
EXPECT_SUCCESS_ENV = "H5READER_EXPECT_ML_SUCCESS"
EXPECT_INPUT_FAILURE_ENV = "H5READER_EXPECT_ML_INPUT_FAILURE"


REQUIRED_RUNTIME_FILES = (
    "model.ts",
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
    return all(
        (helper.parent / name).is_file()
        for name in REQUIRED_RUNTIME_FILES
        if name not in {"model.ts", "manifest.json", "infer.exe"}
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


def _enabled(name: str) -> bool:
    return os.environ.get(name, "").lower() not in {"", "0", "false", "no"}


def _device_preference() -> str:
    preference = os.environ.get(DEVICE_ENV, "").strip().lower()
    return preference if preference in {"cpu", "rocm"} else "auto"


def _find_tree_node(nodes, field: str):
    for node in nodes:
        if node.get("field") == field:
            return node
        child = _find_tree_node(node.get("children", []), field)
        if child is not None:
            return child
    return None


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
    expected_input_failure = _enabled(EXPECT_INPUT_FAILURE_ENV)
    assert ml["inferenceReady"] is not expected_input_failure
    if expected_input_failure:
        assert (
            ml["inferenceError"]
            == "frame 0 required model input larsen_hbond_water_term.npy is absent"
        )
    preference = _device_preference()
    assert ml["devicePreference"] == preference
    assert ml["configuredDevice"] == ("cpu" if preference == "cpu" else "rocm")
    assert ml["cpuFallbackConfigured"] is (preference == "auto")
    assert ml["runtime"] in {"development", "installed"}
    if _expect_stale_dev_fallback():
        assert ml["runtime"] == "installed"
        assert "model.ts" in ml["developmentMissing"]

    manifest = ml["manifest"]
    assert manifest["name"] == "Experimental Shielding ML"
    assert manifest["bundleVersion"] == "F003-R004-v1-reader-runtime"
    assert manifest["bundleDate"] == "2026-07-19"
    assert manifest["inferenceSchemaVersion"] == 2
    assert manifest["target"] == "ORCA total shielding: isotropic 0e plus traceless 2e"

    models = {model["id"]: model for model in manifest["models"]}
    assert set(models) == {"f003_r004"}
    model = models["f003_r004"]
    assert model["modelFile"] == "model.ts"
    assert model["inputPreset"] == "f003_no_mopac_common_sense"
    assert model["training"]["bestEpoch"] == 95
    assert model["training"]["run"] == "R004-F003-no-mopac-common-sense-seed0-full96"
    assert model["training"]["producerCommit"] == "2bb3a5fa52b8a1e158d14405936f008e646ce712"

    assert ml["inputProfile"]["loaded"] is True
    assert ml["inputProfile"]["contract"] == "july_full720_f003_no_mopac_common_sense_v1"
    assert ml["selectedModel"]["id"] == "f003_r004"
    assert ml["selectedModel"]["modelFile"] == "model.ts"
    assert ml["selectedModel"]["reason"] == "july_contract_loaded"


@pytest.mark.skipif(
    not _enabled(EXPECT_SUCCESS_ENV),
    reason="fixture is not the complete F003 static acceptance member",
)
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
        assert state["protein"] == "A0A822IKI2_WT"
        # Independent eager-Python graph/model oracle for the complete July
        # full720 member. This pins coordinates, all feature blocks and masks,
        # categorical IDs, edges, and radial basis across the C++ bridge while
        # allowing CPU/ROCm floating-point noise.
        assert abs(value - 29.752159118652344) < 0.001

        ml = state["experimentalShieldingMl"]
        assert ml["inferenceRunning"] is False
        assert ml["activeDevice"] in {"cpu", "rocm"}
        expected_active_device = os.environ.get(EXPECT_ACTIVE_DEVICE_ENV, "").lower()
        if expected_active_device:
            assert ml["activeDevice"] == expected_active_device
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


@pytest.mark.skipif(
    not _enabled(EXPECT_SUCCESS_ENV),
    reason="fixture is not the complete F003 static acceptance member",
)
def test_experimental_shielding_ml_tensor_reaches_scene_and_inspector(rest):
    response = rest.client.post(
        "/dashboard/metric",
        json={
            "descriptor_id": "ml:experimental_shielding_t2",
            "anchor": {"atom": 16},
            "modes": ["static.tensor"],
        },
    )
    assert response.status_code == 200, response.text
    signal_id = response.json()["id"]

    try:
        selected = rest.client.post("/selection/pick", json={"atom": 16})
        assert selected.status_code == 204, selected.text

        deadline = time.monotonic() + 60.0
        tensor = None
        while time.monotonic() < deadline:
            state = rest.client.get("/ui/state").json()
            candidate = state["experimentalShieldingMl"]["tensorDisplay"]
            if candidate["active"] and candidate["resident"] and candidate["displayed"]:
                tensor = candidate
                break
            time.sleep(0.05)

        assert tensor is not None, "F003 tensor never reached the shared scene glyph"
        assert tensor["descriptorId"] == "ml:experimental_shielding_t2"
        assert tensor["source"] == "Experimental Shielding ML"
        assert tensor["modelId"] == "f003_r004"
        assert tensor["atom"] == 16
        assert tensor["frame"] == 0
        assert abs(tensor["sigmaIsoPpm"] - 29.752159118652344) < 0.001

        expected_t2 = (
            -4.285722732543945,
            -5.985587120056152,
            0.8754682540893555,
            8.84444522857666,
            2.2115061283111572,
        )
        assert len(tensor["t2"]) == len(expected_t2)
        for actual, expected in zip(tensor["t2"], expected_t2, strict=True):
            assert abs(actual - expected) < 0.001

        tree = rest.client.get("/inspector/tree").json()
        group = _find_tree_node(tree, "Shielding tensor (Experimental Shielding ML)")
        assert group is not None, "inspector did not identify the displayed tensor source"
        source = _find_tree_node(group.get("children", []), "Source")
        assert source is not None
        assert source["value"] == "f003_r004"
    finally:
        rest.client.post("/dashboard/metric/remove", json={"id": signal_id})

    tensor = rest.client.get("/ui/state").json()["experimentalShieldingMl"]["tensorDisplay"]
    assert tensor["active"] is False
    assert tensor["displayed"] is False


@pytest.mark.skipif(
    not _enabled(EXPECT_INPUT_FAILURE_ENV),
    reason="fixture is not expected to omit a required F003 input",
)
def test_experimental_shielding_ml_reports_missing_required_input(rest):
    response = rest.client.post(
        "/dashboard/metric",
        json={
            "descriptor_id": "ml:experimental_shielding_iso",
            "anchor": {"atom": 16},
            "modes": ["strip.scalar"],
        },
    )
    assert response.status_code == 409, response.text
    assert "not available: Absent" in response.json()["error"]

    diagnostic = rest.client.get("/ui/state").json()["diagnostic"]
    assert diagnostic["present"] is True
    assert diagnostic["source"] == "ExperimentalShieldingMlStore"
    assert (
        diagnostic["message"]
        == "frame 0 required model input larsen_hbond_water_term.npy is absent"
    )

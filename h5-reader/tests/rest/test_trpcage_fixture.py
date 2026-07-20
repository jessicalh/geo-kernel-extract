"""Trp-cage dense-DFT fixture smoke for the shipped July calcset."""

from __future__ import annotations

import math
import os

import pytest


def _fixture_is_trpcage() -> bool:
    fixture = os.environ.get("H5READER_REST_FIXTURE", "").lower()
    return "trpcage" in fixture or "1l2y" in fixture or "tc5b" in fixture


pytestmark = pytest.mark.skipif(
    not _fixture_is_trpcage(),
    reason="set H5READER_REST_FIXTURE to the Trp-cage .LGS to run",
)


def test_trpcage_dense_dft_fixture_loads_and_probes_csa(rest):
    state = rest.client.get("/ui/state").json()
    assert state["protein"] == "1L2Y_5292"
    assert state["frames"] == 1501

    atoms = rest.client.get("/protein/atoms").json()
    assert atoms["count"] == 304

    expected_sigma_iso = {
        0: 228.92000000000016,
        750: 229.8833333333333,
        1500: 225.06133333333347,
    }
    for frame, expected in expected_sigma_iso.items():
        response = rest.client.post("/frame/set", json={"frame": frame})
        assert response.status_code == 204, response.text
        csa = rest.client.get("/csa", params={"atom": 0, "frame": frame}).json()
        assert csa["atom"] == 0
        assert csa["frame"] == frame
        assert csa["dft_present"] is True
        assert csa["valid"] is True
        assert math.isfinite(csa["sigma_iso"])
        assert abs(csa["sigma_iso"] - expected) < 1e-6
        assert math.isfinite(csa["span"])
        assert math.isfinite(csa["eta"])

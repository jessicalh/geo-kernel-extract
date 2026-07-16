"""Optional chignolin/CLN025 dense-DFT fixture smoke.

This test is intentionally separate from the 1P9J-calibrated REST suite.  It
verifies that the chignolin one-off loads as a normal reader calcset and that
ORCA shielding tensors can be probed from nontrivial DFT frames.
"""

from __future__ import annotations

import math
import os

import pytest


def _fixture_is_chignolin() -> bool:
    fixture = os.environ.get("H5READER_REST_FIXTURE", "").lower()
    return "chignolin" in fixture or "cln025" in fixture or "5awl" in fixture


pytestmark = pytest.mark.skipif(
    not _fixture_is_chignolin(),
    reason="set H5READER_REST_FIXTURE to the chignolin/CLN025 .LGS to run",
)


def test_chignolin_dense_dft_fixture_loads_and_probes_csa(rest):
    atoms = rest.client.get("/protein/atoms").json()
    assert atoms["count"] == 166

    # These frame-0 / mid-trajectory values come straight from the copied
    # CLN025 ORCA outputs.  They anchor both .LGS DFT provenance and tensor
    # parsing without invoking any 1P9J-specific metric expectations.
    expected_sigma_iso = {
        0: 222.17733333333342,
        2500: 231.35400000000004,
    }
    for frame, expected in expected_sigma_iso.items():
        r = rest.client.post("/frame/set", json={"frame": frame})
        assert r.status_code == 204, r.text
        csa = rest.client.get("/csa", params={"atom": 0, "frame": frame}).json()
        assert csa["atom"] == 0
        assert csa["frame"] == frame
        assert csa["dft_present"] is True
        assert csa["valid"] is True
        assert math.isfinite(csa["sigma_iso"])
        assert abs(csa["sigma_iso"] - expected) < 1e-6
        assert math.isfinite(csa["span"])
        assert math.isfinite(csa["eta"])

    assert rest.client.post("/selection/pick", json={"atom": 0}).status_code == 204
    tree = rest.client.get("/inspector/tree").json()
    assert tree and "Atom 0" in tree[0]["field"]

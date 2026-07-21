"""REST coverage for the canonical mainstream overlay controls."""

from __future__ import annotations


CANONICAL_OVERLAYS = ("ribbon", "rings", "butterfly", "nullcone", "bfield")


def test_canonical_overlay_names_toggle_and_stale_shadow_is_not_advertised(rest):
    try:
        for name in CANONICAL_OVERLAYS:
            enabled = rest.client.post("/overlay", json={"name": name, "visible": True})
            assert enabled.status_code == 204, (name, enabled.text)

            disabled = rest.client.post(
                "/overlay", json={"name": name, "visible": False}
            )
            assert disabled.status_code == 204, (name, disabled.text)

        stale = rest.client.post("/overlay", json={"name": "shadow", "visible": True})
        assert stale.status_code == 400
        assert stale.json()["error"] == (
            'unknown overlay "shadow" (ribbon|rings|butterfly|nullcone|bfield)'
        )
    finally:
        # Restore the reader's ordinary startup presentation for later tests.
        rest.client.post("/overlay", json={"name": "ribbon", "visible": True})
        rest.client.post("/overlay", json={"name": "rings", "visible": True})
        for name in ("butterfly", "nullcone", "bfield"):
            rest.client.post("/overlay", json={"name": name, "visible": False})

"""REST coverage for transient resthero composition helpers."""

from __future__ import annotations


def test_ghost_trail_accepts_explicit_frames(rest):
    before = rest.client.get("/frame/current").json()["frame"]
    response = rest.client.post(
        "/resthero/ghost_trail",
        json={
            "atom": 16,
            "frames": [2, 0],
            "axes": "sigma33",
            "hide_selection_marker": True,
        },
    )
    assert response.status_code == 200, response.text
    payload = response.json()
    assert payload["mode"] == "frames"
    assert [g["frame"] for g in payload["ghosts"]] == [2, 0]
    assert payload["reference_frame"] == 2
    assert payload["kept"] == 2
    assert rest.client.get("/frame/current").json()["frame"] == before
    rest.client.post("/resthero/clear")


def test_legacy_butterfly_zero_extent_keeps_minimum_window(rest):
    response = rest.client.post(
        "/resthero/butterfly",
        json={
            "ring": 0,
            "frame": 0,
            "mode": "biot_savart",
            "dim": 16,
            "extent": 0.0,
        },
    )
    assert response.status_code == 200, response.text
    assert response.json()["extent"] == 0.1
    rest.client.post("/resthero/clear")

"""Adversarial REST coverage for mixed selection, filter, and resthero state.

This is intentionally a state-stack test, not a visual baseline. The user path
it protects is the poster/inspection workflow where a focused atom, residue
filter, dashboard signal, null collar, butterfly, point cloud, and ghost tensor
can all be active while frames are walked and screenshots are taken.
"""

from __future__ import annotations

from io import BytesIO


HERO_ATOM = 568
HERO_RING = 4
CROSSING_FROM = 538
CROSSING_TO = 540


def _scene_png(rest):
    response = rest.client.post("/api/screenshot", json={"target": "scene"})
    assert response.status_code == 200, response.text
    assert response.headers["content-type"] == "image/png"
    assert response.content.startswith(b"\x89PNG"), "response is not a PNG"
    return response.content


def _assert_nonblank_png(png_bytes):
    from PIL import Image

    image = Image.open(BytesIO(png_bytes))
    assert image.mode == "RGB"
    width, height = image.size
    assert width >= 200 and height >= 200, f"viewport too small: {image.size}"

    pixels = list(image.getdata())
    sampled = pixels[:: max(1, len(pixels) // 2000)]
    unique_colors = len({tuple(pixel) for pixel in sampled})
    assert unique_colors > 5, (
        f"scene appears blank or collapsed ({unique_colors} sampled colors)"
    )


def _listed_signal(rest, signal_id):
    return next(
        (
            signal
            for signal in rest.client.get("/dashboard/signals").json()
            if signal["id"] == signal_id
        ),
        None,
    )


def test_selection_filter_dashboard_and_resthero_survive_frame_walk(rest):
    rest.client.post("/resthero/clear")
    rest.client.post("/filter", json={"residues": []})
    rest.client.post("/field/null_cone", json={"visible": False})

    added_signal_id = None
    try:
        crossing = rest.client.post(
            "/api/ring/null_crossings",
            json={
                "atom": HERO_ATOM,
                "ring": HERO_RING,
                "start_frame": CROSSING_FROM,
                "end_frame": CROSSING_TO,
                "include_signal_stamps": False,
            },
        )
        assert crossing.status_code == 200, crossing.text
        crossing_payload = crossing.json()
        assert crossing_payload["summary"]["entry_count"] == 1
        entry = crossing_payload["entries"][0]
        atom_residue = entry["atom_identity"]["residue_index"]
        ring_residue = entry["ring_identity"]["parent_residue_index"]

        selected_residues = sorted({atom_residue, ring_residue})
        filter_response = rest.client.post(
            "/filter",
            json={"residues": selected_residues},
        )
        assert filter_response.status_code == 204, filter_response.text

        rest.client.post("/frame/set", json={"frame": CROSSING_FROM})
        rest.client.post(
            "/selection/pick",
            json={"atom": HERO_ATOM, "modifiers": "none"},
        )
        selection = rest.client.get("/selection").json()
        assert selection["focus"] == HERO_ATOM
        assert selection["atoms"] == [HERO_ATOM]

        signal = rest.client.post(
            "/dashboard/metric",
            json={
                "descriptor_id": "h5:sasa_time_series",
                "anchor": {"atom": HERO_ATOM},
                "modes": ["strip.scalar"],
            },
        )
        assert signal.status_code == 200, signal.text
        added_signal_id = signal.json()["id"]
        listed = _listed_signal(rest, added_signal_id)
        assert listed is not None
        assert listed["anchor"]["atom"] == HERO_ATOM

        null_cone = rest.client.post(
            "/field/null_cone",
            json={"visible": True, "opacity": 0.42, "length": 16.0},
        )
        assert null_cone.status_code == 204, null_cone.text

        style = rest.client.post(
            "/resthero/molecule_style",
            json={
                "preset": "scaffold",
                "render_atoms": False,
                "render_bonds": True,
                "bond_radius": 0.035,
                "bond_color": "#2f3439",
            },
        )
        assert style.status_code == 200, style.text
        assert style.json()["will_restore_on_clear"] is True

        ring_field = rest.client.post(
            "/resthero/ring_field",
            json={"ring": HERO_RING},
        )
        assert ring_field.status_code == 200, ring_field.text
        assert ring_field.json()["ring"] == HERO_RING

        butterfly = rest.client.post(
            "/resthero/butterfly",
            json={
                "ring": HERO_RING,
                "frame": CROSSING_TO,
                "dim": 32,
                "threshold_ppm": 0.04,
                "opacity": 0.16,
                "extent": 17.5,
                "peak": 1.0,
                "mode": "biot_savart",
            },
        )
        assert butterfly.status_code == 200, butterfly.text
        butterfly_payload = butterfly.json()
        assert butterfly_payload["ring"] == HERO_RING
        assert butterfly_payload["grid_dim"] == 32
        assert butterfly_payload["shielded_cells"] > 0
        assert butterfly_payload["deshielded_cells"] > 0

        track = rest.client.post(
            "/resthero/atom_track",
            json={
                "atom": HERO_ATOM,
                "ring": HERO_RING,
                "frame_source": "dft",
                "start_frame": 500,
                "end_frame": 580,
                "max_points": 80,
                "reference_frame": CROSSING_TO,
                "current_frame": CROSSING_TO,
                "coordinate_space": "source_ring_local",
                "color_by": "ring_current",
                "color_mode": "signed",
                "point_shape": "sphere",
                "dot_radius_A": 0.055,
                "point_opacity": 0.92,
                "show_halos": False,
                "show_lines": False,
                "color_scale": 0.2,
                "min_color_fraction": 0.18,
            },
        )
        assert track.status_code == 200, track.text
        track_payload = track.json()
        assert track_payload["atom"] == HERO_ATOM
        assert track_payload["ring"] == HERO_RING
        assert track_payload["coordinate_space"] == "source_ring_local"
        assert track_payload["vet"]["ok"] is True
        assert track_payload["kept"] > 10
        assert track_payload["min_intensity"] < 0 < track_payload["max_intensity"]

        ghost = rest.client.post(
            "/resthero/ghost_trail",
            json={
                "atom": HERO_ATOM,
                "frames": [CROSSING_TO, CROSSING_FROM],
                "axes": "sigma33",
                "hide_selection_marker": True,
            },
        )
        assert ghost.status_code == 200, ghost.text
        assert [item["frame"] for item in ghost.json()["ghosts"]] == [
            CROSSING_TO,
            CROSSING_FROM,
        ]

        for frame in (CROSSING_FROM, CROSSING_TO, CROSSING_TO + 2):
            response = rest.client.post("/frame/set", json={"frame": frame})
            assert response.status_code == 204, response.text
            assert rest.client.get("/frame/current").json()["frame"] == frame
            assert rest.client.get("/selection").json()["focus"] == HERO_ATOM
            assert _listed_signal(rest, added_signal_id) is not None

        _assert_nonblank_png(_scene_png(rest))

        cleared = rest.client.post("/resthero/clear")
        assert cleared.status_code == 204, cleared.text
        assert rest.client.get("/selection").json()["focus"] == HERO_ATOM
        assert _listed_signal(rest, added_signal_id) is not None
        _assert_nonblank_png(_scene_png(rest))
    finally:
        rest.client.post("/resthero/clear")
        rest.client.post("/filter", json={"residues": []})
        rest.client.post("/field/null_cone", json={"visible": False})
        if added_signal_id:
            rest.client.post("/dashboard/metric/remove", json={"id": added_signal_id})

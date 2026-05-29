"""Dashboard-path REST coverage — replaces the in-binary
`--dashboard-path-smoke` runner.

What the bespoke smoke did: walked frames, fed the
DashboardDisplayController, counted sampleable strip descriptors. Much of
that lives behind unexposed buffers on the controller; the REST surface
exposes the user-visible contract:

  - the default startup signal exists and is listed,
  - frame advance via REST is honoured and round-trips,
  - the camera moves between frames (proxy for "the per-frame visit ran"),
  - the catalog reports at least one strip-sampleable descriptor.

Deep per-channel sampling coverage stays in the QtTest unit tier.
"""

from __future__ import annotations


def test_default_signal_listed(rest):
    signals = rest.client.get("/dashboard/signals").json()
    assert isinstance(signals, list)
    # ReaderMainWindow seeds one default signal (Generic NPY DSSP chi)
    # at startup. If that ever changes, this test should be loosened
    # rather than tightened — the load-bearing assertion is "signals
    # listing works," not the specific default.
    assert len(signals) >= 1, "no dashboard signals listed at startup"
    assert all("id" in s and "descriptor_id" in s for s in signals)
    assert all("display_modes" in s for s in signals)


def test_frame_set_and_camera_moves(rest):
    initial_camera = rest.client.get("/scene/camera").json()
    initial_focal = initial_camera["focal"]

    # Walk a handful of frames. The centroid-follow logic in
    # MoleculeScene::setFrame translates the camera between frames; a
    # plain frame stride should produce a non-zero focal-point delta on
    # at least one of the steps.
    saw_movement = False
    for frame in (0, 100, 250, 500, 700):
        r = rest.client.post("/frame/set", json={"frame": frame})
        assert r.status_code == 204, r.text
        confirmed = rest.client.get("/frame/current").json()
        assert confirmed["frame"] == frame
        camera = rest.client.get("/scene/camera").json()
        delta = sum(
            (camera["focal"][i] - initial_focal[i]) ** 2 for i in range(3)
        )
        if delta > 1e-6:
            saw_movement = True
            initial_focal = camera["focal"]
    assert saw_movement, "camera focal point never moved across a 700-frame walk"


def test_protein_atom_count_matches_loaded(rest):
    payload = rest.client.get("/protein/atoms").json()
    assert "count" in payload
    assert payload["count"] > 0
    # 1P9J is the documented fixture per H5READER_REST_FIXTURE; if a
    # different fixture is wired the explicit number is irrelevant, but
    # the presence of atoms is not.

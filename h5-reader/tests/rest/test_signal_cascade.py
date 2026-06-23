"""Lightweight cascade verification via the REST signal listing.

The deep signal/panel cascade behaviour (clear-emits-removed-ids,
re-entrancy guard, panel-clear cleans display refs) is fully covered by
the QtTest unit tier in `tests/dashboard_model_tests.cpp`. This file
verifies that the live binary's coordinator wiring is alive: a signal
added via REST is listed, and walking frames does not silently prune it
(which would happen if the cleanup cascade fired incorrectly on innocuous
events like frame-change).

The dashboard opens EMPTY (the pruned reader seeds no default signal), so
each test adds its own signal first -- the end-to-end add → assert listed
→ frame-walk → assert survives → delete the module previously deferred
until POST /dashboard/metric existed (it now does).
"""

from __future__ import annotations


def _add_sasa(rest) -> str:
    r = rest.client.post("/dashboard/metric", json={
        "descriptor_id": "h5:sasa_time_series",
        "anchor": {"atom": 16},
        "modes": ["strip.scalar"],
    })
    assert r.status_code == 200, r.text
    return r.json()["id"]


def test_added_signal_survives_frame_walk(rest):
    sid = _add_sasa(rest)
    try:
        initial_ids = {s["id"] for s in rest.client.get("/dashboard/signals").json()}
        assert sid in initial_ids, "added signal not listed"

        for frame in (0, 100, 250, 500, 700):
            rest.client.post("/frame/set", json={"frame": frame})

        after_ids = {s["id"] for s in rest.client.get("/dashboard/signals").json()}
        assert sid in after_ids, (
            "added signal disappeared during an innocuous frame walk; "
            "this points at the cleanup cascade firing on frame-change events"
        )
    finally:
        rest.client.post("/dashboard/metric/remove", json={"id": sid})


def test_signal_listing_is_well_formed(rest):
    signals = rest.client.get("/dashboard/signals").json()
    assert isinstance(signals, list)
    for s in signals:
        for required in ("id", "label", "descriptor_id", "concept_key",
                          "display_modes", "anchor", "enabled"):
            assert required in s, f"signal entry missing {required}: {s!r}"
        assert isinstance(s["display_modes"], list)
        assert isinstance(s["anchor"], dict)
        assert "kind" in s["anchor"]

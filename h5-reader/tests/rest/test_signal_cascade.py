"""Lightweight cascade verification via the REST signal listing.

The deep signal/panel cascade behaviour (clear-emits-removed-ids,
re-entrancy guard, panel-clear cleans display refs) is fully covered by
the QtTest unit tier in `tests/dashboard_model_tests.cpp`. This file
just verifies that the live binary's coordinator wiring is alive: the
default startup signal is present, and walking frames does not silently
prune it (which would happen if the cleanup cascade fired incorrectly
on innocuous events like frame-change).

Add/delete endpoints are deliberately not exercised here — they're not
yet exposed via REST. When the add-signal endpoint lands, expand this
file with end-to-end add → assert listed → delete → assert pruned.
"""

from __future__ import annotations


def test_default_signal_survives_frame_walk(rest):
    initial = rest.client.get("/dashboard/signals").json()
    initial_ids = {s["id"] for s in initial}
    assert len(initial_ids) >= 1, "no signals listed at startup"

    for frame in (0, 100, 250, 500, 700):
        rest.client.post("/frame/set", json={"frame": frame})

    after = rest.client.get("/dashboard/signals").json()
    after_ids = {s["id"] for s in after}
    missing = initial_ids - after_ids
    assert not missing, (
        f"signal(s) disappeared during innocuous frame walk: {missing}; "
        "this points at the cascade firing on frame-change events"
    )


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

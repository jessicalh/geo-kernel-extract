"""REST coverage for the metric picker re-axed onto MetricTaxonomy.

The picker (SignalDisplayDialog) now presents the catalog as a mechanism ->
concept -> form TREE (DescriptorTreeModel) instead of a flat 188-row table.
These tests drive the headless GUI through /dashboard/picker/* to pin the
invariants the tree must preserve:

  - opening the picker on an atom auto-selects a concrete descriptor LEAF (the
    tree's depth-first first-leaf walk, preferring a displayable one), NOT a
    group / concept header -> GET /dashboard/picker reports a candidate
    descriptor id plus at least one offered display mode;
  - that auto-selected leaf is addable: /dashboard/picker/add resolves the tree
    leaf back to a SignalDescriptor and binds a signal, exactly as the old flat
    table row did.

No picker test existed before the re-axe; this is its regression net.
"""

from __future__ import annotations


def _signal_ids(rest) -> set:
    return {s["id"] for s in rest.client.get("/dashboard/signals").json()}


def test_picker_opens_with_leaf_autoselected(rest):
    """Open the picker focused on a backbone atom; the tree auto-selects a
    concrete, displayable descriptor leaf with at least one offered mode."""
    r = rest.client.post("/dashboard/picker/open", json={"atom": 16})
    assert r.status_code == 200, r.text
    state = r.json()
    assert state.get("open") is True
    # A leaf -- not a group/concept header -- is current: candidate_descriptor
    # is non-null and carries a colon-namespaced id (e.g. "h5:sasa_time_series").
    descriptor = state.get("candidate_descriptor")
    assert descriptor, f"no descriptor leaf auto-selected: {state}"
    assert ":" in descriptor, f"unexpected descriptor id shape: {descriptor!r}"
    # The auto-selected leaf is actionable: at least one display mode is offered.
    modes = state.get("modes", [])
    assert any(m.get("enabled") for m in modes), f"no enabled mode offered: {modes}"


def test_picker_add_binds_a_signal(rest):
    """The auto-selected leaf resolves back to a descriptor and binds: a
    /dashboard/picker/add yields a new entry in /dashboard/signals."""
    before = _signal_ids(rest)
    opened = rest.client.post("/dashboard/picker/open", json={"atom": 16})
    assert opened.status_code == 200, opened.text
    assert opened.json().get("candidate_descriptor"), "precondition: a leaf must be auto-selected"

    added = rest.client.post("/dashboard/picker/add")
    assert added.status_code == 200, added.text
    new_ids = _signal_ids(rest) - before
    try:
        assert new_ids, "picker add did not bind any new signal"
    finally:
        for sid in new_ids:
            rest.client.post("/dashboard/metric/remove", json={"id": sid})

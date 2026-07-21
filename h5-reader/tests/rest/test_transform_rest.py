"""REST coverage for transform state fidelity and round-tripping."""

from __future__ import annotations


SUBSET = [1, 100, 200, 300]


def test_explicit_subset_round_trips_without_becoming_backbone(rest):
    response = rest.client.post(
        "/transform",
        json={"kind": "fit_subset", "reference_frame": 0, "subset_atoms": SUBSET},
    )
    assert response.status_code == 204, response.text

    state = rest.client.get("/transform").json()
    assert state["kind"] == "fit_subset"
    assert state["subset_atoms"] == SUBSET
    assert state["subset_size"] == len(SUBSET)

    response = rest.client.post("/transform", json=state)
    assert response.status_code == 204, response.text
    round_tripped = rest.client.get("/transform").json()
    assert round_tripped["kind"] == "fit_subset"
    assert round_tripped["subset_atoms"] == SUBSET


def test_typed_backbone_subset_reports_backbone_fit(rest):
    response = rest.client.post(
        "/transform",
        json={"kind": "backbone_fit", "reference_frame": 0},
    )
    assert response.status_code == 204, response.text

    state = rest.client.get("/transform").json()
    assert state["kind"] == "backbone_fit"
    assert state["subset_size"] >= 3

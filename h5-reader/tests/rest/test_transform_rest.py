"""REST coverage for transform state fidelity and round-tripping."""

from __future__ import annotations

import numpy as np


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


def test_unframed_csa_rotates_with_displayed_molecule(rest):
    original_transform = rest.client.get("/transform").json()
    atom_count = rest.client.get("/protein/atoms").json()["count"]
    last_frame = rest.client.get("/frame/current").json()["count"] - 1
    atoms = list(range(atom_count))

    def displayed_positions():
        response = rest.client.post("/positions", json={"atoms": atoms, "frame": 0})
        assert response.status_code == 200, response.text
        return np.array([row["position"] for row in response.json()["positions"]])

    def tensor(probe):
        axes = np.array(probe["pas_axes"]).T
        return axes @ np.diag(probe["principal_values"]) @ axes.T

    try:
        response = rest.client.post(
            "/transform", json={"kind": "all_atom_fit", "reference_frame": 0}
        )
        assert response.status_code == 204, response.text
        for atom in atoms:
            response = rest.client.get("/csa", params={"atom": atom, "frame": 0})
            assert response.status_code == 200, response.text
            before = response.json()
            if before["valid"] and not before["framed"]:
                break
        else:
            raise AssertionError("Fixture has no valid unframed DFT tensor")
        positions_before = displayed_positions()

        response = rest.client.post(
            "/transform", json={"kind": "all_atom_fit", "reference_frame": last_frame}
        )
        assert response.status_code == 204, response.text
        positions_after = displayed_positions()
        response = rest.client.get("/csa", params={"atom": atom, "frame": 0})
        assert response.status_code == 200, response.text
        after = response.json()
        assert after["valid"] and not after["framed"]

        # Recover the display rotation independently from the atom coordinates.
        centered_before = positions_before - positions_before.mean(axis=0)
        centered_after = positions_after - positions_after.mean(axis=0)
        left, _, right_transpose = np.linalg.svd(centered_before.T @ centered_after)
        rotation = right_transpose.T @ left.T
        assert np.linalg.det(rotation) > 0.999999
        np.testing.assert_allclose(
            centered_before @ rotation.T, centered_after, atol=1e-8
        )
        expected = rotation @ tensor(before) @ rotation.T
        assert np.linalg.norm(expected - tensor(before)) > 1e-3
        np.testing.assert_allclose(
            after["principal_values"], before["principal_values"], atol=1e-8
        )
        np.testing.assert_allclose(tensor(after), expected, atol=1e-8)
    finally:
        response = rest.client.post("/transform", json=original_transform)
        assert response.status_code == 204, response.text

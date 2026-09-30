"""Check the displayed measurement, using independent coordinate calculations."""

from __future__ import annotations

import os
import time

import numpy as np
import pytest


def _post(rest, path, body=None):
    response = rest.client.post(path, json=body)
    assert response.is_success, response.text
    return response


def _readout(rest):
    return rest.client.get("/ui/state").json()["measurements"]


def _expected(points):
    if len(points) == 2:
        return "Distance", np.linalg.norm(points[1] - points[0]), 3
    if len(points) == 3:
        first, second = points[0] - points[1], points[2] - points[1]
        angle = np.arctan2(
            np.linalg.norm(np.cross(first, second)), np.dot(first, second)
        )
        return "Angle", np.degrees(angle), 1
    # Project the end bonds onto the plane perpendicular to the middle bond.
    axis = points[2] - points[1]
    axis /= np.linalg.norm(axis)
    first, last = points[0] - points[1], points[3] - points[2]
    first -= np.dot(first, axis) * axis
    last -= np.dot(last, axis) * axis
    angle = np.arctan2(np.dot(np.cross(axis, first), last), np.dot(first, last))
    # The extraction and Reader use negative IUPAC torsions.
    return "Dihedral", -np.degrees(angle), 1


@pytest.mark.parametrize("atoms", [[0, 1], [0, 1, 2], [0, 1, 2, 3]])
def test_measurement_matches_displayed_positions_across_frames_and_fits(rest, atoms):
    frames = rest.sampled_frames(3)
    original_transform = rest.client.get("/transform").json()
    try:
        _post(rest, "/selection/atoms", {"atoms": atoms})
        for transform in ("all_atom_fit", "backbone_fit"):
            _post(rest, "/transform", {"kind": transform})
            for frame in frames:
                _post(rest, "/frame/set", {"frame": frame})
                positions = _post(rest, "/positions", {"atoms": atoms, "frame": frame})
                points = np.array(
                    [row["position"] for row in positions.json()["positions"]]
                )
                kind, expected, precision = _expected(points)
                shown = _readout(rest)
                assert shown["frame"] == frame
                assert shown["atoms"] == atoms
                assert shown["kind"] == kind
                assert float(
                    shown["value"].split()[0].rstrip("\u00b0")
                ) == pytest.approx(expected, abs=0.51 * 10**-precision)
                assert len(shown["labels"].splitlines()) == len(atoms)
                assert shown["foreground"], shown
    finally:
        _post(rest, "/transform", original_transform)


def test_repeated_measurement_returns_to_foreground_and_clears(rest):
    _post(rest, "/selection/atoms", {"atoms": [0, 1]})
    assert _readout(rest)["foreground"]
    _post(rest, "/selection/pick", {"atom": 2, "modifiers": "none"})
    assert _readout(rest)["value"] == ""
    _post(rest, "/selection/pick", {"atom": 3, "modifiers": "shift"})
    shown = _readout(rest)
    assert shown["atoms"] == [2, 3]
    assert shown["foreground"], shown
    _post(rest, "/docks/visible", {"visible": False})
    assert not _readout(rest)["visible"]
    _post(rest, "/docks/visible", {"visible": True})
    _post(rest, "/selection/clear")
    shown = _readout(rest)
    assert shown["atoms"] == []
    assert shown["value"] == shown["labels"] == ""


def test_reload_resets_measurement_and_inspector_frame(rest):
    fixture = os.environ["H5READER_REST_FIXTURE"]
    _post(rest, "/frame/set", {"frame": 10})
    _post(rest, "/selection/atoms", {"atoms": [0, 1]})
    _post(rest, "/api/run/load", {"path": fixture})
    shown = _readout(rest)
    assert shown["frame"] == 0
    assert shown["atoms"] == []
    assert shown["value"] == ""
    _post(rest, "/selection/pick", {"atom": 0, "modifiers": "none"})
    state = rest.client.get("/ui/state").json()
    tree = rest.client.get("/inspector/tree").json()
    assert state["currentFrame"] == 0
    assert tree[0]["value"].startswith("frame 1 /"), tree[0]


def test_paused_step_loads_the_displayed_tensor_without_a_probe(rest):
    _post(rest, "/api/run/load", {"path": os.environ["H5READER_REST_FIXTURE"]})
    _post(rest, "/selection/pick", {"atom": 0, "modifiers": "none"})

    for frame in (0, 2):
        _post(rest, "/frame/set", {"frame": frame})
        deadline = time.monotonic() + 30.0
        while True:
            tree = rest.client.get("/inspector/tree").json()
            tensors = [
                node
                for node in tree[0].get("children", [])
                if node["field"] == "Shielding tensor (ORCA DFT)"
            ]
            if tensors:
                assert tree[0]["value"].startswith(f"frame {frame + 1} /")
                isotropic = next(
                    node
                    for node in tensors[0]["children"]
                    if node["field"] == "sigma_iso"
                )
                assert np.isfinite(float(isotropic["value"].split()[0]))
                break
            assert time.monotonic() < deadline, tree
            time.sleep(0.02)

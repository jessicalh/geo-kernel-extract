"""The visible tensor's numerical key uses the same values as its glyph."""

from io import BytesIO
import json
import os
from pathlib import Path
import time

import numpy as np
from PIL import Image
import pytest


def child(node, field):
    return next(item for item in node["children"] if item["field"] == field)


def shown_tensor(rest, field):
    deadline = time.monotonic() + 30
    while True:
        tree = rest.client.get("/inspector/tree").json()
        groups = [item for item in tree[1:] if item["field"] == field]
        if groups and groups[0]["value"] == "Shown":
            return tree[0], groups[0]
        assert time.monotonic() < deadline, tree
        time.sleep(0.02)


def tensor_visibility(rest, shielding, orientation, *, suppressed=False):
    expected = {
        "Shielding tensor (ORCA DFT)": shielding,
        "Bond orientation tensor": orientation,
    }
    deadline = time.monotonic() + 10
    while True:
        groups = {item["field"]: item for item in rest.client.get("/inspector/tree").json()[1:]}
        if all(name in groups and groups[name].get("checked") == visible
               and (groups[name]["value"] == "Shown") == (visible and not suppressed)
               for name, visible in expected.items()):
            return groups
        assert time.monotonic() < deadline, groups
        time.sleep(0.02)


def test_tensor_visibility_is_independent_and_preserves_selection(rest):
    h5py = pytest.importorskip("h5py")
    manifest_path = Path(os.environ["H5READER_REST_FIXTURE"])
    manifest = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
    with h5py.File(manifest_path.parent / manifest["trajectory"]["trajectory_h5"], "r") as trajectory:
        bonds = trajectory["trajectory/reorientational_dynamics"]
        atom = int(bonds["tail_atom"][0])
        bond_atoms = set(bonds["tail_atom"][:]) | set(bonds["head_atom"][:])
        shielding_only_atom = next(index for index in range(max(bond_atoms)) if index not in bond_atoms)
    response = rest.client.post("/selection/pick", json={"atom": atom})
    assert response.is_success, response.text
    original_selection = rest.client.get("/selection").json()
    tensor_visibility(rest, True, True)
    try:
        for shielding, orientation in [(False, True), (True, False), (False, False), (True, True)]:
            for name, visible in [("shielding", shielding), ("orientation", orientation)]:
                response = rest.client.post("/overlay", json={"name": name, "visible": visible})
                assert response.status_code == 204, response.text
            tensor_visibility(rest, shielding, orientation)
            assert rest.client.get("/selection").json() == original_selection
            response = rest.client.post("/frame/set", json={"frame": 2})
            assert response.status_code == 204, response.text
            tensor_visibility(rest, shielding, orientation)
            response = rest.client.post("/selection/clear")
            assert response.is_success, response.text
            response = rest.client.post("/selection/pick", json={"atom": atom})
            assert response.is_success, response.text
            tensor_visibility(rest, shielding, orientation)
        response = rest.client.post("/selection/pick", json={"atom": shielding_only_atom})
        assert response.is_success, response.text
        shown_tensor(rest, "Shielding tensor (ORCA DFT)")
        response = rest.client.post("/overlay", json={"name": "shielding", "visible": False})
        assert response.status_code == 204, response.text
        deadline = time.monotonic() + 10
        while True:
            groups = rest.client.get("/inspector/tree").json()[1:]
            if len(groups) == 1 and groups[0].get("checked") is False and groups[0]["value"] == "":
                break
            assert time.monotonic() < deadline, groups
            time.sleep(0.02)
        state = rest.client.get("/ui/state").json()
        assert state["measurements"]["atoms"] == [shielding_only_atom]
        assert state["measurements"]["visible"]
    finally:
        for name in ["shielding", "orientation"]:
            response = rest.client.post("/overlay", json={"name": name, "visible": True})
            assert response.status_code == 204, response.text


def test_comparison_hero_respects_current_tensor_choices(rest):
    h5py = pytest.importorskip("h5py")
    manifest_path = Path(os.environ["H5READER_REST_FIXTURE"])
    manifest = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
    with h5py.File(manifest_path.parent / manifest["trajectory"]["trajectory_h5"], "r") as trajectory:
        atom = int(trajectory["trajectory/reorientational_dynamics/tail_atom"][0])
    response = rest.client.post("/selection/pick", json={"atom": atom})
    assert response.is_success, response.text
    tensor_visibility(rest, True, True)
    response = rest.client.post("/overlay", json={"name": "shielding", "visible": False})
    assert response.status_code == 204, response.text
    try:
        response = rest.client.post("/resthero/ring_tensor_compare", json={
            "atoms": [atom], "rings": [4], "theta_resolution": 24, "phi_resolution": 24,
        }, timeout=300)
        assert response.status_code == 200, response.text
        tensor_visibility(rest, False, True, suppressed=True)
        for name, visible in [("shielding", True), ("orientation", False)]:
            response = rest.client.post("/overlay", json={"name": name, "visible": visible})
            assert response.status_code == 204, response.text
        response = rest.client.post("/frame/set", json={"frame": 2})
        assert response.status_code == 204, response.text
        tensor_visibility(rest, True, False, suppressed=True)
        response = rest.client.post("/selection/pick", json={"atom": 0})
        assert response.is_success, response.text
        response = rest.client.post("/selection/pick", json={"atom": atom})
        assert response.is_success, response.text
        tensor_visibility(rest, True, False, suppressed=True)
        response = rest.client.post("/resthero/clear")
        assert response.status_code == 204, response.text
        tensor_visibility(rest, True, False)
    finally:
        rest.client.post("/resthero/clear").raise_for_status()
        for name in ["shielding", "orientation"]:
            rest.client.post("/overlay", json={"name": name, "visible": True}).raise_for_status()


def test_shielding_key_matches_glyph_values(rest, tmp_path):
    response = rest.client.post("/selection/pick", json={"atom": 0})
    assert response.is_success, response.text
    for frame in (0, 2):
        response = rest.client.post("/frame/set", json={"frame": frame})
        assert response.is_success, response.text
        tree, tensor = shown_tensor(rest, "Shielding tensor (ORCA DFT)")
        assert all(item["field"] != tensor["field"] for item in tree["children"])
        expected = rest.client.get("/csa", params={"atom": 0, "frame": frame}).json()
        assert expected["valid"]
        mean = expected["sigma_iso"]
        assert float(child(tensor, "sigma_iso")["value"].split()[0]) == pytest.approx(mean, rel=6e-4)
        axes = child(tensor, "Principal values")
        assert axes["value"] == "ppm (from mean)"
        for axis, principal in zip(axes["children"], expected["principal_values"], strict=True):
            value, difference = axis["value"].split(" (")
            assert float(value) == pytest.approx(principal, rel=6e-4, abs=1e-6)
            assert float(difference.removesuffix(")")) == pytest.approx(principal - mean, rel=6e-4, abs=1e-6)
        assert child(tensor, "Glyph size")["value"] == "Normalised"
        assert not any("surface" in node["field"].lower() for node in tensor["children"])
        assert child(child(tensor, "Details"), "Scope")["value"] == "Current frame"

    response = rest.client.post("/diagnostics/screenshot", json={
        "target": "widget", "object_name": "QtAtomInspectorDock"
    })
    assert response.is_success, response.text
    image = Image.open(BytesIO(response.content))
    assert image.width >= 260 and image.height >= 300
    (tmp_path / "tensor-inspector.png").write_bytes(response.content)


def test_bond_orientation_key_matches_stored_tensor(rest, tmp_path):
    h5py = pytest.importorskip("h5py", reason="independent check of the stored bond tensor")
    manifest_path = Path(os.environ["H5READER_REST_FIXTURE"])
    manifest = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
    h5_path = manifest_path.parent / manifest["trajectory"]["trajectory_h5"]
    with h5py.File(h5_path, "r") as trajectory:
        group = trajectory["trajectory/reorientational_dynamics"]
        tail, head = int(group["tail_atom"][0]), int(group["head_atom"][0])
        matrix = group["bond_orientation_tensor"][0]
    eigenvalues = np.linalg.eigvalsh(0.5 * (matrix + matrix.T))[::-1]
    expected_order = (3 * np.sum(eigenvalues**2) - 1) / 2

    response = rest.client.post("/selection/pick", json={"atom": head})
    assert response.is_success, response.text
    head_tree = rest.client.get("/inspector/tree").json()[0]
    head_name = child(child(head_tree, "Identity"), "IUPAC name")["value"]
    response = rest.client.post("/selection/pick", json={"atom": tail})
    assert response.is_success, response.text
    for frame in (0, 2):
        response = rest.client.post("/frame/set", json={"frame": frame})
        assert response.is_success, response.text
        tree, tensor = shown_tensor(rest, "Bond orientation tensor")
        identity = child(tree, "Identity")
        residue, number = child(identity, "Residue")["value"].split(" #")
        chain = child(identity, "Chain")["value"]
        tail_name = child(identity, "IUPAC name")["value"]
        assert child(tensor, "Bond")["value"] == f"{chain}:{residue}{number} {tail_name}-{head_name}"
        for index, expected in enumerate(eigenvalues, start=1):
            actual = float(child(tensor, f"lambda_{index}")["value"])
            assert actual == pytest.approx(expected, rel=6e-4)
        assert float(child(tensor, "S^2 (order parameter)")["value"]) == pytest.approx(expected_order, rel=6e-4, abs=1e-6)
        assert child(tensor, "Average over")["value"] == "Trajectory"
        assert child(tensor, "Main axis follows")["value"] == "Current bond"

    response = rest.client.post("/diagnostics/screenshot", json={
        "target": "widget", "object_name": "QtAtomInspectorDock"
    })
    assert response.is_success, response.text
    (tmp_path / "bond-and-shielding.png").write_bytes(response.content)

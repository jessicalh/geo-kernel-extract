"""Selection context stays inside Atom Info and describes existing state."""

from io import BytesIO
import os

from PIL import Image
import pytest


def post(rest, path, body=None):
    response = rest.client.post(path, json=body)
    assert response.is_success, response.text
    return response


def context(rest):
    return rest.client.get("/ui/state").json()["context"]


@pytest.fixture(autouse=True)
def free_camera(rest):
    post(rest, "/camera/clear")
    yield
    post(rest, "/camera/clear")


def test_hiding_atom_info_keeps_selection_and_camera_across_frames(rest):
    post(rest, "/selection/atoms", {"atoms": [0, 1]})
    post(rest, "/camera/mode", {"mode": "bond", "atoms": [0, 1]})
    shown = context(rest)
    assert shown["visible"] and not shown["window"]
    assert "atom pair" in shown["camera"]
    selected = rest.client.get("/selection").json()
    camera = rest.client.get("/camera/mode").json()

    post(rest, "/docks/visible", {"visible": False})
    try:
        for frame in [1, 2, 3]:
            post(rest, "/frame/set", {"frame": frame})
            assert not context(rest)["visible"]
            assert context(rest)["measurement"]["frame"] == frame
            assert rest.client.get("/selection").json() == selected
            assert rest.client.get("/camera/mode").json() == camera

        post(rest, "/selection/pick", {"atom": 2})
        assert context(rest)["visible"]
        assert not context(rest)["window"]
        assert context(rest)["measurement"]["atoms"] == [2]
        # The lock still refers to the original pair, not the newly picked atom.
        assert rest.client.get("/camera/mode").json() == camera
    finally:
        post(rest, "/docks/visible", {"visible": True})


def test_clear_selection_and_release_camera_are_separate_buttons(rest):
    post(rest, "/selection/atoms", {"atoms": [0, 1]})
    post(rest, "/camera/mode", {"mode": "bond", "atoms": [0, 1]})
    post(rest, "/ui/context", {"action": "clear_selection"})
    cleared = context(rest)
    assert cleared["measurement"]["atoms"] == []
    assert cleared["visible"]  # The camera lock is still worth describing.
    assert rest.client.get("/camera/mode").json()["mode"] == "bond"
    post(rest, "/ui/context", {"action": "release_camera"})
    assert rest.client.get("/camera/mode").json()["mode"] == "free"
    assert context(rest)["visible"]
    assert context(rest)["measurement"]["kind"] == "No atoms selected"


@pytest.mark.parametrize("atoms", [[0], [0, 1], [0, 1, 2], [0, 1, 2, 3]])
def test_selection_context_renders(rest, tmp_path, atoms):
    post(rest, "/selection/atoms", {"atoms": atoms})
    assert rest.client.get("/inspector/tree").json()
    assert context(rest)["measurement"]["atoms"] == atoms
    shot = post(rest, "/diagnostics/screenshot", {
        "target": "widget", "object_name": "QtAtomInspectorDock"
    })
    image = Image.open(BytesIO(shot.content))
    assert 260 <= image.width <= 1200 and image.height >= 300
    assert len(image.getcolors(image.width * image.height)) > 20
    (tmp_path / "selection-context.png").write_bytes(shot.content)


def test_single_pick_opens_the_existing_atom_info_panel(rest):
    post(rest, "/docks/visible", {"visible": False})
    try:
        assert not rest.client.get("/ui/state").json()["atomInfoVisible"]
        post(rest, "/selection/pick", {"atom": 0})
        assert rest.client.get("/ui/state").json()["atomInfoVisible"]
        assert not context(rest)["window"]
        assert rest.client.get("/inspector/tree").json()
    finally:
        post(rest, "/docks/visible", {"visible": True})


def test_geometry_description_names_the_vertex_and_axis(rest):
    post(rest, "/selection/atoms", {"atoms": [0, 1, 2]})
    assert "atom 2" in context(rest)["meaning"]
    assert context(rest)["measurement"]["kind"] == "Angle"
    post(rest, "/selection/pick", {"atom": 3, "modifiers": "shift"})
    assert "atoms 2 and 3" in context(rest)["meaning"]
    assert context(rest)["measurement"]["kind"] == "Dihedral"
    post(rest, "/ui/context", {"action": "clear_selection"})
    assert context(rest)["measurement"]["value"] == ""
    assert len(rest.client.get("/inspector/tree").json()) == 1


def test_context_is_not_a_separate_dismissible_window(rest):
    post(rest, "/selection/clear")
    post(rest, "/docks/visible", {"visible": True})
    assert context(rest)["visible"]
    assert not context(rest)["window"]
    assert context(rest)["measurement"]["kind"] == "No atoms selected"
    assert context(rest)["camera"] == "Camera: free"
    response = rest.client.post("/ui/context", json={"action": "dismiss"})
    assert response.status_code == 400
    assert context(rest)["visible"]


def test_loading_another_molecule_discards_the_old_context(rest):
    second = os.environ.get("H5READER_REST_RELOAD_FIXTURE")
    if not second:
        pytest.skip("a second molecule is not configured")
    original = os.environ["H5READER_REST_FIXTURE"]
    post(rest, "/selection/atoms", {"atoms": [0, 1, 2, 3]})
    try:
        post(rest, "/api/run/load", {"path": second})
        state = context(rest)
        assert not state["window"]
        assert state["measurement"]["atoms"] == []
        assert state["highlight"] == ""
        assert len(rest.client.get("/inspector/tree").json()) == 1
        post(rest, "/selection/pick", {"atom": 0})
        assert context(rest)["visible"]
        assert context(rest)["measurement"]["atoms"] == [0]
    finally:
        post(rest, "/api/run/load", {"path": original})

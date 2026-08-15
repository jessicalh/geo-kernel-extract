"""REST namespace contract for general, resthero, and diagnostics routes."""

from __future__ import annotations


def test_rest_interface_declares_namespaces(rest):
    response = rest.client.get("/api/interface")
    assert response.status_code == 200
    payload = response.json()

    assert payload["version"] == 1
    namespaces = {item["prefix"]: item for item in payload["namespaces"]}
    assert namespaces["/api"]["tier"] == "general"
    assert namespaces["/api"]["audience"] == "user/mcp"
    assert namespaces["/field"]["tier"] == "general"
    assert namespaces["/resthero"]["tier"] == "figure_composition"
    assert "/resthero/clear" in namespaces["/resthero"]["contract"]
    assert namespaces["/catalog"]["tier"] == "diagnostic"
    assert namespaces["/diagnostics"]["tier"] == "diagnostic"
    assert "/" not in namespaces
    assert "compatibility" not in payload

    routes = {(item["method"], item["path"]): item for item in payload["routes"]}
    assert ("GET", "/api/interface") in routes
    assert ("POST", "/api/screenshot") in routes
    assert ("POST", "/api/alignment/export") in routes
    assert ("POST", "/api/ring/null_crossings") in routes
    assert ("POST", "/api/ring/current_face_collar") in routes
    assert ("POST", "/field/null_cone") in routes
    assert ("POST", "/resthero/atom_track") in routes
    assert ("POST", "/diagnostics/screenshot") in routes
    assert ("GET", "/catalog") in routes


def test_ring_api_canonical_routes(rest):
    response = rest.client.post(
        "/api/ring/null_crossings",
        json={
            "atom": 568,
            "ring": 4,
            "start_frame": 538,
            "end_frame": 540,
            "include_signal_stamps": False,
        },
    )
    assert response.status_code == 200, response.text
    payload = response.json()
    assert payload["kind"] == "ring_null_collar_collection"
    assert payload["summary"]["entry_count"] == 1
    assert payload["entries"][0]["from"]["frame"] == 538
    assert payload["entries"][0]["to"]["frame"] == 540

    response = rest.client.post("/api/ring/current_face_collar", json={})
    assert response.status_code == 400
    assert "full_scan" in response.json()["error"]


def test_resthero_and_diagnostics_canonical_routes(rest):
    response = rest.client.post("/resthero/clear")
    assert response.status_code == 204

    response = rest.client.post("/diagnostics/screenshot", json={"target": "bogus"})
    assert response.status_code == 400
    assert "target must" in response.json()["error"]


def test_old_route_spellings_are_not_kept_as_aliases(rest):
    for path in (
        "/rest/interface",
        "/screenshot",
        "/ring/null_crossings",
        "/ring/current_face_collar",
        "/ring_null_crossings",
        "/ring_current_face_collar",
        "/heroshot/clear",
    ):
        response = rest.client.post(path) if path != "/rest/interface" else rest.client.get(path)
        assert response.status_code == 404

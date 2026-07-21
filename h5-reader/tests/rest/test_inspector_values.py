"""Generic Atom Info value-quality checks for any loaded protein fixture."""

from __future__ import annotations


def _walk(nodes):
    for node in nodes:
        yield node
        yield from _walk(node.get("children", []))


def test_inspector_omits_nonfinite_optional_scalars(rest):
    response = rest.client.post(
        "/selection/pick", json={"atom": 16, "modifiers": "none"}
    )
    assert response.status_code == 204, response.text

    tree = rest.client.get("/inspector/tree").json()
    assert tree, "inspector tree empty"
    nonfinite = []
    for node in _walk(tree):
        value = str(node.get("value", "")).strip().lower()
        if (
            value.startswith("nan")
            or value.startswith("inf")
            or value.startswith("-inf")
        ):
            nonfinite.append((node.get("field", ""), node.get("value", "")))
    assert not nonfinite, f"non-finite optional values surfaced: {nonfinite}"

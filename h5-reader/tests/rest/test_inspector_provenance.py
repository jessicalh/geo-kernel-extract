"""Inspector provenance — the curated metric-inventory clarity, asserted via REST.

The Atom Info panel carries each primary metric's provenance + USE/PLACEHOLDER
status as a hover tooltip (METRIC_INVENTORY_20260624). Tooltips do not appear in a
screenshot, so GET /inspector/tree serializes the panel (field / value / tooltip)
and we assert the curated primaries and their provenance labels are present.

The documented fixture (H5READER_REST_FIXTURE = the complete stride-1 1P9J extract)
carries the McConnell producers + ORCA, so the full forward sum is present on a
backbone atom.
"""

from __future__ import annotations


def _flatten(nodes, out):
    for n in nodes:
        out[n["field"]] = n
        if "children" in n:
            _flatten(n["children"], out)
    return out


def _tree_for_atom(rest, atom):
    assert rest.client.post("/selection/pick", json={"atom": atom}).status_code == 204
    tree = rest.client.get("/inspector/tree").json()
    assert isinstance(tree, list) and tree, "inspector tree empty"
    return _flatten(tree, {})


def test_inspector_tree_is_well_formed(rest):
    rows = _tree_for_atom(rest, 16)
    assert any("Atom 16" in f for f in rows), "atom title row missing"
    assert "Identity" in rows
    assert "Position" in rows


def test_forward_sum_primaries_carry_provenance(rest):
    rows = _tree_for_atom(rest, 16)
    assert "Classical forward sum (derived)" in rows, \
        "forward-sum section missing (does the fixture carry the producers?)"
    ring = rows.get("ring contribution")
    assert ring is not None, "ring contribution row missing"
    assert "Biot-Savart" in ring.get("tooltip", ""), ring
    assert "[USE]" in ring.get("tooltip", ""), ring
    mc = rows.get("McConnell contribution")
    assert mc is not None, "McConnell contribution missing (fixture lacks mc_*_bo?)"
    assert "DeltaChi" in mc.get("tooltip", ""), mc
    s0 = rows.get("sigma0 (placeholder baseline)")
    assert s0 is not None, "sigma0 row missing"
    assert "[PLACEHOLDER]" in s0.get("tooltip", ""), s0


def test_efield_efg_promoted_to_primary(rest):
    rows = _tree_for_atom(rest, 16)
    assert "Electric field & EFG" in rows, "promoted E-field/EFG section missing"
    epar = rows.get("signed E_parallel")
    assert epar is not None, "signed E_parallel primary missing"
    assert "Buckingham" in epar.get("tooltip", ""), epar

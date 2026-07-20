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


def _tree_nodes_for_atom(rest, atom):
    assert rest.client.post("/selection/pick", json={"atom": atom}).status_code == 204
    tree = rest.client.get("/inspector/tree").json()
    assert isinstance(tree, list) and tree, "inspector tree empty"
    return tree


def _tree_for_atom(rest, atom):
    tree = _tree_nodes_for_atom(rest, atom)
    return _flatten(tree, {})


def _find_node(nodes, field):
    for n in nodes:
        if n.get("field") == field:
            return n
        child = _find_node(n.get("children", []), field)
        if child is not None:
            return child
    return None


def _child(node, field):
    for child in node.get("children", []):
        if child.get("field") == field:
            return child
    return None


def test_inspector_tree_is_well_formed(rest):
    rows = _tree_for_atom(rest, 16)
    assert any("Atom 16" in f for f in rows), "atom title row missing"
    assert "Identity" in rows
    assert "Position" in rows


def test_forward_sum_primaries_carry_provenance(rest):
    rows = _tree_for_atom(rest, 16)
    assert "Local classical estimate (tentative)" in rows, \
        "local estimate section missing (does the fixture carry the producers?)"
    ring = rows.get("estimated ring contribution")
    assert ring is not None, "ring contribution row missing"
    assert "Biot-Savart" in ring.get("tooltip", ""), ring
    assert "[LOCAL ESTIMATE]" in ring.get("tooltip", ""), ring
    mc = rows.get("estimated McConnell contribution")
    assert mc is not None, "McConnell contribution missing (fixture lacks mc_*_bo?)"
    assert "DeltaChi" in mc.get("tooltip", ""), mc
    assert "[LOCAL ESTIMATE]" in mc.get("tooltip", ""), mc
    assert "sigma0 (placeholder baseline)" not in rows, \
        "placeholder sigma0 must not ship in the inspector"
    assert "estimated sigma_cl (tentative)" not in rows, \
        "absolute sigma_cl should stay hidden while sigma0 is only a placeholder"


def test_efield_efg_promoted_to_primary(rest):
    # The July contract defines E_parallel only for a parent-to-H bond axis.
    # Atom 19 is the backbone amide H of VAL 2; atom 16 is a carbon and must
    # correctly carry NaN for this field.
    rows = _tree_for_atom(rest, 19)
    assert "Electric field & EFG" in rows, "promoted E-field/EFG section missing"
    epar = rows.get("signed E_parallel")
    assert epar is not None, "signed E_parallel primary missing"
    assert "Buckingham" in epar.get("tooltip", ""), epar


def test_tensor_rows_are_nested_summaries(rest):
    tree = _tree_nodes_for_atom(rest, 16)

    shielding = _find_node(tree, "bs_shielding")
    assert shielding is not None, "shielding tensor row missing"
    assert "T0=" in shielding.get("value", ""), shielding
    assert "span=" in shielding.get("value", ""), shielding
    assert "|T2|=" in shielding.get("value", ""), shielding
    assert _child(shielding, "T0 signed iso") is not None, shielding
    assert _child(shielding, "|T2| anisotropy") is not None, shielding
    pas = _child(shielding, "PAS / shape")
    assert pas is not None, shielding
    assert _child(pas, "sigma_11") is not None, pas
    assert _child(pas, "sigma_22") is not None, pas
    assert _child(pas, "sigma_33") is not None, pas
    raw = _child(shielding, "raw irreps")
    assert raw is not None, shielding
    assert _child(raw, "T2 components (library basis)") is not None, raw

    efg = _find_node(tree, "aimnet2_efg")
    assert efg is not None, "EFG tensor row missing"
    assert "Vzz=" in efg.get("value", ""), efg
    assert "eta=" in efg.get("value", ""), efg
    assert "|T2|=" in efg.get("value", ""), efg
    assert _child(efg, "|T2| invariant") is not None, efg
    efg_pas = _child(efg, "PAS / EFG convention")
    assert efg_pas is not None, efg
    assert _child(efg_pas, "Vzz") is not None, efg_pas
    assert _child(efg, "raw T2 components (library basis)") is not None, efg


def test_inspector_does_not_surface_placeholders(rest):
    rows = _tree_for_atom(rest, 16)
    leaks = []
    for field, row in rows.items():
        text = f"{field} {row.get('value', '')} {row.get('tooltip', '')}".lower()
        if "placeholder" in text:
            leaks.append(field)
    assert not leaks, "placeholder text leaked into inspector: " + ", ".join(leaks)

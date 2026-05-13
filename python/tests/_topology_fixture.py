"""Shared fixture helper: write the minimal topology sidecar that the
SDK's strict ``load()`` now requires.

Mirrors the schema declared in ``src/TopologySidecar.cpp``. Used by
``test_larsen_hbond_group.py``, ``test_tripeptide_group.py``, and any
future test that constructs a synthetic extraction directory.

The minimal sidecar is valid but empty: zero bonds, zero rings,
zero ring memberships. That matches the codex first-pass validation
gates — an empty axis is valid as long as it's declared.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


# Byte-for-byte mirrors of the C++ structured dtypes in TopologySidecar.cpp.
_RESIDUES_DTYPE = np.dtype([
    ("residue_index", "<i4"),
    ("chain_id", "|S2"),
    ("residue_number", "<i4"),
    ("insertion_code", "|S1"),
    ("residue_type", "i1"),
    ("amber_residue_3letter", "|S4"),
    ("iupac_residue_3letter", "|S4"),
    ("one_letter", "|S1"),
    ("protonation_variant_index", "i1"),
    ("terminal_state", "i1"),
    ("prev_residue_index", "<i4"),
    ("next_residue_index", "<i4"),
    ("prev_residue_type", "i1"),
    ("next_residue_type", "i1"),
    ("atom_count", "<i4"),
    ("is_proline", "i1"),
    ("is_aromatic", "i1"),
    ("is_titratable", "i1"),
    ("has_amide_h", "i1"),
    ("is_xpro_context", "i1"),
])

_BONDS_DTYPE = np.dtype([
    ("bond_index", "<i4"),
    ("atom_index_a", "<i4"),
    ("atom_index_b", "<i4"),
    ("bond_order", "i1"),
    ("bond_category", "i1"),
    ("is_rotatable", "i1"),
    ("is_aromatic", "i1"),
    ("is_peptide", "i1"),
    ("is_backbone", "i1"),
])

_RINGS_DTYPE = np.dtype([
    ("ring_id", "<i4"),
    ("ring_kind", "i1"),
    ("ring_type_index", "i1"),
    ("atom_count", "i1"),
    ("_pad0", "i1"),
    ("native_axis_index", "<i4"),
    ("parent_residue_index", "<i4"),
    ("parent_residue_number", "<i4"),
    ("fused_partner_ring_id", "<i4"),
])

_RING_MEMBERSHIP_DTYPE = np.dtype([
    ("ring_id", "<i4"),
    ("atom_index", "<i4"),
    ("ring_atom_order", "i1"),
    ("is_vertex", "i1"),
    ("is_substituent", "i1"),
    ("_pad0", "i1"),
])


def write_minimal_topology_sidecar(
    out_dir: Path,
    n_atoms: int,
    n_residues: int = None,
    n_bonds: int = 0,
    n_aromatic_rings: int = 0,
    n_saturated_rings: int = 0,
    protein_id: str = "test_fixture",
) -> None:
    """Write the five files ``TopologySidecar::WriteFeatures`` emits.

    Defaults to an empty topology: zero bonds, zero rings, zero
    ring memberships. The axis-size declarations in the manifest
    still match the protein, so the codex validation invariants
    (atom rows == atom axis size, etc.) hold.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Default residue count: one residue per 4 atoms (matches the
    # fake-extraction convention in existing test_*_group.py fixtures).
    if n_residues is None:
        n_residues = max(1, (n_atoms + 3) // 4)

    # Distribute n_atoms across n_residues -- ceil division so the last
    # residue absorbs the remainder. The manifest's residue axis size and
    # residues.npy atom_count must sum to n_atoms.
    residues = np.zeros(n_residues, dtype=_RESIDUES_DTYPE)
    if n_residues > 0:
        base = n_atoms // n_residues
        rem = n_atoms - base * n_residues
        for i in range(n_residues):
            residues[i]["residue_index"] = i
            residues[i]["residue_number"] = i + 1
            residues[i]["residue_type"] = 0  # AminoAcid enum 0 = ALA in fixtures
            residues[i]["atom_count"] = base + (1 if i < rem else 0)
            residues[i]["prev_residue_index"] = i - 1 if i > 0 else -1
            residues[i]["next_residue_index"] = i + 1 if i + 1 < n_residues else -1
    np.save(out_dir / "residues.npy", residues)

    bonds = np.zeros(n_bonds, dtype=_BONDS_DTYPE)
    for i in range(n_bonds):
        bonds[i]["bond_index"] = i
        bonds[i]["atom_index_a"] = min(i, max(0, n_atoms - 1))
        bonds[i]["atom_index_b"] = min(i + 1, max(0, n_atoms - 1))
    np.save(out_dir / "bonds.npy", bonds)

    rings = np.zeros(n_aromatic_rings + n_saturated_rings, dtype=_RINGS_DTYPE)
    # Aromatic rows come first, saturated rows next. ring_id is the
    # absolute row; native_axis_index is contiguous within each kind.
    for i in range(n_aromatic_rings):
        rings[i]["ring_id"] = i
        rings[i]["ring_kind"] = 0
        rings[i]["native_axis_index"] = i
    for si in range(n_saturated_rings):
        i = n_aromatic_rings + si
        rings[i]["ring_id"] = i
        rings[i]["ring_kind"] = 1
        rings[i]["native_axis_index"] = si
    np.save(out_dir / "rings.npy", rings)

    np.save(out_dir / "ring_membership.npy",
            np.zeros(0, dtype=_RING_MEMBERSHIP_DTYPE))

    manifest = {
        "schema_version": "1.0",
        "extractor": "nmr_extract",
        "generated_at_utc": "1970-01-01T00:00:00Z",
        "protein_id": protein_id,
        "topology": {
            "source": "test_fixture",
            "has_atom_semantic": False,
            "has_ff_atom_types": False,
            "has_ff_mass": False,
        },
        "axis_sizes": {
            "atom": n_atoms,
            "residue": n_residues,
            "bond": n_bonds,
            "aromatic_ring": n_aromatic_rings,
            "saturated_ring": n_saturated_rings,
            "ring": n_aromatic_rings + n_saturated_rings,
            "ring_membership": 0,
        },
        "axis_alignment": {
            "atom": "fixture: row[i] == atom_index i for all atom-axis arrays",
        },
    }
    (out_dir / "extraction_manifest.json").write_text(json.dumps(manifest))

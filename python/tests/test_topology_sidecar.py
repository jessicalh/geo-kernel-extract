"""SDK round-trip tests for the topology sidecar.

Verifies that ``nmr_extract.load()`` exposes the new ``Protein.topology``
group (TopologyGroup) with ``bonds``, ``rings``, ``ring_membership``,
and parsed ``manifest`` attributes, and that the catalog declares
``native_axis`` / ``mechanism`` / ``irreps`` / ``units`` /
``sign_convention`` / ``tensor_rank`` / ``parity`` for every entry.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from nmr_extract import (
    CATALOG,
    Residues,
    Bonds,
    Rings,
    RingMembership,
    ExtractionManifest,
    TopologyGroup,
    load,
)
from _topology_fixture import (
    write_minimal_topology_sidecar,
    _RESIDUES_DTYPE,
    _BONDS_DTYPE,
    _RINGS_DTYPE,
    _RING_MEMBERSHIP_DTYPE,
)


N_ATOMS = 16


# ── Catalog metadata invariants ────────────────────────────────────


class TestCatalogMetadata:
    """Resolves OI-016: every catalog entry declares axis + irreps + units."""

    def test_every_entry_declares_native_axis(self):
        bad = [k for k, s in CATALOG.items() if not s.native_axis]
        assert not bad, f"Entries missing native_axis: {bad}"

    def test_every_entry_declares_mechanism(self):
        bad = [k for k, s in CATALOG.items() if not s.mechanism]
        assert not bad, f"Entries missing mechanism: {bad}"

    def test_every_entry_has_valid_parity(self):
        bad = [k for k, s in CATALOG.items() if s.parity not in ("even", "odd")]
        assert not bad, f"Entries with invalid parity: {bad}"

    def test_tensor_rank_is_in_0_1_2(self):
        bad = [k for k, s in CATALOG.items() if s.tensor_rank not in (0, 1, 2)]
        assert not bad, f"Entries with invalid tensor_rank: {bad}"

    def test_topology_sidecar_entries_registered(self):
        for stem in ("residues", "bonds", "rings", "ring_membership"):
            assert stem in CATALOG, f"CATALOG missing topology entry {stem!r}"
            assert CATALOG[stem].required, (
                f"topology sidecar entry {stem!r} must be required")
            assert CATALOG[stem].mechanism == "topology"

    def test_axis_value_set(self):
        """Native axes are drawn from a closed set declared by the contract."""
        allowed = {
            "atom", "residue", "aromatic_ring", "saturated_ring", "ring",
            "ring_contribution_pair", "bond", "ring_membership",
            "mutation_match_pair", "protein",
        }
        seen = {s.native_axis for s in CATALOG.values()}
        unknown = seen - allowed
        assert not unknown, f"Unknown native_axis values: {unknown}"

    def test_shielding_entries_carry_sign_convention(self):
        """Every shielding tensor entry declares its sign convention."""
        for stem, spec in CATALOG.items():
            if "shielding" in stem and spec.tensor_rank == 2:
                assert spec.sign_convention, (
                    f"shielding entry {stem!r} missing sign_convention")


# ── TopologyGroup load + wrappers ──────────────────────────────────


def _required_identity_npys(out_dir, n_atoms):
    """Minimal required-by-catalog NPYs the loader checks for."""
    np.save(out_dir / "pos.npy", np.zeros((n_atoms, 3), dtype=np.float64))
    np.save(out_dir / "element.npy", np.full(n_atoms, 6, dtype=np.int32))
    np.save(out_dir / "residue_index.npy",
            np.arange(n_atoms, dtype=np.int32) // 4)
    np.save(out_dir / "residue_type.npy",
            np.zeros(n_atoms, dtype=np.int32))
    np.save(out_dir / "ring_contributions.npy",
            np.zeros((0, 59), dtype=np.float64))
    np.save(out_dir / "ring_geometry.npy",
            np.zeros((0, 10), dtype=np.float64))


def _required_calculator_npys(out_dir, n_atoms):
    stems_9 = ["bs_shielding", "hm_shielding", "pq_shielding",
               "disp_shielding", "ringchi_shielding",
               "mc_shielding", "hbond_shielding"]
    for s in stems_9:
        np.save(out_dir / f"{s}.npy", np.zeros((n_atoms, 9), dtype=np.float64))
    for stem in ("bs", "hm", "pq", "disp"):
        np.save(out_dir / f"{stem}_per_type_T0.npy",
                np.zeros((n_atoms, 8), dtype=np.float64))
        np.save(out_dir / f"{stem}_per_type_T2.npy",
                np.zeros((n_atoms, 40), dtype=np.float64))
    np.save(out_dir / "bs_total_B.npy", np.zeros((n_atoms, 3), dtype=np.float64))
    np.save(out_dir / "bs_ring_counts.npy",
            np.zeros((n_atoms, 4), dtype=np.float64))
    np.save(out_dir / "mc_category_T2.npy",
            np.zeros((n_atoms, 25), dtype=np.float64))
    np.save(out_dir / "mc_scalars.npy",
            np.zeros((n_atoms, 6), dtype=np.float64))
    np.save(out_dir / "hbond_scalars.npy",
            np.zeros((n_atoms, 3), dtype=np.float64))
    np.save(out_dir / "dssp_backbone.npy",
            np.zeros((n_atoms, 5), dtype=np.float64))
    np.save(out_dir / "aimnet2_charges.npy",
            np.zeros(n_atoms, dtype=np.float64))
    np.save(out_dir / "aimnet2_aim.npy",
            np.zeros((n_atoms, 256), dtype=np.float32))
    # EFG schema rev 2026-05-18: T2 only (5 components).
    np.save(out_dir / "aimnet2_efg.npy",
            np.zeros((n_atoms, 5), dtype=np.float64))
    np.save(out_dir / "aimnet2_efg_aromatic.npy",
            np.zeros((n_atoms, 5), dtype=np.float64))
    np.save(out_dir / "aimnet2_efg_backbone.npy",
            np.zeros((n_atoms, 5), dtype=np.float64))


class TestTopologyLoad:

    @pytest.fixture
    def fake_extraction(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        # 4 residues, 4 atoms each. Two synthetic bonds, one aromatic
        # ring with 6 members, one saturated ring with 5 members.
        n_residues = N_ATOMS // 4
        residues = np.zeros(n_residues, dtype=_RESIDUES_DTYPE)
        for i in range(n_residues):
            residues[i]["residue_index"] = i
            residues[i]["residue_number"] = i + 1
            residues[i]["atom_count"] = 4
            residues[i]["prev_residue_index"] = i - 1 if i > 0 else -1
            residues[i]["next_residue_index"] = i + 1 if i + 1 < n_residues else -1
        np.save(tmp_path / "residues.npy", residues)

        bonds = np.zeros(2, dtype=_BONDS_DTYPE)
        bonds[0] = (0, 0, 1, 0, 0, 0, 0, 0, 1)   # bond_index, a, b, order=Single, cat=PeptideCO, rotate, arom, peptide, backbone
        bonds[1] = (1, 1, 2, 4, 1, 0, 0, 1, 1)   # peptide bond
        np.save(tmp_path / "bonds.npy", bonds)

        rings = np.zeros(2, dtype=_RINGS_DTYPE)
        rings[0] = (0, 0, 0, 6, 0, 0, 1, 100, -1)  # aromatic PheBenzene
        rings[1] = (1, 1, 8, 5, 0, 0, 2, 200, -1)  # saturated Pro
        np.save(tmp_path / "rings.npy", rings)

        memb = np.zeros(11, dtype=_RING_MEMBERSHIP_DTYPE)
        for k in range(6):
            memb[k] = (0, k, k, 1, 0, 0)
        for k in range(5):
            memb[6 + k] = (1, 6 + k, k, 1, 0, 0)
        np.save(tmp_path / "ring_membership.npy", memb)

        manifest = {
            "schema_version": "1.0",
            "extractor": "nmr_extract",
            "generated_at_utc": "2026-05-13T00:00:00Z",
            "protein_id": "test_fixture",
            "topology": {"has_atom_semantic": True,
                         "has_ff_atom_types": False,
                         "has_ff_mass": False},
            "axis_sizes": {
                "atom": N_ATOMS, "residue": n_residues, "bond": 2,
                "aromatic_ring": 1, "saturated_ring": 1, "ring": 2,
                "ring_membership": 11,
            },
            "axis_alignment": {"atom": "rows align with atom_index"},
        }
        (tmp_path / "extraction_manifest.json").write_text(json.dumps(manifest))
        return tmp_path

    def test_residues_wrapper(self, fake_extraction):
        p = load(fake_extraction)
        r = p.topology.residues
        assert isinstance(r, Residues)
        assert r.n_residues == 4
        # Total atom_count must sum to atom axis size (validation invariant).
        assert int(r.atom_count.sum()) == N_ATOMS
        # First residue has prev = -1 (chain start).
        assert r.prev_residue_index[0] == -1
        # Last residue has next = -1 (chain end).
        assert r.next_residue_index[-1] == -1

    def test_protein_topology_is_a_topology_group(self, fake_extraction):
        p = load(fake_extraction)
        assert isinstance(p.topology, TopologyGroup)

    def test_bonds_wrapper(self, fake_extraction):
        p = load(fake_extraction)
        b = p.topology.bonds
        assert isinstance(b, Bonds)
        assert b.n_bonds == 2
        assert b.atom_index_a[0] == 0 and b.atom_index_b[0] == 1
        assert b.atom_index_a[1] == 1 and b.atom_index_b[1] == 2
        assert b.is_peptide[0] == False  # bond 0 is Single
        assert b.is_peptide[1] == True   # bond 1 is Peptide
        assert b.is_backbone[0] == True
        assert b.is_backbone[1] == True

    def test_rings_wrapper(self, fake_extraction):
        p = load(fake_extraction)
        r = p.topology.rings
        assert isinstance(r, Rings)
        assert r.n_rings == 2
        assert r.is_aromatic[0]
        assert r.is_saturated[1]
        assert r.atom_count[0] == 6
        assert r.atom_count[1] == 5
        assert r.parent_residue_number[0] == 100
        assert r.parent_residue_number[1] == 200
        assert r.fused_partner_ring_id[0] == -1

    def test_ring_membership_wrapper(self, fake_extraction):
        p = load(fake_extraction)
        m = p.topology.ring_membership
        assert isinstance(m, RingMembership)
        assert m.n_rows == 11
        # First six rows reference ring_id 0 (aromatic).
        assert np.all(m.ring_id[:6] == 0)
        # Next five reference ring 1 (saturated).
        assert np.all(m.ring_id[6:] == 1)
        # All rows are vertices in our model.
        assert np.all(m.is_vertex)

    def test_manifest_wrapper(self, fake_extraction):
        p = load(fake_extraction)
        man = p.topology.manifest
        assert isinstance(man, ExtractionManifest)
        assert man.schema_version == "1.0"
        assert man.protein_id == "test_fixture"
        assert man.has_atom_semantic()
        assert man.axis_size("atom") == N_ATOMS
        assert man.axis_size("bond") == 2
        assert man.axis_size("ring") == 2
        assert man.axis_size("ring_membership") == 11


class TestStrictLoad:
    """Loader rejects extractions without the topology sidecar."""

    def test_missing_bonds_npy_fails(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        # No residues / bonds / rings / ring_membership / manifest.
        with pytest.raises(FileNotFoundError):
            load(tmp_path)

    def test_missing_manifest_fails(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        # residues/bonds/rings/ring_membership present; manifest absent.
        write_minimal_topology_sidecar(tmp_path, n_atoms=N_ATOMS)
        (tmp_path / "extraction_manifest.json").unlink()
        with pytest.raises(FileNotFoundError, match="extraction_manifest"):
            load(tmp_path)


class TestValidationInvariants:
    """SDK enforces the codex first-pass validation gates at load time."""

    def test_axis_size_mismatch_raises(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        # Write a topology sidecar with a manifest that lies about axis sizes.
        write_minimal_topology_sidecar(tmp_path, n_atoms=N_ATOMS)
        # Re-emit a corrupted manifest claiming a different atom count.
        manifest_path = tmp_path / "extraction_manifest.json"
        man = json.loads(manifest_path.read_text())
        man["axis_sizes"]["atom"] = N_ATOMS + 7   # lie
        manifest_path.write_text(json.dumps(man))
        with pytest.raises(ValueError, match="atom axis"):
            load(tmp_path)

    def test_bond_endpoint_out_of_range_raises(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        write_minimal_topology_sidecar(tmp_path, n_atoms=N_ATOMS, n_bonds=2)
        # Corrupt bonds.npy with an out-of-range endpoint.
        bonds = np.load(tmp_path / "bonds.npy")
        bonds[0]["atom_index_a"] = N_ATOMS + 100  # invalid
        np.save(tmp_path / "bonds.npy", bonds)
        with pytest.raises(ValueError, match="bonds.npy"):
            load(tmp_path)

    def test_ring_membership_out_of_range_raises(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        write_minimal_topology_sidecar(
            tmp_path, n_atoms=N_ATOMS, n_aromatic_rings=1)
        # Membership rows are empty by default; need to add a corrupt one.
        memb = np.zeros(1, dtype=_RING_MEMBERSHIP_DTYPE)
        memb[0]["ring_id"] = 99  # invalid (only 1 ring exists)
        memb[0]["atom_index"] = 0
        memb[0]["is_vertex"] = 1
        np.save(tmp_path / "ring_membership.npy", memb)
        man_path = tmp_path / "extraction_manifest.json"
        man = json.loads(man_path.read_text())
        man["axis_sizes"]["ring_membership"] = 1
        man_path.write_text(json.dumps(man))
        with pytest.raises(ValueError, match="ring_id"):
            load(tmp_path)

    def test_residue_atom_count_sum_mismatch_raises(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        write_minimal_topology_sidecar(tmp_path, n_atoms=N_ATOMS)
        # Corrupt residues.npy so atom_count sum != N_ATOMS.
        residues = np.load(tmp_path / "residues.npy")
        residues[0]["atom_count"] = 999
        np.save(tmp_path / "residues.npy", residues)
        with pytest.raises(ValueError, match="atom_count sums"):
            load(tmp_path)

    def test_residue_row_identity_mismatch_raises(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        write_minimal_topology_sidecar(tmp_path, n_atoms=N_ATOMS)
        residues = np.load(tmp_path / "residues.npy")
        residues[1]["residue_index"] = 99
        np.save(tmp_path / "residues.npy", residues)
        with pytest.raises(ValueError, match="residue_index"):
            load(tmp_path)

    def test_bond_row_identity_mismatch_raises(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        write_minimal_topology_sidecar(tmp_path, n_atoms=N_ATOMS, n_bonds=2)
        bonds = np.load(tmp_path / "bonds.npy")
        bonds[1]["bond_index"] = 99
        np.save(tmp_path / "bonds.npy", bonds)
        with pytest.raises(ValueError, match="bond_index"):
            load(tmp_path)

    def test_ring_row_identity_mismatch_raises(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        write_minimal_topology_sidecar(
            tmp_path, n_atoms=N_ATOMS, n_aromatic_rings=2)
        rings = np.load(tmp_path / "rings.npy")
        rings[1]["ring_id"] = 99
        np.save(tmp_path / "rings.npy", rings)
        with pytest.raises(ValueError, match="ring_id"):
            load(tmp_path)

    def test_atom_residue_counts_mismatch_raises(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        write_minimal_topology_sidecar(tmp_path, n_atoms=N_ATOMS)
        # Keep all residue references in range but redistribute one atom
        # onto residue 0; residues.npy atom_count no longer matches.
        residue_index = np.load(tmp_path / "residue_index.npy")
        residue_index[-1] = 0
        np.save(tmp_path / "residue_index.npy", residue_index)
        with pytest.raises(ValueError, match="residue_index.npy counts"):
            load(tmp_path)

    def test_atom_residue_out_of_range_raises(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        write_minimal_topology_sidecar(tmp_path, n_atoms=N_ATOMS)
        residue_index = np.load(tmp_path / "residue_index.npy")
        residue_index[0] = 99
        np.save(tmp_path / "residue_index.npy", residue_index)
        with pytest.raises(ValueError, match="residue_index.npy references"):
            load(tmp_path)

    def test_ring_parent_residue_out_of_range_raises(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        write_minimal_topology_sidecar(
            tmp_path, n_atoms=N_ATOMS, n_aromatic_rings=1)
        rings = np.load(tmp_path / "rings.npy")
        rings[0]["parent_residue_index"] = 99
        np.save(tmp_path / "rings.npy", rings)
        with pytest.raises(ValueError, match="parent_residue_index"):
            load(tmp_path)

    def test_ring_native_axis_noncontiguous_raises(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        write_minimal_topology_sidecar(
            tmp_path, n_atoms=N_ATOMS, n_aromatic_rings=2)
        rings = np.load(tmp_path / "rings.npy")
        rings[1]["native_axis_index"] = 3
        np.save(tmp_path / "rings.npy", rings)
        with pytest.raises(ValueError, match="native_axis_index"):
            load(tmp_path)

    def test_ring_atom_count_membership_mismatch_raises(self, tmp_path):
        _required_identity_npys(tmp_path, N_ATOMS)
        _required_calculator_npys(tmp_path, N_ATOMS)
        write_minimal_topology_sidecar(
            tmp_path, n_atoms=N_ATOMS, n_aromatic_rings=1)
        rings = np.load(tmp_path / "rings.npy")
        rings[0]["atom_count"] = 2
        np.save(tmp_path / "rings.npy", rings)

        memb = np.zeros(1, dtype=_RING_MEMBERSHIP_DTYPE)
        memb[0]["ring_id"] = 0
        memb[0]["atom_index"] = 0
        memb[0]["is_vertex"] = 1
        np.save(tmp_path / "ring_membership.npy", memb)
        man_path = tmp_path / "extraction_manifest.json"
        man = json.loads(man_path.read_text())
        man["axis_sizes"]["ring_membership"] = 1
        man_path.write_text(json.dumps(man))

        with pytest.raises(ValueError, match="atom_count does not match"):
            load(tmp_path)


# The 6 new fields on atoms_category_info are exercised end-to-end on
# the C++ side at tests/test_category_info_projection.cpp::WriteFeaturesEmitsNpy
# (header check) and the round-trip from real extraction is implicit:
# accessing ``info.chain_id`` on an NPY without the column raises
# ValueError from numpy structured-array access (loud fail), which is
# the contract.

"""
Loader: builds an AnalysisProtein from an analysis trajectory H5 file.

The H5 is the only input. Self-contained.
"""

from __future__ import annotations

import h5py
import numpy as np
from pathlib import Path
from typing import Union

from .enums import (
    Element, AtomRole, Hybridisation, BondOrder, BondCategory,
    RingTypeIndex, AminoAcid, amino_acid_from_code,
)
from .model import (
    AnalysisProtein, AnalysisAtom, AnalysisResidue, AnalysisRing, AnalysisBond,
    PositionsGroup, RingCurrentGroup, EfgGroup, BondAnisoGroup,
    QuadrupoleGroup, DispersionGroup, HBondGroup, SasaGroup, WaterGroup,
    ChargesGroup, AimNet2EmbeddingGroup, PerRingGroup, RingGeometryGroup,
)


# Map from atomic number (H5 storage) to Element enum.
_ELEMENT_FROM_ATOMIC_NUMBER = {e.value: e for e in Element}


def _safe_element(atomic_number: int) -> Element:
    """Convert atomic number to Element, raising on unknown values."""
    elem = _ELEMENT_FROM_ATOMIC_NUMBER.get(atomic_number)
    if elem is None:
        raise ValueError(
            f"Unknown atomic number {atomic_number} in atoms/element. "
            f"Expected one of: {list(_ELEMENT_FROM_ATOMIC_NUMBER.keys())}"
        )
    return elem


def _safe_enum(enum_cls, value: int):
    """Convert an integer to an IntEnum member, returning the last member on unknown."""
    try:
        return enum_cls(value)
    except ValueError:
        # Return the "Unknown" or last member as fallback
        members = list(enum_cls)
        return members[-1]


def load_analysis(path: Union[str, Path]) -> AnalysisProtein:
    """
    Load an analysis trajectory H5 file and return an AnalysisProtein.

    The returned object keeps the H5 file open for lazy dataset access.
    Use it as a context manager or call .close() when done:

        with load_analysis("path.h5") as protein:
            ...

    Or:
        protein = load_analysis("path.h5")
        ...
        protein.close()
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"H5 file not found: {path}")

    h5 = h5py.File(str(path), "r")

    protein = AnalysisProtein()
    protein._h5 = h5

    # ================================================================
    # Metadata
    # ================================================================

    meta = h5["meta"]
    protein.protein_id = str(meta.attrs["protein_id"])
    protein.n_atoms = int(meta.attrs["n_atoms"])
    protein.n_frames = int(meta.attrs["n_frames"])
    protein.n_residues = int(meta.attrs["n_residues"])
    protein.stride = int(meta.attrs["stride"])
    protein.frame_times = meta["frame_times"][:]
    protein.frame_indices = meta["frame_indices"][:]

    N = protein.n_atoms
    R = protein.n_residues

    # ================================================================
    # Atoms (constant identity)
    # ================================================================

    atoms_grp = h5["atoms"]
    elements_raw = atoms_grp["element"][:]             # (N,) int32 atomic numbers
    residue_indices = atoms_grp["residue_index"][:]    # (N,) int32
    atom_roles = atoms_grp["atom_role"][:]             # (N,) int32
    hybridisations = atoms_grp["hybridisation"][:]     # (N,) int32
    is_backbone = atoms_grp["is_backbone"][:]          # (N,) int8
    is_conjugated = atoms_grp["is_conjugated"][:]      # (N,) int8
    n_bonded = atoms_grp["n_bonded"][:]                # (N,) int32
    graph_dist_ring = atoms_grp["graph_dist_ring"][:]  # (N,) int32
    partial_charges = atoms_grp["partial_charge"][:]   # (N,) float64
    vdw_radii = atoms_grp["vdw_radius"][:]             # (N,) float64

    # PDB atom names -- string dataset
    atom_names_raw = atoms_grp["atom_name"][:]
    atom_names = [
        x.decode() if isinstance(x, bytes) else str(x)
        for x in atom_names_raw
    ]

    # Optional enrichment boolean fields
    def _read_bool_array(name: str) -> np.ndarray:
        if name in atoms_grp:
            return atoms_grp[name][:].astype(bool)
        return np.zeros(N, dtype=bool)

    is_amide_H = _read_bool_array("is_amide_H")
    is_alpha_H = _read_bool_array("is_alpha_H")
    is_methyl = _read_bool_array("is_methyl")
    is_aromatic_H = _read_bool_array("is_aromatic_H")
    is_on_aromatic_residue = _read_bool_array("is_on_aromatic_residue")
    is_hbond_donor = _read_bool_array("is_hbond_donor")
    is_hbond_acceptor = _read_bool_array("is_hbond_acceptor")
    parent_is_sp2 = _read_bool_array("parent_is_sp2")

    # Optional graph topology fields
    def _read_int_array(name: str, default: int = -1) -> np.ndarray:
        if name in atoms_grp:
            return atoms_grp[name][:]
        return np.full(N, default, dtype=np.int32)

    def _read_float_array(name: str, default: float = 0.0) -> np.ndarray:
        if name in atoms_grp:
            return atoms_grp[name][:]
        return np.full(N, default, dtype=np.float64)

    graph_dist_N = _read_int_array("graph_dist_N")
    graph_dist_O = _read_int_array("graph_dist_O")
    n_pi_bonds_3 = _read_int_array("n_pi_bonds_3", default=0)
    eneg_sum_1 = _read_float_array("eneg_sum_1")
    eneg_sum_2 = _read_float_array("eneg_sum_2")
    bfs_to_nearest_ring = _read_int_array("bfs_to_nearest_ring")
    bfs_decay = _read_float_array("bfs_decay")

    # Build atom objects
    atoms = []
    for i in range(N):
        a = AnalysisAtom(protein, i)
        a.element = _safe_element(int(elements_raw[i]))
        a.pdb_name = atom_names[i]
        a.residue_index = int(residue_indices[i])
        a.role = _safe_enum(AtomRole, int(atom_roles[i]))
        a.hybridisation = _safe_enum(Hybridisation, int(hybridisations[i]))
        a.is_backbone = bool(is_backbone[i])
        a.is_conjugated = bool(is_conjugated[i])
        a.n_bonded = int(n_bonded[i])
        a.graph_dist_ring = int(graph_dist_ring[i])
        a.partial_charge = float(partial_charges[i])
        a.vdw_radius = float(vdw_radii[i])
        a.bond_indices = []  # filled below from topology

        # Enrichment booleans
        a.is_amide_H = bool(is_amide_H[i])
        a.is_alpha_H = bool(is_alpha_H[i])
        a.is_methyl = bool(is_methyl[i])
        a.is_aromatic_H = bool(is_aromatic_H[i])
        a.is_on_aromatic_residue = bool(is_on_aromatic_residue[i])
        a.is_hbond_donor = bool(is_hbond_donor[i])
        a.is_hbond_acceptor = bool(is_hbond_acceptor[i])
        a.parent_is_sp2 = bool(parent_is_sp2[i])

        # Graph topology
        a.graph_dist_N = int(graph_dist_N[i])
        a.graph_dist_O = int(graph_dist_O[i])
        a.n_pi_bonds_3 = int(n_pi_bonds_3[i])
        a.eneg_sum_1 = float(eneg_sum_1[i])
        a.eneg_sum_2 = float(eneg_sum_2[i])
        a.bfs_to_nearest_ring = int(bfs_to_nearest_ring[i])
        a.bfs_decay = float(bfs_decay[i])

        atoms.append(a)

    protein._atoms = atoms

    # ================================================================
    # Residues
    # ================================================================

    res_grp = h5["residues"]
    res_names_raw = res_grp["residue_name"][:]
    res_numbers = res_grp["residue_number"][:]
    res_chains_raw = res_grp["chain_id"][:]

    residues = []
    for i in range(R):
        r = AnalysisResidue(protein, i)
        name_str = res_names_raw[i]
        if isinstance(name_str, bytes):
            name_str = name_str.decode()
        r.amino_acid = amino_acid_from_code(str(name_str))
        r.sequence_number = int(res_numbers[i])
        chain_str = res_chains_raw[i]
        if isinstance(chain_str, bytes):
            chain_str = chain_str.decode()
        r.chain_id = str(chain_str)
        r.atom_indices = []
        residues.append(r)

    protein._residues = residues

    # Wire atoms to residues
    for a in atoms:
        if 0 <= a.residue_index < R:
            residues[a.residue_index].atom_indices.append(a.index)

    # ================================================================
    # Topology: bonds
    # ================================================================

    topo = h5["topology"]
    bond_atoms_arr = topo["bond_atoms"][:]          # (B, 2) int32
    bond_category_arr = topo["bond_category"][:]    # (B,) int32
    bond_order_arr = topo["bond_order"][:]          # (B,) int32
    parent_atom_idx = topo["parent_atom_index"][:]  # (N,) int32

    n_bonds = int(topo.attrs["n_bonds"])
    bonds = []
    for i in range(n_bonds):
        b = AnalysisBond(protein, i)
        b.atom_index_a = int(bond_atoms_arr[i, 0])
        b.atom_index_b = int(bond_atoms_arr[i, 1])
        b.order = _safe_enum(BondOrder, int(bond_order_arr[i]))
        b.category = _safe_enum(BondCategory, int(bond_category_arr[i]))
        bonds.append(b)

        # Wire bond indices to atoms
        if b.atom_index_a < N:
            atoms[b.atom_index_a].bond_indices.append(i)
        if b.atom_index_b < N:
            atoms[b.atom_index_b].bond_indices.append(i)

    protein._bonds = bonds

    # Parent atom index (hydrogen -> heavy atom mapping)
    for i in range(N):
        atoms[i].parent_atom_index = int(parent_atom_idx[i])

    # ================================================================
    # Topology: rings (CSR ragged array)
    # ================================================================

    n_rings = int(topo.attrs["n_rings"])
    ring_type_arr = topo["ring_type"][:]            # (n_rings,) int32
    ring_residue_arr = topo["ring_residue"][:]      # (n_rings,) int32
    ring_fused_arr = topo["ring_fused_partner"][:]  # (n_rings,) int32
    ring_offsets = topo["ring_offsets"][:]           # (n_rings+1,) int32
    ring_atom_indices = topo["ring_atom_indices"][:]  # (K,) int32

    rings = []
    for i in range(n_rings):
        r = AnalysisRing(protein, i)
        r.type_index = _safe_enum(RingTypeIndex, int(ring_type_arr[i]))
        r.residue_index = int(ring_residue_arr[i])
        r.fused_partner_index = int(ring_fused_arr[i])
        start = int(ring_offsets[i])
        end = int(ring_offsets[i + 1])
        r.atom_indices = [int(x) for x in ring_atom_indices[start:end]]
        rings.append(r)

    protein._rings = rings

    # ================================================================
    # Physics groups (lazy h5py dataset references)
    # ================================================================

    if "positions" in h5:
        protein.positions = PositionsGroup(h5["positions/xyz"])

    if "ring_current" in h5:
        protein.ring_current = RingCurrentGroup(h5["ring_current"])

    if "efg" in h5:
        protein.efg = EfgGroup(h5["efg"])

    if "bond_aniso" in h5:
        protein.bond_aniso = BondAnisoGroup(h5["bond_aniso"])

    if "quadrupole" in h5:
        protein.quadrupole = QuadrupoleGroup(h5["quadrupole"])

    if "dispersion" in h5:
        protein.dispersion = DispersionGroup(h5["dispersion"])

    if "hbond" in h5:
        protein.hbond = HBondGroup(h5["hbond"])

    if "sasa" in h5:
        protein.sasa = SasaGroup(h5["sasa"])

    if "water" in h5:
        protein.water = WaterGroup(h5["water"])

    if "charges" in h5:
        protein.charges = ChargesGroup(h5["charges"])

    if "aimnet2_embedding" in h5:
        protein.aimnet2_embedding = AimNet2EmbeddingGroup(h5["aimnet2_embedding"])

    if "per_ring" in h5:
        protein.per_ring = PerRingGroup(h5["per_ring"])

    if "ring_geometry" in h5:
        protein.ring_geometry = RingGeometryGroup(h5["ring_geometry"])

    return protein

"""
Python object model for analysis trajectory H5 files.

The design mirrors the C++ Protein/ConformationAtom split:
  - Protein = identity + topology (constant)
  - Conformation = frame index into (T, N, ...) arrays

AnalysisProtein is built once from the H5. Frame data is accessed
as numpy slices through the atom or through physics group accessors.
There are no materialised frame objects.
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Sequence

from .enums import (
    Element, AtomRole, Hybridisation, BondOrder, BondCategory,
    RingTypeIndex, AminoAcid, amino_acid_from_code,
)


# ============================================================================
# SphericalTensor: thin view over a 9-element array
# Layout: [T0, T1_x, T1_y, T1_z, T2_m-2, T2_m-1, T2_m0, T2_m+1, T2_m+2]
# ============================================================================

class SphericalTensorView:
    """
    View over a (..., 9) array giving named access to T0, T1, T2 components.

    Can wrap a single (9,) array or a time series (T, 9) array.
    Never copies data -- all access returns views/slices of the backing array.
    """
    __slots__ = ("_data",)

    def __init__(self, data: np.ndarray):
        if data.shape[-1] != 9:
            raise ValueError(f"SphericalTensor expects last dim=9, got {data.shape}")
        self._data = data

    @property
    def raw(self) -> np.ndarray:
        """The full (..., 9) backing array."""
        return self._data

    @property
    def T0(self) -> np.ndarray:
        """Isotropic component. Shape (...)."""
        return self._data[..., 0]

    @property
    def T1(self) -> np.ndarray:
        """Antisymmetric (pseudovector) components. Shape (..., 3)."""
        return self._data[..., 1:4]

    @property
    def T2(self) -> np.ndarray:
        """Traceless symmetric components (m=-2..+2). Shape (..., 5)."""
        return self._data[..., 4:9]

    @property
    def T2_magnitude(self) -> np.ndarray:
        """L2 norm of T2 components. Shape (...)."""
        return np.sqrt(np.sum(self._data[..., 4:9] ** 2, axis=-1))


# ============================================================================
# Atom: identity only. Position and computed data are in the time series.
# ============================================================================

class AnalysisAtom:
    """
    One atom in the protein. Identity and topology are constant.
    Time-varying data is accessed through the parent protein's arrays
    using this atom's index.

    Usage:
        atom = protein.atom(42)
        atom.element          # Element.C
        atom.pdb_name         # "CA"
        atom.is_backbone      # True

        # Time series -- all frames:
        protein.ring_current.bs_T0_per_type[:, 42, :]   # (T, 8)

        # Or through the atom convenience:
        atom.position(t=0)    # (3,) xyz at frame 0
        atom.position()       # (T, 3) all frames
    """
    __slots__ = (
        "_protein", "index",
        "element", "pdb_name", "residue_index", "role", "hybridisation",
        "is_backbone", "is_conjugated", "n_bonded", "graph_dist_ring",
        "partial_charge", "vdw_radius",
        "parent_atom_index", "bond_indices",
        # Graph topology fields
        "graph_dist_N", "graph_dist_O", "n_pi_bonds_3",
        "eneg_sum_1", "eneg_sum_2",
        "bfs_to_nearest_ring", "bfs_decay",
        # Boolean enrichment flags
        "is_amide_H", "is_alpha_H", "is_methyl",
        "is_aromatic_H", "is_on_aromatic_residue",
        "is_hbond_donor", "is_hbond_acceptor", "parent_is_sp2",
    )

    def __init__(self, protein: AnalysisProtein, index: int):
        self._protein = protein
        self.index = index

    @property
    def residue(self) -> AnalysisResidue:
        """The residue this atom belongs to."""
        return self._protein.residue(self.residue_index)

    @property
    def parent_atom(self) -> Optional[AnalysisAtom]:
        """For hydrogen atoms: the bonded heavy atom. None otherwise."""
        if self.parent_atom_index < 0:
            return None
        return self._protein.atom(self.parent_atom_index)

    @property
    def is_hydrogen(self) -> bool:
        return self.element == Element.H

    @property
    def bonds(self) -> list[AnalysisBond]:
        """All bonds involving this atom."""
        return [self._protein.bond(i) for i in self.bond_indices]

    # --- Time series convenience ---

    def position(self, t: Optional[int] = None) -> np.ndarray:
        """
        Atom position(s) in Angstroms.
        t=None: all frames, shape (T, 3)
        t=int:  one frame, shape (3,)
        """
        xyz = self._protein.positions.xyz
        if t is not None:
            return xyz[t, self.index, :]
        return xyz[:, self.index, :]

    def __repr__(self) -> str:
        return (
            f"AnalysisAtom({self.index}, {self.element.symbol} "
            f"{self.pdb_name!r}, res={self.residue_index})"
        )


# ============================================================================
# Residue
# ============================================================================

class AnalysisResidue:
    """
    One amino acid residue. Identity is constant.
    """
    __slots__ = (
        "_protein", "index",
        "amino_acid", "sequence_number", "chain_id",
        "atom_indices",
    )

    def __init__(self, protein: AnalysisProtein, index: int):
        self._protein = protein
        self.index = index
        self.atom_indices: list[int] = []

    @property
    def atoms(self) -> list[AnalysisAtom]:
        """All atoms belonging to this residue."""
        return [self._protein.atom(i) for i in self.atom_indices]

    @property
    def name(self) -> str:
        """Three-letter code (e.g. 'ALA')."""
        return self.amino_acid.three_letter_code

    def __repr__(self) -> str:
        return (
            f"AnalysisResidue({self.index}, {self.amino_acid.name} "
            f"{self.sequence_number}, chain={self.chain_id!r})"
        )


# ============================================================================
# Ring
# ============================================================================

class AnalysisRing:
    """
    One aromatic ring. Identity is constant. Geometry varies per frame
    and is in ring_geometry/.
    """
    __slots__ = (
        "_protein", "index",
        "type_index", "residue_index", "fused_partner_index",
        "atom_indices",
    )

    def __init__(self, protein: AnalysisProtein, index: int):
        self._protein = protein
        self.index = index
        self.atom_indices: list[int] = []
        self.fused_partner_index: int = -1

    @property
    def residue(self) -> AnalysisResidue:
        return self._protein.residue(self.residue_index)

    @property
    def atoms(self) -> list[AnalysisAtom]:
        return [self._protein.atom(i) for i in self.atom_indices]

    @property
    def is_fused(self) -> bool:
        return self.fused_partner_index >= 0

    @property
    def fused_partner(self) -> Optional[AnalysisRing]:
        if self.fused_partner_index < 0:
            return None
        return self._protein.ring(self.fused_partner_index)

    @property
    def short_name(self) -> str:
        return self.type_index.short_name

    @property
    def ring_size(self) -> int:
        return self.type_index.ring_size

    def __repr__(self) -> str:
        return (
            f"AnalysisRing({self.index}, {self.type_index.short_name}, "
            f"res={self.residue_index}, atoms={self.atom_indices})"
        )


# ============================================================================
# Bond
# ============================================================================

class AnalysisBond:
    """
    One covalent bond. Topology is constant. Length/direction vary per frame.
    """
    __slots__ = (
        "_protein", "index",
        "atom_index_a", "atom_index_b",
        "order", "category",
    )

    def __init__(self, protein: AnalysisProtein, index: int):
        self._protein = protein
        self.index = index

    @property
    def atom_a(self) -> AnalysisAtom:
        return self._protein.atom(self.atom_index_a)

    @property
    def atom_b(self) -> AnalysisAtom:
        return self._protein.atom(self.atom_index_b)

    def length(self, t: int) -> float:
        """Bond length at frame t in Angstroms."""
        pa = self._protein.positions.xyz[t, self.atom_index_a, :]
        pb = self._protein.positions.xyz[t, self.atom_index_b, :]
        return float(np.linalg.norm(pb - pa))

    def __repr__(self) -> str:
        return (
            f"AnalysisBond({self.index}, {self.atom_index_a}-{self.atom_index_b}, "
            f"{self.category.name})"
        )


# ============================================================================
# Physics group accessors: thin wrappers that hold h5py dataset references.
# No data is loaded until you slice.
# ============================================================================

class PositionsGroup:
    """positions/ group. xyz is (T, N, 3) float64, Angstroms, PBC-fixed."""
    __slots__ = ("xyz",)

    def __init__(self, xyz: "h5py.Dataset"):
        self.xyz = xyz


class RingCurrentGroup:
    """ring_current/ group. All Biot-Savart, Haigh-Mallion, and ring susceptibility data."""
    __slots__ = (
        "bs_T0_per_type", "bs_T2_per_type", "bs_shielding",
        "hm_T0_per_type", "hm_T2_per_type", "hm_shielding",
        "rs_shielding",
        "n_rings_3A", "n_rings_5A", "n_rings_8A", "n_rings_12A",
        "mean_ring_dist", "nearest_ring_atom",
        "G_iso_exp_sum", "G_T2_exp_sum", "G_iso_var_8A",
        "total_B_field",
    )

    def __init__(self, grp):
        self.bs_T0_per_type = grp["bs_T0_per_type"]     # (T,N,8)
        self.bs_T2_per_type = grp["bs_T2_per_type"]     # (T,N,8,5)
        self.bs_shielding = grp["bs_shielding"]         # (T,N,9)
        self.hm_T0_per_type = grp["hm_T0_per_type"]     # (T,N,8)
        self.hm_T2_per_type = grp["hm_T2_per_type"]     # (T,N,8,5)
        self.hm_shielding = grp["hm_shielding"]         # (T,N,9)
        self.rs_shielding = grp["rs_shielding"]         # (T,N,9)
        self.n_rings_3A = grp["n_rings_3A"]             # (T,N)
        self.n_rings_5A = grp["n_rings_5A"]             # (T,N)
        self.n_rings_8A = grp["n_rings_8A"]             # (T,N)
        self.n_rings_12A = grp["n_rings_12A"]           # (T,N)
        self.mean_ring_dist = grp["mean_ring_dist"]     # (T,N)
        self.nearest_ring_atom = grp["nearest_ring_atom"]  # (T,N)
        self.G_iso_exp_sum = grp["G_iso_exp_sum"]       # (T,N)
        self.G_T2_exp_sum = grp["G_T2_exp_sum"]         # (T,N,5)
        self.G_iso_var_8A = grp["G_iso_var_8A"]         # (T,N)
        self.total_B_field = grp["total_B_field"]       # (T,N,3)


class EfgGroup:
    """efg/ group. Coulomb EFG, E-fields, APBS, AIMNet2 EFG."""
    __slots__ = (
        "coulomb_total", "coulomb_backbone", "coulomb_aromatic", "coulomb_shielding",
        "E_total", "E_backbone", "E_sidechain", "E_aromatic", "E_solvent",
        "E_magnitude", "E_bond_proj", "E_backbone_frac",
        "apbs_efg", "apbs_efield",
        "aimnet2_total", "aimnet2_backbone", "aimnet2_aromatic", "aimnet2_shielding",
    )

    def __init__(self, grp):
        self.coulomb_total = grp["coulomb_total"]         # (T,N,9)
        self.coulomb_backbone = grp["coulomb_backbone"]   # (T,N,9)
        self.coulomb_aromatic = grp["coulomb_aromatic"]   # (T,N,9)
        self.coulomb_shielding = grp["coulomb_shielding"] # (T,N,9)
        self.E_total = grp["E_total"]                     # (T,N,3)
        self.E_backbone = grp["E_backbone"]               # (T,N,3)
        self.E_sidechain = grp["E_sidechain"]             # (T,N,3)
        self.E_aromatic = grp["E_aromatic"]               # (T,N,3)
        self.E_solvent = grp["E_solvent"]                 # (T,N,3)
        self.E_magnitude = grp["E_magnitude"]             # (T,N)
        self.E_bond_proj = grp["E_bond_proj"]             # (T,N)
        self.E_backbone_frac = grp["E_backbone_frac"]     # (T,N)
        self.apbs_efg = grp["apbs_efg"]                   # (T,N,9)
        self.apbs_efield = grp["apbs_efield"]             # (T,N,3)
        self.aimnet2_total = grp["aimnet2_total"]         # (T,N,9)
        self.aimnet2_backbone = grp["aimnet2_backbone"]   # (T,N,9)
        self.aimnet2_aromatic = grp["aimnet2_aromatic"]   # (T,N,9)
        self.aimnet2_shielding = grp["aimnet2_shielding"] # (T,N,9)


class BondAnisoGroup:
    """bond_aniso/ group. McConnell dipolar anisotropy."""
    __slots__ = (
        "mc_shielding",
        "T2_backbone", "T2_sidechain", "T2_aromatic",
        "T2_CO_nearest", "T2_CN_nearest",
        "co_sum", "cn_sum", "sidechain_sum", "aromatic_sum",
        "co_nearest", "nearest_CO_dist", "nearest_CN_dist",
        "nearest_CO_midpoint", "dir_nearest_CO",
    )

    def __init__(self, grp):
        self.mc_shielding = grp["mc_shielding"]               # (T,N,9)
        self.T2_backbone = grp["T2_backbone"]                 # (T,N,9)
        self.T2_sidechain = grp["T2_sidechain"]               # (T,N,9)
        self.T2_aromatic = grp["T2_aromatic"]                 # (T,N,9)
        self.T2_CO_nearest = grp["T2_CO_nearest"]             # (T,N,9)
        self.T2_CN_nearest = grp["T2_CN_nearest"]             # (T,N,9)
        self.co_sum = grp["co_sum"]                           # (T,N)
        self.cn_sum = grp["cn_sum"]                           # (T,N)
        self.sidechain_sum = grp["sidechain_sum"]             # (T,N)
        self.aromatic_sum = grp["aromatic_sum"]               # (T,N)
        self.co_nearest = grp["co_nearest"]                   # (T,N)
        self.nearest_CO_dist = grp["nearest_CO_dist"]         # (T,N)
        self.nearest_CN_dist = grp["nearest_CN_dist"]         # (T,N)
        self.nearest_CO_midpoint = grp["nearest_CO_midpoint"] # (T,N,3)
        self.dir_nearest_CO = grp["dir_nearest_CO"]           # (T,N,3)


class QuadrupoleGroup:
    """quadrupole/ group. Pi-quadrupole contributions per ring type."""
    __slots__ = ("pq_shielding", "pq_T0_per_type", "pq_T2_per_type")

    def __init__(self, grp):
        self.pq_shielding = grp["pq_shielding"]       # (T,N,9)
        self.pq_T0_per_type = grp["pq_T0_per_type"]   # (T,N,8)
        self.pq_T2_per_type = grp["pq_T2_per_type"]   # (T,N,8,5)


class DispersionGroup:
    """dispersion/ group. London dispersion contributions per ring type."""
    __slots__ = ("disp_shielding", "disp_T0_per_type", "disp_T2_per_type")

    def __init__(self, grp):
        self.disp_shielding = grp["disp_shielding"]       # (T,N,9)
        self.disp_T0_per_type = grp["disp_T0_per_type"]   # (T,N,8)
        self.disp_T2_per_type = grp["disp_T2_per_type"]   # (T,N,8,5)


class HBondGroup:
    """hbond/ group. Hydrogen bond geometry and tensor data."""
    __slots__ = (
        "hbond_shielding", "nearest_spherical",
        "nearest_dist", "nearest_dir", "inv_d3",
        "count_3_5A", "is_donor", "is_acceptor", "is_backbone",
    )

    def __init__(self, grp):
        self.hbond_shielding = grp["hbond_shielding"]     # (T,N,9)
        self.nearest_spherical = grp["nearest_spherical"] # (T,N,9)
        self.nearest_dist = grp["nearest_dist"]           # (T,N)
        self.nearest_dir = grp["nearest_dir"]             # (T,N,3)
        self.inv_d3 = grp["inv_d3"]                       # (T,N)
        self.count_3_5A = grp["count_3_5A"]               # (T,N)
        self.is_donor = grp["is_donor"]                   # (T,N)
        self.is_acceptor = grp["is_acceptor"]             # (T,N)
        self.is_backbone = grp["is_backbone"]             # (T,N)


class SasaGroup:
    """sasa/ group. Solvent-accessible surface area."""
    __slots__ = ("sasa", "normal")

    def __init__(self, grp):
        self.sasa = grp["sasa"]       # (T,N) float64 Angstrom^2
        self.normal = grp["normal"]   # (T,N,3) outward surface normal


class WaterGroup:
    """water/ group. Explicit solvent fields and hydration shell geometry."""
    __slots__ = (
        "efield", "efg", "efield_first", "efg_first",
        "n_first", "n_second",
        "half_shell_asymmetry", "dipole_cos",
        "nearest_ion_dist", "nearest_ion_charge",
        "dipole_vector", "surface_normal",
        "sasa_asymmetry", "sasa_dipole_align", "sasa_dipole_cohere",
        "sasa_first_shell_n",
    )

    def __init__(self, grp):
        self.efield = grp["efield"]                   # (T,N,3)
        self.efg = grp["efg"]                         # (T,N,9)
        self.efield_first = grp["efield_first"]       # (T,N,3)
        self.efg_first = grp["efg_first"]             # (T,N,9)
        self.n_first = grp["n_first"]                 # (T,N) int16
        self.n_second = grp["n_second"]               # (T,N) int16
        self.half_shell_asymmetry = grp["half_shell_asymmetry"]  # (T,N)
        self.dipole_cos = grp["dipole_cos"]           # (T,N)
        self.nearest_ion_dist = grp["nearest_ion_dist"]   # (T,N)
        self.nearest_ion_charge = grp["nearest_ion_charge"] # (T,N)
        self.dipole_vector = grp["dipole_vector"]         # (T,N,3)
        self.surface_normal = grp["surface_normal"]       # (T,N,3)
        self.sasa_asymmetry = grp["sasa_asymmetry"]       # (T,N)
        self.sasa_dipole_align = grp["sasa_dipole_align"] # (T,N)
        self.sasa_dipole_cohere = grp["sasa_dipole_cohere"] # (T,N)
        self.sasa_first_shell_n = grp["sasa_first_shell_n"] # (T,N) int16


class ChargesGroup:
    """charges/ group. Per-frame partial charges."""
    __slots__ = ("aimnet2_charge", "eeq_charge", "eeq_cn")

    def __init__(self, grp):
        self.aimnet2_charge = grp["aimnet2_charge"]  # (T,N) Hirshfeld (e)
        self.eeq_charge = grp["eeq_charge"]          # (T,N) Caldeweyher (e)
        self.eeq_cn = grp["eeq_cn"]                  # (T,N) coordination number


class AimNet2EmbeddingGroup:
    """aimnet2_embedding/ group. Learned electronic structure embedding."""
    __slots__ = ("aim",)

    def __init__(self, grp):
        self.aim = grp["aim"]  # (T,N,256) float64


class PerRingGroup:
    """
    per_ring/ group. K nearest rings per atom, with geometry and T2 tensors.

    All arrays have shape (T, N, K, ...) where K=6 (configurable, read from attrs).

    ring_type: (T, N, K) int8 -- RingTypeIndex of the k-th nearest ring.
               Value -1 means no ring at that slot (atom has fewer than K nearby rings).

    geometry: (T, N, K, F) float64 where F = len(geometry_fields).
              Fields named by geometry_fields (e.g. ["distance", "rho", "z", ...]).

    bs_T2, hm_T2, hm_H_T2, pq_T2, chi_T2, disp_T2: (T, N, K, 5) float64
              T2 components from each calculator for each of the K nearest rings.

    disp_scalar: (T, N, K) float64 -- dispersion scalar per ring.
    """
    __slots__ = (
        "K", "geometry_fields",
        "ring_type", "geometry",
        "bs_T2", "hm_T2", "hm_H_T2", "pq_T2", "chi_T2", "disp_T2",
        "disp_scalar",
    )

    def __init__(self, grp):
        self.K = int(grp.attrs["K"])
        gf_raw = grp["geometry_fields"][:]
        self.geometry_fields = [
            x.decode() if isinstance(x, bytes) else str(x) for x in gf_raw
        ]
        self.ring_type = grp["ring_type"]        # (T,N,K) int8
        self.geometry = grp["geometry"]           # (T,N,K,F)
        self.bs_T2 = grp["bs_T2"]                # (T,N,K,5)
        self.hm_T2 = grp["hm_T2"]                # (T,N,K,5)
        self.hm_H_T2 = grp["hm_H_T2"]           # (T,N,K,5)
        self.pq_T2 = grp["pq_T2"]                # (T,N,K,5)
        self.chi_T2 = grp["chi_T2"]              # (T,N,K,5)
        self.disp_T2 = grp["disp_T2"]            # (T,N,K,5)
        self.disp_scalar = grp["disp_scalar"]     # (T,N,K)

    def geometry_field_index(self, name: str) -> int:
        """Index of a named geometry field (e.g. 'distance', 'rho')."""
        return self.geometry_fields.index(name)


class RingGeometryGroup:
    """
    ring_geometry/ group. Per-ring, per-frame geometry (center, normal, radius, etc.).

    data: (T, n_rings, F) float64 where F = len(fields).
    fields: names of the geometry columns.
    """
    __slots__ = ("data", "fields", "_field_index")

    def __init__(self, grp):
        self.data = grp["data"]  # (T, n_rings, F)
        raw_fields = grp["fields"][:]
        self.fields = [x.decode() if isinstance(x, bytes) else str(x) for x in raw_fields]
        self._field_index = {name: i for i, name in enumerate(self.fields)}

    def field_index(self, name: str) -> int:
        return self._field_index[name]


# ============================================================================
# AnalysisProtein: the top-level object
# ============================================================================

class AnalysisProtein:
    """
    A protein loaded from an analysis trajectory H5 file.

    Topology (atoms, residues, rings, bonds) is constant.
    Time-varying data is accessed through physics group properties,
    each of which holds h5py dataset references for lazy loading.

    The H5 file handle is kept open for the lifetime of this object.
    Call close() or use as a context manager to release it.
    """

    def __init__(self):
        # Set by loader
        self._h5 = None  # h5py.File, kept open

        # Metadata
        self.protein_id: str = ""
        self.n_atoms: int = 0
        self.n_frames: int = 0
        self.n_residues: int = 0
        self.stride: int = 1
        self.frame_times: np.ndarray = np.array([])   # (T,) ps
        self.frame_indices: np.ndarray = np.array([])  # (T,) original XTC frame numbers

        # Topology objects
        self._atoms: list[AnalysisAtom] = []
        self._residues: list[AnalysisResidue] = []
        self._rings: list[AnalysisRing] = []
        self._bonds: list[AnalysisBond] = []

        # Physics groups (set by loader, hold h5py dataset references)
        self.positions: Optional[PositionsGroup] = None
        self.ring_current: Optional[RingCurrentGroup] = None
        self.efg: Optional[EfgGroup] = None
        self.bond_aniso: Optional[BondAnisoGroup] = None
        self.quadrupole: Optional[QuadrupoleGroup] = None
        self.dispersion: Optional[DispersionGroup] = None
        self.hbond: Optional[HBondGroup] = None
        self.sasa: Optional[SasaGroup] = None
        self.water: Optional[WaterGroup] = None
        self.charges: Optional[ChargesGroup] = None
        self.aimnet2_embedding: Optional[AimNet2EmbeddingGroup] = None
        self.per_ring: Optional[PerRingGroup] = None
        self.ring_geometry: Optional[RingGeometryGroup] = None

    # --- Topology access ---

    def atom(self, index: int) -> AnalysisAtom:
        """Get atom by index (0-based)."""
        return self._atoms[index]

    def residue(self, index: int) -> AnalysisResidue:
        """Get residue by index (0-based)."""
        return self._residues[index]

    def ring(self, index: int) -> AnalysisRing:
        """Get ring by index (0-based)."""
        return self._rings[index]

    def bond(self, index: int) -> AnalysisBond:
        """Get bond by index (0-based)."""
        return self._bonds[index]

    @property
    def atoms(self) -> list[AnalysisAtom]:
        return self._atoms

    @property
    def residues(self) -> list[AnalysisResidue]:
        return self._residues

    @property
    def rings(self) -> list[AnalysisRing]:
        return self._rings

    @property
    def bonds(self) -> list[AnalysisBond]:
        return self._bonds

    @property
    def n_rings(self) -> int:
        return len(self._rings)

    @property
    def n_bonds(self) -> int:
        return len(self._bonds)

    # --- Query helpers ---

    def atoms_by_element(self, element: Element) -> list[AnalysisAtom]:
        """All atoms of a given element."""
        return [a for a in self._atoms if a.element == element]

    def atoms_by_residue(self, residue_index: int) -> list[AnalysisAtom]:
        """All atoms in a residue."""
        return self._residues[residue_index].atoms

    def backbone_atoms(self) -> list[AnalysisAtom]:
        """All backbone atoms."""
        return [a for a in self._atoms if a.is_backbone]

    def find_atoms(self, element: Optional[Element] = None,
                   role: Optional[AtomRole] = None,
                   residue_index: Optional[int] = None) -> list[AnalysisAtom]:
        """Filter atoms by element, role, and/or residue."""
        result = self._atoms
        if element is not None:
            result = [a for a in result if a.element == element]
        if role is not None:
            result = [a for a in result if a.role == role]
        if residue_index is not None:
            result = [a for a in result if a.residue_index == residue_index]
        return result

    # --- SphericalTensor convenience ---

    def spherical_tensor(self, dataset_path: str,
                         atom_index: Optional[int] = None,
                         frame: Optional[int] = None) -> SphericalTensorView:
        """
        Wrap a (T, N, 9) dataset as a SphericalTensorView, optionally
        selecting a specific atom and/or frame.

        Examples:
            protein.spherical_tensor("ring_current/bs_shielding")
                -> SphericalTensorView over (T, N, 9)

            protein.spherical_tensor("ring_current/bs_shielding", atom_index=42)
                -> SphericalTensorView over (T, 9)

            protein.spherical_tensor("ring_current/bs_shielding", atom_index=42, frame=0)
                -> SphericalTensorView over (9,)
        """
        ds = self._h5[dataset_path]
        if atom_index is not None and frame is not None:
            data = ds[frame, atom_index, :]
        elif atom_index is not None:
            data = ds[:, atom_index, :]
        elif frame is not None:
            data = ds[frame, :, :]
        else:
            data = ds[:]
        return SphericalTensorView(data)

    # --- Raw H5 access ---

    @property
    def h5(self):
        """The underlying h5py.File handle, for direct access to any dataset."""
        return self._h5

    def dataset(self, path: str):
        """Access any h5py dataset by path (e.g. 'charges/eeq_charge')."""
        return self._h5[path]

    # --- Lifecycle ---

    def close(self):
        """Close the underlying H5 file."""
        if self._h5 is not None:
            self._h5.close()
            self._h5 = None

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def __repr__(self) -> str:
        return (
            f"AnalysisProtein({self.protein_id!r}, "
            f"{self.n_atoms} atoms, {self.n_residues} residues, "
            f"{self.n_rings} rings, {self.n_frames} frames)"
        )

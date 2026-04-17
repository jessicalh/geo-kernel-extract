"""
Python object model for analysis trajectory H5 files.

Design principles:
    1. Protein is identity/topology. Built once, never changes.
    2. Time-varying data is accessed as (T,...) numpy array slices through
       the atom or residue, NOT as materialized per-frame objects.
    3. FrameView is a lightweight index wrapper for "give me everything at frame t".
    4. All enum fields use the typed Python enums, never raw ints.
    5. The H5 datasets are loaded lazily as numpy arrays and held in memory.
       For 1B1V (335 atoms, 626 frames), ring_current alone is ~210 MB.
       That is fine on 64 GB machines. For memory-constrained use, the
       loader can be extended to use h5py dataset objects (lazy I/O) instead.

Naming: classes are prefixed Analysis* to distinguish from the C++ types
that live in the nmr namespace of the engine.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from numpy.typing import NDArray

from .enums import (
    Element,
    AtomRole,
    Hybridisation,
    BondOrder,
    BondCategory,
    RingTypeIndex,
    AminoAcid,
    amino_acid_from_code,
)


# ---------------------------------------------------------------------------
# Per-group time series containers
# ---------------------------------------------------------------------------
# Each holds the (T, N, ...) arrays for one H5 group, attributed to atoms
# by index. The atom object gets a reference to slice into these.
# Groups that may not exist in the H5 (not yet written) are Optional.
# ---------------------------------------------------------------------------

@dataclass
class TimeSeries:
    """
    Raw numpy arrays for one H5 group, exactly as stored.

    The arrays are (T, N, ...) where T = n_frames, N = n_atoms (or R = n_residues
    for per-residue groups). Atom/residue access is by index slicing: arr[:, i, ...].

    Field names match the H5 dataset names within their group.
    """
    _arrays: dict[str, NDArray] = field(default_factory=dict, repr=False)

    def __getattr__(self, name: str) -> NDArray:
        if name.startswith("_"):
            raise AttributeError(name)
        arrays = object.__getattribute__(self, "_arrays")
        if name in arrays:
            return arrays[name]
        raise AttributeError(
            f"Dataset '{name}' not present in this group. "
            f"Available: {sorted(arrays.keys())}"
        )

    def __contains__(self, name: str) -> bool:
        return name in self._arrays

    def datasets(self) -> list[str]:
        """List available dataset names."""
        return sorted(self._arrays.keys())


# Typed group containers — one per H5 group.
# These are thin wrappers that document what datasets to expect.
# The actual data lives in the underlying TimeSeries._arrays dict.

class RingCurrentData(TimeSeries):
    """
    ring_current/ group.

    Expected datasets:
        bs_T0_per_type       (T,N,8)     Biot-Savart isotropic per ring type
        bs_T2_per_type       (T,N,8,5)   Biot-Savart T2 per ring type
        bs_shielding         (T,N,9)     Biot-Savart total SphericalTensor
        hm_T0_per_type       (T,N,8)     Haigh-Mallion isotropic per ring type
        hm_T2_per_type       (T,N,8,5)   Haigh-Mallion T2 per ring type
        hm_shielding         (T,N,9)     Haigh-Mallion total SphericalTensor
        rs_shielding         (T,N,9)     Ring susceptibility
        n_rings_3A           (T,N)       Ring count within 3A
        n_rings_5A           (T,N)       Ring count within 5A
        n_rings_8A           (T,N)       Ring count within 8A
        n_rings_12A          (T,N)       Ring count within 12A
        mean_ring_dist       (T,N)       Mean distance to rings
        nearest_ring_atom    (T,N)       Nearest ring atom distance
        G_iso_exp_sum        (T,N)       Exponential-weighted isotropic sum
        G_T2_exp_sum         (T,N,5)     Exponential-weighted T2 sum
        G_iso_var_8A         (T,N)       Isotropic variance within 8A
        total_B_field        (T,N,3)     Total induced B-field
    """
    pass


class EfgData(TimeSeries):
    """
    efg/ group — electric field gradients.

    Expected datasets: coulomb_total, coulomb_backbone, coulomb_aromatic (T,N,9),
    E_total ... E_solvent (T,N,3), E_magnitude (T,N), E_bond_proj (T,N),
    E_backbone_frac (T,N), apbs_efg (T,N,9), apbs_efield (T,N,3),
    aimnet2_total/backbone/aromatic/shielding (T,N,9).
    """
    pass


class BondAnisoData(TimeSeries):
    """
    bond_aniso/ group — McConnell bond anisotropy.

    Expected datasets: mc_shielding (T,N,9), T2_backbone/sidechain/aromatic/
    CO_nearest/CN_nearest (T,N,9), co_sum/cn_sum/... scalars (T,N),
    nearest_CO_dist (T,N), dir_nearest_CO (T,N,3), etc.
    """
    pass


class QuadrupoleData(TimeSeries):
    """quadrupole/ group — pi-quadrupole shielding."""
    pass


class DispersionData(TimeSeries):
    """dispersion/ group — ring dispersion shielding."""
    pass


class HBondData(TimeSeries):
    """hbond/ group — hydrogen bond fields."""
    pass


class SasaData(TimeSeries):
    """sasa/ group — solvent-accessible surface area."""
    pass


class WaterData(TimeSeries):
    """water/ group — explicit solvent fields and hydration shell geometry."""
    pass


class ChargesData(TimeSeries):
    """charges/ group — AIMNet2 and EEQ charges."""
    pass


class AimNet2EmbeddingData(TimeSeries):
    """aimnet2_embedding/ group — learned electronic embeddings."""
    pass


class PredictionsData(TimeSeries):
    """predictions/ group — ridge and MLP predictions, per-group contributions."""
    pass


class ProjectionsData(TimeSeries):
    """projections/ group — subspace coords, residual projections, Mahalanobis."""
    pass


class DihedralData(TimeSeries):
    """dihedrals/ group — per-RESIDUE (T,R) dihedral angles."""
    pass


class DsspData(TimeSeries):
    """dssp/ group — per-RESIDUE (T,R) secondary structure."""
    pass


# ---------------------------------------------------------------------------
# Atom view — per-atom slice into the time series arrays
# ---------------------------------------------------------------------------

class AtomView:
    """
    Provides time series access for a single atom.

    Usage:
        atom = protein.atoms[42]
        # Full time series (T,) or (T,...):
        t0_all_frames = atom.ring_current.bs_T0_per_type  # shape (T, 8)
        # Single frame:
        t0_frame_5 = atom.ring_current.bs_T0_per_type[5]  # shape (8,)
    """

    def __init__(self, atom_index: int, group: TimeSeries):
        self._idx = atom_index
        self._group = group

    def __getattr__(self, name: str) -> NDArray:
        if name.startswith("_"):
            raise AttributeError(name)
        # Get the (T, N, ...) array and slice to (T, ...)
        full = getattr(self._group, name)
        return full[:, self._idx]

    def datasets(self) -> list[str]:
        return self._group.datasets()


class ResidueView:
    """Per-residue time series access, analogous to AtomView."""

    def __init__(self, res_index: int, group: TimeSeries):
        self._idx = res_index
        self._group = group

    def __getattr__(self, name: str) -> NDArray:
        if name.startswith("_"):
            raise AttributeError(name)
        full = getattr(self._group, name)
        return full[:, self._idx]

    def datasets(self) -> list[str]:
        return self._group.datasets()


# ---------------------------------------------------------------------------
# Topology objects
# ---------------------------------------------------------------------------

class AnalysisAtom:
    """
    One atom in the protein. Identity and topology only — no position.

    Time-varying data is accessed through property-based views that slice
    into the protein-level arrays:
        atom.ring_current.bs_T0_per_type   -> (T, 8) for this atom
        atom.positions                     -> (T, 3) xyz time series
    """

    __slots__ = (
        "index", "element", "pdb_name", "role", "hybridisation",
        "residue_index", "is_backbone", "is_conjugated",
        "n_bonded", "graph_dist_ring", "parent_atom_index",
        "bond_indices", "_protein",
    )

    def __init__(
        self,
        index: int,
        element: Element,
        pdb_name: str,
        role: AtomRole,
        hybridisation: Hybridisation,
        residue_index: int,
        is_backbone: bool,
        is_conjugated: bool,
        n_bonded: int,
        graph_dist_ring: int,
        parent_atom_index: int,  # -1 = no parent (heavy atom)
    ):
        self.index = index
        self.element = element
        self.pdb_name = pdb_name
        self.role = role
        self.hybridisation = hybridisation
        self.residue_index = residue_index
        self.is_backbone = is_backbone
        self.is_conjugated = is_conjugated
        self.n_bonded = n_bonded
        self.graph_dist_ring = graph_dist_ring
        self.parent_atom_index = parent_atom_index
        self.bond_indices: list[int] = []
        self._protein: Optional[AnalysisProtein] = None

    @property
    def residue(self) -> AnalysisResidue:
        """The residue this atom belongs to."""
        assert self._protein is not None
        return self._protein.residues[self.residue_index]

    @property
    def is_hydrogen(self) -> bool:
        return self.element == Element.H

    @property
    def parent_atom(self) -> Optional[AnalysisAtom]:
        """For hydrogens: the nearest bonded heavy atom. None for heavy atoms."""
        if self.parent_atom_index < 0 or self._protein is None:
            return None
        return self._protein.atoms[self.parent_atom_index]

    @property
    def bonds(self) -> list[AnalysisBond]:
        assert self._protein is not None
        return [self._protein.bonds[i] for i in self.bond_indices]

    # --- Time series access ---

    @property
    def positions(self) -> NDArray:
        """(T, 3) xyz positions in Angstroms, PBC-fixed."""
        assert self._protein is not None
        return self._protein._positions[:, self.index, :]

    @property
    def ring_current(self) -> AtomView:
        assert self._protein is not None and self._protein.ring_current is not None
        return AtomView(self.index, self._protein.ring_current)

    @property
    def efg(self) -> AtomView:
        assert self._protein is not None and self._protein.efg is not None
        return AtomView(self.index, self._protein.efg)

    @property
    def bond_aniso(self) -> AtomView:
        assert self._protein is not None and self._protein.bond_aniso is not None
        return AtomView(self.index, self._protein.bond_aniso)

    @property
    def quadrupole(self) -> AtomView:
        assert self._protein is not None and self._protein.quadrupole is not None
        return AtomView(self.index, self._protein.quadrupole)

    @property
    def dispersion(self) -> AtomView:
        assert self._protein is not None and self._protein.dispersion is not None
        return AtomView(self.index, self._protein.dispersion)

    @property
    def hbond(self) -> AtomView:
        assert self._protein is not None and self._protein.hbond is not None
        return AtomView(self.index, self._protein.hbond)

    @property
    def sasa(self) -> AtomView:
        assert self._protein is not None and self._protein.sasa is not None
        return AtomView(self.index, self._protein.sasa)

    @property
    def water(self) -> AtomView:
        assert self._protein is not None and self._protein.water is not None
        return AtomView(self.index, self._protein.water)

    @property
    def charges(self) -> AtomView:
        assert self._protein is not None and self._protein.charges is not None
        return AtomView(self.index, self._protein.charges)

    @property
    def aimnet2_embedding(self) -> AtomView:
        assert self._protein is not None and self._protein.aimnet2_embedding is not None
        return AtomView(self.index, self._protein.aimnet2_embedding)

    @property
    def predictions(self) -> AtomView:
        assert self._protein is not None and self._protein.predictions is not None
        return AtomView(self.index, self._protein.predictions)

    @property
    def projections(self) -> AtomView:
        assert self._protein is not None and self._protein.projections is not None
        return AtomView(self.index, self._protein.projections)

    def __repr__(self) -> str:
        res = self.residue if self._protein else None
        res_str = f"{res.name}{res.sequence_number}" if res else f"res={self.residue_index}"
        return f"AnalysisAtom({self.index}, {self.element.symbol} {self.pdb_name}, {res_str})"


class AnalysisResidue:
    """One residue in the protein. Identity only."""

    __slots__ = (
        "index", "amino_acid", "name", "sequence_number",
        "atom_indices", "_protein",
    )

    def __init__(
        self,
        index: int,
        amino_acid: AminoAcid,
        name: str,
        sequence_number: int,
    ):
        self.index = index
        self.amino_acid = amino_acid
        self.name = name
        self.sequence_number = sequence_number
        self.atom_indices: list[int] = []
        self._protein: Optional[AnalysisProtein] = None

    @property
    def atoms(self) -> list[AnalysisAtom]:
        assert self._protein is not None
        return [self._protein.atoms[i] for i in self.atom_indices]

    @property
    def is_aromatic(self) -> bool:
        return self.amino_acid.is_aromatic

    @property
    def dihedrals(self) -> Optional[ResidueView]:
        """Per-residue dihedral time series (T,) per angle."""
        if self._protein is None or self._protein.dihedrals is None:
            return None
        return ResidueView(self.index, self._protein.dihedrals)

    @property
    def dssp(self) -> Optional[ResidueView]:
        """Per-residue secondary structure time series."""
        if self._protein is None or self._protein.dssp is None:
            return None
        return ResidueView(self.index, self._protein.dssp)

    def __repr__(self) -> str:
        return f"AnalysisResidue({self.index}, {self.name} {self.sequence_number})"


class AnalysisRing:
    """
    One aromatic ring in the protein.

    Ring atom indices are stored as a flat list. The ring type determines
    the expected size (5, 6, or 9 atoms).
    """

    __slots__ = (
        "index", "type_index", "atom_indices", "parent_residue_index",
        "_protein",
    )

    def __init__(
        self,
        index: int,
        type_index: RingTypeIndex,
        atom_indices: list[int],
        parent_residue_index: int,
    ):
        self.index = index
        self.type_index = type_index
        self.atom_indices = atom_indices
        self.parent_residue_index = parent_residue_index
        self._protein: Optional[AnalysisProtein] = None

    @property
    def atoms(self) -> list[AnalysisAtom]:
        assert self._protein is not None
        return [self._protein.atoms[i] for i in self.atom_indices]

    @property
    def residue(self) -> AnalysisResidue:
        assert self._protein is not None
        return self._protein.residues[self.parent_residue_index]

    @property
    def type_name(self) -> str:
        return self.type_index.type_name

    @property
    def ring_size(self) -> int:
        return self.type_index.ring_size

    def center(self, frame: int) -> NDArray:
        """Geometric center of ring atoms at frame t. Returns (3,) array."""
        assert self._protein is not None
        positions = self._protein._positions[frame, self.atom_indices, :]
        return positions.mean(axis=0)

    def centers(self) -> NDArray:
        """Centers for all frames. Returns (T, 3) array."""
        assert self._protein is not None
        positions = self._protein._positions[:, self.atom_indices, :]
        return positions.mean(axis=1)

    def __repr__(self) -> str:
        return f"AnalysisRing({self.index}, {self.type_name}, {len(self.atom_indices)} atoms)"


class AnalysisBond:
    """One covalent bond. Topology only — geometry computed from positions."""

    __slots__ = (
        "index", "atom_index_a", "atom_index_b", "order", "category",
        "_protein",
    )

    def __init__(
        self,
        index: int,
        atom_index_a: int,
        atom_index_b: int,
        order: BondOrder,
        category: BondCategory,
    ):
        self.index = index
        self.atom_index_a = atom_index_a
        self.atom_index_b = atom_index_b
        self.order = order
        self.category = category
        self._protein: Optional[AnalysisProtein] = None

    @property
    def atom_a(self) -> AnalysisAtom:
        assert self._protein is not None
        return self._protein.atoms[self.atom_index_a]

    @property
    def atom_b(self) -> AnalysisAtom:
        assert self._protein is not None
        return self._protein.atoms[self.atom_index_b]

    @property
    def is_backbone(self) -> bool:
        return self.category in (
            BondCategory.PeptideCO,
            BondCategory.PeptideCN,
            BondCategory.BackboneOther,
        )

    @property
    def is_peptide(self) -> bool:
        return self.order == BondOrder.Peptide

    def length(self, frame: int) -> float:
        """Bond length at frame t, in Angstroms."""
        assert self._protein is not None
        pa = self._protein._positions[frame, self.atom_index_a]
        pb = self._protein._positions[frame, self.atom_index_b]
        return float(np.linalg.norm(pb - pa))

    def midpoint(self, frame: int) -> NDArray:
        """Bond midpoint at frame t. Returns (3,) array."""
        assert self._protein is not None
        pa = self._protein._positions[frame, self.atom_index_a]
        pb = self._protein._positions[frame, self.atom_index_b]
        return 0.5 * (pa + pb)

    def __repr__(self) -> str:
        return (
            f"AnalysisBond({self.index}, "
            f"{self.atom_index_a}-{self.atom_index_b}, "
            f"{self.order.name}, {self.category.name})"
        )


# ---------------------------------------------------------------------------
# Frame view — lightweight per-frame accessor
# ---------------------------------------------------------------------------

class FrameAtomView:
    """Single-frame, single-atom accessor. Returned by FrameView[atom_index]."""

    def __init__(self, frame: int, atom_index: int, protein: AnalysisProtein):
        self._t = frame
        self._i = atom_index
        self._protein = protein

    @property
    def position(self) -> NDArray:
        """(3,) xyz at this frame."""
        return self._protein._positions[self._t, self._i]

    def _scalar_or_vector(self, group: Optional[TimeSeries], name: str) -> NDArray:
        if group is None:
            raise AttributeError(f"Group not loaded")
        arr = getattr(group, name)
        return arr[self._t, self._i]

    def ring_current(self, name: str) -> NDArray:
        return self._scalar_or_vector(self._protein.ring_current, name)

    def efg(self, name: str) -> NDArray:
        return self._scalar_or_vector(self._protein.efg, name)


class FrameView:
    """
    Lightweight view of one frame across all atoms.

    Usage:
        fv = protein.frame(42)
        fv.positions           -> (N, 3) all atom positions at frame 42
        fv.ring_current("bs_T0_per_type")  -> (N, 8) all atoms at frame 42
        fv[atom_index]         -> FrameAtomView for one atom
        fv.time_ps             -> simulation time in ps
        fv.frame_index         -> original XTC frame index
    """

    def __init__(self, frame: int, protein: AnalysisProtein):
        self._t = frame
        self._protein = protein

    @property
    def time_ps(self) -> float:
        return float(self._protein.frame_times[self._t])

    @property
    def frame_index(self) -> int:
        return int(self._protein.frame_indices[self._t])

    @property
    def positions(self) -> NDArray:
        """(N, 3) atom positions at this frame."""
        return self._protein._positions[self._t]

    def _group_slice(self, group: Optional[TimeSeries], name: str) -> NDArray:
        if group is None:
            raise AttributeError(f"Group not loaded")
        return getattr(group, name)[self._t]

    def ring_current(self, name: str) -> NDArray:
        return self._group_slice(self._protein.ring_current, name)

    def efg(self, name: str) -> NDArray:
        return self._group_slice(self._protein.efg, name)

    def bond_aniso(self, name: str) -> NDArray:
        return self._group_slice(self._protein.bond_aniso, name)

    def quadrupole(self, name: str) -> NDArray:
        return self._group_slice(self._protein.quadrupole, name)

    def dispersion(self, name: str) -> NDArray:
        return self._group_slice(self._protein.dispersion, name)

    def hbond(self, name: str) -> NDArray:
        return self._group_slice(self._protein.hbond, name)

    def sasa(self, name: str) -> NDArray:
        return self._group_slice(self._protein.sasa, name)

    def water(self, name: str) -> NDArray:
        return self._group_slice(self._protein.water, name)

    def charges(self, name: str) -> NDArray:
        return self._group_slice(self._protein.charges, name)

    def predictions(self, name: str) -> NDArray:
        return self._group_slice(self._protein.predictions, name)

    def projections(self, name: str) -> NDArray:
        return self._group_slice(self._protein.projections, name)

    def __getitem__(self, atom_index: int) -> FrameAtomView:
        return FrameAtomView(self._t, atom_index, self._protein)

    def __repr__(self) -> str:
        return f"FrameView(t={self._t}, time={self.time_ps:.1f} ps, index={self.frame_index})"


# ---------------------------------------------------------------------------
# The protein
# ---------------------------------------------------------------------------

class AnalysisProtein:
    """
    Root object for an analysis trajectory H5 file.

    Built by load_analysis(). Holds:
        - Constant topology: atoms, residues, rings, bonds
        - Time series data: positions + per-group (T, N, ...) arrays
        - Metadata: protein_id, frame times/indices, dimensions

    The protein IS the molecule's identity. Time-varying data is accessed
    through atom properties that slice the shared arrays:

        protein.atoms[42].ring_current.bs_T0_per_type  -> (T, 8) numpy array

    Or through frame views:

        protein.frame(100).positions  -> (N, 3) numpy array
    """

    def __init__(
        self,
        protein_id: str,
        atoms: list[AnalysisAtom],
        residues: list[AnalysisResidue],
        rings: list[AnalysisRing],
        bonds: list[AnalysisBond],
        positions: NDArray,
        frame_times: NDArray,
        frame_indices: NDArray,
    ):
        self.protein_id = protein_id
        self.atoms = atoms
        self.residues = residues
        self.rings = rings
        self.bonds = bonds
        self._positions = positions  # (T, N, 3)
        self.frame_times = frame_times  # (T,)
        self.frame_indices = frame_indices  # (T,)

        # Per-group time series data — populated by the loader.
        # None means the group was not present in the H5.
        self.ring_current: Optional[RingCurrentData] = None
        self.efg: Optional[EfgData] = None
        self.bond_aniso: Optional[BondAnisoData] = None
        self.quadrupole: Optional[QuadrupoleData] = None
        self.dispersion: Optional[DispersionData] = None
        self.hbond: Optional[HBondData] = None
        self.sasa: Optional[SasaData] = None
        self.water: Optional[WaterData] = None
        self.charges: Optional[ChargesData] = None
        self.aimnet2_embedding: Optional[AimNet2EmbeddingData] = None
        self.predictions: Optional[PredictionsData] = None
        self.projections: Optional[ProjectionsData] = None
        self.dihedrals: Optional[DihedralData] = None
        self.dssp: Optional[DsspData] = None

        # Wire back-pointers
        for a in self.atoms:
            a._protein = self
        for r in self.residues:
            r._protein = self
        for ring in self.rings:
            ring._protein = self
        for b in self.bonds:
            b._protein = self

    # --- Dimensions ---

    @property
    def n_atoms(self) -> int:
        return len(self.atoms)

    @property
    def n_residues(self) -> int:
        return len(self.residues)

    @property
    def n_rings(self) -> int:
        return len(self.rings)

    @property
    def n_bonds(self) -> int:
        return len(self.bonds)

    @property
    def n_frames(self) -> int:
        return self._positions.shape[0]

    # --- Frame access ---

    def frame(self, t: int) -> FrameView:
        """Lightweight view of frame t."""
        if t < 0 or t >= self.n_frames:
            raise IndexError(f"Frame {t} out of range [0, {self.n_frames})")
        return FrameView(t, self)

    # --- Atom queries ---

    def atoms_by_element(self, element: Element) -> list[AnalysisAtom]:
        """All atoms of a given element."""
        return [a for a in self.atoms if a.element == element]

    def atoms_by_role(self, role: AtomRole) -> list[AnalysisAtom]:
        """All atoms with a given AtomRole."""
        return [a for a in self.atoms if a.role == role]

    def atoms_by_residue(self, residue_index: int) -> list[AnalysisAtom]:
        """All atoms in a given residue."""
        return self.residues[residue_index].atoms

    def atom_by_name(self, residue_index: int, pdb_name: str) -> Optional[AnalysisAtom]:
        """Find a specific atom by residue and PDB name (e.g., "CA")."""
        for a in self.residues[residue_index].atoms:
            if a.pdb_name == pdb_name:
                return a
        return None

    # --- Ring queries ---

    def rings_by_type(self, ring_type: RingTypeIndex) -> list[AnalysisRing]:
        """All rings of a given type."""
        return [r for r in self.rings if r.type_index == ring_type]

    # --- Positions ---

    @property
    def positions(self) -> NDArray:
        """(T, N, 3) position array."""
        return self._positions

    # --- Available groups ---

    def available_groups(self) -> list[str]:
        """Names of loaded time series groups."""
        groups = []
        for name in [
            "ring_current", "efg", "bond_aniso", "quadrupole", "dispersion",
            "hbond", "sasa", "water", "charges", "aimnet2_embedding",
            "predictions", "projections", "dihedrals", "dssp",
        ]:
            if getattr(self, name) is not None:
                groups.append(name)
        return groups

    # --- String representation ---

    def __repr__(self) -> str:
        groups = self.available_groups()
        return (
            f"AnalysisProtein('{self.protein_id}', "
            f"{self.n_atoms} atoms, {self.n_residues} residues, "
            f"{self.n_rings} rings, {self.n_bonds} bonds, "
            f"{self.n_frames} frames, "
            f"groups=[{', '.join(groups)}])"
        )

    def summary(self) -> str:
        """Human-readable summary of the loaded protein."""
        lines = [
            f"Protein: {self.protein_id}",
            f"  Atoms:    {self.n_atoms}",
            f"  Residues: {self.n_residues}",
            f"  Rings:    {self.n_rings}",
            f"  Bonds:    {self.n_bonds}",
            f"  Frames:   {self.n_frames}",
            f"  Time:     {self.frame_times[0]:.1f} - {self.frame_times[-1]:.1f} ps",
            f"  Groups:   {', '.join(self.available_groups())}",
        ]
        # Element counts
        from collections import Counter
        elem_counts = Counter(a.element for a in self.atoms)
        elem_str = ", ".join(
            f"{e.symbol}={c}" for e, c in sorted(elem_counts.items(), key=lambda x: x[0].value)
        )
        lines.append(f"  Elements: {elem_str}")
        # Ring types
        if self.rings:
            ring_types = Counter(r.type_name for r in self.rings)
            ring_str = ", ".join(f"{t}={c}" for t, c in sorted(ring_types.items()))
            lines.append(f"  Rings:    {ring_str}")
        return "\n".join(lines)

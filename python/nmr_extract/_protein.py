"""Protein dataclass and loader.

    from nmr_extract import load
    p = load("path/to/extraction")
    p.biot_savart.shielding.T2      # (N, 5)
    p.ring_contributions.bs.T2      # (P, 5)
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

from ._tensors import (
    ShieldingTensor,
    EFGTensor,
    VectorField,
    MagneticVectorField,
    PositionField,
    PerRingTypeT0,
    PerRingTypeT1,
    PerRingTypeT2,
    RingCounts,
    McConnellNearFieldCounts,
    CoulombScalars,
    HBondScalars,
    DsspScalars,
    MopacScalars,
    MopacGlobal,
    MopacAtomPopulations,
    MopacAtomicOrbitalPopulations,
    MopacAtomicOrbitalPopulationTotals,
    MopacUniqueBondOrders,
    MopacTopologyBondOrdersFull,
    BondOrders,
    DeltaScalars,
    DeltaAPBS,
    DeltaRingProximity,
    AIMNet2Charges,
    AIMNet2AimEmbedding,
    AIMNet2ChargeResponseGradient,
)
from ._ring import RingContributions, RingGeometry, RingPairGeometry
from ._catalog import CATALOG


# ── Calculator group dataclasses ────────────────────────────────────


@dataclass(frozen=True)
class RingKernelGroup:
    """Packed ring-current kernel plus its per-type decomposition.

    The wrapper name does not imply ppm: consult each array's catalog units
    and ``scaling_contract`` before treating a native kernel as shielding.
    """
    shielding: ShieldingTensor
    per_type_T0: PerRingTypeT0
    per_type_T2: PerRingTypeT2
    per_type_T1: Optional[PerRingTypeT1] = None


@dataclass(frozen=True)
class BiotSavartGroup(RingKernelGroup):
    """Biot-Savart with additional B-field and ring counts."""
    total_B: MagneticVectorField = None
    ring_B_field: MagneticVectorField = None
    ring_B_cylindrical: Optional[np.ndarray] = None
    ring_counts: RingCounts = None


@dataclass(frozen=True)
class HaighMallionGroup(RingKernelGroup):
    """Unscaled Å⁻¹ Haigh-Mallion kernels and effective-field diagnostics.

    The producer applies no ring-current intensity or ppm conversion; the
    physical/model scale is downstream.
    """
    ring_B_field: MagneticVectorField = None


@dataclass(frozen=True)
class PiQuadrupoleGroup:
    """Geometry-only pi-quadrupole kernels; physical scale is learnable.

    The producer applies no quadrupole prefactor. ``local_tensor`` and
    ``local_T2`` are raw project-basis geometry tensors, not e3nn tensors.
    """
    per_type_T0: PerRingTypeT0
    quad_scalar: Optional[np.ndarray] = None
    axial_scalar_per_type_T0: Optional[PerRingTypeT0] = None
    local_tensor: Optional[np.ndarray] = None
    local_T2: Optional[np.ndarray] = None
    local_frame: Optional[np.ndarray] = None
    local_geometry: Optional[np.ndarray] = None

    @property
    def pq_per_type_T0(self) -> PerRingTypeT0:
        """Backward-compatible name for the axial geometry scalar."""
        return self.per_type_T0

    @property
    def piquad_axial_scalar_per_type_T0(self) -> PerRingTypeT0:
        """Preferred honest name, with fallback for pre-alias extractions."""
        return (self.axial_scalar_per_type_T0
                if self.axial_scalar_per_type_T0 is not None
                else self.per_type_T0)


@dataclass(frozen=True)
class RingSusceptibilityGroup:
    """Ring-susceptibility geometry scalar on union ring-neighbour rows.

    A sparse zero has no dedicated calculator mask and is ambiguous between
    a physical zero and a row not accepted by RingSusceptibilityResult.
    """
    scalar: np.ndarray
    per_type_T0: PerRingTypeT0


@dataclass(frozen=True)
class DispersionGroup:
    """Dispersion scalar decomposition by ring type."""
    per_type_T0: PerRingTypeT0
    aromatic_r6_proximity_per_type_T0: Optional[PerRingTypeT0] = None


@dataclass(frozen=True)
class EnrichmentGroup:
    """Per-atom enrichment classifications from EnrichmentResult.

    ``parent_is_sp2`` is a compatibility heuristic (aromatic parent or cached
    backbone C/N), not a general typed-hybridisation query. Flag column 6,
    exposed by ``is_hbond_acceptor``, is only the coarse producer predicate
    "element is N or O"; use ``acceptor_class`` for typed chemistry.
    """
    role: np.ndarray
    hybridisation: np.ndarray
    flags: np.ndarray
    parent_is_sp2: Optional[np.ndarray] = None
    polar_h_kind: Optional[np.ndarray] = None
    planar_group_kind: Optional[np.ndarray] = None
    formal_charge: Optional[np.ndarray] = None
    ring_position: Optional[np.ndarray] = None
    locant: Optional[np.ndarray] = None
    donor_class: Optional[np.ndarray] = None
    acceptor_class: Optional[np.ndarray] = None
    hybridisation_class: Optional[np.ndarray] = None

    @property
    def is_backbone(self) -> np.ndarray:
        return self.flags[:, 0] != 0

    @property
    def is_amide_H(self) -> np.ndarray:
        return self.flags[:, 1] != 0

    @property
    def is_alpha_H(self) -> np.ndarray:
        return self.flags[:, 2] != 0

    @property
    def is_methyl(self) -> np.ndarray:
        return self.flags[:, 3] != 0

    @property
    def is_aromatic_H(self) -> np.ndarray:
        return self.flags[:, 4] != 0

    @property
    def is_hbond_donor(self) -> np.ndarray:
        return self.flags[:, 5] != 0

    @property
    def is_hbond_acceptor(self) -> np.ndarray:
        """Coarse N/O candidate flag, not a chemically typed acceptor mask."""
        return self.flags[:, 6] != 0

    @property
    def is_on_aromatic_residue(self) -> np.ndarray:
        return self.flags[:, 7] != 0


@dataclass(frozen=True)
class ChargeAssignmentGroup:
    """Configured typed ChargeSource table projected onto atom rows.

    The source may be a parameter file, AMBER/runtime-tleap prmtop, or a
    caller-preloaded table; this group does not imply a single force field.
    """
    partial_charge: np.ndarray
    pb_radius: np.ndarray


@dataclass(frozen=True)
class SpatialIndexGroup:
    """Sparse atom-neighbor rows from SpatialIndexResult."""
    neighbors: np.ndarray


@dataclass(frozen=True)
class MolecularGraphGroup:
    """Per-atom typed through-bond graph features."""
    integer: np.ndarray
    floating: np.ndarray
    compatibility: np.ndarray


@dataclass(frozen=True)
class McConnellGroup:
    """Two-channel McConnell source response by source category.

    Values are unscaled Å⁻³ unit-susceptibility geometry responses, not ppm.
    The physical response scale is downstream; the extraction manifest pins
    peptide-C=O shape ratios but declares the axial magnitude learned.

    Each field is a project-native packed SphericalTensor view over one
    ``mc_<cat>_<ch>.npy`` array with component order
    ``T0,T1_x,T1_y,T1_z,T2_m-2..+2``. Convert explicitly before e3nn use.
    """
    peptide_co_fixed: ShieldingTensor
    peptide_co_bo: ShieldingTensor
    peptide_co_rhombic: ShieldingTensor
    peptide_cn_fixed: ShieldingTensor
    peptide_cn_bo: ShieldingTensor
    backbone_other_fixed: ShieldingTensor
    backbone_other_bo: ShieldingTensor
    sidechain_co_fixed: ShieldingTensor
    sidechain_co_bo: ShieldingTensor
    sidechain_other_fixed: ShieldingTensor
    sidechain_other_bo: ShieldingTensor
    disulfide_fixed: ShieldingTensor
    disulfide_bo: ShieldingTensor
    aromatic_zeroed_fixed: ShieldingTensor
    aromatic_zeroed_bo: ShieldingTensor
    backbone_xh_fixed: ShieldingTensor
    backbone_xh_bo: ShieldingTensor
    sidechain_xh_fixed: ShieldingTensor
    sidechain_xh_bo: ShieldingTensor
    s_h_fixed: ShieldingTensor
    s_h_bo: ShieldingTensor
    nearfield_counts: Optional[McConnellNearFieldCounts] = None
    nearest_co_dir: Optional[VectorField] = None
    nearest_co_midpoint: Optional[PositionField] = None
    nearest_co_T2: Optional[ShieldingTensor] = None
    nearest_cn_T2: Optional[ShieldingTensor] = None
    bond_neighbors: Optional[np.ndarray] = None

    @property
    def fixed(self) -> dict[str, ShieldingTensor]:
        return {
            "peptide_co": self.peptide_co_fixed,
            "peptide_cn": self.peptide_cn_fixed,
            "backbone_other": self.backbone_other_fixed,
            "sidechain_co": self.sidechain_co_fixed,
            "sidechain_other": self.sidechain_other_fixed,
            "disulfide": self.disulfide_fixed,
            "aromatic_zeroed": self.aromatic_zeroed_fixed,
            "backbone_xh": self.backbone_xh_fixed,
            "sidechain_xh": self.sidechain_xh_fixed,
            "s_h": self.s_h_fixed,
        }

    @property
    def bo(self) -> dict[str, ShieldingTensor]:
        return {
            "peptide_co": self.peptide_co_bo,
            "peptide_cn": self.peptide_cn_bo,
            "backbone_other": self.backbone_other_bo,
            "sidechain_co": self.sidechain_co_bo,
            "sidechain_other": self.sidechain_other_bo,
            "disulfide": self.disulfide_bo,
            "aromatic_zeroed": self.aromatic_zeroed_bo,
            "backbone_xh": self.backbone_xh_bo,
            "sidechain_xh": self.sidechain_xh_bo,
            "s_h": self.s_h_bo,
        }


@dataclass(frozen=True)
class SidechainCarbonylAnisotropyGroup:
    """Typed side-chain C=O sources, local frames, and tensor audits."""
    source_bonds: np.ndarray
    frame: np.ndarray
    frame_quality: np.ndarray
    fixed_T2: ShieldingTensor
    bo_T2: ShieldingTensor
    scalar_audit: np.ndarray


@dataclass(frozen=True)
class CoulombGroup:
    """Coulomb E-field and bare EFG decompositions.

    E-field wrappers expose only a conditional polar irrep because the
    """
    efg: ShieldingTensor
    E: VectorField
    E_backbone: VectorField
    E_sidechain: VectorField
    E_aromatic: VectorField
    efg_backbone: EFGTensor
    efg_sidechain: EFGTensor
    efg_aromatic: EFGTensor
    scalars: CoulombScalars
    efg_t2: Optional[EFGTensor] = None
    aromatic_E_proj: Optional[np.ndarray] = None
    aromatic_n_src: Optional[np.ndarray] = None
    E_solvent: Optional[VectorField] = None
    efg_solvent: Optional[EFGTensor] = None


@dataclass(frozen=True)
class EeqCoulombGroup:
    """Coulomb E-field and bare EFG decompositions from EEQ charges."""
    efg: ShieldingTensor
    E: VectorField
    E_backbone: VectorField
    E_sidechain: VectorField
    E_aromatic: VectorField
    efg_backbone: EFGTensor
    efg_sidechain: EFGTensor
    efg_aromatic: EFGTensor
    scalars: CoulombScalars
    aromatic_E_proj: Optional[np.ndarray] = None
    aromatic_n_src: Optional[np.ndarray] = None


@dataclass(frozen=True)
class HBondGroup:
    scalars: HBondScalars
    nearest_dir: Optional[VectorField] = None
    flags: Optional[np.ndarray] = None

    @property
    def is_backbone(self) -> Optional[np.ndarray]:
        if self.flags is None:
            return None
        return self.flags[:, 0] != 0

    @property
    def is_donor(self) -> Optional[np.ndarray]:
        if self.flags is None:
            return None
        return self.flags[:, 1] != 0

    @property
    def is_acceptor(self) -> Optional[np.ndarray]:
        if self.flags is None:
            return None
        return self.flags[:, 2] != 0


@dataclass(frozen=True)
class MopacCoreGroup:
    """Legacy-quantized compatibility projection of diskless libmopac data."""
    charges: np.ndarray
    scalars: MopacScalars
    bond_orders: BondOrders
    global_: MopacGlobal
    bond_neighbors: Optional[np.ndarray] = None


@dataclass(frozen=True)
class MopacFullGroup:
    """Additional direct libmopac arrays and compatibility diagnostics.

    This is not every field of the API structs: control-plane failures and
    structural-N/A job outputs are intentionally not scientific NPYs. Core
    compatibility charge/bond/global arrays live in sibling ``MopacCoreGroup``.
    Field-level ``Optional`` annotations preserve loading of older/partial SDK
    schemas; they are not current scientific absence states. A successful
    current ``MopacResult::WriteFeatures`` attempts the complete declared
    core/direct family, while worker/API/dimension/non-finite failure prevents
    the whole MopacResult from attaching rather than selectively zero-filling
    or omitting one direct feature.
    """

    charges_full_precision: Optional[np.ndarray] = None
    bond_orders_full_precision: Optional[BondOrders] = None
    bond_valencies_full_precision: Optional[np.ndarray] = None
    atom_populations: Optional[MopacAtomPopulations] = None
    atomic_orbital_populations: Optional[MopacAtomicOrbitalPopulations] = None
    atomic_orbital_population_totals: Optional[MopacAtomicOrbitalPopulationTotals] = None
    bond_valencies: Optional[np.ndarray] = None
    bond_orders_unique: Optional[MopacUniqueBondOrders] = None
    topology_bond_orders_full: Optional[MopacTopologyBondOrdersFull] = None
    heat_kcal_mol: Optional[np.ndarray] = None
    dipole_debye: Optional[np.ndarray] = None
    dipole_point_charge_debye: Optional[np.ndarray] = None
    dipole_hybridization_debye: Optional[np.ndarray] = None
    bond_index: Optional[np.ndarray] = None
    bond_atom: Optional[np.ndarray] = None
    bond_order: Optional[np.ndarray] = None
    ao_max_orbitals: Optional[np.ndarray] = None
    ao_orbitals_per_atom: Optional[np.ndarray] = None
    atom_ao_density: Optional[np.ndarray] = None
    atomic_orbital_populations_full_precision: Optional[np.ndarray] = None
    bond_ao_density_directed: Optional[np.ndarray] = None
    bond_density_pairs: Optional[np.ndarray] = None
    bond_ao_density: Optional[np.ndarray] = None
    atom_electron_population: Optional[np.ndarray] = None
    atom_s_population: Optional[np.ndarray] = None
    atom_p_population: Optional[np.ndarray] = None
    atom_d_population: Optional[np.ndarray] = None
    lewis_bond_count: Optional[np.ndarray] = None
    lewis_bond_atoms: Optional[np.ndarray] = None
    lmo_energy_levels: Optional[np.ndarray] = None
    lmo_occupied_atom_counts: Optional[np.ndarray] = None
    lmo_occupied_atoms: Optional[np.ndarray] = None
    lmo_occupied_coefficients: Optional[np.ndarray] = None
    lmo_virtual_atom_counts: Optional[np.ndarray] = None
    lmo_virtual_atoms: Optional[np.ndarray] = None
    lmo_virtual_coefficients: Optional[np.ndarray] = None
    lmo_occupied_atom_offsets_native: Optional[np.ndarray] = None
    lmo_virtual_atom_offsets_native: Optional[np.ndarray] = None
    lmo_occupied_coefficient_offsets_native: Optional[np.ndarray] = None
    lmo_virtual_coefficient_offsets_native: Optional[np.ndarray] = None
    lmo_occupied_atom_storage_native: Optional[np.ndarray] = None
    lmo_virtual_atom_storage_native: Optional[np.ndarray] = None
    lmo_occupied_coefficient_storage_native: Optional[np.ndarray] = None
    lmo_virtual_coefficient_storage_native: Optional[np.ndarray] = None
    mozyme_state_dimensions: Optional[np.ndarray] = None


@dataclass(frozen=True)
class MopacCoulombGroup:
    """All-pairs fields from legacy F15.6 diskless-libmopac Coulson charges.

    Cartesian E has a polar law only for fixed charges on a non-exceptional
    """
    efg: ShieldingTensor
    E: VectorField
    E_backbone: VectorField
    E_sidechain: VectorField
    E_aromatic: VectorField
    efg_backbone: EFGTensor
    efg_sidechain: EFGTensor
    efg_aromatic: EFGTensor
    scalars: CoulombScalars
    aromatic_E_proj: Optional[np.ndarray] = None
    aromatic_n_src: Optional[np.ndarray] = None
    E_clamp_mask: Optional[np.ndarray] = None
    E_clamp_scale: Optional[np.ndarray] = None


@dataclass(frozen=True)
class MopacMcConnellGroup:
    """Direct read-back of the MOPAC Wiberg-weighted McConnell projection.

    Tensor fields are unscaled Å⁻³ project-native Full9 values, not ppm.
    The two nearest tensors have a source only where their corresponding
    distance is below the producer's ``NO_DATA_SENTINEL``.
    """
    co_sum: np.ndarray
    cn_sum: np.ndarray
    sidechain_sum: np.ndarray
    aromatic_sum: np.ndarray
    co_nearest: np.ndarray
    nearest_co_dist: np.ndarray
    nearest_cn_dist: np.ndarray
    nearest_co_T2: ShieldingTensor
    nearest_cn_T2: ShieldingTensor
    backbone_total: ShieldingTensor
    sidechain_total: ShieldingTensor
    aromatic_total: ShieldingTensor
    shielding: ShieldingTensor


@dataclass(frozen=True)
class MopacGroup:
    core: MopacCoreGroup
    coulomb: MopacCoulombGroup
    full: Optional[MopacFullGroup] = None
    mcconnell: Optional[MopacMcConnellGroup] = None


@dataclass(frozen=True)
class APBSGroup:
    """APBS reaction-field outputs from the finite-difference PB solve."""
    E: VectorField
    efg: EFGTensor
    phi: Optional[np.ndarray] = None
    E_clamp_mask: Optional[np.ndarray] = None
    E_clamp_scale: Optional[np.ndarray] = None
    nonfinite_sanitizer_mask: Optional[np.ndarray] = None
    E_total_diagnostic: Optional[VectorField] = None
    efg_total_diagnostic: Optional[EFGTensor] = None


@dataclass(frozen=True)
class OrcaGroup:
    total: ShieldingTensor
    diamagnetic: ShieldingTensor
    paramagnetic: ShieldingTensor


@dataclass(frozen=True)
class DeltaGroup:
    """WT - mutant deltas, indexed by WT atom row.

    `shielding` is the total DFT shielding delta (sigma_total = sigma_dia
    + sigma_para). The dia/para components are also stored on each side
    AND as deltas — six tensors total — so analyses can stratify the
    mutation shift by physical mechanism (electron-density change vs.
    orbital-response change) without joining against the per-conformation
    ORCA NPYs in two extraction directories.

    All component tensors are optional; they are present when an ORCA
    NMR output file accompanied both the WT and mutant runs.
    """
    shielding: ShieldingTensor
    scalars: DeltaScalars
    apbs: Optional[DeltaAPBS]
    ring_proximity: DeltaRingProximity
    graph: Optional[np.ndarray] = None

    # WT side
    wt_shielding_diamagnetic: Optional[ShieldingTensor] = None
    wt_shielding_paramagnetic: Optional[ShieldingTensor] = None
    # mut side
    mut_shielding_diamagnetic: Optional[ShieldingTensor] = None
    mut_shielding_paramagnetic: Optional[ShieldingTensor] = None
    # WT - mut for each component
    delta_shielding_diamagnetic: Optional[ShieldingTensor] = None
    delta_shielding_paramagnetic: Optional[ShieldingTensor] = None


# ── Per-atom invariant categorical record ───────────────────────────


class CategoryInfo:
    """Per-atom invariant categorical record from ``atoms_category_info.npy``.

    One structured-dtype row per atom, populated by ``CategoryInfoProjection``
    at one-shot per-protein emission time. Holds every categorical fact
    the C++ side knows about each atom: identity (atom_index, residue_index,
    element), atom names across naming systems (AMBER / IUPAC / BMRB),
    per-residue 3-letter / 1-letter codes (also AMBER / IUPAC / BMRB),
    mechanical identity (locant, branch, di_index, backbone_role), and
    chemistry classification (prochiral, planar group/stereo, polar-H
    kind, ring position, pseudoatom membership, aromatic, formal_charge,
    is_exchangeable). See ``src/CategoryInfoProjection.h`` and
    ``src/SemanticEnums.h`` for the enum value mappings.

    Convenience properties decode the most-used columns; the rest are
    accessible as raw int8 arrays via ``info.data[<field>]``. Do
    stratification with numpy boolean masks:

        info.is_backbone & (info.element == 7)            # backbone N
        info.is_aromatic & (info.element == 6)            # aromatic C
        info.amber_residue_3letter == b'CYX'              # disulfide CYS
        info.data['polar_h_kind'] != 0                    # any polar H
    """

    __slots__ = ("_data",)

    def __init__(self, data: np.ndarray):
        if data.dtype.fields is None:
            raise ValueError(
                "CategoryInfo: expected a numpy structured array; got "
                f"flat dtype {data.dtype}. The C++ writer emits a "
                "structured-dtype NPY; if this fires, the file is "
                "either truncated or the schema diverged from the "
                "loader.")
        self._data = data

    # ── Raw access ──
    @property
    def data(self) -> np.ndarray:
        """Raw numpy structured array. All fields per ``data.dtype.names``."""
        return self._data

    def __len__(self) -> int:
        return len(self._data)

    @property
    def n_atoms(self) -> int:
        return len(self._data)

    # ── Decoded string columns ──
    @property
    def amber_atom_name(self) -> np.ndarray:
        """Per-atom AMBER ff14SB atom name as ``bytes`` (S8)."""
        return self._data["amber_atom_name"]

    @property
    def iupac_atom_name(self) -> np.ndarray:
        return self._data["iupac_atom_name"]

    @property
    def bmrb_atom_name(self) -> np.ndarray:
        return self._data["bmrb_atom_name"]

    @property
    def amber_residue_3letter(self) -> np.ndarray:
        """Variant-specific AMBER residue name (CYX / HID / ASH / ...)."""
        return self._data["amber_residue_3letter"]

    @property
    def iupac_residue_3letter(self) -> np.ndarray:
        """Canonical IUPAC residue name (CYS / HIS / ASP / ...)."""
        return self._data["iupac_residue_3letter"]

    @property
    def bmrb_residue_3letter(self) -> np.ndarray:
        """Canonical BMRB residue name (== IUPAC for the standard 20)."""
        return self._data["bmrb_residue_3letter"]

    @property
    def residue_1letter(self) -> np.ndarray:
        return self._data["residue_1letter"]

    # ── Numeric columns commonly used in stratification ──
    @property
    def atom_index(self) -> np.ndarray:
        return self._data["atom_index"]

    @property
    def residue_index(self) -> np.ndarray:
        return self._data["residue_index"]

    @property
    def element(self) -> np.ndarray:
        """Per-atom atomic number as int8. 1=H, 6=C, 7=N, 8=O, 16=S.

        Matches the convention of the existing ``element.npy``; codex
        can stratify by element with ``info.element == 6`` etc.
        """
        return self._data["element"]

    @property
    def residue_type(self) -> np.ndarray:
        """Per-atom AminoAcid enum index (canonical, NOT variant-aware)."""
        return self._data["residue_type"]

    @property
    def residue_variant_index(self) -> np.ndarray:
        """AMBER protonation variant index per residue. -1 = canonical."""
        return self._data["residue_variant_index"]

    @property
    def terminal_state(self) -> np.ndarray:
        """ResidueTerminalState int: Internal/Nterm{Charged,Neutral}/Cterm{...}."""
        return self._data["terminal_state"]

    @property
    def formal_charge(self) -> np.ndarray:
        return self._data["formal_charge"]

    # ── Convenience boolean masks ──
    @property
    def is_backbone(self) -> np.ndarray:
        """``backbone_role != None`` — N/CA/C/O/H/HA backbone slots."""
        return self._data["backbone_role"] != 0

    @property
    def is_aromatic(self) -> np.ndarray:
        return self._data["aromatic"] != 0

    @property
    def is_exchangeable(self) -> np.ndarray:
        return self._data["is_exchangeable"] != 0

    @property
    def is_polar_h(self) -> np.ndarray:
        """``polar_h_kind != NotPolar`` — atoms classified as exchangeable Hs."""
        return self._data["polar_h_kind"] != 0

    @property
    def in_super_group(self) -> np.ndarray:
        """In a Markley super-aggregate (QG/QD/QH/QR)."""
        return self._data["in_super_group"] != 0

    # ── Topology sidecar fields (2026-05-13 extension) ──
    # Direct dtype access. An old NPY without these columns produces
    # KeyError at access -- the model is canonical, the file is wrong.

    @property
    def chain_id(self) -> np.ndarray:
        """Per-atom chain identifier (``|S2`` bytes; empty if not multi-chain).

        Reflects ``Residue::chain_id`` (Residue.h:40).
        """
        return self._data["chain_id"]

    @property
    def residue_number(self) -> np.ndarray:
        """Per-atom biological residue number (PDB sequence_number).

        Reflects ``Residue::sequence_number`` (Residue.h:39). Independent
        from the dense ``residue_index`` axis -- this is the source-level
        biological numbering that BMRB/RefDB joins key on.
        """
        return self._data["residue_number"]

    @property
    def insertion_code(self) -> np.ndarray:
        """Per-atom PDB insertion code (``|S1``; empty when absent).

        Reflects ``Residue::insertion_code`` (Residue.h:41).
        """
        return self._data["insertion_code"]

    @property
    def parent_atom_index(self) -> np.ndarray:
        """For element 1 (H) atoms, the bonded heavy atom's index.

        All other elements read ``-1``. Reflects
        ``Atom::parent_atom_index`` (Atom.h:29); ``SIZE_MAX`` is mapped
        to ``-1`` at projection time. Hydrogens whose parent wasn't
        assigned (rare; happens for isolated atoms) also read ``-1`` —
        cross-check ``element == 1`` to distinguish unassigned-H from
        heavy-atom-by-convention.
        """
        return self._data["parent_atom_index"]

    @property
    def ff_atom_type_string(self) -> np.ndarray:
        """Per-atom AMBER ff14SB atom type string (``|S4``; empty if no FF data).

        Reflects ``LegacyAmberTopology::AtomtypeString(ai)``
        (LegacyAmberTopology.h:158). Populated when the load path
        supplied FF data (PRMTOP / GROMACS readback); empty for
        PDB-only loads per LegacyAmberInvariants's empty-vector
        convention.
        """
        return self._data["ff_atom_type_string"]

    @property
    def equivalence_class(self) -> np.ndarray:
        """Per-atom RDKit canonical-rank equivalence class.

        Reflects ``AtomSemanticTable::equivalence_class``
        (SemanticEnums.h:877). Populated by the substrate generator
        (``tools/topology/build_semantic_tables.cpp``) from RDKit
        canonical_rank. Zero is the unassigned sentinel; consumers
        should mask via ``> 0`` for equivalence-class stratification.
        """
        return self._data["equivalence_class"]


# ─── Topology-sidecar additive projections ──────────────────────────
#
# Bonds / Rings / RingMembership / ExtractionManifest.
#
# Wrappers over the structured-NPY / JSON sibling artifacts emitted by
# ``src/TopologySidecar.cpp``. Same architectural shape as
# ``CategoryInfo``: thin views on a structured array (or parsed JSON)
# with convenience accessors; no model is built.
#

class Residues:
    """Per-residue record from ``residues.npy``.

    One row per residue. Codex's RequiredTable contract item, plus
    prev/next chain-aware links (-1 at chain boundaries) and the
    Markley-style 1-letter / 3-letter renderings.

    Dtype: ``residue_index`` ``chain_id`` (S2) ``residue_number``
    ``insertion_code`` (S1) ``residue_type`` (AminoAcid enum)
    ``amber_residue_3letter`` ``iupac_residue_3letter`` ``one_letter``
    ``protonation_variant_index`` ``terminal_state``
    ``prev_residue_index`` ``next_residue_index`` ``prev_residue_type``
    ``next_residue_type`` ``atom_count`` ``is_proline`` ``is_aromatic``
    ``is_titratable`` ``has_amide_h``.
    """

    __slots__ = ("_data",)

    def __init__(self, data: np.ndarray):
        if data.dtype.fields is None:
            raise ValueError(
                "Residues: expected a numpy structured array; got "
                f"flat dtype {data.dtype}.")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    def __len__(self) -> int:
        return len(self._data)

    @property
    def n_residues(self) -> int:
        return len(self._data)

    @property
    def residue_index(self) -> np.ndarray:
        return self._data["residue_index"]

    @property
    def chain_id(self) -> np.ndarray:
        return self._data["chain_id"]

    @property
    def residue_number(self) -> np.ndarray:
        return self._data["residue_number"]

    @property
    def insertion_code(self) -> np.ndarray:
        return self._data["insertion_code"]

    @property
    def residue_type(self) -> np.ndarray:
        return self._data["residue_type"]

    @property
    def amber_residue_3letter(self) -> np.ndarray:
        return self._data["amber_residue_3letter"]

    @property
    def iupac_residue_3letter(self) -> np.ndarray:
        return self._data["iupac_residue_3letter"]

    @property
    def one_letter(self) -> np.ndarray:
        return self._data["one_letter"]

    @property
    def protonation_variant_index(self) -> np.ndarray:
        return self._data["protonation_variant_index"]

    @property
    def terminal_state(self) -> np.ndarray:
        return self._data["terminal_state"]

    @property
    def prev_residue_index(self) -> np.ndarray:
        return self._data["prev_residue_index"]

    @property
    def next_residue_index(self) -> np.ndarray:
        return self._data["next_residue_index"]

    @property
    def prev_residue_type(self) -> np.ndarray:
        return self._data["prev_residue_type"]

    @property
    def next_residue_type(self) -> np.ndarray:
        return self._data["next_residue_type"]

    @property
    def atom_count(self) -> np.ndarray:
        return self._data["atom_count"]

    @property
    def is_proline(self) -> np.ndarray:
        return self._data["is_proline"] != 0

    @property
    def is_aromatic(self) -> np.ndarray:
        return self._data["is_aromatic"] != 0

    @property
    def is_titratable(self) -> np.ndarray:
        return self._data["is_titratable"] != 0

    @property
    def has_amide_h(self) -> np.ndarray:
        return self._data["has_amide_h"] != 0

    @property
    def is_xpro_context(self) -> np.ndarray:
        """True when residue ``i+1`` is PRO (X→Pro peptide context).

        Relevant for the ``(i, i+1)`` ω dihedral — X→Pro permits cis
        isomerism that standard non-Pro context does not. False at
        chain ends (``next_residue_index == -1``).
        """
        return self._data["is_xpro_context"] != 0


class Bonds:
    """Per-bond record from ``bonds.npy``.

    One row per bond from ``LegacyAmberTopology::BondList()``. The
    structured dtype carries: ``bond_index``, ``atom_index_a``,
    ``atom_index_b``, ``bond_order`` (BondOrder enum), ``bond_category``
    (BondCategory enum), ``is_rotatable``, ``is_aromatic``,
    ``is_peptide``, ``is_backbone``.

    See ``src/TopologySidecar.cpp`` for the canonical dtype declaration.
    """

    __slots__ = ("_data",)

    def __init__(self, data: np.ndarray):
        if data.dtype.fields is None:
            raise ValueError(
                "Bonds: expected a numpy structured array; got "
                f"flat dtype {data.dtype}.")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    def __len__(self) -> int:
        return len(self._data)

    @property
    def n_bonds(self) -> int:
        return len(self._data)

    @property
    def bond_index(self) -> np.ndarray:
        return self._data["bond_index"]

    @property
    def atom_index_a(self) -> np.ndarray:
        return self._data["atom_index_a"]

    @property
    def atom_index_b(self) -> np.ndarray:
        return self._data["atom_index_b"]

    @property
    def bond_order(self) -> np.ndarray:
        """BondOrder enum: 0=Single 1=Double 2=Triple 3=Aromatic 4=Peptide 5=Unknown."""
        return self._data["bond_order"]

    @property
    def bond_category(self) -> np.ndarray:
        """BondCategory enum: 0=PeptideCO 1=PeptideCN 2=BackboneOther 3=SidechainCO 4=Aromatic 5=Disulfide 6=SidechainOther 7=Unknown."""
        return self._data["bond_category"]

    @property
    def is_rotatable(self) -> np.ndarray:
        return self._data["is_rotatable"] != 0

    @property
    def is_aromatic(self) -> np.ndarray:
        return self._data["is_aromatic"] != 0

    @property
    def is_peptide(self) -> np.ndarray:
        return self._data["is_peptide"] != 0

    @property
    def is_backbone(self) -> np.ndarray:
        return self._data["is_backbone"] != 0


class Rings:
    """Per-ring record from ``rings.npy``.

    One row per ring from ``LegacyAmberTopology::Rings()``. Aromatic
    rings come first (rows ``0..aromatic_ring_count-1``), then
    saturated rings. ``ring_id`` is the absolute row index;
    ``native_axis_index`` is the index within the aromatic-only or
    saturated-only axis (matches ``ring_geometry.npy`` row order for
    aromatic).

    Dtype: ``ring_id`` ``ring_kind`` (0=aromatic, 1=saturated)
    ``ring_type_index`` (RingTypeIndex enum) ``atom_count``
    ``native_axis_index`` ``parent_residue_index``
    ``parent_residue_number`` ``fused_partner_ring_id`` (-1 if not fused).
    """

    __slots__ = ("_data",)

    def __init__(self, data: np.ndarray):
        if data.dtype.fields is None:
            raise ValueError(
                "Rings: expected a numpy structured array; got "
                f"flat dtype {data.dtype}.")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    def __len__(self) -> int:
        return len(self._data)

    @property
    def n_rings(self) -> int:
        return len(self._data)

    @property
    def ring_id(self) -> np.ndarray:
        return self._data["ring_id"]

    @property
    def ring_kind(self) -> np.ndarray:
        """0 = aromatic, 1 = saturated."""
        return self._data["ring_kind"]

    @property
    def is_aromatic(self) -> np.ndarray:
        return self._data["ring_kind"] == 0

    @property
    def is_saturated(self) -> np.ndarray:
        return self._data["ring_kind"] == 1

    @property
    def ring_type_index(self) -> np.ndarray:
        """RingTypeIndex enum: 0=PHE, 1=TYR, 2=TRPbenzene, ..., 8=ProPyrrolidine."""
        return self._data["ring_type_index"]

    @property
    def atom_count(self) -> np.ndarray:
        return self._data["atom_count"]

    @property
    def native_axis_index(self) -> np.ndarray:
        """Index within the aromatic-only or saturated-only axis."""
        return self._data["native_axis_index"]

    @property
    def parent_residue_index(self) -> np.ndarray:
        return self._data["parent_residue_index"]

    @property
    def parent_residue_number(self) -> np.ndarray:
        return self._data["parent_residue_number"]

    @property
    def fused_partner_ring_id(self) -> np.ndarray:
        """Absolute ring_id of the fused partner. -1 if not fused."""
        return self._data["fused_partner_ring_id"]


class RingMembership:
    """Per (ring, ring-vertex-atom) record from ``ring_membership.npy``.

    One row per (ring, vertex atom) pair, in canonical ring walk order.
    Codex's RingMembershipTable contract -- the only acceptable basis
    for ring-to-atom projection of pucker / aromatic-chi2 / ring-current
    contributions.

    Dtype: ``ring_id`` ``atom_index`` ``ring_atom_order`` ``is_vertex``
    ``is_substituent`` (currently always 0 -- reserved for future
    ring-substituent extraction).
    """

    __slots__ = ("_data",)

    def __init__(self, data: np.ndarray):
        if data.dtype.fields is None:
            raise ValueError(
                "RingMembership: expected a numpy structured array; got "
                f"flat dtype {data.dtype}.")
        self._data = data

    @property
    def data(self) -> np.ndarray:
        return self._data

    def __len__(self) -> int:
        return len(self._data)

    @property
    def n_rows(self) -> int:
        return len(self._data)

    @property
    def ring_id(self) -> np.ndarray:
        return self._data["ring_id"]

    @property
    def atom_index(self) -> np.ndarray:
        return self._data["atom_index"]

    @property
    def ring_atom_order(self) -> np.ndarray:
        return self._data["ring_atom_order"]

    @property
    def is_vertex(self) -> np.ndarray:
        return self._data["is_vertex"] != 0


class ExtractionManifest:
    """Parsed ``extraction_manifest.json`` sibling.

    Holds schema version, protein id, generated_at, topology-population
    flags, axis sizes, and the axis-alignment statements. Use as a
    self-describing companion to the NPY tree; codex's first-pass
    validation gates read from here.
    """

    __slots__ = ("_data",)

    def __init__(self, parsed: dict):
        self._data = parsed

    @property
    def data(self) -> dict:
        return self._data

    @property
    def schema_version(self) -> str:
        return self._data.get("schema_version", "")

    @property
    def protein_id(self) -> str:
        return self._data.get("protein_id", "")

    @property
    def generated_at_utc(self) -> str:
        return self._data.get("generated_at_utc", "")

    @property
    def topology(self) -> dict:
        return self._data.get("topology", {})

    @property
    def axis_sizes(self) -> dict:
        return self._data.get("axis_sizes", {})

    @property
    def axis_alignment(self) -> dict:
        return self._data.get("axis_alignment", {})

    @property
    def feature_metadata(self) -> dict:
        return self._data.get("feature_metadata", {})

    def has_atom_semantic(self) -> bool:
        return bool(self.topology.get("has_atom_semantic", False))

    def axis_size(self, axis: str) -> int:
        """Return the row count declared for ``axis``. Raises KeyError if absent."""
        return int(self.axis_sizes[axis])


@dataclass(frozen=True)
class TopologyGroup:
    """Topology sidecar projections from ``TopologySidecar::WriteFeatures``.

    Holds the four structured-NPY tables (residues / bonds / rings /
    ring_membership) and the parsed manifest JSON. Always present on a
    successful ``load()`` -- every field is populated by the loader.
    """
    residues: Optional[Residues] = None
    bonds: Optional[Bonds] = None
    rings: Optional[Rings] = None
    ring_membership: Optional[RingMembership] = None
    manifest: Optional[ExtractionManifest] = None


@dataclass(frozen=True)
class AIMNet2Group:
    """Static outputs from the caller-selected AIMNet2 TorchScript model.

    ``aim`` and ``aim_projection`` are opaque learned-channel arrays, not
    Cartesian/e3nn tensors. The E/EFG fields have their stated polar/rank-2
    laws only conditional on the loaded model emitting scalar charges under a
    rigid transform; the producer checks the model interface and finiteness,
    not arbitrary-``.jpt`` equivariance. See each :data:`CATALOG` entry for
    exact units, source partitions, basis conversion and absence rules.
    """
    charges: AIMNet2Charges
    aim: AIMNet2AimEmbedding
    efg: EFGTensor
    efg_aromatic: EFGTensor
    efg_backbone: EFGTensor
    efg_sidechain: EFGTensor
    E: VectorField
    E_backbone: VectorField
    E_sidechain: VectorField
    E_aromatic: VectorField
    energy_mlp: np.ndarray
    energy_shifted_local: np.ndarray
    energy_terms: np.ndarray
    d3_e_disp_atom: np.ndarray
    d3_cn: np.ndarray
    d3_c6_stats: np.ndarray
    aim_projection: np.ndarray
    # Current successful CLI extractions compute these in a separate Result
    # after the main AIMNet2 Result. They remain Optional in the reader for
    # pre-2026-05-09 files. Field renamed from `polarisability` to
    # `charge_response_gradient` 2026-05-20 (commit 58594f5) — the
    # emission is ∂(Σq²)/∂r, NOT a Buckingham α tensor.
    charge_response_gradient: Optional[AIMNet2ChargeResponseGradient] = None
    charge_response_gradient_scalar: Optional[np.ndarray] = None


@dataclass(frozen=True)
class PlanarGeometryGroup:
    """Planar and near-planar geometry outputs from PlanarGeometryResult.

    Axes are intentionally kept explicit because the result mixes atom,
    residue, aromatic-ring, and saturated-ring quantities:

    - ``pyramidalization``, ``pyramidalization_valid``, and
      ``pyramidalization_center_type`` are per atom, shape ``(N,)``.
    - ``omega_actual``, ``omega_deviation``, and ``omega_is_xpro`` are
      per residue, shape ``(R,)``.
    - ``aromatic_chi2`` is per aromatic ring.
    - ``pucker_Q`` and ``pucker_theta`` are per saturated ring.

    Convenience methods project the residue-level omega arrays back to
    atom rows using ``Protein.residue_index``. Ring-level arrays are not
    projected here because aromatic and saturated ring axes are distinct
    from the existing per-aromatic-ring ``RingGeometry`` table.
    """
    pyramidalization: Optional[np.ndarray] = None
    pyramidalization_valid: Optional[np.ndarray] = None
    pyramidalization_center_type: Optional[np.ndarray] = None
    omega_actual: Optional[np.ndarray] = None
    omega_deviation: Optional[np.ndarray] = None
    omega_is_xpro: Optional[np.ndarray] = None
    aromatic_chi2: Optional[np.ndarray] = None
    pucker_Q: Optional[np.ndarray] = None
    pucker_theta: Optional[np.ndarray] = None

    @staticmethod
    def _per_atom_residue_values(values: Optional[np.ndarray],
                                 residue_index: np.ndarray) -> np.ndarray:
        out = np.full(len(residue_index), np.nan, dtype=np.float64)
        if values is None:
            return out
        vals = np.asarray(values)
        ri = np.asarray(residue_index, dtype=np.intp)
        ok = (ri >= 0) & (ri < len(vals))
        out[ok] = vals[ri[ok]]
        return out

    def omega_actual_per_atom(self, residue_index: np.ndarray) -> np.ndarray:
        """Per-atom omega angle by each atom's parent residue."""
        return self._per_atom_residue_values(self.omega_actual, residue_index)

    def omega_deviation_per_atom(self, residue_index: np.ndarray) -> np.ndarray:
        """Per-atom wrapped omega deviation by each atom's parent residue."""
        return self._per_atom_residue_values(self.omega_deviation, residue_index)

    def omega_is_xpro_per_atom(self, residue_index: np.ndarray) -> np.ndarray:
        """Per-atom X-Pro mask by each atom's parent residue."""
        return self._per_atom_residue_values(self.omega_is_xpro, residue_index)


@dataclass(frozen=True)
class WaterFieldGroup:
    """Explicit water E-field and EFG from full-system trajectory."""
    efield: VectorField             # (N, 3) total water E-field
    efield_first: VectorField       # (N, 3) first-shell E-field
    efg: EFGTensor                  # (N, 5) total water EFG (T2 only)
    efg_first: EFGTensor            # (N, 5) first-shell EFG (T2 only)
    shell_counts: np.ndarray        # (N, 2) [n_first, n_second]
    efield_clamp_mask: Optional[np.ndarray] = None
    efield_clamp_scale: Optional[np.ndarray] = None
    efield_first_clamp_mask: Optional[np.ndarray] = None
    efield_first_clamp_scale: Optional[np.ndarray] = None


@dataclass(frozen=True)
class WaterHBondGeometryGroup:
    """Raw explicit-water H-bond candidates and per-atom summaries."""
    candidates: np.ndarray
    counts: np.ndarray
    nearest: np.ndarray


@dataclass(frozen=True)
class HydrationGroup:
    """Per-atom explicit-solvent hydration shell geometry.

    Columns are exposed-water fraction, mean dipole cosine, nearest-ion
    distance (Å), and charge (e). No shell gives ``0,0``; no ion inside the
    configured cutoff gives ``+inf,0``.
    """
    data: np.ndarray                # (N, 4) [asymmetry, dipole_cos, ion_dist, ion_charge]


@dataclass(frozen=True)
class WaterPolarizationGroup:
    """Water polarisation features using SASA-normal reference frame."""
    data: np.ndarray                # (N, 10) packed columns

    @property
    def dipole_vector(self) -> np.ndarray:
        """Net first-shell water dipole (N, 3)."""
        return self.data[:, 0:3]

    @property
    def surface_normal(self) -> np.ndarray:
        """SASA-derived outward surface normal (N, 3)."""
        return self.data[:, 3:6]

    @property
    def asymmetry(self) -> np.ndarray:
        """Half-shell asymmetry using SASA normal (N,)."""
        return self.data[:, 6]

    @property
    def dipole_alignment(self) -> np.ndarray:
        """cos(net dipole, SASA normal) (N,)."""
        return self.data[:, 7]

    @property
    def coherence(self) -> np.ndarray:
        """Legacy mean net dipole ``|sum d_i| / n`` in e_Angstrom (N,)."""
        return self.data[:, 8]

    @property
    def mean_net_dipole_eA(self) -> np.ndarray:
        """Mean net first-shell dipole ``|sum d_i| / n`` in e_Angstrom (N,)."""
        return self.data[:, 8]

    @property
    def shell_count(self) -> np.ndarray:
        """First-shell water count (N,)."""
        return self.data[:, 9]


@dataclass(frozen=True)
class EeqGroup:
    """Project-local QEq/EEQ-style geometry-dependent charges.

    project-local QEq/EEQ-style charge-equilibration model with error-function coordination number, CN-dependent electronegativity shift, Gaussian self term, and Ohno-Klopman off-diagonal kernel; parameters are in-repo/project-local and are not a validated dftd4/multicharge port.
    """
    charges: np.ndarray             # (N,) partial charges (elementary charges)
    cn: np.ndarray                  # (N,) coordination number
    chi_eff: Optional[np.ndarray] = None
    hardness: Optional[np.ndarray] = None
    coulomb: Optional[EeqCoulombGroup] = None


@dataclass(frozen=True)
class LarsenHBondPairs:
    """Split raw Larsen pair provenance plus the compatibility view.

    ``index`` columns are donor/acceptor atoms, residues, classes,
    disposition, six frame atoms, then five target masks. Dispositions are
    0=missing frame, 1=theta miss, 2=finite valid grid miss, 3=success,
    4=invalid/nonfinite frame, and 5=carboxylate symmetry-filtered sibling.
    ``geometry`` is
    ``[r, theta, rho, any_imputed, imputed_corner_count, frame_valid]``;
    ``frame_valid`` certifies both finite query geometry and a real,
    non-degenerate donor rotation frame (never an identity fallback).
    ``isotropic`` is the four Table-2 terms, diagnostic C-beta, and their
    scientific Table-2 total. ``compatibility`` concatenates 16+6+6 columns.
    """
    index: np.ndarray
    geometry: np.ndarray
    isotropic: np.ndarray
    compatibility: np.ndarray


@dataclass(frozen=True)
class LarsenHBondGroup:
    """Larsen 2015 ProCS15 H-bond term shielding contributions.

    Per Larsen 2015 Eq 5, each H-bond pair contributes four terms (1° +
    2° on each of HB / HαB) plus a water-term offset Δσ_w on amide H
    atoms with no H-bond partner. Tensors are ppm in the protein lab
    frame (rotated from the canonical donor frame at calculator time)
    and emitted as SphericalTensor-packed (T0+T1+T2 = 9 cols).  The signed-rho
    grid lookup is chiral and has a proper-rotation contract only, so these
    packed tensors are deliberately exposed as raw ndarrays rather than O(3)
    tensor wrappers.

    Methods accumulate side-by-side with the scalar-geometry ``HBondGroup``
    — both calculators cover overlapping physics (amide-H / backbone-O
    subset) but use different formulations (kernel × η vs grid lookup).
    Per-atom-type differences are themselves thesis-reportable. See
    ``feedback_methods_accumulate`` memory entry.

    The producer covers both backbone amide-H and Hα donors via spatial
    enumeration, with backbone carbonyl, sidechain carbonyl, hydroxyl, and
    carboxylate acceptor classes. SidechainCarbonyl is explicitly surfaced as
    the frozen BackboneCarbonyl/NMA archive approximation.

    Per-atom convention:

    - ``shielding`` is the sum over all four contribution classes that
      apply at this atom per Larsen 2015 Table 2 dispatch (encoded as
      ``LarsenContribDispatch::Applies`` in the C++ side). Cβ does NOT
      contribute (Table 2 says so) — its zero column here is the
      physics statement, not an absence.
    - The per-class columns hold each contribution separately for
      downstream ML stratification.
    - ``diagnostic_CB`` should be near-zero in production (Larsen Table
      2 says Cβ gets no HB term; non-zero would signal a pipeline bug
      in the parser/loader/rotation path).
    - ``water_term`` is 2.07 ppm isotropic on amide Hs with zero geometric
      H-bond candidates under the Larsen 4.2 Å / theta >= 90 degree gate.
    - ``count`` counts H-bond pairs that contributed to this atom under
      any of the four Table 2 classes — the diagnostic CB does NOT
      increment it.
    """
    shielding: Optional[np.ndarray] = None
    pHB_1: Optional[np.ndarray] = None
    pHB_2: Optional[np.ndarray] = None
    pHaB_1: Optional[np.ndarray] = None
    pHaB_2: Optional[np.ndarray] = None
    diagnostic_CB: Optional[np.ndarray] = None
    water_term: Optional[np.ndarray] = None
    count: Optional[np.ndarray] = None
    corner_imputed: Optional[np.ndarray] = None
    imputed_pair_count: Optional[np.ndarray] = None
    sidechain_carbonyl_pair_count: Optional[np.ndarray] = None
    pairs: Optional[LarsenHBondPairs] = None

    @property
    def pairs_index(self) -> Optional[np.ndarray]:
        return None if self.pairs is None else self.pairs.index

    @property
    def pairs_geometry(self) -> Optional[np.ndarray]:
        return None if self.pairs is None else self.pairs.geometry

    @property
    def pairs_isotropic(self) -> Optional[np.ndarray]:
        return None if self.pairs is None else self.pairs.isotropic

    @property
    def pairs_compatibility(self) -> Optional[np.ndarray]:
        return None if self.pairs is None else self.pairs.compatibility


@dataclass(frozen=True)
class LarsenSidechainDonorAuditGroup:
    """Audit-only sidechain polar-H donors omitted by Larsen Table 2."""
    atoms: np.ndarray
    candidates: np.ndarray


# ── Top-level protein container ─────────────────────────────────────


@dataclass
class Protein:
    """All extracted features for one protein conformation.

    Every field is typed. Optional groups are None when the extractor
    did not produce them (MOPAC, APBS, Orca, Delta).
    """
    protein_id: str
    n_atoms: int

    # Identity
    pos: PositionField
    element: np.ndarray
    residue_type: np.ndarray
    residue_index: np.ndarray

    # Ring calculators
    biot_savart: BiotSavartGroup
    haigh_mallion: HaighMallionGroup
    pi_quadrupole: PiQuadrupoleGroup
    ring_susceptibility: Optional[RingSusceptibilityGroup]
    dispersion: DispersionGroup

    # Per-ring sparse data
    ring_contributions: RingContributions = None
    ring_direction_to_center: Optional[VectorField] = None
    ring_geometry: RingGeometry = None
    ring_pair_geometry: RingPairGeometry = None

    # Foundation calculators
    enrichment: Optional[EnrichmentGroup] = None
    charge_assignment: Optional[ChargeAssignmentGroup] = None
    spatial_index: Optional[SpatialIndexGroup] = None
    molecular_graph: Optional[MolecularGraphGroup] = None

    # Bond calculators
    mcconnell: McConnellGroup = None
    sidechain_carbonyl_anisotropy: Optional[
        SidechainCarbonylAnisotropyGroup] = None
    coulomb: CoulombGroup = None
    hbond: HBondGroup = None
    dssp: DsspScalars = None
    dssp_observed: np.ndarray = None
    dssp_ss8: np.ndarray = None
    dssp_ppii: np.ndarray = None
    dssp_hbond_energy: np.ndarray = None
    dssp_chi: np.ndarray = None
    sasa: np.ndarray = None
    sasa_normal: Optional[VectorField] = None  # (N, 3) outward surface normal
    sasa_fraction: Optional[np.ndarray] = None

    # Optional calculator groups
    mopac: Optional[MopacGroup] = None
    apbs: Optional[APBSGroup] = None
    orca: Optional[OrcaGroup] = None
    delta: Optional[DeltaGroup] = None
    aimnet2: Optional[AIMNet2Group] = None
    planar_geometry: Optional[PlanarGeometryGroup] = None

    # Per-atom invariant categorical record (CategoryInfoProjection).
    category_info: Optional[CategoryInfo] = None

    # Topology sidecar projections (TopologySidecar, 2026-05-13).
    # Always present on a successful load -- bonds/rings/ring_membership
    # NPYs are required and ``extraction_manifest.json`` is required.
    topology: TopologyGroup = None

    # Explicit solvent (trajectory path only)
    water_field: Optional[WaterFieldGroup] = None
    water_hbond_geometry: Optional[WaterHBondGeometryGroup] = None
    hydration: Optional[HydrationGroup] = None
    water_polarization: Optional[WaterPolarizationGroup] = None
    gromacs_energy: Optional[np.ndarray] = None
    # (N,7) equal-share local BondedParameters evaluation; not an EDR
    # per-atom decomposition.
    bonded_energy: Optional[np.ndarray] = None

    # Optional legacy NPY projections from TrajectoryProtein::WriteFeatures.
    # Rich trajectory consumers should prefer trajectory.h5.
    hm_welford: Optional[np.ndarray] = None       # (N,11), scalar rollup
    mc_welford: Optional[np.ndarray] = None       # (N,11), scalar rollup
    sasa_welford: Optional[np.ndarray] = None     # (N,7), finite-grid rollup

    # Geometry-dependent charges
    eeq: Optional[EeqGroup] = None

    larsen_hbond: Optional[LarsenHBondGroup] = None
    larsen_sidechain_donor_audit: Optional[
        LarsenSidechainDonorAuditGroup] = None

    @property
    def relative_sasa(self) -> Optional[np.ndarray]:
        """Alias for the producer-emitted normalized SASA fraction."""
        return self.sasa_fraction


# ── Loader ──────────────────────────────────────────────────────────


def _wrap(spec, data: np.ndarray):
    """Wrap raw array in its typed class, or return as-is for np.ndarray."""
    if spec.wrapper is np.ndarray:
        return data
    return spec.wrapper(data)


def _validate_topology_invariants(tg: TopologyGroup, n_atoms: int,
                                    protein_id: str,
                                    atom_residue_index: np.ndarray | None = None
                                    ) -> None:
    """Cross-check declared axis sizes + reference integrity.

    Codex contract first-pass validation gates: the manifest declares
    axis sizes; the structured-NPY tables must agree, and every
    cross-reference (bond endpoint -> atom, ring membership -> atom +
    ring) must be in-bounds. Raises ``ValueError`` on any violation.

    This is the SDK's role in "force python to use the model": a
    malformed export fails loud at load() before the consumer sees any
    derived analysis.
    """
    man = tg.manifest
    if man is None:
        raise ValueError(
            f"{protein_id}: topology sidecar has no manifest; "
            "cannot validate invariants")

    def check_axis(axis: str, actual: int) -> None:
        declared = int(man.axis_sizes.get(axis, -1))
        if declared != actual:
            raise ValueError(
                f"{protein_id}: manifest declares {axis} axis size "
                f"{declared} but on-disk row count is {actual}")

    check_axis("atom", n_atoms)
    check_axis("residue", tg.residues.n_residues)
    check_axis("bond", tg.bonds.n_bonds)
    check_axis("aromatic_ring",
                 int((tg.rings.ring_kind == 0).sum()))
    check_axis("saturated_ring",
                 int((tg.rings.ring_kind == 1).sum()))
    check_axis("ring", tg.rings.n_rings)
    check_axis("ring_membership", tg.ring_membership.n_rows)

    # Residue's atom_count sums to total atom count (each atom belongs
    # to exactly one residue in our model).
    if tg.residues.n_residues > 0:
        if (tg.residues.atom_count < 0).any():
            raise ValueError(
                f"{protein_id}: residues.npy atom_count contains negative values")
        total_atoms_from_residues = int(tg.residues.atom_count.sum())
        if total_atoms_from_residues != n_atoms:
            raise ValueError(
                f"{protein_id}: residues.npy atom_count sums to "
                f"{total_atoms_from_residues}; atom axis size is {n_atoms}")

    if atom_residue_index is not None:
        ari = np.asarray(atom_residue_index)
        if len(ari) != n_atoms:
            raise ValueError(
                f"{protein_id}: residue_index.npy has {len(ari)} rows; "
                f"atom axis size is {n_atoms}")
        if tg.residues.n_residues == 0 and n_atoms != 0:
            raise ValueError(
                f"{protein_id}: residue axis is empty but atom axis has "
                f"{n_atoms} rows")
        if tg.residues.n_residues > 0:
            bad_atom_residue = (
                (ari < 0) | (ari >= tg.residues.n_residues)
            )
            if bad_atom_residue.any():
                raise ValueError(
                    f"{protein_id}: residue_index.npy references outside "
                    f"the residue axis [0, {tg.residues.n_residues})")
            atom_counts = np.bincount(
                ari.astype(np.intp), minlength=tg.residues.n_residues)
            if not np.array_equal(atom_counts, tg.residues.atom_count):
                raise ValueError(
                    f"{protein_id}: residue_index.npy counts do not match "
                    "residues.npy atom_count")

    def check_row_identity(label: str, values: np.ndarray) -> None:
        expected = np.arange(len(values), dtype=values.dtype)
        if not np.array_equal(values, expected):
            raise ValueError(
                f"{protein_id}: {label} must equal its row index for every row")

    check_row_identity("residues.npy residue_index", tg.residues.residue_index)
    check_row_identity("bonds.npy bond_index", tg.bonds.bond_index)
    check_row_identity("rings.npy ring_id", tg.rings.ring_id)

    # Atom-axis residue_index.npy references the residue axis.
    ri = np.asarray(tg.residues.residue_index)
    if tg.residues.n_residues > 0:
        bad_prev = (
            (tg.residues.prev_residue_index != -1) &
            ((tg.residues.prev_residue_index < 0) |
             (tg.residues.prev_residue_index >= tg.residues.n_residues))
        )
        bad_next = (
            (tg.residues.next_residue_index != -1) &
            ((tg.residues.next_residue_index < 0) |
             (tg.residues.next_residue_index >= tg.residues.n_residues))
        )
        if bad_prev.any() or bad_next.any():
            raise ValueError(
                f"{protein_id}: residues.npy prev/next residue references "
                f"outside the residue axis [0, {tg.residues.n_residues})")
        if np.any(ri < 0):
            raise ValueError(
                f"{protein_id}: residues.npy residue_index contains negative rows")

    # Bond endpoints reference the atom axis.
    bonds = tg.bonds
    if bonds.n_bonds > 0:
        bad_a = ((bonds.atom_index_a < 0) | (bonds.atom_index_a >= n_atoms))
        bad_b = ((bonds.atom_index_b < 0) | (bonds.atom_index_b >= n_atoms))
        if bad_a.any() or bad_b.any():
            raise ValueError(
                f"{protein_id}: bonds.npy carries endpoints outside the "
                f"atom axis [0, {n_atoms})")

    # Ring membership references the atom + ring axes.
    rm = tg.ring_membership
    if rm.n_rows > 0:
        bad_atom = ((rm.atom_index < 0) | (rm.atom_index >= n_atoms))
        bad_ring = ((rm.ring_id < 0) | (rm.ring_id >= tg.rings.n_rings))
        if bad_atom.any():
            raise ValueError(
                f"{protein_id}: ring_membership.npy carries atom_index outside "
                f"the atom axis [0, {n_atoms})")
        if bad_ring.any():
            raise ValueError(
                f"{protein_id}: ring_membership.npy carries ring_id outside "
                f"the ring axis [0, {tg.rings.n_rings})")

    # Fused-ring partner index must reference a valid ring (-1 = none).
    rings = tg.rings
    if rings.n_rings > 0:
        fp = rings.fused_partner_ring_id
        bad_fp = ((fp != -1) & ((fp < 0) | (fp >= rings.n_rings)))
        if bad_fp.any():
            raise ValueError(
                f"{protein_id}: rings.npy has fused_partner_ring_id outside "
                f"[-1, {rings.n_rings})")

        parent = rings.parent_residue_index
        bad_parent = (
            (parent < 0) | (parent >= tg.residues.n_residues)
        )
        if bad_parent.any():
            raise ValueError(
                f"{protein_id}: rings.npy parent_residue_index outside the "
                f"residue axis [0, {tg.residues.n_residues})")

        for kind, axis_name in ((0, "aromatic_ring"), (1, "saturated_ring")):
            mask = rings.ring_kind == kind
            native = np.asarray(rings.native_axis_index[mask])
            expected = np.arange(len(native), dtype=native.dtype)
            if not np.array_equal(native, expected):
                raise ValueError(
                    f"{protein_id}: rings.npy native_axis_index for "
                    f"{axis_name} must be contiguous 0..N-1")

        valid_kinds = (rings.ring_kind == 0) | (rings.ring_kind == 1)
        if not valid_kinds.all():
            raise ValueError(
                f"{protein_id}: rings.npy ring_kind must be 0 aromatic or "
                "1 saturated")

        membership_counts = np.bincount(
            rm.ring_id.astype(np.intp), minlength=rings.n_rings)
        if not np.array_equal(membership_counts, rings.atom_count):
            raise ValueError(
                f"{protein_id}: rings.npy atom_count does not match "
                "ring_membership.npy rows per ring")


def load(path: str | Path) -> Protein:
    """Load an extraction directory into a fully typed Protein.

    Validates every file against the catalog. Warns on unregistered files
    (forward-compatible). Errors on missing required files.
    """
    path = Path(path)
    protein_id = path.name

    # Load all NPY files
    available: dict[str, np.ndarray] = {}
    for npy in path.glob("*.npy"):
        available[npy.stem] = np.load(npy)

    # Warn on unregistered (don't error — forward-compatible)
    unregistered = set(available.keys()) - set(CATALOG.keys())
    if unregistered:
        warnings.warn(
            f"Unregistered NPY files in {protein_id}: {sorted(unregistered)}. "
            f"Update the catalog in nmr_extract._catalog to wrap them.",
            stacklevel=2,
        )

    # Check required files
    for stem, spec in CATALOG.items():
        if spec.required and stem not in available:
            raise FileNotFoundError(
                f"Required file {stem}.npy missing for {protein_id}")

    # Wrap each array
    wrapped: dict[str, object] = {}
    for stem, data in available.items():
        if stem in CATALOG:
            wrapped[stem] = _wrap(CATALOG[stem], data)

    def get(stem: str, default=None):
        return wrapped.get(stem, default)

    n_atoms = available["pos"].shape[0]

    # ── Assemble groups ─────────────────────────────────────────

    biot_savart = BiotSavartGroup(
        shielding=get("bs_shielding"),
        per_type_T0=get("bs_per_type_T0"),
        per_type_T2=get("bs_per_type_T2"),
        per_type_T1=get("bs_per_type_T1"),
        total_B=get("bs_total_B"),
        ring_B_field=get("bs_ring_B_field"),
        ring_B_cylindrical=get("bs_ring_B_cylindrical"),
        ring_counts=get("bs_ring_counts"),
    )
    haigh_mallion = HaighMallionGroup(
        shielding=get("hm_shielding"),
        per_type_T0=get("hm_per_type_T0"),
        per_type_T2=get("hm_per_type_T2"),
        per_type_T1=get("hm_per_type_T1"),
        ring_B_field=get("hm_ring_B_field"),
    )
    pi_quadrupole = PiQuadrupoleGroup(
        per_type_T0=get("pq_per_type_T0"),
        quad_scalar=get("piquad_quad_scalar"),
        axial_scalar_per_type_T0=get("piquad_axial_scalar_per_type_T0"),
        local_tensor=get("piquad_local_tensor"),
        local_T2=get("piquad_local_T2"),
        local_frame=get("piquad_local_frame"),
        local_geometry=get("piquad_local_geometry"),
    )
    ring_susceptibility = RingSusceptibilityGroup(
        scalar=get("ringchi_scalar"),
        per_type_T0=get("ringchi_per_type_T0"),
    )
    dispersion = DispersionGroup(
        per_type_T0=get("disp_per_type_T0"),
        aromatic_r6_proximity_per_type_T0=(
            get("aromatic_r6_proximity_per_type_T0")
            if "aromatic_r6_proximity_per_type_T0" in available else None
        ),
    )
    mcconnell = McConnellGroup(
        peptide_co_fixed=get("mc_peptide_co_fixed"),
        peptide_co_bo=get("mc_peptide_co_bo"),
        peptide_co_rhombic=get("mc_peptide_co_rhombic"),
        peptide_cn_fixed=get("mc_peptide_cn_fixed"),
        peptide_cn_bo=get("mc_peptide_cn_bo"),
        backbone_other_fixed=get("mc_backbone_other_fixed"),
        backbone_other_bo=get("mc_backbone_other_bo"),
        sidechain_co_fixed=get("mc_sidechain_co_fixed"),
        sidechain_co_bo=get("mc_sidechain_co_bo"),
        sidechain_other_fixed=get("mc_sidechain_other_fixed"),
        sidechain_other_bo=get("mc_sidechain_other_bo"),
        disulfide_fixed=get("mc_disulfide_fixed"),
        disulfide_bo=get("mc_disulfide_bo"),
        aromatic_zeroed_fixed=get("mc_aromatic_zeroed_fixed"),
        aromatic_zeroed_bo=get("mc_aromatic_zeroed_bo"),
        backbone_xh_fixed=get("mc_backbone_xh_fixed"),
        backbone_xh_bo=get("mc_backbone_xh_bo"),
        sidechain_xh_fixed=get("mc_sidechain_xh_fixed"),
        sidechain_xh_bo=get("mc_sidechain_xh_bo"),
        s_h_fixed=get("mc_s_h_fixed"),
        s_h_bo=get("mc_s_h_bo"),
        nearfield_counts=get("mc_nearfield_counts"),
        nearest_co_dir=get("mc_nearest_co_dir"),
        nearest_co_midpoint=get("mc_nearest_co_midpoint"),
        nearest_co_T2=get("mc_nearest_co_T2"),
        nearest_cn_T2=get("mc_nearest_cn_T2"),
        bond_neighbors=get("mc_bond_neighbors")
            if "mc_bond_neighbors" in available else None,
    )
    sidechain_carbonyl_anisotropy = None
    if "sidechain_co_source_bonds" in available:
        sidechain_carbonyl_anisotropy = SidechainCarbonylAnisotropyGroup(
            source_bonds=get("sidechain_co_source_bonds"),
            frame=get("sidechain_co_frame"),
            frame_quality=get("sidechain_co_frame_quality"),
            fixed_T2=get("sidechain_co_fixed_T2"),
            bo_T2=get("sidechain_co_bo_T2"),
            scalar_audit=get("sidechain_co_scalar_audit"),
        )
    coulomb = CoulombGroup(
        efg=get("coulomb_efg"),
        E=get("coulomb_E"),
        E_backbone=get("coulomb_E_backbone"),
        E_sidechain=get("coulomb_E_sidechain"),
        E_aromatic=get("coulomb_E_aromatic"),
        efg_backbone=get("coulomb_efg_backbone"),
        efg_sidechain=get("coulomb_efg_sidechain"),
        efg_aromatic=get("coulomb_efg_aromatic"),
        scalars=get("coulomb_scalars"),
        efg_t2=get("coulomb_efg_t2"),
        aromatic_E_proj=get("coulomb_aromatic_E_proj"),
        aromatic_n_src=get("coulomb_aromatic_n_src"),
        E_solvent=get("coulomb_E_solvent")
            if "coulomb_E_solvent" in available else None,
        efg_solvent=get("coulomb_efg_solvent")
            if "coulomb_efg_solvent" in available else None,
    )
    hbond = HBondGroup(
        scalars=get("hbond_scalars"),
        nearest_dir=get("hbond_nearest_dir"),
        flags=get("hbond_flags"),
    )

    enrichment = None
    if "enrichment_role" in available:
        enrichment = EnrichmentGroup(
            role=get("enrichment_role"),
            hybridisation=get("enrichment_hybridisation"),
            flags=get("enrichment_flags"),
            parent_is_sp2=get("enrichment_parent_is_sp2")
                if "enrichment_parent_is_sp2" in available else None,
            polar_h_kind=get("semantic_polar_h_kind")
                if "semantic_polar_h_kind" in available else None,
            planar_group_kind=get("semantic_planar_group_kind")
                if "semantic_planar_group_kind" in available else None,
            formal_charge=get("semantic_formal_charge")
                if "semantic_formal_charge" in available else None,
            ring_position=get("semantic_ring_position")
                if "semantic_ring_position" in available else None,
            locant=get("semantic_locant")
                if "semantic_locant" in available else None,
            donor_class=get("enrichment_donor_class")
                if "enrichment_donor_class" in available else None,
            acceptor_class=get("enrichment_acceptor_class")
                if "enrichment_acceptor_class" in available else None,
            hybridisation_class=get("enrichment_hybridisation_class")
                if "enrichment_hybridisation_class" in available else None,
        )

    charge_assignment = None
    if "ff_partial_charge" in available:
        charge_assignment = ChargeAssignmentGroup(
            partial_charge=get("ff_partial_charge"),
            pb_radius=get("ff_pb_radius"),
        )

    spatial_index = None
    if "spatial_neighbors" in available:
        spatial_index = SpatialIndexGroup(neighbors=get("spatial_neighbors"))

    molecular_graph = None
    if "molecular_graph_int" in available:
        molecular_graph = MolecularGraphGroup(
            integer=get("molecular_graph_int"),
            floating=get("molecular_graph_float"),
            compatibility=get("molecular_graph"),
        )

    # MOPAC (optional)
    mopac = None
    if "mopac_charges" in available:
        mopac_full_stems = {
            "mopac_charges_full_precision",
            "mopac_bond_orders_full_precision",
            "mopac_bond_valencies_full_precision",
            "mopac_atom_populations",
            "mopac_atomic_orbital_populations",
            "mopac_atomic_orbital_population_totals",
            "mopac_bond_valencies",
            "mopac_bond_orders_unique",
            "mopac_topology_bond_orders_full",
            "mopac_heat_kcal_mol",
            "mopac_dipole_debye",
            "mopac_dipole_point_charge_debye",
            "mopac_dipole_hybridization_debye",
            "mopac_bond_index",
            "mopac_bond_atom",
            "mopac_bond_order",
            "mopac_ao_max_orbitals",
            "mopac_ao_orbitals_per_atom",
            "mopac_atom_ao_density",
            "mopac_atomic_orbital_populations_full_precision",
            "mopac_bond_ao_density_directed",
            "mopac_bond_density_pairs",
            "mopac_bond_ao_density",
            "mopac_atom_electron_population",
            "mopac_atom_s_population",
            "mopac_atom_p_population",
            "mopac_atom_d_population",
            "mopac_lewis_bond_count",
            "mopac_lewis_bond_atoms",
            "mopac_lmo_energy_levels",
            "mopac_lmo_occupied_atom_counts",
            "mopac_lmo_occupied_atoms",
            "mopac_lmo_occupied_coefficients",
            "mopac_lmo_virtual_atom_counts",
            "mopac_lmo_virtual_atoms",
            "mopac_lmo_virtual_coefficients",
            "mopac_lmo_occupied_atom_offsets_native",
            "mopac_lmo_virtual_atom_offsets_native",
            "mopac_lmo_occupied_coefficient_offsets_native",
            "mopac_lmo_virtual_coefficient_offsets_native",
            "mopac_lmo_occupied_atom_storage_native",
            "mopac_lmo_virtual_atom_storage_native",
            "mopac_lmo_occupied_coefficient_storage_native",
            "mopac_lmo_virtual_coefficient_storage_native",
            "mopac_mozyme_state_dimensions",
        }
        mopac_full = None
        if any(stem in available for stem in mopac_full_stems):
            def mopac_optional(stem: str):
                return get(stem) if stem in available else None

            mopac_full = MopacFullGroup(
                charges_full_precision=mopac_optional("mopac_charges_full_precision"),
                bond_orders_full_precision=mopac_optional("mopac_bond_orders_full_precision"),
                bond_valencies_full_precision=mopac_optional("mopac_bond_valencies_full_precision"),
                atom_populations=mopac_optional("mopac_atom_populations"),
                atomic_orbital_populations=mopac_optional("mopac_atomic_orbital_populations"),
                atomic_orbital_population_totals=mopac_optional("mopac_atomic_orbital_population_totals"),
                bond_valencies=mopac_optional("mopac_bond_valencies"),
                bond_orders_unique=mopac_optional("mopac_bond_orders_unique"),
                topology_bond_orders_full=mopac_optional("mopac_topology_bond_orders_full"),
                heat_kcal_mol=mopac_optional("mopac_heat_kcal_mol"),
                dipole_debye=mopac_optional("mopac_dipole_debye"),
                dipole_point_charge_debye=mopac_optional("mopac_dipole_point_charge_debye"),
                dipole_hybridization_debye=mopac_optional("mopac_dipole_hybridization_debye"),
                bond_index=mopac_optional("mopac_bond_index"),
                bond_atom=mopac_optional("mopac_bond_atom"),
                bond_order=mopac_optional("mopac_bond_order"),
                ao_max_orbitals=mopac_optional("mopac_ao_max_orbitals"),
                ao_orbitals_per_atom=mopac_optional("mopac_ao_orbitals_per_atom"),
                atom_ao_density=mopac_optional("mopac_atom_ao_density"),
                atomic_orbital_populations_full_precision=mopac_optional(
                    "mopac_atomic_orbital_populations_full_precision"),
                bond_ao_density_directed=mopac_optional("mopac_bond_ao_density_directed"),
                bond_density_pairs=mopac_optional("mopac_bond_density_pairs"),
                bond_ao_density=mopac_optional("mopac_bond_ao_density"),
                atom_electron_population=mopac_optional("mopac_atom_electron_population"),
                atom_s_population=mopac_optional("mopac_atom_s_population"),
                atom_p_population=mopac_optional("mopac_atom_p_population"),
                atom_d_population=mopac_optional("mopac_atom_d_population"),
                lewis_bond_count=mopac_optional("mopac_lewis_bond_count"),
                lewis_bond_atoms=mopac_optional("mopac_lewis_bond_atoms"),
                lmo_energy_levels=mopac_optional("mopac_lmo_energy_levels"),
                lmo_occupied_atom_counts=mopac_optional("mopac_lmo_occupied_atom_counts"),
                lmo_occupied_atoms=mopac_optional("mopac_lmo_occupied_atoms"),
                lmo_occupied_coefficients=mopac_optional("mopac_lmo_occupied_coefficients"),
                lmo_virtual_atom_counts=mopac_optional("mopac_lmo_virtual_atom_counts"),
                lmo_virtual_atoms=mopac_optional("mopac_lmo_virtual_atoms"),
                lmo_virtual_coefficients=mopac_optional("mopac_lmo_virtual_coefficients"),
                lmo_occupied_atom_offsets_native=mopac_optional("mopac_lmo_occupied_atom_offsets_native"),
                lmo_virtual_atom_offsets_native=mopac_optional("mopac_lmo_virtual_atom_offsets_native"),
                lmo_occupied_coefficient_offsets_native=mopac_optional("mopac_lmo_occupied_coefficient_offsets_native"),
                lmo_virtual_coefficient_offsets_native=mopac_optional("mopac_lmo_virtual_coefficient_offsets_native"),
                lmo_occupied_atom_storage_native=mopac_optional("mopac_lmo_occupied_atom_storage_native"),
                lmo_virtual_atom_storage_native=mopac_optional("mopac_lmo_virtual_atom_storage_native"),
                lmo_occupied_coefficient_storage_native=mopac_optional("mopac_lmo_occupied_coefficient_storage_native"),
                lmo_virtual_coefficient_storage_native=mopac_optional("mopac_lmo_virtual_coefficient_storage_native"),
                mozyme_state_dimensions=mopac_optional("mopac_mozyme_state_dimensions"),
            )
        mopac = MopacGroup(
            core=MopacCoreGroup(
                charges=get("mopac_charges"),
                scalars=get("mopac_scalars"),
                bond_orders=get("mopac_bond_orders"),
                global_=get("mopac_global"),
                bond_neighbors=get("mopac_bond_neighbors")
                    if "mopac_bond_neighbors" in available else None,
            ),
            coulomb=MopacCoulombGroup(
                efg=get("mopac_coulomb_efg"),
                E=get("mopac_coulomb_E"),
                E_backbone=get("mopac_coulomb_E_backbone"),
                E_sidechain=get("mopac_coulomb_E_sidechain"),
                E_aromatic=get("mopac_coulomb_E_aromatic"),
                efg_backbone=get("mopac_coulomb_efg_backbone"),
                efg_sidechain=get("mopac_coulomb_efg_sidechain"),
                efg_aromatic=get("mopac_coulomb_efg_aromatic"),
                scalars=get("mopac_coulomb_scalars"),
                aromatic_E_proj=get("mopac_coulomb_aromatic_E_proj"),
                aromatic_n_src=get("mopac_coulomb_aromatic_n_src"),
                E_clamp_mask=get("mopac_coulomb_E_clamp_mask"),
                E_clamp_scale=get("mopac_coulomb_E_clamp_scale"),
            ),
            full=mopac_full,
            mcconnell=MopacMcConnellGroup(
                co_sum=get("mopac_mc_co_sum"),
                cn_sum=get("mopac_mc_cn_sum"),
                sidechain_sum=get("mopac_mc_sidechain_sum"),
                aromatic_sum=get("mopac_mc_aromatic_sum"),
                co_nearest=get("mopac_mc_co_nearest"),
                nearest_co_dist=get("mopac_mc_nearest_co_dist"),
                nearest_cn_dist=get("mopac_mc_nearest_cn_dist"),
                nearest_co_T2=get("mopac_mc_nearest_co_T2"),
                nearest_cn_T2=get("mopac_mc_nearest_cn_T2"),
                backbone_total=get("mopac_mc_backbone_total"),
                sidechain_total=get("mopac_mc_sidechain_total"),
                aromatic_total=get("mopac_mc_aromatic_total"),
                shielding=get("mopac_mc_shielding"),
            ) if "mopac_mc_shielding" in available else None,
        )

    # APBS (optional)
    apbs = None
    if "apbs_E" in available:
        apbs = APBSGroup(
            E=get("apbs_E"),
            efg=get("apbs_efg"),
            phi=get("apbs_phi") if "apbs_phi" in available else None,
            E_clamp_mask=get("apbs_E_clamp_mask")
                if "apbs_E_clamp_mask" in available else None,
            E_clamp_scale=get("apbs_E_clamp_scale")
                if "apbs_E_clamp_scale" in available else None,
            nonfinite_sanitizer_mask=get("apbs_nonfinite_sanitizer_mask")
                if "apbs_nonfinite_sanitizer_mask" in available else None,
            E_total_diagnostic=get("apbs_E_total_diagnostic")
                if "apbs_E_total_diagnostic" in available else None,
            efg_total_diagnostic=get("apbs_efg_total_diagnostic")
                if "apbs_efg_total_diagnostic" in available else None,
        )

    # Orca DFT (optional)
    orca = None
    if "orca_total" in available:
        orca = OrcaGroup(
            total=get("orca_total"),
            diamagnetic=get("orca_diamagnetic"),
            paramagnetic=get("orca_paramagnetic"),
        )

    # Delta (optional)
    delta = None
    if "delta_shielding" in available:
        delta = DeltaGroup(
            shielding=get("delta_shielding"),
            scalars=get("delta_scalars"),
            apbs=get("delta_apbs"),
            ring_proximity=get("delta_ring_proximity"),
            graph=get("delta_graph"),
            wt_shielding_diamagnetic=get("wt_shielding_diamagnetic"),
            wt_shielding_paramagnetic=get("wt_shielding_paramagnetic"),
            mut_shielding_diamagnetic=get("mut_shielding_diamagnetic"),
            mut_shielding_paramagnetic=get("mut_shielding_paramagnetic"),
            delta_shielding_diamagnetic=get("delta_shielding_diamagnetic"),
            delta_shielding_paramagnetic=get("delta_shielding_paramagnetic"),
        )

    # Per-atom invariant categorical record (optional; CategoryInfoProjection
    # emits it whenever LegacyAmber substrate is populated, which is every
    # production load path). The catalog declares wrapper=np.ndarray to
    # avoid a circular import; we wrap to CategoryInfo here.
    category_info = None
    if "atoms_category_info" in available:
        category_info = CategoryInfo(available["atoms_category_info"])

    # Topology sidecar (TopologySidecar). bonds.npy / rings.npy /
    # ring_membership.npy are required NPYs (declared in CATALOG);
    # the missing-required check above already failed if any were
    # absent. extraction_manifest.json is required separately.
    manifest_path = path / "extraction_manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(
            f"Required topology sidecar extraction_manifest.json missing for {protein_id}")
    import json
    with open(manifest_path) as f:
        manifest_obj = ExtractionManifest(json.load(f))
    topology_group = TopologyGroup(
        residues=Residues(available["residues"]),
        bonds=Bonds(available["bonds"]),
        rings=Rings(available["rings"]),
        ring_membership=RingMembership(available["ring_membership"]),
        manifest=manifest_obj,
    )
    _validate_topology_invariants(
        topology_group, n_atoms, protein_id, available.get("residue_index"))
    pair_geometry_raw = np.asarray(available["ring_pair_geometry"])
    n_aromatic = int((topology_group.rings.ring_kind == 0).sum())
    expected_pairs = n_aromatic * (n_aromatic - 1) // 2
    if pair_geometry_raw.shape != (expected_pairs, 13):
        raise ValueError(
            f"{protein_id}: ring_pair_geometry.npy has shape "
            f"{pair_geometry_raw.shape}; expected ({expected_pairs}, 13)")
    if expected_pairs:
        expected_ab = np.array(
            [(a, b) for a in range(n_aromatic)
                    for b in range(a + 1, n_aromatic)], dtype=np.intp)
        actual_ab = pair_geometry_raw[:, :2].astype(np.intp)
        if not np.array_equal(actual_ab, expected_ab):
            raise ValueError(
                f"{protein_id}: ring_pair_geometry.npy rows must be "
                "lexicographic i<j over the aromatic-ring axis")

    # AIMNet2 production arrays are required by the catalog; this guard keeps
    # legacy/pre-contract extractions from constructing a partial group.
    aimnet2 = None
    if "aimnet2_charges" in available:
        aimnet2 = AIMNet2Group(
            charges=get("aimnet2_charges"),
            aim=get("aimnet2_aim"),
            efg=get("aimnet2_efg"),
            efg_aromatic=get("aimnet2_efg_aromatic"),
            efg_backbone=get("aimnet2_efg_backbone"),
            efg_sidechain=get("aimnet2_efg_sidechain"),
            E=get("aimnet2_E"),
            E_backbone=get("aimnet2_E_backbone"),
            E_sidechain=get("aimnet2_E_sidechain"),
            E_aromatic=get("aimnet2_E_aromatic"),
            energy_mlp=get("aimnet2_energy_mlp"),
            energy_shifted_local=get("aimnet2_energy_shifted_local"),
            energy_terms=get("aimnet2_energy_terms"),
            d3_e_disp_atom=get("aimnet2_d3_e_disp_atom"),
            d3_cn=get("aimnet2_d3_cn"),
            d3_c6_stats=get("aimnet2_d3_c6_stats"),
            aim_projection=get("aimnet2_aim_projection"),
            charge_response_gradient=get("aimnet2_charge_response_gradient")
                if "aimnet2_charge_response_gradient" in available else None,
            charge_response_gradient_scalar=get("aimnet2_charge_response_gradient_scalar")
                if "aimnet2_charge_response_gradient_scalar" in available else None,
        )

    # Planar geometry (optional in pre-2026-05-09 outputs and in unusual
    # fixtures where the topology substrate is not populated).
    planar_geometry = None
    planar_stems = {
        "pyramidalization",
        "pyramidalization_valid",
        "pyramidalization_center_type",
        "omega_actual",
        "omega_deviation",
        "omega_is_xpro",
        "aromatic_chi2",
        "pucker_Q",
        "pucker_theta",
    }
    if any(stem in available for stem in planar_stems):
        planar_geometry = PlanarGeometryGroup(
            pyramidalization=get("pyramidalization")
                if "pyramidalization" in available else None,
            pyramidalization_valid=get("pyramidalization_valid")
                if "pyramidalization_valid" in available else None,
            pyramidalization_center_type=get("pyramidalization_center_type")
                if "pyramidalization_center_type" in available else None,
            omega_actual=get("omega_actual")
                if "omega_actual" in available else None,
            omega_deviation=get("omega_deviation")
                if "omega_deviation" in available else None,
            omega_is_xpro=get("omega_is_xpro")
                if "omega_is_xpro" in available else None,
            aromatic_chi2=get("aromatic_chi2")
                if "aromatic_chi2" in available else None,
            pucker_Q=get("pucker_Q")
                if "pucker_Q" in available else None,
            pucker_theta=get("pucker_theta")
                if "pucker_theta" in available else None,
        )

    # Water field (trajectory path — optional)
    water_field = None
    if "water_efield" in available:
        water_field = WaterFieldGroup(
            efield=get("water_efield"),
            efield_first=get("water_efield_first"),
            efg=get("water_efg"),
            efg_first=get("water_efg_first"),
            shell_counts=get("water_shell_counts"),
            efield_clamp_mask=get("water_efield_clamp_mask")
                if "water_efield_clamp_mask" in available else None,
            efield_clamp_scale=get("water_efield_clamp_scale")
                if "water_efield_clamp_scale" in available else None,
            efield_first_clamp_mask=get("water_efield_first_clamp_mask")
                if "water_efield_first_clamp_mask" in available else None,
            efield_first_clamp_scale=get("water_efield_first_clamp_scale")
                if "water_efield_first_clamp_scale" in available else None,
        )

    water_hbond_geometry = None
    if "water_hbond_candidates" in available:
        water_hbond_geometry = WaterHBondGeometryGroup(
            candidates=get("water_hbond_candidates"),
            counts=get("water_hbond_counts"),
            nearest=get("water_hbond_nearest"),
        )

    # Hydration shell (trajectory path — optional)
    hydration = None
    if "hydration_shell" in available:
        hydration = HydrationGroup(data=get("hydration_shell"))

    # Water polarisation — SASA-normal reference frame (trajectory path — optional)
    water_polarization = None
    if "water_polarization" in available:
        water_polarization = WaterPolarizationGroup(data=get("water_polarization"))

    # project-local QEq/EEQ-style charge-equilibration model with error-function coordination number, CN-dependent electronegativity shift, Gaussian self term, and Ohno-Klopman off-diagonal kernel; parameters are in-repo/project-local and are not a validated dftd4/multicharge port.
    # Optional: attach the group only when the producer emitted EEQ arrays.
    eeq = None
    if "eeq_charges" in available:
        eeq_coulomb = None
        if "eeq_coulomb_efg" in available:
            eeq_coulomb = EeqCoulombGroup(
                efg=get("eeq_coulomb_efg"),
                E=get("eeq_coulomb_E"),
                E_backbone=get("eeq_coulomb_E_backbone"),
                E_sidechain=get("eeq_coulomb_E_sidechain"),
                E_aromatic=get("eeq_coulomb_E_aromatic"),
                efg_backbone=get("eeq_coulomb_efg_backbone"),
                efg_sidechain=get("eeq_coulomb_efg_sidechain"),
                efg_aromatic=get("eeq_coulomb_efg_aromatic"),
                scalars=get("eeq_coulomb_scalars"),
                aromatic_E_proj=get("eeq_coulomb_aromatic_E_proj")
                    if "eeq_coulomb_aromatic_E_proj" in available else None,
                aromatic_n_src=get("eeq_coulomb_aromatic_n_src")
                    if "eeq_coulomb_aromatic_n_src" in available else None,
            )
        eeq = EeqGroup(
            charges=get("eeq_charges"),
            cn=get("eeq_cn"),
            chi_eff=get("eeq_chi_eff") if "eeq_chi_eff" in available else None,
            hardness=get("eeq_hardness") if "eeq_hardness" in available else None,
            coulomb=eeq_coulomb,
        )

    # Larsen H-bond term shielding (Larsen 2015): spatial enumeration of
    # backbone amide-H and Halpha donors against all typed acceptor classes.
    # Group attached if any larsen_hbond_* NPY is present.
    larsen_hbond = None
    larsen_hbond_stems = {
        "larsen_hbond_shielding",
        "larsen_hbond_1pHB_shielding",
        "larsen_hbond_2pHB_shielding",
        "larsen_hbond_1pHaB_shielding",
        "larsen_hbond_2pHaB_shielding",
        "larsen_hbond_diagnostic_CB_shielding",
        "larsen_hbond_water_term",
        "larsen_hbond_count",
        "larsen_corner_imputed",
        "larsen_imputed_pair_count",
        "larsen_sidechain_carbonyl_pair_count",
        "larsen_hbond_pairs_index",
        "larsen_hbond_pairs_geometry",
        "larsen_hbond_pairs_isotropic",
        "larsen_hbond_pairs",
    }
    if any(stem in available for stem in larsen_hbond_stems):
        larsen_hbond = LarsenHBondGroup(
            shielding=get("larsen_hbond_shielding")
                if "larsen_hbond_shielding" in available else None,
            pHB_1=get("larsen_hbond_1pHB_shielding")
                if "larsen_hbond_1pHB_shielding" in available else None,
            pHB_2=get("larsen_hbond_2pHB_shielding")
                if "larsen_hbond_2pHB_shielding" in available else None,
            pHaB_1=get("larsen_hbond_1pHaB_shielding")
                if "larsen_hbond_1pHaB_shielding" in available else None,
            pHaB_2=get("larsen_hbond_2pHaB_shielding")
                if "larsen_hbond_2pHaB_shielding" in available else None,
            diagnostic_CB=get("larsen_hbond_diagnostic_CB_shielding")
                if "larsen_hbond_diagnostic_CB_shielding" in available else None,
            water_term=get("larsen_hbond_water_term")
                if "larsen_hbond_water_term" in available else None,
            count=get("larsen_hbond_count")
                if "larsen_hbond_count" in available else None,
            corner_imputed=get("larsen_corner_imputed")
                if "larsen_corner_imputed" in available else None,
            imputed_pair_count=get("larsen_imputed_pair_count")
                if "larsen_imputed_pair_count" in available else None,
            sidechain_carbonyl_pair_count=get(
                "larsen_sidechain_carbonyl_pair_count")
                if "larsen_sidechain_carbonyl_pair_count" in available else None,
            pairs=LarsenHBondPairs(
                index=get("larsen_hbond_pairs_index"),
                geometry=get("larsen_hbond_pairs_geometry"),
                isotropic=get("larsen_hbond_pairs_isotropic"),
                compatibility=get("larsen_hbond_pairs"),
            ) if "larsen_hbond_pairs_index" in available else None,
        )

    larsen_sidechain_donor_audit = None
    if "larsen_sidechain_donor_atoms" in available:
        larsen_sidechain_donor_audit = LarsenSidechainDonorAuditGroup(
            atoms=get("larsen_sidechain_donor_atoms"),
            candidates=get("larsen_sidechain_donor_candidates"),
        )

    return Protein(
        protein_id=protein_id,
        n_atoms=n_atoms,
        pos=get("pos"),
        element=get("element"),
        residue_type=get("residue_type"),
        residue_index=get("residue_index"),
        biot_savart=biot_savart,
        haigh_mallion=haigh_mallion,
        pi_quadrupole=pi_quadrupole,
        ring_susceptibility=ring_susceptibility,
        dispersion=dispersion,
        ring_contributions=get("ring_contributions"),
        ring_direction_to_center=get("ring_direction_to_center"),
        ring_geometry=get("ring_geometry"),
        ring_pair_geometry=get("ring_pair_geometry"),
        enrichment=enrichment,
        charge_assignment=charge_assignment,
        spatial_index=spatial_index,
        molecular_graph=molecular_graph,
        mcconnell=mcconnell,
        sidechain_carbonyl_anisotropy=sidechain_carbonyl_anisotropy,
        coulomb=coulomb,
        hbond=hbond,
        dssp=get("dssp_backbone"),
        dssp_observed=get("dssp_observed"),
        dssp_ss8=get("dssp_ss8"),
        dssp_ppii=get("dssp_ppii"),
        dssp_hbond_energy=get("dssp_hbond_energy"),
        dssp_chi=get("dssp_chi"),
        sasa=get("atom_sasa"),
        sasa_normal=get("sasa_normal"),
        sasa_fraction=get("atom_sasa_fraction"),
        mopac=mopac,
        apbs=apbs,
        orca=orca,
        delta=delta,
        aimnet2=aimnet2,
        planar_geometry=planar_geometry,
        water_field=water_field,
        water_hbond_geometry=water_hbond_geometry,
        hydration=hydration,
        water_polarization=water_polarization,
        gromacs_energy=get("gromacs_energy"),
        bonded_energy=get("bonded_energy"),
        hm_welford=get("hm_welford"),
        mc_welford=get("mc_welford"),
        sasa_welford=get("sasa_welford"),
        eeq=eeq,
        category_info=category_info,
        topology=topology_group,
        larsen_hbond=larsen_hbond,
        larsen_sidechain_donor_audit=larsen_sidechain_donor_audit,
    )

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
    PerRingTypeT0,
    PerRingTypeT2,
    PerBondCategoryT2,
    RingCounts,
    McConnellScalars,
    CoulombScalars,
    HBondScalars,
    DsspScalars,
    MopacScalars,
    MopacGlobal,
    BondOrders,
    DeltaScalars,
    DeltaAPBS,
    DeltaRingProximity,
    AIMNet2Charges,
    AIMNet2AimEmbedding,
    AIMNet2Polarisability,
)
from ._ring import RingContributions, RingGeometry
from ._catalog import CATALOG


# ── Calculator group dataclasses ────────────────────────────────────


@dataclass(frozen=True)
class RingKernelGroup:
    """Shielding + per-type decomposition for a ring current calculator."""
    shielding: ShieldingTensor
    per_type_T0: PerRingTypeT0
    per_type_T2: PerRingTypeT2


@dataclass(frozen=True)
class BiotSavartGroup(RingKernelGroup):
    """Biot-Savart with additional B-field and ring counts."""
    total_B: VectorField = None
    ring_counts: RingCounts = None


@dataclass(frozen=True)
class McConnellGroup:
    """McConnell bond anisotropy: shielding + T2 per bond category + scalars."""
    shielding: ShieldingTensor
    category_T2: PerBondCategoryT2
    scalars: McConnellScalars


@dataclass(frozen=True)
class CoulombGroup:
    """Coulomb E-field: shielding + E-field vector + EFG decompositions."""
    shielding: ShieldingTensor
    E: VectorField
    efg_backbone: EFGTensor
    efg_aromatic: EFGTensor
    scalars: CoulombScalars


@dataclass(frozen=True)
class HBondGroup:
    shielding: ShieldingTensor
    scalars: HBondScalars


@dataclass(frozen=True)
class MopacCoreGroup:
    charges: np.ndarray
    scalars: MopacScalars
    bond_orders: BondOrders
    global_: MopacGlobal


@dataclass(frozen=True)
class MopacCoulombGroup:
    shielding: ShieldingTensor
    E: VectorField
    efg_backbone: EFGTensor
    efg_aromatic: EFGTensor
    scalars: CoulombScalars


@dataclass(frozen=True)
class MopacMcConnellGroup:
    shielding: ShieldingTensor
    category_T2: PerBondCategoryT2
    scalars: McConnellScalars


@dataclass(frozen=True)
class MopacGroup:
    core: MopacCoreGroup
    coulomb: MopacCoulombGroup
    mcconnell: MopacMcConnellGroup


@dataclass(frozen=True)
class APBSGroup:
    E: VectorField
    efg: EFGTensor


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
    charges: AIMNet2Charges
    aim: AIMNet2AimEmbedding
    efg: EFGTensor
    efg_aromatic: EFGTensor
    efg_backbone: EFGTensor
    # Polarisability fields are present only when the extraction was run
    # with AIMNet2 model loaded after the 2026-05-09 always-on promotion
    # of AIMNet2PolarisabilityResult. Old outputs leave these as None.
    polarisability: Optional[AIMNet2Polarisability] = None
    polarisability_scalar: Optional[np.ndarray] = None


@dataclass(frozen=True)
class PlanarGeometryGroup:
    """Planar and near-planar geometry outputs from PlanarGeometryResult.

    Axes are intentionally kept explicit because the result mixes atom,
    residue, aromatic-ring, and saturated-ring quantities:

    - ``pyramidalization`` is per atom, shape ``(N,)``.
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
    efg: EFGTensor                  # (N, 9) total water EFG
    efg_first: EFGTensor            # (N, 9) first-shell EFG
    shell_counts: np.ndarray        # (N, 2) [n_first, n_second]


@dataclass(frozen=True)
class HydrationGroup:
    """Per-atom hydration shell geometry."""
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
        """Dipole coherence ``|sum d_i| / n`` (N,)."""
        return self.data[:, 8]

    @property
    def shell_count(self) -> np.ndarray:
        """First-shell water count (N,)."""
        return self.data[:, 9]


@dataclass(frozen=True)
class EeqGroup:
    """EEQ geometry-dependent charges (Caldeweyher et al. 2019)."""
    charges: np.ndarray             # (N,) partial charges (elementary charges)
    cn: np.ndarray                  # (N,) coordination number


@dataclass(frozen=True)
class TripeptideGroup:
    """ProCS15 (Larsen 2015) tripeptide DFT shielding lookup.

    σ_BB^i (the central-residue backbone-shielding tensor) and
    Δσ_BB^{i±1} (the cap-side neighbour correction per Larsen Eq 3) are
    emitted by the C++ ``TripeptideBackboneShieldingResult`` and
    ``TripeptideNeighborShieldingResult`` calculators when the
    extraction is run against a host that has the ``tensorcs15``
    PostgreSQL replica available (DSN configured under
    ``[databases].tensorcs15`` in the runtime TOML).

    All fields are optional; the group is attached to a Protein only
    when at least one of the seven NPY files is present on disk.

    Per-atom convention:
    - ``bb_shielding`` carries the tensor on backbone N/CA/C/O/H/HA and
      central-residue sidechain atoms, NaN-filled elsewhere.
    - ``neighbor_shielding`` carries the summed Δσ_{i-1} + Δσ_{i+1}.
    - ``bb_residual_vec`` is the central-residue match residual
      (aligned_dft - protein_position) as a Vec3, NaN where the residue
      had no central match.
    - ``neighbor_residual_vec_prev/_next`` carry the per-direction cap
      residuals. NaN distinguishes "the i-1 (or i+1) direction did not
      contribute" from a coincidentally-zero residual.
    - ``bb_method_tag`` encodes the DFT engine that produced the row
      (1=OPBE Gaussian per Larsen, 2=PBE ORCA per the SER regen). 0
      means no match.
    """
    bb_shielding: Optional[ShieldingTensor] = None
    bb_residual_vec: Optional[VectorField] = None
    bb_match_distance: Optional[np.ndarray] = None
    bb_method_tag: Optional[np.ndarray] = None
    neighbor_shielding: Optional[ShieldingTensor] = None
    neighbor_residual_vec_prev: Optional[VectorField] = None
    neighbor_residual_vec_next: Optional[VectorField] = None


@dataclass(frozen=True)
class LarsenHBondGroup:
    """Larsen 2015 ProCS15 H-bond term shielding contributions.

    Per Larsen 2015 Eq 5, each H-bond pair contributes four terms (1° +
    2° on each of HB / HαB) plus a water-term offset Δσ_w on amide H
    atoms with no H-bond partner. Tensors are ppm in the protein lab
    frame (rotated from the canonical donor frame at calculator time)
    and emitted as SphericalTensor-packed (T0+T1+T2 = 9 cols).

    Methods accumulate side-by-side with the kernel-form ``HBondGroup``
    — both calculators cover overlapping physics (amide-H / backbone-O
    subset) but use different formulations (kernel × η vs grid lookup).
    Per-atom-type differences are themselves thesis-reportable. See
    ``feedback_methods_accumulate`` memory entry.

    Phase 1 (landed) covers amide-H donor / backbone-O acceptor only.
    Phase 2 adds Hα donors (via spatial search) and sidechain
    acceptors. The Phase 1 fields are populated; the Phase 2 fields
    (1pHaB / 2pHaB) are present-but-zero until Phase 2 lands.

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
    - ``water_term`` is 2.07 ppm isotropic on amide Hs that DSSP saw as
      solvent-exposed. A C-term-acceptor H-bond that our grid path
      couldn't process is NOT spuriously assigned the water term —
      the DSSP-paired bookkeeping is separate from the grid-paired
      bookkeeping (codex M2).
    - ``count`` counts H-bond pairs that contributed to this atom under
      any of the four Table 2 classes — the diagnostic CB does NOT
      increment it.
    """
    shielding: Optional[ShieldingTensor] = None
    pHB_1: Optional[ShieldingTensor] = None
    pHB_2: Optional[ShieldingTensor] = None
    pHaB_1: Optional[ShieldingTensor] = None
    pHaB_2: Optional[ShieldingTensor] = None
    diagnostic_CB: Optional[ShieldingTensor] = None
    water_term: Optional[np.ndarray] = None
    count: Optional[np.ndarray] = None


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
    pos: VectorField
    element: np.ndarray
    residue_type: np.ndarray
    residue_index: np.ndarray

    # Ring calculators
    biot_savart: BiotSavartGroup
    haigh_mallion: RingKernelGroup
    pi_quadrupole: RingKernelGroup
    dispersion: RingKernelGroup
    ring_susceptibility: ShieldingTensor

    # Per-ring sparse data
    ring_contributions: RingContributions = None
    ring_geometry: RingGeometry = None

    # Bond calculators
    mcconnell: McConnellGroup = None
    coulomb: CoulombGroup = None
    hbond: HBondGroup = None
    dssp: DsspScalars = None
    dssp_ss8: np.ndarray = None
    dssp_hbond_energy: np.ndarray = None
    dssp_chi: np.ndarray = None
    sasa: np.ndarray = None
    sasa_normal: Optional[VectorField] = None  # (N, 3) outward surface normal

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
    hydration: Optional[HydrationGroup] = None
    water_polarization: Optional[WaterPolarizationGroup] = None
    gromacs_energy: Optional[np.ndarray] = None

    # Geometry-dependent charges
    eeq: Optional[EeqGroup] = None

    # Tripeptide DFT shielding (ProCS15 / Larsen 2015) — emitted when
    # the extractor was run with the tensorcs15 Postgres DSN configured.
    tripeptide: Optional[TripeptideGroup] = None
    larsen_hbond: Optional[LarsenHBondGroup] = None


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
        total_B=get("bs_total_B"),
        ring_counts=get("bs_ring_counts"),
    )
    haigh_mallion = RingKernelGroup(
        shielding=get("hm_shielding"),
        per_type_T0=get("hm_per_type_T0"),
        per_type_T2=get("hm_per_type_T2"),
    )
    pi_quadrupole = RingKernelGroup(
        shielding=get("pq_shielding"),
        per_type_T0=get("pq_per_type_T0"),
        per_type_T2=get("pq_per_type_T2"),
    )
    dispersion = RingKernelGroup(
        shielding=get("disp_shielding"),
        per_type_T0=get("disp_per_type_T0"),
        per_type_T2=get("disp_per_type_T2"),
    )
    mcconnell = McConnellGroup(
        shielding=get("mc_shielding"),
        category_T2=get("mc_category_T2"),
        scalars=get("mc_scalars"),
    )
    coulomb = CoulombGroup(
        shielding=get("coulomb_shielding"),
        E=get("coulomb_E"),
        efg_backbone=get("coulomb_efg_backbone"),
        efg_aromatic=get("coulomb_efg_aromatic"),
        scalars=get("coulomb_scalars"),
    )
    hbond = HBondGroup(
        shielding=get("hbond_shielding"),
        scalars=get("hbond_scalars"),
    )

    # MOPAC (optional)
    mopac = None
    if "mopac_charges" in available:
        mopac = MopacGroup(
            core=MopacCoreGroup(
                charges=get("mopac_charges"),
                scalars=get("mopac_scalars"),
                bond_orders=get("mopac_bond_orders"),
                global_=get("mopac_global"),
            ),
            coulomb=MopacCoulombGroup(
                shielding=get("mopac_coulomb_shielding"),
                E=get("mopac_coulomb_E"),
                efg_backbone=get("mopac_coulomb_efg_backbone"),
                efg_aromatic=get("mopac_coulomb_efg_aromatic"),
                scalars=get("mopac_coulomb_scalars"),
            ),
            mcconnell=MopacMcConnellGroup(
                shielding=get("mopac_mc_shielding"),
                category_T2=get("mopac_mc_category_T2"),
                scalars=get("mopac_mc_scalars"),
            ),
        )

    # APBS (optional)
    apbs = None
    if "apbs_E" in available:
        apbs = APBSGroup(E=get("apbs_E"), efg=get("apbs_efg"))

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

    # AIMNet2 (optional)
    aimnet2 = None
    if "aimnet2_charges" in available:
        aimnet2 = AIMNet2Group(
            charges=get("aimnet2_charges"),
            aim=get("aimnet2_aim"),
            efg=get("aimnet2_efg"),
            efg_aromatic=get("aimnet2_efg_aromatic"),
            efg_backbone=get("aimnet2_efg_backbone"),
            polarisability=get("aimnet2_polarisability")
                if "aimnet2_polarisability" in available else None,
            polarisability_scalar=get("aimnet2_polarisability_scalar")
                if "aimnet2_polarisability_scalar" in available else None,
        )

    # Planar geometry (optional in pre-2026-05-09 outputs and in unusual
    # fixtures where the topology substrate is not populated).
    planar_geometry = None
    planar_stems = {
        "pyramidalization",
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
        )

    # Hydration shell (trajectory path — optional)
    hydration = None
    if "hydration_shell" in available:
        hydration = HydrationGroup(data=get("hydration_shell"))

    # Water polarisation — SASA-normal reference frame (trajectory path — optional)
    water_polarization = None
    if "water_polarization" in available:
        water_polarization = WaterPolarizationGroup(data=get("water_polarization"))

    # EEQ charges (Caldeweyher 2019 — optional)
    eeq = None
    if "eeq_charges" in available:
        eeq = EeqGroup(
            charges=get("eeq_charges"),
            cn=get("eeq_cn"),
        )

    # Tripeptide DFT shielding (ProCS15 / Larsen 2015) — attached when
    # any of the seven tripeptide NPYs is present. Individual fields
    # are wrapped (or left None) based on per-stem availability so a
    # partial output (e.g., BB calculator ran but Neighbor did not) is
    # still consumable.
    tripeptide = None
    tripeptide_stems = {
        "tripeptide_bb_shielding",
        "tripeptide_bb_residual_vec",
        "tripeptide_bb_match_distance",
        "tripeptide_bb_method_tag",
        "tripeptide_neighbor_shielding",
        "tripeptide_neighbor_residual_vec_prev",
        "tripeptide_neighbor_residual_vec_next",
    }
    if any(stem in available for stem in tripeptide_stems):
        tripeptide = TripeptideGroup(
            bb_shielding=get("tripeptide_bb_shielding")
                if "tripeptide_bb_shielding" in available else None,
            bb_residual_vec=get("tripeptide_bb_residual_vec")
                if "tripeptide_bb_residual_vec" in available else None,
            bb_match_distance=get("tripeptide_bb_match_distance")
                if "tripeptide_bb_match_distance" in available else None,
            bb_method_tag=get("tripeptide_bb_method_tag")
                if "tripeptide_bb_method_tag" in available else None,
            neighbor_shielding=get("tripeptide_neighbor_shielding")
                if "tripeptide_neighbor_shielding" in available else None,
            neighbor_residual_vec_prev=get("tripeptide_neighbor_residual_vec_prev")
                if "tripeptide_neighbor_residual_vec_prev" in available else None,
            neighbor_residual_vec_next=get("tripeptide_neighbor_residual_vec_next")
                if "tripeptide_neighbor_residual_vec_next" in available else None,
        )

    # Larsen H-bond term shielding (Larsen 2015) — Phase 1 covers
    # amide-H / backbone-O subset; Phase 2 adds Hα + sidechain. Group
    # attached if any larsen_hbond_* NPY is present.
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
        dispersion=dispersion,
        ring_susceptibility=get("ringchi_shielding"),
        ring_contributions=get("ring_contributions"),
        ring_geometry=get("ring_geometry"),
        mcconnell=mcconnell,
        coulomb=coulomb,
        hbond=hbond,
        dssp=get("dssp_backbone"),
        dssp_ss8=get("dssp_ss8"),
        dssp_hbond_energy=get("dssp_hbond_energy"),
        dssp_chi=get("dssp_chi"),
        sasa=get("atom_sasa"),
        sasa_normal=get("sasa_normal"),
        mopac=mopac,
        apbs=apbs,
        orca=orca,
        delta=delta,
        aimnet2=aimnet2,
        planar_geometry=planar_geometry,
        water_field=water_field,
        hydration=hydration,
        water_polarization=water_polarization,
        gromacs_energy=get("gromacs_energy"),
        eeq=eeq,
        category_info=category_info,
        topology=topology_group,
        tripeptide=tripeptide,
        larsen_hbond=larsen_hbond,
    )

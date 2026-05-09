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

    # Per-atom invariant categorical record (CategoryInfoProjection).
    category_info: Optional[CategoryInfo] = None

    # Explicit solvent (trajectory path only)
    water_field: Optional[WaterFieldGroup] = None
    hydration: Optional[HydrationGroup] = None
    water_polarization: Optional[WaterPolarizationGroup] = None
    gromacs_energy: Optional[np.ndarray] = None

    # Geometry-dependent charges
    eeq: Optional[EeqGroup] = None


# ── Loader ──────────────────────────────────────────────────────────


def _wrap(spec, data: np.ndarray):
    """Wrap raw array in its typed class, or return as-is for np.ndarray."""
    if spec.wrapper is np.ndarray:
        return data
    return spec.wrapper(data)


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
        water_field=water_field,
        hydration=hydration,
        water_polarization=water_polarization,
        gromacs_energy=get("gromacs_energy"),
        eeq=eeq,
        category_info=category_info,
    )

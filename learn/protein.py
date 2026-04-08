"""
Typed protein loader for calibration features.

Reads a protein feature directory (51 NPY arrays from --mutant extraction)
and returns a ProteinFeatures object with every array wrapped in its
correct type, grouped by physical calculator.

    from learn.protein import load_protein

    p = load_protein("learn/features/CalibrationExtractionTest/A0A7C5FAR6")
    p.biot_savart.shielding.T2          # (543, 5)
    p.mcconnell.category_T2.for_category("CO_nearest")  # (543, 5)
    p.delta.scalars.matched_mask        # (543,) bool
    p.mopac.bond_orders.to_dense(p.n_atoms)  # (543, 543)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

from features import (
    BondOrders,
    CoulombScalars,
    DeltaAPBS,
    DeltaRingProximity,
    DeltaScalars,
    DsspBackbone,
    EFGTensor,
    HBondScalars,
    McConnellScalars,
    MopacAtomScalars,
    MopacGlobal,
    PerBondCategoryT2,
    PerRingTypeT0,
    PerRingTypeT2,
    RingCounts,
    ShieldingTensor,
    VectorField,
)


# ── Calculator group dataclasses ─────────────────────────────────────

@dataclass
class BiotSavart:
    """Biot-Savart ring current calculator outputs."""
    shielding: ShieldingTensor
    per_type_T0: PerRingTypeT0
    per_type_T2: PerRingTypeT2
    ring_counts: RingCounts
    total_B: VectorField


@dataclass
class McConnell:
    """McConnell bond anisotropy calculator outputs."""
    shielding: ShieldingTensor
    category_T2: PerBondCategoryT2
    scalars: McConnellScalars


@dataclass
class Coulomb:
    """Coulomb E-field calculator outputs (ff14SB charges)."""
    shielding: ShieldingTensor
    E: VectorField
    efg_backbone: EFGTensor
    efg_aromatic: EFGTensor
    scalars: CoulombScalars


@dataclass
class HBond:
    """Hydrogen bond calculator outputs."""
    shielding: ShieldingTensor
    scalars: HBondScalars


@dataclass
class MopacCore:
    """MOPAC electronic structure (PM7)."""
    charges: np.ndarray          # (N,) Mulliken charges
    scalars: MopacAtomScalars    # (N, 3) charge, s_pop, p_pop
    bond_orders: BondOrders      # (M, 3) sparse Wiberg
    global_: MopacGlobal         # (4,) graph-level


@dataclass
class MopacCoulomb:
    """Coulomb calculator re-evaluated with MOPAC charges."""
    shielding: ShieldingTensor
    E: VectorField
    efg_backbone: EFGTensor
    efg_aromatic: EFGTensor
    scalars: CoulombScalars


@dataclass
class MopacMcConnell:
    """McConnell calculator weighted by MOPAC Wiberg bond orders."""
    shielding: ShieldingTensor
    category_T2: PerBondCategoryT2
    scalars: McConnellScalars


@dataclass
class Mopac:
    """All MOPAC-derived results."""
    core: MopacCore
    coulomb: MopacCoulomb
    mcconnell: MopacMcConnell


@dataclass
class Orca:
    """DFT reference shielding tensors."""
    total: ShieldingTensor
    diamagnetic: ShieldingTensor
    paramagnetic: ShieldingTensor


@dataclass
class APBS:
    """APBS Poisson-Boltzmann electrostatics."""
    E: VectorField
    efg: EFGTensor


@dataclass
class Delta:
    """Mutation delta (WT - ALA) results."""
    shielding: ShieldingTensor
    scalars: DeltaScalars
    apbs: Optional[DeltaAPBS]
    ring_proximity: DeltaRingProximity


# ── Top-level protein container ──────────────────────────────────────

@dataclass
class ProteinFeatures:
    """All features for one calibration protein.

    Every field is typed — no raw arrays escape except identity data.
    Use .n_atoms for the atom count.  Optional groups are None when
    the extractor did not produce them.
    """
    protein_id: str
    n_atoms: int

    # Identity
    pos: np.ndarray              # (N, 3) float64  Angstroms
    element: np.ndarray          # (N,) int32
    residue_type: np.ndarray     # (N,) int32
    residue_index: np.ndarray    # (N,) int32

    # Classical calculators
    biot_savart: BiotSavart
    mcconnell: McConnell
    coulomb: Coulomb
    hbond: HBond
    haigh_mallion: "HaighMallion"
    pople_karplus: "PopleKarplus"
    dispersion: "Dispersion"
    ring_susceptibility: ShieldingTensor
    dssp: DsspBackbone

    # MOPAC-derived
    mopac: Optional[Mopac] = None

    # APBS electrostatics
    apbs: Optional[APBS] = None

    # DFT reference
    orca: Optional[Orca] = None

    # Mutation delta
    delta: Optional[Delta] = None


@dataclass
class HaighMallion:
    """Haigh-Mallion ring current calculator outputs."""
    shielding: ShieldingTensor
    per_type_T0: PerRingTypeT0
    per_type_T2: PerRingTypeT2


@dataclass
class PopleKarplus:
    """Pople-Karplus ring current calculator outputs."""
    shielding: ShieldingTensor
    per_type_T0: PerRingTypeT0
    per_type_T2: PerRingTypeT2


@dataclass
class Dispersion:
    """Dispersion shielding calculator outputs."""
    shielding: ShieldingTensor
    per_type_T0: PerRingTypeT0
    per_type_T2: PerRingTypeT2


# ── Registry: filename → (group_path, wrapper_factory) ───────────────
#
# Each entry maps an NPY stem to:
#   - dotted path into ProteinFeatures (used by the loader)
#   - wrapper class or callable that wraps the raw ndarray
#   - whether the file is required (True) or optional (False)

_REG: list[tuple[str, str, type, bool]] = [
    # Identity
    ("pos",              "pos",                               np.ndarray, True),
    ("element",          "element",                           np.ndarray, True),
    ("residue_type",     "residue_type",                      np.ndarray, True),
    ("residue_index",    "residue_index",                     np.ndarray, True),

    # Biot-Savart
    ("bs_shielding",     "biot_savart.shielding",             ShieldingTensor, True),
    ("bs_per_type_T0",   "biot_savart.per_type_T0",           PerRingTypeT0,   True),
    ("bs_per_type_T2",   "biot_savart.per_type_T2",           PerRingTypeT2,   True),
    ("bs_ring_counts",   "biot_savart.ring_counts",           RingCounts,      True),
    ("bs_total_B",       "biot_savart.total_B",               VectorField,     True),

    # McConnell
    ("mc_shielding",     "mcconnell.shielding",               ShieldingTensor, True),
    ("mc_category_T2",   "mcconnell.category_T2",             PerBondCategoryT2, True),
    ("mc_scalars",       "mcconnell.scalars",                 McConnellScalars, True),

    # Coulomb
    ("coulomb_shielding",       "coulomb.shielding",          ShieldingTensor, True),
    ("coulomb_E",               "coulomb.E",                  VectorField,     True),
    ("coulomb_efg_backbone",    "coulomb.efg_backbone",       EFGTensor,       True),
    ("coulomb_efg_aromatic",    "coulomb.efg_aromatic",       EFGTensor,       True),
    ("coulomb_scalars",         "coulomb.scalars",            CoulombScalars,  True),

    # H-bond
    ("hbond_shielding",  "hbond.shielding",                   ShieldingTensor, True),
    ("hbond_scalars",    "hbond.scalars",                     HBondScalars,    True),

    # Haigh-Mallion
    ("hm_shielding",     "haigh_mallion.shielding",           ShieldingTensor, True),
    ("hm_per_type_T0",   "haigh_mallion.per_type_T0",         PerRingTypeT0,   True),
    ("hm_per_type_T2",   "haigh_mallion.per_type_T2",         PerRingTypeT2,   True),

    # Pople-Karplus
    ("pq_shielding",     "pople_karplus.shielding",           ShieldingTensor, True),
    ("pq_per_type_T0",   "pople_karplus.per_type_T0",         PerRingTypeT0,   True),
    ("pq_per_type_T2",   "pople_karplus.per_type_T2",         PerRingTypeT2,   True),

    # Dispersion
    ("disp_shielding",   "dispersion.shielding",              ShieldingTensor, True),
    ("disp_per_type_T0", "dispersion.per_type_T0",            PerRingTypeT0,   True),
    ("disp_per_type_T2", "dispersion.per_type_T2",            PerRingTypeT2,   True),

    # Ring susceptibility
    ("ringchi_shielding", "ring_susceptibility",              ShieldingTensor, True),

    # DSSP
    ("dssp_backbone",    "dssp",                              DsspBackbone,    True),

    # MOPAC core
    ("mopac_charges",    "mopac.core.charges",                np.ndarray,      False),
    ("mopac_scalars",    "mopac.core.scalars",                MopacAtomScalars, False),
    ("mopac_bond_orders","mopac.core.bond_orders",            BondOrders,      False),
    ("mopac_global",     "mopac.core.global_",                MopacGlobal,     False),

    # MOPAC Coulomb
    ("mopac_coulomb_shielding",    "mopac.coulomb.shielding",     ShieldingTensor, False),
    ("mopac_coulomb_E",            "mopac.coulomb.E",             VectorField,     False),
    ("mopac_coulomb_efg_backbone", "mopac.coulomb.efg_backbone",  EFGTensor,       False),
    ("mopac_coulomb_efg_aromatic", "mopac.coulomb.efg_aromatic",  EFGTensor,       False),
    ("mopac_coulomb_scalars",      "mopac.coulomb.scalars",       CoulombScalars,  False),

    # MOPAC McConnell
    ("mopac_mc_shielding",     "mopac.mcconnell.shielding",      ShieldingTensor,     False),
    ("mopac_mc_category_T2",   "mopac.mcconnell.category_T2",    PerBondCategoryT2,   False),
    ("mopac_mc_scalars",       "mopac.mcconnell.scalars",        McConnellScalars,    False),

    # APBS
    ("apbs_E",           "apbs.E",                            VectorField,     False),
    ("apbs_efg",         "apbs.efg",                          EFGTensor,       False),

    # Orca DFT reference
    ("orca_total",       "orca.total",                        ShieldingTensor, False),
    ("orca_diamagnetic", "orca.diamagnetic",                  ShieldingTensor, False),
    ("orca_paramagnetic","orca.paramagnetic",                 ShieldingTensor, False),

    # Delta (mutation)
    ("delta_shielding",       "delta.shielding",              ShieldingTensor,      False),
    ("delta_scalars",         "delta.scalars",                DeltaScalars,         False),
    ("delta_apbs",            "delta.apbs",                   DeltaAPBS,            False),
    ("delta_ring_proximity",  "delta.ring_proximity",         DeltaRingProximity,   False),
]

# Build lookup dict
_REGISTRY: dict[str, tuple[str, type, bool]] = {
    stem: (path, cls, req) for stem, path, cls, req in _REG
}


# ── Loader ───────────────────────────────────────────────────────────

def _wrap(cls: type, data: np.ndarray):
    """Wrap raw array in its typed class, or return as-is for np.ndarray."""
    if cls is np.ndarray:
        return data
    return cls(data)


def _build_group(group_cls: type, attrs: dict):
    """Construct a dataclass from a dict of {field_name: value}."""
    return group_cls(**attrs)


def load_protein(path: str | Path) -> ProteinFeatures:
    """Load a protein feature directory into a fully typed ProteinFeatures.

    Args:
        path: directory containing *.npy files from nmr_extract --mutant

    Raises:
        FileNotFoundError: if a required file is missing
        ValueError: if an array has unexpected shape
    """
    path = Path(path)
    protein_id = path.name

    # Load all available arrays
    available: dict[str, np.ndarray] = {}
    for npy in path.glob("*.npy"):
        available[npy.stem] = np.load(npy)

    # Check for unregistered files (data on the floor)
    unregistered = set(available.keys()) - set(_REGISTRY.keys())
    if unregistered:
        raise ValueError(
            f"Unregistered NPY files in {protein_id}: {sorted(unregistered)}. "
            f"Add them to the registry in protein.py")

    # Check required files
    for stem, (dotpath, cls, required) in _REGISTRY.items():
        if required and stem not in available:
            raise FileNotFoundError(
                f"Required file {stem}.npy missing for {protein_id}")

    # Wrap each array in its typed class
    wrapped: dict[str, object] = {}
    for stem, data in available.items():
        dotpath, cls, _ = _REGISTRY[stem]
        wrapped[dotpath] = _wrap(cls, data)

    # Determine n_atoms from pos
    n_atoms = available["pos"].shape[0]

    # Assemble groups bottom-up
    def get(dotpath: str, default=None):
        return wrapped.get(dotpath, default)

    # Classical calculators (always present)
    biot_savart = BiotSavart(
        shielding=get("biot_savart.shielding"),
        per_type_T0=get("biot_savart.per_type_T0"),
        per_type_T2=get("biot_savart.per_type_T2"),
        ring_counts=get("biot_savart.ring_counts"),
        total_B=get("biot_savart.total_B"),
    )
    mcconnell = McConnell(
        shielding=get("mcconnell.shielding"),
        category_T2=get("mcconnell.category_T2"),
        scalars=get("mcconnell.scalars"),
    )
    coulomb = Coulomb(
        shielding=get("coulomb.shielding"),
        E=get("coulomb.E"),
        efg_backbone=get("coulomb.efg_backbone"),
        efg_aromatic=get("coulomb.efg_aromatic"),
        scalars=get("coulomb.scalars"),
    )
    hbond = HBond(
        shielding=get("hbond.shielding"),
        scalars=get("hbond.scalars"),
    )
    haigh_mallion = HaighMallion(
        shielding=get("haigh_mallion.shielding"),
        per_type_T0=get("haigh_mallion.per_type_T0"),
        per_type_T2=get("haigh_mallion.per_type_T2"),
    )
    pople_karplus = PopleKarplus(
        shielding=get("pople_karplus.shielding"),
        per_type_T0=get("pople_karplus.per_type_T0"),
        per_type_T2=get("pople_karplus.per_type_T2"),
    )
    dispersion = Dispersion(
        shielding=get("dispersion.shielding"),
        per_type_T0=get("dispersion.per_type_T0"),
        per_type_T2=get("dispersion.per_type_T2"),
    )

    # MOPAC (optional)
    mopac = None
    if "mopac_charges" in available:
        mopac_core = MopacCore(
            charges=get("mopac.core.charges"),
            scalars=get("mopac.core.scalars"),
            bond_orders=get("mopac.core.bond_orders"),
            global_=get("mopac.core.global_"),
        )
        mopac_coulomb = MopacCoulomb(
            shielding=get("mopac.coulomb.shielding"),
            E=get("mopac.coulomb.E"),
            efg_backbone=get("mopac.coulomb.efg_backbone"),
            efg_aromatic=get("mopac.coulomb.efg_aromatic"),
            scalars=get("mopac.coulomb.scalars"),
        )
        mopac_mc = MopacMcConnell(
            shielding=get("mopac.mcconnell.shielding"),
            category_T2=get("mopac.mcconnell.category_T2"),
            scalars=get("mopac.mcconnell.scalars"),
        )
        mopac = Mopac(core=mopac_core, coulomb=mopac_coulomb, mcconnell=mopac_mc)

    # APBS (optional)
    apbs = None
    if "apbs_E" in available:
        apbs = APBS(E=get("apbs.E"), efg=get("apbs.efg"))

    # Orca DFT (optional)
    orca = None
    if "orca_total" in available:
        orca = Orca(
            total=get("orca.total"),
            diamagnetic=get("orca.diamagnetic"),
            paramagnetic=get("orca.paramagnetic"),
        )

    # Delta (optional)
    delta = None
    if "delta_shielding" in available:
        delta = Delta(
            shielding=get("delta.shielding"),
            scalars=get("delta.scalars"),
            apbs=get("delta.apbs"),
            ring_proximity=get("delta.ring_proximity"),
        )

    return ProteinFeatures(
        protein_id=protein_id,
        n_atoms=n_atoms,
        pos=get("pos"),
        element=get("element"),
        residue_type=get("residue_type"),
        residue_index=get("residue_index"),
        biot_savart=biot_savart,
        mcconnell=mcconnell,
        coulomb=coulomb,
        hbond=hbond,
        haigh_mallion=haigh_mallion,
        pople_karplus=pople_karplus,
        dispersion=dispersion,
        ring_susceptibility=get("ring_susceptibility"),
        dssp=get("dssp"),
        mopac=mopac,
        apbs=apbs,
        orca=orca,
        delta=delta,
    )


def list_proteins(run_dir: str | Path) -> list[str]:
    """List protein IDs with extracted features in a run directory."""
    run_dir = Path(run_dir)
    return sorted(d.name for d in run_dir.iterdir()
                  if d.is_dir() and (d / "pos.npy").exists())

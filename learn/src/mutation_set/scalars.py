"""Scalar feature assembly from nmr_extract SDK.

Each feature group is a named function returning (name, array).
Positions are tracked dynamically — no hardcoded column indices.
One-hot columns identified by name, not position.

The final scalar vector includes kernel normalization scales so the
MLP knows the magnitude that per-protein normalization stripped.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from nmr_extract import RingType
from .config import Config


# ── Helpers ──────────────────────────────────────────────────────────

ELEMENTS = [1, 6, 7, 8, 16]  # H, C, N, O, S


def _one_hot(values, n_classes: int) -> np.ndarray:
    M = len(values)
    oh = np.zeros((M, n_classes), dtype=np.float64)
    for i, v in enumerate(values):
        if 0 <= v < n_classes:
            oh[i, v] = 1.0
    return oh


def _sentinel_to_inv(arr, sentinel: float, epsilon: float) -> np.ndarray:
    return np.where(arr > sentinel, 0.0, 1.0 / (arr + epsilon))


# ── Feature groups ───────────────────────────────────────────────────

@dataclass
class ScalarBlock:
    """A named group of scalar features."""
    name: str
    data: np.ndarray       # (M, width)
    categorical: bool = False  # skip z-score if True

    @property
    def width(self) -> int:
        return self.data.shape[1]


def _identity_features(p, idx, M) -> list[ScalarBlock]:
    """Element one-hot (5), residue type one-hot (20), ring counts (4)."""
    elem = p.element[idx]
    elem_oh = np.zeros((M, len(ELEMENTS)), dtype=np.float64)
    for i, e in enumerate(ELEMENTS):
        elem_oh[:, i] = (elem == e).astype(float)

    return [
        ScalarBlock("element", elem_oh, categorical=True),
        ScalarBlock("residue_type", _one_hot(p.residue_type[idx], 20),
                    categorical=True),
        ScalarBlock("ring_counts", p.biot_savart.ring_counts.data[idx] / 10.0),
    ]


def _calculator_scalars(p, idx, M, cfg) -> list[ScalarBlock]:
    """McConnell (6), Coulomb (4), HBond (3) scalar summaries."""
    sent, eps = cfg.sentinel_distance, cfg.inverse_epsilon

    mc = p.mcconnell.scalars.data[idx].copy()
    mc[:, 4] = _sentinel_to_inv(mc[:, 4], sent, eps)
    mc[:, 5] = _sentinel_to_inv(mc[:, 5], sent, eps)

    hb = p.hbond.scalars.data[idx].copy()
    hb[:, 0] = _sentinel_to_inv(hb[:, 0], sent, eps)

    return [
        ScalarBlock("mcconnell", mc),
        ScalarBlock("coulomb", p.coulomb.scalars.data[idx]),
        ScalarBlock("hbond", hb),
    ]


def _delta_features(p, idx, M, cfg) -> list[ScalarBlock]:
    """Delta metadata (4): 1/nearest_ring_dist, charge deltas, match distance."""
    ds = p.delta.scalars
    delta = np.column_stack([
        _sentinel_to_inv(ds.nearest_removed_ring_dist[idx],
                         cfg.sentinel_distance, cfg.inverse_epsilon),
        ds.delta_partial_charge[idx],
        ds.delta_mopac_charge[idx],
        ds.match_distance[idx],
    ])
    return [ScalarBlock("delta", delta)]


def _mopac_charge(p, idx, M) -> list[ScalarBlock]:
    """MOPAC Mulliken charge (1)."""
    q = np.zeros((M, 1), dtype=np.float64)
    if p.mopac:
        q[:, 0] = p.mopac.core.charges[idx]
    return [ScalarBlock("mopac_charge", q)]


def _t1_magnitudes(p, idx, M) -> list[ScalarBlock]:
    """T1 (antisymmetric) magnitudes from all calculators (10)."""
    shieldings = [
        p.biot_savart.shielding, p.haigh_mallion.shielding,
        p.pi_quadrupole.shielding, p.dispersion.shielding,
        p.mcconnell.shielding, p.coulomb.shielding,
        p.hbond.shielding, p.ring_susceptibility,
    ]
    t1 = np.zeros((M, 10), dtype=np.float64)
    for i, s in enumerate(shieldings):
        t1[:, i] = np.linalg.norm(s.T1[idx], axis=-1)
    if p.mopac:
        t1[:, 8] = np.linalg.norm(p.mopac.coulomb.shielding.T1[idx], axis=-1)
        t1[:, 9] = np.linalg.norm(p.mopac.mcconnell.shielding.T1[idx], axis=-1)
    return [ScalarBlock("t1_magnitudes", t1)]


def _ring_proximity(p, idx, M, cfg) -> list[ScalarBlock]:
    """Per-ring geometry + dispersion scalars for top-K nearest rings.

    Per ring: [mcconnell_factor, exp_decay, 1/dist, z, rho, theta,
               disp_scalar, disp_contacts] + n_rings count.

    disp_scalar encodes vertex proximity asymmetry: same center-distance
    but different disp_scalar means the atom sees an edge vs a face.
    disp_contacts is the effective solid angle (how many vertices in range).
    """
    K = cfg.top_k_rings
    COLS = 8  # 6 geometry + disp_scalar + disp_contacts
    rp = np.zeros((M, K * COLS + 1), dtype=np.float64)

    rc = p.ring_contributions
    if rc is not None and rc.n_pairs > 0:
        for j in range(M):
            atom_rc = rc.for_atom(idx[j])
            if atom_rc.n_pairs == 0:
                continue
            order = np.argsort(atom_rc.distance)
            for ki, ri in enumerate(order[:K]):
                base = ki * COLS
                rp[j, base]     = atom_rc.mcconnell_factor[ri]
                rp[j, base + 1] = atom_rc.exp_decay[ri]
                rp[j, base + 2] = 1.0 / (atom_rc.distance[ri] + cfg.inverse_epsilon)
                rp[j, base + 3] = atom_rc.z[ri]
                rp[j, base + 4] = atom_rc.rho[ri]
                rp[j, base + 5] = atom_rc.theta[ri]
                rp[j, base + 6] = atom_rc.disp_scalar[ri]
                rp[j, base + 7] = atom_rc.disp_contacts[ri]
            rp[j, -1] = atom_rc.n_pairs

    return [ScalarBlock("ring_proximity", rp)]


def _bond_orders(p, idx, M, cfg) -> list[ScalarBlock]:
    """Per-atom bond order features (4): max, count, total, n_aromatic."""
    bo = np.zeros((M, 4), dtype=np.float64)
    if p.mopac:
        dense = p.mopac.core.bond_orders.to_dense(p.n_atoms)[idx]
        bo[:, 0] = dense.max(axis=1)
        bo[:, 1] = (dense > cfg.bond_order_floor).sum(axis=1).astype(float)
        bo[:, 2] = dense.sum(axis=1)
        bo[:, 3] = (dense > cfg.bond_order_pi).sum(axis=1).astype(float)
    return [ScalarBlock("bond_orders", bo)]


def _mutation_identity(p, idx, M) -> list[ScalarBlock]:
    """Which aromatic residue types are present (4 one-hot)."""
    mutid = np.zeros((M, 4), dtype=np.float64)
    active = np.abs(p.biot_savart.per_type_T0.data[idx]).sum(axis=0) > 1e-10
    if active[RingType.PHE]:       mutid[:, 0] = 1.0
    if active[RingType.TYR]:       mutid[:, 1] = 1.0
    if any(active[t] for t in [RingType.TRP_benzene,
           RingType.TRP_pyrrole, RingType.TRP_perimeter]):
                                   mutid[:, 2] = 1.0
    if any(active[t] for t in [RingType.HIS, RingType.HID, RingType.HIE]):
                                   mutid[:, 3] = 1.0
    return [ScalarBlock("mutation_id", mutid, categorical=True)]


def _per_type_t0(p, idx, M) -> list[ScalarBlock]:
    """Per-ring-type T0 from BS/HM/Disp (24 total)."""
    return [
        ScalarBlock("t0_bs", p.biot_savart.per_type_T0.data[idx]),
        ScalarBlock("t0_hm", p.haigh_mallion.per_type_T0.data[idx]),
        ScalarBlock("t0_disp", p.dispersion.per_type_T0.data[idx]),
    ]


def _mopac_electronic(p, idx, M) -> list[ScalarBlock]:
    """MOPAC s/p orbital populations (2)."""
    elec = np.zeros((M, 2), dtype=np.float64)
    if p.mopac:
        elec[:, 0] = p.mopac.core.scalars.s_pop[idx]
        elec[:, 1] = p.mopac.core.scalars.p_pop[idx]
    return [ScalarBlock("mopac_electronic", elec)]


def _t0_magnitudes(p, idx, M) -> list[ScalarBlock]:
    """T0 isotropic magnitudes (6): |MC|, |Coulomb|, |BS|, |HBond|, |MopacCoul|, |MopacMC|."""
    t0 = np.zeros((M, 6), dtype=np.float64)
    t0[:, 0] = np.abs(p.mcconnell.shielding.isotropic[idx])
    t0[:, 1] = np.abs(p.coulomb.shielding.isotropic[idx])
    t0[:, 2] = np.abs(p.biot_savart.shielding.isotropic[idx])
    t0[:, 3] = np.abs(p.hbond.shielding.isotropic[idx])
    if p.mopac:
        t0[:, 4] = np.abs(p.mopac.coulomb.shielding.isotropic[idx])
        t0[:, 5] = np.abs(p.mopac.mcconnell.shielding.isotropic[idx])
    return [ScalarBlock("t0_magnitudes", t0)]


def _mopac_calculator_scalars(p, idx, M, cfg) -> list[ScalarBlock]:
    """MOPAC Coulomb (4) and McConnell (6) scalars."""
    coul = np.zeros((M, 4), dtype=np.float64)
    mc = np.zeros((M, 6), dtype=np.float64)
    if p.mopac:
        coul = p.mopac.coulomb.scalars.data[idx]
        mc = p.mopac.mcconnell.scalars.data[idx].copy()
        mc[:, 4] = _sentinel_to_inv(mc[:, 4], cfg.sentinel_distance,
                                    cfg.inverse_epsilon)
        mc[:, 5] = _sentinel_to_inv(mc[:, 5], cfg.sentinel_distance,
                                    cfg.inverse_epsilon)
    return [
        ScalarBlock("mopac_coulomb", coul),
        ScalarBlock("mopac_mcconnell", mc),
    ]


# ── Assembly ─────────────────────────────────────────────────────────

@dataclass(frozen=True)
class ScalarLayout:
    """Tracks scalar block names, widths, and categorical flags."""
    names: list[str]
    widths: list[int]
    categorical_mask: list[bool]
    total: int

    @property
    def categorical_columns(self) -> list[int]:
        """Column indices that should skip z-score normalization."""
        cols = []
        offset = 0
        for width, is_cat in zip(self.widths, self.categorical_mask):
            if is_cat:
                cols.extend(range(offset, offset + width))
            offset += width
        return cols


def assemble_scalars(protein, idx: np.ndarray, kernel_scales: np.ndarray,
                     cfg: Config) -> tuple[np.ndarray, ScalarLayout]:
    """Build scalar feature vector and its layout.

    Returns (scalars_array, layout) where layout tracks names and
    which columns are categorical (skip z-score).
    """
    M = len(idx)
    p = protein
    sc = cfg.scalars

    blocks: list[ScalarBlock] = []
    blocks += _identity_features(p, idx, M)
    blocks += _calculator_scalars(p, idx, M, sc)
    blocks += _delta_features(p, idx, M, sc)
    blocks += _mopac_charge(p, idx, M)
    blocks += _t1_magnitudes(p, idx, M)
    blocks += _ring_proximity(p, idx, M, sc)
    blocks += _bond_orders(p, idx, M, sc)
    blocks += _mutation_identity(p, idx, M)
    blocks += _per_type_t0(p, idx, M)
    blocks += _mopac_electronic(p, idx, M)
    blocks += _t0_magnitudes(p, idx, M)
    blocks += _mopac_calculator_scalars(p, idx, M, sc)
    blocks.append(ScalarBlock("kernel_scales", np.tile(kernel_scales, (M, 1))))

    scalars = np.concatenate([b.data for b in blocks], axis=1)
    layout = ScalarLayout(
        names=[b.name for b in blocks],
        widths=[b.width for b in blocks],
        categorical_mask=[b.categorical for b in blocks],
        total=scalars.shape[1],
    )
    return scalars, layout

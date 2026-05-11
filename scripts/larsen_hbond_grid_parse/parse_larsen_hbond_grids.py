#!/usr/bin/env python3
"""Parse Larsen 2015 ProCS15 H-bond DFT grids from
`/mnt/expansion/larsen_archive/hydrogenbondnmrlogs.tar`.

Streams through 6 nested archives (NMA|ALA donor × NMA|HOMe|acetate
acceptor). Per Gaussian log:
  - extracts atom positions from Standard orientation block;
  - identifies donor H + anchor + acceptor O + anchor by chemical
    perception (one-time per molecule, cached);
  - extracts per-atom GIAO shielding tensors;
  - computes (rOH, θ, ρ) from atom positions (not filename);
  - rotates donor-side readout tensors into a canonical donor frame
    so all grid points share a common tensor basis;
  - bins into a 3D grid keyed on actual (r, θ, ρ).

Per archive, emits:
  - `<DONOR><ACCEPTOR>_grid.npz` — packed tensors + axes + readout
    atom names + canonical donor frame anchor positions.
  - `<DONOR><ACCEPTOR>_meta.json` — human-readable metadata.

Reference subtraction (parser provides Δσ = σ_grid − σ_ref):
  - Reference σ per readout atom = average over the largest r-edge
    grid points (r ≥ r_max − 1 grid step). Proxy for free-monomer
    DFT (not in our ERDA archive). Documented limitation.

Usage:
  python3 parse_larsen_hbond_grids.py \\
      --tar /mnt/expansion/larsen_archive/hydrogenbondnmrlogs.tar \\
      --out data/larsen_hbond_grids/

Specific archive:
  python3 parse_larsen_hbond_grids.py --archive ALANMA --out ...
"""

from __future__ import annotations

import argparse
import bz2
import io
import json
import math
import re
import sys
import tarfile
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, Optional

import numpy as np


# ---------------------------------------------------------------------------
# Bond cutoffs (reused from scripts/perceive_larsen_tripeptide.py shape)
# ---------------------------------------------------------------------------

BOND_CUTOFFS_A: dict[frozenset, float] = {
    frozenset({"C", "C"}): 1.62,
    frozenset({"C", "N"}): 1.62,
    frozenset({"C", "O"}): 1.62,
    frozenset({"N", "N"}): 1.62,
    frozenset({"N", "O"}): 1.62,
    frozenset({"O", "O"}): 1.50,
    frozenset({"C", "H"}): 1.20,
    frozenset({"N", "H"}): 1.20,
    frozenset({"O", "H"}): 1.20,
}

ELEMENT_BY_Z = {1: "H", 6: "C", 7: "N", 8: "O", 16: "S"}


# ---------------------------------------------------------------------------
# Archive configuration — donor / acceptor kinds per archive name
# ---------------------------------------------------------------------------

# DonorKind in {"NMA", "ALA"}
# AcceptorKind in {"NMA", "HOMe", "Acetate"}
ARCHIVE_CONFIG = {
    "NMANMA": ("NMA", "NMA"),
    "NMACOH": ("NMA", "HOMe"),
    "NMACOO": ("NMA", "Acetate"),
    "ALANMA": ("ALA", "NMA"),
    "ALACOH": ("ALA", "HOMe"),
    "ALACOO": ("ALA", "Acetate"),
}


# Readout atom names per donor — these are the atoms Larsen 2015 Table 2
# specifies as receiving 1° (donor-side) and 2° (acceptor-side) contributions.
DONOR_READOUTS = {
    # ALA donor's residue i atoms (Δσ_1°HαB readouts)
    "ALA": ["N", "CA", "CB", "C", "HA", "HN"],
    # NMA donor's atoms (Δσ_1°HB readouts; mapped onto residue i)
    # The NMA donor has: amide H (HN), amide N, the two methyls' Cs
    # mapped onto Cα (the methyl on the N side, since N–CH3 ~ N–Cα) and
    # C' (the methyl on the C=O side, since C=O has the methyl on the
    # acyl side). HA is averaged over 3 methyl Hs.
    "NMA": ["N", "CA", "C", "HA", "HN"],
}

ACCEPTOR_READOUTS = {
    # NMA acceptor's atoms (Δσ_2° readouts — representing residue i+1)
    # The acceptor NMA's amide N+H are read as i+1's N+HN; one methyl's
    # 3 Hs averaged to represent i+1's HA.
    "NMA": ["N", "HN", "HA"],  # Hα = methyl-3H-average
    "HOMe": [],   # No 2° readout — no i+1 mapping for methanol
    "Acetate": [],  # No 2° readout — no i+1 mapping for acetate
}


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def vec(a, b):
    return [b[0] - a[0], b[1] - a[1], b[2] - a[2]]


def norm(v):
    return math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


def normalize(v):
    n = norm(v)
    return [v[0] / n, v[1] / n, v[2] / n]


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a, b):
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]


def distance(p, q):
    return norm(vec(p, q))


def angle_deg(a, b, c):
    """Angle at vertex b, formed by a-b-c."""
    ba = vec(b, a)
    bc = vec(b, c)
    cos_a = dot(ba, bc) / (norm(ba) * norm(bc))
    cos_a = max(-1.0, min(1.0, cos_a))
    return math.degrees(math.acos(cos_a))


def dihedral_deg(a, b, c, d):
    """Dihedral a-b-c-d (IUPAC convention)."""
    b1 = vec(a, b)
    b2 = vec(b, c)
    b3 = vec(c, d)
    n1 = cross(b1, b2)
    n2 = cross(b2, b3)
    m1 = cross(n1, normalize(b2))
    x = dot(n1, n2)
    y = dot(m1, n2)
    return math.degrees(math.atan2(y, x))


# ---------------------------------------------------------------------------
# Log parsing
# ---------------------------------------------------------------------------

_RE_ATOM_GEOM = re.compile(
    r"^\s*(\d+)\s+(\d+)\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$"
)

# "  1  N    Isotropic =   150.1855   Anisotropy =   121.5041"
_RE_TENSOR_HEADER = re.compile(
    r"^\s*(\d+)\s+([A-Z][a-z]?)\s+Isotropic\s*=\s*(-?\d+\.\d+)\s+Anisotropy\s*=\s*(-?\d+\.\d+)\s*$"
)

# "   XX=    90.2658   YX=   -79.9446   ZX=   -56.4689"
_RE_TENSOR_ROW1 = re.compile(
    r"^\s*XX=\s*(-?\d+\.\d+)\s+YX=\s*(-?\d+\.\d+)\s+ZX=\s*(-?\d+\.\d+)\s*$"
)
_RE_TENSOR_ROW2 = re.compile(
    r"^\s*XY=\s*(-?\d+\.\d+)\s+YY=\s*(-?\d+\.\d+)\s+ZY=\s*(-?\d+\.\d+)\s*$"
)
_RE_TENSOR_ROW3 = re.compile(
    r"^\s*XZ=\s*(-?\d+\.\d+)\s+YZ=\s*(-?\d+\.\d+)\s+ZZ=\s*(-?\d+\.\d+)\s*$"
)


@dataclass
class LogAtom:
    idx: int  # 1-based
    Z: int
    element: str
    pos: list[float]  # Å, in Standard orientation frame
    iso: float = 0.0
    anis: float = 0.0
    sigma: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))


def parse_log(text: str) -> list[LogAtom]:
    """Parse a Gaussian log: returns atoms with positions + σ tensors.

    Reads the FIRST Standard orientation block, then the GIAO shielding
    tensor block. Tensor block contains per-atom 3x3 σ in PPM in the
    same Cartesian frame as Standard orientation.
    """
    atoms_by_idx: dict[int, LogAtom] = {}

    # Pass 1: Standard orientation (positions).
    capture = False
    sep_count = 0
    for line in text.splitlines():
        if "Standard orientation" in line:
            capture = True
            sep_count = 0
            atoms_by_idx.clear()
            continue
        if capture:
            if re.match(r"^\s*-+\s*$", line):
                sep_count += 1
                if sep_count >= 3:
                    break
                continue
            m = _RE_ATOM_GEOM.match(line)
            if m:
                idx = int(m.group(1))
                Z = int(m.group(2))
                x, y, z = float(m.group(3)), float(m.group(4)), float(m.group(5))
                atoms_by_idx[idx] = LogAtom(
                    idx=idx, Z=Z, element=ELEMENT_BY_Z.get(Z, "X"),
                    pos=[x, y, z],
                )
    if not atoms_by_idx:
        raise ValueError("Standard orientation block not found")

    # Pass 2: SCF GIAO Magnetic shielding tensor block.
    state = "scan"  # scan -> after-header (filling 3 rows)
    pending: Optional[LogAtom] = None
    row = 0
    sigma = np.zeros((3, 3))
    in_tensor_block = False
    for line in text.splitlines():
        if "SCF GIAO Magnetic shielding tensor" in line:
            in_tensor_block = True
            continue
        if not in_tensor_block:
            continue
        if state == "scan":
            m = _RE_TENSOR_HEADER.match(line)
            if m:
                idx = int(m.group(1))
                if idx not in atoms_by_idx:
                    raise ValueError(f"tensor for unknown atom idx {idx}")
                pending = atoms_by_idx[idx]
                pending.iso = float(m.group(3))
                pending.anis = float(m.group(4))
                sigma = np.zeros((3, 3))
                row = 0
                state = "row1"
                continue
        elif state == "row1":
            m = _RE_TENSOR_ROW1.match(line)
            if m:
                sigma[0, 0] = float(m.group(1))  # XX
                sigma[1, 0] = float(m.group(2))  # YX
                sigma[2, 0] = float(m.group(3))  # ZX
                state = "row2"
                continue
        elif state == "row2":
            m = _RE_TENSOR_ROW2.match(line)
            if m:
                sigma[0, 1] = float(m.group(1))  # XY
                sigma[1, 1] = float(m.group(2))  # YY
                sigma[2, 1] = float(m.group(3))  # ZY
                state = "row3"
                continue
        elif state == "row3":
            m = _RE_TENSOR_ROW3.match(line)
            if m:
                sigma[0, 2] = float(m.group(1))  # XZ
                sigma[1, 2] = float(m.group(2))  # YZ
                sigma[2, 2] = float(m.group(3))  # ZZ
                if pending is not None:
                    pending.sigma = sigma.copy()
                pending = None
                state = "scan"
                continue

    return list(atoms_by_idx.values())


# ---------------------------------------------------------------------------
# Bond graph + chemical perception
# ---------------------------------------------------------------------------

def build_bond_graph(atoms: list[LogAtom]) -> list[set[int]]:
    """0-indexed adjacency list."""
    n = len(atoms)
    adj: list[set[int]] = [set() for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            cutoff = BOND_CUTOFFS_A.get(
                frozenset({atoms[i].element, atoms[j].element}), 0.0
            )
            if cutoff <= 0:
                continue
            if distance(atoms[i].pos, atoms[j].pos) < cutoff:
                adj[i].add(j)
                adj[j].add(i)
    return adj


def cluster_components(adj: list[set[int]], n: int) -> list[set[int]]:
    seen = set()
    clusters = []
    for i in range(n):
        if i in seen:
            continue
        stack = [i]
        comp = set()
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            comp.add(x)
            stack.extend(adj[x])
        clusters.append(comp)
    return clusters


@dataclass
class CanonicalSystem:
    """Canonical atom indices (0-based into the log's atom list) for a single
    donor-acceptor system. Identified ONCE per archive from a reference log,
    then reused across all grid points (since Larsen froze monomer geometries).
    """
    # Donor side
    donor_H: int
    donor_anchor: int          # Cα for ALA, amide N for NMA
    donor_third: int           # third atom for canonical-frame definition
    donor_readout: dict[str, int]  # readout name → atom index

    # Acceptor side
    acceptor_O: int
    acceptor_C: int            # carbonyl C (or hydroxyl C for HOMe)
    acceptor_third: int        # third atom for dihedral reference
    acceptor_readout: dict[str, int]  # readout name → atom index (empty if no 2°)


def perceive_canonical(
    atoms: list[LogAtom],
    donor_kind: str,
    acceptor_kind: str,
) -> CanonicalSystem:
    """Identify canonical atoms for the donor + acceptor system.

    Bond-graph perception. Works on a single log; result is shared across
    all grid points in the archive (monomer geometries frozen).
    """
    n = len(atoms)
    elem = [a.element for a in atoms]
    pos = [a.pos for a in atoms]
    adj = build_bond_graph(atoms)
    clusters = cluster_components(adj, n)
    if len(clusters) != 2:
        raise ValueError(f"expected 2 molecules, found {len(clusters)} clusters")

    # Find donor cluster.
    # - ALA donor (Ac-A-NMe = 6 C) is unambiguously larger than NMA acceptor (3 C).
    # - NMA donor + NMA acceptor: BOTH clusters have 3 C. Disambiguate by
    #   H-bond proximity: the donor cluster is the one whose amide N-H is
    #   close to an O in the other cluster.
    # - NMA donor + HOMe acceptor: NMA donor has 3 C, HOMe has 1 C. Donor is larger.
    # - NMA donor + Acetate acceptor: NMA has 3 C, acetate has 2 C. Donor is larger.
    if donor_kind == "ALA":
        clusters_sorted = sorted(clusters, key=len, reverse=True)
        donor_cluster = clusters_sorted[0]
        acceptor_cluster = clusters_sorted[1]
    elif donor_kind == "NMA":
        sizes = [len(c) for c in clusters]
        if sizes[0] != sizes[1]:
            # Asymmetric: donor is the larger cluster.
            clusters_sorted = sorted(clusters, key=len, reverse=True)
            donor_cluster = clusters_sorted[0]
            acceptor_cluster = clusters_sorted[1]
        else:
            # NMA-NMA symmetric: pick by H-bond proximity.
            # For each cluster, find its amide N-H. The donor is the cluster
            # whose amide H is closer to an O in the OTHER cluster.
            def find_amide_H(cluster):
                for n_idx in cluster:
                    if elem[n_idx] != "N":
                        continue
                    nbs = adj[n_idx]
                    n_h = sum(1 for x in nbs if elem[x] == "H")
                    n_c = sum(1 for x in nbs if elem[x] == "C")
                    if n_h == 1 and n_c == 2:
                        for x in nbs:
                            if elem[x] == "H":
                                return x
                return None
            c1, c2 = clusters
            h1 = find_amide_H(c1)
            h2 = find_amide_H(c2)
            o_in_c1 = [i for i in c1 if elem[i] == "O"]
            o_in_c2 = [i for i in c2 if elem[i] == "O"]
            d1 = min(distance(pos[h1], pos[o]) for o in o_in_c2) if h1 else float("inf")
            d2 = min(distance(pos[h2], pos[o]) for o in o_in_c1) if h2 else float("inf")
            if d1 < d2:
                donor_cluster, acceptor_cluster = c1, c2
            else:
                donor_cluster, acceptor_cluster = c2, c1

    # --- Donor side identification ---
    donor_H = None
    donor_anchor = None
    donor_third = None
    donor_readout: dict[str, int] = {}

    if donor_kind == "ALA":
        # Find Cα: C with neighbours (1 N, 2 C, 1 H).
        for c in donor_cluster:
            if elem[c] != "C":
                continue
            nbs = adj[c]
            n_h = sum(1 for x in nbs if elem[x] == "H")
            n_n = sum(1 for x in nbs if elem[x] == "N")
            n_c = sum(1 for x in nbs if elem[x] == "C")
            if n_h == 1 and n_n == 1 and n_c == 2:
                donor_anchor = c  # Cα
                donor_readout["CA"] = c
                for x in nbs:
                    if elem[x] == "H":
                        donor_H = x  # Hα
                        donor_readout["HA"] = x
                    elif elem[x] == "N":
                        donor_third = x  # alanine N
                        donor_readout["N"] = x
                # Identify Cβ (methyl C, 3 H + 1 C neighbours, the bonded C is Cα)
                # and C' (carbonyl C, has =O and -N neighbours)
                for x in nbs:
                    if elem[x] != "C":
                        continue
                    x_nbs = adj[x]
                    x_h = sum(1 for y in x_nbs if elem[y] == "H")
                    x_o = sum(1 for y in x_nbs if elem[y] == "O")
                    x_n = sum(1 for y in x_nbs if elem[y] == "N")
                    if x_h == 3:
                        donor_readout["CB"] = x
                    elif x_o == 1 and x_n == 1:
                        donor_readout["C"] = x
                break
        # Identify HN: H bonded to the alanine N (donor_third).
        if donor_third is not None:
            for x in adj[donor_third]:
                if elem[x] == "H":
                    donor_readout["HN"] = x
                    break
    elif donor_kind == "NMA":
        # Find amide N: N with neighbours (1 H, 2 C).
        for n_idx in donor_cluster:
            if elem[n_idx] != "N":
                continue
            nbs = adj[n_idx]
            n_h = sum(1 for x in nbs if elem[x] == "H")
            n_c = sum(1 for x in nbs if elem[x] == "C")
            if n_h == 1 and n_c == 2:
                donor_anchor = n_idx  # amide N
                donor_readout["N"] = n_idx
                for x in nbs:
                    if elem[x] == "H":
                        donor_H = x
                        donor_readout["HN"] = x
                # Among the 2 C neighbours: one is carbonyl C' (has =O),
                # the other is the methyl C (becomes Hα-stand-in via 3 methyl Hs).
                # Per the Δσ_1°HB convention: the methyl on the N side is
                # the "Cα" stand-in (N-CH3 ~ N-Cα).
                for x in nbs:
                    if elem[x] != "C":
                        continue
                    x_nbs = adj[x]
                    if any(elem[y] == "O" for y in x_nbs):
                        donor_readout["C"] = x  # carbonyl C
                        donor_third = x  # third frame atom: the carbonyl C
                    else:
                        donor_readout["CA"] = x  # methyl-as-Cα-stand-in
                # HA stand-in: average of 3 methyl Hs (handled at extraction time)
                if "CA" in donor_readout:
                    methyl_Hs = [y for y in adj[donor_readout["CA"]] if elem[y] == "H"]
                    donor_readout["HA"] = methyl_Hs[0]  # placeholder; we'll
                                                         # average at extraction
                    donor_readout["_HA_average"] = tuple(methyl_Hs)
                break

    if donor_H is None or donor_anchor is None or donor_third is None:
        raise ValueError(f"failed to identify donor anchors (donor_kind={donor_kind})")

    # --- Acceptor side identification ---
    acceptor_O = None
    acceptor_C = None
    acceptor_third = None
    acceptor_readout: dict[str, int] = {}

    if acceptor_kind == "NMA":
        # NMA acceptor: 1 O on a C=O; find amide N + H + the methyl on N-side
        # The acceptor O is the only O in this cluster.
        os_in = [i for i in acceptor_cluster if elem[i] == "O"]
        if len(os_in) != 1:
            raise ValueError(f"NMA acceptor expected 1 O, got {len(os_in)}")
        acceptor_O = os_in[0]
        # Acceptor carbonyl C: the C bonded to that O.
        for x in adj[acceptor_O]:
            if elem[x] == "C":
                acceptor_C = x
                break
        # Acceptor third: the N bonded to that C (for Hα..O=C-N dihedral).
        # Larsen's convention "Hα..O=C-N or Hα..O=C-C" — pick N preferentially.
        ns_on_C = [x for x in adj[acceptor_C] if elem[x] == "N"]
        cs_on_C = [x for x in adj[acceptor_C] if elem[x] == "C" and x != acceptor_O]
        if ns_on_C:
            acceptor_third = ns_on_C[0]
        elif cs_on_C:
            acceptor_third = cs_on_C[0]
        # 2° readouts (acceptor's amide group → represents i+1 backbone)
        amide_N = ns_on_C[0] if ns_on_C else None
        if amide_N is not None:
            acceptor_readout["N"] = amide_N
            for x in adj[amide_N]:
                if elem[x] == "H":
                    acceptor_readout["HN"] = x
                # The methyl on the N (representing i+1's Cα methyl, Hα stand-in)
                elif elem[x] == "C":
                    methyl_Hs = [y for y in adj[x] if elem[y] == "H"]
                    if len(methyl_Hs) == 3:
                        acceptor_readout["HA"] = methyl_Hs[0]
                        acceptor_readout["_HA_average"] = tuple(methyl_Hs)
    elif acceptor_kind == "HOMe":
        # Methanol acceptor: HO-CH3. Find the O with neighbours (1 H, 1 C).
        for x in acceptor_cluster:
            if elem[x] != "O":
                continue
            nbs = adj[x]
            n_h = sum(1 for y in nbs if elem[y] == "H")
            n_c = sum(1 for y in nbs if elem[y] == "C")
            if n_h == 1 and n_c == 1:
                acceptor_O = x
                for y in nbs:
                    if elem[y] == "C":
                        acceptor_C = y
                # Acceptor third: the H on the O (for Hα..O-H dihedral
                # per Larsen's "Hα..O-C(..)H^O" convention).
                for y in nbs:
                    if elem[y] == "H":
                        acceptor_third = y
                break
    elif acceptor_kind == "Acetate":
        # Acetate (CH3-COO⁻): 2 O atoms on a single C. Pick the O nearer
        # the donor H (closest inter-cluster H-O contact).
        os_in = [i for i in acceptor_cluster if elem[i] == "O"]
        if len(os_in) != 2:
            raise ValueError(f"Acetate acceptor expected 2 O, got {len(os_in)}")
        # Closest to donor H
        best = None
        for o in os_in:
            d = distance(pos[donor_H], pos[o])
            if best is None or d < best[0]:
                best = (d, o)
        acceptor_O = best[1]
        for x in adj[acceptor_O]:
            if elem[x] == "C":
                acceptor_C = x
                break
        # Acceptor third: the OTHER O on the C (defines C-O reference).
        # Larsen says "Hα..O=C-C or Hα..O=C-O"-like; the symmetric COO has
        # one C and one O on the carboxylate C; pick the OTHER O.
        for x in adj[acceptor_C]:
            if x != acceptor_O and elem[x] == "O":
                acceptor_third = x
                break
        if acceptor_third is None:
            # Fallback: pick the methyl C on the carboxylate C.
            for x in adj[acceptor_C]:
                if elem[x] == "C":
                    acceptor_third = x
                    break

    if acceptor_O is None or acceptor_C is None or acceptor_third is None:
        raise ValueError(f"failed to identify acceptor anchors (acceptor_kind={acceptor_kind})")

    return CanonicalSystem(
        donor_H=donor_H,
        donor_anchor=donor_anchor,
        donor_third=donor_third,
        donor_readout=donor_readout,
        acceptor_O=acceptor_O,
        acceptor_C=acceptor_C,
        acceptor_third=acceptor_third,
        acceptor_readout=acceptor_readout,
    )


# ---------------------------------------------------------------------------
# Canonical donor frame + tensor rotation
# ---------------------------------------------------------------------------

def canonical_donor_frame(
    donor_H_pos: list[float],
    donor_anchor_pos: list[float],
    donor_third_pos: list[float],
) -> np.ndarray:
    """Return 3x3 rotation matrix R such that R @ v_in_log_frame gives v
    in the canonical donor frame.

    Canonical frame:
      - Origin at donor H (translation removed; we only return R).
      - z-axis: donor_anchor → donor_H normalized.
      - x-axis: in the donor_third → donor_H direction projected
        orthogonal to z.
      - y-axis: z × x.
    """
    z = normalize(vec(donor_anchor_pos, donor_H_pos))
    # third → H direction
    v3 = vec(donor_third_pos, donor_H_pos)
    # Remove z component
    proj = dot(v3, z)
    x_raw = [v3[i] - proj * z[i] for i in range(3)]
    x = normalize(x_raw)
    y = cross(z, x)
    # Rotation matrix: rows are the canonical basis vectors expressed in log frame.
    # R @ v_log = v_canonical, where v_canonical's components are
    # (v_log · x, v_log · y, v_log · z) — so R has rows [x; y; z].
    return np.array([x, y, z])


def rotate_tensor(sigma_log: np.ndarray, R: np.ndarray) -> np.ndarray:
    """sigma_canonical = R @ sigma_log @ R.T"""
    return R @ sigma_log @ R.T


# ---------------------------------------------------------------------------
# Per-log processing
# ---------------------------------------------------------------------------

@dataclass
class GridPoint:
    r: float
    theta: float
    rho: float
    donor_sigma: dict[str, np.ndarray]      # readout name → σ in canonical frame
    acceptor_sigma: dict[str, np.ndarray]   # readout name → σ in canonical frame
    log_name: str = ""


def process_log(
    text: str,
    canon: CanonicalSystem,
    donor_readouts: list[str],
    acceptor_readouts: list[str],
) -> GridPoint:
    """Parse one log; return its grid point.

    Geometry is computed from the log's atom positions; tensors are rotated
    into the canonical donor frame so the grid has a common tensor basis.
    """
    atoms = parse_log(text)
    pos = [a.pos for a in atoms]

    # Verify canonical atom indices are present
    if canon.donor_H >= len(atoms) or canon.acceptor_O >= len(atoms):
        raise ValueError("atom index out of range")

    r = distance(pos[canon.donor_H], pos[canon.acceptor_O])
    theta = angle_deg(pos[canon.donor_H], pos[canon.acceptor_O], pos[canon.acceptor_C])
    rho = dihedral_deg(
        pos[canon.donor_H],
        pos[canon.acceptor_O],
        pos[canon.acceptor_C],
        pos[canon.acceptor_third],
    )

    # Build canonical donor frame rotation
    R = canonical_donor_frame(
        pos[canon.donor_H],
        pos[canon.donor_anchor],
        pos[canon.donor_third],
    )

    # Extract donor-side readouts in canonical frame
    donor_sigma: dict[str, np.ndarray] = {}
    for name in donor_readouts:
        if name not in canon.donor_readout:
            continue  # not all donors have all readouts
        i = canon.donor_readout[name]
        sigma = atoms[i].sigma
        # Special case: HA average for NMA donor (3 methyl Hs)
        if name == "HA" and "_HA_average" in canon.donor_readout:
            hs = canon.donor_readout["_HA_average"]
            avg = np.zeros((3, 3))
            for h_i in hs:
                avg += atoms[h_i].sigma
            sigma = avg / len(hs)
        donor_sigma[name] = rotate_tensor(sigma, R)

    # Extract acceptor-side readouts in canonical frame
    acceptor_sigma: dict[str, np.ndarray] = {}
    for name in acceptor_readouts:
        if name not in canon.acceptor_readout:
            continue
        i = canon.acceptor_readout[name]
        sigma = atoms[i].sigma
        if name == "HA" and "_HA_average" in canon.acceptor_readout:
            hs = canon.acceptor_readout["_HA_average"]
            avg = np.zeros((3, 3))
            for h_i in hs:
                avg += atoms[h_i].sigma
            sigma = avg / len(hs)
        acceptor_sigma[name] = rotate_tensor(sigma, R)

    return GridPoint(
        r=r, theta=theta, rho=rho,
        donor_sigma=donor_sigma,
        acceptor_sigma=acceptor_sigma,
    )


# ---------------------------------------------------------------------------
# Archive streaming
# ---------------------------------------------------------------------------

def iter_archive_logs(
    outer_tar_path: Path,
    archive_stem: str,
) -> Iterator[tuple[str, str]]:
    """Yield (log_filename, log_text) from the inner archive.

    `archive_stem` is one of ALANMA, ALACOO, ALACOH, NMANMA, NMACOH, NMACOO.
    """
    inner_path = f"hydrogenbonds/{archive_stem}nmrlog.tar.bz2"
    with tarfile.open(outer_tar_path, "r") as outer:
        inner_member = outer.getmember(inner_path)
        inner_f = outer.extractfile(inner_member)
        if inner_f is None:
            raise RuntimeError(f"could not extract {inner_path}")
        # The inner file is bz2-compressed tar
        inner_data = bz2.decompress(inner_f.read())
        inner_buf = io.BytesIO(inner_data)
        with tarfile.open(fileobj=inner_buf, mode="r:") as inner_tar:
            for member in inner_tar:
                if not member.isfile():
                    continue
                if not member.name.endswith(".log"):
                    continue
                f = inner_tar.extractfile(member)
                if f is None:
                    continue
                text = f.read().decode("utf-8", errors="replace")
                yield Path(member.name).name, text


# ---------------------------------------------------------------------------
# Grid assembly + emission
# ---------------------------------------------------------------------------

@dataclass
class ArchiveGrid:
    archive_stem: str
    donor_kind: str
    acceptor_kind: str
    canon: CanonicalSystem
    points: list[GridPoint] = field(default_factory=list)
    donor_readout_names: list[str] = field(default_factory=list)
    acceptor_readout_names: list[str] = field(default_factory=list)


def build_archive_grid(
    tar_path: Path,
    archive_stem: str,
    max_logs: Optional[int] = None,
) -> ArchiveGrid:
    donor_kind, acceptor_kind = ARCHIVE_CONFIG[archive_stem]
    donor_readouts = DONOR_READOUTS[donor_kind]
    acceptor_readouts = ACCEPTOR_READOUTS[acceptor_kind]

    canon: Optional[CanonicalSystem] = None
    points: list[GridPoint] = []
    seen_keys: set[tuple[float, float, float]] = set()
    n_processed = 0
    n_skipped_duplicate = 0
    n_failed_dft = 0  # logs with no Standard orientation block (Larsen-skipped)

    for log_name, log_text in iter_archive_logs(tar_path, archive_stem):
        if max_logs is not None and n_processed >= max_logs:
            break
        # Skip Gaussian-failed calculations (Larsen 2015 notes these were
        # interpolated in his pipeline; we drop them and let runtime trilinear
        # interpolation fill the gaps).
        if "Standard orientation" not in log_text:
            n_failed_dft += 1
            continue
        try:
            if canon is None:
                # One-time perception on first successful log
                atoms = parse_log(log_text)
                canon = perceive_canonical(atoms, donor_kind, acceptor_kind)
                print(f"  [{archive_stem}] canonical anchors identified: "
                      f"donor_H={canon.donor_H+1}, donor_anchor={canon.donor_anchor+1}, "
                      f"acceptor_O={canon.acceptor_O+1}, acceptor_C={canon.acceptor_C+1}")

            gp = process_log(log_text, canon, donor_readouts, acceptor_readouts)
            gp.log_name = log_name
            # Round to grid resolution to detect duplicates
            key = (round(gp.r, 3), round(gp.theta, 2), round(gp.rho, 2))
            if key in seen_keys:
                n_skipped_duplicate += 1
                continue
            seen_keys.add(key)
            points.append(gp)
            n_processed += 1
            if n_processed % 500 == 0:
                print(f"  [{archive_stem}] processed {n_processed} logs...")
        except Exception as e:
            print(f"  [{archive_stem}] FAIL on {log_name}: {e}", file=sys.stderr)
            raise

    print(f"  [{archive_stem}] DONE: {n_processed} grid points, "
          f"{n_skipped_duplicate} duplicates skipped, "
          f"{n_failed_dft} Larsen-failed DFT logs skipped")

    return ArchiveGrid(
        archive_stem=archive_stem,
        donor_kind=donor_kind,
        acceptor_kind=acceptor_kind,
        canon=canon,
        points=points,
        donor_readout_names=donor_readouts,
        acceptor_readout_names=acceptor_readouts,
    )


def compute_reference_and_subtract(grid: ArchiveGrid) -> None:
    """Reference σ per readout = average over r-max grid edge.

    For each readout atom (donor + acceptor): compute mean σ across grid
    points where r is at or above (max_r - one_grid_step). Subtract this
    reference from every grid point's σ → grid stores Δσ.
    """
    if not grid.points:
        return
    max_r = max(p.r for p in grid.points)
    rs = sorted({round(p.r, 3) for p in grid.points})
    grid_step_r = rs[1] - rs[0] if len(rs) >= 2 else 0.2
    r_threshold = max_r - grid_step_r * 0.5  # only the largest-r layer
    edge_points = [p for p in grid.points if p.r >= r_threshold]
    if not edge_points:
        raise ValueError(f"no points at r-max edge for {grid.archive_stem}")

    # Reference per readout atom
    donor_ref: dict[str, np.ndarray] = {}
    for name in grid.donor_readout_names:
        if name not in edge_points[0].donor_sigma:
            continue
        avg = np.mean([p.donor_sigma[name] for p in edge_points if name in p.donor_sigma], axis=0)
        donor_ref[name] = avg
    acceptor_ref: dict[str, np.ndarray] = {}
    for name in grid.acceptor_readout_names:
        if not edge_points[0].acceptor_sigma:
            continue
        if name not in edge_points[0].acceptor_sigma:
            continue
        avg = np.mean([p.acceptor_sigma[name] for p in edge_points if name in p.acceptor_sigma], axis=0)
        acceptor_ref[name] = avg

    # Subtract from every point
    for p in grid.points:
        for name, ref in donor_ref.items():
            if name in p.donor_sigma:
                p.donor_sigma[name] = p.donor_sigma[name] - ref
        for name, ref in acceptor_ref.items():
            if name in p.acceptor_sigma:
                p.acceptor_sigma[name] = p.acceptor_sigma[name] - ref

    print(f"  [{grid.archive_stem}] reference subtracted from "
          f"{len(edge_points)} r-max edge points (r ≥ {r_threshold:.3f} Å)")


def emit_grid(grid: ArchiveGrid, out_dir: Path) -> None:
    """Write <stem>_grid.npz and <stem>_meta.json."""
    out_dir.mkdir(parents=True, exist_ok=True)

    n_points = len(grid.points)
    rs = np.array([p.r for p in grid.points])
    thetas = np.array([p.theta for p in grid.points])
    rhos = np.array([p.rho for p in grid.points])

    save_arrays: dict[str, np.ndarray] = {
        "r": rs.astype(np.float64),
        "theta": thetas.astype(np.float64),
        "rho": rhos.astype(np.float64),
    }

    for name in grid.donor_readout_names:
        if not any(name in p.donor_sigma for p in grid.points):
            continue
        arr = np.array([p.donor_sigma.get(name, np.full((3, 3), np.nan))
                        for p in grid.points])
        save_arrays[f"donor_{name}"] = arr.astype(np.float64)

    for name in grid.acceptor_readout_names:
        if not any(name in p.acceptor_sigma for p in grid.points):
            continue
        arr = np.array([p.acceptor_sigma.get(name, np.full((3, 3), np.nan))
                        for p in grid.points])
        save_arrays[f"acceptor_{name}"] = arr.astype(np.float64)

    npz_path = out_dir / f"{grid.archive_stem}_grid.npz"
    np.savez_compressed(npz_path, **save_arrays)
    print(f"  [{grid.archive_stem}] wrote {npz_path} ({n_points} points)")

    meta = {
        "archive_stem": grid.archive_stem,
        "donor_kind": grid.donor_kind,
        "acceptor_kind": grid.acceptor_kind,
        "n_points": n_points,
        "axes": {
            "r_min": float(rs.min()), "r_max": float(rs.max()),
            "theta_min": float(thetas.min()), "theta_max": float(thetas.max()),
            "rho_min": float(rhos.min()), "rho_max": float(rhos.max()),
        },
        "donor_readouts": [k.replace("donor_", "")
                           for k in save_arrays if k.startswith("donor_")],
        "acceptor_readouts": [k.replace("acceptor_", "")
                              for k in save_arrays if k.startswith("acceptor_")],
        "canonical_indices_1based": {
            "donor_H": grid.canon.donor_H + 1,
            "donor_anchor": grid.canon.donor_anchor + 1,
            "donor_third": grid.canon.donor_third + 1,
            "acceptor_O": grid.canon.acceptor_O + 1,
            "acceptor_C": grid.canon.acceptor_C + 1,
            "acceptor_third": grid.canon.acceptor_third + 1,
            "donor_readout_atoms": {k: v + 1 for k, v in grid.canon.donor_readout.items()
                                    if isinstance(v, int)},
            "acceptor_readout_atoms": {k: v + 1 for k, v in grid.canon.acceptor_readout.items()
                                       if isinstance(v, int)},
        },
        "reference_subtraction_method": "r-max-edge average (proxy for free-monomer DFT)",
    }
    meta_path = out_dir / f"{grid.archive_stem}_meta.json"
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)
    print(f"  [{grid.archive_stem}] wrote {meta_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--tar",
        default="/mnt/expansion/larsen_archive/hydrogenbondnmrlogs.tar",
        type=Path,
    )
    ap.add_argument(
        "--out",
        default=Path("data/larsen_hbond_grids"),
        type=Path,
    )
    ap.add_argument(
        "--archive",
        default=None,
        help="Specific archive stem (ALANMA, ALACOO, ALACOH, NMANMA, NMACOH, NMACOO). "
             "If omitted, processes all 6.",
    )
    ap.add_argument(
        "--max-logs",
        type=int,
        default=None,
        help="Limit per archive (for smoke testing).",
    )
    args = ap.parse_args()

    archives_to_process = (
        [args.archive] if args.archive else list(ARCHIVE_CONFIG.keys())
    )
    for stem in archives_to_process:
        if stem not in ARCHIVE_CONFIG:
            print(f"unknown archive: {stem}", file=sys.stderr)
            continue
        print(f"\n=== {stem} ===")
        grid = build_archive_grid(args.tar, stem, max_logs=args.max_logs)
        compute_reference_and_subtract(grid)
        emit_grid(grid, args.out)


if __name__ == "__main__":
    main()

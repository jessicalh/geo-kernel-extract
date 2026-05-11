#!/usr/bin/env python3
"""Bond-graph perception of a Larsen ProCS15 tripeptide DFT record.

This is the Python proof-of-concept that validates the perception
algorithm before the C++ port to `src/LarsenResidue.{h,cpp}`. It
reads `data/topology/tripeptide_orderings_dump.json` (one row per
`(tripeptide, frame_type)`) and for each row:

  1. Builds the bond graph from per-atom (position, element) using
     element-pair-specific distance cutoffs.
  2. Identifies the 4 peptide amides (C=O···N motif).
  3. Segments the molecule into the 5 pieces (ACE, NCapAla, Central,
     CCapAla, NME) by cutting the amides and finding connected
     components.
  4. Orders the pieces by walking the amide chain ACE → NME.
  5. Per piece, matches each atom against canonical AminoAcidType
     chemistry by element + bond-neighbour multiset (graph signature)
     iteration.
  6. Reports the per-atom typed identity (canonical name).

Failure modes are reported per row; no row should fail in a clean DB.
The algorithm here is the spec the C++ port implements.

Cross-check against the original Larsen Gaussian log for AAA is
available via --source-log.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
# Bond perception
# ---------------------------------------------------------------------------

BOND_CUTOFFS_A: dict[frozenset[str], float] = {
    # Heavy-heavy: covers single / double / aromatic / partial.
    frozenset({"C", "C"}): 1.85,
    frozenset({"C", "N"}): 1.85,
    frozenset({"C", "O"}): 1.85,
    frozenset({"N", "N"}): 1.85,
    frozenset({"N", "O"}): 1.85,
    frozenset({"O", "O"}): 1.60,  # only peroxide hits; carboxylate O atoms
                                   # are ~2.3 Å apart, won't bond
    # Heavy-S.
    frozenset({"C", "S"}): 2.10,
    frozenset({"N", "S"}): 2.10,
    frozenset({"O", "S"}): 2.10,
    frozenset({"S", "S"}): 2.30,  # disulfide
    # X-H.
    frozenset({"C", "H"}): 1.30,
    frozenset({"N", "H"}): 1.30,
    frozenset({"O", "H"}): 1.30,
    frozenset({"S", "H"}): 1.50,
}


def vlen(a: list[float], b: list[float]) -> float:
    dx, dy, dz = a[0] - b[0], a[1] - b[1], a[2] - b[2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def build_bond_graph(atoms: list[dict]) -> list[set[int]]:
    """Return adjacency list (per-atom set of 0-indexed neighbours)."""
    n = len(atoms)
    adj: list[set[int]] = [set() for _ in range(n)]
    pos = [[a["x"], a["y"], a["z"]] for a in atoms]
    elem = [a["element"] for a in atoms]
    for i in range(n):
        for j in range(i + 1, n):
            cutoff = BOND_CUTOFFS_A.get(frozenset({elem[i], elem[j]}), 0.0)
            if cutoff <= 0:
                continue
            if vlen(pos[i], pos[j]) < cutoff:
                adj[i].add(j); adj[j].add(i)
    return adj


# ---------------------------------------------------------------------------
# Amide identification + piece segmentation
# ---------------------------------------------------------------------------

def find_amide_bonds(atoms: list[dict], adj: list[set[int]]) -> list[tuple[int, int]]:
    """Return list of (carbonyl_C_idx, amide_N_idx) pairs. 0-indexed.

    Conservative C=O cutoff: 1.30 Å is too tight for some optimised
    geometries (planar amide C=O can sit at 1.20-1.24 Å but with
    sp2 contributions widening to ~1.28). Widen to 1.32 Å for the
    detection, then re-filter by sp2 geometry (the carbonyl C must
    have exactly 3 neighbours).
    """
    elem = [a["element"] for a in atoms]
    pos = [[a["x"], a["y"], a["z"]] for a in atoms]
    out: list[tuple[int, int]] = []
    for i in range(len(atoms)):
        if elem[i] != "C":
            continue
        if len(adj[i]) != 3:
            continue  # sp2 carbonyl C has 3 neighbours (=O, -N or -OH, -X)
        o_double = any(elem[j] == "O" and vlen(pos[i], pos[j]) < 1.32
                       for j in adj[i])
        if not o_double:
            continue
        for j in adj[i]:
            if elem[j] != "N":
                continue
            d = vlen(pos[i], pos[j])
            if 1.20 < d < 1.50:  # amide N-C partial-double widened
                out.append((i, j))
    return out


def connected_components(adj: list[set[int]], n: int) -> list[list[int]]:
    visited = [False] * n
    out: list[list[int]] = []
    for start in range(n):
        if visited[start]:
            continue
        stack = [start]; comp: list[int] = []
        while stack:
            cur = stack.pop()
            if visited[cur]:
                continue
            visited[cur] = True; comp.append(cur)
            stack.extend(adj[cur])
        out.append(sorted(comp))
    return out


def order_pieces_along_chain(
    pieces: list[list[int]],
    amides: list[tuple[int, int]],
    adj: list[set[int]],
) -> list[list[int]]:
    """Return pieces in linear ACE → NCapAla → Central → CCapAla → NME order
    by walking the directed amide chain (C-piece → N-piece)."""
    piece_of = [-1] * sum(len(p) for p in pieces)
    for pi, comp in enumerate(pieces):
        for ai in comp:
            piece_of[ai] = pi
    # Directed edges: piece(C) -> piece(N) per amide.
    out_edges: dict[int, int] = {}
    in_edges: dict[int, int] = {}
    for c, nj in amides:
        out_edges[piece_of[c]] = piece_of[nj]
        in_edges[piece_of[nj]] = piece_of[c]
    # ACE is the piece with an outgoing amide but no incoming amide.
    start = None
    for pi in range(len(pieces)):
        if pi in out_edges and pi not in in_edges:
            start = pi; break
    if start is None:
        raise RuntimeError("no chain start found (cyclic or disconnected)")
    chain = [start]
    cur = start
    while cur in out_edges:
        nxt = out_edges[cur]; chain.append(nxt); cur = nxt
    return [pieces[pi] for pi in chain]


# ---------------------------------------------------------------------------
# Canonical chemistry: copied minimally from src/AminoAcidType.cpp.
# ---------------------------------------------------------------------------

BB = [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"), ("H", "H"), ("HA", "H")]
BB_PRO = [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"), ("HA", "H")]
BB_GLY = [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"), ("H", "H"),
          ("HA2", "H"), ("HA3", "H")]

# Per-residue canonical atom list (name, element). Bonds derived from
# chi-angle defs + extras below.
CANONICAL: dict[str, list[tuple[str, str]]] = {
    "ALA": BB + [("CB","C"),("HB1","H"),("HB2","H"),("HB3","H")],
    "ARG": BB + [("CB","C"),("HB2","H"),("HB3","H"),("CG","C"),("HG2","H"),
                  ("HG3","H"),("CD","C"),("HD2","H"),("HD3","H"),("NE","N"),
                  ("HE","H"),("CZ","C"),("NH1","N"),("HH11","H"),("HH12","H"),
                  ("NH2","N"),("HH21","H"),("HH22","H")],
    "ASN": BB + [("CB","C"),("HB2","H"),("HB3","H"),("CG","C"),("OD1","O"),
                  ("ND2","N"),("HD21","H"),("HD22","H")],
    "ASP": BB + [("CB","C"),("HB2","H"),("HB3","H"),("CG","C"),("OD1","O"),
                  ("OD2","O")],
    "CYS": BB + [("CB","C"),("HB2","H"),("HB3","H"),("SG","S"),("HG","H")],
    "GLN": BB + [("CB","C"),("HB2","H"),("HB3","H"),("CG","C"),("HG2","H"),
                  ("HG3","H"),("CD","C"),("OE1","O"),("NE2","N"),
                  ("HE21","H"),("HE22","H")],
    "GLU": BB + [("CB","C"),("HB2","H"),("HB3","H"),("CG","C"),("HG2","H"),
                  ("HG3","H"),("CD","C"),("OE1","O"),("OE2","O")],
    "GLY": BB_GLY,
    # HIS canonical here is the HIE tautomer (H on N-epsilon2, none on
    # N-delta1). HID and HIP variants are handled by `try_his_variants`
    # — see PerceiveHisVariants below.
    "HIS": BB + [("CB","C"),("HB2","H"),("HB3","H"),("CG","C"),("ND1","N"),
                  ("CD2","C"),("HD2","H"),("CE1","C"),("HE1","H"),
                  ("NE2","N"),("HE2","H")],
    "ILE": BB + [("CB","C"),("HB","H"),("CG1","C"),("HG12","H"),("HG13","H"),
                  ("CG2","C"),("HG21","H"),("HG22","H"),("HG23","H"),
                  ("CD1","C"),("HD11","H"),("HD12","H"),("HD13","H")],
    "LEU": BB + [("CB","C"),("HB2","H"),("HB3","H"),("CG","C"),("HG","H"),
                  ("CD1","C"),("HD11","H"),("HD12","H"),("HD13","H"),
                  ("CD2","C"),("HD21","H"),("HD22","H"),("HD23","H")],
    "LYS": BB + [("CB","C"),("HB2","H"),("HB3","H"),("CG","C"),("HG2","H"),
                  ("HG3","H"),("CD","C"),("HD2","H"),("HD3","H"),
                  ("CE","C"),("HE2","H"),("HE3","H"),("NZ","N"),
                  ("HZ1","H"),("HZ2","H"),("HZ3","H")],
    "MET": BB + [("CB","C"),("HB2","H"),("HB3","H"),("CG","C"),("HG2","H"),
                  ("HG3","H"),("SD","S"),("CE","C"),("HE1","H"),("HE2","H"),
                  ("HE3","H")],
    "PHE": BB + [("CB","C"),("HB2","H"),("HB3","H"),("CG","C"),("CD1","C"),
                  ("HD1","H"),("CD2","C"),("HD2","H"),("CE1","C"),("HE1","H"),
                  ("CE2","C"),("HE2","H"),("CZ","C"),("HZ","H")],
    "PRO": BB_PRO + [("CB","C"),("HB2","H"),("HB3","H"),("CG","C"),("HG2","H"),
                      ("HG3","H"),("CD","C"),("HD2","H"),("HD3","H")],
    "SER": BB + [("CB","C"),("HB2","H"),("HB3","H"),("OG","O"),("HG","H")],
    "THR": BB + [("CB","C"),("HB","H"),("OG1","O"),("HG1","H"),
                  ("CG2","C"),("HG21","H"),("HG22","H"),("HG23","H")],
    "TRP": BB + [("CB","C"),("HB2","H"),("HB3","H"),("CG","C"),("CD1","C"),
                  ("HD1","H"),("CD2","C"),("NE1","N"),("HE1","H"),
                  ("CE2","C"),("CE3","C"),("HE3","H"),("CZ2","C"),
                  ("HZ2","H"),("CZ3","C"),("HZ3","H"),("CH2","C"),("HH2","H")],
    "TYR": BB + [("CB","C"),("HB2","H"),("HB3","H"),("CG","C"),("CD1","C"),
                  ("HD1","H"),("CD2","C"),("HD2","H"),("CE1","C"),("HE1","H"),
                  ("CE2","C"),("HE2","H"),("CZ","C"),("OH","O"),("HH","H")],
    "VAL": BB + [("CB","C"),("HB","H"),("CG1","C"),("HG11","H"),("HG12","H"),
                  ("HG13","H"),("CG2","C"),("HG21","H"),("HG22","H"),
                  ("HG23","H")],
}

HIS_VARIANTS: dict[str, list[tuple[str, str]]] = {
    "HID": BB + [("CB","C"),("HB2","H"),("HB3","H"),("CG","C"),("ND1","N"),
                  ("HD1","H"),("CD2","C"),("HD2","H"),("CE1","C"),("HE1","H"),
                  ("NE2","N")],
    "HIE": BB + [("CB","C"),("HB2","H"),("HB3","H"),("CG","C"),("ND1","N"),
                  ("CD2","C"),("HD2","H"),("CE1","C"),("HE1","H"),
                  ("NE2","N"),("HE2","H")],
    "HIP": BB + [("CB","C"),("HB2","H"),("HB3","H"),("CG","C"),("ND1","N"),
                  ("HD1","H"),("CD2","C"),("HD2","H"),("CE1","C"),("HE1","H"),
                  ("NE2","N"),("HE2","H")],
}


CHI: dict[str, list[tuple[str, ...]]] = {
    "ARG": [("N","CA","CB","CG"),("CA","CB","CG","CD"),("CB","CG","CD","NE"),
            ("CG","CD","NE","CZ")],
    "ASN": [("N","CA","CB","CG"),("CA","CB","CG","OD1")],
    "ASP": [("N","CA","CB","CG"),("CA","CB","CG","OD1")],
    "CYS": [("N","CA","CB","SG")],
    "GLN": [("N","CA","CB","CG"),("CA","CB","CG","CD"),("CB","CG","CD","OE1")],
    "GLU": [("N","CA","CB","CG"),("CA","CB","CG","CD"),("CB","CG","CD","OE1")],
    "HIS": [("N","CA","CB","CG"),("CA","CB","CG","ND1")],
    "ILE": [("N","CA","CB","CG1"),("CA","CB","CG1","CD1")],
    "LEU": [("N","CA","CB","CG"),("CA","CB","CG","CD1")],
    "LYS": [("N","CA","CB","CG"),("CA","CB","CG","CD"),("CB","CG","CD","CE"),
            ("CG","CD","CE","NZ")],
    "MET": [("N","CA","CB","CG"),("CA","CB","CG","SD"),("CB","CG","SD","CE")],
    "PHE": [("N","CA","CB","CG"),("CA","CB","CG","CD1")],
    "PRO": [("N","CA","CB","CG"),("CA","CB","CG","CD")],
    "SER": [("N","CA","CB","OG")],
    "THR": [("N","CA","CB","OG1")],
    "TRP": [("N","CA","CB","CG"),("CA","CB","CG","CD1")],
    "TYR": [("N","CA","CB","CG"),("CA","CB","CG","CD1")],
    "VAL": [("N","CA","CB","CG1")],
}

EXTRA_BONDS: dict[str, list[tuple[str, str]]] = {
    "ARG": [("CZ","NH1"),("CZ","NH2")],
    "ASN": [("CG","ND2")],
    "ASP": [("CG","OD2")],
    "GLN": [("CD","NE2")],
    "GLU": [("CD","OE2")],
    "HIS": [("CG","CD2"),("CD2","NE2"),("NE2","CE1"),("CE1","ND1")],
    "ILE": [("CB","CG2")],
    "LEU": [("CG","CD2")],
    "PHE": [("CG","CD2"),("CD1","CE1"),("CD2","CE2"),("CE1","CZ"),
            ("CE2","CZ")],
    "PRO": [("CD","N")],
    "THR": [("CB","CG2")],
    "TRP": [("CG","CD2"),("CD1","NE1"),("NE1","CE2"),("CD2","CE2"),
            ("CD2","CE3"),("CE2","CZ2"),("CE3","CZ3"),("CZ2","CH2"),
            ("CZ3","CH2")],
    "TYR": [("CG","CD2"),("CD1","CE1"),("CD2","CE2"),("CE1","CZ"),
            ("CE2","CZ"),("CZ","OH")],
    "VAL": [("CB","CG2")],
}

# ACE: -C(=O)-CH3. Atom names from AMBER ACE residue template:
#   C (carbonyl), O (carbonyl), CH3 (methyl), HH31, HH32, HH33 (3 H on methyl).
ACE_ATOMS = [("C", "C"), ("O", "O"), ("CH3", "C"),
             ("HH31", "H"), ("HH32", "H"), ("HH33", "H")]
ACE_BONDS = [("C", "O"), ("C", "CH3"), ("CH3", "HH31"),
             ("CH3", "HH32"), ("CH3", "HH33")]

# NME: -NH-CH3.
NME_ATOMS = [("N", "N"), ("H", "H"), ("CH3", "C"),
             ("HH31", "H"), ("HH32", "H"), ("HH33", "H")]
NME_BONDS = [("N", "H"), ("N", "CH3"), ("CH3", "HH31"),
             ("CH3", "HH32"), ("CH3", "HH33")]


def canonical_bonds_of_residue(residue: str) -> list[tuple[str, str]]:
    """Build canonical bond list for a residue from BB + chi + extras."""
    atoms = {n for n, _ in CANONICAL[residue]}
    bonds: list[tuple[str, str]] = []

    def add(a: str, b: str) -> None:
        if a in atoms and b in atoms:
            bonds.append(tuple(sorted([a, b])))

    # Backbone.
    add("N","CA"); add("CA","C"); add("C","O")
    if "H"  in atoms: add("N", "H")
    if "HA" in atoms: add("CA","HA")
    for h in ("HA2","HA3"):
        if h in atoms: add("CA", h)
    if "CB" in atoms: add("CA","CB")

    # Chi-chain bonds.
    for chi in CHI.get(residue, []):
        for i in range(len(chi) - 1):
            add(chi[i], chi[i+1])

    # Extras.
    for a, b in EXTRA_BONDS.get(residue, []):
        add(a, b)

    # Hydrogen attachments by locant convention.
    for atom_name, elem in CANONICAL[residue]:
        if elem != "H": continue
        if atom_name in {"H","HA","HA2","HA3","HE","HG","HZ","HH","HG1"}:
            continue  # already attached above or handled below
        # HBx → CB
        if atom_name.startswith("HB"): add(atom_name, "CB"); continue
        # HG1x → CG1; HG2x → CG2; HG3x → CG3 (none in standard 20)
        if atom_name.startswith("HG1") and "CG1" in atoms:
            add(atom_name, "CG1"); continue
        if atom_name.startswith("HG2") and "CG2" in atoms:
            add(atom_name, "CG2"); continue
        if atom_name.startswith("HG2") and "CG" in atoms:
            add(atom_name, "CG"); continue
        if atom_name.startswith("HG3"):
            add(atom_name, "CG" if "CG" in atoms else "CG3"); continue
        if atom_name.startswith("HD1") and "CD1" in atoms:
            add(atom_name, "CD1"); continue
        if atom_name.startswith("HD2") and "CD2" in atoms:
            add(atom_name, "CD2"); continue
        if atom_name.startswith("HD2") and "ND2" in atoms:
            add(atom_name, "ND2"); continue           # ASN HD21/HD22
        if atom_name.startswith("HD") and "CD" in atoms:
            add(atom_name, "CD"); continue
        if atom_name.startswith("HE1") and "CE1" in atoms:
            add(atom_name, "CE1"); continue
        if atom_name.startswith("HE1") and "NE1" in atoms:
            add(atom_name, "NE1"); continue           # TRP HE1
        if atom_name.startswith("HE2") and "CE2" in atoms:
            add(atom_name, "CE2"); continue
        if atom_name.startswith("HE2") and "NE2" in atoms:
            add(atom_name, "NE2"); continue           # GLN HE21/HE22, HIE HE2
        if atom_name.startswith("HE3") and "CE3" in atoms:
            add(atom_name, "CE3"); continue
        if atom_name.startswith("HE") and "CE" in atoms:
            add(atom_name, "CE"); continue
        if atom_name.startswith("HZ1") and "NZ" in atoms:
            add(atom_name, "NZ"); continue
        if atom_name.startswith("HZ2") and "NZ" in atoms:
            add(atom_name, "NZ"); continue
        if atom_name.startswith("HZ3") and "NZ" in atoms:
            add(atom_name, "NZ"); continue
        if atom_name.startswith("HZ2") and "CZ2" in atoms:
            add(atom_name, "CZ2"); continue
        if atom_name.startswith("HZ3") and "CZ3" in atoms:
            add(atom_name, "CZ3"); continue
        if atom_name == "HZ": add(atom_name, "CZ"); continue
        if atom_name.startswith("HH1"): add(atom_name, "NH1"); continue
        if atom_name.startswith("HH2") and "NH2" in atoms:
            add(atom_name, "NH2"); continue
        if atom_name.startswith("HH2") and "CH2" in atoms:
            add(atom_name, "CH2"); continue
        if atom_name == "HH": add(atom_name, "OH"); continue

    # Polar Hs not caught above (because they were filtered out at the
    # top of the loop to avoid mis-attaching them to a C/N parent).
    if "HE" in atoms and "NE" in atoms: add("HE", "NE")
    if "HG" in atoms:
        if "OG" in atoms: add("HG", "OG")
        elif "SG" in atoms: add("HG", "SG")
        elif "CG" in atoms: add("HG", "CG")  # LEU HG
    if "HG1" in atoms and "OG1" in atoms: add("HG1", "OG1")
    if "HZ" in atoms and "CZ" in atoms: add("HZ", "CZ")   # PHE
    if "HH" in atoms and "OH" in atoms: add("HH", "OH")   # TYR

    # De-duplicate.
    return sorted(set(bonds))


def canonical_adj(residue: str) -> dict[str, set[str]]:
    """Adjacency lookup for canonical residue atoms."""
    out: dict[str, set[str]] = {n: set() for n, _ in CANONICAL[residue]}
    for a, b in canonical_bonds_of_residue(residue):
        out[a].add(b); out[b].add(a)
    return out


def canonical_adj_for_cap(cap_atoms: list[tuple[str, str]],
                          cap_bonds: list[tuple[str, str]]
                          ) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {n: set() for n, _ in cap_atoms}
    for a, b in cap_bonds:
        out[a].add(b); out[b].add(a)
    return out


def canonical_adj_for_variant(variant_name: str,
                               variant_atoms: list[tuple[str, str]]
                               ) -> dict[str, set[str]]:
    """Adjacency for HIS variants (HID/HIE/HIP)."""
    atoms = {n for n, _ in variant_atoms}
    bonds: list[tuple[str, str]] = []

    def add(a: str, b: str) -> None:
        if a in atoms and b in atoms:
            bonds.append(tuple(sorted([a, b])))

    add("N","CA"); add("CA","C"); add("C","O")
    if "H"  in atoms: add("N","H")
    if "HA" in atoms: add("CA","HA")
    add("CA","CB")
    # Chi1: N-CA-CB-CG  (same as base HIS)
    add("CB","CG")
    # Chi2: CA-CB-CG-ND1 (HIS chi2 anchor)
    add("CG","ND1")
    # HIS imidazole ring (5-membered).
    add("CG","CD2"); add("CD2","NE2"); add("NE2","CE1"); add("CE1","ND1")
    # HBs.
    add("HB2","CB"); add("HB3","CB")
    # CD2 has HD2 (in all variants).
    add("HD2","CD2"); add("HE1","CE1")
    # Variant-specific.
    if "HD1" in atoms: add("HD1","ND1")
    if "HE2" in atoms: add("HE2","NE2")

    out: dict[str, set[str]] = {n: set() for n, _ in variant_atoms}
    for a, b in bonds:
        out[a].add(b); out[b].add(a)
    return out


# ---------------------------------------------------------------------------
# Graph isomorphism (residue-scale; canonical is rigid)
# ---------------------------------------------------------------------------

def canonical_wl_sigs(canon_adj: dict[str, set[str]],
                       atom_elem: dict[str, str],
                       rounds: int = 3) -> dict[str, int]:
    """K=3 Weisfeiler-Lehman per-atom signature ids for the canonical
    chemistry. After K=3, all chemically-distinct atoms in the standard
    20 residues are in singleton signature classes; the only residual
    multi-atom classes are graph-automorphic pairs (PHE/TYR CD1↔CD2,
    ARG NH1↔NH2 etc.) which spatial tiebreak at match time resolves.
    """
    cur: dict[str, int] = {}
    for n, e in atom_elem.items():
        cur[n] = hash((e, 0))
    for _ in range(rounds):
        nxt: dict[str, int] = {}
        for n in atom_elem:
            nbr = tuple(sorted((atom_elem[m], cur[m] & 0x7fffffff)
                                for m in canon_adj[n]))
            nxt[n] = hash((atom_elem[n], len(canon_adj[n]), nbr))
        cur = nxt
    return cur


def perceived_wl_sigs(piece_atoms: list[int],
                       sub_adj: list[set[int]],
                       elements: list[str],
                       rounds: int = 3) -> dict[int, int]:
    """K=3 WL signatures for the perceived atoms."""
    cur: dict[int, int] = {i: hash((elements[i], 0)) for i in piece_atoms}
    for _ in range(rounds):
        nxt: dict[int, int] = {}
        for i in piece_atoms:
            nbr = tuple(sorted((elements[j], cur[j] & 0x7fffffff)
                                for j in sub_adj[i]))
            nxt[i] = hash((elements[i], len(sub_adj[i]), nbr))
        cur = nxt
    return cur


def match_piece(
    piece_atoms: list[int],
    full_elements: list[str],
    full_adj: list[set[int]],
    canon_atoms: list[tuple[str, str]],
    canon_adj: dict[str, set[str]],
    debug: bool = False,
    debug_label: str = "",
) -> Optional[tuple[dict[int, str], dict[int, bool]]]:
    """Match perceived atoms (0-indexed in full record) to canonical
    atom names. Returns ``(name_by_perceived, ambiguous_by_perceived)``
    or None on failure.

    ``ambiguous_by_perceived[p]`` is True iff the canonical name for
    perceived atom ``p`` was chosen from a WL signature class of size
    ≥ 2 — i.e., the graph isomorphism could not distinguish the atom
    from a sibling (PHE/TYR CD1↔CD2, ARG NH1↔NH2, ASN HD21↔HD22,
    methyl-Hs whose DiastereotopicIndex collapses). Singleton classes
    (the common case after K=3 WL, including chemistry-distinct
    branches like ILE CG1/CG2) flag False — the canonical
    BranchAddress + DiastereotopicIndex are deterministic.

    The flag is consumed downstream by the assembler: ambiguous=True
    relaxes the typed-identity equality to (Element, Locant,
    BackboneRole) and resolves the within-class assignment by
    nearest-spatial. ambiguous=False keeps strict identity binding.
    """
    canon_elem = {n: e for n, e in canon_atoms}
    # Build adjacency restricted to this piece (in the full graph).
    in_piece = set(piece_atoms)
    sub_adj: list[set[int]] = [set() for _ in range(len(full_adj))]
    for i in piece_atoms:
        for j in full_adj[i]:
            if j in in_piece:
                sub_adj[i].add(j)

    canon_sig = canonical_wl_sigs(canon_adj, canon_elem)
    perceived_sig = perceived_wl_sigs(piece_atoms, sub_adj, full_elements)

    # Group canonical atoms by WL signature.
    canon_by_sig: dict[int, list[str]] = {}
    for name, s in canon_sig.items():
        canon_by_sig.setdefault(s, []).append(name)
    perceived_by_sig: dict[int, list[int]] = {}
    for ai, s in perceived_sig.items():
        perceived_by_sig.setdefault(s, []).append(ai)

    if debug:
        print(f"  [debug {debug_label}] canonical sigs:", file=sys.stderr)
        for sig, names in sorted(canon_by_sig.items()):
            print(f"    {sig} -> {names}", file=sys.stderr)
        print(f"  [debug {debug_label}] perceived sigs:", file=sys.stderr)
        for sig, idxs in sorted(perceived_by_sig.items()):
            print(f"    {sig} -> {idxs}", file=sys.stderr)

    if set(canon_by_sig.keys()) != set(perceived_by_sig.keys()):
        only_canon = set(canon_by_sig.keys()) - set(perceived_by_sig.keys())
        only_perceived = set(perceived_by_sig.keys()) - set(canon_by_sig.keys())
        if debug:
            print(f"  [debug {debug_label}] signature sets differ", file=sys.stderr)
            print(f"    only canonical: {only_canon}", file=sys.stderr)
            print(f"    only perceived: {only_perceived}", file=sys.stderr)
        return None  # signature sets differ → not isomorphic

    names_out: dict[int, str] = {}
    ambig_out: dict[int, bool] = {}
    for sig in canon_by_sig:
        cs = canon_by_sig[sig]; ps = perceived_by_sig[sig]
        if len(cs) != len(ps):
            if debug:
                print(f"  [debug {debug_label}] sig {sig} cardinality "
                      f"canonical={len(cs)} perceived={len(ps)}",
                      file=sys.stderr)
            return None  # cardinality mismatch within signature class
        ambiguous = len(cs) > 1
        # Map perceived atoms to canonical atoms 1:1 in arbitrary order
        # within the equivalence class.
        for p_idx, c_name in zip(ps, cs):
            names_out[p_idx] = c_name
            ambig_out[p_idx] = ambiguous

    return names_out, ambig_out


# ---------------------------------------------------------------------------
# Top-level perception
# ---------------------------------------------------------------------------

ONE_LETTER_TO_THREE = {
    "A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","Q":"GLN","E":"GLU",
    "G":"GLY","H":"HIS","I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE",
    "P":"PRO","S":"SER","T":"THR","W":"TRP","Y":"TYR","V":"VAL",
}


def perceive(atoms: list[dict], expected_central_one_letter: str,
              debug: bool = False) -> dict:
    """Run the full perception pipeline. Returns a result dict with
    'ok' (bool), 'reason' (str if not ok), and 'pieces' (list of dicts
    with 'kind' and per-atom assignments)."""
    elements = [a["element"] for a in atoms]

    adj = build_bond_graph(atoms)
    amides_all = find_amide_bonds(atoms, adj)

    # ASN/GLN sidechain primary amides are also detected. Filter to
    # peptide amides: a peptide amide N has at least one heavy neighbour
    # other than the carbonyl C (the CA of the next residue, or CD for
    # PRO). A sidechain primary amide N (ASN ND2, GLN NE2) has only
    # H neighbours besides the carbonyl C.
    amides: list[tuple[int, int]] = []
    for c, n in amides_all:
        n_heavy_other = sum(
            1 for j in adj[n] if j != c and elements[j] != "H"
        )
        if n_heavy_other >= 1:
            amides.append((c, n))

    if len(amides) != 4:
        return {"ok": False,
                "reason": f"expected 4 peptide amides, got {len(amides)} "
                          f"(all-amides={len(amides_all)})"}

    # Cut amides → connected components.
    cut_adj = [s.copy() for s in adj]
    for c, nj in amides:
        cut_adj[c].discard(nj); cut_adj[nj].discard(c)
    comps = connected_components(cut_adj, len(atoms))
    if len(comps) != 5:
        return {"ok": False, "reason": f"expected 5 pieces, got {len(comps)} "
                                        f"(sizes={[len(c) for c in comps]})"}

    ordered = order_pieces_along_chain(comps, amides, adj)
    # ordered[0] = ACE, [1] = NCapAla, [2] = Central, [3] = CCapAla, [4] = NME

    kinds = ["ACE", "NCapAla", "Central", "CCapAla", "NME"]
    expected_central_three = ONE_LETTER_TO_THREE[expected_central_one_letter]

    pieces_out: list[dict] = []
    for kind, piece in zip(kinds, ordered):
        # Initialise canon_atoms/canon_adj and the dual outputs from
        # match_piece (names + per-perceived-atom ambiguous flag).
        canon_atoms: list[tuple[str, str]] = []
        canon_adj: dict[str, set[str]] = {}
        names_map: Optional[dict[int, str]] = None
        ambig_map: Optional[dict[int, bool]] = None

        if kind == "ACE":
            canon_atoms = ACE_ATOMS; canon_adj = canonical_adj_for_cap(ACE_ATOMS, ACE_BONDS)
        elif kind == "NME":
            canon_atoms = NME_ATOMS; canon_adj = canonical_adj_for_cap(NME_ATOMS, NME_BONDS)
        elif kind in ("NCapAla", "CCapAla"):
            canon_atoms = CANONICAL["ALA"]; canon_adj = canonical_adj("ALA")
        elif expected_central_three == "HIS":
            # HIS Central: try each HIS variant; pick the one matching.
            variant_match: Optional[tuple[
                str, list[tuple[str, str]], dict[str, set[str]],
                dict[int, str], dict[int, bool]]] = None
            for variant_name, variant_atoms in HIS_VARIANTS.items():
                if len(piece) != len(variant_atoms):
                    continue
                var_canon_adj = canonical_adj_for_variant(variant_name,
                                                           variant_atoms)
                m = match_piece(piece, elements, adj, variant_atoms,
                                var_canon_adj, debug=debug,
                                debug_label=f"Central/{variant_name}")
                if m is not None:
                    variant_match = (variant_name, variant_atoms,
                                      var_canon_adj, m[0], m[1]); break
            if variant_match is None:
                return {"ok": False,
                        "reason": f"Central: no HIS variant matched "
                                  f"perceived_n={len(piece)}"}
            _, canon_atoms, canon_adj, names_map, ambig_map = variant_match
        else:  # Central, not HIS
            canon_atoms = CANONICAL[expected_central_three]
            canon_adj   = canonical_adj(expected_central_three)
            if len(piece) != len(canon_atoms):
                return {"ok": False,
                        "reason": f"{kind}: atom-count mismatch perceived={len(piece)} canonical={len(canon_atoms)}"}
            res = match_piece(piece, elements, adj, canon_atoms,
                              canon_adj, debug=debug,
                              debug_label=kind)
            if res is None:
                return {"ok": False,
                        "reason": f"{kind}: graph-isomorphism failed"}
            names_map, ambig_map = res

        # For non-Central pieces, do the count + match steps now.
        if kind != "Central":
            if len(piece) != len(canon_atoms):
                return {"ok": False,
                        "reason": f"{kind}: atom-count mismatch perceived={len(piece)} canonical={len(canon_atoms)}"}
            res = match_piece(piece, elements, adj, canon_atoms,
                              canon_adj, debug=debug,
                              debug_label=kind)
            if res is None:
                return {"ok": False,
                        "reason": f"{kind}: graph-isomorphism failed"}
            names_map, ambig_map = res

        # mypy: names_map / ambig_map are populated on every reachable path above.
        assert names_map is not None and ambig_map is not None
        pieces_out.append({
            "kind": kind,
            "residue": "ALA" if kind in ("NCapAla","CCapAla")
                       else "ACE" if kind == "ACE"
                       else "NME" if kind == "NME"
                       else expected_central_three,
            "atoms": [{"dft_atom_idx": atoms[i]["atom_idx"],
                       "element": atoms[i]["element"],
                       "canonical_name": names_map[i],
                       "canonical_assignment_ambiguous": ambig_map[i]}
                      for i in piece],
        })

    return {"ok": True, "pieces": pieces_out}


# ---------------------------------------------------------------------------
# Source-log parser (optional cross-check)
# ---------------------------------------------------------------------------

def parse_standard_orientation_block(path: Path) -> list[dict]:
    """Parse the LAST Standard orientation block from a Gaussian log.

    Returns list of {atom_idx, element, x, y, z}.
    """
    text = path.read_text()
    # Find all Standard orientation blocks; we want the last one
    # (the optimised geometry).
    matches = list(re.finditer(r"Standard orientation:", text))
    if not matches:
        raise RuntimeError(f"no Standard orientation block in {path}")
    start = matches[-1].end()
    block = text[start:]
    # Skip past the Standard orientation header block (3 lines + 2
    # dashed separators) to the atom rows.
    #     "Standard orientation:"   <- already past
    #     "---------------------"   <- first dashed
    #     "Center  Atomic  Atomic   Coordinates (Angstroms)"
    #     "Number  Number   Type      X     Y     Z"
    #     "---------------------"   <- second dashed
    #     <atom rows>
    #     "---------------------"   <- closing dashed
    lines = block.splitlines()
    i = 0
    while i < len(lines) and not lines[i].lstrip().startswith("---"):
        i += 1
    i += 1  # past first ---
    i += 2  # past 2 header text lines
    i += 1  # past second ---
    Z_TO_ELEM = {1:"H", 6:"C", 7:"N", 8:"O", 16:"S"}
    out: list[dict] = []
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("---"): break
        parts = line.split()
        if len(parts) < 6: break
        idx = int(parts[0]); z = int(parts[1])
        x = float(parts[3]); y = float(parts[4]); z_coord = float(parts[5])
        out.append({"atom_idx": idx, "element": Z_TO_ELEM[z],
                    "x": x, "y": y, "z": z_coord})
        i += 1
    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dump", type=Path, required=True,
                    help="Path to tripeptide_orderings_dump.json")
    ap.add_argument("--source-log", type=Path, default=None,
                    help="Optional: parse this Gaussian log and run perception "
                         "on the Standard orientation block instead of the DB row")
    ap.add_argument("--source-log-central", type=str, default="A",
                    help="One-letter code for the central residue in the source log")
    ap.add_argument("--verbose", action="store_true",
                    help="Print per-atom assignments")
    ap.add_argument("--debug", type=str, default=None,
                    help="Print signature-class debug info for the given "
                         "tripeptide (e.g., AFA, AHA)")
    args = ap.parse_args()

    if args.source_log:
        atoms = parse_standard_orientation_block(args.source_log)
        print(f"parsed {len(atoms)} atoms from {args.source_log}", file=sys.stderr)
        result = perceive(atoms, args.source_log_central)
        if not result["ok"]:
            print(f"FAIL: {result['reason']}", file=sys.stderr); return 1
        print(f"PERCEIVED from source log {args.source_log}:", file=sys.stderr)
        for piece in result["pieces"]:
            print(f"  {piece['kind']:8} ({piece['residue']:3}) "
                  f"n={len(piece['atoms'])}", file=sys.stderr)
            if args.verbose:
                for a in piece["atoms"]:
                    print(f"      idx={a['dft_atom_idx']:3} {a['element']} "
                          f"→ {a['canonical_name']}", file=sys.stderr)
        return 0

    with args.dump.open() as fh:
        rows = json.load(fh)

    all_ok = True
    for row in sorted(rows, key=lambda r: (r["tripeptide"], r["frame_type"])):
        r_letter = row["central_residue"]
        debug = (row["tripeptide"] == args.debug)
        result = perceive(row["atoms"], r_letter, debug=debug)
        tag = "OK " if result["ok"] else "FAIL"
        ratio = f"n={row['n_atoms']:>3} calc_id={row['calc_id']}"
        if result["ok"]:
            print(f"  {tag}  {row['tripeptide']:3} {row['frame_type']:32} {ratio}",
                  file=sys.stderr)
            if args.verbose:
                for piece in result["pieces"]:
                    print(f"      {piece['kind']:8} ({piece['residue']:3}) "
                          f"n={len(piece['atoms'])}", file=sys.stderr)
                    for a in piece["atoms"]:
                        print(f"        idx={a['dft_atom_idx']:3} "
                              f"{a['element']} → {a['canonical_name']}",
                              file=sys.stderr)
        else:
            all_ok = False
            print(f"  {tag}  {row['tripeptide']:3} {row['frame_type']:32} {ratio}  "
                  f"reason={result['reason']}", file=sys.stderr)

    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())

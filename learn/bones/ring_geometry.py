"""
Ring geometry extraction from atomic positions.

Identifies aromatic rings from residue type and atom element patterns,
computes ring center, normal, and per-atom cylindrical coordinates (ρ, z, θ).

This is a Python stopgap. The C++ code has proper ring detection via
CovalentTopology and SVD-based normals. A dedicated geometry extractor
should replace this when we move past the extract-train-explore cycle.
"""

import numpy as np
from dataclasses import dataclass
from typing import List, Optional


# Aromatic residue types (matching C++ AminoAcid enum order)
# 0=ALA,1=ARG,2=ASN,3=ASP,4=CYS,5=GLN,6=GLU,7=GLY,
# 8=HIS,9=ILE,10=LEU,11=LYS,12=MET,13=PHE,14=PRO,15=SER,
# 16=THR,17=TRP,18=TYR,19=VAL
AROMATIC_RESIDUE_TYPES = {8, 13, 17, 18}  # HIS, PHE, TRP, TYR

# Ring carbon/nitrogen counts by residue type
# PHE: 6 ring C
# TYR: 6 ring C
# TRP: 9 ring atoms (8C + 1N, two fused rings)
# HIS: 5 ring atoms (3C + 2N)


@dataclass
class Ring:
    """One aromatic ring with geometry."""
    residue_index: int
    residue_type: int
    center: np.ndarray      # (3,)
    normal: np.ndarray       # (3,) unit vector
    radius: float            # mean distance from center to vertices
    vertices: np.ndarray     # (n_vertices, 3)
    ring_type: str           # "PHE", "TYR", "TRP6", "TRP5", "HIS", etc.


@dataclass
class AtomRingCoords:
    """Cylindrical coordinates of one atom relative to one ring."""
    ring_idx: int
    ring_type: str
    rho: float        # radial distance in ring plane (Å)
    z: float          # signed height above ring plane (Å)
    distance: float   # total distance to ring center (Å)
    theta: float      # cone angle from normal: atan2(rho, |z|)


def _fit_ring_normal(vertices: np.ndarray) -> np.ndarray:
    """SVD-based ring normal (same approach as C++ GeometryResult)."""
    centroid = vertices.mean(axis=0)
    centered = vertices - centroid
    _, _, Vt = np.linalg.svd(centered)
    normal = Vt[2]  # smallest singular value = normal direction
    # Convention: normal points toward positive z (arbitrary but consistent)
    if normal[2] < 0:
        normal = -normal
    return normal / np.linalg.norm(normal)


def find_rings(pos: np.ndarray, element: np.ndarray,
               residue_index: np.ndarray, residue_type: np.ndarray) -> List[Ring]:
    """Find aromatic rings from atom data.

    Identifies ring atoms as non-hydrogen atoms in aromatic residues
    that are carbon or nitrogen (ring atoms, not sidechain CH2 etc.).
    Groups by residue and fits ring planes.
    """
    rings = []
    n_atoms = len(pos)

    # Group atoms by residue
    residue_atoms = {}
    for i in range(n_atoms):
        ri = residue_index[i]
        if ri not in residue_atoms:
            residue_atoms[ri] = []
        residue_atoms[ri].append(i)

    for ri, atom_indices in residue_atoms.items():
        rt = residue_type[atom_indices[0]]
        if rt not in AROMATIC_RESIDUE_TYPES:
            continue

        # Get ring atoms: C and N only (element 6, 7), no H
        ring_atoms = [i for i in atom_indices
                      if element[i] in (6, 7)]

        if rt == 13:  # PHE: 6-membered ring (6 C)
            ring_c = [i for i in ring_atoms if element[i] == 6]
            # PHE has 6 ring C + 1 Cβ + 1 Cα. Ring C are further from backbone.
            # Simple heuristic: take the 6 C furthest from Cα
            if len(ring_c) >= 6:
                verts = pos[ring_c]
                center = verts.mean(axis=0)
                # Sort by distance to center, take closest 6 (ring, not Cα/Cβ)
                dists = np.linalg.norm(verts - center, axis=1)
                ring_idx = np.argsort(dists)[:6]
                verts = verts[ring_idx]
                normal = _fit_ring_normal(verts)
                center = verts.mean(axis=0)
                radius = np.linalg.norm(verts - center, axis=1).mean()
                rings.append(Ring(ri, rt, center, normal, radius, verts, "PHE"))

        elif rt == 18:  # TYR: same as PHE but with OH
            ring_c = [i for i in ring_atoms if element[i] == 6]
            if len(ring_c) >= 6:
                verts = pos[ring_c]
                center = verts.mean(axis=0)
                dists = np.linalg.norm(verts - center, axis=1)
                ring_idx = np.argsort(dists)[:6]
                verts = verts[ring_idx]
                normal = _fit_ring_normal(verts)
                center = verts.mean(axis=0)
                radius = np.linalg.norm(verts - center, axis=1).mean()
                rings.append(Ring(ri, rt, center, normal, radius, verts, "TYR"))

        elif rt == 8:  # HIS: 5-membered ring (3C + 2N)
            if len(ring_atoms) >= 5:
                verts = pos[ring_atoms[:5]]
                normal = _fit_ring_normal(verts)
                center = verts.mean(axis=0)
                radius = np.linalg.norm(verts - center, axis=1).mean()
                rings.append(Ring(ri, rt, center, normal, radius, verts, "HIS"))

        elif rt == 17:  # TRP: two fused rings
            # TRP has 9 ring atoms. We treat it as one 9-membered system
            # (matching C++ TRP9 ring type for the fused system)
            if len(ring_atoms) >= 9:
                verts = pos[ring_atoms[:9]]
                normal = _fit_ring_normal(verts)
                center = verts.mean(axis=0)
                radius = np.linalg.norm(verts - center, axis=1).mean()
                rings.append(Ring(ri, rt, center, normal, radius, verts, "TRP"))

    return rings


def atom_ring_coordinates(pos: np.ndarray, rings: List[Ring],
                          max_dist: float = 15.0) -> List[List[AtomRingCoords]]:
    """Compute cylindrical coordinates of each atom relative to each ring.

    Returns list of lists: result[atom_index] = [AtomRingCoords, ...]
    """
    n_atoms = len(pos)
    result = [[] for _ in range(n_atoms)]

    for ri, ring in enumerate(rings):
        # Vector from ring center to each atom
        d = pos - ring.center  # (N, 3)
        distances = np.linalg.norm(d, axis=1)

        # Only process atoms within max_dist
        nearby = np.where(distances < max_dist)[0]

        for ai in nearby:
            z = np.dot(d[ai], ring.normal)  # height above plane (signed)
            in_plane = d[ai] - z * ring.normal
            rho = np.linalg.norm(in_plane)  # radial distance in plane
            theta = np.arctan2(rho, abs(z))  # cone angle from normal

            result[ai].append(AtomRingCoords(
                ring_idx=ri,
                ring_type=ring.ring_type,
                rho=rho,
                z=z,
                distance=distances[ai],
                theta=theta,
            ))

    return result

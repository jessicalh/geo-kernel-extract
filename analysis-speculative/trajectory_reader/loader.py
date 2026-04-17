"""
Loader: reads an analysis trajectory H5 file and returns an AnalysisProtein.

The H5 is the only input. No external config, no TOML, no manifest.
Groups that don't exist in the H5 are silently skipped (the corresponding
attribute on AnalysisProtein will be None).
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import h5py
import numpy as np

from .enums import (
    Element,
    AtomRole,
    Hybridisation,
    BondOrder,
    BondCategory,
    RingTypeIndex,
    amino_acid_from_code,
)
from .model import (
    AnalysisProtein,
    AnalysisAtom,
    AnalysisResidue,
    AnalysisRing,
    AnalysisBond,
    RingCurrentData,
    EfgData,
    BondAnisoData,
    QuadrupoleData,
    DispersionData,
    HBondData,
    SasaData,
    WaterData,
    ChargesData,
    AimNet2EmbeddingData,
    PredictionsData,
    ProjectionsData,
    DihedralData,
    DsspData,
    TimeSeries,
)


def _read_group(h5: h5py.File, group_name: str, cls: type[TimeSeries]) -> Optional[TimeSeries]:
    """
    Read all datasets in an H5 group into a TimeSeries subclass.

    Returns None if the group does not exist.
    Reads all datasets eagerly into numpy arrays.
    """
    if group_name not in h5:
        return None
    grp = h5[group_name]
    arrays = {}
    for name in grp:
        ds = grp[name]
        if isinstance(ds, h5py.Dataset):
            arrays[name] = ds[()]
    if not arrays:
        return None
    return cls(_arrays=arrays)


def _decode_string(val: bytes | str) -> str:
    """Decode a value that may be bytes or str from h5py."""
    if isinstance(val, bytes):
        return val.decode("utf-8")
    return str(val)


def load_analysis(path: str | Path) -> AnalysisProtein:
    """
    Load an analysis trajectory H5 file into an AnalysisProtein.

    Parameters
    ----------
    path : str or Path
        Path to the md_analysis.h5 file.

    Returns
    -------
    AnalysisProtein
        Fully populated protein with topology and time series data.

    Raises
    ------
    FileNotFoundError
        If the H5 file does not exist.
    KeyError
        If required groups (meta, atoms, residues, positions, topology) are missing.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"H5 file not found: {path}")

    with h5py.File(path, "r") as h5:
        # ---------------------------------------------------------------
        # Meta
        # ---------------------------------------------------------------
        protein_id = _read_meta_string(h5, "protein_id")
        frame_times = h5["meta/frame_times"][()]
        frame_indices = h5["meta/frame_indices"][()]
        n_atoms = int(h5["meta"].attrs.get("n_atoms", 0))
        n_residues = int(h5["meta"].attrs.get("n_residues", 0))

        # ---------------------------------------------------------------
        # Positions (T, N, 3)
        # ---------------------------------------------------------------
        positions = h5["positions/xyz"][()]

        # ---------------------------------------------------------------
        # Atoms
        # ---------------------------------------------------------------
        atom_elements = h5["atoms/element"][()]
        atom_names = h5["atoms/atom_name"][()]
        atom_roles = h5["atoms/atom_role"][()]
        atom_hyb = h5["atoms/hybridisation"][()]
        atom_residue_idx = h5["atoms/residue_index"][()]
        atom_is_backbone = h5["atoms/is_backbone"][()]
        atom_is_conjugated = h5["atoms/is_conjugated"][()]
        atom_n_bonded = h5["atoms/n_bonded"][()]
        atom_graph_dist_ring = h5["atoms/graph_dist_ring"][()]

        # parent_atom_index: stored as int32, -1 = no parent (heavy atom)
        parent_atom_idx = h5["topology/parent_atom_index"][()]

        atoms: list[AnalysisAtom] = []
        for i in range(len(atom_elements)):
            elem_val = int(atom_elements[i])
            try:
                element = Element(elem_val)
            except ValueError:
                element = Element.Unknown

            role_val = int(atom_roles[i])
            try:
                role = AtomRole(role_val)
            except ValueError:
                role = AtomRole.Unknown

            hyb_val = int(atom_hyb[i])
            try:
                hyb = Hybridisation(hyb_val)
            except ValueError:
                hyb = Hybridisation.Unassigned

            atoms.append(AnalysisAtom(
                index=i,
                element=element,
                pdb_name=_decode_string(atom_names[i]),
                role=role,
                hybridisation=hyb,
                residue_index=int(atom_residue_idx[i]),
                is_backbone=bool(atom_is_backbone[i]),
                is_conjugated=bool(atom_is_conjugated[i]),
                n_bonded=int(atom_n_bonded[i]),
                graph_dist_ring=int(atom_graph_dist_ring[i]),
                parent_atom_index=int(parent_atom_idx[i]),
            ))

        # ---------------------------------------------------------------
        # Residues
        # ---------------------------------------------------------------
        res_names = h5["residues/residue_name"][()]
        res_numbers = h5["residues/residue_number"][()]

        residues: list[AnalysisResidue] = []
        for i in range(len(res_names)):
            name = _decode_string(res_names[i])
            residues.append(AnalysisResidue(
                index=i,
                amino_acid=amino_acid_from_code(name),
                name=name,
                sequence_number=int(res_numbers[i]),
            ))

        # Wire atoms -> residues (build atom_indices per residue)
        for atom in atoms:
            if 0 <= atom.residue_index < len(residues):
                residues[atom.residue_index].atom_indices.append(atom.index)

        # ---------------------------------------------------------------
        # Bonds
        # ---------------------------------------------------------------
        bonds: list[AnalysisBond] = []
        if "topology/bond_atoms" in h5:
            bond_atoms = h5["topology/bond_atoms"][()]
            bond_orders = h5["topology/bond_order"][()]
            bond_categories = h5["topology/bond_category"][()]

            for i in range(len(bond_atoms)):
                a_idx = int(bond_atoms[i, 0])
                b_idx = int(bond_atoms[i, 1])

                order_val = int(bond_orders[i])
                try:
                    order = BondOrder(order_val)
                except ValueError:
                    order = BondOrder.Unknown

                cat_val = int(bond_categories[i])
                try:
                    category = BondCategory(cat_val)
                except ValueError:
                    category = BondCategory.Unknown

                bonds.append(AnalysisBond(
                    index=i,
                    atom_index_a=a_idx,
                    atom_index_b=b_idx,
                    order=order,
                    category=category,
                ))

            # Wire bonds -> atoms
            for bond in bonds:
                if 0 <= bond.atom_index_a < len(atoms):
                    atoms[bond.atom_index_a].bond_indices.append(bond.index)
                if 0 <= bond.atom_index_b < len(atoms):
                    atoms[bond.atom_index_b].bond_indices.append(bond.index)

        # ---------------------------------------------------------------
        # Rings
        # ---------------------------------------------------------------
        rings: list[AnalysisRing] = []
        if "topology/ring_type" in h5:
            ring_types = h5["topology/ring_type"][()]
            ring_residues = h5["topology/ring_residue"][()]
            ring_atom_indices_flat = h5["topology/ring_atom_indices"][()]
            ring_offsets = h5["topology/ring_offsets"][()]

            n_rings = len(ring_types)
            for i in range(n_rings):
                start = int(ring_offsets[i])
                end = int(ring_offsets[i + 1]) if i + 1 < len(ring_offsets) else len(ring_atom_indices_flat)

                type_val = int(ring_types[i])
                try:
                    ring_type = RingTypeIndex(type_val)
                except ValueError:
                    ring_type = RingTypeIndex.PheBenzene

                atom_idxs = [int(x) for x in ring_atom_indices_flat[start:end]]

                rings.append(AnalysisRing(
                    index=i,
                    type_index=ring_type,
                    atom_indices=atom_idxs,
                    parent_residue_index=int(ring_residues[i]),
                ))

        # ---------------------------------------------------------------
        # Build the protein
        # ---------------------------------------------------------------
        protein = AnalysisProtein(
            protein_id=protein_id,
            atoms=atoms,
            residues=residues,
            rings=rings,
            bonds=bonds,
            positions=positions,
            frame_times=frame_times,
            frame_indices=frame_indices,
        )

        # ---------------------------------------------------------------
        # Load time series groups (skip missing ones)
        # ---------------------------------------------------------------
        protein.ring_current = _read_group(h5, "ring_current", RingCurrentData)
        protein.efg = _read_group(h5, "efg", EfgData)
        protein.bond_aniso = _read_group(h5, "bond_aniso", BondAnisoData)
        protein.quadrupole = _read_group(h5, "quadrupole", QuadrupoleData)
        protein.dispersion = _read_group(h5, "dispersion", DispersionData)
        protein.hbond = _read_group(h5, "hbond", HBondData)
        protein.sasa = _read_group(h5, "sasa", SasaData)
        protein.water = _read_group(h5, "water", WaterData)
        protein.charges = _read_group(h5, "charges", ChargesData)
        protein.aimnet2_embedding = _read_group(h5, "aimnet2_embedding", AimNet2EmbeddingData)
        protein.predictions = _read_group(h5, "predictions", PredictionsData)
        protein.projections = _read_group(h5, "projections", ProjectionsData)
        protein.dihedrals = _read_group(h5, "dihedrals", DihedralData)
        protein.dssp = _read_group(h5, "dssp", DsspData)

    return protein


def _read_meta_string(h5: h5py.File, key: str) -> str:
    """
    Read a metadata string, checking both dataset and attribute locations.

    The spec says protein_id is a scalar dataset under meta/, but the actual
    H5 stores it as an attribute on the meta/ group. We check both.
    """
    # Check as dataset first
    ds_path = f"meta/{key}"
    if ds_path in h5:
        val = h5[ds_path][()]
        if isinstance(val, bytes):
            return val.decode("utf-8")
        return str(val)
    # Check as attribute on meta group
    if "meta" in h5 and key in h5["meta"].attrs:
        val = h5["meta"].attrs[key]
        if isinstance(val, bytes):
            return val.decode("utf-8")
        return str(val)
    # Fallback: derive from filename
    return Path(h5.filename).stem

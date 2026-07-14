"""Synthetic SDK coverage for the complete diskless-MOPAC NPY group."""

from __future__ import annotations

from dataclasses import fields

import numpy as np

from nmr_extract import load
from nmr_extract._protein import MopacFullGroup
from _topology_fixture import (
    write_minimal_topology_sidecar,
    write_required_sdk_npys,
)


N = 3
W = 4
E = 7
U = 2
O = 4
V = 2


F64_SHAPES = {
    "mopac_charges": (N,),
    "mopac_charges_full_precision": (N,),
    "mopac_scalars": (N, 4),
    "mopac_bond_orders": (U, 3),
    "mopac_bond_orders_full_precision": (U, 3),
    "mopac_bond_neighbors": (2 * U, 4),
    "mopac_global": (4,),
    "mopac_atom_populations": (N, 12),
    "mopac_atomic_orbital_populations": (N, 9),
    "mopac_atomic_orbital_population_totals": (N, 3),
    "mopac_bond_valencies": (N,),
    "mopac_bond_valencies_full_precision": (N,),
    "mopac_bond_orders_unique": (U, 8),
    "mopac_topology_bond_orders_full": (U, 8),
    "mopac_heat_kcal_mol": (1,),
    "mopac_dipole_debye": (3,),
    "mopac_dipole_point_charge_debye": (3,),
    "mopac_dipole_hybridization_debye": (3,),
    "mopac_bond_order": (E,),
    "mopac_atom_ao_density": (N, W, W),
    "mopac_atomic_orbital_populations_full_precision": (N, W),
    "mopac_bond_ao_density_directed": (E, W, W),
    "mopac_bond_ao_density": (U, W, W),
    "mopac_atom_electron_population": (N,),
    "mopac_atom_s_population": (N,),
    "mopac_atom_p_population": (N,),
    "mopac_atom_d_population": (N,),
    "mopac_lmo_energy_levels": (O + V,),
    "mopac_lmo_occupied_coefficients": (24,),
    "mopac_lmo_virtual_coefficients": (12,),
    "mopac_lmo_occupied_coefficient_storage_native": (1104,),
    "mopac_lmo_virtual_coefficient_storage_native": (1052,),
}


I32_SHAPES = {
    "mopac_bond_index": (N + 1,),
    "mopac_bond_atom": (E,),
    "mopac_ao_max_orbitals": (1,),
    "mopac_ao_orbitals_per_atom": (N,),
    "mopac_bond_density_pairs": (U, 2),
    "mopac_lewis_bond_count": (N,),
    "mopac_lewis_bond_atoms": (N, 9),
    "mopac_lmo_occupied_atom_counts": (O,),
    "mopac_lmo_occupied_atoms": (12,),
    "mopac_lmo_virtual_atom_counts": (V,),
    "mopac_lmo_virtual_atoms": (6,),
    "mopac_lmo_occupied_atom_offsets_native": (O,),
    "mopac_lmo_virtual_atom_offsets_native": (V,),
    "mopac_lmo_occupied_coefficient_offsets_native": (O,),
    "mopac_lmo_virtual_coefficient_offsets_native": (V,),
    "mopac_lmo_occupied_atom_storage_native": (1052,),
    "mopac_lmo_virtual_atom_storage_native": (1026,),
    "mopac_mozyme_state_dimensions": (7,),
}


def _save_sentinel(path, stem: str, shape, dtype, sentinel: int) -> None:
    np.save(path / f"{stem}.npy", np.full(shape, sentinel, dtype=dtype))


def test_complete_direct_mopac_group_loads_all_50_arrays(tmp_path):
    write_required_sdk_npys(tmp_path, N, n_residues=1)
    write_minimal_topology_sidecar(tmp_path, N, n_residues=1)

    sentinel_by_stem = {}
    for sentinel, (stem, shape) in enumerate(F64_SHAPES.items(), start=1):
        _save_sentinel(tmp_path, stem, shape, np.float64, sentinel)
        sentinel_by_stem[stem] = sentinel
    for sentinel, (stem, shape) in enumerate(I32_SHAPES.items(), start=101):
        _save_sentinel(tmp_path, stem, shape, np.int32, sentinel)
        sentinel_by_stem[stem] = sentinel

    assert len(sentinel_by_stem) == 50
    protein = load(tmp_path)
    assert protein.mopac is not None
    assert protein.mopac.full is not None

    full = protein.mopac.full
    assert isinstance(full, MopacFullGroup)
    assert len(fields(MopacFullGroup)) == 45
    for field in fields(MopacFullGroup):
        value = getattr(full, field.name)
        assert value is not None, field.name
        array = value.data if hasattr(value, "data") else value
        stem = f"mopac_{field.name}"
        assert np.asarray(array).flat[0] == sentinel_by_stem[stem], stem

    np.testing.assert_array_equal(
        protein.mopac.core.charges,
        np.full((N,), sentinel_by_stem["mopac_charges"], dtype=np.float64),
    )
    assert protein.mopac.core.scalars.data.shape == (N, 4)
    assert protein.mopac.core.bond_orders.data.shape == (U, 3)
    assert protein.mopac.core.bond_neighbors.shape == (2 * U, 4)
    assert protein.mopac.core.global_.data.shape == (4,)

#!/usr/bin/env python3
"""Strict NumPy/schema acceptance check for the production water probe."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from nmr_extract._catalog import CATALOG


def legacy_printed(value: float, precision: int) -> float:
    """Value recovered after the legacy parser read a Fortran F field."""
    return float(f"{float(value):.{precision}f}")


def legacy_compact_bond_rows(
        bond_index: np.ndarray,
        bond_atom: np.ndarray,
        bond_order: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Reproduce print_bonds_compact's physical first row per atom."""
    atom_count = len(bond_index) - 1
    rows: list[tuple[int, int, float]] = []
    directed_sums = np.zeros(atom_count, dtype=np.float64)
    for atom in range(atom_count):
        candidates = [
            (int(bond_atom[entry]), float(bond_order[entry]))
            for entry in range(int(bond_index[atom]),
                               int(bond_index[atom + 1]))
            if int(bond_atom[entry]) != atom
        ]
        # bonds_for_MOZYME supplied ascending atom indices before the Fortran
        # selection sort. Its inflated running maximum deliberately preserves
        # near-tie order.
        candidates.sort(key=lambda item: item[0])
        for start in range(len(candidates)):
            selected = start
            running = 0.0
            found = False
            for candidate in range(start, len(candidates)):
                value = candidates[candidate][1]
                if value > running:
                    selected = candidate
                    running = value * (1.0 + 1.0e-4)
                    found = True
            assert found
            candidates[start], candidates[selected] = (
                candidates[selected], candidates[start])

        for other, value in candidates[:6]:
            printed = legacy_printed(value, 3)
            if printed > 0.01:
                rows.append((atom, other, printed))
                directed_sums[atom] += printed

    return np.asarray(rows, dtype=np.float64).reshape((-1, 3)), directed_sums


def unordered_pair_projection(rows: np.ndarray) -> np.ndarray:
    """Legacy parser's maximum over the retained directed compact rows."""
    pairs: dict[tuple[int, int], float] = {}
    for atom_a, atom_b, order in rows:
        pair = tuple(sorted((int(atom_a), int(atom_b))))
        pairs[pair] = max(pairs.get(pair, -np.inf), float(order))
    return np.asarray(
        [(atom_a, atom_b, pairs[(atom_a, atom_b)])
         for atom_a, atom_b in sorted(pairs)],
        dtype=np.float64,
    ).reshape((-1, 3))


def rows_sorted_by_atoms(rows: np.ndarray) -> np.ndarray:
    """Canonicalize row sets whose production order is intentionally legacy."""
    if len(rows) == 0:
        return rows.copy()
    return rows[np.lexsort((rows[:, 1], rows[:, 0]))]


def main() -> None:
    root = Path(sys.argv[1])
    probe_archive = Path(sys.argv[2])
    probe_input = Path(sys.argv[3])
    n, width, entries = 3, 4, 7
    api_unique, compact_unique, topology = 2, 2, 2
    occupied, virtual = 4, 2

    f64_shapes = {
        "mopac_charges": (n,),
        "mopac_charges_full_precision": (n,),
        "mopac_scalars": (n, 4),
        "mopac_bond_orders": (compact_unique, 3),
        "mopac_bond_orders_full_precision": (api_unique, 3),
        "mopac_bond_neighbors": (2 * compact_unique, 4),
        "mopac_global": (4,),
        "mopac_atom_populations": (n, 12),
        "mopac_atomic_orbital_populations": (n, 9),
        "mopac_atomic_orbital_population_totals": (n, 3),
        "mopac_bond_valencies": (n,),
        "mopac_bond_valencies_full_precision": (n,),
        "mopac_bond_orders_unique": (api_unique, 8),
        "mopac_topology_bond_orders_full": (topology, 8),
        "mopac_heat_kcal_mol": (1,),
        "mopac_dipole_debye": (3,),
        "mopac_dipole_point_charge_debye": (3,),
        "mopac_dipole_hybridization_debye": (3,),
        "mopac_bond_order": (entries,),
        "mopac_atom_ao_density": (n, width, width),
        "mopac_atomic_orbital_populations_full_precision": (n, width),
        "mopac_bond_ao_density_directed": (entries, width, width),
        "mopac_bond_ao_density": (api_unique, width, width),
        "mopac_atom_electron_population": (n,),
        "mopac_atom_s_population": (n,),
        "mopac_atom_p_population": (n,),
        "mopac_atom_d_population": (n,),
        "mopac_lmo_energy_levels": (occupied + virtual,),
        "mopac_lmo_occupied_coefficients": (24,),
        "mopac_lmo_virtual_coefficients": (12,),
    }
    i32_shapes = {
        "mopac_bond_index": (n + 1,),
        "mopac_bond_atom": (entries,),
        "mopac_ao_max_orbitals": (1,),
        "mopac_ao_orbitals_per_atom": (n,),
        "mopac_bond_density_pairs": (api_unique, 2),
        "mopac_lewis_bond_count": (n,),
        "mopac_lewis_bond_atoms": (n, 9),
        "mopac_lmo_occupied_atom_counts": (occupied,),
        "mopac_lmo_occupied_atoms": (12,),
        "mopac_lmo_virtual_atom_counts": (virtual,),
        "mopac_lmo_virtual_atoms": (6,),
        "mopac_lmo_occupied_atom_offsets_native": (occupied,),
        "mopac_lmo_virtual_atom_offsets_native": (virtual,),
        "mopac_lmo_occupied_coefficient_offsets_native": (occupied,),
        "mopac_lmo_virtual_coefficient_offsets_native": (virtual,),
        "mopac_mozyme_state_dimensions": (7,),
    }

    state = np.load(root / "mopac_mozyme_state_dimensions.npy", allow_pickle=False)
    np.testing.assert_array_equal(state, [3, 4, 2, 1052, 1026, 1104, 1052])
    _, _, _, icocc_dim, icvir_dim, cocc_dim, cvir_dim = map(int, state)
    i32_shapes.update({
        "mopac_lmo_occupied_atom_storage_native": (icocc_dim,),
        "mopac_lmo_virtual_atom_storage_native": (icvir_dim,),
    })
    f64_shapes.update({
        "mopac_lmo_occupied_coefficient_storage_native": (cocc_dim,),
        "mopac_lmo_virtual_coefficient_storage_native": (cvir_dim,),
    })

    expected = set(f64_shapes) | set(i32_shapes)
    actual = {path.stem for path in root.glob("mopac_*.npy")}
    assert actual == expected, (sorted(expected - actual), sorted(actual - expected))
    assert len(expected) == 50

    arrays: dict[str, np.ndarray] = {}
    for stem, shape in f64_shapes.items():
        value = np.load(root / f"{stem}.npy", allow_pickle=False)
        assert value.dtype == np.dtype("<f8"), (stem, value.dtype)
        assert value.shape == shape, (stem, value.shape, shape)
        arrays[stem] = value
    for stem, shape in i32_shapes.items():
        value = np.load(root / f"{stem}.npy", allow_pickle=False)
        assert value.dtype == np.dtype("<i4"), (stem, value.dtype)
        assert value.shape == shape, (stem, value.shape, shape)
        arrays[stem] = value

    for stem, value in arrays.items():
        assert stem in CATALOG, stem
        spec = CATALOG[stem]
        if spec.cols is not None:
            assert value.shape[-1] == spec.cols, (stem, value.shape, spec.cols)

    # Independent forcing reference: this NPZ was generated immediately
    # before validation by the authoritative mopac2 mopac_feature_probe.py
    # on the same checked-in Cartesian input and exact pinned library.
    with np.load(probe_archive, allow_pickle=False) as probe:
        summary = json.loads(str(probe["summary_json"].item()))
        assert Path(summary["structure"]).resolve() == probe_input.resolve()
        assert Path(summary["library"]).resolve() == Path(
            "/shared/2026Thesis/mopac2/local/lib/libmopac.so.2"
        ).resolve()
        assert summary["atoms"] == n and summary["charge"] == 0

        float_pairs = {
            "mopac_charges_full_precision": "charges",
            "mopac_bond_valencies_full_precision": "bond_valencies",
            "mopac_bond_order": "bond_order",
            "mopac_heat_kcal_mol": "heat_kcal_mol",
            "mopac_dipole_debye": "dipole_debye",
            "mopac_dipole_point_charge_debye":
                "mopac_dipole_point_charge_debye",
            "mopac_dipole_hybridization_debye":
                "mopac_dipole_hybridization_debye",
            "mopac_atom_ao_density": "mopac_atom_ao_density",
            "mopac_atomic_orbital_populations_full_precision":
                "mopac_atomic_orbital_populations",
            "mopac_bond_ao_density": "mopac_bond_ao_density",
            "mopac_atom_electron_population":
                "mopac_atom_electron_population",
            "mopac_atom_s_population": "mopac_atom_s_population",
            "mopac_atom_p_population": "mopac_atom_p_population",
            "mopac_atom_d_population": "mopac_atom_d_population",
            "mopac_lmo_energy_levels": "mopac_lmo_energy_levels",
            "mopac_lmo_occupied_coefficients":
                "mopac_lmo_occupied_coefficients",
            "mopac_lmo_virtual_coefficients":
                "mopac_lmo_virtual_coefficients",
        }
        for production_name, probe_name in float_pairs.items():
            np.testing.assert_allclose(
                arrays[production_name], np.asarray(probe[probe_name]).reshape(
                    arrays[production_name].shape),
                rtol=0.0, atol=1e-12,
                err_msg=f"production/probe mismatch for {production_name}",
            )

        int_pairs = {
            "mopac_bond_index": "bond_index",
            "mopac_bond_atom": "bond_atom",
            "mopac_ao_orbitals_per_atom": "mopac_ao_orbitals_per_atom",
            "mopac_bond_density_pairs": "mopac_bond_density_pairs",
            "mopac_lewis_bond_count": "mopac_lewis_bond_count",
            "mopac_lewis_bond_atoms": "mopac_lewis_bond_atoms",
            "mopac_lmo_occupied_atom_counts":
                "mopac_lmo_occupied_atom_counts",
            "mopac_lmo_occupied_atoms": "mopac_lmo_occupied_atoms",
            "mopac_lmo_virtual_atom_counts":
                "mopac_lmo_virtual_atom_counts",
            "mopac_lmo_virtual_atoms": "mopac_lmo_virtual_atoms",
        }
        for production_name, probe_name in int_pairs.items():
            np.testing.assert_array_equal(
                arrays[production_name], probe[probe_name],
                err_msg=f"production/probe mismatch for {production_name}",
            )
        np.testing.assert_allclose(
            arrays["mopac_bond_orders_full_precision"],
            probe["bond_orders_unique"],
            rtol=0.0, atol=1e-12,
        )
        probe_ao_populations = probe["mopac_atomic_orbital_populations"]
        for atom, live_width in enumerate(
                arrays["mopac_ao_orbitals_per_atom"]):
            live_width = int(live_width)
            np.testing.assert_allclose(
                arrays["mopac_atomic_orbital_populations"][
                    atom, :live_width],
                [legacy_printed(value, 5)
                 for value in probe_ao_populations[atom, :live_width]],
                rtol=0.0, atol=0.0,
            )

    expected_full_precision_charges = np.array([
        -0.6481523853995341,
        0.3240182106742622,
        0.3241341747298975,
    ])
    expected_full_precision_valencies = np.array([
        1.788263447550035,
        0.8950121991514495,
        0.8949370367721683,
    ])
    expected_pairs = np.array([[0, 1], [0, 2]], dtype=np.int32)
    expected_full_precision_orders = np.array([
        0.8941693043878478,
        0.8940941420085942,
    ])
    np.testing.assert_allclose(
        arrays["mopac_charges_full_precision"],
        expected_full_precision_charges, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(
        arrays["mopac_bond_valencies_full_precision"],
        expected_full_precision_valencies, rtol=0.0, atol=1e-12)
    np.testing.assert_array_equal(
        arrays["mopac_bond_orders_full_precision"][:, :2], expected_pairs)
    np.testing.assert_allclose(
        arrays["mopac_bond_orders_full_precision"][:, 2],
        expected_full_precision_orders, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(arrays["mopac_heat_kcal_mol"],
                               [-57.79011922141217], rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(
        arrays["mopac_dipole_debye"],
        [1.313778525459614, 1.6987678419003842, 5.064050067329224e-05],
        rtol=0.0, atol=1e-12,
    )
    np.testing.assert_allclose(
        arrays["mopac_dipole_point_charge_debye"] +
        arrays["mopac_dipole_hybridization_debye"],
        arrays["mopac_dipole_debye"], rtol=0.0, atol=1e-12,
    )

    # Every retained compatibility projection reproduces the numeric value
    # that the old parser recovered from its fixed-precision text field.
    expected_compatibility_charges = np.array([
        legacy_printed(value, 6)
        for value in expected_full_precision_charges
    ])
    np.testing.assert_array_equal(arrays["mopac_charges"],
                                  expected_compatibility_charges)
    expected_compatibility_s = np.array([
        legacy_printed(value, 5)
        for value in arrays["mopac_atom_s_population"]
    ])
    expected_compatibility_p = np.array([
        legacy_printed(value, 5)
        for value in arrays["mopac_atom_p_population"]
    ])
    np.testing.assert_allclose(
        arrays["mopac_scalars"][:, :3],
        np.column_stack((expected_compatibility_charges,
                         expected_compatibility_s,
                         expected_compatibility_p)),
        rtol=0.0, atol=0.0,
    )
    expected_compatibility_global = np.array([
        legacy_printed(arrays["mopac_heat_kcal_mol"][0], 5),
        *(legacy_printed(value, 3)
          for value in arrays["mopac_dipole_debye"]),
    ])
    np.testing.assert_allclose(
        arrays["mopac_global"], expected_compatibility_global,
        rtol=0.0, atol=0.0,
    )

    bond_index = arrays["mopac_bond_index"]
    bond_atom = arrays["mopac_bond_atom"]
    bond_order = arrays["mopac_bond_order"]
    assert bond_index[0] == 0 and bond_index[-1] == entries
    assert np.all(np.diff(bond_index) >= 0)
    raw_pairs: dict[tuple[int, int], list[float]] = {}
    diagonal = np.empty(n)
    for atom in range(n):
        for entry in range(int(bond_index[atom]), int(bond_index[atom + 1])):
            other = int(bond_atom[entry])
            if atom == other:
                diagonal[atom] = bond_order[entry]
            else:
                raw_pairs.setdefault(tuple(sorted((atom, other))), []).append(
                    float(bond_order[entry]))
    np.testing.assert_allclose(diagonal, expected_full_precision_valencies,
                               rtol=0.0, atol=1e-12)
    assert set(raw_pairs) == {(0, 1), (0, 2)}
    assert all(len(values) == 2 for values in raw_pairs.values())
    for values in raw_pairs.values():
        np.testing.assert_allclose(values, values[0], rtol=0.0, atol=1e-12)

    compact_rows, compact_directed_sums = legacy_compact_bond_rows(
        bond_index, bond_atom, bond_order)
    compact_pair_rows = unordered_pair_projection(compact_rows)
    np.testing.assert_allclose(
        rows_sorted_by_atoms(arrays["mopac_bond_orders"]),
        compact_pair_rows, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(arrays["mopac_scalars"][:, 3],
                               compact_directed_sums,
                               rtol=0.0, atol=0.0)

    np.testing.assert_array_equal(arrays["mopac_bond_density_pairs"], expected_pairs)
    np.testing.assert_array_equal(
        arrays["mopac_bond_orders_full_precision"][:, :2], expected_pairs)
    np.testing.assert_allclose(
        arrays["mopac_bond_orders_full_precision"][:, 2],
        expected_full_precision_orders,
        rtol=0.0, atol=1e-12,
    )
    reconstructed = np.sum(arrays["mopac_bond_ao_density"] ** 2, axis=(1, 2))
    np.testing.assert_allclose(
        reconstructed,
        arrays["mopac_bond_orders_full_precision"][:, 2],
        rtol=0.0, atol=1e-12)

    unique_rows = arrays["mopac_bond_orders_unique"]
    np.testing.assert_array_equal(unique_rows[:, :2], expected_pairs)
    expected_unique_orders = np.array([
        [max(legacy_printed(value, 3) for value in raw_pairs[tuple(pair)]),
         np.mean([legacy_printed(value, 3)
                  for value in raw_pairs[tuple(pair)]])]
        for pair in expected_pairs
    ])
    np.testing.assert_allclose(unique_rows[:, 2:4], expected_unique_orders,
                               rtol=0.0, atol=0.0)
    np.testing.assert_array_equal(unique_rows[:, 7], [0.0, 1.0])

    topology_rows = arrays["mopac_topology_bond_orders_full"]
    np.testing.assert_array_equal(topology_rows[:, 0], np.arange(topology))
    np.testing.assert_array_equal(topology_rows[:, 1:3], expected_pairs)
    np.testing.assert_allclose(topology_rows[:, 3],
                               expected_unique_orders[:, 0],
                               rtol=0.0, atol=0.0)
    np.testing.assert_array_equal(topology_rows[:, 4], 1.0)
    np.testing.assert_array_equal(topology_rows[:, 5], np.arange(api_unique))
    np.testing.assert_array_equal(topology_rows[:, 6], 0.0)

    topology_index = {
        tuple(sorted(map(int, row[1:3]))): int(row[0])
        for row in topology_rows
    }
    neighbor_rows = arrays["mopac_bond_neighbors"]
    expected_neighbor_rows = np.asarray([
        (atom, other, order,
         topology_index[tuple(sorted((int(atom), int(other))))])
        for atom_a, atom_b, order in compact_pair_rows
        for atom, other in ((atom_a, atom_b), (atom_b, atom_a))
    ], dtype=np.float64)
    np.testing.assert_allclose(
        rows_sorted_by_atoms(neighbor_rows),
        rows_sorted_by_atoms(expected_neighbor_rows),
        rtol=0.0, atol=0.0)
    for atom in range(n):
        atom_rows = neighbor_rows[neighbor_rows[:, 0] == atom]
        assert np.all(np.diff(atom_rows[:, 2]) <= 0.0)

    ao_widths = arrays["mopac_ao_orbitals_per_atom"]
    np.testing.assert_array_equal(ao_widths, [4, 1, 1])
    assert arrays["mopac_ao_max_orbitals"][0] == width
    atom_density = arrays["mopac_atom_ao_density"]
    ao_pop = arrays["mopac_atomic_orbital_populations"]
    ao_pop_full = arrays[
        "mopac_atomic_orbital_populations_full_precision"]
    diagonal_pop = np.diagonal(atom_density, axis1=1, axis2=2)
    np.testing.assert_allclose(ao_pop_full, diagonal_pop,
                               rtol=0.0, atol=0.0)
    for atom, live_width in enumerate(ao_widths):
        live_width = int(live_width)
        expected_printed_ao = np.array([
            legacy_printed(value, 5)
            for value in diagonal_pop[atom, :live_width]
        ])
        np.testing.assert_allclose(ao_pop[atom, :live_width],
                                   expected_printed_ao,
                                   rtol=0.0, atol=0.0)
        assert np.isnan(ao_pop[atom, live_width:]).all()
    np.testing.assert_allclose(arrays["mopac_atom_electron_population"],
                               diagonal_pop.sum(axis=1), rtol=0.0, atol=1e-15)
    np.testing.assert_allclose(arrays["mopac_atom_s_population"],
                               diagonal_pop[:, 0], rtol=0.0, atol=0.0)
    np.testing.assert_allclose(arrays["mopac_atom_p_population"],
                               diagonal_pop[:, 1:4].sum(axis=1), rtol=0.0, atol=1e-15)
    np.testing.assert_allclose(arrays["mopac_atom_d_population"], 0.0,
                               rtol=0.0, atol=0.0)
    compatibility_shells = arrays["mopac_atomic_orbital_population_totals"]
    for atom, live_width in enumerate(ao_widths):
        live_width = int(live_width)
        np.testing.assert_allclose(compatibility_shells[atom, 0],
                                   ao_pop[atom, 0],
                                   rtol=0.0, atol=0.0)
        if live_width >= 4:
            np.testing.assert_allclose(
                compatibility_shells[atom, 1],
                np.sum(ao_pop[atom, 1:4]),
                rtol=0.0, atol=0.0)
        else:
            assert np.isnan(compatibility_shells[atom, 1])
        if live_width >= 9:
            np.testing.assert_allclose(
                compatibility_shells[atom, 2],
                np.sum(ao_pop[atom, 4:9]),
                rtol=0.0, atol=0.0)
        else:
            assert np.isnan(compatibility_shells[atom, 2])

    atom_pop = arrays["mopac_atom_populations"]
    expected_compatibility_electrons = np.array([
        legacy_printed(value, 4)
        for value in arrays["mopac_atom_electron_population"]
    ])
    np.testing.assert_allclose(atom_pop[:, :3], np.column_stack((
        expected_compatibility_charges,
        expected_compatibility_electrons,
        expected_compatibility_s,
    )), rtol=0.0, atol=0.0)
    np.testing.assert_allclose(atom_pop[0, 3], expected_compatibility_p[0],
                               rtol=0.0, atol=0.0)
    assert np.isnan(atom_pop[1:, 3]).all()
    assert np.isnan(atom_pop[:, 4]).all()
    expected_compatibility_valencies = np.array([
        legacy_printed(value, 3)
        for value in expected_full_precision_valencies
    ])
    np.testing.assert_allclose(arrays["mopac_bond_valencies"],
                               expected_compatibility_valencies,
                               rtol=0.0, atol=0.0)
    np.testing.assert_allclose(atom_pop[:, 10],
                               expected_compatibility_valencies,
                               rtol=0.0, atol=0.0)
    np.testing.assert_allclose(atom_pop[:, 11], compact_directed_sums,
                               rtol=0.0, atol=0.0)
    assert np.isnan(atom_pop[:, 5:10]).all()
    assert np.isnan(arrays["mopac_bond_orders_unique"][:, 4:7]).all()
    assert np.isnan(arrays["mopac_topology_bond_orders_full"][:, 7]).all()

    lewis_count = arrays["mopac_lewis_bond_count"]
    lewis_atoms = arrays["mopac_lewis_bond_atoms"]
    for atom, count in enumerate(lewis_count):
        count = int(count)
        assert np.all((0 <= lewis_atoms[atom, :count]) &
                      (lewis_atoms[atom, :count] < n))
        assert np.all(lewis_atoms[atom, count:] == -1)

    def check_lmos(kind: str, lmo_count: int) -> None:
        counts = arrays[f"mopac_lmo_{kind}_atom_counts"]
        packed_atoms = arrays[f"mopac_lmo_{kind}_atoms"]
        packed_coefficients = arrays[f"mopac_lmo_{kind}_coefficients"]
        atom_offsets = arrays[f"mopac_lmo_{kind}_atom_offsets_native"]
        coefficient_offsets = arrays[
            f"mopac_lmo_{kind}_coefficient_offsets_native"]
        native_atoms = arrays[f"mopac_lmo_{kind}_atom_storage_native"]
        native_coefficients = arrays[
            f"mopac_lmo_{kind}_coefficient_storage_native"]
        assert len(counts) == lmo_count
        atom_cursor = 0
        coefficient_cursor = 0
        for count, atom_offset, coefficient_offset in zip(
                counts, atom_offsets, coefficient_offsets, strict=True):
            count = int(count)
            atom_offset = int(atom_offset)
            coefficient_offset = int(coefficient_offset)
            atoms = native_atoms[atom_offset:atom_offset + count] - 1
            np.testing.assert_array_equal(
                packed_atoms[atom_cursor:atom_cursor + count], atoms)
            coefficient_count = int(np.sum(ao_widths[atoms]))
            np.testing.assert_allclose(
                packed_coefficients[
                    coefficient_cursor:coefficient_cursor + coefficient_count],
                native_coefficients[
                    coefficient_offset:coefficient_offset + coefficient_count],
                rtol=0.0, atol=0.0,
            )
            atom_cursor += count
            coefficient_cursor += coefficient_count
        assert atom_cursor == len(packed_atoms)
        assert coefficient_cursor == len(packed_coefficients)

    check_lmos("occupied", occupied)
    check_lmos("virtual", virtual)


if __name__ == "__main__":
    main()

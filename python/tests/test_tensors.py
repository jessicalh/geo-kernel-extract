import inspect

import numpy as np
import torch
from e3nn.o3 import Irreps, spherical_harmonics

from nmr_extract import (
    CATALOG,
    EFGTensor,
    MopacAtomicOrbitalPopulations,
    MopacAtomicOrbitalPopulationTotals,
    PerBondCategoryT2,
    ShieldingTensor,
    project_full9_to_e3nn,
    project_t2_to_e3nn,
)
from nmr_extract._protein import WaterFieldGroup
from nmr_extract._tensors import DeltaAPBS


def test_project_t2_to_e3nn_basis_vectors_and_norms():
    expected = np.array([
        [0.0, 1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, -0.5, 0.0, np.sqrt(3.0) / 2.0],
        [1.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, -np.sqrt(3.0) / 2.0, 0.0, -0.5],
    ])
    got = project_t2_to_e3nn(np.eye(5))
    np.testing.assert_allclose(got, expected)
    np.testing.assert_allclose(np.linalg.norm(got, axis=-1), 1.0)


def test_project_t2_to_e3nn_pinned_cartesian_quadratics():
    vectors = [
        (
            np.array([0, 0, -1/np.sqrt(6), 0, 1/np.sqrt(2)]),
            np.array([0, 0, -1/np.sqrt(6), 0, -1/np.sqrt(2)]),
        ),
        (
            np.array([0, 0, -1/np.sqrt(6), 0, -1/np.sqrt(2)]),
            np.array([0, 0, 2/np.sqrt(6), 0, 0]),
        ),
        (
            np.array([0, 0, 2/np.sqrt(6), 0, 0]),
            np.array([0, 0, -1/np.sqrt(6), 0, 1/np.sqrt(2)]),
        ),
    ]
    for project, e3nn in vectors:
        converted = project_t2_to_e3nn(project)
        np.testing.assert_allclose(converted, e3nn, atol=1e-15)
        np.testing.assert_allclose(np.linalg.norm(converted), np.linalg.norm(project))


def test_project_t2_to_e3nn_matches_live_e3nn_norm_spherical_harmonics():
    rng = np.random.default_rng(20260709)
    unit = rng.normal(size=(12, 3))
    unit /= np.linalg.norm(unit, axis=1, keepdims=True)
    x, y, z = unit.T

    project_t2 = np.column_stack([
        np.sqrt(2.0) * x * y,
        np.sqrt(2.0) * y * z,
        np.sqrt(3.0 / 2.0) * (z * z - 1.0 / 3.0),
        np.sqrt(2.0) * x * z,
        (x * x - y * y) / np.sqrt(2.0),
    ])
    e3nn_2e = spherical_harmonics(
        2,
        torch.as_tensor(unit, dtype=torch.float64),
        normalize=True,
        normalization="norm",
    ).numpy()

    expected = np.sqrt(2.0 / 3.0) * e3nn_2e
    np.testing.assert_allclose(project_t2_to_e3nn(project_t2), expected, atol=1e-14)

    full9 = np.zeros((unit.shape[0], 9), dtype=np.float64)
    full9[:, :4] = rng.normal(size=(unit.shape[0], 4))
    full9[:, 4:9] = project_t2
    converted = project_full9_to_e3nn(full9)
    np.testing.assert_allclose(converted[:, :4], full9[:, :4])
    np.testing.assert_allclose(converted[:, 4:9], expected, atol=1e-14)


def test_project_full9_to_e3nn_preserves_t0_t1_and_does_not_mutate():
    full = np.arange(18, dtype=np.float64).reshape(2, 9)
    original = full.copy()
    converted = project_full9_to_e3nn(full)
    np.testing.assert_array_equal(full, original)
    np.testing.assert_array_equal(converted[:, :4], full[:, :4])
    np.testing.assert_allclose(converted[:, 4:9], project_t2_to_e3nn(full[:, 4:9]))


def test_raw_wrappers_convert_explicitly_to_e3nn():
    st = ShieldingTensor(np.zeros((2, 9), dtype=np.float64))
    assert not hasattr(st, "irreps")
    assert st.tensor_basis == "project_native_full9_spherical_tensor_v1"
    assert st.to_e3nn().irreps == Irreps("1x0e + 1x1e + 1x2e")
    assert st.to_e3nn_T2().irreps == Irreps("1x2e")

    efg = EFGTensor(np.zeros((2, 5), dtype=np.float64))
    assert not hasattr(efg, "irreps")
    assert efg.tensor_basis == "project_native_t2_isometric_real_tesseral_v1"
    assert efg.to_e3nn().irreps == Irreps("1x2e")

    block = PerBondCategoryT2(np.zeros((2, 25), dtype=np.float64))
    assert not hasattr(block, "irreps")
    assert block.component_order == "T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"
    assert block.to_e3nn().irreps == Irreps("5x2e")


def test_mopac_ao_population_wrappers_and_totals_metadata():
    raw = MopacAtomicOrbitalPopulations(np.array([
        [1.0, 0.1, 0.2, 0.3, 0.01, 0.02, 0.03, 0.04, 0.05],
    ]))
    totals = MopacAtomicOrbitalPopulationTotals(np.array([[1.0, 0.6, 0.15]]))
    assert raw.model_scalar_source is False
    assert raw.diagnostic_frame_dependent is True
    assert totals.model_scalar_source is True
    assert totals.diagnostic_frame_dependent is False
    np.testing.assert_array_equal(raw.data[0], [
        1.0, 0.1, 0.2, 0.3, 0.01, 0.02, 0.03, 0.04, 0.05,
    ])
    np.testing.assert_allclose(totals.data[0], [1.0, 0.6, 0.15])


def test_delta_apbs_t2_accessor_and_catalog_metadata():
    data = np.arange(24, dtype=np.float64).reshape(2, 12)
    apbs = DeltaAPBS(data)
    np.testing.assert_array_equal(apbs.delta_E.data, data[:, :3])
    np.testing.assert_array_equal(apbs.delta_efg.data, data[:, 3:12])
    np.testing.assert_array_equal(apbs.delta_efg_full9.data, data[:, 3:12])
    np.testing.assert_array_equal(apbs.delta_efg_t2.data, data[:, 7:12])
    np.testing.assert_array_equal(apbs.delta_efg_t2.data, apbs.delta_efg.T2)
    spec = CATALOG["delta_apbs"]
    assert spec.structural_zero_components == "T0,T1_x,T1_y,T1_z"
    assert "T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2" in spec.tensor_component_order


def test_catalog_tensor_metadata_and_water_doc_shape():
    assert CATALOG["water_efg"].cols == 5
    assert CATALOG["water_efg_first"].cols == 5
    assert "(N, 9)" not in inspect.getsource(WaterFieldGroup)

    assert CATALOG["coulomb_efg"].cols == 9
    assert CATALOG["coulomb_efg"].structural_zero_components == "T0,T1_x,T1_y,T1_z"
    assert CATALOG["eeq_coulomb_efg"].cols == 9
    assert CATALOG["eeq_coulomb_efg"].structural_zero_components == \
        "T0,T1_x,T1_y,T1_z"
    assert CATALOG["coulomb_efg_t2"].cols == 5
    assert CATALOG["coulomb_efg_t2"].wrapper is EFGTensor
    assert CATALOG["mopac_coulomb_efg"].cols == 9
    assert CATALOG["mopac_atomic_orbital_population_totals"].wrapper is \
        MopacAtomicOrbitalPopulationTotals

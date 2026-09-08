"""Worked nmr_extract-equivalent Biot-Savart calculation through Reader REST."""

from __future__ import annotations

import math
import os

import pytest


pytestmark = pytest.mark.skipif(
    "1p9j" not in os.environ.get("H5READER_REST_FIXTURE", "").lower(),
    reason="fixed producer comparison uses the complete 1P9J fixture",
)


def _assert_vector_close(actual, expected, *, absolute_tolerance):
    assert len(actual) == len(expected)
    for observed, reference in zip(actual, expected, strict=True):
        assert math.isclose(
            observed, reference, rel_tol=1.0e-12, abs_tol=absolute_tolerance
        )


def test_one_ring_matches_nmr_extract_and_exposes_the_segment_sum(rest):
    response = rest.client.post(
        "/api/ring/biot_savart",
        json={"frame": 0, "ring": 1, "atom": 348, "include_segments": True},
    )
    assert response.status_code == 200, response.text
    payload = response.json()

    assert payload["kind"] == "biot_savart_ring_current"
    assert payload["frame"] == 0
    assert payload["ring"] == 1
    assert payload["target"]["atom"] == 348
    assert payload["ring_identity"]["type_name"] == "PHE"
    assert payload["sample"]["status"] == "evaluated"
    assert payload["sample"]["intensity_nA_per_T"] == 1.0

    geometry = payload["ring_geometry"]
    _assert_vector_close(
        geometry["center_A"],
        [81.8303108215332, 62.7879524230957, 33.71971845626831],
        absolute_tolerance=1.0e-12,
    )
    _assert_vector_close(
        geometry["normal"],
        [-0.678949302748889, -0.2018531824122542, -0.7058917318164616],
        absolute_tolerance=1.0e-12,
    )
    assert math.isclose(
        geometry["mean_radius_A"], 1.3926680043082793,
        rel_tol=1.0e-12, abs_tol=1.0e-12,
    )
    assert geometry["lobe_offset_A"] == 0.64

    expected_field = [
        1.1482225029724714e-08,
        -1.4160182596153648e-08,
        -2.215760591963646e-08,
    ]
    _assert_vector_close(
        payload["sample"]["induced_B_T_per_T"],
        expected_field,
        absolute_tolerance=1.0e-20,
    )
    for row in range(3):
        for column in range(3):
            expected = -expected_field[row] * geometry["normal"][column] * 1.0e6
            assert math.isclose(
                payload["sample"]["shielding_cartesian_ppm"][row][column],
                expected,
                rel_tol=1.0e-12,
                abs_tol=1.0e-12,
            )
    _assert_vector_close(
        payload["sample"]["shielding_spherical"]["T1"],
        [
            -0.002761486273060498,
            -0.011574549400530265,
            0.005965884881939523,
        ],
        absolute_tolerance=1.0e-12,
    )
    assert math.isclose(
        payload["sample"]["shielding_spherical"]["T0"],
        -0.0035677666860512193,
        rel_tol=1.0e-12,
        abs_tol=1.0e-12,
    )
    _assert_vector_close(
        payload["sample"]["shielding_spherical"]["T2"],
        [
            -0.005159279072948812,
            -0.010230519258260715,
            -0.014786472364342114,
            -0.004906390069361093,
            0.007533605165426149,
        ],
        absolute_tolerance=1.0e-12,
    )

    segments = payload["sample"]["segments"]
    assert len(segments) == 12
    summed_field = [
        sum(segment["induced_B_T_per_T"][axis] for segment in segments)
        for axis in range(3)
    ]
    _assert_vector_close(summed_field, expected_field, absolute_tolerance=1.0e-20)

    scaled_response = rest.client.post(
        "/api/ring/biot_savart",
        json={
            "frame": 0,
            "ring": 1,
            "atom": 348,
            "intensity_nA_per_T": -12.0,
            "include_segments": False,
        },
    )
    assert scaled_response.status_code == 200, scaled_response.text
    scaled_sample = scaled_response.json()["sample"]
    assert scaled_sample["intensity_nA_per_T"] == -12.0
    _assert_vector_close(
        scaled_sample["induced_B_T_per_T"],
        [-12.0 * component for component in expected_field],
        absolute_tolerance=1.0e-19,
    )
    assert math.isclose(
        scaled_sample["shielding_spherical"]["T0"],
        -12.0 * payload["sample"]["shielding_spherical"]["T0"],
        rel_tol=1.0e-12,
        abs_tol=1.0e-12,
    )


def test_point_inside_mean_ring_radius_is_not_mistaken_for_a_singularity(rest):
    reference = rest.client.post(
        "/api/ring/biot_savart",
        json={"frame": 0, "ring": 1, "atom": 348, "include_segments": False},
    )
    assert reference.status_code == 200, reference.text
    geometry = reference.json()["ring_geometry"]
    point = [
        geometry["center_A"][axis] + 0.5 * geometry["normal"][axis]
        for axis in range(3)
    ]

    response = rest.client.post(
        "/api/ring/biot_savart",
        json={"frame": 0, "ring": 1, "point_A": point, "include_segments": False},
    )
    assert response.status_code == 200, response.text
    sample = response.json()["sample"]
    assert sample["status"] == "evaluated"
    assert sample["distance_to_ring_center_A"] < geometry["mean_radius_A"]
    assert any(abs(component) > 0.0 for component in sample["induced_B_T_per_T"])

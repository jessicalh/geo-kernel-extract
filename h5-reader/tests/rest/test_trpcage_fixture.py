"""Trp-cage dense-DFT fixture smoke for the shipped July calcset."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import math
import os
import time

import httpx
import pytest


def _fixture_is_trpcage() -> bool:
    fixture = os.environ.get("H5READER_REST_FIXTURE", "").lower()
    return "trpcage" in fixture or "1l2y" in fixture or "tc5b" in fixture


pytestmark = pytest.mark.skipif(
    not _fixture_is_trpcage(),
    reason="set H5READER_REST_FIXTURE to the Trp-cage .LGS to run",
)


def test_trpcage_dense_dft_fixture_loads_and_probes_csa(rest):
    state = rest.client.get("/ui/state").json()
    assert state["protein"] == "1L2Y_5292"
    assert state["frames"] == 1501

    atoms = rest.client.get("/protein/atoms").json()
    assert atoms["count"] == 304

    expected_sigma_iso = {
        0: 228.92000000000016,
        750: 229.8833333333333,
        1500: 225.06133333333347,
    }
    for frame, expected in expected_sigma_iso.items():
        response = rest.client.post("/frame/set", json={"frame": frame})
        assert response.status_code == 204, response.text
        csa = rest.client.get("/csa", params={"atom": 0, "frame": frame}).json()
        assert csa["atom"] == 0
        assert csa["frame"] == frame
        assert csa["dft_present"] is True
        assert csa["valid"] is True
        assert math.isfinite(csa["sigma_iso"])
        assert abs(csa["sigma_iso"] - expected) < 1e-6
        assert math.isfinite(csa["span"])
        assert math.isfinite(csa["eta"])


def test_trpcage_candidate_a_ring_system_cloud_matches_frozen_analysis(rest):
    surface_response = rest.client.post(
        "/resthero/butterfly",
        json={
            "rings": [1, 2],
            "frame": 750,
            "mode": "candidate_a",
            "dim": 32,
            "threshold_ppm": 0.5,
            "opacity": 0.12,
            "extent": 0.0,
            "show_source_loops": True,
            "source_loop_tube_radius_A": 0.018,
        },
    )
    assert surface_response.status_code == 200, surface_response.text
    surface = surface_response.json()
    assert surface["rings"] == [1, 2]
    assert surface["mode"] == "candidate_a"
    assert surface["frame"] == 750
    assert surface["extent"] == 0.0
    assert surface["min_T0"] < -0.5
    assert surface["max_T0"] > 0.5
    assert surface["shielded_cells"] > 0
    assert surface["deshielded_cells"] > 0
    assert surface["source_loops_visible"] is True
    assert surface["source_loop_count"] == 4
    assert surface["source_loop_actor_count"] == 2
    assert surface["surface_quantity"] == "candidate_shielding_T0_ppm"
    assert surface["blue_surface_T0_ppm"] == 0.5
    assert surface["blue_surface_predicted_shift_delta_ppm"] == -0.5
    assert surface["red_surface_T0_ppm"] == -0.5
    assert surface["red_surface_predicted_shift_delta_ppm"] == 0.5
    assert len(surface["circular_sources"]) == 2
    assert [source["ring"] for source in surface["circular_sources"]] == [1, 2]
    assert [source["loop_radius_A"] for source in surface["circular_sources"]] == [
        1.39,
        1.182,
    ]
    assert [source["loop_offset_A"] for source in surface["circular_sources"]] == [
        0.64,
        0.64,
    ]
    for source in surface["circular_sources"]:
        assert source["plane_rms_A"] < 0.05
        assert source["total_current_nA_per_T"] < 0.0
        assert math.isclose(
            source["current_per_loop_nA_per_T"],
            0.5 * source["total_current_nA_per_T"],
            abs_tol=1e-12,
        )

    response = rest.client.post(
        "/resthero/ring_system_cloud",
        json={
            "atoms": [180, 275, 272, 273],
            "rings": [1, 2],
            "reference_ring": 1,
            "reference_frame": 750,
            "start_frame": 0,
            "end_frame": 1500,
            "step": 1,
            "max_frames": 1501,
            "show_points": True,
        },
    )
    assert response.status_code == 200, response.text
    payload = response.json()
    assert payload["model"] == "candidate_a_circular_two_loop"
    assert payload["coordinate_space"] == "source_ring_local"
    assert payload["ring_plane_fit"] == "least_squares_svd"
    assert payload["ring_plane_sign"] == "aligned_to_winding_normal"
    assert payload["color_quantity"] == "predicted_shift_delta_ppm"
    assert payload["rings"] == [1, 2]
    assert payload["reference_ring"] == 1
    assert payload["frame_count"] == 1501
    assert payload["sample_count"] == 6004
    assert payload["max_local_roundtrip_delta_A"] < 1e-12
    assert payload["actor_count"] == 1
    assert payload["color_scale_mode"] == "auto"
    assert payload["requested_color_scale"] == 0.0
    assert payload["color_scale"] > 5.0

    expected = {
        180: (2.8925452931342166, 0.08080587623455238, 4.944671935394664, 1501),
        272: (0.7806484364587049, 0.07299702946563898, 1.484135500342942, 1501),
        273: (2.201855426620336, -0.1974103051381922, 5.401362752847593, 1500),
        275: (2.4972384196286392, 0.058734549850106454, 4.521378642194338, 1501),
    }
    observed = {item["atom"]: item for item in payload["atoms"]}
    assert observed.keys() == expected.keys()
    expected_names = {
        180: ("GLY", 11, "HA2"),
        272: ("PRO", 18, "HB2"),
        273: ("PRO", 18, "HB3"),
        275: ("PRO", 18, "HA"),
    }
    for atom, (mean, minimum, maximum, upfield_frames) in expected.items():
        item = observed[atom]
        residue, residue_number, atom_name = expected_names[atom]
        assert item["identity"]["residue_label_iupac"] == residue
        assert item["identity"]["residue_number"] == residue_number
        assert item["identity"]["atom_label_iupac"] == atom_name
        assert item["frames"] == 1501
        assert math.isclose(item["mean_candidate_shielding_T0_ppm"], mean, abs_tol=1e-9)
        assert math.isclose(
            item["min_candidate_shielding_T0_ppm"], minimum, abs_tol=1e-9
        )
        assert math.isclose(
            item["max_candidate_shielding_T0_ppm"], maximum, abs_tol=1e-9
        )
        assert item["upfield_frames"] == upfield_frames

    coordinate_response = rest.client.post(
        "/resthero/ring_system_cloud",
        json={
            "atoms": [180],
            "rings": [1, 2],
            "reference_ring": 1,
            "reference_frame": 750,
            "start_frame": 0,
            "end_frame": 0,
            "echo_samples": True,
            "show_points": True,
        },
    )
    assert coordinate_response.status_code == 200, coordinate_response.text
    coordinate_payload = coordinate_response.json()
    assert coordinate_payload["actor_count"] == 1
    local_position = coordinate_payload["samples"][0]["ring_local_position"]
    expected_local_position = (
        -0.06601516077142318,
        0.5437654448108414,
        3.6338207545275916,
    )
    for observed_value, expected_value in zip(local_position, expected_local_position):
        assert math.isclose(observed_value, expected_value, abs_tol=1e-9)

    cleared = rest.client.post("/resthero/clear")
    assert cleared.status_code == 204, cleared.text


def test_trpcage_candidate_a_tensor_comparison_matches_frozen_analysis(rest):
    frame_response = rest.client.post("/frame/set", json={"frame": 237})
    assert frame_response.status_code == 204, frame_response.text

    request_body = {
        "atoms": [180, 275, 272, 273],
        "rings": [1, 2],
        "scale_A_per_ppm": 0.42,
        "theta_resolution": 24,
        "phi_resolution": 24,
    }

    def run_comparison():
        with httpx.Client(base_url=rest.base_url, timeout=60.0) as client:
            return client.post("/resthero/ring_tensor_compare", json=request_body)

    with ThreadPoolExecutor(max_workers=1) as executor:
        result = executor.submit(run_comparison)
        with httpx.Client(base_url=rest.base_url, timeout=5.0) as observer:
            deadline = time.monotonic() + 5.0
            while True:
                active = observer.post("/resthero/ring_tensor_compare", json={})
                if active.status_code == 409:
                    break
                assert active.status_code == 400, active.text
                assert not result.done(), (
                    "tensor comparison finished before its active state was observed"
                )
                assert time.monotonic() < deadline, (
                    "tensor comparison did not enter its active state"
                )
                time.sleep(0.01)

            health_start = time.monotonic()
            health = observer.get("/health")
            health_elapsed = time.monotonic() - health_start
            assert health.status_code == 200, health.text
            assert health.json()["ok"] is True
            assert health_elapsed < 2.0

        response = result.result(timeout=60.0)

    assert response.status_code == 200, response.text
    payload = response.json()
    assert payload["model"] == "candidate_a_circular_two_loop"
    assert payload["reference"] == "whole_protein_orca_total_shielding"
    assert payload["comparison"] == (
        "residue_N_CA_C_local_symmetric_traceless_time_demeaned"
    )
    assert payload["frame"] == 237
    assert payload["mean_frame_count"] == 1501
    assert payload["sample_count"] == 4
    assert payload["rings"] == [1, 2]
    assert [item["ring"] for item in payload["ring_identities"]] == [1, 2]
    assert [
        item["identity"]["parent_residue_number"] for item in payload["ring_identities"]
    ] == [6, 6]
    assert [item["identity"]["type_name"] for item in payload["ring_identities"]] == [
        "TRP6",
        "TRP5",
    ]
    assert payload["candidate_representation"] == ("solid_signed_quadratic_surface")
    assert payload["orca_representation"] == "wire_signed_quadratic_surface"
    assert payload["minimum_radius_A"] == 0.015
    assert payload["radial_equation"] == ("r(n)=max(minimum_radius,scale*abs(n^T T n))")

    expected = {
        180: (0.8645527514059487, 0.8471206226213741, 2.1371272368444987),
        275: (0.7501685511575085, 0.6838929258507388, 1.3904509002649705),
        272: (-0.06623557859023582, -0.7114257183233731, -0.11714440183333874),
        273: (0.7414491605299572, 0.6043055023170284, 1.7754607137586014),
    }
    observed = {item["atom"]: item for item in payload["atoms"]}
    assert observed.keys() == expected.keys()

    expected_gly11_candidate = [
        [-0.8413822913565538, 0.08300698290198638, -0.8158685895035199],
        [0.08300698290198638, 1.8190865671950482, -0.8455866639973766],
        [-0.8158685895035196, -0.8455866639973775, -0.9777042758384935],
    ]
    observed_gly11_candidate = observed[180]["candidate_local_residual_T2"]
    for observed_row, expected_row in zip(
        observed_gly11_candidate, expected_gly11_candidate
    ):
        assert observed_row == pytest.approx(expected_row, abs=1e-9)

    for atom, (pooled_cosine, live_cosine, rmse_reduction) in expected.items():
        item = observed[atom]
        assert math.isclose(
            item["pooled_local_demeaned_T2_cosine"], pooled_cosine, abs_tol=1e-9
        )
        assert math.isclose(
            item["live_frame_local_T2_cosine"], live_cosine, abs_tol=1e-9
        )
        assert math.isclose(
            item["demeaned_fixed_candidate_rmse_reduction_ppm"],
            rmse_reduction,
            abs_tol=1e-9,
        )

    cleared = rest.client.post("/resthero/clear")
    assert cleared.status_code == 204, cleared.text


def test_trpcage_candidate_a_two_state_tensor_comparison(rest):
    near_frames = [1, 2, 3, 4, 5, 6, 7, 9, 11, 12, 13, 14]
    far_frames = [0, 8, 10, 17, 19, 20, 33, 35, 42, 54, 56, 64]

    incomplete = rest.client.post(
        "/resthero/ring_tensor_compare",
        json={
            "atoms": [273],
            "rings": [1, 2],
            "state_a_frames": near_frames,
        },
    )
    assert incomplete.status_code == 400, incomplete.text

    frame_response = rest.client.post("/frame/set", json={"frame": 755})
    assert frame_response.status_code == 204, frame_response.text
    response = rest.client.post(
        "/resthero/ring_tensor_compare",
        json={
            "atoms": [273],
            "rings": [1, 2],
            "state_a_label": "down_pucker_HB3_near",
            "state_b_label": "up_pucker_HB3_far",
            "state_a_frames": near_frames,
            "state_b_frames": far_frames,
            "normal_reference_ring": 3,
            "theta_resolution": 24,
            "phi_resolution": 24,
        },
    )
    assert response.status_code == 200, response.text
    payload = response.json()
    assert payload["comparison"] == (
        "residue_N_CA_C_local_symmetric_traceless_state_b_minus_state_a"
    )
    assert payload["reference_ring"] == 1
    assert payload["state_a_label"] == "down_pucker_HB3_near"
    assert payload["state_b_label"] == "up_pucker_HB3_far"
    assert payload["state_a_frame_count"] == len(near_frames)
    assert payload["state_b_frame_count"] == len(far_frames)
    assert payload["unassigned_frame_count"] == 1501 - 24
    assert payload["normal_reference_ring"] == 3
    assert payload["reference_ring_normal_averaging"] == "state_a"
    assert payload["reference_ring_normal_frame_count"] == len(near_frames)

    atom = payload["atoms"][0]
    assert atom["atom"] == 273
    assert atom["identity"]["atom_label_iupac"] == "HB3"
    assert len(atom["candidate_state_b_minus_state_a_local_T2"]) == 3
    assert len(atom["orca_state_b_minus_state_a_local_T2"]) == 3
    assert -1.0 <= atom["state_b_minus_state_a_local_T2_cosine"] <= 1.0
    assert 0.0 <= atom["dominant_axis_separation_degrees"] <= 90.0
    assert (
        0.0
        <= atom["candidate_dominant_axis_from_reference_ring_normal_degrees"]
        <= 90.0
    )
    assert 0.0 <= atom["orca_dominant_axis_from_reference_ring_normal_degrees"] <= 90.0
    assert math.isfinite(
        atom["candidate_state_b_minus_state_a_predicted_shift_delta_ppm"]
    )
    assert math.isfinite(atom["orca_state_b_minus_state_a_shift_like_ppm"])
    expected_candidate = [
        [-1.3842560160068536, -1.6438210746872404, 2.3294999614069596],
        [-1.6438210746872404, 1.1377650122491263, 1.519120659097952],
        [2.3294999614069596, 1.519120659097952, 0.2464910037577268],
    ]
    expected_orca = [
        [-5.6436151094178095, -5.7636623352155105, 3.8982311702059866],
        [-5.7636623352155105, 5.855574044532208, 1.324078282624608],
        [3.8982311702059866, 1.324078282624608, -0.2119589351143964],
    ]
    for observed_row, expected_row in zip(
        atom["candidate_state_b_minus_state_a_local_T2"], expected_candidate
    ):
        assert observed_row == pytest.approx(expected_row, abs=1e-9)
    for observed_row, expected_row in zip(
        atom["orca_state_b_minus_state_a_local_T2"], expected_orca
    ):
        assert observed_row == pytest.approx(expected_row, abs=1e-9)
    assert atom[
        "candidate_state_b_minus_state_a_predicted_shift_delta_ppm"
    ] == pytest.approx(1.6526703462058592, abs=1e-9)
    assert atom["orca_state_b_minus_state_a_shift_like_ppm"] == pytest.approx(
        1.8028333333333464, abs=1e-9
    )

    director_response = rest.client.post(
        "/resthero/ring_tensor_compare",
        json={
            "atoms": [273],
            "rings": [1, 2],
            "state_a_frames": near_frames,
            "state_b_frames": far_frames,
            "normal_reference_ring": 3,
            "representation": "director",
            "show_reference_ring_normal": True,
            "director_half_length_A": 1.4,
            "director_radius_A": 0.02,
            "theta_resolution": 24,
            "phi_resolution": 24,
        },
    )
    assert director_response.status_code == 200, director_response.text
    director = director_response.json()
    assert director["representation"] == "director"
    assert director["candidate_representation"] == "dominant_unoriented_director"
    assert director["orca_representation"] == "dominant_unoriented_director"
    assert director["director_definition"] == (
        "eigenvector_of_largest_absolute_eigenvalue"
    )
    assert director["director_length_encodes"] == "nothing_stylistic_only"
    assert director["reference_ring_normal_visible"] is True
    assert director["director_half_length_A"] == 1.4
    assert director["candidate_director_half_length_A"] == 1.4
    assert director["orca_director_half_length_A"] == pytest.approx(1.568)
    assert director["reference_ring_normal_half_length_A"] == pytest.approx(1.736)
    assert director["director_radius_A"] == 0.02

    director_atom = director["atoms"][0]
    candidate_axis = director_atom["candidate_dominant_director_display"]
    orca_axis = director_atom["orca_dominant_director_display"]
    ring_normal = director_atom["reference_ring_mean_normal_display"]
    for vector in (candidate_axis, orca_axis, ring_normal):
        assert math.sqrt(
            sum(component * component for component in vector)
        ) == pytest.approx(1.0, abs=1e-12)

    def acute_angle_degrees(first, second):
        cosine = abs(sum(a * b for a, b in zip(first, second)))
        return math.degrees(math.acos(min(1.0, max(-1.0, cosine))))

    assert acute_angle_degrees(candidate_axis, orca_axis) == pytest.approx(
        director_atom["dominant_axis_separation_degrees"], abs=1e-10
    )
    assert acute_angle_degrees(candidate_axis, ring_normal) == pytest.approx(
        director_atom["candidate_dominant_axis_from_reference_ring_normal_degrees"],
        abs=1e-10,
    )
    assert acute_angle_degrees(orca_axis, ring_normal) == pytest.approx(
        director_atom["orca_dominant_axis_from_reference_ring_normal_degrees"],
        abs=1e-10,
    )

    cleared = rest.client.post("/resthero/clear")
    assert cleared.status_code == 204, cleared.text


def test_trpcage_shutdown_during_active_tensor_comparison_is_clean(rest):
    request_body = {
        "atoms": [180, 275, 272, 273],
        "rings": [1, 2],
        "theta_resolution": 24,
        "phi_resolution": 24,
    }

    def run_comparison():
        with httpx.Client(base_url=rest.base_url, timeout=60.0) as client:
            return client.post("/resthero/ring_tensor_compare", json=request_body)

    with ThreadPoolExecutor(max_workers=1) as executor:
        result = executor.submit(run_comparison)
        with httpx.Client(base_url=rest.base_url, timeout=5.0) as observer:
            deadline = time.monotonic() + 5.0
            while True:
                active = observer.post("/resthero/ring_tensor_compare", json={})
                if active.status_code == 409:
                    break
                assert active.status_code == 400, active.text
                assert not result.done()
                assert time.monotonic() < deadline
                time.sleep(0.01)

            with httpx.Client(base_url=rest.base_url, timeout=5.0) as shutdown_client:
                shutdown = shutdown_client.post("/shutdown")
                assert shutdown.status_code == 204, shutdown.text

        rest.process.wait(timeout=10.0)
        assert rest.process.returncode == 0
        try:
            comparison = result.result(timeout=5.0)
        except httpx.HTTPError:
            pass
        else:
            assert comparison.status_code == 200, comparison.text

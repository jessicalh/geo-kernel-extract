"""REST contract for the ring-current weak-signal receiver."""

from __future__ import annotations


def test_ring_current_face_collar_requires_filter_or_explicit_full_scan(rest):
    response = rest.client.post("/ring/current_face_collar", json={})
    assert response.status_code == 400
    assert "full_scan" in response.json()["error"]


def test_ring_current_face_collar_fits_orca_to_expected_relationship(rest):
    response = rest.client.post(
        "/ring/current_face_collar",
        json={
            "atom": 568,
            "ring": 4,
            "start_frame": 500,
            "end_frame": 580,
            "min_samples": 6,
            "max_entries": 1,
            "include_samples": True,
        },
    )
    assert response.status_code == 200, response.text
    payload = response.json()

    assert payload["kind"] == "ring_current_face_collar"
    assert "expected_relationship_value" in payload["receiver"]
    assert "ORCA_component = intercept + scale" in payload["receiver"]["fit_model"]
    assert "finite Johnson-Bovey" in payload["receiver"]["biot_savart_relationship"]
    assert payload["summary"]["complete"]
    assert not payload["summary"]["truncated_by_max_entries"]
    assert payload["summary"]["entry_count"] == 1
    assert payload["summary"]["paths_considered"] == 1
    assert payload["summary"]["paths_rejected_for_hard_crossing"] == 0
    assert payload["summary"]["paths_rejected_for_weak_lobes"] == 0
    assert payload["options"]["min_samples_per_lobe"] == 3
    assert payload["options"]["min_expected_relationship_span"] == 0.02
    assert payload["options"]["min_abs_lobe_expected_value"] == 0.005

    entry = payload["entries"][0]
    assert entry["atom"] == 568
    assert entry["ring"] == 4
    assert entry["hard_lobe_crossing"]
    assert entry["positive_template_samples"] > 0
    assert entry["negative_template_samples"] > 0
    assert entry["positive_template_samples"] >= payload["options"]["min_samples_per_lobe"]
    assert entry["negative_template_samples"] >= payload["options"]["min_samples_per_lobe"]
    assert entry["template_sign_changes"] > 1
    assert entry["min_expected_relationship_value"] < 0.0
    assert entry["max_expected_relationship_value"] > 0.0
    assert entry["expected_relationship_span"] >= payload["options"]["min_expected_relationship_span"]
    assert abs(entry["min_expected_relationship_value"]) >= payload["options"]["min_abs_lobe_expected_value"]
    assert abs(entry["max_expected_relationship_value"]) >= payload["options"]["min_abs_lobe_expected_value"]
    assert entry["sample_count"] == 41
    assert len(entry["samples"]) == entry["sample_count"]

    samples = entry["samples"]
    values = [sample["expected_relationship_value"] for sample in samples]
    assert min(values) == entry["min_expected_relationship_value"]
    assert max(values) == entry["max_expected_relationship_value"]
    assert all("orca" in sample for sample in samples)
    assert all("biot_savart" in sample for sample in samples)

    fit = entry["fits"]["orca_total_T0"]
    assert fit["valid"]
    assert fit["sample_count"] == entry["sample_count"]
    assert fit["null_shift_count"] == entry["sample_count"] - 1
    assert len(fit["null_r2"]) == entry["sample_count"] - 1
    assert fit["scale"] > 0.0
    assert fit["r2"] > 0.5
    assert fit["r2"] > fit["null_median_r2"]
    assert fit["null_ge_real_fraction"] < 0.25

    distance_only = entry["confound_fits"]["distance_only_orca_total_T0"]
    angular_only = entry["confound_fits"]["angular_only_orca_total_T0"]
    assert distance_only["valid"]
    assert angular_only["valid"]
    assert fit["r2"] > distance_only["r2"]
    assert fit["r2"] > angular_only["r2"]

    predictor_fit = entry["predictor_diagnostics"]["expected_relationship_value_vs_biot_savart_T0"]
    assert predictor_fit["valid"]
    assert predictor_fit["sample_count"] == entry["sample_count"]

    para = entry["fits"]["orca_paramagnetic_T0"]
    assert para["valid"]
    assert para["sample_count"] == entry["sample_count"]

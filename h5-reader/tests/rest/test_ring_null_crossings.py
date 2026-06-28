"""REST contract for the operational ring-null collar collection."""

from __future__ import annotations

import math


def test_ring_null_collar_requires_filter_or_explicit_full_scan(rest):
    response = rest.client.post("/api/ring/null_crossings", json={})
    assert response.status_code == 400
    assert "full_scan" in response.json()["error"]


def test_ring_null_collar_collects_adjacent_dft_orca_pair(rest):
    response = rest.client.post(
        "/api/ring/null_crossings",
        json={"atom": 568, "ring": 4, "start_frame": 538, "end_frame": 540},
    )
    assert response.status_code == 200, response.text
    payload = response.json()

    assert payload["kind"] == "ring_null_collar_collection"
    assert payload["statistic"]["name"] == "ring_null_margin_A"
    assert payload["statistic"]["definition"] == "radial_A - sqrt(2) * abs(axial_A)"
    assert payload["summary"]["complete"]
    assert payload["summary"]["dft_pairs_scanned"] == 1
    assert payload["summary"]["entry_count"] == 1
    assert payload["summary"]["signal_stamp_count"] == 2
    assert payload["options"]["include_signal_stamps"]
    assert payload["options"]["stamp_radius_dft"] == 2
    assert "score" not in payload
    assert "samples" not in payload
    assert "dft_change_set" not in payload

    entry = payload["entries"][0]
    assert entry["atom"] == 568
    assert entry["ring"] == 4
    assert entry["atom_identity"]["atom_label_amber"]
    assert entry["atom_identity"]["residue_number"] > 0
    assert entry["atom_identity"]["residue_label_amber"]
    assert entry["ring_identity"]["type_name"]
    assert entry["ring_identity"]["kind"] == "aromatic"
    assert entry["ring_identity"]["parent_residue_number"] > 0
    assert entry["from"]["frame"] == 538
    assert entry["to"]["frame"] == 540
    assert entry["from"]["null"]["side"] == "axial"
    assert entry["to"]["null"]["side"] == "equatorial"
    assert payload["collar"]["mode"] == "zero_width_surface_crossing"
    assert payload["collar"]["physical_width_A"] == 0.0
    assert "not a finite physical aperture" in payload["collar"]["width_semantics"]

    event_frame = entry["event_frame"]
    assert event_frame["phase_coordinate"] == "signed_null_margin_A"
    assert 0.0 <= event_frame["zero_fraction"] <= 1.0
    assert math.isclose(
        event_frame["signed_null_margin_step_A"],
        event_frame["to_signed_null_margin_A"] - event_frame["from_signed_null_margin_A"],
        rel_tol=0.0,
        abs_tol=1e-9,
    )

    motion = entry["motion"]
    world = motion["world_vector_A"]
    expected_motion = math.sqrt(sum(component * component for component in world))
    assert math.isclose(motion["distance_A"], expected_motion, rel_tol=0.0, abs_tol=1e-9)

    stamps = entry["signal_stamps"]
    assert [stamp["frame"] for stamp in stamps] == [538, 540]
    assert stamps[0]["dft_ordinal_offset_from_zero"] < 0.0
    assert stamps[1]["dft_ordinal_offset_from_zero"] > 0.0
    for stamp in stamps:
        assert stamp["orca"]["present"]
        assert stamp["snapshot_present"]
        assert stamp["mopac"]["present"]
        assert set(stamp["mopac"]) >= {
            "present",
            "charge",
            "core",
            "coulomb_E",
            "coulomb_scalars",
            "coulomb_shielding",
            "coulomb_efg_backbone",
            "coulomb_efg_aromatic",
            "mc_shielding",
            "mc_category_T2",
            "mc_scalars",
        }
        assert stamp["mopac"]["charge"]["present"]
        assert stamp["mopac"]["core"]["present"]
        assert stamp["mopac"]["coulomb_E"]["present"]
        assert stamp["mopac"]["coulomb_scalars"]["present"]
        assert stamp["mopac"]["coulomb_shielding"]["present"]
        assert stamp["mopac"]["coulomb_efg_backbone"]["present"]
        assert stamp["mopac"]["coulomb_efg_aromatic"]["present"]

    for snapshot in (entry["from"], entry["to"]):
        null = snapshot["null"]
        expected = null["radial_A"] - math.sqrt(2.0) * null["abs_axial_A"]
        assert math.isclose(null["null_margin_A"], expected, rel_tol=0.0, abs_tol=1e-9)

        orca = snapshot["orca"]
        assert orca["element"] == "H"
        assert set(orca) >= {
            "total",
            "diamagnetic",
            "paramagnetic",
            "total_raw",
            "diamagnetic_raw",
            "paramagnetic_raw",
            "orca_coord_A",
        }
        assert len(orca["total_raw"]) == 3
        assert len(orca["total_raw"][0]) == 3
        assert set(orca["total"]) >= {"T0", "T1", "T2", "T2_magnitude"}
        assert snapshot["shape"]["total"]["valid"]

    delta_iso = entry["to"]["orca"]["total"]["T0"] - entry["from"]["orca"]["total"]["T0"]
    assert math.isclose(delta_iso, -2.421, abs_tol=1e-6)

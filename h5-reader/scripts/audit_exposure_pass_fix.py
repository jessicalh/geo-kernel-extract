#!/usr/bin/env python3
"""Proof-contract audits for the exposure-pass blocker fixes.

Uses Python's csv module intentionally; these files contain quoted fields and
must not be parsed with comma splitting.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import tempfile
from pathlib import Path
from typing import Any


def finite_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        parsed = float(value)
    except ValueError:
        return None
    return parsed if math.isfinite(parsed) else None


def int_value(value: str | None) -> int | None:
    if value is None or value == "":
        return None
    try:
        return int(value)
    except ValueError:
        return None


def audit_cross_axis(path: Path) -> dict[str, Any]:
    result: dict[str, Any] = {
        "path": str(path),
        "rows": 0,
        "zero_population_rows": 0,
        "zero_population_crosswalk_missing_rows": 0,
    }
    if not path.exists():
        result["skipped"] = "missing"
        return result

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"{path}: missing CSV header")
        population_column = (
            "cohort_cell_n_proteins"
            if "cohort_cell_n_proteins" in reader.fieldnames
            else "n_proteins"
        )
        if population_column not in reader.fieldnames or "comparability_reason" not in reader.fieldnames:
            raise ValueError(f"{path}: missing population or comparability columns")
        result["population_column"] = population_column
        for row in reader:
            result["rows"] += 1
            population = int_value(row.get(population_column))
            if population == 0:
                result["zero_population_rows"] += 1
                label = row.get("comparability_reason", "")
                if label.startswith("axis1_crosswalk_missing"):
                    result["zero_population_crosswalk_missing_rows"] += 1
    return result


def audit_bounded_sigma(path: Path, angle_threshold_deg: float) -> dict[str, Any]:
    result: dict[str, Any] = {
        "path": str(path),
        "rows": 0,
        "finite_pas_rows": 0,
        "pas_ordering_violations": 0,
        "pas_trace_violations": 0,
        "max_pas_trace_abs_error": 0.0,
        "finite_sp2_pyramidalization_rows": 0,
        "sp2_pyramidalization_angle_violations": 0,
        "max_abs_sp2_pyramidalization_angle_deg": 0.0,
        "angle_threshold_deg": angle_threshold_deg,
    }
    if not path.exists():
        result["skipped"] = "missing"
        return result

    required = {
        "sigma_iso",
        "pas_delta11",
        "pas_delta22",
        "pas_delta33",
        "hyb",
        "frame_kind",
        "pyramidalization_oop_angle_deg",
        "pyramidalization_support_class",
    }
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"{path}: missing CSV header")
        missing = required.difference(reader.fieldnames)
        if missing:
            raise ValueError(f"{path}: missing columns: {sorted(missing)}")
        for row in reader:
            result["rows"] += 1
            sigma_iso = finite_float(row.get("sigma_iso"))
            d11 = finite_float(row.get("pas_delta11"))
            d22 = finite_float(row.get("pas_delta22"))
            d33 = finite_float(row.get("pas_delta33"))
            if sigma_iso is not None and d11 is not None and d22 is not None and d33 is not None:
                result["finite_pas_rows"] += 1
                if not (d11 + 1e-10 >= d22 and d22 + 1e-10 >= d33):
                    result["pas_ordering_violations"] += 1
                trace_abs = abs((d11 + d22 + d33) - 3.0 * sigma_iso)
                scale = max(1.0, abs(d11) + abs(d22) + abs(d33) + abs(3.0 * sigma_iso))
                result["max_pas_trace_abs_error"] = max(
                    result["max_pas_trace_abs_error"],
                    trace_abs,
                )
                if trace_abs > max(1e-8, 1e-10 * scale):
                    result["pas_trace_violations"] += 1

            angle = finite_float(row.get("pyramidalization_oop_angle_deg"))
            support = row.get("pyramidalization_support_class", "")
            is_sp2_or_aromatic = (
                row.get("hyb", "").lower() == "sp2"
                or row.get("frame_kind", "") == "aromatic_ring_local"
            )
            if is_sp2_or_aromatic and support != "insufficient" and angle is not None:
                result["finite_sp2_pyramidalization_rows"] += 1
                abs_angle = abs(angle)
                result["max_abs_sp2_pyramidalization_angle_deg"] = max(
                    result["max_abs_sp2_pyramidalization_angle_deg"],
                    abs_angle,
                )
                if abs_angle > angle_threshold_deg:
                    result["sp2_pyramidalization_angle_violations"] += 1
    return result


def audit_paths(paths: list[Path], angle_threshold_deg: float) -> dict[str, Any]:
    cross_axis = []
    bounded = []
    for base in paths:
        cross_axis.append(audit_cross_axis(base / "cross_axis_overlay.csv"))
        bounded.append(audit_bounded_sigma(base / "bounded_sigma.csv", angle_threshold_deg))

    summary = {
        "cross_axis": cross_axis,
        "bounded_sigma": bounded,
        "failures": {
            "zero_population_crosswalk_missing_rows": sum(
                item.get("zero_population_crosswalk_missing_rows", 0)
                for item in cross_axis
                if item.get("skipped") is None
            ),
            "pas_ordering_violations": sum(
                item.get("pas_ordering_violations", 0)
                for item in bounded
                if item.get("skipped") is None
            ),
            "pas_trace_violations": sum(
                item.get("pas_trace_violations", 0)
                for item in bounded
                if item.get("skipped") is None
            ),
            "sp2_pyramidalization_angle_violations": sum(
                item.get("sp2_pyramidalization_angle_violations", 0)
                for item in bounded
                if item.get("skipped") is None
            ),
        },
    }
    return summary


def has_failures(summary: dict[str, Any]) -> bool:
    failures = summary["failures"]
    return any(value != 0 for value in failures.values())


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def run_self_test(angle_threshold_deg: float) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        bad = Path(tmp) / "bad"
        write_csv(
            bad / "cross_axis_overlay.csv",
            [
                {
                    "cohort_cell_n_proteins": "0",
                    "n_proteins": "0",
                    "comparability_reason": "axis1_crosswalk_missing_populated_unjoined",
                }
            ],
        )
        write_csv(
            bad / "bounded_sigma.csv",
            [
                {
                    "sigma_iso": "2.0",
                    "pas_delta11": "1.0",
                    "pas_delta22": "2.0",
                    "pas_delta33": "3.0",
                    "hyb": "sp2",
                    "frame_kind": "aromatic_ring_local",
                    "pyramidalization_oop_angle_deg": str(angle_threshold_deg + 1.0),
                    "pyramidalization_support_class": "insufficient",
                },
                {
                    "sigma_iso": "2.0",
                    "pas_delta11": "3.0",
                    "pas_delta22": "2.0",
                    "pas_delta33": "2.0",
                    "hyb": "sp2",
                    "frame_kind": "aromatic_ring_local",
                    "pyramidalization_oop_angle_deg": str(angle_threshold_deg + 1.0),
                    "pyramidalization_support_class": "full",
                },
            ],
        )
        summary = audit_paths([bad], angle_threshold_deg)
        failures = summary["failures"]
        expected = {
            "zero_population_crosswalk_missing_rows",
            "pas_ordering_violations",
            "pas_trace_violations",
            "sp2_pyramidalization_angle_violations",
        }
        missing = [name for name in expected if failures.get(name, 0) == 0]
        if missing:
            raise SystemExit(f"self-test failed to trigger: {missing}\n{json.dumps(summary, indent=2)}")
    print("self_test_ok")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*", type=Path, help="smoke output directories")
    parser.add_argument("--angle-threshold-deg", type=float, default=12.0)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()

    if args.self_test:
        run_self_test(args.angle_threshold_deg)
        return 0
    if not args.paths:
        parser.error("at least one smoke output directory is required unless --self-test is used")

    summary = audit_paths(args.paths, args.angle_threshold_deg)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 1 if has_failures(summary) else 0


if __name__ == "__main__":
    raise SystemExit(main())

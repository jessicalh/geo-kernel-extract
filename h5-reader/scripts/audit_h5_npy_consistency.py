#!/usr/bin/env python3
"""Compare July trajectory HDF5 products with authoritative frame NPYs."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Callable

import h5py
import numpy as np


RTOL = 1.0e-13
ATOL = 1.0e-13


@dataclass(frozen=True)
class FieldMapping:
    name: str
    h5_dataset: str
    npy_name: str
    select_npy: Callable[[np.ndarray], np.ndarray] = lambda values: values


DIRECT_MAPPINGS = (
    FieldMapping("positions", "positions/xyz", "pos.npy"),
    FieldMapping("Biot-Savart tensor", "bs_shielding_time_series/xyz", "bs_shielding.npy"),
    FieldMapping("Haigh-Mallion tensor", "hm_shielding_time_series/xyz", "hm_shielding.npy"),
    FieldMapping(
        "MOPAC Coulomb EFG T2",
        "mopac_coulomb_efg_time_series/t2",
        "mopac_coulomb_efg.npy",
        lambda values: values[:, 4:9],
    ),
    FieldMapping(
        "MOPAC McConnell tensor",
        "mopac_mc_shielding_time_series/xyz",
        "mopac_mc_shielding.npy",
    ),
    FieldMapping("SASA", "sasa_time_series/sasa", "atom_sasa.npy"),
    FieldMapping("AIMNet charges", "aimnet2_charge_time_series/charge", "aimnet2_charges.npy"),
    FieldMapping(
        "AIMNet embedding",
        "aimnet2_embedding_time_series/embedding",
        "aimnet2_aim.npy",
    ),
    FieldMapping(
        "AIMNet charge-response vector",
        "aimnet2_charge_response_gradient_time_series/charge_response_gradient_vector",
        "aimnet2_charge_response_gradient.npy",
    ),
    FieldMapping(
        "AIMNet charge-response scalar",
        "aimnet2_charge_response_gradient_time_series/charge_response_gradient_scalar",
        "aimnet2_charge_response_gradient_scalar.npy",
    ),
    FieldMapping("APBS electric field", "apbs_efield_time_series/xyz", "apbs_E.npy"),
    FieldMapping("APBS EFG T2", "apbs_efg_time_series/t2", "apbs_efg.npy"),
    FieldMapping(
        "bonded energy total",
        "bonded_energy_time_series/total",
        "bonded_energy.npy",
        lambda values: values[:, 7],
    ),
    FieldMapping("water electric field", "water_field_time_series/efield", "water_efield.npy"),
    FieldMapping(
        "first-shell water electric field",
        "water_field_time_series/efield_first",
        "water_efield_first.npy",
    ),
    FieldMapping("water EFG T2", "water_field_time_series/efg", "water_efg.npy"),
    FieldMapping(
        "first-shell water EFG T2",
        "water_field_time_series/efg_first",
        "water_efg_first.npy",
    ),
    FieldMapping(
        "first-shell water count",
        "water_field_time_series/n_first",
        "water_shell_counts.npy",
        lambda values: values[:, 0],
    ),
    FieldMapping(
        "second-shell water count",
        "water_field_time_series/n_second",
        "water_shell_counts.npy",
        lambda values: values[:, 1],
    ),
    FieldMapping(
        "nearest-ion distance",
        "hydration_shell_time_series/nearest_ion_distance",
        "hydration_shell.npy",
        lambda values: values[:, 2],
    ),
    FieldMapping(
        "nearest-ion charge",
        "hydration_shell_time_series/nearest_ion_charge",
        "hydration_shell.npy",
        lambda values: values[:, 3],
    ),
    FieldMapping(
        "mean water-dipole cosine",
        "hydration_shell_time_series/mean_water_dipole_cos",
        "hydration_shell.npy",
        lambda values: values[:, 1],
    ),
    FieldMapping(
        "hydration-shell asymmetry",
        "hydration_shell_time_series/half_shell_asymmetry",
        "hydration_shell.npy",
        lambda values: values[:, 0],
    ),
    FieldMapping(
        "water dipole vector",
        "hydration_geometry_time_series/dipole_vector",
        "water_polarization.npy",
        lambda values: values[:, 0:3],
    ),
    FieldMapping(
        "water surface normal",
        "hydration_geometry_time_series/surface_normal",
        "water_polarization.npy",
        lambda values: values[:, 3:6],
    ),
    FieldMapping(
        "water half-shell asymmetry",
        "hydration_geometry_time_series/half_shell_asymmetry",
        "water_polarization.npy",
        lambda values: values[:, 6],
    ),
    FieldMapping(
        "water dipole alignment",
        "hydration_geometry_time_series/dipole_alignment",
        "water_polarization.npy",
        lambda values: values[:, 7],
    ),
    FieldMapping(
        "water dipole coherence",
        "hydration_geometry_time_series/dipole_coherence",
        "water_polarization.npy",
        lambda values: values[:, 8],
    ),
    FieldMapping(
        "water first-shell count",
        "hydration_geometry_time_series/first_shell_count",
        "water_polarization.npy",
        lambda values: values[:, 9],
    ),
)


MCCONNELL_FIXED_FIELDS = (
    "mc_peptide_co_fixed.npy",
    "mc_peptide_cn_fixed.npy",
    "mc_backbone_other_fixed.npy",
    "mc_sidechain_co_fixed.npy",
    "mc_sidechain_other_fixed.npy",
    "mc_disulfide_fixed.npy",
    "mc_aromatic_fixed.npy",
    "mc_backbone_xh_fixed.npy",
    "mc_sidechain_xh_fixed.npy",
    "mc_s_h_fixed.npy",
)


RETIRED_GROUPS = (
    "mopac_vs_ff14sb_reconciliation",
    "tripeptide_bb_shielding_time_series",
    "tripeptide_neighbor_shielding_time_series",
    "tripeptide_bb_residual_vec_time_series",
    "tripeptide_neighbor_residual_vec_prev_time_series",
    "tripeptide_neighbor_residual_vec_next_time_series",
    "tripeptide_bb_method_tag_time_series",
)


EXPECTED_KERNEL_CHANNELS = (
    "bs_T0",
    "bs_absT2",
    "hm_T0",
    "hm_absT2",
    "mc_T0",
    "mc_absT2",
    "apbs_absT2",
)


def extraction_root(path: Path) -> Path:
    if (path / "trajectory.h5").is_file():
        return path
    nested = path / "extraction"
    if (nested / "trajectory.h5").is_file():
        return nested
    raise FileNotFoundError(f"no trajectory.h5 under {path}")


def decode_strings(values: np.ndarray) -> tuple[str, ...]:
    return tuple(
        value.decode("utf-8") if isinstance(value, bytes) else str(value)
        for value in values
    )


def selected_frame_positions(frame_count: int, all_frames: bool) -> list[int]:
    if all_frames:
        return list(range(frame_count))
    return sorted({0, frame_count // 2, frame_count - 1})


def group_for_dataset(trajectory: h5py.Group, dataset_path: str) -> h5py.Group:
    group_name = dataset_path.split("/", maxsplit=1)[0]
    return trajectory[group_name]


def compare_values(actual: np.ndarray, expected: np.ndarray) -> str:
    if actual.shape != expected.shape:
        return f"shape mismatch: H5 {actual.shape}, NPY {expected.shape}"
    if np.array_equal(actual, expected, equal_nan=True):
        return "exact"
    if np.allclose(actual, expected, rtol=RTOL, atol=ATOL, equal_nan=True):
        return "close"
    finite = np.isfinite(actual) & np.isfinite(expected)
    max_error = float(np.max(np.abs(actual[finite] - expected[finite]))) if finite.any() else float("nan")
    return f"value mismatch: max_abs_error={max_error:.17g}"


def validate_group_frame_axis(
    trajectory: h5py.Group,
    group_name: str,
    frame_indices: np.ndarray,
    frame_times: np.ndarray,
) -> list[str]:
    errors: list[str] = []
    group = trajectory[group_name]
    if "frame_indices" in group and not np.array_equal(group["frame_indices"][:], frame_indices):
        errors.append(f"{group_name}: frame_indices differ from /trajectory/frames")
    if "frame_times" in group and not np.array_equal(group["frame_times"][:], frame_times):
        errors.append(f"{group_name}: frame_times differ from /trajectory/frames")
    return errors


def audit_run(path: Path, all_frames: bool) -> bool:
    root = extraction_root(path)
    h5_path = root / "trajectory.h5"
    npy_root = root / "npys"
    errors: list[str] = []
    exact = 0
    close = 0
    masked = 0

    with h5py.File(h5_path, "r") as handle:
        trajectory = handle["trajectory"]
        frame_indices = trajectory["frames/original_index"][:]
        frame_times = trajectory["frames/time_ps"][:]
        frame_positions = selected_frame_positions(len(frame_indices), all_frames)

        for retired in RETIRED_GROUPS:
            if retired in trajectory:
                errors.append(f"retired H5 group is present: /trajectory/{retired}")

        for kernel_group_name in ("kernel_dynamics", "kernel_coherence"):
            if kernel_group_name not in trajectory:
                errors.append(f"missing H5 group: /trajectory/{kernel_group_name}")
                continue
            names = decode_strings(trajectory[kernel_group_name]["channel_names"][:])
            if names != EXPECTED_KERNEL_CHANNELS:
                errors.append(
                    f"{kernel_group_name}: channels {names!r}, expected {EXPECTED_KERNEL_CHANNELS!r}"
                )

        canonical = "mopac_coulomb_efg_time_series"
        legacy = "mopac_coulomb_shielding_time_series"
        if canonical not in trajectory:
            errors.append(f"missing canonical H5 group: /trajectory/{canonical}")
        if legacy in trajectory:
            alias = trajectory[legacy].attrs.get("alias_of", "")
            if alias != f"/trajectory/{canonical}":
                errors.append(f"{legacy}: alias_of={alias!r}")

        checked_groups: set[str] = set()
        for mapping in DIRECT_MAPPINGS:
            group_name = mapping.h5_dataset.split("/", maxsplit=1)[0]
            if group_name not in trajectory:
                errors.append(f"{mapping.name}: missing H5 group {group_name}")
                continue
            if mapping.h5_dataset not in trajectory:
                errors.append(f"{mapping.name}: missing H5 dataset {mapping.h5_dataset}")
                continue
            if group_name not in checked_groups:
                errors.extend(
                    validate_group_frame_axis(
                        trajectory, group_name, frame_indices, frame_times
                    )
                )
                checked_groups.add(group_name)

            dataset = trajectory[mapping.h5_dataset]
            group = group_for_dataset(trajectory, mapping.h5_dataset)
            source_mask = group.get("source_attached_per_frame")
            for frame_position in frame_positions:
                original_index = int(frame_indices[frame_position])
                frame_dir = npy_root / f"frame_{original_index:06d}"
                npy_path = frame_dir / mapping.npy_name
                if source_mask is not None and int(source_mask[frame_position]) == 0:
                    masked += 1
                    continue
                if not npy_path.is_file():
                    errors.append(f"{mapping.name}: missing {npy_path}")
                    continue
                actual = dataset[:, frame_position, ...]
                expected = mapping.select_npy(np.load(npy_path, allow_pickle=False))
                result = compare_values(actual, expected)
                if result == "exact":
                    exact += 1
                elif result == "close":
                    close += 1
                else:
                    errors.append(
                        f"{mapping.name}, frame {original_index}: {result}"
                    )

        mc_dataset = trajectory["mc_shielding_time_series/xyz"]
        for frame_position in frame_positions:
            original_index = int(frame_indices[frame_position])
            frame_dir = npy_root / f"frame_{original_index:06d}"
            components = [
                np.load(frame_dir / name, allow_pickle=False)
                for name in MCCONNELL_FIXED_FIELDS
            ]
            expected = np.sum(np.stack(components, axis=0), axis=0)
            result = compare_values(mc_dataset[:, frame_position, :], expected)
            if result == "exact":
                exact += 1
            elif result == "close":
                close += 1
            else:
                errors.append(
                    f"McConnell fixed aggregate, frame {original_index}: {result}"
                )

        if canonical in trajectory and legacy in trajectory:
            for frame_position in frame_positions:
                result = compare_values(
                    trajectory[f"{canonical}/t2"][:, frame_position, :],
                    trajectory[f"{legacy}/t2"][:, frame_position, :],
                )
                if result != "exact":
                    errors.append(
                        f"deprecated MOPAC alias, frame "
                        f"{int(frame_indices[frame_position])}: {result}"
                    )

    mode = "all" if all_frames else "first/middle/last"
    print(f"{root}:")
    print(
        f"  frames checked: {len(frame_positions)}/{len(frame_indices)} ({mode}); "
        f"field-frame comparisons exact={exact}, close={close}, masked={masked}"
    )
    if errors:
        for error in errors:
            print(f"  ERROR: {error}")
        return False
    print("  PASS")
    return True


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "runs",
        nargs="+",
        type=Path,
        help="extraction directories, or one-off directories containing extraction/",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="compare every frame instead of first, middle, and last",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    passed = True
    for run in args.runs:
        try:
            passed = audit_run(run, args.all) and passed
        except Exception as error:
            passed = False
            print(f"{run}: ERROR: {error}", file=sys.stderr)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())

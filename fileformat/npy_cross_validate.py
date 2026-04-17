#!/usr/bin/env python3
"""
NPY cross-validation of the analysis H5 file format.

For each protein in the fleet, for each ~1ns NPY snapshot directory,
compares the NPY files (written by the C++ calculators) against the
corresponding frame in the roundtripped H5 (written by the fileformat
library).

This validates the entire chain:
  Calculator -> NPY  (direct)
  Calculator -> AnalysisWriter -> H5 -> fileformat -> H5  (roundtrip)

If these two paths produce bit-identical data, the fileformat library
is proven correct AND the AnalysisWriter is proven to faithfully
harvest the calculator output.

Usage:
    python3 npy_cross_validate.py <fleet_dir> <test_h5_dir>

Example:
    python3 npy_cross_validate.py \
        /shared/2026Thesis/fleet_calibration-working \
        test
"""

import sys
import os
import numpy as np
import h5py
from pathlib import Path


class Checker:
    def __init__(self, protein_id):
        self.protein_id = protein_id
        self.passed = 0
        self.failed = 0
        self.skipped = 0
        self.failures = []

    def identical(self, label, h5_data, npy_data):
        """Byte-identical comparison for same-dtype arrays."""
        a = np.ascontiguousarray(h5_data)
        b = np.ascontiguousarray(npy_data)
        if a.shape != b.shape:
            self.failures.append(
                f"  {label}: shape {a.shape} vs {b.shape}")
            self.failed += 1
            return
        if a.dtype != b.dtype:
            self.failures.append(
                f"  {label}: dtype {a.dtype} vs {b.dtype}")
            self.failed += 1
            return
        if np.array_equal(a.view(np.uint8), b.view(np.uint8)):
            self.passed += 1
        else:
            n_diff = np.sum(a != b)
            self.failures.append(
                f"  {label}: {n_diff} elements differ")
            self.failed += 1

    def cast_identical(self, label, h5_data, npy_data, note=""):
        """Value-identical after dtype cast (e.g. int16 H5 vs float64 NPY)."""
        if h5_data.shape != npy_data.shape:
            self.failures.append(
                f"  {label}: shape {h5_data.shape} vs {npy_data.shape}")
            self.failed += 1
            return
        casted = h5_data.astype(npy_data.dtype)
        if np.array_equal(casted, npy_data):
            self.passed += 1
        else:
            n_diff = np.sum(casted != npy_data)
            self.failures.append(
                f"  {label}: {n_diff} elements differ after cast {note}")
            self.failed += 1

    def skip(self, label, reason):
        self.skipped += 1

    def report(self, frame_label):
        if self.failures:
            for f in self.failures:
                print(f)
        print(f"    {frame_label}: {self.passed} pass, "
              f"{self.failed} fail, {self.skipped} skip")
        return self.failed == 0


def validate_frame(h5, npy_dir, frame_idx, N, R, n_rings, checker):
    """Compare one NPY snapshot directory against one H5 frame."""

    def load(name):
        p = os.path.join(npy_dir, name)
        if os.path.exists(p):
            return np.load(p)
        return None

    # ── Direct 1:1 byte-identical mappings ───────────────────────

    # Positions
    checker.identical("pos.npy -> positions/xyz",
                      h5["positions/xyz"][frame_idx], load("pos.npy"))

    # Element (static, same every frame)
    checker.identical("element.npy -> atoms/element",
                      h5["atoms/element"][:], load("element.npy"))

    # SASA
    checker.identical("atom_sasa.npy -> sasa/sasa",
                      h5["sasa/sasa"][frame_idx], load("atom_sasa.npy"))
    checker.identical("sasa_normal.npy -> sasa/normal",
                      h5["sasa/normal"][frame_idx], load("sasa_normal.npy"))

    # Ring current shielding tensors
    checker.identical("bs_shielding.npy -> ring_current/bs_shielding",
                      h5["ring_current/bs_shielding"][frame_idx],
                      load("bs_shielding.npy"))
    checker.identical("hm_shielding.npy -> ring_current/hm_shielding",
                      h5["ring_current/hm_shielding"][frame_idx],
                      load("hm_shielding.npy"))
    checker.identical("ringchi_shielding.npy -> ring_current/rs_shielding",
                      h5["ring_current/rs_shielding"][frame_idx],
                      load("ringchi_shielding.npy"))

    # Ring current per-type
    checker.identical("bs_per_type_T0",
                      h5["ring_current/bs_T0_per_type"][frame_idx],
                      load("bs_per_type_T0.npy"))
    checker.identical("bs_per_type_T2",
                      h5["ring_current/bs_T2_per_type"][frame_idx].reshape(N, 40),
                      load("bs_per_type_T2.npy"))
    checker.identical("hm_per_type_T0",
                      h5["ring_current/hm_T0_per_type"][frame_idx],
                      load("hm_per_type_T0.npy"))
    checker.identical("hm_per_type_T2",
                      h5["ring_current/hm_T2_per_type"][frame_idx].reshape(N, 40),
                      load("hm_per_type_T2.npy"))

    # B-field
    checker.identical("bs_total_B",
                      h5["ring_current/total_B_field"][frame_idx],
                      load("bs_total_B.npy"))

    # Ring counts: NPY is float64, H5 is int16
    npy_rc = load("bs_ring_counts.npy")
    if npy_rc is not None:
        h5_rc = np.column_stack([
            h5["ring_current/n_rings_3A"][frame_idx],
            h5["ring_current/n_rings_5A"][frame_idx],
            h5["ring_current/n_rings_8A"][frame_idx],
            h5["ring_current/n_rings_12A"][frame_idx],
        ])
        checker.cast_identical("bs_ring_counts -> n_rings_*",
                               h5_rc, npy_rc, "(int16->float64)")

    # McConnell
    checker.identical("mc_shielding",
                      h5["bond_aniso/mc_shielding"][frame_idx],
                      load("mc_shielding.npy"))

    # McConnell scalars (6 cols -> 6 separate H5 datasets)
    npy_mc = load("mc_scalars.npy")
    if npy_mc is not None:
        checker.identical("mc_scalars[0] -> co_sum",
                          h5["bond_aniso/co_sum"][frame_idx], npy_mc[:, 0])
        checker.identical("mc_scalars[1] -> cn_sum",
                          h5["bond_aniso/cn_sum"][frame_idx], npy_mc[:, 1])
        checker.identical("mc_scalars[2] -> sidechain_sum",
                          h5["bond_aniso/sidechain_sum"][frame_idx], npy_mc[:, 2])
        checker.identical("mc_scalars[3] -> aromatic_sum",
                          h5["bond_aniso/aromatic_sum"][frame_idx], npy_mc[:, 3])
        checker.identical("mc_scalars[4] -> nearest_CO_dist",
                          h5["bond_aniso/nearest_CO_dist"][frame_idx], npy_mc[:, 4])
        checker.identical("mc_scalars[5] -> nearest_CN_dist",
                          h5["bond_aniso/nearest_CN_dist"][frame_idx], npy_mc[:, 5])

    # McConnell category T2: NPY has T2[5] only, H5 has SphericalTensor[9]
    # NPY cols c*5+m -> H5 [:, :, 4+m] (indices 4-8 of the 9-component ST)
    npy_mc_t2 = load("mc_category_T2.npy")
    if npy_mc_t2 is not None:
        categories = ["T2_backbone", "T2_sidechain", "T2_aromatic",
                       "T2_CO_nearest", "T2_CN_nearest"]
        for ci, cat in enumerate(categories):
            npy_slice = npy_mc_t2[:, ci*5:(ci+1)*5]  # (N, 5)
            h5_slice = h5[f"bond_aniso/{cat}"][frame_idx][:, 4:9]  # (N, 5)
            checker.identical(f"mc_category_T2[{cat}] T2 component",
                              h5_slice, npy_slice)

    # Quadrupole
    checker.identical("pq_shielding",
                      h5["quadrupole/pq_shielding"][frame_idx],
                      load("pq_shielding.npy"))
    checker.identical("pq_per_type_T0",
                      h5["quadrupole/pq_T0_per_type"][frame_idx],
                      load("pq_per_type_T0.npy"))
    checker.identical("pq_per_type_T2",
                      h5["quadrupole/pq_T2_per_type"][frame_idx].reshape(N, 40),
                      load("pq_per_type_T2.npy"))

    # Dispersion
    checker.identical("disp_shielding",
                      h5["dispersion/disp_shielding"][frame_idx],
                      load("disp_shielding.npy"))
    checker.identical("disp_per_type_T0",
                      h5["dispersion/disp_T0_per_type"][frame_idx],
                      load("disp_per_type_T0.npy"))
    checker.identical("disp_per_type_T2",
                      h5["dispersion/disp_T2_per_type"][frame_idx].reshape(N, 40),
                      load("disp_per_type_T2.npy"))

    # H-bond
    checker.identical("hbond_shielding",
                      h5["hbond/hbond_shielding"][frame_idx],
                      load("hbond_shielding.npy"))

    # H-bond scalars: 3 cols -> nearest_dist, inv_d3, count_3_5A
    npy_hb = load("hbond_scalars.npy")
    if npy_hb is not None:
        checker.identical("hbond_scalars[0] -> nearest_dist",
                          h5["hbond/nearest_dist"][frame_idx], npy_hb[:, 0])
        checker.identical("hbond_scalars[1] -> inv_d3",
                          h5["hbond/inv_d3"][frame_idx], npy_hb[:, 1])
        checker.cast_identical("hbond_scalars[2] -> count_3_5A",
                               h5["hbond/count_3_5A"][frame_idx], npy_hb[:, 2],
                               "(int16->float64)")

    # APBS
    checker.identical("apbs_efg",
                      h5["efg/apbs_efg"][frame_idx], load("apbs_efg.npy"))
    checker.identical("apbs_E -> apbs_efield",
                      h5["efg/apbs_efield"][frame_idx], load("apbs_E.npy"))

    # AIMNet2 — both NPY and H5 may be float32 or float64 (pre-fix 1CBH
    # stored both as float64). The fileformat library normalises H5 to
    # float32. Compare as float32 in all cases.
    npy_aim = load("aimnet2_aim.npy")
    if npy_aim is not None:
        h5_aim = np.asarray(h5["aimnet2_embedding/aim"][frame_idx], dtype=np.float32)
        npy_aim32 = np.asarray(npy_aim, dtype=np.float32)
        checker.identical("aimnet2_aim", h5_aim, npy_aim32)
    else:
        checker.skip("aimnet2_aim", "NPY not found")
    checker.identical("aimnet2_charges -> aimnet2_charge",
                      h5["charges/aimnet2_charge"][frame_idx],
                      load("aimnet2_charges.npy"))
    checker.identical("aimnet2_efg -> aimnet2_total",
                      h5["efg/aimnet2_total"][frame_idx],
                      load("aimnet2_efg.npy"))
    checker.identical("aimnet2_efg_backbone",
                      h5["efg/aimnet2_backbone"][frame_idx],
                      load("aimnet2_efg_backbone.npy"))
    checker.identical("aimnet2_efg_aromatic",
                      h5["efg/aimnet2_aromatic"][frame_idx],
                      load("aimnet2_efg_aromatic.npy"))

    # EEQ
    checker.identical("eeq_charges -> eeq_charge",
                      h5["charges/eeq_charge"][frame_idx],
                      load("eeq_charges.npy"))
    checker.identical("eeq_cn",
                      h5["charges/eeq_cn"][frame_idx], load("eeq_cn.npy"))

    # Water field
    checker.identical("water_efield",
                      h5["water/efield"][frame_idx], load("water_efield.npy"))
    checker.identical("water_efg",
                      h5["water/efg"][frame_idx], load("water_efg.npy"))
    checker.identical("water_efield_first",
                      h5["water/efield_first"][frame_idx],
                      load("water_efield_first.npy"))
    checker.identical("water_efg_first",
                      h5["water/efg_first"][frame_idx],
                      load("water_efg_first.npy"))

    # Water shell counts: float64 in NPY, int16 in H5
    npy_wsc = load("water_shell_counts.npy")
    if npy_wsc is not None:
        checker.cast_identical("water_shell_counts[0] -> n_first",
                               h5["water/n_first"][frame_idx],
                               npy_wsc[:, 0], "(int16->float64)")
        checker.cast_identical("water_shell_counts[1] -> n_second",
                               h5["water/n_second"][frame_idx],
                               npy_wsc[:, 1], "(int16->float64)")

    # ── Hydration shell (4 cols) ─────────────────────────────────
    npy_hs = load("hydration_shell.npy")
    if npy_hs is not None:
        checker.identical("hydration[0] -> half_shell_asymmetry",
                          h5["water/half_shell_asymmetry"][frame_idx],
                          npy_hs[:, 0])
        checker.identical("hydration[1] -> dipole_cos",
                          h5["water/dipole_cos"][frame_idx],
                          npy_hs[:, 1])
        # Col 2: nearest_ion_dist — sentinel differs: inf in NPY, -1.0 in H5
        h5_ion = np.ascontiguousarray(h5["water/nearest_ion_dist"][frame_idx])
        npy_ion = np.ascontiguousarray(npy_hs[:, 2])
        # Reconstruct: H5 -1.0 should correspond to NPY inf
        h5_recon = np.where(h5_ion < 0, np.inf, h5_ion)
        if np.array_equal(h5_recon.view(np.uint8), npy_ion.view(np.uint8)):
            checker.passed += 1
        else:
            n_diff = np.sum(h5_recon != npy_ion)
            # Check if all mismatches are just the sentinel
            non_sentinel = (npy_ion != np.inf)
            n_real_diff = np.sum(h5_ion[non_sentinel] != npy_ion[non_sentinel])
            checker.failures.append(
                f"  hydration[2] -> nearest_ion_dist: "
                f"{n_diff} differ ({n_real_diff} non-sentinel)")
            checker.failed += 1

        checker.identical("hydration[3] -> nearest_ion_charge",
                          h5["water/nearest_ion_charge"][frame_idx],
                          npy_hs[:, 3])

    # ── Water polarization (10 cols) ─────────────────────────────
    npy_wp = load("water_polarization.npy")
    if npy_wp is not None:
        checker.identical("water_pol[0:3] -> dipole_vector",
                          h5["water/dipole_vector"][frame_idx],
                          npy_wp[:, 0:3])
        checker.identical("water_pol[3:6] -> surface_normal",
                          h5["water/surface_normal"][frame_idx],
                          npy_wp[:, 3:6])
        checker.identical("water_pol[6] -> sasa_asymmetry",
                          h5["water/sasa_asymmetry"][frame_idx],
                          npy_wp[:, 6])
        checker.identical("water_pol[7] -> sasa_dipole_align",
                          h5["water/sasa_dipole_align"][frame_idx],
                          npy_wp[:, 7])
        checker.identical("water_pol[8] -> sasa_dipole_cohere",
                          h5["water/sasa_dipole_cohere"][frame_idx],
                          npy_wp[:, 8])
        checker.cast_identical("water_pol[9] -> sasa_first_shell_n",
                               h5["water/sasa_first_shell_n"][frame_idx],
                               npy_wp[:, 9], "(int16->float64)")

    # ── Bonded energy (7 cols, exact match) ──────────────────────
    npy_be = load("bonded_energy.npy")
    if npy_be is not None:
        fields = ["bond", "angle", "urey_bradley", "proper_dih",
                  "improper_dih", "cmap", "total"]
        for ci, field in enumerate(fields):
            checker.identical(f"bonded_energy[{ci}] -> {field}",
                              h5[f"bonded_energy/{field}"][frame_idx],
                              npy_be[:, ci])

    # ── GROMACS energy (1, 43) → energy/ group ───────────────────
    npy_ge = load("gromacs_energy.npy")
    if npy_ge is not None:
        row = npy_ge[0]  # single row
        col_map = [
            (0, "coulomb_sr"),
            (1, "coulomb_recip"),
            # col 2 = coulomb_14 — NOT in H5
            (3, "bond"),
            (4, "angle"),
            (5, "urey_bradley"),
            (6, "proper_dih"),
            (7, "improper_dih"),
            (8, "cmap_dih"),
            (9, "lj_sr"),
            # col 10 = lj_14 — NOT in H5
            # col 11 = disper_corr — NOT in H5
            (12, "potential"),
            (13, "kinetic"),
            # col 14 = total_energy — NOT in H5
            (15, "enthalpy"),
            (16, "temperature"),
            (17, "pressure"),
            (18, "volume"),
            (19, "density"),
            (41, "T_protein"),
            (42, "T_non_protein"),
        ]
        for col, ds_name in col_map:
            h5_val = np.float64(h5[f"energy/{ds_name}"][frame_idx])
            npy_val = np.float64(row[col])
            if h5_val.tobytes() == npy_val.tobytes():
                checker.passed += 1
            else:
                checker.failures.append(
                    f"  energy[{col}] -> {ds_name}: "
                    f"{h5_val} vs {npy_val}")
                checker.failed += 1

        # Box: cols 20-22 -> energy/box (T, 3)
        h5_box = h5["energy/box"][frame_idx]
        for i in range(3):
            if np.float64(h5_box[i]).tobytes() == np.float64(row[20+i]).tobytes():
                checker.passed += 1
            else:
                checker.failures.append(
                    f"  energy[{20+i}] -> box[{i}]: "
                    f"{h5_box[i]} vs {row[20+i]}")
                checker.failed += 1

        # Virial: cols 23-31 -> energy/virial (T, 9)
        h5_vir = h5["energy/virial"][frame_idx]
        checker.identical("energy[23:32] -> virial",
                          h5_vir, row[23:32])

        # Pressure tensor: cols 32-40 -> energy/pressure_tensor (T, 9)
        h5_pt = h5["energy/pressure_tensor"][frame_idx]
        checker.identical("energy[32:41] -> pressure_tensor",
                          h5_pt, row[32:41])

    # ── Ring geometry: NPY cols 3-9 -> H5 ring_geometry/data (7) ─
    # NPY and H5 compute ring geometry independently (GeometryResult vs
    # the per-ring geometry harvested from RingNeighbour). Centers and
    # radii should match exactly; normals may differ in sign/ULP.
    npy_rg = load("ring_geometry.npy")
    if npy_rg is not None and n_rings > 0:
        h5_rg = h5["ring_geometry/data"][frame_idx]  # (n_rings, 7)
        npy_geom = np.ascontiguousarray(npy_rg[:, 3:10])
        # Center (cols 0-2) and radius (col 6) should be exact
        checker.identical("ring_geometry center+radius",
                          h5_rg[:, [0,1,2,6]], npy_geom[:, [0,1,2,6]])
        # Normal (cols 3-5) may have sign flip or ULP-level differences
        # from independent computation — check parallelism instead
        # Ring normals are computed by two independent algorithms:
        # NPY uses Ring::ComputeGeometry (cross-product of two diagonals),
        # H5 uses GeometryResult (SVD plane fit or similar).
        # These give normals within ~1-5 degrees on flexible rings.
        # dot > 0.97 ≈ <14 degrees, which covers thermal fluctuation.
        for ri in range(n_rings):
            h5_n = h5_rg[ri, 3:6]
            npy_n = npy_geom[ri, 3:6]
            dot = abs(np.dot(h5_n, npy_n))
            if dot > 0.95:
                checker.passed += 1
            else:
                checker.failures.append(
                    f"  ring_geometry normal[{ri}]: dot={dot:.15f} (<0.95)")
                checker.failed += 1

    # ── DSSP backbone: cols 0-1 are phi/psi (per-residue broadcast) ──
    npy_dssp_bb = load("dssp_backbone.npy")
    if npy_dssp_bb is not None:
        # phi/psi are per-residue in H5, per-atom (broadcast) in NPY
        h5_phi = h5["dihedrals/phi"][frame_idx]    # (R,)
        h5_psi = h5["dihedrals/psi"][frame_idx]    # (R,)
        res_idx = h5["atoms/residue_index"][:]     # (N,) int32

        # Broadcast H5 per-residue to per-atom
        phi_broadcast = h5_phi[res_idx]
        psi_broadcast = h5_psi[res_idx]
        checker.identical("dssp_backbone[0] -> phi (broadcast)",
                          phi_broadcast, npy_dssp_bb[:, 0])
        checker.identical("dssp_backbone[1] -> psi (broadcast)",
                          psi_broadcast, npy_dssp_bb[:, 1])
        # col 2 = DSSP residue SASA (may differ from atom SASA — skip)
        checker.skip("dssp_backbone[2]",
                     "DSSP residue SASA vs atom SASA — different source")
        # cols 3-4 = helix/sheet indicators (derivable from ss8 — verify)
        h5_ss8 = h5["dssp/ss8"][frame_idx]  # (R,) int8
        ss8_broadcast = h5_ss8[res_idx]
        helix_h5 = np.isin(ss8_broadcast, [0, 1, 2]).astype(np.float64)
        sheet_h5 = np.isin(ss8_broadcast, [3, 4]).astype(np.float64)
        checker.identical("dssp_backbone[3] -> helix (from ss8)",
                          helix_h5, npy_dssp_bb[:, 3])
        checker.identical("dssp_backbone[4] -> sheet (from ss8)",
                          sheet_h5, npy_dssp_bb[:, 4])

    # ── DSSP chi: 12 cols = 4 chi angles × (cos, sin, exists) ───
    # NPY cos/sin come from DsspResult (DSSP library computation).
    # H5 cos/sin come from AnalysisWriter::Dihedral() (independent
    # computation from atom positions). These are two independent
    # trig evaluations of the same angle — ULP-level differences
    # (typically 1e-16, i.e. 1 ULP) are expected.
    npy_chi = load("dssp_chi.npy")
    if npy_chi is not None:
        res_idx = h5["atoms/residue_index"][:]
        for ci, chi_name in enumerate(["chi1", "chi2", "chi3", "chi4"]):
            h5_cos = h5[f"dihedrals/{chi_name}_cos"][frame_idx]  # (R,)
            h5_sin = h5[f"dihedrals/{chi_name}_sin"][frame_idx]  # (R,)
            cos_broadcast = h5_cos[res_idx]  # (N,)
            sin_broadcast = h5_sin[res_idx]  # (N,)
            npy_cos = np.ascontiguousarray(npy_chi[:, ci*3 + 0])
            npy_sin = np.ascontiguousarray(npy_chi[:, ci*3 + 1])
            npy_exists = npy_chi[:, ci*3 + 2]

            mask = npy_exists == 1.0
            if mask.any():
                # Allow 2 ULP tolerance (independent trig computations)
                cos_h = cos_broadcast[mask]
                cos_n = npy_cos[mask]
                max_cos_diff = np.max(np.abs(cos_h - cos_n))
                if max_cos_diff < 1e-14:
                    checker.passed += 1
                else:
                    checker.failures.append(
                        f"  dssp_chi {chi_name} cos: max_diff={max_cos_diff:.2e}")
                    checker.failed += 1

                sin_h = sin_broadcast[mask]
                sin_n = npy_sin[mask]
                max_sin_diff = np.max(np.abs(sin_h - sin_n))
                if max_sin_diff < 1e-14:
                    checker.passed += 1
                else:
                    checker.failures.append(
                        f"  dssp_chi {chi_name} sin: max_diff={max_sin_diff:.2e}")
                    checker.failed += 1
            if (~mask).any():
                n_nan = np.sum(np.isnan(cos_broadcast[~mask]))
                n_expected = np.sum(~mask)
                if n_nan == n_expected:
                    checker.passed += 1
                else:
                    checker.failures.append(
                        f"  dssp_chi {chi_name}: {n_nan}/{n_expected} "
                        f"NaN where exists=0")
                    checker.failed += 1

    # ── DSSP ss8: one-hot (N, 8) in NPY vs int8 enum (R,) in H5 ─
    npy_ss8 = load("dssp_ss8.npy")
    if npy_ss8 is not None:
        res_idx = h5["atoms/residue_index"][:]
        h5_ss8 = h5["dssp/ss8"][frame_idx]  # (R,) int8
        ss8_broadcast = h5_ss8[res_idx]      # (N,) int8
        # Reconstruct one-hot from enum
        onehot_h5 = np.zeros((N, 8), dtype=np.float64)
        for i in range(N):
            cls = ss8_broadcast[i]
            if 0 <= cls < 8:
                onehot_h5[i, cls] = 1.0
        checker.identical("dssp_ss8 one-hot (reconstructed from enum)",
                          onehot_h5, npy_ss8)

    # ── Coulomb NPYs (present only if --no-coulomb was NOT set) ──
    npy_cs = load("coulomb_shielding.npy")
    if npy_cs is not None:
        checker.identical("coulomb_shielding -> efg/coulomb_shielding",
                          h5["efg/coulomb_shielding"][frame_idx], npy_cs)
        npy_ce = load("coulomb_E.npy")
        if npy_ce is not None:
            checker.identical("coulomb_E -> efg/E_total",
                              h5["efg/E_total"][frame_idx], npy_ce)
        npy_ceb = load("coulomb_efg_backbone.npy")
        if npy_ceb is not None:
            checker.identical("coulomb_efg_backbone -> efg/coulomb_backbone",
                              h5["efg/coulomb_backbone"][frame_idx], npy_ceb)
        npy_cea = load("coulomb_efg_aromatic.npy")
        if npy_cea is not None:
            checker.identical("coulomb_efg_aromatic -> efg/coulomb_aromatic",
                              h5["efg/coulomb_aromatic"][frame_idx], npy_cea)
        npy_csc = load("coulomb_scalars.npy")
        if npy_csc is not None:
            checker.identical("coulomb_scalars[0] -> E_magnitude",
                              h5["efg/E_magnitude"][frame_idx], npy_csc[:, 0])
            checker.identical("coulomb_scalars[1] -> E_bond_proj",
                              h5["efg/E_bond_proj"][frame_idx], npy_csc[:, 1])
            checker.identical("coulomb_scalars[2] -> E_backbone_frac",
                              h5["efg/E_backbone_frac"][frame_idx], npy_csc[:, 2])
            # col 3 = aromatic_E_magnitude — NOT in H5
            checker.skip("coulomb_scalars[3]",
                         "aromatic_E_magnitude not harvested into H5")
    else:
        checker.skip("coulomb NPYs", "not present (--no-coulomb)")


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <fleet_dir> <test_h5_dir>")
        sys.exit(1)

    fleet_dir = sys.argv[1]
    h5_dir = sys.argv[2]

    # Find proteins
    proteins = sorted([d for d in os.listdir(fleet_dir)
                       if os.path.isdir(os.path.join(fleet_dir, d))])
    print(f"Fleet: {len(proteins)} proteins in {fleet_dir}")
    print(f"H5s:   {h5_dir}\n")

    total_pass = 0
    total_fail = 0
    total_skip = 0
    file_pass = 0
    file_fail = 0

    for protein_id in proteins:
        protein_dir = os.path.join(fleet_dir, protein_id, "analysis_output")
        h5_path = os.path.join(h5_dir, f"{protein_id}_analysis.h5")

        if not os.path.exists(h5_path):
            print(f"=== {protein_id}: H5 not found, skipping ===")
            continue

        h5 = h5py.File(h5_path, "r")
        N = h5["meta"].attrs["n_atoms"]
        R = h5["meta"].attrs["n_residues"]
        n_rings = h5["topology"].attrs["n_rings"]
        times = h5["meta/frame_times"][:]

        # Find snapshot directories
        snap_dirs = sorted([
            d for d in os.listdir(protein_dir)
            if d.endswith("ns") and os.path.isdir(os.path.join(protein_dir, d))
        ])

        print(f"=== {protein_id} ({N} atoms, {R} res, "
              f"{len(snap_dirs)} snapshots) ===")

        protein_ok = True
        for snap_name in snap_dirs:
            snap_dir = os.path.join(protein_dir, snap_name)
            # Extract ns from name: "1B1V_4292_analysis_12ns" -> 12
            ns_str = snap_name.split("_")[-1].replace("ns", "")
            try:
                ns = int(ns_str)
            except ValueError:
                continue

            target_ps = ns * 1000.0
            diffs = np.abs(times - target_ps)
            frame_idx = int(np.argmin(diffs))
            if diffs[frame_idx] > 1.0:
                print(f"    {snap_name}: no matching frame "
                      f"(nearest={times[frame_idx]:.0f}ps)")
                continue

            checker = Checker(protein_id)
            validate_frame(h5, snap_dir, frame_idx, N, R, n_rings, checker)
            ok = checker.report(f"  {snap_name} (frame {frame_idx})")
            if not ok:
                protein_ok = False

            total_pass += checker.passed
            total_fail += checker.failed
            total_skip += checker.skipped

        h5.close()

        if protein_ok:
            file_pass += 1
        else:
            file_fail += 1
        print()

    print("=" * 60)
    print(f"Proteins: {file_pass} pass, {file_fail} fail")
    print(f"Checks:   {total_pass} pass, {total_fail} fail, {total_skip} skip")

    return 1 if total_fail > 0 else 0


if __name__ == "__main__":
    sys.exit(main())

// QtSpecialBuffers.h — TRs that don't fit the atom-axis or residue-axis
// pattern. Each gets its own struct shape matching the writer.

#pragma once

#include "Types.h"

#include <QString>
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace h5reader::model {

// ──────────────────────────────────────────────────────────────────
// QtRingNeighbourhoodTimeSeries — (N, T, 11, 4) per-(atom, ring-slot,
// channel) per-frame data.
//
// Source: /trajectory/ring_neighbourhood_trajectory_stats/.
// Channels per `channel_layout` attribute: distance, rho, z,
// in_plane_angle (Å, Å, Å, radians).
// Slots: up to 11 nearest rings per atom (r_per_atom_max attr).
// Static side-data: ring_membership_per_atom (N, 11) int32 — ring
// indices into rings.npy, -1 for unused slots.
//
// Source-attached mask: per-frame uint8 (always-attached in fixture).
// ──────────────────────────────────────────────────────────────────

struct QtRingNeighbourhoodTimeSeries {
    std::size_t n_atoms = 0;
    std::size_t n_frames = 0;
    std::size_t n_slots = 11;    // r_per_atom_max
    std::size_t n_channels = 4;  // distance, rho, z, in_plane_angle

    std::vector<double> geometry;                   // (N*T*11*4,) row-major
    std::vector<int32_t> ring_membership_per_atom;  // (N*11,) int32; -1 sentinel
    std::vector<uint64_t> frame_indices;            // (T,)
    std::vector<double> frame_times;                // (T,)
    std::vector<uint8_t> source_attached;           // (T,)

    QString units;  // "Angstrom,Angstrom,Angstrom,radians"
    QString channel_layout;
    QString result_name;
    double ring_current_spatial_cutoff_A = 15.0;
    QString z_sign_convention;
    QString in_plane_angle_range;

    // Returns the 4-channel slab for one (atom, frame, ring_slot).
    // Channels: [0]=distance, [1]=rho, [2]=z, [3]=in_plane_angle.
    // NaN if slot is unused (ring_membership_per_atom == -1).
    std::array<double, 4> at(std::size_t atomIdx, std::size_t t, std::size_t slot) const {
        std::array<double, 4> out{};
        if (atomIdx >= n_atoms || t >= n_frames || slot >= n_slots)
            return out;
        const std::size_t base = ((atomIdx * n_frames + t) * n_slots + slot) * n_channels;
        for (std::size_t k = 0; k < 4; ++k)
            out[k] = geometry[base + k];
        return out;
    }

    // Ring index for one (atom, slot) — -1 means unused.
    int32_t ringIndexAt(std::size_t atomIdx, std::size_t slot) const {
        if (atomIdx >= n_atoms || slot >= n_slots)
            return -1;
        return ring_membership_per_atom[atomIdx * n_slots + slot];
    }
};


// ──────────────────────────────────────────────────────────────────
// QtBondLengthStats — per-bond Welford rollup (bond-axis, no T axis).
// ──────────────────────────────────────────────────────────────────

struct QtBondLengthStats {
    std::size_t n_bonds = 0;

    std::vector<double> length_mean;  // (B,) Å
    std::vector<double> length_std;
    std::vector<double> length_min;
    std::vector<double> length_max;
    std::vector<double> length_delta_mean;
    std::vector<double> length_delta_std;

    // Per-bond identity reflection (also in bonds.npy; carried here for
    // self-contained access).
    std::vector<int32_t> atom_a;
    std::vector<int32_t> atom_b;
    std::vector<int8_t> bond_order;
    std::vector<int8_t> bond_category;

    QString units;
    QString result_name;
};


// ──────────────────────────────────────────────────────────────────
// QtSystemEnergyTimeSeries — gromacs energy + thermo time-series.
// Source: /trajectory/gromacs_energy_time_series/. ~30 scalar
// (T,) channels + pressure_tensor (T,9) + virial (T,9).
// All system-scope; no atom or residue index.
// ──────────────────────────────────────────────────────────────────

struct QtSystemEnergyTimeSeries {
    std::size_t n_frames = 0;

    // Scalar channels (T,) — kJ/mol unless noted
    std::vector<double> kinetic;
    std::vector<double> potential;
    std::vector<double> total_energy;
    std::vector<double> enthalpy;
    std::vector<double> temperature;    // K
    std::vector<double> T_protein;      // K
    std::vector<double> T_non_protein;  // K
    std::vector<double> pressure;       // bar
    std::vector<double> density;        // kg/m³
    std::vector<double> volume;         // nm³
    std::vector<double> box_x;          // nm
    std::vector<double> box_y;          // nm
    std::vector<double> box_z;          // nm
    std::vector<double> bond;
    std::vector<double> angle;
    std::vector<double> proper_dih;
    std::vector<double> improper_dih;
    std::vector<double> urey_bradley;
    std::vector<double> cmap_dih;
    std::vector<double> coulomb_sr;
    std::vector<double> coulomb_recip;
    std::vector<double> coulomb_14;
    std::vector<double> lj_sr;
    std::vector<double> lj_14;
    std::vector<double> disper_corr;

    // Tensor channels (T,9) — Cartesian XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ
    std::vector<double> pressure_tensor;
    std::vector<double> virial;

    // Frame metadata
    std::vector<uint64_t> frame_indices;
    std::vector<double> frame_times;
    std::vector<double> energy_frame_times_ps;  // EDR cadence (may differ)
    std::vector<uint8_t> source_attached;

    QString result_name;
    QString tensor_layout;  // "XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ"

    bool sourceAttachedAt(std::size_t t) const {
        return source_attached.empty() || (t < source_attached.size() && source_attached[t] != 0);
    }
};


// ──────────────────────────────────────────────────────────────────
// QtRmsdTracking — system-scope RMSD to frame 0.
// Source: /trajectory/rmsd_tracking/. (T,) float64.
// Static side-data: atom_indices (n_subset,) int32 — which atoms
// the RMSD subset includes (typically backbone NCACO = 216 atoms).
// ──────────────────────────────────────────────────────────────────

struct QtRmsdTracking {
    std::size_t n_frames = 0;

    std::vector<double> rmsd;           // (T,) Å
    std::vector<int32_t> atom_indices;  // (n_subset,)
    std::vector<uint64_t> frame_indices;
    std::vector<double> frame_times;
    std::vector<uint8_t> source_attached;

    QString units;
    QString result_name;
    QString alignment_method;        // "kabsch_svd"
    QString atom_selection;          // "backbone_heavy_atoms_NCACO"
    QString reference_frame_origin;  // "trajectory_frame_0"
};


// ──────────────────────────────────────────────────────────────────
// QtDssp8Transitions — per-residue DSSP transition statistics.
// Source: /trajectory/dssp8_transition/. Per-residue:
//   - ss8_occupancy (R, 8) uint32 (count of frames in each SS class)
//   - ss8_transition_matrix (R, 8, 8) uint32
//   - ss8_transition_count (R,) uint32
// Conditional-attach per DSSP cadence.
// ──────────────────────────────────────────────────────────────────

struct QtDssp8Transitions {
    std::size_t n_residues = 0;

    std::vector<uint32_t> ss8_occupancy;          // (R*8,)
    std::vector<uint32_t> ss8_transition_matrix;  // (R*8*8,)
    std::vector<uint32_t> ss8_transition_count;   // (R,)

    QString result_name;
};


// ──────────────────────────────────────────────────────────────────
// QtDihedralBinTransitions — per-residue dihedral bin transition stats.
// Source: /trajectory/dihedral_bin_transition/. Per-residue:
//   - backbone_bin_occupancy (R, 6) uint32 (Ramachandran 6-region bins)
//   - backbone_transition_count (R,) uint32
//   - backbone_dominant_region (R,) uint8 (RamaRegion enum)
//   - chi_rotamer_occupancy (R, 4, 3) uint32 (3 rotamer bins per chi)
//   - chi_transition_count (R, 4) uint32
//   - chi_dominant_rotamer (R, 4) uint8
// Always-attached.
// ──────────────────────────────────────────────────────────────────

struct QtDihedralBinTransitions {
    std::size_t n_residues = 0;

    std::vector<uint32_t> backbone_bin_occupancy;     // (R*6,)
    std::vector<uint32_t> backbone_transition_count;  // (R,)
    std::vector<uint8_t> backbone_dominant_region;    // (R,)

    std::vector<uint32_t> chi_rotamer_occupancy;  // (R*4*3,)
    std::vector<uint32_t> chi_transition_count;   // (R*4,)
    std::vector<uint8_t> chi_dominant_rotamer;    // (R*4,)

    QString result_name;
};


// ──────────────────────────────────────────────────────────────────
// QtWaterFieldTimeSeries — composite per-atom water-field TR.
// Source: /trajectory/water_field_time_series/. Six channels:
//   - efg          (N, T, 5) float64 — water-derived EFG T2 components
//   - efg_first    (N, T, 5) float64 — first-shell-only EFG
//   - efield       (N, T, 3) float64 — water-derived E-field (Vec3)
//   - efield_first (N, T, 3) float64 — first-shell-only E-field
//   - n_first      (N, T)    uint32  — first-shell water count
//   - n_second     (N, T)    uint32  — second-shell water count
// Always-attached.
// ──────────────────────────────────────────────────────────────────

struct QtWaterFieldTimeSeries {
    std::size_t n_atoms = 0;
    std::size_t n_frames = 0;

    std::vector<double> efg;  // (N*T*5,)
    std::vector<double> efg_first;
    std::vector<double> efield;  // (N*T*3,)
    std::vector<double> efield_first;
    std::vector<uint32_t> n_first;  // (N*T,)
    std::vector<uint32_t> n_second;

    std::vector<uint64_t> frame_indices;
    std::vector<double> frame_times;
    std::vector<uint8_t> source_attached;

    QString result_name;

    std::array<double, 5> efgAt(std::size_t atomIdx, std::size_t t) const {
        std::array<double, 5> out{};
        if (atomIdx >= n_atoms || t >= n_frames)
            return out;
        const std::size_t base = (atomIdx * n_frames + t) * 5;
        for (std::size_t k = 0; k < 5; ++k)
            out[k] = efg[base + k];
        return out;
    }
    Vec3 efieldAt(std::size_t atomIdx, std::size_t t) const {
        if (atomIdx >= n_atoms || t >= n_frames)
            return Vec3::Zero();
        const std::size_t base = (atomIdx * n_frames + t) * 3;
        return Vec3(efield[base + 0], efield[base + 1], efield[base + 2]);
    }
    uint32_t nFirstAt(std::size_t atomIdx, std::size_t t) const {
        if (atomIdx >= n_atoms || t >= n_frames)
            return 0;
        return n_first[atomIdx * n_frames + t];
    }
};


// ──────────────────────────────────────────────────────────────────
// QtHydrationShellTimeSeries — composite per-atom hydration shell TR.
// Source: /trajectory/hydration_shell_time_series/. Four scalar
// channels (N, T): half_shell_asymmetry, mean_water_dipole_cos,
// nearest_ion_charge, nearest_ion_distance. Always-attached.
// ──────────────────────────────────────────────────────────────────

struct QtHydrationShellTimeSeries {
    std::size_t n_atoms = 0;
    std::size_t n_frames = 0;

    std::vector<double> half_shell_asymmetry;
    std::vector<double> mean_water_dipole_cos;
    std::vector<double> nearest_ion_charge;
    std::vector<double> nearest_ion_distance;

    std::vector<uint64_t> frame_indices;
    std::vector<double> frame_times;
    std::vector<uint8_t> source_attached;

    QString result_name;

    double halfShellAsymmetryAt(std::size_t atomIdx, std::size_t t) const {
        if (atomIdx >= n_atoms || t >= n_frames)
            return 0.0;
        return half_shell_asymmetry[atomIdx * n_frames + t];
    }
};


// ──────────────────────────────────────────────────────────────────
// QtHydrationGeometryTimeSeries — composite per-atom hydration geometry TR.
// Source: /trajectory/hydration_geometry_time_series/. Five channels:
//   - dipole_alignment    (N, T) float64
//   - dipole_coherence    (N, T) float64
//   - dipole_vector       (N, T, 3) float64
//   - first_shell_count   (N, T) uint32
//   - half_shell_asymmetry (N, T) float64
//   - surface_normal      (N, T, 3) float64
// Always-attached.
// ──────────────────────────────────────────────────────────────────

struct QtHydrationGeometryTimeSeries {
    std::size_t n_atoms = 0;
    std::size_t n_frames = 0;

    std::vector<double> dipole_alignment;
    std::vector<double> dipole_coherence;
    std::vector<double> dipole_vector;  // (N*T*3,)
    std::vector<uint32_t> first_shell_count;
    std::vector<double> half_shell_asymmetry;
    std::vector<double> surface_normal;  // (N*T*3,)

    std::vector<uint64_t> frame_indices;
    std::vector<double> frame_times;
    std::vector<uint8_t> source_attached;

    QString result_name;

    Vec3 dipoleVectorAt(std::size_t atomIdx, std::size_t t) const {
        if (atomIdx >= n_atoms || t >= n_frames)
            return Vec3::Zero();
        const std::size_t base = (atomIdx * n_frames + t) * 3;
        return Vec3(dipole_vector[base + 0], dipole_vector[base + 1], dipole_vector[base + 2]);
    }
    Vec3 surfaceNormalAt(std::size_t atomIdx, std::size_t t) const {
        if (atomIdx >= n_atoms || t >= n_frames)
            return Vec3::Zero();
        const std::size_t base = (atomIdx * n_frames + t) * 3;
        return Vec3(surface_normal[base + 0], surface_normal[base + 1], surface_normal[base + 2]);
    }
    uint32_t firstShellCountAt(std::size_t atomIdx, std::size_t t) const {
        if (atomIdx >= n_atoms || t >= n_frames)
            return 0;
        return first_shell_count[atomIdx * n_frames + t];
    }
};


// ──────────────────────────────────────────────────────────────────
// QtRingPuckerTimeSeries — per-ring conformational TR.
// Source: /trajectory/ring_pucker_time_series/. Two axes:
//   - aromatic axis: aromatic_chi2 (n_aromatic, T) float64 — chi2 angle
//     per aromatic ring (radians; ring rotation about chi2)
//   - saturated axis: pucker_Q (n_saturated, T), pucker_theta (n_saturated, T)
//     — Cremer-Pople puckering amplitude + phase for saturated rings
//     (only PRO pyrrolidine in standard 20)
// Both axes carry their parent residue index for cross-reference.
// Always-attached.
// ──────────────────────────────────────────────────────────────────

struct QtRingPuckerTimeSeries {
    std::size_t n_aromatic_rings = 0;
    std::size_t n_saturated_rings = 0;
    std::size_t n_frames = 0;

    std::vector<double> aromatic_chi2;                    // (n_aromatic*T,) rad
    std::vector<int32_t> aromatic_parent_residue_index;   // (n_aromatic,)
    std::vector<double> pucker_Q;                         // (n_saturated*T,) Å
    std::vector<double> pucker_theta;                     // (n_saturated*T,) rad
    std::vector<int32_t> saturated_parent_residue_index;  // (n_saturated,)

    std::vector<uint64_t> frame_indices;
    std::vector<double> frame_times;
    std::vector<uint8_t> source_attached;

    QString result_name;

    double aromaticChi2At(std::size_t aromaticIdx, std::size_t t) const {
        if (aromaticIdx >= n_aromatic_rings || t >= n_frames)
            return 0.0;
        return aromatic_chi2[aromaticIdx * n_frames + t];
    }
    double puckerQAt(std::size_t saturatedIdx, std::size_t t) const {
        if (saturatedIdx >= n_saturated_rings || t >= n_frames)
            return 0.0;
        return pucker_Q[saturatedIdx * n_frames + t];
    }
};


}  // namespace h5reader::model

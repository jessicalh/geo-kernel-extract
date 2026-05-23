#pragma once
//
// TrajectoryH5: typed read boundary for the per-TR-emits-its-own-group
// trajectory.h5 format emitted by `nmr_extract --trajectory`.
//
// One HDF5 group per attached TrajectoryResult under
// `/trajectory/<name>/`, with per-TR schema documented at
// OBJECT_MODEL.md "H5 file layout".
//
// Shape: the constructor opens the file, validates the structural
// minimum (`/atoms`, `/trajectory/frames`, root attrs), eagerly reads
// the bounded slices the viewer needs (Welford rollups at the atom
// axis + frame-0 slabs of each per-atom time series), then closes the
// file. Owns no HighFive handles after construction.
//
// Each per-TR accessor returns `std::optional<>`; absent groups are
// normal (sparse TR sets — May-8 fixtures have ~6 of ~24 TR groups).
// There is no generic `Get(name)` adapter surface: per the PATTERNS
// discipline, the typed object answers typed questions.
//
// Construction throws `HighFive::Exception` on structural failure;
// ComputeWorker is responsible for catching and converting to a log
// line (the HighFive boundary lives entirely in this class).
//
// Note on T1 basis (deferred): per OBJECT_MODEL.md, runtime
// `SphericalTensor.T1` carries Cartesian (v_x, v_y, v_z) while the
// file's `irrep_layout` attribute names m-basis. T0 and |T2| are
// well-defined regardless. The MVP exposes only T0 and |T2| for
// shielding time series; T1 surfacing waits on a basis resolution.

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

class TrajectoryH5 {
public:
    // Mean + std rollup pair, the headline scalars per metric.
    struct WelfordPair {
        double mean = 0.0;
        double std = 0.0;
    };

    // Tensor-source Welford rollup: T0 (isotropic) + |T2| (rms of the
    // 5 m-components, computed as sqrt(sum_{m=-2..+2} T2_m^2)). Shared
    // shape across all shielding calculators (BS, HM, McConnell,
    // PiQuad, RingChi, Disp, HBond kernel-form). The accessor method
    // name carries the physics source; the row type carries only the
    // common shape.
    struct ShieldingWelfordRow {
        WelfordPair t0;
        WelfordPair t2magnitude;  // mean/std of |T2| samples
    };

    // Scalar-source Welford rollups — distinct channel names per
    // physics so the call site reads naturally.
    struct SasaWelfordRow {
        WelfordPair sasa;
    };
    struct EeqWelfordRow {
        WelfordPair charge;
    };
    struct HBondCountWelfordRow {
        WelfordPair count;
    };

    // Tensor-source frame-0 slab. T0 + |T2| only; T1 deferred until
    // the file-vs-runtime basis question is resolved (the file emits
    // m-basis under irrep_layout, runtime SphericalTensor.T1 carries
    // Cartesian — see header doc above). T2_magnitude is
    // sqrt(sum_{m=-2..+2} T2_m^2), same convention as the Welford
    // t2magnitude channel above.
    struct ShieldingFrame0Row {
        double T0 = 0.0;
        double T2_magnitude = 0.0;
    };

    // Cartesian Vec3 frame-0 row (APBS efield etc.). Deliberately
    // NOT `nmr::Vec3` (which is `Eigen::Vector3d`): keeping `Types.h`
    // out of this header avoids pulling Eigen into every viewer
    // translation unit that includes TrajectoryH5.h. The inspector
    // and REST sites compute their own magnitudes from .x/.y/.z;
    // there's nothing to convert.
    struct Vec3F0 {
        double x = 0.0, y = 0.0, z = 0.0;
    };

    // Throws HighFive::Exception on structural failure.
    explicit TrajectoryH5(const std::string& path);

    // --- Top-level ----------------------------------------------------
    std::size_t AtomCount() const { return n_atoms_; }
    std::size_t FrameCount() const { return n_frames_; }
    const std::string& ProteinId() const { return protein_id_; }

    // --- /atoms identity (for cross-check against library Protein) ---
    int ElementAt(std::size_t i) const { return atom_element_[i]; }
    const std::string& AtomNameAt(std::size_t i) const { return atom_name_[i]; }
    std::size_t ResidueIndexAt(std::size_t i) const { return residue_index_[i]; }

    // --- /trajectory/frames ------------------------------------------
    double FrameTimePs(std::size_t t) const { return frame_times_[t]; }

    // --- /trajectory attrs -------------------------------------------
    const std::string& XtcPath() const { return xtc_path_; }
    const std::string& TprPath() const { return tpr_path_; }
    const std::string& EdrPath() const { return edr_path_; }

    // --- Per-Welford rollups at the atom axis ------------------------
    std::optional<ShieldingWelfordRow> BsWelford(std::size_t atom_idx) const;
    std::optional<ShieldingWelfordRow> HmWelford(std::size_t atom_idx) const;
    std::optional<ShieldingWelfordRow> McWelford(std::size_t atom_idx) const;
    std::optional<SasaWelfordRow> SasaWelford(std::size_t atom_idx) const;
    std::optional<EeqWelfordRow> EeqWelford(std::size_t atom_idx) const;
    std::optional<HBondCountWelfordRow> HBondCountWelford(std::size_t atom_idx) const;

    // --- Per-TS frame-0 slabs ----------------------------------------
    // Tensor shielding (xyz[N, T, 9] frame-0 → T0 + |T2|).
    std::optional<ShieldingFrame0Row> BsShieldingFrame0(std::size_t atom_idx) const;
    std::optional<ShieldingFrame0Row> HmShieldingFrame0(std::size_t atom_idx) const;
    std::optional<ShieldingFrame0Row> McShieldingFrame0(std::size_t atom_idx) const;
    std::optional<ShieldingFrame0Row> PiQuadShieldingFrame0(std::size_t atom_idx) const;
    std::optional<ShieldingFrame0Row> RingChiShieldingFrame0(std::size_t atom_idx) const;
    std::optional<ShieldingFrame0Row> DispShieldingFrame0(std::size_t atom_idx) const;
    std::optional<ShieldingFrame0Row> HBondShieldingFrame0(std::size_t atom_idx) const;

    // Scalar TS (xyz[N, T] frame 0).
    std::optional<double> SasaFrame0(std::size_t atom_idx) const;
    std::optional<double> Aimnet2ChargeFrame0(std::size_t atom_idx) const;

    // Vec3 TS (xyz[N, T, 3] frame 0).
    std::optional<Vec3F0> ApbsEfieldFrame0(std::size_t atom_idx) const;

    // --- Sparse-set introspection ------------------------------------
    bool HasGroup(const std::string& name) const;
    const std::vector<std::string>& GroupsPresent() const { return groups_present_; }

private:
    // Top-level
    std::size_t n_atoms_ = 0;
    std::size_t n_frames_ = 0;
    std::string protein_id_;

    // /atoms
    std::vector<int> atom_element_;
    std::vector<std::string> atom_name_;
    std::vector<std::size_t> residue_index_;

    // /trajectory/frames
    std::vector<double> frame_times_;

    // /trajectory attrs (optional — populated when present)
    std::string xtc_path_, tpr_path_, edr_path_;

    // Inventory: every direct child of /trajectory present in this file.
    std::vector<std::string> groups_present_;

    // Per-Welford rollup buffers. outer optional<> = group absent;
    // inner vector sized n_atoms_ when present.
    std::optional<std::vector<ShieldingWelfordRow>> bs_welford_;
    std::optional<std::vector<ShieldingWelfordRow>> hm_welford_;
    std::optional<std::vector<ShieldingWelfordRow>> mc_welford_;
    std::optional<std::vector<SasaWelfordRow>> sasa_welford_;
    std::optional<std::vector<EeqWelfordRow>> eeq_welford_;
    std::optional<std::vector<HBondCountWelfordRow>> hbond_count_welford_;

    // Per-TS frame-0 slabs (atom-axis, N-sized).
    std::optional<std::vector<ShieldingFrame0Row>> bs_shielding_f0_;
    std::optional<std::vector<ShieldingFrame0Row>> hm_shielding_f0_;
    std::optional<std::vector<ShieldingFrame0Row>> mc_shielding_f0_;
    std::optional<std::vector<ShieldingFrame0Row>> piquad_shielding_f0_;
    std::optional<std::vector<ShieldingFrame0Row>> ringchi_shielding_f0_;
    std::optional<std::vector<ShieldingFrame0Row>> disp_shielding_f0_;
    std::optional<std::vector<ShieldingFrame0Row>> hbond_shielding_f0_;
    std::optional<std::vector<double>> sasa_f0_;
    std::optional<std::vector<double>> aimnet2_charge_f0_;
    std::optional<std::vector<Vec3F0>> apbs_efield_f0_;
};

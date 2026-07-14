#include "TripeptideNeighborShieldingTimeSeriesTrajectoryResult.h"
#include "TripeptideNeighborShieldingResult.h"
#include "TrajectoryProtein.h"
#include "ProteinConformation.h"
#include "ConformationAtom.h"
#include "OperationLog.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <limits>
#include <typeinfo>

namespace nmr {


std::unique_ptr<TripeptideNeighborShieldingTimeSeriesTrajectoryResult>
TripeptideNeighborShieldingTimeSeriesTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<TripeptideNeighborShieldingTimeSeriesTrajectoryResult>();
    r->per_atom_shielding_.assign(tp.AtomCount(),
                                  std::vector<SphericalTensor>{});
    r->per_atom_has_match_.assign(tp.AtomCount(),
                                  std::vector<std::uint8_t>{});
    return r;
}


// ── Compute ──────────────────────────────────────────────────────
//
// Append this frame's per-atom SphericalTensor to the growing buffers.
// Records frame_idx + time_ps so WriteH5Group emits the frame list
// that downstream consumers need to align columns with time.

void TripeptideNeighborShieldingTimeSeriesTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;
    // Record both frame-level source presence and atom-level contribution
    // applicability. The latter distinguishes an unmatched zero-default
    // tensor from a matched contribution whose physical sum is zero.
    const bool source_attached = force_source_present_for_testing_
        || conf.HasResult<TripeptideNeighborShieldingResult>();
    source_present_per_frame_.push_back(source_attached ? 1u : 0u);

    const std::size_t N = conf.AtomCount();
    for (std::size_t i = 0; i < N; ++i) {
        const ConformationAtom& atom = conf.AtomAt(i);
        per_atom_shielding_[i].push_back(
            atom.tripeptide_neighbor_shielding_spherical);
        per_atom_has_match_[i].push_back(
            atom.tripeptide_neighbor_has_match ? 1u : 0u);
    }
    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    ++n_frames_;
}


// ── Finalize ─────────────────────────────────────────────────────
//
// Transfer growing per-atom buffers into a contiguous atom-major
// DenseBuffer<SphericalTensor> owned by TrajectoryProtein.

// Idempotent: the per-atom bounds check (src.size() != n_frames_)
// makes a second call short-circuit naturally — after the first
// Finalize swaps per_atom_shielding_[i] to empty, src.size() is 0,
// the bounds check skips that atom, the buffer comes out empty, and
// we don't AdoptDenseBuffer. The data flow makes idempotency fall
// out without a "have I run yet" state flag.
void TripeptideNeighborShieldingTimeSeriesTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                                    Trajectory& traj) {
    (void)traj;
    const std::size_t N = tp.AtomCount();

    auto buffer =
        std::make_unique<DenseBuffer<SphericalTensor>>(N, n_frames_);
    std::size_t atoms_written = 0;
    for (std::size_t i = 0; i < N; ++i) {
        const auto& src = per_atom_shielding_[i];
        if (src.size() != n_frames_) continue;  // skip mismatched / post-swap
        for (std::size_t f = 0; f < n_frames_; ++f) {
            buffer->At(i, f) = src[f];
        }
        std::vector<SphericalTensor>().swap(per_atom_shielding_[i]);
        ++atoms_written;
    }

    if (atoms_written == 0) return;  // idempotent no-op on second call

    tp.AdoptDenseBuffer<SphericalTensor>(
        std::move(buffer),
        std::type_index(typeid(TripeptideNeighborShieldingTimeSeriesTrajectoryResult)));

    finalized_ = true;

    OperationLog::Info(LogCalcOther,
        "TripeptideNeighborShieldingTimeSeriesTrajectoryResult::Finalize",
        "transferred (" + std::to_string(N) + " atoms x " +
        std::to_string(n_frames_) +
        " frames) SphericalTensor time-series to tp dense buffer");
}


// ── WriteH5Group ─────────────────────────────────────────────────
//
// Flat (N · T · 9) double array via explicit component access on each
// SphericalTensor — .T0 / .T1[k] / .T2[k] — so the layout is
// independent of any assumption about struct packing. The 9-component
// trailing axis follows SphericalTensor::PackFull9:
//
//   index 0: T0   (l=0, m=0)
//   index 1: T1_x   (Cartesian Levi-Civita dual)
//   index 2: T1_y
//   index 3: T1_z
//   index 4: T2_{-2}
//   index 5: T2_{-1}
//   index 6: T2_{0}
//   index 7: T2_{+1}
//   index 8: T2_{+2}
//
// The payload order above is explicit; group attributes also record the
// proper-rotation-only transformation scope of the chiral lookup.

void TripeptideNeighborShieldingTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    auto* buffer = const_cast<TrajectoryProtein&>(tp)
        .GetDenseBuffer<SphericalTensor>(std::type_index(
            typeid(TripeptideNeighborShieldingTimeSeriesTrajectoryResult)));
    if (!buffer) {
        OperationLog::Warn(
            "TripeptideNeighborShieldingTimeSeriesTrajectoryResult::WriteH5Group",
            "no dense buffer present (Finalize not called?)");
        return;
    }

    // "Absent, not faked" — if no frame had the source-present flag,
    // skip emission.
    std::size_t source_present_count = 0;
    for (auto v : source_present_per_frame_)
        if (v) ++source_present_count;
    if (source_present_count == 0) {
        OperationLog::Warn(
            "TripeptideNeighborShieldingTimeSeriesTrajectoryResult::WriteH5Group",
            "TripeptideNeighborShieldingResult was not attached in any of "
            + std::to_string(source_present_per_frame_.size()) +
            " frames; skipping /trajectory/tripeptide_neighbor_shielding_time_series/ "
            "emission per 'absent, not faked' discipline.");
        return;
    }

    const std::size_t N = buffer->AtomCount();
    const std::size_t T = buffer->StridePerAtom();

    auto grp = file.createGroup("/trajectory/tripeptide_neighbor_shielding_time_series");

    // Provenance header.
    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_atoms",     N);
    grp.createAttribute("n_frames",    T);
    grp.createAttribute("finalized",   finalized_);

    // Schema metadata.
    grp.createAttribute("irrep_layout",
        std::string("T0,T1_x,T1_y,T1_z,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"));
    grp.createAttribute("normalization", std::string("isometric_real_sph"));
    grp.createAttribute("parity",        std::string("mixed"));
    grp.createAttribute("coordinate_frame",
        std::string("conformation_cartesian_xyz"));
    grp.createAttribute("tensor_basis",
        std::string("project_native_full9_spherical_tensor_v1"));
    grp.createAttribute("tensor_component_order", std::string(
        "T0,T1_x,T1_y,T1_z,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"));
    grp.createAttribute("tensor_frame",
        std::string("conformation_cartesian_xyz"));
    grp.createAttribute("tensor_t1_semantics", std::string(
        "Cartesian Levi-Civita dual x,y,z (not real-Y1m): "
        "a=((T_yz-T_zy)/2,(T_zx-T_xz)/2,(T_xy-T_yx)/2); "
        "axial a'=det(R) R a; generically nonzero"));
    grp.createAttribute("tensor_t1_structural_zero", false);
    grp.createAttribute("tensor_structural_zero_components",
        std::string("none"));
    grp.createAttribute("e3nn_export", std::string(
        "explicit project-basis to e3nn conversion required before use"));
    grp.createAttribute("normalization_scope", std::string(
        "xyz tensor payload: T2 uses isometric real-tesseral "
        "normalization; T1 uses the tensor_t1_semantics convention"));
    grp.createAttribute("transformation", std::string(
        "even_rank2 under proper rotations: T'=R T R^T; typed tripeptide "
        "lookup/Kabsch alignment is L-amino-acid chirality-conditioned and "
        "has no improper-transform contract"));
    grp.createAttribute("units",         std::string("ppm"));

    // Flat (N, T, 9) via explicit component access. NaN-fill rows where
    // the source is absent or this atom had no neighboring tripeptide
    // contribution, preserving a matched physical zero as zero.
    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> flat(N * T * 9);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t t = 0; t < T; ++t) {
            const std::size_t base = (i * T + t) * 9;
            const bool atom_matched = i < per_atom_has_match_.size()
                && t < per_atom_has_match_[i].size()
                && per_atom_has_match_[i][t] != 0;
            if (t >= source_present_per_frame_.size()
                || source_present_per_frame_[t] == 0
                || !atom_matched) {
                for (std::size_t k = 0; k < 9; ++k) flat[base + k] = kNaN;
                continue;
            }
            const SphericalTensor& st = buffer->At(i, t);
            st.PackFull9(&flat[base]);
        }
    }

    std::vector<std::size_t> dims = {N, T, std::size_t(9)};
    HighFive::DataSpace space(dims);
    auto ds = grp.createDataSet<double>("xyz", space);
    ds.write_raw(flat.data());

    grp.createDataSet("frame_indices", frame_indices_);
    grp.createDataSet("frame_times",   frame_times_);

    // Provenance mask.
    grp.createDataSet("source_attached_per_frame", source_present_per_frame_);
}

}  // namespace nmr

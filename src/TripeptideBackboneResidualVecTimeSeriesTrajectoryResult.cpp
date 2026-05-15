#include "TripeptideBackboneResidualVecTimeSeriesTrajectoryResult.h"
#include "TripeptideBackboneShieldingResult.h"
#include "TrajectoryProtein.h"
#include "ProteinConformation.h"
#include "ConformationAtom.h"
#include "OperationLog.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <typeinfo>

namespace nmr {


std::unique_ptr<TripeptideBackboneResidualVecTimeSeriesTrajectoryResult>
TripeptideBackboneResidualVecTimeSeriesTrajectoryResult::Create(
        const TrajectoryProtein& tp) {
    auto r = std::make_unique<
        TripeptideBackboneResidualVecTimeSeriesTrajectoryResult>();
    r->per_atom_residual_.assign(tp.AtomCount(), std::vector<Vec3>{});
    return r;
}


// ── Compute ──────────────────────────────────────────────────────
//
// Append this frame's per-atom Vec3 residual to the growing buffers.
// Records frame_idx + time_ps so WriteH5Group emits the frame list
// that downstream consumers need to align columns with time.

void TripeptideBackboneResidualVecTimeSeriesTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;
    const std::size_t N = conf.AtomCount();
    for (std::size_t i = 0; i < N; ++i) {
        per_atom_residual_[i].push_back(
            conf.AtomAt(i).tripeptide_bb_residual_vec);
    }
    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    ++n_frames_;
}


// ── Finalize ─────────────────────────────────────────────────────
//
// Transfer growing per-atom buffers into a contiguous atom-major
// DenseBuffer<Vec3> owned by TrajectoryProtein.
//
// Idempotent: the per-atom bounds check (src.size() != n_frames_)
// makes a second call short-circuit naturally — after the first
// Finalize swaps per_atom_residual_[i] to empty, src.size() is 0,
// the bounds check skips that atom, the buffer comes out empty, and
// we don't AdoptDenseBuffer. The data flow makes idempotency fall
// out without a "have I run yet" state flag (see
// feedback_bounds_check_over_state_flag).

void TripeptideBackboneResidualVecTimeSeriesTrajectoryResult::Finalize(
        TrajectoryProtein& tp, Trajectory& traj) {
    (void)traj;
    const std::size_t N = tp.AtomCount();

    auto buffer = std::make_unique<DenseBuffer<Vec3>>(N, n_frames_);
    std::size_t atoms_written = 0;
    for (std::size_t i = 0; i < N; ++i) {
        const auto& src = per_atom_residual_[i];
        if (src.size() != n_frames_) continue;  // skip mismatched / post-swap
        for (std::size_t f = 0; f < n_frames_; ++f) {
            buffer->At(i, f) = src[f];
        }
        std::vector<Vec3>().swap(per_atom_residual_[i]);
        ++atoms_written;
    }

    if (atoms_written == 0) return;  // idempotent no-op on second call

    tp.AdoptDenseBuffer<Vec3>(
        std::move(buffer),
        std::type_index(typeid(
            TripeptideBackboneResidualVecTimeSeriesTrajectoryResult)));

    finalized_ = true;

    OperationLog::Info(LogCalcOther,
        "TripeptideBackboneResidualVecTimeSeriesTrajectoryResult::Finalize",
        "transferred (" + std::to_string(N) + " atoms x " +
        std::to_string(n_frames_) +
        " frames) Vec3 residual time-series to tp dense buffer");
}


// ── WriteH5Group ─────────────────────────────────────────────────
//
// Flat (N · T · 3) double array via explicit .x()/.y()/.z() component
// access on each Vec3. No reinterpret_cast on the buffer — Eigen's
// Matrix<double,3,1> is 3 contiguous doubles in practice, but that's
// not a contract the library exposes. Explicit component access is
// the same pattern PositionsTimeSeriesTrajectoryResult uses.
//
// The three Cartesian doubles per atom per frame are x, y, z in the
// protein's lab frame — same frame as ConformationAtom positions.
// The (irrep_layout="x,y,z", normalization="cartesian") attribute
// pair tells downstream e3nn consumers to apply the Cartesian →
// real-spherical change of basis themselves (rather than the
// SphericalTensor TRs which emit pre-decomposed (T0, T1_m, T2_m)).

void TripeptideBackboneResidualVecTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    auto* buffer = const_cast<TrajectoryProtein&>(tp)
        .GetDenseBuffer<Vec3>(std::type_index(
            typeid(TripeptideBackboneResidualVecTimeSeriesTrajectoryResult)));
    if (!buffer) {
        OperationLog::Warn(
            "TripeptideBackboneResidualVecTimeSeriesTrajectoryResult::WriteH5Group",
            "no dense buffer present (Finalize not called?)");
        return;
    }

    const std::size_t N = buffer->AtomCount();
    const std::size_t T = buffer->StridePerAtom();

    auto grp = file.createGroup(
        "/trajectory/tripeptide_bb_residual_vec_time_series");

    // Provenance header.
    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_atoms",     N);
    grp.createAttribute("n_frames",    T);
    grp.createAttribute("finalized",   finalized_);

    // e3nn-consumable metadata. Cartesian xyz layout, polar vector
    // (parity 1o), units Å. The "angstrom" string keeps the H5
    // attribute byte-clean (ASCII); downstream consumers map it to Å.
    grp.createAttribute("irrep_layout",  std::string("x,y,z"));
    grp.createAttribute("normalization", std::string("cartesian"));
    grp.createAttribute("parity",        std::string("1o"));
    grp.createAttribute("units",         std::string("angstrom"));

    // Flat (N, T, 3) via explicit component access. Atom-major:
    // [atom_0_frame_0_xyz, atom_0_frame_1_xyz, ..., atom_1_frame_0_xyz, ...].
    std::vector<double> flat(N * T * 3);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t t = 0; t < T; ++t) {
            const Vec3& v = buffer->At(i, t);
            const std::size_t base = (i * T + t) * 3;
            flat[base + 0] = v.x();
            flat[base + 1] = v.y();
            flat[base + 2] = v.z();
        }
    }

    std::vector<std::size_t> dims = {N, T, std::size_t(3)};
    HighFive::DataSpace space(dims);
    auto ds = grp.createDataSet<double>("xyz", space);
    ds.write_raw(flat.data());

    grp.createDataSet("frame_indices", frame_indices_);
    grp.createDataSet("frame_times",   frame_times_);
}

}  // namespace nmr

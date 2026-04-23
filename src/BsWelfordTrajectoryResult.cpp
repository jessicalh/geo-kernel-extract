#include "BsWelfordTrajectoryResult.h"
#include "BiotSavartResult.h"
#include "TrajectoryProtein.h"
#include "TrajectoryAtom.h"
#include "TrajectoryMoments.h"
#include "ProteinConformation.h"
#include "ConformationAtom.h"
#include "Types.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <cmath>
#include <limits>
#include <typeinfo>

namespace nmr {


std::vector<std::type_index> BsWelfordTrajectoryResult::Dependencies() const {
    return { std::type_index(typeid(BiotSavartResult)) };
}


std::unique_ptr<BsWelfordTrajectoryResult>
BsWelfordTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto result = std::make_unique<BsWelfordTrajectoryResult>();
    const size_t N = tp.AtomCount();
    result->prev_t0_.assign(N, 0.0);
    result->prev_valid_.assign(N, false);
    return result;
}


// ── Compute ──────────────────────────────────────────────────────
//
// Welford online update per atom (mean / M2), tracking min/max with
// frame index, plus a frame-to-frame delta tracker. All writes are
// to TrajectoryAtom fields owned by this result (singleton per type,
// one writer per field).

// ── ConformationAtom → TrajectoryAtom bridge ─────────────────────
//
// This is the canonical pattern. Every per-atom TrajectoryResult
// follows the same shape:
//
//   READS (per-frame, from ConformationAtom on the current conf):
//     conf.AtomAt(i).bs_shielding_contribution.T0        (ppm, scalar)
//     conf.AtomAt(i).bs_shielding_contribution.T2Magnitude()  (ppm)
//
//   WRITES (running, onto TrajectoryAtom owned by this result):
//     ta.bs_t0_mean / bs_t0_m2 / bs_t0_min / bs_t0_max
//     ta.bs_t0_min_frame / bs_t0_max_frame / bs_n_frames
//     ta.bs_t2mag_{mean,m2,min,max,min_frame,max_frame}
//     ta.bs_t0_delta_{mean,m2,min,max,n}
//
// Singleton-per-writer is enforced by TrajectoryProtein::AttachResult
// (one instance per type). Each field above has exactly one writer:
// THIS class. No other TrajectoryResult touches these fields.
//
// Required per-frame precondition: BiotSavartResult must have run on
// `conf` before this Compute fires. Declared via Dependencies() +
// validated in Trajectory::Run Phase 2 against the RunConfiguration's
// required ConformationResult set.

void BsWelfordTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)traj; (void)time_ps;
    const size_t N = tp.AtomCount();

    for (size_t i = 0; i < N; ++i) {
        // ── READ from ConformationAtom (this frame's per-atom state) ──
        const auto& ca = conf.AtomAt(i);
        const double t0    = ca.bs_shielding_contribution.T0;
        const double t2mag = ca.bs_shielding_contribution.T2Magnitude();

        // ── WRITE through TrajectoryAtom (running across frames) ─────
        TrajectoryAtom& ta = tp.MutableAtomAt(i);
        const size_t n_new = ta.bs_n_frames + 1;

        MomentsUpdate(ta.bs_t0_mean, ta.bs_t0_m2, t0, n_new);
        MinMaxUpdate(ta.bs_t0_min, ta.bs_t0_max,
                     ta.bs_t0_min_frame, ta.bs_t0_max_frame,
                     t0, frame_idx);

        MomentsUpdate(ta.bs_t2mag_mean, ta.bs_t2mag_m2, t2mag, n_new);
        MinMaxUpdate(ta.bs_t2mag_min, ta.bs_t2mag_max,
                     ta.bs_t2mag_min_frame, ta.bs_t2mag_max_frame,
                     t2mag, frame_idx);

        ta.bs_n_frames = n_new;

        // ── Frame-to-frame T0 delta ──────────────────────────────────
        // Skip the first frame per atom; delta needs a prior value.
        if (prev_valid_[i]) {
            const double delta_x = t0 - prev_t0_[i];
            const size_t dn_new  = ta.bs_t0_delta_n + 1;
            MomentsUpdate(ta.bs_t0_delta_mean, ta.bs_t0_delta_m2,
                          delta_x, dn_new);
            MinMaxUpdate(ta.bs_t0_delta_min, ta.bs_t0_delta_max, delta_x);
            ta.bs_t0_delta_n = dn_new;
        }
        prev_t0_[i] = t0;
        prev_valid_[i] = true;
    }

    ++n_frames_;
}


// ── Finalize ─────────────────────────────────────────────────────
//
// Convert Welford M2 → unbiased standard deviation (M2 / (n-1)).
// For n ≤ 1, std is 0 (single-sample or empty).

void BsWelfordTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                         Trajectory& traj) {
    (void)traj;
    const size_t N = tp.AtomCount();
    for (size_t i = 0; i < N; ++i) {
        TrajectoryAtom& ta = tp.MutableAtomAt(i);
        ta.bs_t0_std       = MomentsStd(ta.bs_t0_m2,       ta.bs_n_frames);
        ta.bs_t2mag_std    = MomentsStd(ta.bs_t2mag_m2,    ta.bs_n_frames);
        ta.bs_t0_delta_std = MomentsStd(ta.bs_t0_delta_m2, ta.bs_t0_delta_n);
    }

    finalized_ = true;

    OperationLog::Info(LogCalcOther, "BsWelfordTrajectoryResult::Finalize",
        "finalized across " + std::to_string(n_frames_) + " frames, " +
        std::to_string(N) + " atoms");
}


// ── WriteFeatures (NPY) ──────────────────────────────────────────
//
// bs_welford.npy: (N, 11) float64. Columns:
//   0: bs_t0_mean
//   1: bs_t0_std
//   2: bs_t0_min
//   3: bs_t0_max
//   4: bs_t2mag_mean
//   5: bs_t2mag_std
//   6: bs_t2mag_min
//   7: bs_t2mag_max
//   8: bs_t0_delta_mean
//   9: bs_t0_delta_std
//  10: n_frames (as double; downstream casts back to int)

int BsWelfordTrajectoryResult::WriteFeatures(
        const TrajectoryProtein& tp,
        const std::string& output_dir) const {
    const size_t N = tp.AtomCount();
    constexpr size_t K = 11;

    std::vector<double> data(N * K);
    for (size_t i = 0; i < N; ++i) {
        const auto& ta = tp.AtomAt(i);
        data[i * K + 0]  = ta.bs_t0_mean;
        data[i * K + 1]  = ta.bs_t0_std;
        data[i * K + 2]  = ta.bs_t0_min;
        data[i * K + 3]  = ta.bs_t0_max;
        data[i * K + 4]  = ta.bs_t2mag_mean;
        data[i * K + 5]  = ta.bs_t2mag_std;
        data[i * K + 6]  = ta.bs_t2mag_min;
        data[i * K + 7]  = ta.bs_t2mag_max;
        data[i * K + 8]  = ta.bs_t0_delta_mean;
        data[i * K + 9]  = ta.bs_t0_delta_std;
        data[i * K + 10] = static_cast<double>(ta.bs_n_frames);
    }

    NpyWriter::WriteFloat64(output_dir + "/bs_welford.npy",
                            data.data(), N, K);
    return 1;
}


// ── WriteH5Group ─────────────────────────────────────────────────
//
// /trajectory/bs_welford/ with a provenance header + one dataset per
// column. The provenance (name + frame count) is emitted as group
// attributes; downstream readers can identify which Result produced
// this group without consulting an external registry. This is the
// pattern all *TrajectoryResult::WriteH5Group emissions follow —
// rollups emit name + n_frames + any time-range attributes; time-series
// Results additionally emit full frame_indices/frame_times datasets
// alongside their dense buffers.
//
// Separate datasets per column (rather than one 2D block) keep the
// schema tolerant to per-column additions without affecting existing
// consumers.

void BsWelfordTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    const size_t N = tp.AtomCount();

    std::vector<double> t0_mean(N), t0_std(N), t0_min(N), t0_max(N);
    std::vector<double> t2mag_mean(N), t2mag_std(N), t2mag_min(N), t2mag_max(N);
    std::vector<double> t0_delta_mean(N), t0_delta_std(N);
    std::vector<size_t> n_frames(N);

    for (size_t i = 0; i < N; ++i) {
        const auto& ta = tp.AtomAt(i);
        t0_mean[i]       = ta.bs_t0_mean;
        t0_std[i]        = ta.bs_t0_std;
        t0_min[i]        = ta.bs_t0_min;
        t0_max[i]        = ta.bs_t0_max;
        t2mag_mean[i]    = ta.bs_t2mag_mean;
        t2mag_std[i]     = ta.bs_t2mag_std;
        t2mag_min[i]     = ta.bs_t2mag_min;
        t2mag_max[i]     = ta.bs_t2mag_max;
        t0_delta_mean[i] = ta.bs_t0_delta_mean;
        t0_delta_std[i]  = ta.bs_t0_delta_std;
        n_frames[i]      = ta.bs_n_frames;
    }

    auto grp = file.createGroup("/trajectory/bs_welford");

    // Provenance header.
    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_frames",    n_frames_);
    grp.createAttribute("finalized",   finalized_);

    grp.createDataSet("t0_mean",       t0_mean);
    grp.createDataSet("t0_std",        t0_std);
    grp.createDataSet("t0_min",        t0_min);
    grp.createDataSet("t0_max",        t0_max);
    grp.createDataSet("t2mag_mean",    t2mag_mean);
    grp.createDataSet("t2mag_std",     t2mag_std);
    grp.createDataSet("t2mag_min",     t2mag_min);
    grp.createDataSet("t2mag_max",     t2mag_max);
    grp.createDataSet("t0_delta_mean", t0_delta_mean);
    grp.createDataSet("t0_delta_std",  t0_delta_std);
    grp.createDataSet("n_frames_per_atom", n_frames);
}

}  // namespace nmr

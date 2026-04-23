#include "BsAnomalousAtomMarkerTrajectoryResult.h"
#include "BsWelfordTrajectoryResult.h"
#include "BiotSavartResult.h"
#include "TrajectoryProtein.h"
#include "TrajectoryAtom.h"
#include "ProteinConformation.h"
#include "ConformationAtom.h"
#include "AtomEvent.h"
#include "OperationLog.h"

#include <cmath>
#include <sstream>
#include <typeinfo>

namespace nmr {


std::vector<std::type_index>
BsAnomalousAtomMarkerTrajectoryResult::Dependencies() const {
    return {
        std::type_index(typeid(BsWelfordTrajectoryResult)),
        std::type_index(typeid(BiotSavartResult)),
    };
}


std::unique_ptr<BsAnomalousAtomMarkerTrajectoryResult>
BsAnomalousAtomMarkerTrajectoryResult::Create(const TrajectoryProtein& tp) {
    (void)tp;
    return std::make_unique<BsAnomalousAtomMarkerTrajectoryResult>();
}


// ── Compute ──────────────────────────────────────────────────────
//
// This is the canonical cross-Result read. For each atom:
//
//   READS:
//     conf.AtomAt(i).bs_shielding_contribution.T0      (this frame,
//                                                       via BiotSavartResult)
//     tp.AtomAt(i).bs_t0_mean                          (running,
//                                                       via BsWelfordTrajectoryResult)
//     tp.AtomAt(i).bs_t0_m2
//     tp.AtomAt(i).bs_n_frames
//
//   WRITES (on detected anomaly):
//     tp.MutableAtomAt(i).events.Push(...)             (this Result
//                                                       is the ONLY
//                                                       writer to this
//                                                       atom's event bag
//                                                       here — future
//                                                       emitters can push
//                                                       their own kinds)
//
// The three fields we read from tp.AtomAt(i) are owned (singleton-
// per-type + one-writer-per-field) by BsWelfordTrajectoryResult. We
// read them, we don't modify them — reading another Result's stashed
// state is what the "open buffer" discipline supports.

void BsAnomalousAtomMarkerTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)traj;
    const std::size_t N = tp.AtomCount();

    for (std::size_t i = 0; i < N; ++i) {
        const TrajectoryAtom& ta_const = tp.AtomAt(i);

        // CROSS-RESULT READ: bs_n_frames, bs_t0_m2, bs_t0_mean are
        // owned by BsWelfordTrajectoryResult. Attach order places
        // BsWelford first, so these reflect this frame's update.
        // See this class's header for the rationale; rename-changes
        // require BsWelford's writer-side marker to stay in sync.
        const std::size_t n = ta_const.bs_n_frames;
        if (n < MIN_BURN_IN_FRAMES) continue;

        const double m2 = ta_const.bs_t0_m2;
        const double std_est = std::sqrt(m2 / static_cast<double>(n - 1));
        if (std_est < 1e-10) continue;  // constant signal — no meaningful z

        const double current = conf.AtomAt(i).bs_shielding_contribution.T0;
        const double mean    = ta_const.bs_t0_mean;   // CROSS-RESULT READ (BsWelford)
        const double z       = (current - mean) / std_est;

        if (std::fabs(z) <= Z_THRESHOLD) continue;

        // Two kinds, one emitter: high vs low z-score. A downstream
        // consumer queries ta.events.ByKind<BsAnomalyHighT0>() or
        // .ByKind<BsAnomalyLowT0>() to pick one tail.
        const std::type_index kind =
            (z > 0.0)
                ? std::type_index(typeid(BsAnomalyHighT0))
                : std::type_index(typeid(BsAnomalyLowT0));

        std::ostringstream z_ss;
        z_ss << z;
        std::ostringstream mean_ss;
        mean_ss << mean;
        std::ostringstream curr_ss;
        curr_ss << current;

        TrajectoryAtom& ta_mut = tp.MutableAtomAt(i);
        ta_mut.events.Push(AtomEvent(
            std::type_index(typeid(BsAnomalousAtomMarkerTrajectoryResult)),
            kind,
            frame_idx,
            time_ps,
            {
                {"z_sigma",      z_ss.str()},
                {"current_t0",   curr_ss.str()},
                {"running_mean", mean_ss.str()},
                {"running_n",    std::to_string(n)},
            }));
        ++n_events_;
        if (z > 0.0) ++n_high_events_; else ++n_low_events_;
    }
}

}  // namespace nmr

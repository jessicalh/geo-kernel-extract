// BsWelfordTrajectoryResult.cpp
//
// Per WIP_OBJECT_MODEL.md §4 worked example. Bodies of Compute and
// Finalize are the spec's literal pseudocode.

#include "BsWelfordTrajectoryResult.h"

#include <cmath>

#include "TrajectoryAtom.h"
#include "TrajectoryProtein.h"

std::vector<std::type_index> BsWelfordTrajectoryResult::Dependencies() const {
    // Per-frame dependency: BiotSavartResult must run each frame.
    // typeid(BiotSavartResult) — BiotSavartResult is forward-declared
    // in Stubs.h.
    return { std::type_index(typeid(BiotSavartResult)) };
}

std::unique_ptr<BsWelfordTrajectoryResult>
BsWelfordTrajectoryResult::Create(const TrajectoryProtein& /*tp*/) {
    return std::make_unique<BsWelfordTrajectoryResult>();
}

void BsWelfordTrajectoryResult::Compute(
    const ProteinConformation& conf,
    TrajectoryProtein& tp,
    std::size_t /*frame_idx*/,
    double /*time_ps*/)
{
    // Read from frame's ConformationAtoms; write to tp's TrajectoryAtoms.
    // Welford update in place, field-by-field. Per §4 concrete Compute.
    const std::size_t N = tp.AtomCount();
    for (std::size_t i = 0; i < N; ++i) {
        const double x = conf.AtomAt(i).bs_shielding_contribution.T0;
        TrajectoryAtom& ta = tp.MutableAtomAt(i);

        // Running min/max.
        if (x < ta.bs_t0_min) ta.bs_t0_min = x;
        if (x > ta.bs_t0_max) ta.bs_t0_max = x;

        // Welford update. Fields are read, updated, written back.
        // After this frame, ta.bs_t0_mean is the mean over frames [0..idx].
        const std::size_t n = ta.bs_n_frames + 1;
        const double delta = x - ta.bs_t0_mean;
        const double new_mean = ta.bs_t0_mean + delta / static_cast<double>(n);
        const double new_m2 = ta.bs_t0_m2 + delta * (x - new_mean);
        ta.bs_t0_mean = new_mean;
        ta.bs_t0_m2 = new_m2;
        ta.bs_n_frames = n;
    }
    ++n_frames_;
}

void BsWelfordTrajectoryResult::Finalize(TrajectoryProtein& tp) {
    const std::size_t N = tp.AtomCount();
    for (std::size_t i = 0; i < N; ++i) {
        TrajectoryAtom& ta = tp.MutableAtomAt(i);
        if (ta.bs_n_frames > 1) {
            ta.bs_t0_std = std::sqrt(
                ta.bs_t0_m2 / static_cast<double>(ta.bs_n_frames - 1));
        } else {
            ta.bs_t0_std = 0.0;
        }
    }
    finalized_ = true;
}

double BsWelfordTrajectoryResult::MeanAtAtom(std::size_t atom_idx) const {
    // Intentionally do not take a TrajectoryProtein here — callers reach
    // tp.AtomAt(i).bs_t0_mean directly. Keeping the convenience accessor
    // as documented in §4 example would require either holding a back-
    // pointer to tp (not in the spec) or taking tp as an argument. The
    // spec shows `double MeanAtAtom(size_t atom_idx) const;` without a
    // tp parameter; we preserve the signature and return 0.0 here as
    // a sketch — a real implementation would pair this with a back-
    // pointer design the spec does not commit to. See README.
    (void)atom_idx;
    return 0.0;
}

double BsWelfordTrajectoryResult::StdAtAtom(std::size_t atom_idx) const {
    (void)atom_idx;
    return 0.0;  // same shape as MeanAtAtom; post-finalize field on tp.
}

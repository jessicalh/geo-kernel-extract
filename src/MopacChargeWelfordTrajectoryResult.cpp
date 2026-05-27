#include "MopacChargeWelfordTrajectoryResult.h"

#include "MopacResult.h"
#include "ConformationAtom.h"
#include "OperationLog.h"
#include "ProteinConformation.h"
#include "TrajectoryAtom.h"
#include "TrajectoryMoments.h"  // WelfordUpdate / WelfordFinalize
#include "TrajectoryProtein.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <cstddef>
#include <string>

namespace nmr {

std::unique_ptr<MopacChargeWelfordTrajectoryResult>
MopacChargeWelfordTrajectoryResult::Create(const TrajectoryProtein& tp) {
    (void)tp;
    return std::make_unique<MopacChargeWelfordTrajectoryResult>();
}


// ── Compute ──────────────────────────────────────────────────────
//
// Sparse-cadence gate: MopacResult attaches only when MOPAC ran this
// frame (TimedAttach in OperationRunner, not Require). HasResult
// false → skip Welford update + record mask=0. The per-atom
// `mopac_charge_welford.n_frames` counts only the frames the atom
// actually contributed a sample, which is the divisor WelfordFinalize
// needs for unbiased std.

void MopacChargeWelfordTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)traj;

    const bool source_present = conf.HasResult<MopacResult>();
    if (source_present) {
        const std::size_t N = conf.AtomCount();
        for (std::size_t i = 0; i < N; ++i) {
            auto& ws = tp.MutableAtomAt(i).mopac_charge_welford;
            const double q = conf.AtomAt(i).mopac_charge;
            const std::size_t n_new = ws.n_frames + 1;
            WelfordUpdate(ws.charge, q, n_new, frame_idx);
            ws.n_frames = n_new;
        }
        ++source_attached_count_;
    }
    // Sparse cadence is normal — don't log every absent frame.

    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    source_attached_per_frame_.push_back(source_present ? 1u : 0u);
    ++n_frames_;
}


// ── Finalize ─────────────────────────────────────────────────────
//
// Per canonical Welford TR pattern: call WelfordFinalize per atom to
// derive std from the accumulated frames. Only
// touch the per-atom Welford slots if at least one frame attached
// the source — same defensive pattern as AIMNet2ChargeResponseGradient
// Welford (MutableAtomAt accessing un-touched slots when no Compute
// ever updated them is technically safe with default-init but the
// gate keeps the test path symmetric to WriteH5Group's skip path).

void MopacChargeWelfordTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                  Trajectory& traj) {
    (void)traj;
    if (source_attached_count_ > 0) {
        const std::size_t N = tp.AtomCount();
        for (std::size_t i = 0; i < N; ++i) {
            auto& ws = tp.MutableAtomAt(i).mopac_charge_welford;
            WelfordFinalize(ws.charge, ws.n_frames);
        }
    }
    finalized_ = true;
    OperationLog::Info(LogCalcOther,
        "MopacChargeWelfordTrajectoryResult::Finalize",
        "finalized across " + std::to_string(n_frames_) +
        " frames; MOPAC ran on " +
        std::to_string(source_attached_count_) + " frames");
}


// ── WriteH5Group ─────────────────────────────────────────────────
//
// Skip the group entirely when MOPAC didn't run in ANY frame —
// canonical "absent, not faked" (OBJECT_MODEL "Conditional-attach TR
// discipline"). Otherwise emit the full canonical Welford row:
// mean / std / m2 / min / max / min_frame / max_frame + n_per_atom +
// frame metadata + source_attached_per_frame mask.

void MopacChargeWelfordTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp, HighFive::File& file) const {
    const std::size_t N = tp.AtomCount();
    const std::size_t T = n_frames_;

    if (source_attached_count_ == 0) {
        OperationLog::Warn(
            "MopacChargeWelfordTrajectoryResult::WriteH5Group",
            "MOPAC ran on 0/" + std::to_string(T) +
            " frames; skipping /trajectory/mopac_charge_welford/ "
            "per canonical 'absent, not faked' discipline.");
        return;
    }

    auto grp = file.createGroup("/trajectory/mopac_charge_welford");

    grp.createAttribute("result_name",            Name());
    grp.createAttribute("n_atoms",                N);
    grp.createAttribute("n_frames",               T);
    grp.createAttribute("source_attached_count",  source_attached_count_);
    grp.createAttribute("finalized",              finalized_);
    grp.createAttribute("units",                  std::string("e"));
    grp.createAttribute("source", std::string(
        "MopacResult.mopac_charge (Mulliken charge from PM7+MOZYME, "
        "units e). Per-atom Welford rollup; emits canonical row "
        "mean + std + m2 + min/max + min_frame/max_frame + n_per_atom. "
        "Minimum-viable v0 — no delta variants."));
    grp.createAttribute("source_attached_policy", std::string(
        "conditional -- MopacResult attaches sparsely per the Mopac "
        "cadence (OperationRunner.cpp:138-148, Attach gated by "
        "!opts.skip_mopac AND non-null Compute return; NOT "
        "RequireConformationResult). Compute's HasResult<MopacResult>() "
        "gate skips Welford update + records mask=0 on absent frames. "
        "WriteH5Group skips the entire group when source_attached_count==0."));

    // Per-dataset units attrs: e for mean/std/min/max; e^2 for m2;
    // frame_index for *_frame; frame_count for n_per_atom (codex H3
    // 2026-05-20 convention).
    std::vector<double>        mean(N), std_(N), m2(N), min_(N), max_(N);
    std::vector<std::uint64_t> min_frame(N), max_frame(N), n_per_atom(N);
    for (std::size_t i = 0; i < N; ++i) {
        const auto& w = tp.AtomAt(i).mopac_charge_welford.charge;
        mean[i]      = w.mean;
        std_[i]      = w.std;
        m2[i]        = w.m2;
        min_[i]      = w.min;
        max_[i]      = w.max;
        min_frame[i] = static_cast<std::uint64_t>(w.min_frame);
        max_frame[i] = static_cast<std::uint64_t>(w.max_frame);
        n_per_atom[i]
            = static_cast<std::uint64_t>(tp.AtomAt(i).mopac_charge_welford.n_frames);
    }

    auto with_units = [](HighFive::DataSet ds, const std::string& u) {
        ds.createAttribute("units", u);
    };
    with_units(grp.createDataSet("charge_mean",      mean),      std::string("e"));
    with_units(grp.createDataSet("charge_std",       std_),      std::string("e"));
    with_units(grp.createDataSet("charge_m2",        m2),        std::string("e^2"));
    with_units(grp.createDataSet("charge_min",       min_),      std::string("e"));
    with_units(grp.createDataSet("charge_max",       max_),      std::string("e"));
    with_units(grp.createDataSet("charge_min_frame", min_frame), std::string("frame_index"));
    with_units(grp.createDataSet("charge_max_frame", max_frame), std::string("frame_index"));
    with_units(grp.createDataSet("n_per_atom",       n_per_atom), std::string("frame_count"));

    grp.createDataSet("frame_indices",            frame_indices_);
    grp.createDataSet("frame_times",              frame_times_)
       .createAttribute("units", std::string("ps"));
    grp.createDataSet("source_attached_per_frame", source_attached_per_frame_);

    OperationLog::Info(LogCalcOther,
        "MopacChargeWelfordTrajectoryResult::WriteH5Group",
        "wrote /trajectory/mopac_charge_welford with " +
        std::to_string(N) + " atoms (" +
        std::to_string(source_attached_count_) + "/" + std::to_string(T) +
        " MOPAC-attached frames)");
}

}  // namespace nmr

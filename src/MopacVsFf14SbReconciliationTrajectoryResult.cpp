#include "MopacVsFf14SbReconciliationTrajectoryResult.h"

#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "CoulombResult.h"
#include "MopacCoulombResult.h"
#include "OperationLog.h"
#include "ProteinConformation.h"
#include "TrajectoryProtein.h"
#include "Types.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <cmath>
#include <limits>

namespace nmr {

namespace {

// |cos(T2_a, T2_b)| in the 5-vector T2 subspace. Returns NaN when
// either |T2| < threshold (cosine undefined).
double AbsCosT2(const SphericalTensor& a,
                const SphericalTensor& b,
                double threshold) {
    const double mag_a = a.T2Magnitude();
    const double mag_b = b.T2Magnitude();
    if (mag_a < threshold || mag_b < threshold) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    double dot = 0.0;
    for (std::size_t k = 0; k < 5; ++k) {
        dot += a.T2[k] * b.T2[k];
    }
    return std::abs(dot / (mag_a * mag_b));
}

}  // namespace


std::unique_ptr<MopacVsFf14SbReconciliationTrajectoryResult>
MopacVsFf14SbReconciliationTrajectoryResult::Create(
        const TrajectoryProtein& tp) {
    auto r = std::make_unique<
        MopacVsFf14SbReconciliationTrajectoryResult>();
    r->per_atom_abs_cos_.assign(tp.AtomCount(), std::vector<double>{});
    // Cache threshold once — same value MopacCoulombResult uses for
    // its near-zero E-field guard at line 244.
    r->zero_magnitude_threshold_ =
        CalculatorConfig::Get("near_zero_vector_norm_threshold");
    return r;
}


// ── Compute ──────────────────────────────────────────────────────
//
// Cross-source gate: BOTH MopacCoulombResult AND CoulombResult must
// be attached this frame for cosine to be defined. Either-absent →
// NaN-fill all atoms + source_attached_per_frame=0. When both attached,
// per-atom |cos|; per-atom NaN where either-side |T2| < threshold
// (undefined cosine).

void MopacVsFf14SbReconciliationTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;

    const bool both_present = conf.HasResult<MopacCoulombResult>()
                           && conf.HasResult<CoulombResult>();
    const double nan_d = std::numeric_limits<double>::quiet_NaN();

    const std::size_t N = conf.AtomCount();
    if (!both_present) {
        for (std::size_t i = 0; i < N; ++i) {
            per_atom_abs_cos_[i].push_back(nan_d);
        }
    } else {
        for (std::size_t i = 0; i < N; ++i) {
            const auto& ca = conf.AtomAt(i);
            per_atom_abs_cos_[i].push_back(
                AbsCosT2(ca.mopac_coulomb_shielding_contribution,
                         ca.coulomb_shielding_contribution,
                         zero_magnitude_threshold_));
        }
        ++source_attached_count_;
    }

    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    source_attached_per_frame_.push_back(both_present ? 1u : 0u);
    ++n_frames_;
}


// ── Finalize ─────────────────────────────────────────────────────

void MopacVsFf14SbReconciliationTrajectoryResult::Finalize(
        TrajectoryProtein& tp, Trajectory& traj) {
    (void)tp; (void)traj;
    finalized_ = true;
    OperationLog::Info(LogCalcOther,
        "MopacVsFf14SbReconciliationTrajectoryResult::Finalize",
        "finalized across " + std::to_string(n_frames_) +
        " frames; both MopacCoulomb + Coulomb attached on " +
        std::to_string(source_attached_count_) + " frames");
}


// ── WriteH5Group ─────────────────────────────────────────────────
//
// Skip the group entirely when no frame had both sources attached.
// Otherwise emit (N, T) double of |cos| values (with NaN where
// the cosine was undefined — see source attr).

void MopacVsFf14SbReconciliationTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    if (source_attached_count_ == 0) {
        OperationLog::Warn(
            "MopacVsFf14SbReconciliationTrajectoryResult::WriteH5Group",
            "both MopacCoulomb + Coulomb attached on 0/" +
            std::to_string(n_frames_) + " frames; skipping "
            "/trajectory/mopac_vs_ff14sb_reconciliation/ per canonical "
            "'absent, not faked' discipline.");
        return;
    }

    const std::size_t N = per_atom_abs_cos_.size();
    const std::size_t T = n_frames_;
    (void)tp;  // tp.AtomCount() may not match if any atoms were skipped;
               // we rely on the Compute-time push counts.

    auto grp = file.createGroup("/trajectory/mopac_vs_ff14sb_reconciliation");

    grp.createAttribute("result_name",            Name());
    grp.createAttribute("n_atoms",                N);
    grp.createAttribute("n_frames",               T);
    grp.createAttribute("source_attached_count",  source_attached_count_);
    grp.createAttribute("finalized",              finalized_);
    grp.createAttribute("parity",                 std::string("0e"));
    grp.createAttribute("units",                  std::string("dimensionless"));
    grp.createAttribute("sources", std::string(
        "MopacCoulombResult.mopac_coulomb_shielding_contribution + "
        "CoulombResult.coulomb_shielding_contribution. Per-atom-per-frame "
        "|cos| in the T2 5-vector subspace (both sources emit T2-only "
        "shielding contributions; absolute value because T2 has sign "
        "ambiguity from the eigenvector convention). |cos| ∈ [0, 1] "
        "where 1 = aligned tensor orientations, 0 = perpendicular."));
    grp.createAttribute("source_attached_policy", std::string(
        "conditional -- requires BOTH MopacCoulombResult (TimedAttach at "
        "OperationRunner.cpp:183) AND CoulombResult (TimedAttach at "
        "OperationRunner.cpp:178). Either-absent → all-N NaN cells + "
        "source_attached_per_frame=0. NO frame had both → WriteH5Group "
        "skips the group entirely per 'absent, not faked'."));
    grp.createAttribute("zero_magnitude_threshold", zero_magnitude_threshold_);
    grp.createAttribute("zero_magnitude_units",     std::string("ppm-like (matches "
        "Coulomb shielding contribution units; same threshold as MopacCoulombResult "
        "near_zero_vector_norm_threshold guard)"));

    // (N, T) flat. SDK readers use isfinite() to distinguish real
    // measurements from undefined-cosine cells.
    std::vector<double> flat(N * T);
    for (std::size_t i = 0; i < N; ++i) {
        const auto& atom_frames = per_atom_abs_cos_[i];
        const std::size_t T_atom = atom_frames.size();
        for (std::size_t t = 0; t < T; ++t) {
            flat[i * T + t] = (t < T_atom)
                ? atom_frames[t]
                : std::numeric_limits<double>::quiet_NaN();
        }
    }
    std::vector<std::size_t> dims = {N, T};
    HighFive::DataSpace space(dims);
    auto ds = grp.createDataSet<double>("abs_cos_t2", space);
    ds.write_raw(flat.data());

    grp.createDataSet("frame_indices",            frame_indices_);
    grp.createDataSet("frame_times",              frame_times_)
       .createAttribute("units", std::string("ps"));
    grp.createDataSet("source_attached_per_frame", source_attached_per_frame_);

    OperationLog::Info(LogCalcOther,
        "MopacVsFf14SbReconciliationTrajectoryResult::WriteH5Group",
        "wrote /trajectory/mopac_vs_ff14sb_reconciliation with " +
        std::to_string(N) + " atoms (" +
        std::to_string(source_attached_count_) + "/" + std::to_string(T) +
        " both-attached frames)");
}

}  // namespace nmr

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

std::unique_ptr<MopacVsFf14SbReconciliationTrajectoryResult>
MopacVsFf14SbReconciliationTrajectoryResult::Create(
        const TrajectoryProtein& tp) {
    auto r = std::make_unique<
        MopacVsFf14SbReconciliationTrajectoryResult>();
    r->per_atom_cos_.assign(tp.AtomCount(), std::vector<double>{});
    // Cache magnitude floor at Create. EFG-scale threshold (V/Å²) —
    // calibrated for the signal magnitude, NOT the project-wide
    // direction-vector floor (1e-10) which would admit FP-noise-
    // dominated atoms. Decision 2026-05-21 per math adversarial
    // review H1.
    r->magnitude_floor_ =
        CalculatorConfig::Get("coulomb_efg_t2_magnitude_floor");
    return r;
}


// ── Compute ──────────────────────────────────────────────────────
//
// Cross-source gate: BOTH MopacCoulombResult AND CoulombResult must
// be attached this frame for cosine to be defined. Either-absent →
// NaN-fill all atoms + source_attached_per_frame=0. When both attached,
// per-atom signed cosine via SphericalTensor::T2CosineWith;
// per-atom NaN where either-side |T2| < magnitude_floor (undefined
// cosine).

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
            per_atom_cos_[i].push_back(nan_d);
        }
    } else {
        for (std::size_t i = 0; i < N; ++i) {
            const auto& ca = conf.AtomAt(i);
            // Signed cos in [-1, 1] via the canonical SphericalTensor
            // method (returns NaN when either |T2| < magnitude_floor).
            per_atom_cos_[i].push_back(
                ca.mopac_coulomb_shielding_contribution.T2CosineWith(
                    ca.coulomb_shielding_contribution,
                    magnitude_floor_));
        }
        ++source_attached_count_;
    }

    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    source_attached_per_frame_.push_back(both_present ? 1u : 0u);
    ++n_frames_;
}


// ── Finalize ─────────────────────────────────────────────────────
//
// Idempotent: no state mutation beyond finalized_ flag and the info
// log (calling twice is harmless — second call sets the same flag,
// emits another log line). per_atom_cos_ buffers stay alive for
// WriteH5Group to read.

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
// Otherwise emit (N, T) double of signed cos values (with NaN where
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

    // N from the authoritative TrajectoryProtein axis. The internal
    // per_atom_cos_ buffer is sized to tp.AtomCount() at Create and
    // pushed unconditionally per frame per atom in Compute, so the
    // sizes match by structural invariant.
    const std::size_t N = tp.AtomCount();
    const std::size_t T = n_frames_;

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
        "CoulombResult.coulomb_shielding_contribution (both bare T2 "
        "EFG kernels in V/Å², no γ multiplication at extraction). "
        "Per-atom-per-frame signed cos in the T2 5-vector subspace via "
        "SphericalTensor::T2CosineWith (Frobenius inner product / "
        "magnitude product; isometric_real_sph normalization "
        "preserves the matrix Frobenius product). cos ∈ [-1, 1]: "
        "+1 = aligned tensor orientations, -1 = opposite-polarisation, "
        "0 = orthogonal. Calibration ridge MUST see the SIGNED cos "
        "(not |cos|) to expose chemistry-driven sign disagreement at "
        "qualitatively distinctive groups (decision 2026-05-21 per "
        "science adversarial review M1)."));
    grp.createAttribute("source_attached_policy", std::string(
        "conditional -- requires BOTH MopacCoulombResult (TimedAttach at "
        "OperationRunner.cpp:183) AND CoulombResult (TimedAttach at "
        "OperationRunner.cpp:178). Either-absent → all-N NaN cells + "
        "source_attached_per_frame=0. NO frame had both → WriteH5Group "
        "skips the group entirely per 'absent, not faked'."));
    grp.createAttribute("magnitude_floor",        magnitude_floor_);
    grp.createAttribute("magnitude_floor_units",  std::string("V/Å^2"));
    grp.createAttribute("magnitude_floor_source", std::string(
        "CalculatorConfig::Get(\"coulomb_efg_t2_magnitude_floor\") — "
        "EFG-scale floor on |T2| for cosine well-definedness. Not the "
        "project-wide near_zero_vector_norm_threshold (1e-10) which is "
        "calibrated for direction-vector normalization, NOT EFG "
        "magnitudes; using that would let FP-noise-dominated atoms "
        "(remote-from-charge with |T2| ~ 1e-8 V/Å²) leak into the "
        "calibration distribution as spurious |cos|≈1 tails. "
        "Decision 2026-05-21 per math adversarial review H1."));

    // (N, T) flat. SDK readers use isfinite() to distinguish real
    // signed-cosine measurements from undefined-cosine cells.
    std::vector<double> flat(N * T);
    for (std::size_t i = 0; i < N; ++i) {
        const auto& atom_frames = per_atom_cos_[i];
        for (std::size_t t = 0; t < T; ++t) {
            flat[i * T + t] = atom_frames[t];
        }
    }
    std::vector<std::size_t> dims = {N, T};
    HighFive::DataSpace space(dims);
    auto ds = grp.createDataSet<double>("cos_t2", space);
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

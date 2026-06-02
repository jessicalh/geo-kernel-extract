#include "Dssp8TransitionTrajectoryResult.h"
#include "Dssp8TimeSeriesTrajectoryResult.h"  // Ss8Code helper

#include "DsspResult.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "TrajectoryProtein.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <cstdint>
#include <string>
#include <typeinfo>

namespace nmr {

namespace {

template <typename Arr>
std::uint8_t Argmax8OrUnassigned(const Arr& arr) {
    std::uint32_t best = 0;
    std::uint8_t  idx  = Dssp8TransitionTrajectoryResult::kSSUnassigned;
    for (std::size_t i = 0; i < arr.size(); ++i) {
        if (arr[i] > best) {
            best = arr[i];
            idx  = static_cast<std::uint8_t>(i);
        }
    }
    return idx;
}

}  // anonymous namespace


std::unique_ptr<Dssp8TransitionTrajectoryResult>
Dssp8TransitionTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<Dssp8TransitionTrajectoryResult>();
    const std::size_t R = tp.ProteinRef().ResidueCount();
    r->prev_ss_.assign(R, kSSUnassigned);
    r->ss8_transition_count_.assign(R, 0);
    r->ss8_occupancy_.assign(R, {});
    r->ss8_transition_matrix_.assign(R, {});
    r->n_frames_observed_.assign(R, 0);
    r->ss8_dominant_.assign(R, kSSUnassigned);
    return r;
}


void Dssp8TransitionTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;

    const bool source_attached = force_source_present_for_testing_
        || conf.HasResult<DsspResult>();
    source_attached_per_frame_.push_back(source_attached ? 1u : 0u);

    const std::size_t R = conf.ProteinRef().ResidueCount();

    if (source_attached) {
        const auto& dssp_residues = conf.Result<DsspResult>().AllResidues();
        const std::size_t dssp_R = dssp_residues.size();
        for (std::size_t ri = 0; ri < R; ++ri) {
            std::uint8_t cur = kSSUnassigned;
            // Codex review 2026-05-19: only count a residue as observed
            // when the DSSP row actually mapped (observed=true),
            // NOT just when ri < dssp_R. Otherwise unmapped residues
            // (caps, insertion-code mismatches) silently look like coil
            // and inflate transition counts to/from C.
            if (ri < dssp_R && dssp_residues[ri].observed) {
                cur = Dssp8TimeSeriesTrajectoryResult::Ss8Code(
                    dssp_residues[ri].secondary_structure);
            }
            if (cur != kSSUnassigned) {
                ++ss8_occupancy_[ri][cur];
                ++n_frames_observed_[ri];
                const std::uint8_t prev = prev_ss_[ri];
                if (prev != kSSUnassigned && prev != cur) {
                    ++ss8_transition_count_[ri];
                    ++ss8_transition_matrix_[ri][prev][cur];
                }
            }
            prev_ss_[ri] = cur;
        }
    } else {
        // Source absent — break the transition chain. All residues'
        // prev_ss becomes unassigned so the next consecutive-observed
        // pair starts fresh.
        for (std::size_t ri = 0; ri < R; ++ri) {
            prev_ss_[ri] = kSSUnassigned;
        }
    }

    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    ++n_frames_;
}


void Dssp8TransitionTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                Trajectory& traj) {
    (void)tp; (void)traj;
    const std::size_t R = ss8_occupancy_.size();
    for (std::size_t ri = 0; ri < R; ++ri) {
        ss8_dominant_[ri] = Argmax8OrUnassigned(ss8_occupancy_[ri]);
    }
    finalized_ = true;

    OperationLog::Info(LogCalcOther,
        "Dssp8TransitionTrajectoryResult::Finalize",
        "finalized across " + std::to_string(n_frames_) +
        " frames, " + std::to_string(R) + " residues.");
}


void Dssp8TransitionTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    (void)tp;
    const std::size_t R = ss8_occupancy_.size();
    const std::size_t T = n_frames_;

    std::size_t source_attached_count = 0;
    for (auto v : source_attached_per_frame_) if (v) ++source_attached_count;
    if (source_attached_count == 0) {
        OperationLog::Warn(
            "Dssp8TransitionTrajectoryResult::WriteH5Group",
            "DsspResult was not attached in any of " +
            std::to_string(source_attached_per_frame_.size()) +
            " frames; skipping /trajectory/dssp8_transition/ emission.");
        return;
    }

    auto grp = file.createGroup("/trajectory/dssp8_transition");

    grp.createAttribute("result_name",  Name());
    grp.createAttribute("n_residues",   R);
    grp.createAttribute("n_frames",     T);
    grp.createAttribute("finalized",    finalized_);
    grp.createAttribute("source_attached_count", source_attached_count);
    grp.createAttribute("ss8_count",    kSSCount);

    grp.createAttribute("ss8_legend", std::string(
        "H=0 (alpha helix), G=1 (3_10 helix), I=2 (pi helix), "
        "E=3 (extended strand), B=4 (beta bridge), T=5 (turn), "
        "S=6 (bend), C=7 (coil); 255 = no observation (dominant only)."));
    grp.createAttribute("transition_matrix_layout", std::string(
        "M[ri, prev, curr] = count of consecutive observed-pair "
        "transitions (prev -> curr) at residue ri. Self-edges (prev == "
        "curr) are NOT counted; the diagonal is identically zero. "
        "Sum over the (prev, curr) off-diagonal equals "
        "ss8_transition_count[ri]."));
    grp.createAttribute("cadence_caveat", std::string(
        "Transition counts are stride-sensitive. At <= 100 ps cadence, "
        "fast SS fluctuations are captured; slow secondary-structure "
        "rearrangements (helix unfolding events ~ns-us) require finer "
        "or longer-trajectory cadence to resolve. Cadence is set by "
        "--stride (default 1 = every TRR frame); coarser strides miss "
        "progressively faster SS dynamics (e.g. stride=300 / 3 ns/frame "
        "misses most transitions)."));
    grp.createAttribute("transition_gate", std::string(
        "Both prev and curr frame must have an observed SS code "
        "(DsspResult attached AND code != unassigned). Source-absent "
        "frames break the chain: prev resets to unassigned so the "
        "next observed pair starts a fresh transition count."));
    grp.createAttribute("source", std::string(
        "DsspResult (libdssp via Joosten 2011 / Kabsch-Sander 1983)."));
    grp.createAttribute("source_attached_policy", std::string(
        "conditional -- DsspResult attaches only when "
        "PerFrameRunOptions::skip_dssp == false. Reader contract: "
        "consult source_attached_per_frame; n_frames_observed[ri] "
        "captures the per-residue effective sample size."));
    grp.createAttribute("residue_axis", std::string("protein_residue_index"));

    // ── Per-frame (T,) ───────────────────────────────────────────────
    grp.createDataSet("frame_indices", frame_indices_)
       .createAttribute("units", std::string("frame_index"));
    grp.createDataSet("frame_times",   frame_times_)
       .createAttribute("units", std::string("ps"));
    grp.createDataSet("source_attached_per_frame", source_attached_per_frame_)
       .createAttribute("units", std::string("dimensionless"));

    // ── Per-residue stats (R,) ───────────────────────────────────────
    grp.createDataSet("ss8_transition_count", ss8_transition_count_)
       .createAttribute("units", std::string("transitions"));
    grp.createDataSet("ss8_dominant", ss8_dominant_)
       .createAttribute("units", std::string("category"));
    grp.createDataSet("n_frames_observed", n_frames_observed_)
       .createAttribute("units", std::string("frame_count"));

    // ── ss8_occupancy (R, 8) uint32 ──────────────────────────────────
    {
        std::vector<std::uint32_t> flat(R * kSSCount);
        for (std::size_t ri = 0; ri < R; ++ri) {
            for (std::size_t s = 0; s < kSSCount; ++s) {
                flat[ri * kSSCount + s] = ss8_occupancy_[ri][s];
            }
        }
        const std::vector<std::size_t> dims = {R, kSSCount};
        HighFive::DataSpace space(dims);
        auto ds = grp.createDataSet<std::uint32_t>("ss8_occupancy", space);
        ds.write_raw(flat.data());
        ds.createAttribute("units", std::string("frame_count"));
    }

    // ── ss8_transition_matrix (R, 8, 8) uint32 ───────────────────────
    {
        std::vector<std::uint32_t> flat(R * kSSCount * kSSCount);
        for (std::size_t ri = 0; ri < R; ++ri) {
            for (std::size_t prev = 0; prev < kSSCount; ++prev) {
                for (std::size_t curr = 0; curr < kSSCount; ++curr) {
                    flat[(ri * kSSCount + prev) * kSSCount + curr] =
                        ss8_transition_matrix_[ri][prev][curr];
                }
            }
        }
        const std::vector<std::size_t> dims = {R, kSSCount, kSSCount};
        HighFive::DataSpace space(dims);
        auto ds = grp.createDataSet<std::uint32_t>(
            "ss8_transition_matrix", space);
        ds.write_raw(flat.data());
        ds.createAttribute("units", std::string("transitions"));
        ds.createAttribute("axis_2_label", std::string("prev_ss"));
        ds.createAttribute("axis_3_label", std::string("curr_ss"));
    }

    OperationLog::Info(LogCalcOther,
        "Dssp8TransitionTrajectoryResult::WriteH5Group",
        "wrote /trajectory/dssp8_transition with " +
        std::to_string(R) + " residues x " + std::to_string(T) +
        " frames (source attached in " +
        std::to_string(source_attached_count) + "/" +
        std::to_string(source_attached_per_frame_.size()) + " frames)");
}


}  // namespace nmr

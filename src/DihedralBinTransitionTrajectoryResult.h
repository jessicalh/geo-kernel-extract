#pragma once
//
// DihedralBinTransitionTrajectoryResult: per-residue rotamer + Rama
// region transition counters and bin occupancy, accumulated across all
// frames of a run.
//
// AV companion to `DihedralTimeSeriesTrajectoryResult` (FO per-frame
// timeline). DihedralTS emits raw phi/psi/omega/chi(R,T) — this TR
// summarises the bin-classified dynamics: how many transitions per
// residue, what fraction of frames each bin was occupied, which bin
// dominated, how many frames contributed.
//
// Why both? DihedralTS is the per-frame source of truth (movie
// playback + downstream re-bin against any Rama/rotamer convention).
// BinTransition is the cheap-to-consume summary (transition rates per
// nanosecond, dominant-state classification, dynamics quality).
//
// Independent compute: this TR queries positions + AminoAcidType.
// chi_angles + Protein::BackbonePredecessor / BackboneSuccessor each
// frame. Does NOT cross-read DihedralTS (which would couple attach
// order and require DihedralTS to expose mid-stream per-frame fields).
// Same Dihedral/RamachandranBin/BinForChi logic as DihedralTS +
// ChiRotamerSelection; copies kept small per PATTERNS Lesson 10.
//
// Pure AV pattern: Compute updates running counters in place, no
// DenseBuffer, no per-frame buffer. Per-residue state lives internally
// to the TR (per-residue prev-bin state + per-residue accumulators),
// not on a TrajectoryResidue (no such type exists; per PATTERNS §13
// trajectory atom is the per-atom store; residue-scope state lives on
// the owning Result, same as BondLengthStats's per-bond state).
//
// Emission: /trajectory/dihedral_bin_transition/
//   Per-residue stats (R,):
//     backbone_transition_count   uint32
//     backbone_dominant_region    uint8
//     n_frames_observed           uint32  (phi+psi both finite)
//
//   Per-residue per backbone bin (R, 6):
//     backbone_bin_occupancy      uint32  frame count per region
//                                 (bin indexing matches DihedralTS:
//                                  0=unassigned, 1=αR, 2=β, 3=αL,
//                                  4=PPII, 5=other)
//
//   Per-residue per chi (R, 4):
//     chi_transition_count        uint32
//     chi_dominant_rotamer        uint8   (0=g+, 1=t, 2=g-, 255=N/A)
//     chi_n_frames_observed       uint32  (chi[k] finite)
//
//   Per-residue per chi per rotamer bin (R, 4, 3):
//     chi_rotamer_occupancy       uint32  frame count per (chi, bin)
//                                 (bin 0=g+, 1=t, 2=g-)
//
//   Per-frame metadata:
//     frame_indices, frame_times, source_attached_per_frame
//
// Transition gate: BOTH prev and curr frame must be observed (both
// finite phi/psi for the backbone bin; both finite chi[k] for chi
// bin). Wrap-tolerant via bin labels, not raw angle deltas.
//
// Conventions match DihedralTS verbatim (legend, boundaries,
// resolution order, chi rotamer bin endpoints). Documented as group
// attrs so consumers see the same vocabulary across both groups.
//
// Source-attached gate: positions always present at tp.Seed time;
// source_attached_per_frame emitted as all-1 for SDK uniformity
// (OBJECT_MODEL.md "Conditional-attach TR discipline").
//

#include "TrajectoryResult.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class DihedralBinTransitionTrajectoryResult : public TrajectoryResult {
public:
    // 6-region literature Rama grid (matches DihedralTS exactly):
    //   0 = unassigned (phi or psi NaN), 1 = αR, 2 = β, 3 = αL,
    //   4 = PPII, 5 = other.
    static constexpr std::size_t kBackboneBinCount = 6;
    static constexpr std::size_t kChiCount = 4;
    static constexpr std::size_t kChiBinCount = 3;  // g+, t, g-

    static constexpr std::uint8_t kBinUnassigned   = 0;
    static constexpr std::uint8_t kBinAlphaR       = 1;
    static constexpr std::uint8_t kBinBeta         = 2;
    static constexpr std::uint8_t kBinAlphaL       = 3;
    static constexpr std::uint8_t kBinPPII         = 4;
    static constexpr std::uint8_t kBinOther        = 5;

    static constexpr std::uint8_t kChiBinGplus     = 0;
    static constexpr std::uint8_t kChiBinTrans     = 1;
    static constexpr std::uint8_t kChiBinGminus    = 2;
    static constexpr std::uint8_t kChiBinUnassigned = 255;

    std::string Name() const override {
        return "DihedralBinTransitionTrajectoryResult";
    }

    // No ConformationResult dependency — reads positions + AminoAcidType
    // chi_angles (via Residue.chi[k]) + Protein::BackbonePredecessor /
    // BackboneSuccessor each frame.
    std::vector<std::type_index> Dependencies() const override {
        return {};
    }

    static std::unique_ptr<DihedralBinTransitionTrajectoryResult>
    Create(const TrajectoryProtein& tp);

    void Compute(const ProteinConformation& conf,
                 TrajectoryProtein& tp,
                 Trajectory& traj,
                 std::size_t frame_idx,
                 double time_ps) override;

    void Finalize(TrajectoryProtein& tp, Trajectory& traj) override;

    void WriteH5Group(const TrajectoryProtein& tp,
                      HighFive::File& file) const override;

    std::size_t NumFrames() const { return n_frames_; }

private:
    // Per-residue previous-frame bin state. kBinUnassigned == 0 marks
    // "no observed previous frame" (either this is frame 0 for the
    // residue, or last frame's bin was unassigned). This is the
    // consecutive-frame gate for transition counting.
    std::vector<std::uint8_t> prev_backbone_bin_;
    std::vector<std::array<std::uint8_t, kChiCount>> prev_chi_bin_;

    // Per-residue running accumulators (R-sized).
    std::vector<std::uint32_t> backbone_transition_count_;
    std::vector<std::array<std::uint32_t, kBackboneBinCount>>
        backbone_bin_occupancy_;
    std::vector<std::uint32_t> n_frames_observed_;

    std::vector<std::array<std::uint32_t, kChiCount>>
        chi_transition_count_;
    std::vector<std::array<std::uint32_t, kChiCount>>
        chi_n_frames_observed_;
    // (R, 4, 3): per residue, per chi, per bin.
    std::vector<
        std::array<std::array<std::uint32_t, kChiBinCount>, kChiCount>>
        chi_rotamer_occupancy_;

    // Finalize-derived per-residue dominants.
    std::vector<std::uint8_t> backbone_dominant_region_;
    std::vector<std::array<std::uint8_t, kChiCount>>
        chi_dominant_rotamer_;

    // Per-frame metadata.
    std::vector<std::size_t>  frame_indices_;
    std::vector<double>       frame_times_;
    std::vector<std::uint8_t> source_attached_per_frame_;

    std::size_t n_frames_  = 0;
    bool        finalized_ = false;
};

}  // namespace nmr

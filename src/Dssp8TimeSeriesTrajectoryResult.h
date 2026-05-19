#pragma once
//
// Dssp8TimeSeriesTrajectoryResult: per-residue per-frame DSSP secondary
// structure + DSSP H-bond partners and energies.
//
// Reads `DsspResult` per frame; emits the 8-state SS code + the four
// H-bond slots (2 acceptors + 2 donors per residue, each with a
// partner residue index and a kcal/mol energy). Phi/psi/SASA from
// DsspResult are NOT mirrored here -- they live in their own canonical
// homes (DihedralTimeSeriesTrajectoryResult owns phi/psi; SasaResult
// owns SASA).
//
// Conditional-source TR: DsspResult attaches only when
// PerFrameRunOptions::skip_dssp is false. Per the 2026-05-15
// conditional-attach discipline (OBJECT_MODEL.md), this TR uses the
// source-attached gate -- when no frame had DsspResult attached, the
// H5 group is skipped entirely. Float data NaN-filled on absent frames;
// int data (ss8_code, partner indices) uses out-of-range sentinels
// (255 for ss8_code, -1 for partner index).
//
// Emission: /trajectory/dssp8_time_series/
//   Per-residue per-frame (R, T):
//     ss8_code             uint8   H=0, G=1, I=2, E=3, B=4, T=5,
//                                  S=6, C=7; 255 = no observation
//     hbond_acceptor_partner (R, T, 2) int32   partner residue index
//                                              or -1 if no partner
//     hbond_acceptor_energy  (R, T, 2) float64 kcal/mol; NaN if no partner
//     hbond_donor_partner    (R, T, 2) int32
//     hbond_donor_energy     (R, T, 2) float64 kcal/mol
//
//   Per-atom lookup:
//     residue_index_per_atom (N,) int32 -- atom_i -> residue_i broadcast
//                                          for the viewer / SDK atom-axis
//                                          consumer.
//
//   Per-frame metadata:
//     frame_indices, frame_times, source_attached_per_frame
//
// Convention attrs:
//   ss8_legend  = "H,G,I,E,B,T,S,C" (codes 0..7; 255 = no observation)
//   ss8_meaning = "alpha helix, 3_10 helix, pi helix, extended strand,
//                   beta bridge, turn, bend, coil"
//   hbond_partner_sentinel = "-1 means no partner"
//   hbond_energy_units = "kcal/mol"
//   hbond_energy_absent_sentinel = "NaN"
//   source = "DsspResult (libdssp via Joosten 2011 / Kabsch-Sander 1983)"
//   source_attached_policy = "conditional" -- attaches only when
//     PerFrameRunOptions::skip_dssp is false; consult
//     source_attached_per_frame mask before downstream stats.
//
// AV-pattern caveat: this TR is FO (per-frame buffers accumulated
// during Compute, flattened at WriteH5Group). DSSP transitions /
// occupancy stats live in the AV companion
// `Dssp8TransitionTrajectoryResult`.
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

class Dssp8TimeSeriesTrajectoryResult : public TrajectoryResult {
public:
    static constexpr std::uint8_t kSSUnassigned = 255;
    static constexpr std::int32_t kNoPartner = -1;

    // SS char -> uint8 code. H/G/I/E/B/T/S/C → 0..7. Anything else
    // (DSSP returns ' ' for unassigned, mapped to 'C' upstream by
    // DsspResult::Compute; 'C' itself for coil) maps to 7.
    static std::uint8_t Ss8Code(char ss) {
        switch (ss) {
            case 'H': return 0;
            case 'G': return 1;
            case 'I': return 2;
            case 'E': return 3;
            case 'B': return 4;
            case 'T': return 5;
            case 'S': return 6;
            default:  return 7;
        }
    }

    std::string Name() const override {
        return "Dssp8TimeSeriesTrajectoryResult";
    }

    // No declared dependency: DsspResult attaches conditionally (when
    // skip_dssp is false). We read whatever's there per frame and
    // record the source-attached mask. Declaring it would force
    // Phase-4 attach validation to fail on legitimate no-DSSP fleet
    // runs.
    std::vector<std::type_index> Dependencies() const override {
        return {};
    }

    static std::unique_ptr<Dssp8TimeSeriesTrajectoryResult>
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

    // Test-only: bypass the per-frame `conf.HasResult<DsspResult>()`
    // check. Used by synthetic unit tests that bypass OperationRunner.
    //
    // SAFETY (codex review 2026-05-19): When the flag is true, Compute
    // takes the source-attached branch and calls `conf.Result<DsspResult>()`
    // — which THROWS if DsspResult was not actually attached. So the
    // test caller MUST attach DsspResult to every ProteinConformation
    // it passes to Compute before setting this flag. (TripeptideBB-
    // and Larsen-family TRs read a per-atom field instead of
    // .Result<>() and are safe; THIS TR reads the Result object, so
    // the constraint is real.)
    void ForceSourcePresentForTesting() {
        force_source_present_for_testing_ = true;
    }

private:
    // Per-residue per-frame growing buffers. R outer x T inner.
    std::vector<std::vector<std::uint8_t>>  ss8_code_;
    std::vector<std::vector<std::array<std::int32_t, 2>>> hbond_acceptor_partner_;
    std::vector<std::vector<std::array<double, 2>>>       hbond_acceptor_energy_;
    std::vector<std::vector<std::array<std::int32_t, 2>>> hbond_donor_partner_;
    std::vector<std::vector<std::array<double, 2>>>       hbond_donor_energy_;

    // Static atom -> residue lookup for the SDK / viewer broadcast.
    std::vector<std::int32_t> residue_index_per_atom_;

    // Per-frame metadata.
    std::vector<std::size_t>  frame_indices_;
    std::vector<double>       frame_times_;
    std::vector<std::uint8_t> source_attached_per_frame_;
    bool force_source_present_for_testing_ = false;

    std::size_t n_frames_  = 0;
    bool        finalized_ = false;
};

}  // namespace nmr

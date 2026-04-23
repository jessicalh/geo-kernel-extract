#pragma once
//
// BsT0AutocorrelationTrajectoryResult: per-atom normalized
// autocorrelation function of the BiotSavart T0 shielding contribution
// at lags k = 0, 1, ..., N_LAGS − 1. Worked example of:
//
//   1. Per-atom buffered history during Compute. The Result keeps a
//      `std::vector<std::vector<double>>` — one growing buffer per
//      atom, Finalize then operates on the full series. All
//      accumulator state lives INSIDE the Result (never on
//      TrajectoryAtom).
//
//   2. Finalize-heavy computation. Compute appends one scalar per
//      atom per frame; Finalize walks each atom's history and
//      produces the classic biased autocorrelation
//
//          ρ(k) = C(k) / C(0),
//          C(k) = (1/N) · Σ_{t=0..N-k-1} (x_t − ⟨x⟩)(x_{t+k} − ⟨x⟩)
//
//      This estimator guarantees |ρ(k)| ≤ 1 at all finite N, unlike
//      the naive ratio-of-running-sums variant which drifts outside
//      [−1, +1] when the right-tail mean differs from the full-range
//      mean. Canonical shape for per-atom-per-lag TrajectoryResults
//      (spectral density, cross-correlation diagnostics, relaxation
//      times) — clone this file and swap the scalar source.
//
//   3. Per-atom × lag DenseBuffer<double>. Third concrete T in the
//      worked-example set (Vec3 via Positions, SphericalTensor via
//      BsShielding, double here). Emission uses explicit indexed
//      access — no reinterpret cast, same discipline as Positions.
//
// Memory note: the full history buffer is O(atoms × n_frames × 8 B).
// At the 685-protein fleet scale with stride 2 on 25 ns runs
// (≈ 625 frames × ~4000 atoms) that's ~20 MB per run — comfortable.
// If a later session needs bounded memory for μs-scale runs, the
// refactor is: replace the per-atom history vector with a circular
// window + per-lag running left/right tail sums, then Finalize
// becomes a closed-form formula over running state. The current
// exemplar prioritizes correctness + clarity; the bounded-memory
// variant is strictly a swap-in optimization.
//
// Emission: /trajectory/bs_t0_autocorrelation/
//   rho             (N, N_LAGS) float64   — autocorrelation ρ(k)
//   lag_frames      (N_LAGS,)   uint64    — k = 0..N_LAGS-1
//   lag_times_ps    (N_LAGS,)   float64   — k × median-frame-interval
//   attrs: result_name, n_atoms, n_lags, n_frames, finalized,
//          sample_interval_ps, units="dimensionless"
//

#include "TrajectoryResult.h"

#include <cstddef>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class BsT0AutocorrelationTrajectoryResult : public TrajectoryResult {
public:
    // Number of lags computed. 120 samples at stride-2 / 20-ps frame
    // cadence ≈ 4.8 ns of lag — enough to see the fall-off of typical
    // shielding autocorrelations. Set at compile time for simplicity;
    // future Results can template over N_LAGS.
    static constexpr std::size_t N_LAGS = 120;

    std::string Name() const override {
        return "BsT0AutocorrelationTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<BsT0AutocorrelationTrajectoryResult> Create(
        const TrajectoryProtein& tp);

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
    // Per-atom growing buffer of bs_t0 values, one entry per frame.
    // At Finalize this is walked to compute ρ(k); the vectors are
    // swapped out and released as the dense buffer is populated.
    std::vector<std::vector<double>> per_atom_history_;

    // Per-Compute recorded frame times. The Result records its own
    // provenance rather than reading from traj.FrameTimes() — that
    // lets manually-orchestrated callers (tests, diagnostics) drive
    // the Result without threading frame metadata through Trajectory.
    std::vector<double> frame_times_;

    std::size_t n_frames_ = 0;
    double sample_interval_ps_ = 0.0;  // median Δt from frame_times_ at Finalize
    bool finalized_ = false;
};

}  // namespace nmr

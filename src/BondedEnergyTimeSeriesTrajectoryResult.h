#pragma once
//
// BondedEnergyTimeSeriesTrajectoryResult: per-atom per-frame timeline of
// the GROMACS CHARMM36m bonded-energy decomposition computed by
// BondedEnergyResult. Seven channels per atom per frame: bond, angle,
// Urey-Bradley, proper dihedral, improper dihedral, CMAP, total.
//
// Shape: per-atom, multi-channel. Internal storage is seven parallel
// std::vector<std::vector<double>> indexed [atom][frame]; H5 emits one
// (N, T) dataset per channel. No DenseBuffer adoption — bonded-energy
// channels are a leaf consumer (no cross-TR read of per-frame breakdown).
//
// Use case: ML feature input for per-atom strain-energy environment
// across the trajectory; pairs with GromacsEnergyTimeSeries (system-
// scalar) which gives the global thermodynamic state for the same frame.
// Both together unlock low-energy-state filtering AND per-atom strain
// features for downstream ridge / e3nn.
//
// Export-everything-upstream (PATTERNS Lesson 25): all 7 channels emitted
// even though `total` is the running sum of the other 6. Calibration may
// weight the decomposition separately from the total.
//
// Emission:
//
//   /trajectory/bonded_energy_time_series/
//     bond           (N, T)  float64  kJ/mol
//     angle          (N, T)  float64  kJ/mol
//     urey_bradley   (N, T)  float64  kJ/mol
//     proper_dih     (N, T)  float64  kJ/mol
//     improper_dih   (N, T)  float64  kJ/mol
//     cmap_dih       (N, T)  float64  kJ/mol     (matches GromacsEnergy.cmap_dih)
//     total          (N, T)  float64  kJ/mol
//     frame_indices  (T,)    uint64
//     frame_times    (T,)    float64  ps
//     attrs: result_name, n_atoms, n_frames, finalized,
//            units = "kJ/mol",
//            split_convention = "evenly_among_participating_atoms"
//
// The split_convention attribute pins how the source BondedEnergyResult
// distributes interaction energy across the 2..5 atoms each interaction
// touches: even share = energy / n_atoms_in_interaction. Downstream that
// wants whole-system totals can sum bond[:, t] across atoms.
//
// Dependencies: BondedEnergyResult (conformation-scope, unconditionally
// attached in PerFrameExtractionSet — no source-attached gate needed).
//

#include "TrajectoryResult.h"

#include <cstddef>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class BondedEnergyTimeSeriesTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "BondedEnergyTimeSeriesTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<BondedEnergyTimeSeriesTrajectoryResult> Create(
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
    // Per-atom growing buffers, one per channel. After Finalize each
    // bond_[i] / angle_[i] / ... has size n_frames_.
    std::vector<std::vector<double>> bond_;
    std::vector<std::vector<double>> angle_;
    std::vector<std::vector<double>> urey_bradley_;
    std::vector<std::vector<double>> proper_dih_;
    std::vector<std::vector<double>> improper_dih_;
    std::vector<std::vector<double>> cmap_;
    std::vector<std::vector<double>> total_;
    std::vector<std::size_t>         frame_indices_;
    std::vector<double>              frame_times_;
    std::size_t                      n_frames_ = 0;
    bool                             finalized_ = false;
};

}  // namespace nmr

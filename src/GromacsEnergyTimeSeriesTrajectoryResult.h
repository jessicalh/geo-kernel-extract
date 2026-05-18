#pragma once
//
// GromacsEnergyTimeSeriesTrajectoryResult: per-frame timeline of all
// system-scalar GROMACS energy / thermodynamic / box / virial / pressure
// terms read from the .edr energy file. The trajectory-scope rollup of
// the per-frame GromacsEnergyResult; one row per frame, one column per
// channel from the GromacsEnergy struct.
//
// Shape: PER-FRAME SYSTEM SCALARS — not per-atom. This TR breaks the
// per-atom DenseBuffer pattern that the BS / HM / Sasa / APBS time-series
// TRs follow. The data being rolled up is whole-system thermodynamic
// state, not per-atom shielding. Storage is a flat std::vector<GromacsEnergy>
// indexed by frame.
//
// Use case: this is the gate for "low-energy-state" ML model framings.
// Once emitted, downstream selection of the bottom-N% lowest-energy frames
// for training is a one-liner:
//
//   energy = traj.energy.gromacs.total_energy[:]
//   low_idx = np.argsort(energy)[:int(0.1 * len(energy))]
//
// Export-everything-upstream (PATTERNS Lesson 25): emit every column of
// the GromacsEnergy struct as its own dataset, including the 3×3 virial
// and pressure tensors. Storage at thesis scale is negligible; removed
// columns are forever, added ones are free.
//
// Emission:
//
//   /trajectory/gromacs_energy_time_series/
//     coulomb_sr, coulomb_recip, coulomb_14            (T,)  float64  kJ/mol
//     bond, angle, urey_bradley                        (T,)  float64  kJ/mol
//     proper_dih, improper_dih, cmap_dih               (T,)  float64  kJ/mol
//     lj_sr, lj_14, disper_corr                        (T,)  float64  kJ/mol
//     potential, kinetic, total_energy, enthalpy       (T,)  float64  kJ/mol
//     temperature, T_protein, T_non_protein            (T,)  float64  K
//     pressure                                         (T,)  float64  bar
//     volume                                           (T,)  float64  nm^3
//     density                                          (T,)  float64  kg/m^3
//     box_x, box_y, box_z                              (T,)  float64  nm
//     virial                                           (T, 9) float64 kJ/mol
//     pressure_tensor                                  (T, 9) float64 bar
//     frame_indices                                    (T,)  uint64
//     frame_times                                      (T,)  float64  ps
//     attrs: result_name, n_frames, finalized, units = "kJ/mol"
//            tensor_layout = "XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ"
//
// Dependencies: GromacsEnergyResult (conformation-scope, unconditionally
// attached in PerFrameExtractionSet — no source-attached gate needed).
//

#include "TrajectoryResult.h"
#include "GromacsEnergyResult.h"

#include <cstddef>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class GromacsEnergyTimeSeriesTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "GromacsEnergyTimeSeriesTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<GromacsEnergyTimeSeriesTrajectoryResult> Create(
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
    // Per-frame GromacsEnergy snapshots, appended in Compute. Each entry
    // is a copy of the frame's GromacsEnergyResult::Energy() — system
    // scalars, NOT per-atom data.
    std::vector<GromacsEnergy> per_frame_energy_;
    std::vector<std::size_t>   frame_indices_;
    std::vector<double>        frame_times_;
    std::size_t                n_frames_ = 0;
    bool                       finalized_ = false;
};

}  // namespace nmr

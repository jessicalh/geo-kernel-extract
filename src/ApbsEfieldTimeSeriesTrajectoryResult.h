#pragma once
//
// ApbsEfieldTimeSeriesTrajectoryResult: per-atom per-frame time
// series of the canonical APBS reaction E-field
// (ConformationAtom::apbs_efield, V/Angstrom). FO dense-buffer
// pattern with Vec3 payload (clones PositionsTimeSeriesTrajectoryResult
// shape against the apbs_efield field). A defensive source gate records
// absence explicitly with NaN payloads and a per-frame mask.
//
// Emission:
//
//   /trajectory/apbs_efield_time_series/
//     xyz            (N, T, 3)  float64  — E_x, E_y, E_z (V/A)
//     clamp_mask     (N, T)     uint8
//     clamp_scale    (N, T)     float64
//     frame_indices  (T,)       uint64
//     frame_times    (T,)       float64
//     source_attached_per_frame (T,) uint8
//     attrs:
//       result_name    = "ApbsEfieldTimeSeriesTrajectoryResult"
//       irrep_layout   = "x,y,z"
//       normalization  = "cartesian"
//       parity         = "1o"
//       units          = "V/Angstrom"
//       n_atoms, n_frames, finalized
//
// Parity "1o": E-field is a polar (true) vector, parity-odd under
// inversion. Same parity convention as PositionsTimeSeries.
//

#include "DenseBuffer.h"
#include "TrajectoryResult.h"
#include "Types.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class ApbsEfieldTimeSeriesTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "ApbsEfieldTimeSeriesTrajectoryResult";
    }

    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<ApbsEfieldTimeSeriesTrajectoryResult> Create(
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
    std::vector<std::vector<Vec3>> per_atom_efield_;
    std::vector<std::vector<std::uint8_t>> per_atom_clamp_mask_;
    std::vector<std::vector<double>> per_atom_clamp_scale_;
    std::vector<std::uint8_t> clamp_mask_flat_;
    std::vector<double> clamp_scale_flat_;
    std::vector<std::uint8_t> source_attached_per_frame_;
    std::vector<std::array<std::uint64_t, 3>> grid_dims_per_frame_;
    std::vector<std::array<double, 3>> grid_lengths_A_per_frame_;
    std::vector<std::array<double, 3>> grid_origin_A_per_frame_;
    std::vector<std::array<double, 3>> grid_spacing_A_per_frame_;
    double apbs_manual_grid_padding_A_ =
        std::numeric_limits<double>::quiet_NaN();
    double apbs_manual_grid_min_dim_A_ =
        std::numeric_limits<double>::quiet_NaN();
    double apbs_temperature_K_ =
        std::numeric_limits<double>::quiet_NaN();
    double apbs_thermal_voltage_V_ =
        std::numeric_limits<double>::quiet_NaN();
    double efield_clamp_threshold_ =
        std::numeric_limits<double>::quiet_NaN();
    std::vector<std::size_t> frame_indices_;
    std::vector<double> frame_times_;
    std::size_t n_frames_ = 0;
    bool finalized_ = false;
};

}  // namespace nmr

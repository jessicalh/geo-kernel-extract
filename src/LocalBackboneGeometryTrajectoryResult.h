#pragma once
//
// LocalBackboneGeometryTrajectoryResult: residue-local backbone valence
// geometry over the frames actually dispatched by Trajectory::Run.
//
// H5-only trajectory result. It owns no ConformationAtom or TrajectoryAtom
// fields and emits no static/per-frame NPY arrays. Per-residue buffers are
// residue-major and retain one value for every dispatched frame:
//
//   /trajectory/local_backbone_geometry/
//     tau_N_CA_C          (R,T)   radians
//     angle_N_CA_CB       (R,T)   radians
//     angle_CB_CA_C       (R,T)   radians
//     angle_Cprev_N_CA    (R,T)   radians
//     angle_CA_C_Nnext    (R,T)   radians
//     cb_deviation        (R,T)   Angstrom
//     cb_local_vector     (R,T,3) Angstrom, observed_CB - ideal_CB in
//                                  global Cartesian coordinates
//
// Static masks are topology-derived once at Create. Backbone adjacency is
// exclusively Protein::BackbonePredecessor/BackboneSuccessor; residue-index
// arithmetic is deliberately absent. Positions are unconditional trajectory
// substrate, so source_attached_per_frame is all one.

#include "TrajectoryResult.h"
#include "Types.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class LocalBackboneGeometryTrajectoryResult : public TrajectoryResult {
public:
    std::string Name() const override {
        return "LocalBackboneGeometryTrajectoryResult";
    }

    // Positions and invariant typed topology are available after Seed.
    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<LocalBackboneGeometryTrajectoryResult>
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
    std::vector<std::vector<double>> tau_n_ca_c_;
    std::vector<std::vector<double>> angle_n_ca_cb_;
    std::vector<std::vector<double>> angle_cb_ca_c_;
    std::vector<std::vector<double>> angle_cprev_n_ca_;
    std::vector<std::vector<double>> angle_ca_c_nnext_;
    std::vector<std::vector<double>> cb_deviation_;
    std::vector<std::vector<Vec3>> cb_local_vector_;

    // Static per-residue topology/type masks.
    std::vector<std::uint8_t> has_n_ca_c_;
    std::vector<std::uint8_t> has_cb_;
    std::vector<std::uint8_t> has_prev_c_;
    std::vector<std::uint8_t> has_next_n_;
    std::vector<std::uint8_t> is_glycine_;
    std::vector<std::uint8_t> is_proline_;

    // Dispatched-frame metadata; T is exactly these vectors' length.
    std::vector<std::size_t> frame_indices_;
    std::vector<double> frame_times_;
    std::vector<std::uint8_t> source_attached_per_frame_;

    std::size_t n_frames_ = 0;
    bool finalized_ = false;
};

}  // namespace nmr

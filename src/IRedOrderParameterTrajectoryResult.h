#pragma once
//
// IRedOrderParameterTrajectoryResult: isotropic Reorientational Eigenmode
// Dynamics (iRED) order parameters for the backbone amide N-H vectors,
// across the whole protein. Prompers & Bruschweiler 2002 (J. Am. Chem.
// Soc. 124, 4522). iRED removes overall tumbling by isotropic averaging
// rather than by superposing onto a reference, so it has no
// reference-structure ambiguity -- the robust whole-protein S^2 and the
// reference-free cross-check to the superposition-based S^2 in
// ReorientationalDynamicsTrajectoryResult.
//
// Method: form the isotropically-averaged covariance of the rank-2
// orientations of the M interaction vectors (lab frame),
//
//   M_ij = < P2( u_i(t) . u_j(t) ) >,    M_ii = 1,
//
// and diagonalise it, M = sum_m lambda_m |m><m|. For a rank-2 (l=2)
// interaction, isotropic overall tumbling spans the FIVE LARGEST
// eigenmodes; the per-vector order parameter is the projection ONTO those
// tumbling modes:
//
//   S^2_i = sum_{m in 5 largest} lambda_m * <i|m>^2
//         = 1 - sum_{m in the rest} lambda_m * <i|m>^2     (since M_ii = 1)
//
// (Rigid body: internal modes carry no weight -> S^2 -> 1. This is the
// codex-corrected orientation; the inverse 1 - sum(5 largest) would give
// S^2 = 0 for a rigid body.) The eigenvalue gap lambda5/lambda6 diagnoses
// whether the 5-mode overall/internal split is clean.
//
// Scope: v1 is the amide N-H set (res.H present, non-Pro) -- the canonical
// iRED application and the standard 15N relaxation probe. ReorientationalDynamics
// carries the broader per-vector set (Calpha-Halpha, C=O, methyl axis,
// aromatic C-H). iRED over all of those at once is a documented follow-up.
//
// Lifecycle: FO. The M x M matrix accumulates online each frame (no
// history); Finalize divides by the frame count and eigendecomposes.
// Cost O(M^2) per frame + O(M^3) once; M ~ residue count, trivial.
//
// Emission /trajectory/ired_order_parameters/:
//   s2_ired           (M,) float64   per-vector iRED order parameter
//   eigenvalues       (M,) float64   descending; first 5 are overall tumbling
//   residue_index     (M,) int32     parent residue of each N-H vector
//   n_atom, h_atom    (M,) int32     the N and H atom indices defining each vector
//   attrs: result_name, n_vectors, n_frames, finalized, n_tumbling_modes=5,
//          separability_gap (lambda5/lambda6), reference="none_isotropic",
//          frame="lab"
//

#include "TrajectoryResult.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class IRedOrderParameterTrajectoryResult : public TrajectoryResult {
public:
    // Number of rank-2 overall-tumbling eigenmodes (2l+1 for l=2).
    static constexpr std::size_t N_TUMBLING_MODES = 5;

    std::string Name() const override {
        return "IRedOrderParameterTrajectoryResult";
    }

    // No declared dependency: positions and Residue.N/H are present after
    // tp.Seed (PATTERNS.md 15).
    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<IRedOrderParameterTrajectoryResult> Create(
        const TrajectoryProtein& tp);

    void Compute(const ProteinConformation& conf,
                 TrajectoryProtein& tp,
                 Trajectory& traj,
                 std::size_t frame_idx,
                 double time_ps) override;

    void Finalize(TrajectoryProtein& tp, Trajectory& traj) override;

    void WriteH5Group(const TrajectoryProtein& tp,
                      HighFive::File& file) const override;

    std::size_t NumVectors() const { return n_vectors_; }

private:
    // Per-vector identity (amide N-H), captured at Create.
    std::vector<std::size_t>  n_atom_;        // backbone N index
    std::vector<std::size_t>  h_atom_;        // amide H index
    std::vector<std::int32_t> residue_index_; // parent residue

    // M x M covariance accumulator, row-major flat, summed each frame.
    std::vector<double> m_accum_;

    // Finalized.
    std::vector<double> s2_ired_;       // M
    std::vector<double> eigenvalues_;   // M, descending
    double separability_gap_ = 0.0;     // lambda5 / lambda6 (NaN if M <= 5)

    std::size_t n_vectors_ = 0;
    std::size_t n_frames_  = 0;
    bool finalized_ = false;
};

}  // namespace nmr

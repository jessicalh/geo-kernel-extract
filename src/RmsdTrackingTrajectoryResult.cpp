#include "RmsdTrackingTrajectoryResult.h"

#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "TrajectoryProtein.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <Eigen/Dense>
#include <Eigen/SVD>

#include <cmath>
#include <limits>
#include <typeinfo>

namespace nmr {

namespace {

// Kabsch alignment of N points: returns the RMSD between
// `current` (rotated to best-fit) and `reference` (the target).
// Both inputs are assembled as 3 × M Eigen matrices with one point per
// column; we use P = current_centered, Q = reference_centered, H = P Q^T.
//
// RMSD = sqrt(mean(|R*(p_i - p_centroid) + q_centroid - q_i|^2))
//      = sqrt(mean(|R*(p_i - p_centroid) - (q_i - q_centroid)|^2))
//
// Degenerate cases:
//   - M < 3:           RMSD undefined (rotation underdetermined).
//                      Returns NaN.
//   - Collapsed point clouds are not special-cased; after centroid
//                      subtraction, identical zero-variance clouds give
//                      zero centered residual.
//
double KabschRmsd(const std::vector<Vec3>& current,
                  const std::vector<Vec3>& reference) {
    const std::size_t M = current.size();
    if (reference.size() != M) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (M < 3) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    // Centroids.
    Vec3 p_cen = Vec3::Zero();
    Vec3 q_cen = Vec3::Zero();
    for (std::size_t i = 0; i < M; ++i) {
        p_cen += current[i];
        q_cen += reference[i];
    }
    p_cen /= static_cast<double>(M);
    q_cen /= static_cast<double>(M);

    // Centered matrices: P (3, M), Q (3, M).
    Eigen::Matrix<double, 3, Eigen::Dynamic> P(3, M);
    Eigen::Matrix<double, 3, Eigen::Dynamic> Q(3, M);
    for (std::size_t i = 0; i < M; ++i) {
        P.col(i) = current[i] - p_cen;
        Q.col(i) = reference[i] - q_cen;
    }

    // Cross-covariance H = P * Q^T (3 x 3).
    const Mat3 H = P * Q.transpose();

    // SVD with reflection correction. Canonical Kabsch uses
    // sign(det(V*Uᵀ)) — not the raw determinant — because the
    // product is theoretically orthogonal with |det| = 1, so the
    // reflection guard only needs the sign. Passing the raw det
    // (≈ ±1 + O(ε)) into the diagonal would inject anisotropic
    // O(ε) scaling and break strict orthogonality of R.
    Eigen::JacobiSVD<Mat3> svd(H,
        Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Mat3& U = svd.matrixU();
    const Mat3& V = svd.matrixV();
    const double det = (V * U.transpose()).determinant();
    Eigen::DiagonalMatrix<double, 3> D(1.0, 1.0,
        (det < 0.0) ? -1.0 : 1.0);
    const Mat3 R = V * D * U.transpose();

    // RMSD after alignment.
    double sum_sq = 0.0;
    for (std::size_t i = 0; i < M; ++i) {
        const Vec3 aligned = R * P.col(i);  // already centered
        const Vec3 ref_centered = Q.col(i);
        sum_sq += (aligned - ref_centered).squaredNorm();
    }
    return std::sqrt(sum_sq / static_cast<double>(M));
}

}  // namespace

std::unique_ptr<RmsdTrackingTrajectoryResult>
RmsdTrackingTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<RmsdTrackingTrajectoryResult>();

    // Build alignment set from typed Residue backbone slots: N, CA, C, O.
    // The canonical heavy-atom backbone set per IUPAC + GROMACS
    // `g_rms -fit backbone`; HA is intentionally excluded since it is
    // not part of the peptide backbone proper, and any HA motion is
    // captured implicitly via CA.
    // Per residue, include each slot independently — ACE caps lack
    // N/CA but have C/O; NME caps lack C/O but have N. The mix is
    // fine — RMSD is over heavy backbone atoms wherever they are.
    const Protein& protein = tp.ProteinRef();
    const std::size_t n_res = protein.ResidueCount();
    for (std::size_t ri = 0; ri < n_res; ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        if (res.N  != Residue::NONE) r->atom_indices_.push_back(res.N);
        if (res.CA != Residue::NONE) r->atom_indices_.push_back(res.CA);
        if (res.C  != Residue::NONE) r->atom_indices_.push_back(res.C);
        if (res.O  != Residue::NONE) r->atom_indices_.push_back(res.O);
    }

    OperationLog::Info(
        "RmsdTrackingTrajectoryResult::Create",
        "alignment set: " + std::to_string(r->atom_indices_.size()) +
        " backbone heavy atoms (N/CA/C/O) across " +
        std::to_string(n_res) + " residues");

    return r;
}

void RmsdTrackingTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;

    const std::size_t M = atom_indices_.size();

    // First frame: capture the reference geometry.
    if (!reference_captured_) {
        reference_positions_.reserve(M);
        for (std::size_t ai : atom_indices_) {
            reference_positions_.push_back(conf.PositionAt(ai));
        }
        reference_captured_ = true;
        // Frame-0 RMSD vs itself is 0 by definition.
        rmsd_.push_back(0.0);
    } else {
        std::vector<Vec3> current;
        current.reserve(M);
        for (std::size_t ai : atom_indices_) {
            current.push_back(conf.PositionAt(ai));
        }
        rmsd_.push_back(KabschRmsd(current, reference_positions_));
    }

    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    source_attached_per_frame_.push_back(1u);
    ++n_frames_;
}

void RmsdTrackingTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                              Trajectory& traj) {
    (void)tp; (void)traj;
    // AV pattern: data is already accumulated. Idempotent.
    finalized_ = true;
}

double RmsdTrackingTrajectoryResult::LatestRmsd() const {
    if (rmsd_.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return rmsd_.back();
}

double RmsdTrackingTrajectoryResult::RmsdAtSampleIndex(
        std::size_t sample_idx) const {
    if (sample_idx >= rmsd_.size()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return rmsd_[sample_idx];
}

void RmsdTrackingTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    (void)tp;

    auto grp = file.createGroup("/trajectory/rmsd_tracking");
    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_atoms", atom_indices_.size());
    grp.createAttribute("n_frames", n_frames_);
    grp.createAttribute("finalized", finalized_);
    grp.createAttribute("alignment_method", std::string("kabsch_svd"));
    grp.createAttribute("atom_selection",
        std::string("backbone_heavy_atoms_NCACO"));
    grp.createAttribute("reference_frame_origin",
        std::string("trajectory_frame_0"));
    grp.createAttribute("units", std::string("Angstrom"));
    grp.createAttribute("source_attached_policy",
        std::string("always_attached -- positions present at tp.Seed; "
                    "reference geometry captured at first Compute call "
                    "(frame 0)"));
    grp.createAttribute("rmsd_frame_0_convention",
        std::string("0.0 exactly -- frame 0 is its own reference"));

    grp.createDataSet("rmsd", rmsd_);

    // atom_indices for alignment set provenance + SDK / analysis
    // bridge to the protein atom axis.
    std::vector<std::int32_t> ai32(atom_indices_.size());
    for (std::size_t i = 0; i < atom_indices_.size(); ++i) {
        ai32[i] = static_cast<std::int32_t>(atom_indices_[i]);
    }
    grp.createDataSet("atom_indices", ai32);

    grp.createDataSet("frame_indices", frame_indices_);
    grp.createDataSet("frame_times", frame_times_);
    grp.createDataSet("source_attached_per_frame",
                       source_attached_per_frame_);
}

}  // namespace nmr

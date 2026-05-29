#include "IRedOrderParameterTrajectoryResult.h"

#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "TrajectoryProtein.h"
#include "TrajectorySpectral.h"   // P2
#include "Types.h"                // Vec3

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <limits>

namespace nmr {


std::unique_ptr<IRedOrderParameterTrajectoryResult>
IRedOrderParameterTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<IRedOrderParameterTrajectoryResult>();
    const Protein& protein = tp.ProteinRef();
    const std::size_t R = protein.ResidueCount();

    // Amide N-H vectors: residues whose backbone N and amide H are both
    // present (Residue.H is NONE for Pro and for residues without a
    // resolved amide proton).
    for (std::size_t ri = 0; ri < R; ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        if (res.N != Residue::NONE && res.H != Residue::NONE) {
            r->n_atom_.push_back(res.N);
            r->h_atom_.push_back(res.H);
            r->residue_index_.push_back(static_cast<std::int32_t>(ri));
        }
    }
    r->n_vectors_ = r->n_atom_.size();
    r->m_accum_.assign(r->n_vectors_ * r->n_vectors_, 0.0);

    OperationLog::Info("IRedOrderParameterTrajectoryResult::Create",
        "amide N-H iRED set: " + std::to_string(r->n_vectors_) +
        " vectors across " + std::to_string(R) + " residues");
    return r;
}


// ── Compute ──────────────────────────────────────────────────────
//
// Accumulate the lab-frame rank-2 covariance M_ij += P2(u_i . u_j) over
// all vector pairs this frame. Online (no history); the matrix is the
// running sum, divided by the frame count at Finalize.

void IRedOrderParameterTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj; (void)frame_idx; (void)time_ps;
    const std::size_t M = n_vectors_;
    ++n_frames_;
    if (M == 0) return;

    std::vector<Vec3> u(M);
    for (std::size_t i = 0; i < M; ++i) {
        const Vec3 d = conf.PositionAt(h_atom_[i]) - conf.PositionAt(n_atom_[i]);
        const double nrm = d.norm();
        if (nrm > 1e-9) u[i] = d / nrm;       // N-H ~1 A, never 0
        else            u[i] = Vec3::Zero();
    }
    for (std::size_t i = 0; i < M; ++i) {
        m_accum_[i * M + i] += 1.0;          // P2(u_i . u_i) = 1 exactly
        for (std::size_t j = i + 1; j < M; ++j) {
            const double p = P2(u[i].dot(u[j]));
            m_accum_[i * M + j] += p;
            m_accum_[j * M + i] += p;
        }
    }
}


// ── Finalize ──────────────────────────────────────────────────────
//
// M_ij /= T, then symmetric eigendecomposition. Eigen sorts eigenvalues
// ASCENDING, so the five largest (overall tumbling) are the last five.
// S^2_i = sum over those five of lambda_m * <i|m>^2 -- the projection onto
// the tumbling modes (rigid body -> S^2 = 1). separability_gap = lambda5 /
// lambda6 diagnoses the overall/internal split.

void IRedOrderParameterTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                  Trajectory& traj) {
    (void)tp; (void)traj;
    if (finalized_) return;
    const std::size_t M = n_vectors_;
    const double nan_val = std::nan("");
    s2_ired_.assign(M, nan_val);
    eigenvalues_.assign(M, nan_val);
    separability_gap_ = nan_val;

    if (M == 0 || n_frames_ < 2) { finalized_ = true; return; }

    Eigen::MatrixXd A(M, M);
    const double T = static_cast<double>(n_frames_);
    for (std::size_t i = 0; i < M; ++i)
        for (std::size_t j = 0; j < M; ++j)
            A(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) =
                m_accum_[i * M + j] / T;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
    if (es.info() != Eigen::Success) {
        OperationLog::Warn("IRedOrderParameterTrajectoryResult::Finalize",
                           "eigendecomposition did not converge");
        finalized_ = true;
        return;
    }
    const Eigen::VectorXd& evals = es.eigenvalues();      // ascending
    const Eigen::MatrixXd& evecs = es.eigenvectors();     // columns = modes

    const std::size_t ntumb = std::min<std::size_t>(N_TUMBLING_MODES, M);
    for (std::size_t i = 0; i < M; ++i) {
        double s = 0.0;
        for (std::size_t t = 0; t < ntumb; ++t) {
            const Eigen::Index m = static_cast<Eigen::Index>(M - 1 - t);
            const double comp = evecs(static_cast<Eigen::Index>(i), m);
            s += evals[m] * comp * comp;
        }
        s2_ired_[i] = s;
    }
    for (std::size_t r = 0; r < M; ++r)
        eigenvalues_[r] = evals[static_cast<Eigen::Index>(M - 1 - r)];  // descending

    if (M > N_TUMBLING_MODES) {
        const double l5 = evals[static_cast<Eigen::Index>(M - 5)];
        const double l6 = evals[static_cast<Eigen::Index>(M - 6)];
        if (l6 > 1e-15)      separability_gap_ = l5 / l6;
        else if (l5 > 1e-15) separability_gap_ =
                                 std::numeric_limits<double>::infinity();  // clean rank-5 split
        else                 separability_gap_ = nan_val;
    }

    std::vector<double>().swap(m_accum_);
    finalized_ = true;
    OperationLog::Info(LogCalcOther,
        "IRedOrderParameterTrajectoryResult::Finalize",
        "iRED over " + std::to_string(M) + " amide N-H vectors, " +
        std::to_string(n_frames_) + " frames; separability_gap = " +
        std::to_string(separability_gap_));
}


// ── WriteH5Group ──────────────────────────────────────────────────

void IRedOrderParameterTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    (void)tp;
    if (!finalized_) {
        OperationLog::Warn("IRedOrderParameterTrajectoryResult::WriteH5Group",
                           "Finalize not called; nothing to write");
        return;
    }
    if (n_vectors_ == 0) {
        OperationLog::Warn("IRedOrderParameterTrajectoryResult::WriteH5Group",
                           "no amide N-H vectors; group omitted");
        return;
    }
    const std::size_t M = n_vectors_;

    auto grp = file.createGroup("/trajectory/ired_order_parameters");
    grp.createAttribute("result_name",      Name());
    grp.createAttribute("n_vectors",        M);
    grp.createAttribute("n_frames",         n_frames_);
    grp.createAttribute("finalized",        finalized_);
    grp.createAttribute("n_tumbling_modes",
                        static_cast<std::size_t>(N_TUMBLING_MODES));
    grp.createAttribute("separability_gap", separability_gap_);
    grp.createAttribute("vector_set",       std::string("amide_N_H"));
    grp.createAttribute("reference",        std::string("none_isotropic"));
    grp.createAttribute("frame",            std::string("lab"));
    grp.createAttribute("s2_definition", std::string(
        "S^2_i = sum over the 5 largest eigenmodes of lambda_m * <i|m>^2 "
        "(projection onto the overall-tumbling modes; rigid -> 1). "
        "Prompers & Bruschweiler 2002, J. Am. Chem. Soc. 124, 4522."));
    grp.createAttribute("separability_note", std::string(
        "lambda5/lambda6: a large gap means the 5-mode overall/internal "
        "split is clean; near 1 means overall and internal motion mix "
        "(common for flexible/multi-domain systems). NaN if n_vectors <= 5."));

    grp.createDataSet("s2_ired", s2_ired_)
       .createAttribute("units", std::string("dimensionless"));
    grp.createDataSet("eigenvalues", eigenvalues_)
       .createAttribute("note", std::string("descending; first 5 are overall tumbling"));
    grp.createDataSet("residue_index", residue_index_)
       .createAttribute("units", std::string("residue_index"));

    std::vector<std::int32_t> n_i32(M), h_i32(M);
    for (std::size_t i = 0; i < M; ++i) {
        n_i32[i] = static_cast<std::int32_t>(n_atom_[i]);
        h_i32[i] = static_cast<std::int32_t>(h_atom_[i]);
    }
    grp.createDataSet("n_atom", n_i32)
       .createAttribute("note", std::string("backbone N atom index of each N-H vector"));
    grp.createDataSet("h_atom", h_i32)
       .createAttribute("note", std::string("amide H atom index of each N-H vector"));
}

}  // namespace nmr

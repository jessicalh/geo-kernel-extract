// TransformedConformation — implementation.

#include "TransformedConformation.h"

#include "QtAtom.h"
#include "QtProtein.h"

#include "../diagnostics/ThreadGuard.h"

#include <Eigen/SVD>

#include <QLoggingCategory>

#include <algorithm>
#include <cstddef>
#include <utility>

namespace h5reader::model {

namespace {
Q_LOGGING_CATEGORY(cXform, "h5reader.transform")
}  // namespace

TransformedConformation::TransformedConformation(Conformation* inner, QObject* parent)
    : Conformation(inner ? inner->protein() : nullptr),
      inner_(inner) {
    setObjectName(QStringLiteral("TransformedConformation"));
    if (parent)
        setParent(parent);
}

TransformedConformation::~TransformedConformation() = default;

std::size_t TransformedConformation::frameCount() const {
    return inner_ ? inner_->frameCount() : 0;
}

double TransformedConformation::timePicoseconds(std::size_t frame) const {
    return inner_ ? inner_->timePicoseconds(frame) : 0.0;
}

std::size_t TransformedConformation::originalFrameIndex(std::size_t frame) const {
    return inner_ ? inner_->originalFrameIndex(frame) : frame;
}

const TrajectoryConformation* TransformedConformation::asTrajectory() const {
    return inner_ ? inner_->asTrajectory() : nullptr;
}

Vec3 TransformedConformation::atomPosition(std::size_t frame, std::size_t atomIdx) const {
    if (!inner_)
        return Vec3::Zero();
    if (mode_ == Mode::Identity)
        return inner_->atomPosition(frame, atomIdx);

    // Look up or compute the transform for this frame.
    auto it = transformCache_.find(frame);
    if (it == transformCache_.end()) {
        const Transform3D tform = computeTransform(frame);
        it = transformCache_.emplace(frame, tform).first;
    }
    const Transform3D& t = it->second;
    const Vec3 raw = inner_->atomPosition(frame, atomIdx);
    return t.R * raw + t.T;
}

std::shared_ptr<const QtConformationSnapshot>
TransformedConformation::loadSnapshot(std::size_t frame) {
    // Snapshots are full-fidelity per-frame source data — atom positions
    // PLUS calculator NPYs. We do NOT decorate snapshots: consumers that
    // need calculator data already read from the snapshot directly, and
    // applying our 3x3 transform to a snapshot's Pos column would diverge
    // from the inner conformation's atomPosition seam. Forward unchanged.
    return inner_ ? inner_->snapshot(frame) : nullptr;
}

void TransformedConformation::setMode(Mode mode,
                                       std::size_t referenceFrame,
                                       std::vector<std::size_t> subsetAtoms) {
    ASSERT_THREAD(this);

    mode_ = mode;
    referenceFrame_ = referenceFrame;
    subsetAtoms_ = std::move(subsetAtoms);

    // Bump the generation counter (semantic marker, primarily for logs /
    // future diagnostics) and wipe the per-frame transform cache so the
    // next atomPosition() lookup recomputes.
    ++generation_;
    transformCache_.clear();
    referencePositions_.clear();

    // Pre-load reference positions for the fit modes so each per-frame
    // Kabsch call doesn't re-scan the inner conformation. For
    // FitReference we need every atom; for FitSubset just the subset.
    if (inner_ && protein_) {
        const std::size_t atomCount = protein_->atomCount();
        if (mode_ == Mode::FitReference) {
            referencePositions_.reserve(atomCount);
            for (std::size_t a = 0; a < atomCount; ++a)
                referencePositions_.push_back(inner_->atomPosition(referenceFrame_, a));
        } else if (mode_ == Mode::FitSubset) {
            referencePositions_.reserve(subsetAtoms_.size());
            for (std::size_t a : subsetAtoms_) {
                if (a >= atomCount) {
                    qCWarning(cXform).noquote()
                        << "FitSubset atom" << a << "out of range; clamping subset to in-range atoms";
                    continue;
                }
                referencePositions_.push_back(inner_->atomPosition(referenceFrame_, a));
            }
        }
    }

    const char* modeName = "identity";
    switch (mode_) {
        case Mode::Identity:     modeName = "identity"; break;
        case Mode::CenterCom:    modeName = "center_com"; break;
        case Mode::FitReference: modeName = "fit_reference"; break;
        case Mode::FitSubset:    modeName = "fit_subset"; break;
    }
    qCInfo(cXform).noquote()
        << "mode set to" << modeName
        << "| ref_frame=" << static_cast<qlonglong>(referenceFrame_)
        << "| subset_size=" << static_cast<qlonglong>(subsetAtoms_.size())
        << "| generation=" << static_cast<qlonglong>(generation_);

    emit transformChanged();
}

std::vector<std::size_t>
TransformedConformation::BackboneSubset(const QtProtein& protein) {
    std::vector<std::size_t> out;
    out.reserve(protein.atomCount() / 4);  // rough estimate (N, CA, C, O per residue)
    for (std::size_t i = 0; i < protein.atomCount(); ++i) {
        if (protein.atom(i).IsBackbone())
            out.push_back(i);
    }
    return out;
}

TransformedConformation::Transform3D
TransformedConformation::computeTransform(std::size_t frame) const {
    Transform3D out;  // identity by default

    if (!inner_ || !protein_)
        return out;
    const std::size_t atomCount = protein_->atomCount();
    if (atomCount == 0)
        return out;

    switch (mode_) {
        case Mode::Identity:
            return out;

        case Mode::CenterCom: {
            // T = -COM(frame); R = identity. Equal-weight COM over all atoms
            // (mass-weighted COM would need atomic masses; equal-weight is
            // the simplest and addresses the "protein drifts across the box"
            // case the centroid-follow code was trying to handle).
            Vec3 com = Vec3::Zero();
            for (std::size_t a = 0; a < atomCount; ++a)
                com += inner_->atomPosition(frame, a);
            com /= static_cast<double>(atomCount);
            out.T = -com;
            return out;
        }

        case Mode::FitReference: {
            // Kabsch over ALL atoms. The reference positions were cached
            // at setMode() time; pull the current frame's positions and
            // hand the two equal-length vectors to KabschFit.
            if (referencePositions_.size() != atomCount)
                return out;  // cache wasn't built; safest is identity
            std::vector<Vec3> current;
            current.reserve(atomCount);
            for (std::size_t a = 0; a < atomCount; ++a)
                current.push_back(inner_->atomPosition(frame, a));
            return KabschFit(current, referencePositions_);
        }

        case Mode::FitSubset: {
            // Kabsch over the subset indices. Reference positions are the
            // subset atoms at referenceFrame_, current is the subset at
            // `frame`. Use the smaller of (subset.size(), reference.size())
            // to defend against subset/reference shape mismatches.
            const std::size_t n = std::min(subsetAtoms_.size(), referencePositions_.size());
            if (n < 3)
                return out;  // underdetermined
            std::vector<Vec3> current;
            current.reserve(n);
            for (std::size_t i = 0; i < n; ++i) {
                const std::size_t a = subsetAtoms_[i];
                if (a >= atomCount)
                    continue;
                current.push_back(inner_->atomPosition(frame, a));
            }
            // KabschFit is robust to mismatched length: it requires the
            // SAME length on both sides — clip the reference to current.size()
            // (the same atoms in the same order; per-atom in subsetAtoms_).
            std::vector<Vec3> reference(referencePositions_.begin(),
                                         referencePositions_.begin() + current.size());
            return KabschFit(current, reference);
        }
    }
    return out;
}

// Kabsch algorithm (W. Kabsch, Acta Cryst. A32, 922-923, 1976) — the
// classical least-squares rigid-body fit. Output: R (3x3), T (3) such
// that for each pair (current[i], reference[i]) the sum of squared
// distances || R * current[i] + T - reference[i] ||^2 is minimised.
//
// Steps:
//   1. centroid_c = mean(current), centroid_r = mean(reference)
//   2. P = current - centroid_c, Q = reference - centroid_r (N x 3 each)
//   3. H = P^T * Q (3 x 3 covariance)
//   4. SVD: H = U * S * V^T
//   5. d = sign(det(V * U^T))   -- determinant correction to avoid
//                                  reflection
//   6. R = V * diag(1, 1, d) * U^T
//   7. T = centroid_r - R * centroid_c
//
// We use Eigen::JacobiSVD with FullU/V; Eigen's SVD ordering puts
// singular values in descending order so U/V columns are already
// sorted by importance. Both inputs are assumed to be the SAME atom
// set in the SAME order at the SAME atom-count.
TransformedConformation::Transform3D
TransformedConformation::KabschFit(const std::vector<Vec3>& current,
                                    const std::vector<Vec3>& reference) {
    Transform3D out;
    const std::size_t n = current.size();
    if (n != reference.size() || n < 3)
        return out;

    // Centroids
    Vec3 cc = Vec3::Zero();
    Vec3 cr = Vec3::Zero();
    for (std::size_t i = 0; i < n; ++i) {
        cc += current[i];
        cr += reference[i];
    }
    cc /= static_cast<double>(n);
    cr /= static_cast<double>(n);

    // Centred coordinate matrices P, Q (3 x N — we use the columnar
    // layout because Eigen's default storage order favours columns).
    Eigen::MatrixXd P(3, n);
    Eigen::MatrixXd Q(3, n);
    for (std::size_t i = 0; i < n; ++i) {
        P.col(static_cast<Eigen::Index>(i)) = current[i] - cc;
        Q.col(static_cast<Eigen::Index>(i)) = reference[i] - cr;
    }

    // Covariance H = P * Q^T (3 x 3).
    Eigen::Matrix3d H = P * Q.transpose();

    // SVD with full U + V so we can compose R = V * diag(1, 1, d) * U^T.
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Eigen::Matrix3d& U = svd.matrixU();
    const Eigen::Matrix3d& V = svd.matrixV();

    // Determinant correction: if det(V * U^T) < 0 we have a reflection
    // (improper rotation); flip the sign on the last column of V to
    // get a proper rotation. This is the canonical Kabsch fix-up.
    Eigen::Matrix3d D = Eigen::Matrix3d::Identity();
    const double d = (V * U.transpose()).determinant();
    D(2, 2) = (d < 0.0) ? -1.0 : 1.0;

    const Eigen::Matrix3d R = V * D * U.transpose();
    const Eigen::Vector3d T = cr - R * cc;

    out.R = R;
    out.T = T;
    return out;
}

}  // namespace h5reader::model

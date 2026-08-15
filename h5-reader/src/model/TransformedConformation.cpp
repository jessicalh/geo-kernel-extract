// TransformedConformation — implementation.

#include "TransformedConformation.h"

#include "QtAtom.h"
#include "QtProtein.h"
#include "ScientificAlignment.h"

#include "../app/FitTargetMath.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"

#include <Eigen/Geometry>

#include <QElapsedTimer>
#include <QLoggingCategory>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace h5reader::model {

namespace {
Q_LOGGING_CATEGORY(cXform, "h5reader.transform")

constexpr int kMeanReferenceMaxIterations = 20;
constexpr double kMeanReferenceEpsAngstrom = 1e-4;

Eigen::Quaterniond NormalisedQuaternion(const Mat3& rotation) {
    constexpr double kMinNorm = 1e-12;
    Eigen::Quaterniond q(rotation);
    if (q.norm() < kMinNorm)
        return Eigen::Quaterniond::Identity();
    q.normalize();
    return q;
}

Vec3 Centroid(const std::vector<Vec3>& positions) {
    Vec3 c = Vec3::Zero();
    if (positions.empty())
        return c;
    for (const Vec3& p : positions)
        c += p;
    return c / static_cast<double>(positions.size());
}

std::vector<Vec3> DemeanedCopy(const std::vector<Vec3>& positions, const Vec3& centroid) {
    std::vector<Vec3> out;
    out.reserve(positions.size());
    for (const Vec3& p : positions)
        out.push_back(p - centroid);
    return out;
}

double RmsDifference(const std::vector<Vec3>& a, const std::vector<Vec3>& b) {
    if (a.empty() || a.size() != b.size())
        return 0.0;
    double ss = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
        ss += (a[i] - b[i]).squaredNorm();
    return std::sqrt(ss / static_cast<double>(a.size()));
}
}  // namespace

TransformedConformation::TransformedConformation(Conformation* inner, QObject* parent)
    : Conformation(inner ? inner->protein() : nullptr),
      inner_(inner) {
    CENSUS_REGISTER(this);
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

std::optional<std::size_t>
TransformedConformation::frameRowForOriginalIndex(std::size_t originalFrame) const {
    return inner_ ? inner_->frameRowForOriginalIndex(originalFrame) : std::nullopt;
}

const TrajectoryConformation* TransformedConformation::asTrajectory() const {
    return inner_ ? inner_->asTrajectory() : nullptr;
}

Vec3 TransformedConformation::atomPosition(std::size_t frame, std::size_t atomIdx) const {
    ASSERT_THREAD(this);
    if (!inner_)
        return Vec3::Zero();

    Transform3D fallback;
    const Transform3D* t = nullptr;
    if (frame < transformCache_.size()) {
        t = &transformCache_[frame];
    } else {
        fallback = computeRawTransform(frame);
        t = &fallback;
    }
    const Vec3 raw = inner_->atomPosition(frame, atomIdx);
    return t->R * raw + t->T;
}

Mat3 TransformedConformation::displayRotation(std::size_t frame) const {
    ASSERT_THREAD(this);
    if (!inner_)
        return Mat3::Identity();
    if (frame < transformCache_.size())
        return transformCache_[frame].R;
    return computeRawTransform(frame).R;
}

std::shared_ptr<const QtConformationSnapshot>
TransformedConformation::loadSnapshot(std::size_t frame) {
    // Snapshots are full-fidelity per-frame source data — atom positions
    // PLUS calculator NPYs. We do NOT decorate snapshots: consumers that
    // need calculator data already read from the snapshot directly, and
    // applying our 3x3 transform to a snapshot's Pos column would diverge
    // from the inner conformation's atomPosition seam. Forward unchanged —
    // but TRIGGER the inner's (synchronous) load, do not merely read its
    // resident slot: the bare accessor returned a snapshot only when the inner
    // already held this exact frame, so NPY strips sampling frames the display
    // never loads saw a stale slot. requestSnapshot() makes inner->snapshot()
    // non-null for the requested frame.
    if (!inner_)
        return nullptr;
    inner_->requestSnapshot(frame);
    return inner_->snapshot(frame);
}

void TransformedConformation::setMode(Mode mode,
                                       std::size_t referenceFrame,
                                       std::vector<std::size_t> subsetAtoms) {
    ASSERT_THREAD(this);

    if (mode == Mode::ScientificAlignment) {
        qCWarning(cXform).noquote()
            << "scientific alignment must be installed with setScientificAlignment";
        return;
    }

    mode_ = mode;
    referenceFrame_ = referenceFrame;
    subsetAtoms_.clear();
    fitAtomIndices_.clear();

    if (inner_) {
        const std::size_t frames = inner_->frameCount();
        if (frames > 0 && referenceFrame_ >= frames) {
            qCWarning(cXform).noquote()
                << "reference frame" << static_cast<qlonglong>(referenceFrame_)
                << "out of range; clamping to" << static_cast<qlonglong>(frames - 1);
            referenceFrame_ = frames - 1;
        }
    }

    if (protein_) {
        const std::size_t atomCount = protein_->atomCount();
        if (mode_ == Mode::FitReference) {
            fitAtomIndices_.reserve(atomCount);
            for (std::size_t a = 0; a < atomCount; ++a)
                fitAtomIndices_.push_back(a);
        } else {
            subsetAtoms_.reserve(subsetAtoms.size());
            for (std::size_t a : subsetAtoms) {
                if (a >= atomCount) {
                    qCWarning(cXform).noquote()
                        << "FitSubset atom" << static_cast<qlonglong>(a)
                        << "out of range; dropping from subset";
                    continue;
                }
                subsetAtoms_.push_back(a);
            }
            fitAtomIndices_ = subsetAtoms_;
        }
    }

    // Bump the generation counter (semantic marker, primarily for logs /
    // future diagnostics) and wipe the transform sequence before rebuilding.
    ++generation_;
    transformCache_.clear();
    referencePositions_.clear();
    referenceCentroid_ = Vec3::Zero();
    rebuildReferenceMean();

    rebuildTransformCache();

    const char* modeName = nullptr;
    switch (mode_) {
        case Mode::FitReference: modeName = "fit_reference"; break;
        case Mode::FitSubset:    modeName = "fit_subset"; break;
        case Mode::ScientificAlignment: modeName = "scientific_alignment"; break;
    }
    qCInfo(cXform).noquote()
        << "mode set to" << modeName
        << "| ref_frame=" << static_cast<qlonglong>(referenceFrame_)
        << "| subset_size=" << static_cast<qlonglong>(subsetAtoms_.size())
        << "| smoothing_window=" << stabilisationWindow_
        << "| transforms=" << static_cast<qlonglong>(transformCache_.size())
        << "| generation=" << static_cast<qlonglong>(generation_);

    emit transformChanged();
}

bool TransformedConformation::setScientificAlignment(
    const ScientificAlignmentResult& alignment,
    QString* error) {
    ASSERT_THREAD(this);

    auto fail = [error](const QString& message) {
        if (error)
            *error = message;
        return false;
    };
    if (!inner_ || !protein_)
        return fail(QStringLiteral("no conformation is available"));
    if (!alignment.ok || !alignment.mean.converged)
        return fail(QStringLiteral("scientific alignment is not complete and converged"));
    if (alignment.policy.seedFrame >= inner_->frameCount()
        || !std::isfinite(alignment.policy.rotationTolerance)
        || !(alignment.policy.rotationTolerance > 0.0)) {
        return fail(QStringLiteral("scientific alignment policy is invalid"));
    }
    if (alignment.frameFits.size() != inner_->frameCount()) {
        return fail(QStringLiteral(
            "scientific transform count does not match the loaded frame count"));
    }
    if (alignment.fitAtomIndices.size() < 3
        || alignment.referencePositions.size() != alignment.fitAtomIndices.size()) {
        return fail(QStringLiteral("scientific fit selection or reference is incomplete"));
    }

    std::vector<bool> seen(protein_->atomCount(), false);
    for (std::size_t atom : alignment.fitAtomIndices) {
        if (atom >= protein_->atomCount() || seen[atom])
            return fail(QStringLiteral("scientific fit selection is invalid"));
        seen[atom] = true;
    }
    for (const Vec3& reference : alignment.referencePositions) {
        if (!reference.allFinite())
            return fail(QStringLiteral(
                "scientific reference contains nonfinite coordinates"));
    }

    std::vector<Transform3D> transforms;
    transforms.reserve(alignment.frameFits.size());
    for (const ScientificFrameFit& fit : alignment.frameFits) {
        const double orthogonalityError =
            (fit.rotation * fit.rotation.transpose() - Mat3::Identity()).norm();
        if (!fit.valid() || !fit.rotation.allFinite()
            || !fit.translation.allFinite()
            || std::abs(fit.rotation.determinant() - 1.0)
                > alignment.policy.rotationTolerance
            || orthogonalityError > alignment.policy.rotationTolerance) {
            return fail(QStringLiteral(
                "scientific transform sequence contains an invalid frame"));
        }
        transforms.push_back({fit.rotation, fit.translation});
    }

    mode_ = Mode::ScientificAlignment;
    referenceFrame_ = alignment.policy.seedFrame;
    subsetAtoms_ = alignment.fitAtomIndices;
    fitAtomIndices_ = alignment.fitAtomIndices;
    stabilisationWindow_ = 0;
    referencePositions_ = alignment.referencePositions;
    referenceCentroid_ = alignment.mean.centroidAnchor;
    transformCache_ = std::move(transforms);
    ++generation_;

    qCInfo(cXform).noquote()
        << "exact scientific alignment installed"
        << "| frames=" << static_cast<qulonglong>(transformCache_.size())
        << "| fit_atoms=" << static_cast<qulonglong>(fitAtomIndices_.size())
        << "| seed_frame=" << static_cast<qulonglong>(referenceFrame_)
        << "| smoothing_window=0"
        << "| generation=" << static_cast<qulonglong>(generation_);
    emit transformChanged();
    return true;
}

void TransformedConformation::setStabilisationWindow(int halfWidth) {
    ASSERT_THREAD(this);

    if (mode_ == Mode::ScientificAlignment) {
        if (halfWidth != 0) {
            qCWarning(cXform).noquote()
                << "smoothing is unavailable for exact scientific alignment";
        }
        return;
    }

    const int clamped = std::max(0, halfWidth);
    if (clamped == stabilisationWindow_)
        return;

    stabilisationWindow_ = clamped;
    ++generation_;
    rebuildTransformCache();

    qCInfo(cXform).noquote()
        << "stabilisation window set to" << stabilisationWindow_
        << "| transforms=" << static_cast<qlonglong>(transformCache_.size())
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

std::vector<Vec3> TransformedConformation::fitPositions(std::size_t frame) const {
    std::vector<Vec3> positions;
    if (!inner_ || !protein_)
        return positions;

    const std::size_t atomCount = protein_->atomCount();
    positions.reserve(fitAtomIndices_.size());
    for (std::size_t a : fitAtomIndices_) {
        if (a >= atomCount) {
            qCWarning(cXform).noquote()
                << "fit atom" << static_cast<qlonglong>(a)
                << "out of range at frame" << static_cast<qlonglong>(frame)
                << "; dropping from fit";
            continue;
        }
        positions.push_back(inner_->atomPosition(frame, a));
    }
    return positions;
}

void TransformedConformation::rebuildReferenceMean() {
    referencePositions_.clear();
    referenceCentroid_ = Vec3::Zero();

    if (!inner_ || !protein_)
        return;

    const std::size_t frames = inner_->frameCount();
    const std::size_t n = fitAtomIndices_.size();
    if (frames == 0 || n == 0)
        return;

    QElapsedTimer timer;
    timer.start();

    const std::vector<Vec3> seed = fitPositions(referenceFrame_);
    if (seed.size() != n) {
        qCWarning(cXform).noquote()
            << "iterative mean reference skipped; fit atom count mismatch"
            << "expected=" << static_cast<qlonglong>(n)
            << "actual=" << static_cast<qlonglong>(seed.size());
        return;
    }

    referenceCentroid_ = Centroid(seed);
    referencePositions_ = DemeanedCopy(seed, referenceCentroid_);

    int iterations = 0;
    double delta = 0.0;
    if (n >= 3) {
        for (int iter = 0; iter < kMeanReferenceMaxIterations; ++iter) {
            std::vector<Vec3> accum(n, Vec3::Zero());
            std::size_t usedFrames = 0;

            for (std::size_t frame = 0; frame < frames; ++frame) {
                const std::vector<Vec3> current = fitPositions(frame);
                if (current.size() != n) {
                    qCWarning(cXform).noquote()
                        << "iterative mean dropped frame" << static_cast<qlonglong>(frame)
                        << "because fit atom count changed"
                        << "expected=" << static_cast<qlonglong>(n)
                        << "actual=" << static_cast<qlonglong>(current.size());
                    continue;
                }

                const Vec3 cc = Centroid(current);
                const std::vector<Vec3> currentDemeaned = DemeanedCopy(current, cc);
                // v2: pass robust/weighted atom weights into KabschFit here.
                const Transform3D fit = KabschFit(currentDemeaned, referencePositions_);
                for (std::size_t i = 0; i < n; ++i)
                    accum[i] += fit.R * currentDemeaned[i];
                ++usedFrames;
            }

            if (usedFrames == 0) {
                qCWarning(cXform).noquote()
                    << "iterative mean reference failed; no frames contributed";
                referencePositions_.clear();
                referenceCentroid_ = Vec3::Zero();
                return;
            }

            for (Vec3& p : accum)
                p /= static_cast<double>(usedFrames);
            const Vec3 avgCentroid = Centroid(accum);
            std::vector<Vec3> avg = DemeanedCopy(accum, avgCentroid);

            delta = RmsDifference(avg, referencePositions_);
            referencePositions_ = std::move(avg);
            iterations = iter + 1;
            if (delta < kMeanReferenceEpsAngstrom)
                break;
        }
    }

    if (n >= 3 && iterations >= kMeanReferenceMaxIterations
            && delta >= kMeanReferenceEpsAngstrom) {
        qCWarning(cXform).noquote()
            << "iterative mean reference did NOT converge | delta_A=" << delta
            << "| eps_A=" << kMeanReferenceEpsAngstrom
            << "| iterations=" << iterations
            << "— shipping last average (rigid + centroid-pinned, still valid)";
    }

    qCInfo(cXform).noquote()
        << "iterative mean reference built"
        << "| fit_atoms=" << static_cast<qlonglong>(n)
        << "| frames=" << static_cast<qlonglong>(frames)
        << "| iterations=" << iterations
        << "| delta_A=" << delta
        << "| eps_A=" << kMeanReferenceEpsAngstrom
        << "| max_iter=" << kMeanReferenceMaxIterations
        << "| anchor_frame=" << static_cast<qlonglong>(referenceFrame_)
        << "| ms=" << timer.elapsed();
}

TransformedConformation::FrameFit
TransformedConformation::computeRawFrameFit(std::size_t frame) const {
    FrameFit out;  // identity by default

    if (!inner_ || referencePositions_.empty())
        return out;

    const std::vector<Vec3> current = fitPositions(frame);
    const std::size_t n = referencePositions_.size();
    if (current.size() != n || n == 0)
        return out;

    out.currentCentroid = Centroid(current);
    if (n >= 3)
        out.transform.R = KabschFit(current, referencePositions_).R;
    out.transform.T = referenceCentroid_ - out.transform.R * out.currentCentroid;
    return out;
}

TransformedConformation::Transform3D
TransformedConformation::computeRawTransform(std::size_t frame) const {
    return computeRawFrameFit(frame).transform;
}

void TransformedConformation::rebuildTransformCache() {
    transformCache_.clear();

    if (!inner_)
        return;
    const std::size_t frames = inner_->frameCount();
    if (frames == 0)
        return;

    QElapsedTimer timer;
    timer.start();
    const std::vector<FrameFit> raw = computeRawTransformSequence();
    transformCache_ = smoothTransformSequence(raw);
    qCInfo(cXform).noquote()
        << "transform cache rebuilt | frames=" << static_cast<qlonglong>(frames)
        << "| window=" << stabilisationWindow_
        << "| ms=" << timer.elapsed();
}

std::vector<TransformedConformation::FrameFit>
TransformedConformation::computeRawTransformSequence() const {
    std::vector<FrameFit> raw;
    if (!inner_)
        return raw;

    const std::size_t frames = inner_->frameCount();
    raw.reserve(frames);
    for (std::size_t frame = 0; frame < frames; ++frame)
        raw.push_back(computeRawFrameFit(frame));
    return raw;
}

std::vector<TransformedConformation::Transform3D>
TransformedConformation::smoothTransformSequence(const std::vector<FrameFit>& raw) const {
    if (raw.empty())
        return {};

    constexpr double kMinQuaternionNorm = 1e-12;
    const std::size_t frames = raw.size();
    std::vector<Transform3D> smoothed(frames);

    if (stabilisationWindow_ <= 0) {
        for (std::size_t frame = 0; frame < frames; ++frame) {
            Transform3D out = raw[frame].transform;
            out.T = referenceCentroid_ - out.R * raw[frame].currentCentroid;
            smoothed[frame] = out;
        }
        return smoothed;
    }

    const std::size_t halfWidth = static_cast<std::size_t>(stabilisationWindow_);

    for (std::size_t frame = 0; frame < frames; ++frame) {
        const std::size_t begin = (frame > halfWidth) ? frame - halfWidth : 0;
        const std::size_t end = std::min(frames - 1, frame + halfWidth);

        const Eigen::Quaterniond qRef = NormalisedQuaternion(raw[frame].transform.R);
        Eigen::Vector4d qSum = Eigen::Vector4d::Zero();

        for (std::size_t k = begin; k <= end; ++k) {
            Eigen::Quaterniond q = NormalisedQuaternion(raw[k].transform.R);
            if (q.dot(qRef) < 0.0)
                q.coeffs() *= -1.0;
            qSum += q.coeffs();
        }

        Transform3D out;
        if (qSum.norm() < kMinQuaternionNorm) {
            out.R = raw[frame].transform.R;
        } else {
            Eigen::Quaterniond qMean;
            qMean.coeffs() = qSum;
            qMean.normalize();
            out.R = qMean.toRotationMatrix();
        }
        out.T = referenceCentroid_ - out.R * raw[frame].currentCentroid;
        smoothed[frame] = out;
    }
    return smoothed;
}

// Kabsch fit — delegates to h5reader::math::ComputeSubsetTransform in
// FitTargetMath.h. The free function is the canonical implementation
// (Codex finding #6) and owns the degeneracy policy (rank-degenerate
// → std::nullopt per Codex finding #4). Both the camera path
// (CameraComposer::writeSubset) and the data path (this) now share one
// failure semantics: if the fit is degenerate, freeze on identity
// rotation with translation-only centroid alignment. This kills the
// divergent failure modes the prior duplicate implementation had —
// camera path nullopt, data path silent-identity-with-bad-T.
//
// Freeze policy when degenerate: R = identity, T = (cr - cc). The atom
// positions become "centred on the reference centroid but otherwise
// unrotated" for that frame; visually this means the molecule continues
// to display from its current orientation without a sudden re-orient
// from a numerically-arbitrary SVD null-space basis. The next frame
// reattempts the fit; if conditioning improves, the rotation
// re-engages smoothly.
TransformedConformation::Transform3D
TransformedConformation::KabschFit(const std::vector<Vec3>& current,
                                    const std::vector<Vec3>& reference) {
    Transform3D out;  // identity by default
    auto transform = h5reader::math::ComputeSubsetTransform(current, reference);
    if (!transform) {
        // Degenerate input (n < 3, rank-deficient, or det validation
        // failed). Freeze rotation; if we have at least one matched pair,
        // still align centroids so the camera/positions land on the
        // expected anchor and the next frame can try again without a
        // visible jump.
        if (current.size() >= 1 && current.size() == reference.size()) {
            Vec3 cc = Vec3::Zero();
            Vec3 cr = Vec3::Zero();
            for (std::size_t i = 0; i < current.size(); ++i) {
                cc += current[i];
                cr += reference[i];
            }
            cc /= static_cast<double>(current.size());
            cr /= static_cast<double>(current.size());
            out.T = cr - cc;  // R already identity
        }
        return out;
    }
    out.R = transform->R;
    out.T = transform->T;
    return out;
}

}  // namespace h5reader::model

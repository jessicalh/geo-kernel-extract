#include "ScientificAlignment.h"

#include <Eigen/SVD>

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace h5reader::model {

namespace {

bool IsFinite(const Vec3& value) {
    return value.allFinite();
}

bool AllFinite(const std::vector<Vec3>& values) {
    return std::all_of(values.begin(), values.end(), IsFinite);
}

Vec3 Centroid(const std::vector<Vec3>& values) {
    Vec3 centroid = Vec3::Zero();
    for (const Vec3& value : values)
        centroid += value;
    if (!values.empty())
        centroid /= static_cast<double>(values.size());
    return centroid;
}

double RmsDifference(const std::vector<Vec3>& left,
                     const std::vector<Vec3>& right) {
    if (left.empty() || left.size() != right.size())
        return std::numeric_limits<double>::infinity();

    double squaredSum = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i)
        squaredSum += (left[i] - right[i]).squaredNorm();
    return std::sqrt(squaredSum / static_cast<double>(left.size()));
}

std::vector<Vec3> SelectFrame(const ScientificPositionTable& positions,
                              std::size_t frame,
                              const std::vector<std::size_t>& atomIndices) {
    std::vector<Vec3> selected;
    selected.reserve(atomIndices.size());
    for (std::size_t atom : atomIndices)
        selected.push_back(positions.at(frame, atom));
    return selected;
}

QString FitFailureMessage(const char* stage,
                          std::size_t frame,
                          ScientificFitStatus status) {
    return QStringLiteral("%1 fit failed at frame %2: %3")
        .arg(QString::fromLatin1(stage))
        .arg(static_cast<qulonglong>(frame))
        .arg(QString::fromLatin1(NameForScientificFitStatus(status)));
}

}  // namespace

const char* NameForScientificFitStatus(ScientificFitStatus status) {
    switch (status) {
        case ScientificFitStatus::Valid: return "valid";
        case ScientificFitStatus::TooFewPoints: return "too_few_points";
        case ScientificFitStatus::InputSizeMismatch: return "input_size_mismatch";
        case ScientificFitStatus::NonFiniteInput: return "nonfinite_input";
        case ScientificFitStatus::RankBelowTwo: return "rank_below_two";
        case ScientificFitStatus::InvalidRotation: return "invalid_rotation";
        case ScientificFitStatus::NonFiniteOutput: return "nonfinite_output";
    }
    return "unknown";
}

bool ScientificPositionTable::hasExpectedSize() const {
    if (frameCount == 0 || atomCount == 0)
        return false;
    if (frameCount > std::numeric_limits<std::size_t>::max() / atomCount)
        return false;
    return values.size() == frameCount * atomCount;
}

const Vec3& ScientificPositionTable::at(std::size_t frame,
                                        std::size_t atom) const {
    return values[frame * atomCount + atom];
}

ScientificFrameFit ComputeScientificKabsch(
    const std::vector<Vec3>& current,
    const std::vector<Vec3>& reference,
    const ScientificAlignmentPolicy& policy) {
    ScientificFrameFit result;
    if (current.size() != reference.size()) {
        result.status = ScientificFitStatus::InputSizeMismatch;
        return result;
    }
    if (current.size() < 3) {
        result.status = ScientificFitStatus::TooFewPoints;
        return result;
    }
    if (!AllFinite(current) || !AllFinite(reference)) {
        result.status = ScientificFitStatus::NonFiniteInput;
        return result;
    }

    const Vec3 currentCentroid = Centroid(current);
    const Vec3 referenceCentroid = Centroid(reference);
    Eigen::MatrixXd currentCentered(3, static_cast<Eigen::Index>(current.size()));
    Eigen::MatrixXd referenceCentered(3, static_cast<Eigen::Index>(reference.size()));
    for (std::size_t i = 0; i < current.size(); ++i) {
        currentCentered.col(static_cast<Eigen::Index>(i)) = current[i] - currentCentroid;
        referenceCentered.col(static_cast<Eigen::Index>(i)) = reference[i] - referenceCentroid;
    }

    const Mat3 covariance = currentCentered * referenceCentered.transpose();
    Eigen::JacobiSVD<Mat3> svd(covariance, Eigen::ComputeFullU | Eigen::ComputeFullV);
    const Mat3& u = svd.matrixU();
    const Mat3& v = svd.matrixV();
    result.singularValues = svd.singularValues();
    if (!result.singularValues.allFinite()) {
        result.status = ScientificFitStatus::NonFiniteOutput;
        return result;
    }

    const double largest = result.singularValues(0);
    result.rankThreshold = std::max(policy.rankAbsoluteTolerance,
                                    policy.rankRelativeTolerance * largest);
    result.numericalRank = 0;
    for (int i = 0; i < 3; ++i) {
        if (result.singularValues(i) > result.rankThreshold)
            ++result.numericalRank;
    }
    if (result.numericalRank < 2) {
        result.status = ScientificFitStatus::RankBelowTwo;
        return result;
    }

    Mat3 determinantCorrection = Mat3::Identity();
    if ((v * u.transpose()).determinant() < 0.0)
        determinantCorrection(2, 2) = -1.0;
    result.rotation = v * determinantCorrection * u.transpose();
    result.translation = referenceCentroid - result.rotation * currentCentroid;

    const double determinant = result.rotation.determinant();
    const double orthogonalityError =
        (result.rotation * result.rotation.transpose() - Mat3::Identity()).norm();
    if (!result.rotation.allFinite()
        || std::abs(determinant - 1.0) > policy.rotationTolerance
        || orthogonalityError > policy.rotationTolerance) {
        result.status = ScientificFitStatus::InvalidRotation;
        return result;
    }
    if (!result.translation.allFinite()) {
        result.status = ScientificFitStatus::NonFiniteOutput;
        return result;
    }

    double squaredResidual = 0.0;
    for (std::size_t i = 0; i < current.size(); ++i) {
        const Vec3 residual = result.rotation * current[i]
            + result.translation - reference[i];
        squaredResidual += residual.squaredNorm();
    }
    result.rmsdAngstrom =
        std::sqrt(squaredResidual / static_cast<double>(current.size()));
    if (!std::isfinite(result.rmsdAngstrom)) {
        result.status = ScientificFitStatus::NonFiniteOutput;
        return result;
    }

    result.status = ScientificFitStatus::Valid;
    return result;
}

ScientificAlignmentResult BuildScientificAlignment(
    const ScientificPositionTable& positions,
    const std::vector<std::size_t>& fitAtomIndices,
    const ScientificAlignmentPolicy& policy) {
    ScientificAlignmentResult result;
    result.policy = policy;
    result.fitAtomIndices = fitAtomIndices;

    if (!positions.hasExpectedSize()) {
        result.error = QStringLiteral("position table dimensions are inconsistent");
        return result;
    }
    if (policy.seedFrame >= positions.frameCount) {
        result.error = QStringLiteral("seed frame is outside the position table");
        return result;
    }
    if (policy.maxIterations <= 0
        || !std::isfinite(policy.convergenceToleranceAngstrom)
        || !(policy.convergenceToleranceAngstrom > 0.0)
        || !std::isfinite(policy.rankRelativeTolerance)
        || !(policy.rankRelativeTolerance >= 0.0)
        || !std::isfinite(policy.rankAbsoluteTolerance)
        || !(policy.rankAbsoluteTolerance >= 0.0)
        || !std::isfinite(policy.rotationTolerance)
        || !(policy.rotationTolerance > 0.0)) {
        result.error = QStringLiteral("scientific alignment policy is invalid");
        return result;
    }
    if (fitAtomIndices.size() < 3) {
        result.error = QStringLiteral("scientific alignment needs at least three fit atoms");
        return result;
    }

    std::vector<bool> seen(positions.atomCount, false);
    for (std::size_t atom : fitAtomIndices) {
        if (atom >= positions.atomCount) {
            result.error = QStringLiteral("fit atom index %1 is outside the atom axis")
                .arg(static_cast<qulonglong>(atom));
            return result;
        }
        if (seen[atom]) {
            result.error = QStringLiteral("fit atom index %1 is duplicated")
                .arg(static_cast<qulonglong>(atom));
            return result;
        }
        seen[atom] = true;
    }

    std::vector<Vec3> reference =
        SelectFrame(positions, policy.seedFrame, fitAtomIndices);
    if (!AllFinite(reference)) {
        result.error = QStringLiteral("seed reference contains a nonfinite position");
        return result;
    }
    const Vec3 centroidAnchor = Centroid(reference);
    result.mean.centroidAnchor = centroidAnchor;

    double finalDelta = std::numeric_limits<double>::infinity();
    for (int iteration = 0; iteration < policy.maxIterations; ++iteration) {
        std::vector<Vec3> accumulated(reference.size(), Vec3::Zero());

        for (std::size_t frame = 0; frame < positions.frameCount; ++frame) {
            const std::vector<Vec3> current =
                SelectFrame(positions, frame, fitAtomIndices);
            const ScientificFrameFit fit =
                ComputeScientificKabsch(current, reference, policy);
            if (!fit.valid()) {
                ++result.mean.referenceBuildDegeneracyCount;
                result.error = FitFailureMessage("reference-build", frame, fit.status);
                return result;
            }
            for (std::size_t i = 0; i < current.size(); ++i)
                accumulated[i] += fit.rotation * current[i] + fit.translation;
        }

        const double inverseFrameCount = 1.0 / static_cast<double>(positions.frameCount);
        for (Vec3& value : accumulated)
            value *= inverseFrameCount;

        // Preserve the seed frame's centroid as the explicit world anchor.
        const Vec3 accumulatedCentroid = Centroid(accumulated);
        for (Vec3& value : accumulated)
            value += centroidAnchor - accumulatedCentroid;

        finalDelta = RmsDifference(accumulated, reference);
        if (!std::isfinite(finalDelta)) {
            result.error = QStringLiteral("iterative mean produced a nonfinite delta");
            return result;
        }
        reference = std::move(accumulated);
        result.mean.iterations = iteration + 1;
        result.mean.finalDeltaAngstrom = finalDelta;
        if (finalDelta <= policy.convergenceToleranceAngstrom) {
            result.mean.converged = true;
            break;
        }
    }

    if (!result.mean.converged) {
        result.error = QStringLiteral(
            "iterative mean did not converge after %1 iterations (delta %2 A, tolerance %3 A)")
            .arg(policy.maxIterations)
            .arg(finalDelta, 0, 'g', 17)
            .arg(policy.convergenceToleranceAngstrom, 0, 'g', 17);
        return result;
    }

    result.referencePositions = reference;
    result.frameFits.reserve(positions.frameCount);
    for (std::size_t frame = 0; frame < positions.frameCount; ++frame) {
        const std::vector<Vec3> current =
            SelectFrame(positions, frame, fitAtomIndices);
        ScientificFrameFit fit = ComputeScientificKabsch(current, reference, policy);
        if (!fit.valid()) {
            result.error = FitFailureMessage("final", frame, fit.status);
            result.frameFits.push_back(std::move(fit));
            return result;
        }
        result.frameFits.push_back(std::move(fit));
    }

    result.ok = true;
    return result;
}

}  // namespace h5reader::model

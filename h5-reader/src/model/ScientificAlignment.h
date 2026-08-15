#pragma once

#include "Types.h"

#include <QString>

#include <cstddef>
#include <cstdint>
#include <vector>

namespace h5reader::model {

enum class ScientificFitStatus : std::uint8_t {
    Valid = 0,
    TooFewPoints = 1,
    InputSizeMismatch = 2,
    NonFiniteInput = 3,
    RankBelowTwo = 4,
    InvalidRotation = 5,
    NonFiniteOutput = 6,
};

const char* NameForScientificFitStatus(ScientificFitStatus status);

struct ScientificAlignmentPolicy {
    std::size_t seedFrame = 0;
    int maxIterations = 20;
    double convergenceToleranceAngstrom = 1e-4;

    // Numerical rank is the number of singular values greater than
    // max(rankAbsoluteTolerance, rankRelativeTolerance * sigma_max).
    // A scientific rigid fit needs rank >= 2: non-collinear planar input is
    // sufficient, while a line or a point is not.
    double rankRelativeTolerance = 1e-12;
    double rankAbsoluteTolerance = 1e-12;
    double rotationTolerance = 1e-8;
};

// Dense raw positions in frame-major order: values[(frame * atoms) + atom].
// Keeping this table independent of Conformation makes the numerical work
// directly testable and makes the H5 read boundary explicit in the exporter.
struct ScientificPositionTable {
    std::size_t frameCount = 0;
    std::size_t atomCount = 0;
    std::vector<Vec3> values;

    bool hasExpectedSize() const;
    const Vec3& at(std::size_t frame, std::size_t atom) const;
};

struct ScientificFrameFit {
    ScientificFitStatus status = ScientificFitStatus::TooFewPoints;
    Mat3 rotation = Mat3::Identity();
    Vec3 translation = Vec3::Zero();
    Vec3 singularValues = Vec3::Zero();
    int numericalRank = 0;
    double rankThreshold = 0.0;
    double rmsdAngstrom = 0.0;

    bool valid() const { return status == ScientificFitStatus::Valid; }
};

struct ScientificMeanDiagnostics {
    bool converged = false;
    int iterations = 0;
    double finalDeltaAngstrom = 0.0;
    std::size_t referenceBuildDegeneracyCount = 0;
    Vec3 centroidAnchor = Vec3::Zero();
};

struct ScientificAlignmentResult {
    bool ok = false;
    QString error;
    ScientificAlignmentPolicy policy;
    std::vector<std::size_t> fitAtomIndices;
    std::vector<Vec3> referencePositions;
    std::vector<ScientificFrameFit> frameFits;
    ScientificMeanDiagnostics mean;
};

ScientificFrameFit ComputeScientificKabsch(
    const std::vector<Vec3>& current,
    const std::vector<Vec3>& reference,
    const ScientificAlignmentPolicy& policy = {});

ScientificAlignmentResult BuildScientificAlignment(
    const ScientificPositionTable& positions,
    const std::vector<std::size_t>& fitAtomIndices,
    const ScientificAlignmentPolicy& policy = {});

}  // namespace h5reader::model

#pragma once

#include "ConformationGeometry.h"
#include "DftShielding.h"

#include "../io/CalcsetManifest.h"

#include <QString>

#include <cstddef>
#include <optional>
#include <vector>

namespace h5reader::model {

class Conformation;
class QtProtein;

struct RingCurrentPathAnalysisOptions {
    std::optional<std::size_t> atom;
    std::optional<std::size_t> ring;
    std::optional<int> startFrame;
    std::optional<int> endFrame;
    double surfaceToleranceA = 1e-9;
    double predictorZeroToleranceA3 = 1e-12;
    bool includeSaturatedRings = false;
    int minSamples = 6;
    int minSamplesPerLobe = 3;
    double minPredictorSpanA3 = 0.02;
    double minAbsLobePredictorA3 = 0.005;
    int maxPaths = 25;
    int nullShiftCount = 64;
};

struct RingCurrentPathSample {
    int frameIndex = -1;
    double timePs = 0.0;
    RingNullMeasurement geometry;
    double pointDipoleGeometryA3 = 0.0;
    double inverseDistanceCubedA3 = 0.0;
    double angularFactor = 0.0;
    SphericalTensor biotSavart;
    DftAtomShielding orca;
};

struct RingCurrentLinearFit {
    bool valid = false;
    int sampleCount = 0;
    double intercept = 0.0;
    double scale = 0.0;
    double r2 = 0.0;
    double correlation = 0.0;
    double sse = 0.0;
    double sst = 0.0;
    int nullShiftCount = 0;
    double nullMedianR2 = 0.0;
    double nullMaxR2 = 0.0;
    double nullGeRealFraction = 0.0;
    std::vector<double> nullR2;
};

struct RingCurrentPathResult {
    std::size_t atom = 0;
    std::size_t ring = 0;
    std::vector<RingCurrentPathSample> samples;

    double minPointDipoleGeometryA3 = 0.0;
    double maxPointDipoleGeometryA3 = 0.0;
    double pointDipoleGeometrySpanA3 = 0.0;
    int positivePredictorSamples = 0;
    int negativePredictorSamples = 0;
    int predictorSignChanges = 0;
    bool hardLobeCrossing = false;

    RingCurrentLinearFit orcaTotalT0VsPointDipole;
    RingCurrentLinearFit orcaDiamagneticT0VsPointDipole;
    RingCurrentLinearFit orcaParamagneticT0VsPointDipole;
    RingCurrentLinearFit orcaTotalT2MagnitudeVsPointDipole;
    RingCurrentLinearFit orcaTotalT0VsBiotSavartT0;
    RingCurrentLinearFit biotSavartT0VsPointDipole;
    RingCurrentLinearFit orcaTotalT0VsDistanceOnly;
    RingCurrentLinearFit orcaTotalT0VsAngleOnly;
};

struct RingCurrentPathAnalysisSummary {
    bool complete = false;
    int dftFramesDeclared = 0;
    int dftFramesLoaded = 0;
    int dftFramesSkipped = 0;
    int atomsScanned = 0;
    int ringsScanned = 0;
    int pathsConsidered = 0;
    int pathsRejectedForSamples = 0;
    int pathsRejectedForHardCrossing = 0;
    int pathsRejectedForWeakLobes = 0;
    int pathCount = 0;
    bool truncatedByMaxPaths = false;
};

// Analyse atom/ring paths that sample both signs of the point-dipole geometry,
// then compare that geometry and a finite-ring Biot-Savart calculation with the
// observed ORCA shielding trace.
class RingCurrentPathAnalysis {
public:
    explicit RingCurrentPathAnalysis(const RingCurrentPathAnalysisOptions& options = {});

    // q = angularFactor / distance_A^3; both return inverse cubic angstroms.
    static double pointDipoleGeometryA3(const RingNullMeasurement& geometry);
    static double inverseDistanceCubedA3(const RingNullMeasurement& geometry);

    bool run(const QtProtein& protein,
             const Conformation& conformation,
             const std::vector<h5reader::io::DftFrame>& dftFrames,
             QString* error = nullptr);

    const RingCurrentPathAnalysisOptions& options() const { return options_; }
    const RingCurrentPathAnalysisSummary& summary() const { return summary_; }
    const std::vector<RingCurrentPathResult>& paths() const { return paths_; }

private:
    RingCurrentPathAnalysisOptions options_;
    RingCurrentPathAnalysisSummary summary_;
    std::vector<RingCurrentPathResult> paths_;
};

}  // namespace h5reader::model

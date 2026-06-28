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

struct RingCurrentFaceCollarOptions {
    std::optional<std::size_t> atom;
    std::optional<std::size_t> ring;
    std::optional<int> startFrame;
    std::optional<int> endFrame;
    double surfaceToleranceA = 1e-9;
    double templateZeroTolerance = 1e-12;
    bool includeSaturatedRings = false;
    int minSamples = 6;
    int minSamplesPerLobe = 3;
    double minExpectedRelationshipSpan = 0.02;
    double minAbsLobeExpectedValue = 0.005;
    int maxEntries = 25;
    int nullShiftCount = 64;
};

struct RingCurrentFaceSample {
    int frameIndex = -1;
    double timePs = 0.0;
    RingNullMeasurement geometry;
    double expectedRelationshipValue = 0.0;
    double distanceOnlyValue = 0.0;
    double angularOnlyValue = 0.0;
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

struct RingCurrentFaceEntry {
    std::size_t atom = 0;
    std::size_t ring = 0;
    std::vector<RingCurrentFaceSample> samples;

    double minTemplate = 0.0;
    double maxTemplate = 0.0;
    double templateSpan = 0.0;
    int positiveTemplateSamples = 0;
    int negativeTemplateSamples = 0;
    int templateSignChanges = 0;
    bool hardLobeCrossing = false;

    RingCurrentLinearFit orcaTotalT0;
    RingCurrentLinearFit orcaDiamagneticT0;
    RingCurrentLinearFit orcaParamagneticT0;
    RingCurrentLinearFit orcaTotalT2Magnitude;
    RingCurrentLinearFit biotSavartOrcaTotalT0;
    RingCurrentLinearFit expectedRelationshipBiotSavartT0;
    RingCurrentLinearFit distanceOnlyOrcaTotalT0;
    RingCurrentLinearFit angularOnlyOrcaTotalT0;
};

struct RingCurrentFaceCollarSummary {
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
    int entryCount = 0;
    bool truncatedByMaxEntries = false;
};

// A simple receiver: detect atom/ring paths whose ring-current template crosses
// both lobes, stash the ORCA samples, then fit the template to the observed QM
// trace. Collection and evaluation stay explicit so the stash can be audited.
class RingCurrentFaceCollar {
public:
    explicit RingCurrentFaceCollar(RingCurrentFaceCollarOptions options = {});

    bool collect(const QtProtein& protein,
                 Conformation& conformation,
                 const std::vector<h5reader::io::DftFrame>& dftFrames,
                 QString* error = nullptr);

    const RingCurrentFaceCollarOptions& options() const { return options_; }
    const RingCurrentFaceCollarSummary& summary() const { return summary_; }
    const std::vector<RingCurrentFaceEntry>& entries() const { return entries_; }

private:
    RingCurrentFaceCollarOptions options_;
    RingCurrentFaceCollarSummary summary_;
    std::vector<RingCurrentFaceEntry> entries_;
};

}  // namespace h5reader::model

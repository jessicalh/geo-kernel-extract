#include "RingCurrentPathAnalysis.h"

#include "Conformation.h"
#include "QtProtein.h"

#include "../calculators/QtBiotSavartCalc.h"
#include "../io/DftShieldingLoader.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <utility>

namespace h5reader::model {
namespace {

struct LoadedDftFrame {
    std::shared_ptr<const DftShieldingFrame> shielding;
};

struct DftCandidate {
    const h5reader::io::DftFrame* declared = nullptr;
    std::size_t frameRow = 0;
};

std::vector<std::size_t> atomScanList(const QtProtein& protein,
                                      const RingCurrentPathAnalysisOptions& options) {
    std::vector<std::size_t> out;
    if (options.atom) {
        out.push_back(*options.atom);
        return out;
    }
    out.reserve(protein.atomCount());
    for (std::size_t i = 0; i < protein.atomCount(); ++i)
        out.push_back(i);
    return out;
}

std::vector<std::size_t> ringScanList(const QtProtein& protein,
                                      const RingCurrentPathAnalysisOptions& options) {
    std::vector<std::size_t> out;
    if (options.ring) {
        out.push_back(*options.ring);
        return out;
    }
    out.reserve(protein.ringCount());
    for (std::size_t i = 0; i < protein.ringCount(); ++i) {
        if (!options.includeSaturatedRings && !protein.ring(i).IsAromatic())
            continue;
        out.push_back(i);
    }
    return out;
}

std::optional<LoadedDftFrame> loadDftFrame(const h5reader::io::DftFrame& declared,
                                           const QtProtein& protein) {
    QString metaError;
    if (!declared.LoadMeta(&metaError))
        return std::nullopt;

    std::shared_ptr<const DftShieldingFrame> frame =
        h5reader::io::DftShieldingLoader::LoadAndValidate(declared.meta_json_abspath,
                                                          &protein);
    if (!frame || !frame->valid || frame->atoms.size() < protein.atomCount())
        return std::nullopt;

    return LoadedDftFrame{std::move(frame)};
}

int signOf(double value, double tolerance) {
    if (!std::isfinite(value) || std::abs(value) <= tolerance)
        return 0;
    return value > 0.0 ? 1 : -1;
}

double t2Magnitude(const SphericalTensor& tensor) {
    return tensor.T2Magnitude();
}

using SampleSignal = double (*)(const RingCurrentPathSample&);

struct FitSignals {
    SampleSignal expected = nullptr;
    SampleSignal observed = nullptr;
};

struct RegressionResult {
    double intercept = 0.0;
    double scale = 0.0;
    double r2 = 0.0;
    double correlation = 0.0;
    double sse = 0.0;
    double sst = 0.0;
};

RingCurrentLinearFit fitLinear(const std::vector<RingCurrentPathSample>& samples,
                               FitSignals channels,
                               int requestedNullShifts) {
    RingCurrentLinearFit out;
    out.sampleCount = static_cast<int>(samples.size());
    const std::size_t n = samples.size();
    if (n < 3 || !channels.expected || !channels.observed)
        return out;

    std::vector<double> x;
    std::vector<double> y;
    x.reserve(n);
    y.reserve(n);
    for (const RingCurrentPathSample& sample : samples) {
        const double xi = channels.expected(sample);
        const double yi = channels.observed(sample);
        if (!std::isfinite(xi) || !std::isfinite(yi))
            return out;
        x.push_back(xi);
        y.push_back(yi);
    }

    const auto regress = [](const std::vector<double>& xs,
                            const std::vector<double>& ys)
                            -> std::optional<RegressionResult> {
        const std::size_t count = xs.size();
        if (count < 3 || count != ys.size())
            return std::nullopt;

        const double meanX = std::accumulate(xs.begin(), xs.end(), 0.0) /
                             static_cast<double>(count);
        const double meanY = std::accumulate(ys.begin(), ys.end(), 0.0) /
                             static_cast<double>(count);
        double sxx = 0.0;
        double syy = 0.0;
        double sxy = 0.0;
        for (std::size_t i = 0; i < count; ++i) {
            const double dx = xs[i] - meanX;
            const double dy = ys[i] - meanY;
            sxx += dx * dx;
            syy += dy * dy;
            sxy += dx * dy;
        }
        if (sxx <= 1e-30 || syy <= 1e-30)
            return std::nullopt;

        RegressionResult result;
        result.scale = sxy / sxx;
        result.intercept = meanY - result.scale * meanX;
        for (std::size_t i = 0; i < count; ++i) {
            const double residual = ys[i] - (result.intercept + result.scale * xs[i]);
            result.sse += residual * residual;
        }
        result.sst = syy;
        result.r2 = 1.0 - (result.sse / result.sst);
        result.correlation = sxy / std::sqrt(sxx * syy);
        if (!std::isfinite(result.r2) || !std::isfinite(result.correlation))
            return std::nullopt;
        return result;
    };

    const std::optional<RegressionResult> fit = regress(x, y);
    if (!fit) {
        return out;
    }
    out.intercept = fit->intercept;
    out.scale = fit->scale;
    out.r2 = fit->r2;
    out.correlation = fit->correlation;
    out.sse = fit->sse;
    out.sst = fit->sst;
    out.valid = true;

    const int maxShifts = std::max(0, std::min(requestedNullShifts, static_cast<int>(n) - 1));
    out.nullR2.reserve(static_cast<std::size_t>(maxShifts));
    for (int shift = 1; shift <= maxShifts; ++shift) {
        std::vector<double> shiftedY;
        shiftedY.reserve(n);
        for (std::size_t i = 0; i < n; ++i)
            shiftedY.push_back(y[(i + static_cast<std::size_t>(shift)) % n]);

        if (const std::optional<RegressionResult> shiftedFit = regress(x, shiftedY))
            out.nullR2.push_back(shiftedFit->r2);
    }

    out.nullShiftCount = static_cast<int>(out.nullR2.size());
    if (!out.nullR2.empty()) {
        std::vector<double> sorted = out.nullR2;
        std::sort(sorted.begin(), sorted.end());
        out.nullMedianR2 = sorted[sorted.size() / 2];
        out.nullMaxR2 = sorted.back();
        int geReal = 0;
        for (double r2 : out.nullR2) {
            if (r2 >= out.r2)
                ++geReal;
        }
        out.nullGeRealFraction =
            static_cast<double>(geReal + 1) / static_cast<double>(out.nullR2.size() + 1);
    }

    return out;
}

double orcaTotalT0(const RingCurrentPathSample& sample) {
    return sample.orca.total.T0;
}

double pointDipoleValue(const RingCurrentPathSample& sample) {
    return sample.pointDipoleGeometryA3;
}

double distanceValue(const RingCurrentPathSample& sample) {
    return sample.inverseDistanceCubedA3;
}

double angularValue(const RingCurrentPathSample& sample) {
    return sample.angularFactor;
}

double biotSavartT0(const RingCurrentPathSample& sample) {
    return sample.biotSavart.T0;
}

double orcaDiaT0(const RingCurrentPathSample& sample) {
    return sample.orca.dia.T0;
}

double orcaParaT0(const RingCurrentPathSample& sample) {
    return sample.orca.para.T0;
}

double orcaTotalT2Magnitude(const RingCurrentPathSample& sample) {
    return t2Magnitude(sample.orca.total);
}

}  // namespace

RingCurrentPathAnalysis::RingCurrentPathAnalysis(const RingCurrentPathAnalysisOptions& options)
    : options_(options) {}

double RingCurrentPathAnalysis::pointDipoleGeometryA3(
    const RingNullMeasurement& geometry) {
    if (!geometry.valid || geometry.distanceA <= 1e-9)
        return 0.0;
    return geometry.angularFactor /
           (geometry.distanceA * geometry.distanceA * geometry.distanceA);
}

double RingCurrentPathAnalysis::inverseDistanceCubedA3(
    const RingNullMeasurement& geometry) {
    if (!geometry.valid || geometry.distanceA <= 1e-9)
        return 0.0;
    return 1.0 /
           (geometry.distanceA * geometry.distanceA * geometry.distanceA);
}

bool RingCurrentPathAnalysis::run(
    const QtProtein& protein,
    const Conformation& conformation,
    const std::vector<h5reader::io::DftFrame>& dftFrames,
    QString* error) {
    summary_ = {};
    paths_.clear();
    summary_.dftFramesDeclared = static_cast<int>(dftFrames.size());

    if (protein.atomCount() == 0 || protein.ringCount() == 0) {
        if (error)
            *error = QStringLiteral("protein has no atoms or rings");
        return false;
    }
    if (conformation.frameCount() == 0) {
        if (error)
            *error = QStringLiteral("conformation has no frames");
        return false;
    }
    if (options_.atom && *options_.atom >= protein.atomCount()) {
        if (error)
            *error = QStringLiteral("atom out of range");
        return false;
    }
    if (options_.ring && *options_.ring >= protein.ringCount()) {
        if (error)
            *error = QStringLiteral("ring out of range");
        return false;
    }
    if (!std::isfinite(options_.surfaceToleranceA) || options_.surfaceToleranceA < 0.0) {
        if (error)
            *error = QStringLiteral("surfaceToleranceA must be finite and >= 0");
        return false;
    }
    if (!std::isfinite(options_.predictorZeroToleranceA3) ||
        options_.predictorZeroToleranceA3 < 0.0) {
        if (error)
            *error = QStringLiteral("predictorZeroToleranceA3 must be finite and >= 0");
        return false;
    }
    if (options_.minSamples < 3) {
        if (error)
            *error = QStringLiteral("minSamples must be >= 3");
        return false;
    }
    if (options_.minSamplesPerLobe < 1) {
        if (error)
            *error = QStringLiteral("minSamplesPerLobe must be >= 1");
        return false;
    }
    if (!std::isfinite(options_.minPredictorSpanA3) ||
        options_.minPredictorSpanA3 < 0.0 ||
        !std::isfinite(options_.minAbsLobePredictorA3) ||
        options_.minAbsLobePredictorA3 < 0.0) {
        if (error)
            *error = QStringLiteral("minimum lobe thresholds must be finite and >= 0");
        return false;
    }
    if (options_.maxPaths < 0 || options_.nullShiftCount < 0) {
        if (error)
            *error = QStringLiteral("maxPaths and nullShiftCount must be >= 0");
        return false;
    }

    int firstOriginalFrame =
        static_cast<int>(conformation.originalFrameIndex(0));
    int lastOriginalFrame = firstOriginalFrame;
    for (std::size_t row = 1; row < conformation.frameCount(); ++row) {
        const int original =
            static_cast<int>(conformation.originalFrameIndex(row));
        firstOriginalFrame = std::min(firstOriginalFrame, original);
        lastOriginalFrame = std::max(lastOriginalFrame, original);
    }
    const int start = options_.startFrame.value_or(firstOriginalFrame);
    const int end = options_.endFrame.value_or(lastOriginalFrame);
    if (start < firstOriginalFrame || end > lastOriginalFrame || end < start) {
        if (error)
            *error = QStringLiteral("frame range out of range");
        return false;
    }

    std::vector<DftCandidate> candidates;
    candidates.reserve(dftFrames.size());
    for (const h5reader::io::DftFrame& declared : dftFrames) {
        if (declared.frame_index < start || declared.frame_index > end)
            continue;
        if (declared.frame_index < 0) {
            ++summary_.dftFramesSkipped;
            continue;
        }
        const std::optional<std::size_t> frameRow =
            conformation.frameRowForOriginalIndex(
                static_cast<std::size_t>(declared.frame_index));
        if (!frameRow) {
            ++summary_.dftFramesSkipped;
            continue;
        }
        candidates.push_back({&declared, *frameRow});
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const DftCandidate& a, const DftCandidate& b) {
        return a.declared->frame_index < b.declared->frame_index;
    });

    std::vector<DftCandidate> uniqueCandidates;
    uniqueCandidates.reserve(candidates.size());
    std::optional<int> lastCandidateFrame;
    for (const DftCandidate& candidate : candidates) {
        if (lastCandidateFrame &&
            *lastCandidateFrame == candidate.declared->frame_index)
            continue;
        lastCandidateFrame = candidate.declared->frame_index;
        uniqueCandidates.push_back(candidate);
    }
    candidates = std::move(uniqueCandidates);

    const std::vector<std::size_t> atoms = atomScanList(protein, options_);
    const std::vector<std::size_t> rings = ringScanList(protein, options_);
    summary_.atomsScanned = static_cast<int>(atoms.size());
    summary_.ringsScanned = static_cast<int>(rings.size());
    const int totalPaths =
        static_cast<int>(atoms.size() * rings.size());
    if (atoms.empty() || rings.empty() || candidates.empty()) {
        summary_.complete = true;
        return true;
    }

    std::vector<std::optional<LoadedDftFrame>> loadedFrameCache(candidates.size());
    std::vector<bool> dftAttempted(candidates.size(), false);
    auto loadedDftAt = [&](std::size_t frameOrdinal) -> const LoadedDftFrame* {
        if (frameOrdinal >= candidates.size())
            return nullptr;
        if (!dftAttempted[frameOrdinal]) {
            dftAttempted[frameOrdinal] = true;
            std::optional<LoadedDftFrame> loaded =
                loadDftFrame(*candidates[frameOrdinal].declared, protein);
            if (!loaded) {
                ++summary_.dftFramesSkipped;
            } else {
                ++summary_.dftFramesLoaded;
                loadedFrameCache[frameOrdinal] = std::move(*loaded);
            }
        }
        if (!loadedFrameCache[frameOrdinal])
            return nullptr;
        return &*loadedFrameCache[frameOrdinal];
    };

    for (std::size_t ring : rings) {
        const QtRing& qtRing = protein.ring(ring);
        std::vector<RingGeometry> ringGeometries;
        ringGeometries.reserve(candidates.size());
        std::vector<std::vector<Vec3>> ringVertices;
        ringVertices.reserve(candidates.size());
        for (const DftCandidate& candidate : candidates) {
            ringVertices.push_back(
                RingVertices(conformation, ring, candidate.frameRow));
            ringGeometries.push_back(FitRingGeometry(ringVertices.back()));
        }

        for (std::size_t atom : atoms) {
            ++summary_.pathsConsidered;

            RingCurrentPathResult path;
            path.atom = atom;
            path.ring = ring;
            path.minPointDipoleGeometryA3 = std::numeric_limits<double>::infinity();
            path.maxPointDipoleGeometryA3 = -std::numeric_limits<double>::infinity();
            std::vector<std::size_t> sampleFrameOrdinals;

            int previousSign = 0;
            for (std::size_t frameOrdinal = 0; frameOrdinal < candidates.size(); ++frameOrdinal) {
                const DftCandidate& candidate = candidates[frameOrdinal];

                RingCurrentPathSample sample;
                sample.frameIndex = candidate.declared->frame_index;
                sample.timePs = candidate.declared->framePs();
                sample.geometry =
                    MeasureRingNull(ringGeometries[frameOrdinal],
                                    conformation.atomPosition(candidate.frameRow, atom),
                                    options_.surfaceToleranceA);
                if (!sample.geometry.valid)
                    continue;
                sample.pointDipoleGeometryA3 =
                    RingCurrentPathAnalysis::pointDipoleGeometryA3(sample.geometry);
                sample.inverseDistanceCubedA3 =
                    RingCurrentPathAnalysis::inverseDistanceCubedA3(sample.geometry);
                sample.angularFactor = sample.geometry.angularFactor;
                if (!std::isfinite(sample.pointDipoleGeometryA3) ||
                    !std::isfinite(sample.inverseDistanceCubedA3) ||
                    !std::isfinite(sample.angularFactor))
                    continue;

                path.minPointDipoleGeometryA3 =
                    std::min(path.minPointDipoleGeometryA3, sample.pointDipoleGeometryA3);
                path.maxPointDipoleGeometryA3 =
                    std::max(path.maxPointDipoleGeometryA3, sample.pointDipoleGeometryA3);
                const int currentSign =
                    signOf(sample.pointDipoleGeometryA3, options_.predictorZeroToleranceA3);
                if (currentSign > 0)
                    ++path.positivePredictorSamples;
                else if (currentSign < 0)
                    ++path.negativePredictorSamples;
                if (previousSign != 0 && currentSign != 0 && previousSign != currentSign)
                    ++path.predictorSignChanges;
                if (currentSign != 0)
                    previousSign = currentSign;

                sampleFrameOrdinals.push_back(frameOrdinal);
                path.samples.push_back(std::move(sample));
            }

            if (static_cast<int>(path.samples.size()) < options_.minSamples) {
                ++summary_.pathsRejectedForSamples;
                continue;
            }

            path.pointDipoleGeometrySpanA3 =
                path.maxPointDipoleGeometryA3 - path.minPointDipoleGeometryA3;
            path.hardLobeCrossing =
                path.positivePredictorSamples >= options_.minSamplesPerLobe &&
                path.negativePredictorSamples >= options_.minSamplesPerLobe &&
                path.predictorSignChanges > 0;
            if (!path.hardLobeCrossing) {
                ++summary_.pathsRejectedForHardCrossing;
                continue;
            }
            if (path.pointDipoleGeometrySpanA3 < options_.minPredictorSpanA3 ||
                path.maxPointDipoleGeometryA3 < options_.minAbsLobePredictorA3 ||
                -path.minPointDipoleGeometryA3 < options_.minAbsLobePredictorA3) {
                ++summary_.pathsRejectedForWeakLobes;
                continue;
            }

            bool hasDftForAllSamples = true;
            for (std::size_t i = 0; i < path.samples.size(); ++i) {
                const std::size_t frameOrdinal = sampleFrameOrdinals[i];
                const LoadedDftFrame* frame = loadedDftAt(frameOrdinal);
                if (!frame || !frame->shielding || atom >= frame->shielding->atoms.size()) {
                    hasDftForAllSamples = false;
                    break;
                }
                RingCurrentPathSample& sample = path.samples[i];
                sample.orca = frame->shielding->atoms[atom];
                sample.biotSavart = h5reader::calculators::EvaluateShielding(
                    sample.geometry.atomPosition,
                    ringGeometries[frameOrdinal],
                    ringVertices[frameOrdinal],
                    qtRing.JohnsonBoveyLobeOffset(),
                    qtRing.LiteratureIntensity());
                if (!std::isfinite(sample.biotSavart.T0)) {
                    hasDftForAllSamples = false;
                    break;
                }
            }
            if (!hasDftForAllSamples) {
                ++summary_.pathsRejectedForSamples;
                continue;
            }

            path.orcaTotalT0VsPointDipole =
                fitLinear(path.samples, {pointDipoleValue, orcaTotalT0}, options_.nullShiftCount);
            path.orcaDiamagneticT0VsPointDipole =
                fitLinear(path.samples, {pointDipoleValue, orcaDiaT0}, options_.nullShiftCount);
            path.orcaParamagneticT0VsPointDipole =
                fitLinear(path.samples, {pointDipoleValue, orcaParaT0}, options_.nullShiftCount);
            path.orcaTotalT2MagnitudeVsPointDipole =
                fitLinear(path.samples, {pointDipoleValue, orcaTotalT2Magnitude},
                          options_.nullShiftCount);
            path.orcaTotalT0VsBiotSavartT0 =
                fitLinear(path.samples, {biotSavartT0, orcaTotalT0},
                          options_.nullShiftCount);
            path.biotSavartT0VsPointDipole =
                fitLinear(path.samples, {pointDipoleValue, biotSavartT0},
                          options_.nullShiftCount);
            path.orcaTotalT0VsDistanceOnly =
                fitLinear(path.samples, {distanceValue, orcaTotalT0}, options_.nullShiftCount);
            path.orcaTotalT0VsAngleOnly =
                fitLinear(path.samples, {angularValue, orcaTotalT0}, options_.nullShiftCount);

            paths_.push_back(std::move(path));
            if (options_.maxPaths > 0 &&
                static_cast<int>(paths_.size()) >= options_.maxPaths) {
                summary_.pathCount = static_cast<int>(paths_.size());
                summary_.truncatedByMaxPaths = summary_.pathsConsidered < totalPaths;
                summary_.complete = !summary_.truncatedByMaxPaths;
                return true;
            }
        }
    }

    summary_.pathCount = static_cast<int>(paths_.size());
    summary_.complete = true;
    return true;
}

}  // namespace h5reader::model

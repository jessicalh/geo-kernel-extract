// ring_current_case_tests -- C++ vet for the TYR24 HD2 / PHE5 poster case.
//
// This intentionally runs through the reader's own loader, topology,
// ConformationGeometry, DftShieldingLoader, and RingCurrentPathAnalysis. The goal
// is not a second pretty chart; it is a C++ assertion that the poster point cloud
// and regression are supported by the same machinery the app uses.

#include "io/QtProteinLoader.h"

#include "model/Conformation.h"
#include "model/ConformationGeometry.h"
#include "model/QtProtein.h"
#include "model/RingCurrentPathAnalysis.h"

#include <QFileInfo>
#include <QVector>
#include <QtTest>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <vector>

using h5reader::io::QtProteinLoader;
using h5reader::model::Conformation;
using h5reader::model::QtProtein;
using h5reader::model::RingCurrentPathAnalysis;
using h5reader::model::RingCurrentPathAnalysisOptions;
using h5reader::model::RingCurrentPathSample;
using h5reader::model::RingLocalFrame;
using h5reader::model::RingNullMeasurement;
using h5reader::model::Vec3;

namespace {

constexpr std::size_t kAtomTyr24Hd2 = 348;
constexpr std::size_t kRingPhe5 = 1;
constexpr std::size_t kReferenceFrame = 856;
constexpr double kColorScaleForPoster = 0.055;

QString fixturePath() {
    const QString env = qEnvironmentVariable("H5READER_RINGCURRENT_VET_FIXTURE");
    if (!env.isEmpty())
        return env;
    return QStringLiteral(
        "C:/projects/reader-data/1p9j-calibration-with-dft/"
        "1p9j-calibration-dense-mopac-live-orca.LGS");
}

bool explicitFixtureWasRequested() {
    return qEnvironmentVariableIsSet("H5READER_RINGCURRENT_VET_FIXTURE");
}

void requireNear(double actual, double expected, double tolerance, const char* label) {
    const double delta = std::abs(actual - expected);
    QVERIFY2(delta <= tolerance,
             qPrintable(QStringLiteral("%1 actual=%2 expected=%3 delta=%4 tolerance=%5")
                            .arg(QString::fromLatin1(label))
                            .arg(actual, 0, 'g', 17)
                            .arg(expected, 0, 'g', 17)
                            .arg(delta, 0, 'g', 17)
                            .arg(tolerance, 0, 'g', 17)));
}

double expectedValue(const RingNullMeasurement& m) {
    return RingCurrentPathAnalysis::pointDipoleGeometryA3(m);
}

double pearson(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.size() < 3)
        return std::numeric_limits<double>::quiet_NaN();
    const double meanX = std::accumulate(x.begin(), x.end(), 0.0) /
                         static_cast<double>(x.size());
    const double meanY = std::accumulate(y.begin(), y.end(), 0.0) /
                         static_cast<double>(y.size());
    double sxx = 0.0;
    double syy = 0.0;
    double sxy = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        const double dx = x[i] - meanX;
        const double dy = y[i] - meanY;
        sxx += dx * dx;
        syy += dy * dy;
        sxy += dx * dy;
    }
    if (sxx <= 0.0 || syy <= 0.0)
        return std::numeric_limits<double>::quiet_NaN();
    return sxy / std::sqrt(sxx * syy);
}

}  // namespace

class RingCurrentCaseTests : public QObject {
    Q_OBJECT

private slots:
    void tyr24Hd2Phe5FullDftCase();
};

void RingCurrentCaseTests::tyr24Hd2Phe5FullDftCase() {
    const QString fixture = fixturePath();
    if (!QFileInfo::exists(fixture)) {
        if (explicitFixtureWasRequested()) {
            QFAIL(qPrintable(QStringLiteral(
                "H5READER_RINGCURRENT_VET_FIXTURE does not exist: %1").arg(fixture)));
        }
        QSKIP(qPrintable(QStringLiteral("1P9J local DFT fixture not present: %1").arg(fixture)));
    }

    auto loaded = QtProteinLoader::LoadRunPath(fixture);
    QVERIFY2(loaded.ok, qPrintable(loaded.error));
    QVERIFY(loaded.protein != nullptr);
    QVERIFY(loaded.conformation != nullptr);
    QVERIFY(loaded.manifest.dft.has_value());

    const QtProtein& protein = *loaded.protein;
    const Conformation& conf = *loaded.conformation;
    QCOMPARE(protein.atomCount(), std::size_t{846});
    QCOMPARE(conf.frameCount(), std::size_t{1501});
    QCOMPARE(loaded.manifest.dft->frames.size(), std::size_t{751});

    const auto& atom = protein.atom(kAtomTyr24Hd2);
    QCOMPARE(atom.atomIndex, static_cast<int32_t>(kAtomTyr24Hd2));
    QCOMPARE(atom.parentAtomIndex, int32_t{347});
    QVERIFY(atom.IsAromaticRingHydrogen());
    QCOMPARE(protein.atomNames(kAtomTyr24Hd2).iupac, QStringLiteral("HD2"));

    const auto& ring = protein.ring(kRingPhe5);
    QCOMPARE(ring.ringId, int32_t{1});
    QCOMPARE(ring.parentResidueNumber, int32_t{5});
    QCOMPARE(QString::fromLatin1(ring.TypeName()), QStringLiteral("PHE"));
    const QVector<int32_t> expectedRingWalk{70, 71, 73, 75, 77, 79};
    QCOMPARE(static_cast<int>(ring.atomIndices.size()), expectedRingWalk.size());
    for (int i = 0; i < expectedRingWalk.size(); ++i)
        QCOMPARE(ring.atomIndices[static_cast<std::size_t>(i)], expectedRingWalk[i]);

    RingCurrentPathAnalysisOptions options;
    options.atom = kAtomTyr24Hd2;
    options.ring = kRingPhe5;
    options.startFrame = 0;
    options.endFrame = 1500;
    options.minSamples = 700;
    options.minSamplesPerLobe = 50;
    options.minPredictorSpanA3 = 0.15;
    options.minAbsLobePredictorA3 = 0.015;
    options.maxPaths = 1;
    options.nullShiftCount = 750;

    RingCurrentPathAnalysis analysis(options);
    QString collectError;
    QVERIFY2(analysis.run(protein, conf, loaded.manifest.dft->frames, &collectError),
             qPrintable(collectError));

    const auto& summary = analysis.summary();
    QVERIFY(summary.complete);
    QVERIFY(!summary.truncatedByMaxPaths);
    QCOMPARE(summary.dftFramesDeclared, 751);
    QCOMPARE(summary.dftFramesLoaded, 751);
    QCOMPARE(summary.dftFramesSkipped, 0);
    QCOMPARE(summary.pathsConsidered, 1);
    QCOMPARE(summary.pathCount, 1);
    QCOMPARE(analysis.paths().size(), std::size_t{1});

    const auto& entry = analysis.paths().front();
    QCOMPARE(entry.atom, kAtomTyr24Hd2);
    QCOMPARE(entry.ring, kRingPhe5);
    QVERIFY(entry.hardLobeCrossing);
    QCOMPARE(entry.samples.size(), std::size_t{751});
    QCOMPARE(entry.positivePredictorSamples, 546);
    QCOMPARE(entry.negativePredictorSamples, 205);
    QVERIFY(entry.predictorSignChanges > 0);
    requireNear(entry.minPointDipoleGeometryA3, -0.0184166372550467, 1e-9,
                "min expected relationship");
    requireNear(entry.maxPointDipoleGeometryA3, 0.1535686443542942, 1e-9,
                "max expected relationship");
    QVERIFY(entry.pointDipoleGeometrySpanA3 > 0.17);

    const auto& fit = entry.orcaTotalT0VsPointDipole;
    QVERIFY(fit.valid);
    QCOMPARE(fit.sampleCount, 751);
    QVERIFY(fit.scale > 0.0);
    requireNear(fit.scale, 25.85396853, 1e-5, "ORCA total T0 scale");
    requireNear(fit.intercept, 24.45256753, 1e-5, "ORCA total T0 intercept");
    requireNear(fit.correlation, 0.9224617929, 1e-6, "ORCA total T0 correlation");
    requireNear(fit.r2, 0.8509357593, 1e-6, "ORCA total T0 r2");
    QCOMPARE(fit.nullShiftCount, 750);
    QVERIFY(fit.nullMedianR2 < 0.07);
    QVERIFY(fit.nullMaxR2 < 0.45);
    QVERIFY(fit.nullGeRealFraction < 0.002);

    QVERIFY(entry.orcaTotalT0VsDistanceOnly.valid);
    QVERIFY(entry.orcaTotalT0VsAngleOnly.valid);
    QVERIFY(fit.r2 > entry.orcaTotalT0VsDistanceOnly.r2 + 0.15);
    QVERIFY(fit.r2 > entry.orcaTotalT0VsAngleOnly.r2 + 0.30);
    QVERIFY(entry.biotSavartT0VsPointDipole.valid);

    int positive = 0;
    int negative = 0;
    int saturatedByPosterScale = 0;
    int signAgreement = 0;
    double sigmaPositiveSum = 0.0;
    double sigmaNegativeSum = 0.0;
    std::vector<double> absFactor;
    std::vector<double> absCenteredSigma;
    absFactor.reserve(entry.samples.size());
    absCenteredSigma.reserve(entry.samples.size());

    for (const RingCurrentPathSample& sample : entry.samples) {
        const double expected = sample.pointDipoleGeometryA3;
        const double sigma = sample.orca.total.T0;
        if (expected > 0.0) {
            ++positive;
            sigmaPositiveSum += sigma;
        } else if (expected < 0.0) {
            ++negative;
            sigmaNegativeSum += sigma;
        }
        if (std::abs(expected) > kColorScaleForPoster)
            ++saturatedByPosterScale;

        const double centered = sigma - fit.intercept;
        if ((centered > 0.0 && expected > 0.0) ||
            (centered < 0.0 && expected < 0.0)) {
            ++signAgreement;
        }
        absFactor.push_back(std::abs(expected));
        absCenteredSigma.push_back(std::abs(centered));
    }
    QCOMPARE(positive, 546);
    QCOMPARE(negative, 205);
    QCOMPARE(saturatedByPosterScale, 119);
    const double sigmaMeanPositive = sigmaPositiveSum / static_cast<double>(positive);
    const double sigmaMeanNegative = sigmaNegativeSum / static_cast<double>(negative);
    requireNear(sigmaMeanPositive, 25.23483333, 1e-5, "positive lobe sigma mean");
    requireNear(sigmaMeanNegative, 24.10895122, 1e-5, "negative lobe sigma mean");
    QVERIFY(sigmaMeanPositive - sigmaMeanNegative > 1.0);
    QVERIFY(static_cast<double>(signAgreement) / static_cast<double>(entry.samples.size()) > 0.84);
    QVERIFY(pearson(absFactor, absCenteredSigma) > 0.92);

    const RingLocalFrame reference =
        h5reader::model::RingLocalFrameAt(conf, kRingPhe5, kReferenceFrame);
    QVERIFY(reference.valid);
    double maxKernelDelta = 0.0;
    double maxDistanceDelta = 0.0;
    double maxNullMarginDelta = 0.0;
    double maxRoundtripDelta = 0.0;
    for (const RingCurrentPathSample& sample : entry.samples) {
        const std::size_t frame = static_cast<std::size_t>(sample.frameIndex);
        const RingLocalFrame source =
            h5reader::model::RingLocalFrameAt(conf, kRingPhe5, frame);
        QVERIFY(source.valid);
        const Vec3 sourcePosition = conf.atomPosition(frame, kAtomTyr24Hd2);
        const Vec3 local = h5reader::model::ToRingLocal(source, sourcePosition);
        const Vec3 drawn = h5reader::model::FromRingLocal(reference, local);
        const RingNullMeasurement sourceMeasure =
            h5reader::model::MeasureRingNull(source.geometry, sourcePosition);
        const RingNullMeasurement projectedMeasure =
            h5reader::model::MeasureRingNull(reference.geometry, drawn);
        QVERIFY(sourceMeasure.valid);
        QVERIFY(projectedMeasure.valid);

        const Vec3 roundtrip = h5reader::model::ToRingLocal(reference, drawn);
        maxKernelDelta = std::max(maxKernelDelta,
                                  std::abs(expectedValue(sourceMeasure) -
                                           expectedValue(projectedMeasure)));
        maxDistanceDelta = std::max(maxDistanceDelta,
                                    std::abs(sourceMeasure.distanceA -
                                             projectedMeasure.distanceA));
        maxNullMarginDelta = std::max(maxNullMarginDelta,
                                      std::abs(sourceMeasure.nullMarginA -
                                               projectedMeasure.nullMarginA));
        maxRoundtripDelta = std::max(maxRoundtripDelta, (roundtrip - local).norm());
    }
    QVERIFY(maxKernelDelta < 1e-12);
    QVERIFY(maxDistanceDelta < 1e-12);
    QVERIFY(maxNullMarginDelta < 1e-12);
    QVERIFY(maxRoundtripDelta < 1e-12);

    qInfo().noquote()
        << "TYR24_HD2_PHE5_CXX_VET"
        << "samples=" << entry.samples.size()
        << "positive=" << positive
        << "negative=" << negative
        << "r=" << fit.correlation
        << "r2=" << fit.r2
        << "scale=" << fit.scale
        << "null_ge_real_fraction=" << fit.nullGeRealFraction
        << "max_kernel_delta=" << maxKernelDelta;
}

QTEST_GUILESS_MAIN(RingCurrentCaseTests)

#include "ring_current_case_tests.moc"

// h5reader_occupancy_shells_math_tests — QtTest binary for the pure density /
// HDR math behind the per-atom occupation-probability envelope shells
// (app/OccupancyShellsMath.h). Pure functions over model::Vec3/Mat3; no VTK, no
// Qt widgets, no I/O — links Qt6::Test + Eigen3 only.
//
// Covers the checklist the two pre-build reviews asked for:
//   - autocorrelation window + statistical inefficiency g + n_eff clamp
//   - anisotropic bandwidth: isotropic scaling + thin-axis eigenvalue floor
//   - grid coarsen for a very mobile atom, and the floor-AFTER-coarsen ordering
//   - HDR levels validated against an analytic isotropic Gaussian (chi-square)
//   - HDR tie handling: the crossing includes the whole tie-shell, not a voxel
//   - rigid / too-few-frames / degenerate skips

#include "app/OccupancyShellsMath.h"

#include <QObject>
#include <QtTest>

#include <algorithm>
#include <cmath>
#include <vector>

using namespace h5reader;

namespace {

bool nearly(double a, double b, double tol = 1e-9) {
    return std::abs(a - b) <= tol;
}
bool relClose(double a, double b, double rel) {
    return std::abs(a - b) <= rel * std::max(std::abs(a), std::abs(b)) + 1e-12;
}
double minEig(const model::Mat3& M) {
    Eigen::SelfAdjointEigenSolver<model::Mat3> es(M);
    return es.eigenvalues()(0);
}
double maxEigVal(const model::Mat3& M) {
    Eigen::SelfAdjointEigenSolver<model::Mat3> es(M);
    return es.eigenvalues()(2);
}

}  // namespace

class OccupancyShellsMathTests : public QObject {
    Q_OBJECT

private slots:
    void testMotionStats_meanCovarianceRmsf();
    void testMotionStats_uncorrelatedNeffNearT();
    void testMotionStats_correlatedNeffBelowT();
    void testMotionStats_degenerateZeroVariance();

    void testBandwidth_isotropicScaling();
    void testBandwidth_floorsThinAxis();

    void testPlanGrid_normalKeepsTargetSpacing();
    void testPlanGrid_hugeSpreadCoarsensWithinCap();

    void testHdr_isotropicGaussianMatchesChiSquare();
    void testHdr_tiesBracketWholeShell();
    void testHdr_nestingOnGaussian();

    void testComputeOccupancy_endToEndMobile();
    void testComputeOccupancy_rigidSkipped();
    void testComputeOccupancy_tooFewFrames();
    void testComputeOccupancy_coarsenRefloorOrdering();
};

// ---- MotionStats -------------------------------------------------------

void OccupancyShellsMathTests::testMotionStats_meanCovarianceRmsf() {
    // Six symmetric offsets about a known mean: ±3 x, ±2 y, ±1 z.
    // Unbiased (1/(T-1)=1/5) variances: 2*9/5=3.6, 2*4/5=1.6, 2*1/5=0.4.
    const model::Vec3 base(10.0, 20.0, 30.0);
    std::vector<model::Vec3> pos = {
        base + model::Vec3( 3, 0, 0), base + model::Vec3(-3, 0, 0),
        base + model::Vec3( 0, 2, 0), base + model::Vec3( 0,-2, 0),
        base + model::Vec3( 0, 0, 1), base + model::Vec3( 0, 0,-1),
    };
    const auto s = math::computeMotionStats(pos);
    QVERIFY(s.valid);
    QVERIFY(!s.degenerate);
    QVERIFY(nearly(s.mean.x(), 10.0));
    QVERIFY(nearly(s.mean.y(), 20.0));
    QVERIFY(nearly(s.mean.z(), 30.0));
    QVERIFY(nearly(s.covariance(0, 0), 3.6, 1e-9));
    QVERIFY(nearly(s.covariance(1, 1), 1.6, 1e-9));
    QVERIFY(nearly(s.covariance(2, 2), 0.4, 1e-9));
    QVERIFY(nearly(s.covariance(0, 1), 0.0, 1e-9));
    QVERIFY(nearly(s.rmsf, std::sqrt(3.6 + 1.6 + 0.4), 1e-9));
}

void OccupancyShellsMathTests::testMotionStats_uncorrelatedNeffNearT() {
    // Alternating ±d about the mean: C(1) is negative, so the window closes
    // at lag 1 with no positive terms -> g = 1, n_eff = T.
    const std::size_t T = 50;
    std::vector<model::Vec3> pos;
    pos.reserve(T);
    for (std::size_t t = 0; t < T; ++t)
        pos.emplace_back((t % 2 == 0) ? model::Vec3(1, 0, 0) : model::Vec3(-1, 0, 0));
    const auto s = math::computeMotionStats(pos);
    QVERIFY(s.valid);
    QVERIFY(nearly(s.statineff, 1.0, 1e-9));
    QVERIFY(nearly(s.nEff, static_cast<double>(T), 1e-9));
}

void OccupancyShellsMathTests::testMotionStats_correlatedNeffBelowT() {
    // A slow sinusoid is strongly autocorrelated: n_eff must drop well below T.
    const std::size_t T = 100;
    std::vector<model::Vec3> pos;
    pos.reserve(T);
    for (std::size_t t = 0; t < T; ++t)
        pos.emplace_back(5.0 * std::sin(2.0 * M_PI * static_cast<double>(t) / T), 0.0, 0.0);
    const auto s = math::computeMotionStats(pos);
    QVERIFY(s.valid);
    QVERIFY(s.statineff > 1.0);
    QVERIFY(s.nEff >= 1.0);
    QVERIFY(s.nEff <= static_cast<double>(T));
    QVERIFY(s.nEff < 0.5 * static_cast<double>(T));  // clearly correlated
}

void OccupancyShellsMathTests::testMotionStats_degenerateZeroVariance() {
    std::vector<model::Vec3> pos(20, model::Vec3(1.0, 2.0, 3.0));
    const auto s = math::computeMotionStats(pos);
    QVERIFY(s.valid);
    QVERIFY(s.degenerate);
    QVERIFY(nearly(s.rmsf, 0.0, 1e-12));
}

// ---- Bandwidth matrix --------------------------------------------------

void OccupancyShellsMathTests::testBandwidth_isotropicScaling() {
    // cov = 4 I, n_eff = 100. Data eigenvalue (scale*4) sits above the floor,
    // so H = scale*4 * I with the Silverman 3-D factor.
    const model::Mat3 cov = 4.0 * model::Mat3::Identity();
    const double nEff = 100.0, voxel = 0.4, floorF = 1.5;
    const double scale = std::pow(0.8, 2.0 / 7.0) * std::pow(nEff, -2.0 / 7.0);
    const double expected = scale * 4.0;
    QVERIFY(expected > floorF * voxel * floorF * voxel);  // not floored
    const model::Mat3 H = math::bandwidthMatrix(cov, nEff, voxel, floorF);
    QVERIFY(nearly(H(0, 0), expected, 1e-9));
    QVERIFY(nearly(H(0, 1), 0.0, 1e-9));
    QVERIFY(nearly(minEig(H), expected, 1e-9));
    QVERIFY(nearly(maxEigVal(H), expected, 1e-9));
}

void OccupancyShellsMathTests::testBandwidth_floorsThinAxis() {
    // Planar motion: var in x,y; ~0 in z. The thin z axis must be floored to
    // (floorFactor*voxel)^2 while x,y stay data-driven.
    model::Mat3 cov = model::Mat3::Zero();
    cov(0, 0) = 4.0;
    cov(1, 1) = 4.0;
    cov(2, 2) = 1e-8;
    const double nEff = 100.0, voxel = 0.4, floorF = 1.5;
    const double floor = (floorF * voxel) * (floorF * voxel);  // 0.36
    const double scale = std::pow(0.8, 2.0 / 7.0) * std::pow(nEff, -2.0 / 7.0);
    const model::Mat3 H = math::bandwidthMatrix(cov, nEff, voxel, floorF);
    QVERIFY(nearly(minEig(H), floor, 1e-9));           // thin axis floored
    QVERIFY(nearly(maxEigVal(H), scale * 4.0, 1e-9));  // fat axis data-driven
}

// ---- Grid planning -----------------------------------------------------

void OccupancyShellsMathTests::testPlanGrid_normalKeepsTargetSpacing() {
    std::vector<model::Vec3> pos = {
        model::Vec3( 2,  2,  2), model::Vec3(-2, -2, -2),
        model::Vec3( 2, -2,  2), model::Vec3(-2,  2, -2),
    };
    const model::Mat3 H = 0.25 * model::Mat3::Identity();
    const auto g = math::planGrid(pos, H, /*voxelTarget=*/0.4, /*maxDim=*/96, /*marginSigma=*/3.5);
    QVERIFY(nearly(g.spacing, 0.4, 1e-9));
    for (int a = 0; a < 3; ++a) {
        QVERIFY(g.dims[a] >= 2);
        QVERIFY(g.dims[a] <= 96);
    }
}

void OccupancyShellsMathTests::testPlanGrid_hugeSpreadCoarsensWithinCap() {
    std::vector<model::Vec3> pos = {
        model::Vec3( 50,  50,  50), model::Vec3(-50, -50, -50),
        model::Vec3( 50, -50,  50), model::Vec3(-50,  50, -50),
    };
    const model::Mat3 H = 0.25 * model::Mat3::Identity();
    const auto g = math::planGrid(pos, H, 0.4, 96, 3.5);
    QVERIFY(g.spacing > 0.4);  // coarsened
    for (int a = 0; a < 3; ++a) QVERIFY(g.dims[a] <= 96);
}

// ---- HDR levels --------------------------------------------------------

// Build an analytic isotropic-Gaussian density field on a fine grid centred at
// the origin: peak value 1 at centre, value exp(-r^2/2sigma^2) elsewhere.
static math::DensityField makeGaussianField(double sigma, double voxel, double halfWidth) {
    math::DensityField f;
    const int n = static_cast<int>(std::round(2.0 * halfWidth / voxel)) + 1;
    f.grid.spacing = voxel;
    f.grid.origin = model::Vec3(-halfWidth, -halfWidth, -halfWidth);
    f.grid.dims = {n, n, n};
    f.values.assign(static_cast<std::size_t>(n) * n * n, 0.0);
    const double s2 = sigma * sigma;
    for (int k = 0; k < n; ++k)
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i) {
                const model::Vec3 p = f.grid.pointAt(i, j, k);
                f.values[static_cast<std::size_t>(f.grid.index(i, j, k))] =
                    std::exp(-0.5 * p.squaredNorm() / s2);
            }
    f.valid = true;
    return f;
}

void OccupancyShellsMathTests::testHdr_isotropicGaussianMatchesChiSquare() {
    // For a 3-D isotropic Gaussian, the highest-density region holding mass p is
    // a ball of radius r_p = sigma*sqrt(chi2_3^{-1}(p)); the density on that
    // shell is peak*exp(-chi2/2). Verify both the enclosed mass bracket and the
    // iso-level against the analytic chi-square values.
    const math::DensityField f = makeGaussianField(/*sigma=*/2.0, /*voxel=*/0.2, /*half=*/8.0);
    double total = 0.0;
    for (double v : f.values) total += v;
    total *= f.grid.voxelVolume();

    const auto levels = math::hdrLevels(f, {0.5, 0.9});
    QCOMPARE(static_cast<int>(levels.size()), 2);

    const double chi2_50 = 2.3660, chi2_90 = 6.2514;  // chi-square_3 quantiles
    const double peak = 1.0;
    const double iso50 = peak * std::exp(-0.5 * chi2_50);  // 0.3064
    const double iso90 = peak * std::exp(-0.5 * chi2_90);  // 0.0439

    for (const auto& lvl : levels) {
        QVERIFY(lvl.valid);
        // Honest bracket: mass(>c) < frac*total <= mass(>=c).
        QVERIFY(lvl.strictlyAboveMass < lvl.fraction * total + 1e-9);
        QVERIFY(lvl.includedMass >= lvl.fraction * total - 1e-9);
        // No gross overshoot on a fine grid.
        QVERIFY(lvl.includedMass <= (lvl.fraction + 0.04) * total);
    }
    QVERIFY(relClose(levels[0].isoValue, iso50, 0.10));
    QVERIFY(relClose(levels[1].isoValue, iso90, 0.10));
}

void OccupancyShellsMathTests::testHdr_tiesBracketWholeShell() {
    // Ten points of equal density 1.0, rest zero. For frac=0.45 the crossing
    // lands inside the tie-shell: c_X = 1.0, mass(>1.0)=0 (< target), and
    // mass(>=1.0) = the WHOLE shell (10 voxels), not one.
    math::DensityField f;
    f.grid.spacing = 1.0;
    f.grid.origin = model::Vec3::Zero();
    f.grid.dims = {100, 1, 1};
    f.values.assign(100, 0.0);
    for (int i = 0; i < 10; ++i) f.values[i] = 1.0;
    f.valid = true;

    const double dV = f.grid.voxelVolume();
    const double total = 10.0 * dV;
    const auto levels = math::hdrLevels(f, {0.45});
    QCOMPARE(static_cast<int>(levels.size()), 1);
    QVERIFY(nearly(levels[0].isoValue, 1.0, 1e-12));
    QVERIFY(nearly(levels[0].strictlyAboveMass, 0.0, 1e-12));   // nothing above 1.0
    QVERIFY(nearly(levels[0].includedMass, total, 1e-12));      // whole tie-shell
    // Bracket holds: 0 < 0.45*total <= total.
    QVERIFY(levels[0].strictlyAboveMass < 0.45 * total);
    QVERIFY(levels[0].includedMass >= 0.45 * total);
}

void OccupancyShellsMathTests::testHdr_nestingOnGaussian() {
    const math::DensityField f = makeGaussianField(2.0, 0.25, 8.0);
    const auto levels = math::hdrLevels(f, {0.5, 0.9});
    QCOMPARE(static_cast<int>(levels.size()), 2);
    // Inner (50%) sits at a HIGHER density level than outer (90%), and encloses
    // LESS mass -> the 50% region is strictly inside the 90% region.
    QVERIFY(levels[0].isoValue > levels[1].isoValue);
    QVERIFY(levels[0].includedMass < levels[1].includedMass);
}

// ---- Orchestrator ------------------------------------------------------

void OccupancyShellsMathTests::testComputeOccupancy_endToEndMobile() {
    const std::size_t T = 200;
    std::vector<model::Vec3> pos;
    pos.reserve(T);
    for (std::size_t t = 0; t < T; ++t) {
        const double th = 2.0 * M_PI * static_cast<double>(t) / T;
        pos.emplace_back(3.0 * std::cos(th), 2.0 * std::sin(th), 1.0 * std::sin(3.0 * th));
    }
    const auto r = math::computeOccupancy(pos);
    QVERIFY(r.computed);
    QVERIFY(r.stats.rmsf > 0.8);
    QCOMPARE(static_cast<int>(r.shells.size()), 2);
    QVERIFY(r.shells[0].valid && r.shells[1].valid);
    QVERIFY(r.shells[0].isoValue > r.shells[1].isoValue);          // 50% inside 90%
    QVERIFY(r.shells[0].includedMass < r.shells[1].includedMass);
    double total = 0.0;
    for (double v : r.field.values) total += v;
    total *= r.field.grid.voxelVolume();
    QVERIFY(nearly(total, 1.0, 1e-6));  // normalised field
    for (const auto& lvl : r.shells) {
        QVERIFY(lvl.strictlyAboveMass < lvl.fraction * total + 1e-9);
        QVERIFY(lvl.includedMass >= lvl.fraction * total - 1e-9);
    }
}

void OccupancyShellsMathTests::testComputeOccupancy_rigidSkipped() {
    std::vector<model::Vec3> pos;
    for (std::size_t t = 0; t < 20; ++t)
        pos.emplace_back((t % 2 == 0) ? model::Vec3(0.05, 0, 0) : model::Vec3(-0.05, 0, 0));
    const auto r = math::computeOccupancy(pos);
    QVERIFY(!r.computed);
    QVERIFY(r.note.find("rigid") != std::string::npos);
}

void OccupancyShellsMathTests::testComputeOccupancy_tooFewFrames() {
    std::vector<model::Vec3> pos = {
        model::Vec3(0, 0, 0), model::Vec3(1, 0, 0), model::Vec3(0, 1, 0),
        model::Vec3(0, 0, 1), model::Vec3(1, 1, 1),
    };
    const auto r = math::computeOccupancy(pos);  // 5 < minFrames (8)
    QVERIFY(!r.computed);
    QVERIFY(r.note.find("few") != std::string::npos);
}

void OccupancyShellsMathTests::testComputeOccupancy_coarsenRefloorOrdering() {
    // Huge x/y spread forces the grid to coarsen; z is essentially flat. The
    // kernel floor must be recomputed at the FINAL (coarsened) spacing, so the
    // smallest eigenvalue of H == (1.5 * final_spacing)^2 — NOT (1.5*0.4)^2.
    const std::size_t T = 60;
    std::vector<model::Vec3> pos;
    pos.reserve(T);
    for (std::size_t t = 0; t < T; ++t) {
        const double th = 2.0 * M_PI * static_cast<double>(t) / T;
        pos.emplace_back(40.0 * std::cos(th), 40.0 * std::sin(1.3 * th), 0.003 * std::sin(th));
    }
    const auto r = math::computeOccupancy(pos);
    QVERIFY(r.computed);
    QVERIFY(r.field.grid.spacing > 0.5);  // coarsened beyond the 0.4 target

    const double expectedFloor =
        (1.5 * r.field.grid.spacing) * (1.5 * r.field.grid.spacing);
    QVERIFY(relClose(minEig(r.bandwidth), expectedFloor, 0.02));
    // Sanity: had the floor used the 0.4 target it would be ~0.36, far off.
    QVERIFY(!relClose(minEig(r.bandwidth), (1.5 * 0.4) * (1.5 * 0.4), 0.02));
}

QTEST_GUILESS_MAIN(OccupancyShellsMathTests)

#include "occupancy_shells_math_tests.moc"

// OccupancyShellsMath — the correctness core of the per-atom occupation-
// probability envelope shells, factored out as pure functions so it can be
// unit-tested directly (Opus/codex review: the math, not the rendering, is
// where a quiet error does the damage). It happens to need no VTK — but that is
// INCIDENTAL, not the goal. The goal is best practice, which here means "use VTK
// where VTK is the right tool, compute directly where it is not":
//   - Where VTK is right, the overlay uses it: it contours this field with
//     vtkContourFilter into translucent actors (the proven QtFieldGridOverlay
//     pattern). That is NOT reimplemented here.
//   - We compute directly only where VTK has no fitting primitive: the
//     ANISOTROPIC-COVARIANCE Gaussian kernel has no vtkGaussianSplatter
//     equivalent (splatter is uniaxial — one normal + eccentricity — and cannot
//     carry a general 3x3 covariance), and exact mass-bracket HDR levels beat
//     vtkImageAccumulate's histogram binning. Hand-rolling those is the right
//     tool, not VTK-avoidance.
// The grid uses vtkImageData's x-fastest point order (see GridSpec::index) so
// the overlay hands the field off with a single std::copy, no re-indexing.
// Modelled on FitTargetMath.h / PlaneFrameMath.h (tested with Qt6::Test + Eigen).
//
// The feature (spec notes/OCCUPANCY_SHELLS_SPEC_2026-06-13.md): given a focused
// atom's positions across all frames of an MD trajectory (already backbone-
// aligned by the caller via the Conformation seam), estimate where it sits as a
// kernel-density field, then extract nested isosurface levels enclosing 50% and
// 90% of the occupation mass (highest-density regions). The VTK overlay layers
// on top of this; ALL the correctness lives here, hence the separation.
//
// Pipeline (each step is a free function so it can be tested in isolation):
//   computeMotionStats  : positions      -> mean, covariance, RMSF, n_eff
//   bandwidthMatrix      : cov,n_eff,voxel-> anisotropic Gaussian kernel H
//   planGrid             : positions,H    -> a padded, dim-capped voxel grid
//   evaluateDensity      : positions,H,grid-> normalized point-sampled density
//   hdrLevels            : density,frac    -> iso-levels + honest bracketed mass
//   computeOccupancy     : the orchestrator wiring all of the above + skips
//
// Design choices that came out of two adversarial reviews (Opus + codex), with
// the rationale kept here so it is not re-litigated:
//   - ANISOTROPIC kernel (H proportional to the sample covariance), NOT an
//     isotropic scalar bandwidth: a near-planar / near-linear motion would make
//     a geometric-mean scale collapse to the voxel floor and the shell would be
//     a grid artifact rather than data. The covariance kernel matches the motion
//     shape; the eigenvalue floor only lifts the *thin* axes to the voxel scale.
//   - n_eff from the autocorrelation (statistical inefficiency g), NOT the raw
//     frame count: MD frames are temporally correlated, so T overstates the
//     information and would over-sharpen the kernel.
//   - HDR level from cumulative mass so "90%" MEANS 90% of the occupation mass.
//     On a discrete grid the guarantee is a BRACKET, not exact: lattice ties on
//     a symmetric Gaussian mean the crossing can include a whole tie-shell, so
//     we report the actual included mass rather than claiming "X% within a
//     voxel" (codex). The contoured isosurface (downstream, in VTK) interpolates
//     the boundary, so its geometric enclosure runs slightly under the
//     point-set mass and tightens as voxel -> 0 (acceptance test territory).

#pragma once

#include "../model/Types.h"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace h5reader::math {

// ===========================================================================
// Config — every threshold is here, documented, and derived rather than magic.
// ===========================================================================
struct OccupancyConfig {
    // Target voxel edge (Angstrom). The actual spacing may be coarsened to
    // satisfy maxDim for a very mobile atom (see planGrid).
    double voxelTarget = 0.4;
    // Per-axis grid point cap. A disordered terminus can wander 15-25 A; at
    // 0.4 A that is 60-80 points, so 96 leaves headroom while bounding the
    // worst case (the HDR sort is O(V log V)). Exceeding it coarsens the voxel.
    int maxDim = 96;
    // Grid pad AND kernel splat cutoff, in units of sqrt(largest eigenvalue of
    // H). 3.5 sigma encloses ~99.95% of a 1-D Gaussian; the 90% HDR radius of a
    // 3-D Gaussian is ~2.5 sigma, so 3.5 comfortably contains the outer shell of
    // even the outermost frames without clipping.
    double marginSigma = 3.5;
    // Kernel eigenvalue floor, in voxels (per principal axis): the smallest
    // kernel sigma is lifted to floorVoxelFactor * spacing so a thin/degenerate
    // motion axis is rendered at grid resolution rather than aliasing. 1.5
    // voxels keeps the kernel just-resolvable on the lattice.
    double floorVoxelFactor = 1.5;
    // Rigid-atom skip: if RMSF < rigidRmsfFactor * voxelTarget the shell would
    // be dominated by the kernel floor (an artifact, not motion), so we skip the
    // KDE and let the caller draw only a marker. Tied to the voxel/floor scale
    // (2 voxels), not a free magic number.
    double rigidRmsfFactor = 2.0;
    // Mass fractions to extract, ascending (inner first). 50% then 90%.
    std::vector<double> massFractions = {0.5, 0.9};
    // Need at least this many frames to estimate a covariance + autocorrelation.
    std::size_t minFrames = 8;
};

// ===========================================================================
// Per-atom motion statistics.
// ===========================================================================
struct MotionStats {
    bool valid = false;          // false => too few frames
    bool degenerate = false;     // zero variance (atom never moved)
    std::size_t frames = 0;
    model::Vec3 mean = model::Vec3::Zero();
    model::Mat3 covariance = model::Mat3::Zero();  // unbiased, 1/(T-1)
    double rmsf = 0.0;           // sqrt(trace(covariance)) — total positional spread
    double statineff = 1.0;      // g = 1 + 2*sum C(k), windowed; clamped >= 1
    double nEff = 0.0;           // T / g, clamped to [1, T]
};

// Compute mean, unbiased covariance, RMSF, and the correlation-corrected
// effective sample size. The autocorrelation uses the scalar fluctuation
// dot-product correlation with the biased (1/T) estimator (Sokal: the biased
// estimator has lower variance and is the standard choice for tau estimation),
// windowed at the first non-positive lag (initial-positive-sequence rule).
inline MotionStats computeMotionStats(const std::vector<model::Vec3>& pos) {
    MotionStats s;
    s.frames = pos.size();
    const std::size_t T = pos.size();
    if (T < 2)
        return s;  // valid stays false

    model::Vec3 mean = model::Vec3::Zero();
    for (const auto& p : pos) mean += p;
    mean /= static_cast<double>(T);
    s.mean = mean;

    // Unbiased covariance for the kernel; scalar variance (biased, 1/T) for the
    // autocorrelation normalisation gamma(0).
    model::Mat3 cov = model::Mat3::Zero();
    double var0 = 0.0;  // gamma(0) = (1/T) sum |dr|^2
    std::vector<model::Vec3> dr(T);
    for (std::size_t t = 0; t < T; ++t) {
        dr[t] = pos[t] - mean;
        cov += dr[t] * dr[t].transpose();
        var0 += dr[t].squaredNorm();
    }
    cov /= static_cast<double>(T - 1);
    var0 /= static_cast<double>(T);
    s.covariance = cov;
    s.rmsf = std::sqrt(std::max(0.0, cov.trace()));

    constexpr double kVarEps = 1e-18;  // below this the atom never moved
    if (var0 <= kVarEps) {
        s.degenerate = true;
        s.valid = true;
        s.nEff = static_cast<double>(T);
        return s;
    }

    // g = 1 + 2 * sum_{k>=1} C(k), windowed at the first non-positive C(k).
    // gamma(k) = (1/T) sum_{t=0}^{T-1-k} dr_t . dr_{t+k}  (biased).
    double sumC = 0.0;
    for (std::size_t k = 1; k < T; ++k) {
        double g_k = 0.0;
        for (std::size_t t = 0; t + k < T; ++t)
            g_k += dr[t].dot(dr[t + k]);
        g_k /= static_cast<double>(T);
        const double Ck = g_k / var0;
        if (Ck <= 0.0)
            break;  // initial-positive-sequence window
        sumC += Ck;
    }
    double g = 1.0 + 2.0 * sumC;
    if (g < 1.0) g = 1.0;  // statistical inefficiency cannot be < 1
    s.statineff = g;
    double nEff = static_cast<double>(T) / g;
    nEff = std::clamp(nEff, 1.0, static_cast<double>(T));
    s.nEff = nEff;
    s.valid = true;
    return s;
}

// ===========================================================================
// Anisotropic Gaussian bandwidth matrix H (kernel COVARIANCE).
// ===========================================================================
// Multivariate Silverman/Scott rule of thumb in 3-D:
//   H = (4/(d+2))^(2/(d+4)) * n^(-2/(d+4)) * Sigma     with d = 3
//     = (4/5)^(2/7) * n_eff^(-2/7) * Sigma
// then each eigenvalue is floored to (floorFactor*voxel)^2 so a thin/degenerate
// principal axis is rendered at grid resolution rather than aliasing. The floor
// depends on the voxel size, so the orchestrator recomputes H after any grid
// coarsening (codex: floor-after-coarsen ordering).
inline model::Mat3 bandwidthMatrix(const model::Mat3& covariance,
                                   double nEff, double voxel,
                                   double floorVoxelFactor) {
    const double n = std::max(1.0, nEff);
    const double scale = std::pow(4.0 / 5.0, 2.0 / 7.0) * std::pow(n, -2.0 / 7.0);
    model::Mat3 H = scale * covariance;
    // Symmetrise against round-off, then eigen-floor.
    H = 0.5 * (H + H.transpose());
    Eigen::SelfAdjointEigenSolver<model::Mat3> es(H);
    model::Vec3 lambda = es.eigenvalues();          // ascending, may have tiny <0
    const model::Mat3 V = es.eigenvectors();
    const double floor = (floorVoxelFactor * voxel) * (floorVoxelFactor * voxel);
    for (int i = 0; i < 3; ++i)
        lambda(i) = std::max(lambda(i), floor);
    return V * lambda.asDiagonal() * V.transpose();
}

// Largest eigenvalue of a symmetric PSD 3x3 (for the splat cutoff / grid pad).
inline double maxEigenvalue(const model::Mat3& M) {
    Eigen::SelfAdjointEigenSolver<model::Mat3> es(M);
    return es.eigenvalues()(2);  // ascending -> last is largest
}

// ===========================================================================
// Voxel grid (point-sampled). dims are POINT counts per axis; the scalar field
// lives on grid POINTS (matching the vtkImageData point-data contour
// downstream), x-fastest linear order.
// ===========================================================================
struct GridSpec {
    model::Vec3 origin = model::Vec3::Zero();  // world coord of point (0,0,0)
    double spacing = 1.0;                       // isotropic voxel edge (A)
    std::array<int, 3> dims = {0, 0, 0};        // POINT counts per axis

    long pointCount() const {
        return static_cast<long>(dims[0]) * dims[1] * dims[2];
    }
    long index(int i, int j, int k) const {
        return static_cast<long>(i) + static_cast<long>(dims[0]) *
               (static_cast<long>(j) + static_cast<long>(dims[1]) * k);
    }
    model::Vec3 pointAt(int i, int j, int k) const {
        return origin + spacing * model::Vec3(i, j, k);
    }
    double voxelVolume() const { return spacing * spacing * spacing; }
};

// Plan a grid covering the position bbox + marginSigma*sqrt(lambda_max(H)),
// at voxelTarget spacing, coarsening the (isotropic) spacing if any axis would
// exceed maxDim. origin is the padded min corner.
inline GridSpec planGrid(const std::vector<model::Vec3>& pos, const model::Mat3& H,
                         double voxelTarget, int maxDim, double marginSigma) {
    GridSpec g;
    if (pos.empty()) return g;

    model::Vec3 lo = pos.front(), hi = pos.front();
    for (const auto& p : pos) {
        lo = lo.cwiseMin(p);
        hi = hi.cwiseMax(p);
    }
    const double margin = marginSigma * std::sqrt(std::max(0.0, maxEigenvalue(H)));
    lo.array() -= margin;
    hi.array() += margin;
    const model::Vec3 extent = hi - lo;

    // Isotropic spacing: start at target, coarsen so the longest axis fits in
    // maxDim points (maxDim-1 intervals).
    double spacing = voxelTarget;
    const double maxExtent = extent.maxCoeff();
    const double minSpacingForCap = maxExtent / static_cast<double>(std::max(1, maxDim - 1));
    if (minSpacingForCap > spacing)
        spacing = minSpacingForCap;

    g.origin = lo;
    g.spacing = spacing;
    for (int a = 0; a < 3; ++a) {
        int n = static_cast<int>(std::ceil(extent(a) / spacing)) + 1;
        n = std::clamp(n, 2, maxDim);
        g.dims[a] = n;
    }
    return g;
}

// ===========================================================================
// Density field — sum of anisotropic Gaussians, normalised so sum(values)*dV=1.
// ===========================================================================
struct DensityField {
    GridSpec grid;
    std::vector<double> values;  // size = grid.pointCount(), x-fastest
    bool valid = false;
};

inline DensityField evaluateDensity(const std::vector<model::Vec3>& pos,
                                    const model::Mat3& H, const GridSpec& grid) {
    DensityField f;
    f.grid = grid;
    const long N = grid.pointCount();
    if (N <= 0 || pos.empty()) return f;
    f.values.assign(static_cast<std::size_t>(N), 0.0);

    const model::Mat3 Hinv = H.inverse();
    const double cutoff = 3.5 * std::sqrt(std::max(1e-12, maxEigenvalue(H)));
    const double spacing = grid.spacing;
    const model::Vec3& origin = grid.origin;

    // Per-frame: splat over the grid points within the cutoff box of that
    // frame's position. The per-frame Gaussian normalisation constant is the
    // same for every frame (shared H), so it cancels under the final
    // normalisation and is omitted here.
    for (const auto& p : pos) {
        std::array<int, 3> ilo, ihi;
        for (int a = 0; a < 3; ++a) {
            const double rel_lo = (p(a) - cutoff - origin(a)) / spacing;
            const double rel_hi = (p(a) + cutoff - origin(a)) / spacing;
            ilo[a] = std::clamp(static_cast<int>(std::floor(rel_lo)), 0, grid.dims[a] - 1);
            ihi[a] = std::clamp(static_cast<int>(std::ceil(rel_hi)), 0, grid.dims[a] - 1);
        }
        for (int k = ilo[2]; k <= ihi[2]; ++k)
            for (int j = ilo[1]; j <= ihi[1]; ++j)
                for (int i = ilo[0]; i <= ihi[0]; ++i) {
                    const model::Vec3 d = grid.pointAt(i, j, k) - p;
                    const double m = d.dot(Hinv * d);  // Mahalanobis^2
                    f.values[static_cast<std::size_t>(grid.index(i, j, k))] +=
                        std::exp(-0.5 * m);
                }
    }

    double total = 0.0;
    for (double v : f.values) total += v;
    total *= grid.voxelVolume();
    if (total <= 0.0) return f;  // valid stays false
    const double inv = 1.0 / total;
    for (double& v : f.values) v *= inv;  // now sum(values)*dV == 1
    f.valid = true;
    return f;
}

// ===========================================================================
// Highest-density-region levels.
// ===========================================================================
struct ShellLevel {
    double fraction = 0.0;          // requested mass fraction (e.g. 0.9)
    double isoValue = 0.0;          // c_X : contour at this density
    double includedMass = 0.0;      // mass(rho >= c_X)  — bracket UPPER (>= frac)
    double strictlyAboveMass = 0.0; // mass(rho >  c_X)  — bracket LOWER (<  frac)
    int components = 0;             // disconnected components (filled downstream)
    bool valid = false;
};

// For each requested fraction, find the density level c_X whose superlevel set
// {rho >= c_X} holds at least that fraction of the mass, and report the honest
// bracket mass(rho>c_X) < frac <= mass(rho>=c_X). Ties at c_X (common on a
// symmetric Gaussian lattice) are handled by summing the whole tie-shell.
inline std::vector<ShellLevel> hdrLevels(const DensityField& field,
                                         const std::vector<double>& fractions) {
    std::vector<ShellLevel> out;
    if (!field.valid || field.values.empty()) return out;
    const double dV = field.grid.voxelVolume();

    std::vector<double> sorted = field.values;
    std::sort(sorted.begin(), sorted.end(), std::greater<double>());
    double total = 0.0;
    for (double v : sorted) total += v;
    total *= dV;
    if (total <= 0.0) return out;

    for (double frac : fractions) {
        ShellLevel lvl;
        lvl.fraction = frac;
        const double target = frac * total;
        double cum = 0.0;
        double cX = sorted.back();  // default to the minimum (encloses all)
        for (double v : sorted) {
            cum += v * dV;
            if (cum >= target) { cX = v; break; }
        }
        // Re-sum including the full tie-shell at cX (codex: the crossing can be
        // a whole shell of equal-density points, not one voxel).
        double inc = 0.0, above = 0.0;
        for (double v : field.values) {
            if (v >= cX) inc += v;
            if (v > cX)  above += v;
        }
        lvl.isoValue = cX;
        lvl.includedMass = inc * dV;
        lvl.strictlyAboveMass = above * dV;
        lvl.valid = cX > 0.0;
        out.push_back(lvl);
    }
    return out;
}

// ===========================================================================
// Orchestrator.
// ===========================================================================
struct OccupancyResult {
    bool computed = false;       // false => skipped (see note)
    std::string note;            // skip reason / annotation
    MotionStats stats;
    model::Mat3 bandwidth = model::Mat3::Zero();  // final kernel H, floored at the FINAL voxel
    DensityField field;
    std::vector<ShellLevel> shells;  // ascending fraction (inner first)
};

inline OccupancyResult computeOccupancy(const std::vector<model::Vec3>& pos,
                                        const OccupancyConfig& cfg = {}) {
    OccupancyResult r;
    if (pos.size() < cfg.minFrames) {
        r.note = "too few frames";
        return r;
    }
    r.stats = computeMotionStats(pos);
    if (!r.stats.valid) {
        r.note = "motion stats invalid";
        return r;
    }
    if (r.stats.degenerate) {
        r.note = "degenerate: atom never moves (zero variance)";
        return r;
    }
    if (r.stats.rmsf < cfg.rigidRmsfFactor * cfg.voxelTarget) {
        r.note = "rigid: RMSF below resolution threshold";
        return r;
    }

    // H with the target voxel, plan the grid (may coarsen), then RE-floor H with
    // the FINAL spacing (the floor is voxel-tied). The data-driven eigenvalues
    // dominate the floor for any non-rigid atom, so the grid extent is stable.
    const model::Mat3 H0 = bandwidthMatrix(r.stats.covariance, r.stats.nEff,
                                           cfg.voxelTarget, cfg.floorVoxelFactor);
    GridSpec grid = planGrid(pos, H0, cfg.voxelTarget, cfg.maxDim, cfg.marginSigma);
    const model::Mat3 H = bandwidthMatrix(r.stats.covariance, r.stats.nEff,
                                          grid.spacing, cfg.floorVoxelFactor);

    r.bandwidth = H;
    r.field = evaluateDensity(pos, H, grid);
    if (!r.field.valid) {
        r.note = "density evaluation failed (zero mass)";
        return r;
    }
    r.shells = hdrLevels(r.field, cfg.massFractions);
    r.computed = true;
    if (grid.spacing > cfg.voxelTarget * 1.0001)
        r.note = "voxel coarsened to fit grid cap";
    return r;
}

}  // namespace h5reader::math

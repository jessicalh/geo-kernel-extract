// CsaShape -- principal-axis-system (PAS) / Haeberlen decomposition of a
// chemical-shielding-anisotropy (CSA) tensor, for the viewer's standalone
// "look at a DFT that hasn't had a rediscover pass yet" path.
//
// Pure Eigen + the shared model Vec3/Mat3 types; no Qt, no VTK -- by design,
// so it is headless-testable from tests/scene_math_tests exactly like
// TensorGlyphMath.h and FitTargetMath.h.
//
// The ALGORITHM is grabbed faithfully from AnalysisAtom.cpp::foldCsa (the
// rediscover engine on the nmr-shielding-analysisatom worktree) so the numbers
// match the extractor -- no functionality lost. But the SHAPE is ours: the
// engine folds this per step over a whole trajectory into an ephemeral atom
// object; the viewer wants one atom, one frame, on demand and interactively,
// so this is a clean pure function, NOT a port of the engine's per-step-series
// / fold machinery. The extractor's way is not the viewer's way. Computed live
// on a loaded DFT tensor, a bare extractor run -- one of the thousands an
// advisor or curious person might browse -- shows an honest, correctly-
// oriented shielding tensor with real shape numbers, no rediscover pass
// required. Staked into src/model so the extractor's AnalysisAtom could later
// call this same code in place of its inline copy when both sides stabilise.
//
// Conventions -- match the engine EXACTLY (a silent reorder corrupts every
// downstream component relationship):
//   - Input is the raw 3x3 lab-frame shielding tensor (ppm). Only the
//     symmetric part carries the CSA shape, so it is symmetrised first; the
//     antisymmetric (T1) part is intentionally dropped here.
//   - principal_values are ASCENDING: sigma11 <= sigma22 <= sigma33, and
//     pas_axes' columns follow that order (SelfAdjointEigenSolver order).
//   - haeberlen_values / haeberlen_axes are REORDERED to the Haeberlen
//     convention |sigma_zz - iso| >= |sigma_xx - iso| >= |sigma_yy - iso|.
//   - span  = sigma33 - sigma11.
//   - skew  = 3 (sigma22 - iso) / span   (sigma22 is the MIDDLE ascending value).
//   - eta in [0,1] = (sigma_yy - sigma_xx) / (sigma_zz - iso).
//   - prevPasAxes (optional) carries cross-frame sign/handedness continuity so
//     the glyph axes do not flip frame-to-frame on a tiny molecular twitch.
//     Thread the returned pas_axes back in as prevPasAxes for the next frame.

#pragma once

#include "Types.h"  // Vec3, Mat3 (Eigen)

#include <Eigen/Eigenvalues>

#include <array>
#include <cmath>
#include <limits>
#include <optional>

namespace h5reader::model {

struct CsaShape {
    bool valid = false;
    // Ascending eigenvalues sigma11 <= sigma22 <= sigma33 (ppm); columns of
    // pas_axes follow this same order (lab-frame director per column).
    Vec3 principal_values = Vec3::Constant(std::numeric_limits<double>::quiet_NaN());
    Mat3 pas_axes = Mat3::Constant(std::numeric_limits<double>::quiet_NaN());
    // Haeberlen-ordered values/axes (zz farthest from iso, then xx, then yy).
    Vec3 haeberlen_values = Vec3::Constant(std::numeric_limits<double>::quiet_NaN());
    Mat3 haeberlen_axes = Mat3::Constant(std::numeric_limits<double>::quiet_NaN());
    double sigma_iso = std::numeric_limits<double>::quiet_NaN();
    double eta = std::numeric_limits<double>::quiet_NaN();    // asymmetry, [0,1]
    double span = std::numeric_limits<double>::quiet_NaN();   // sigma33 - sigma11
    double skew = std::numeric_limits<double>::quiet_NaN();   // 3 (sigma22 - iso)/span
};

// Decompose a raw lab-frame shielding tensor into its PAS + Haeberlen shape.
// Returns an invalid CsaShape (valid == false) on any degeneracy: a failed
// eigensolve, a non-finite tensor, a near-isotropic tensor (|span| ~ 0), or a
// vanishing eta denominator. prevPasAxes, when supplied, enforces sign and
// handedness continuity with the previous frame's axes.
inline CsaShape ComputeCsaShape(const Mat3& rawTensor,
                                const Mat3* prevPasAxes = nullptr) {
    CsaShape out;

    const Mat3 sym = 0.5 * (rawTensor + rawTensor.transpose());
    Eigen::SelfAdjointEigenSolver<Mat3> solver(sym);
    if (solver.info() != Eigen::Success)
        return out;

    Vec3 principal = solver.eigenvalues();    // ascending
    Mat3 pas = solver.eigenvectors();         // columns follow `principal`
    if (!principal.allFinite() || !pas.allFinite())
        return out;

    // Cross-frame sign continuity: flip any column that points opposite to the
    // previous frame's corresponding axis.
    if (prevPasAxes) {
        for (int c = 0; c < 3; ++c) {
            if (prevPasAxes->col(c).dot(pas.col(c)) < 0.0)
                pas.col(c) *= -1.0;
        }
    }
    // Keep a right-handed frame. If the sign flips left it left-handed, flip
    // the least-continuous column (or the last column when there is no prior).
    if (pas.determinant() < 0.0) {
        if (prevPasAxes) {
            int flipCol = 0;
            double smallest = std::numeric_limits<double>::infinity();
            for (int c = 0; c < 3; ++c) {
                const double continuity = std::abs(prevPasAxes->col(c).dot(pas.col(c)));
                if (continuity < smallest) {
                    smallest = continuity;
                    flipCol = c;
                }
            }
            pas.col(flipCol) *= -1.0;
        } else {
            pas.col(2) *= -1.0;
        }
    }

    const double s11 = principal[0];
    const double s22 = principal[1];
    const double s33 = principal[2];
    const double iso = (s11 + s22 + s33) / 3.0;
    const double span = s33 - s11;
    if (!(std::abs(span) > 1e-12))
        return out;

    // Haeberlen ordering: zz is the principal value FARTHEST from isotropic;
    // of the remaining two, xx is the farther-from-iso, yy the closer.
    int zz = 0;
    double farthest = -1.0;
    for (int i = 0; i < 3; ++i) {
        const double dist = std::abs(principal[i] - iso);
        if (dist > farthest) {
            farthest = dist;
            zz = i;
        }
    }
    std::array<int, 2> rem{};
    int k = 0;
    for (int i = 0; i < 3; ++i)
        if (i != zz) rem[static_cast<std::size_t>(k++)] = i;
    int xx = rem[0];
    int yy = rem[1];
    if (std::abs(principal[xx] - iso) < std::abs(principal[yy] - iso))
        std::swap(xx, yy);

    const double denomEta = principal[zz] - iso;
    if (!(std::abs(denomEta) > 1e-12))
        return out;
    double eta = (principal[yy] - principal[xx]) / denomEta;
    if (eta < 0.0 && eta > -1e-12) eta = 0.0;
    if (eta > 1.0 && eta < 1.0 + 1e-12) eta = 1.0;
    if (!(eta >= 0.0 && eta <= 1.0))
        return out;

    out.principal_values = principal;
    out.pas_axes = pas;
    out.sigma_iso = iso;
    out.span = span;
    out.skew = 3.0 * (s22 - iso) / span;
    out.eta = eta;
    out.haeberlen_values = Vec3(principal[xx], principal[yy], principal[zz]);
    out.haeberlen_axes.col(0) = pas.col(xx);
    out.haeberlen_axes.col(1) = pas.col(yy);
    out.haeberlen_axes.col(2) = pas.col(zz);
    out.valid = std::isfinite(out.sigma_iso) && std::isfinite(out.eta)
                && std::isfinite(out.span) && std::isfinite(out.skew)
                && out.haeberlen_axes.allFinite();
    return out;
}

}  // namespace h5reader::model

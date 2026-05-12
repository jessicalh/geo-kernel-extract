// Structure tests for LarsenHBondGrid: trilinear lookup with periodic
// ρ wrap, out-of-range handling, archive selection from (donor,
// acceptor) class, FP tolerance at axis bounds, wrapped ρ in diagnostics,
// validity mask propagation, and the canonical donor frame builder.
//
// Fixture skips if RuntimeEnvironment::LarsenHBondGridDir() is empty
// (test machine without the dense grids configured).

#include "TestEnvironment.h"
#include <gtest/gtest.h>

#include "LarsenHBondGrid.h"
#include "RuntimeEnvironment.h"
#include "Session.h"
#include "errors.h"

#include <cmath>
#include <filesystem>
#include <string>

namespace fs = std::filesystem;

using namespace nmr;


class LarsenHBondGridTest : public ::testing::Test {
protected:
    void SetUp() override {
        const std::string& dir = RuntimeEnvironment::LarsenHBondGridDir();
        if (dir.empty()) {
            GTEST_SKIP() << "RuntimeEnvironment::LarsenHBondGridDir() empty "
                            "(no [paths].larsen_hbond_grids in TOML).";
        }
        if (!fs::exists(dir)) {
            GTEST_SKIP() << "data_dir does not exist: " << dir;
        }
        // Confirm the 6 dense.h5 files exist
        const char* stems[6] = {
            "NMANMA", "NMACOH", "NMACOO",
            "ALANMA", "ALACOH", "ALACOO",
        };
        for (auto stem : stems) {
            fs::path p = fs::path(dir) / (std::string(stem) + "_dense.h5");
            if (!fs::exists(p)) {
                GTEST_SKIP() << "missing grid file: " << p.string();
            }
        }
        ASSERT_EQ(session.LoadLarsenHBondGrid(), kOk)
            << session.LastError();
        ASSERT_TRUE(session.HasLarsenHBondGrid());
        grid = session.LarsenHBondGridPtr();
        ASSERT_NE(grid, nullptr);
    }

    Session session;
    const LarsenHBondGrid* grid = nullptr;
};


// All 6 archive-mapping pairs should produce a hit at a centrally-
// located geometry (r = 2.5 Å, θ = 150°, ρ = 0°) — well inside every
// grid's bounds.
TEST_F(LarsenHBondGridTest, AllArchiveMappingsHit) {
    struct Case {
        HBondDonorClass donor;
        HBondAcceptorClass acceptor;
        const char* label;
    };
    const Case cases[] = {
        {HBondDonorClass::AmideHydrogen, HBondAcceptorClass::BackboneCarbonyl,  "NMA→NMA"},
        {HBondDonorClass::AmideHydrogen, HBondAcceptorClass::SidechainCarbonyl, "NMA→sidechain(NMA approx)"},
        {HBondDonorClass::AmideHydrogen, HBondAcceptorClass::HydroxylOxygen,    "NMA→COH"},
        {HBondDonorClass::AmideHydrogen, HBondAcceptorClass::CarboxylateOxygen, "NMA→COO"},
        {HBondDonorClass::AlphaHydrogen, HBondAcceptorClass::BackboneCarbonyl,  "ALA→NMA"},
        {HBondDonorClass::AlphaHydrogen, HBondAcceptorClass::SidechainCarbonyl, "ALA→sidechain(NMA approx)"},
        {HBondDonorClass::AlphaHydrogen, HBondAcceptorClass::HydroxylOxygen,    "ALA→COH"},
        {HBondDonorClass::AlphaHydrogen, HBondAcceptorClass::CarboxylateOxygen, "ALA→COO"},
    };
    LarsenHBondGeometry geom{2.5, 150.0, 0.0};
    for (const auto& c : cases) {
        auto rec = grid->QueryNearest(c.donor, c.acceptor, geom);
        EXPECT_TRUE(rec.IsHit()) << "miss on " << c.label;
        EXPECT_TRUE(std::isfinite(rec.donor_N.norm()))
            << "non-finite donor_N on " << c.label;
        EXPECT_TRUE(std::isfinite(rec.donor_HA.norm()))
            << "non-finite donor_HA on " << c.label;
    }
}


// ALA donor archives have donor_CB; NMA donor archives don't.
TEST_F(LarsenHBondGridTest, DonorCBPresenceMatchesArchive) {
    LarsenHBondGeometry geom{2.5, 150.0, 0.0};
    auto ala_rec = grid->QueryNearest(
        HBondDonorClass::AlphaHydrogen,
        HBondAcceptorClass::BackboneCarbonyl,
        geom);
    EXPECT_TRUE(ala_rec.has_donor_CB)
        << "ALA donor archive must have CB readout";

    auto nma_rec = grid->QueryNearest(
        HBondDonorClass::AmideHydrogen,
        HBondAcceptorClass::BackboneCarbonyl,
        geom);
    EXPECT_FALSE(nma_rec.has_donor_CB)
        << "NMA donor archive must not have CB readout";
}


// Backbone-carbonyl acceptor archives have acceptor readouts; HOMe
// and acetate acceptor archives do not (Larsen 2015 does not define
// 2° terms for those acceptors).
TEST_F(LarsenHBondGridTest, AcceptorReadoutsOnlyForBackboneCarbonyl) {
    LarsenHBondGeometry geom{2.5, 150.0, 0.0};
    auto nma_acceptor_rec = grid->QueryNearest(
        HBondDonorClass::AlphaHydrogen,
        HBondAcceptorClass::BackboneCarbonyl,
        geom);
    EXPECT_TRUE(nma_acceptor_rec.has_acceptor_readouts);

    auto hydroxyl_rec = grid->QueryNearest(
        HBondDonorClass::AlphaHydrogen,
        HBondAcceptorClass::HydroxylOxygen,
        geom);
    EXPECT_FALSE(hydroxyl_rec.has_acceptor_readouts);

    auto carboxylate_rec = grid->QueryNearest(
        HBondDonorClass::AlphaHydrogen,
        HBondAcceptorClass::CarboxylateOxygen,
        geom);
    EXPECT_FALSE(carboxylate_rec.has_acceptor_readouts);
}


// Out-of-range r returns is_hit=false.
TEST_F(LarsenHBondGridTest, OutOfRangeR) {
    // r = 10 Å is well above every archive's r-axis maximum (4.0 Å).
    auto rec = grid->QueryNearest(
        HBondDonorClass::AlphaHydrogen,
        HBondAcceptorClass::BackboneCarbonyl,
        LarsenHBondGeometry{10.0, 150.0, 0.0});
    EXPECT_FALSE(rec.IsHit());
    // r = 0.5 Å is below every archive's minimum (1.5 Å NMA donor /
    // 1.6 Å ALA donor on the dense-grid axes).
    auto rec2 = grid->QueryNearest(
        HBondDonorClass::AlphaHydrogen,
        HBondAcceptorClass::BackboneCarbonyl,
        LarsenHBondGeometry{0.5, 150.0, 0.0});
    EXPECT_FALSE(rec2.IsHit());
}


// Out-of-range θ returns is_hit=false.
TEST_F(LarsenHBondGridTest, OutOfRangeTheta) {
    auto rec = grid->QueryNearest(
        HBondDonorClass::AlphaHydrogen,
        HBondAcceptorClass::BackboneCarbonyl,
        LarsenHBondGeometry{2.5, 45.0, 0.0});
    EXPECT_FALSE(rec.IsHit());
}


// F5: small FP overshoot at the axis bound should still hit.
// theta = 180.0000000001 should NOT silently miss the bin.
TEST_F(LarsenHBondGridTest, AxisBoundFPTolerance) {
    auto rec_at_180 = grid->QueryNearest(
        HBondDonorClass::AlphaHydrogen,
        HBondAcceptorClass::BackboneCarbonyl,
        LarsenHBondGeometry{2.5, 180.0, 0.0});
    EXPECT_TRUE(rec_at_180.IsHit());

    auto rec_just_past = grid->QueryNearest(
        HBondDonorClass::AlphaHydrogen,
        HBondAcceptorClass::BackboneCarbonyl,
        LarsenHBondGeometry{2.5, 180.0000000001, 0.0});
    EXPECT_TRUE(rec_just_past.IsHit())
        << "axis bound tolerance must absorb FP round-trip overshoot";
}


// ρ wraps periodically: ρ = 185° should give the same result as
// ρ = -175° (modulo numerical precision).
TEST_F(LarsenHBondGridTest, RhoWrapsPeriodically) {
    auto rec_plus  = grid->QueryNearest(
        HBondDonorClass::AlphaHydrogen,
        HBondAcceptorClass::BackboneCarbonyl,
        LarsenHBondGeometry{2.5, 150.0, 185.0});
    auto rec_minus = grid->QueryNearest(
        HBondDonorClass::AlphaHydrogen,
        HBondAcceptorClass::BackboneCarbonyl,
        LarsenHBondGeometry{2.5, 150.0, -175.0});
    ASSERT_TRUE(rec_plus.IsHit());
    ASSERT_TRUE(rec_minus.IsHit());
    double err = (rec_plus.donor_HA - rec_minus.donor_HA).norm();
    EXPECT_LT(err, 1e-4)
        << "ρ=185° should equal ρ=-175° via periodic wrap; got err="
        << err << " on donor_HA";
}


// F6: wrapped ρ is reported in rec.rho_deg (NOT the raw input).
TEST_F(LarsenHBondGridTest, WrappedRhoReportedInRecord) {
    auto rec = grid->QueryNearest(
        HBondDonorClass::AlphaHydrogen,
        HBondAcceptorClass::BackboneCarbonyl,
        LarsenHBondGeometry{2.5, 150.0, 185.0});
    ASSERT_TRUE(rec.IsHit());
    EXPECT_NEAR(rec.rho_deg, -175.0, 1e-6)
        << "rec.rho_deg should be the canonical wrapped value (-175°), "
           "not the raw input (185°)";

    auto rec2 = grid->QueryNearest(
        HBondDonorClass::AlphaHydrogen,
        HBondAcceptorClass::BackboneCarbonyl,
        LarsenHBondGeometry{2.5, 150.0, -185.0});
    EXPECT_NEAR(rec2.rho_deg, 175.0, 1e-6)
        << "rec.rho_deg should canonicalise -185° to +175°";
}


// Interpolated value at midpoint between two grid cells should differ
// from both endpoints by a finite, non-zero amount when the function
// varies.
TEST_F(LarsenHBondGridTest, InterpolationProducesFiniteOutput) {
    auto rec = grid->QueryNearest(
        HBondDonorClass::AlphaHydrogen,
        HBondAcceptorClass::BackboneCarbonyl,
        LarsenHBondGeometry{2.1, 125.0, 37.5});
    ASSERT_TRUE(rec.IsHit());
    EXPECT_TRUE(std::isfinite(rec.donor_HA.norm()));
    EXPECT_TRUE(std::isfinite(rec.donor_N.norm()));
    EXPECT_TRUE(std::isfinite(rec.donor_CA.norm()));
}


// any_corner_imputed propagates from the loaded mask. ALA donor
// archives have ~3-4% imputed nominal bins (Larsen-failed DFT calcs
// at close-approach geometries); scanning a coarse grid of queries
// should encounter at least one imputed cell.
TEST_F(LarsenHBondGridTest, ImputedCornersAreReported) {
    int n_imputed = 0;
    int n_total = 0;
    for (double r = 1.8; r <= 4.0; r += 0.2) {
        for (double th = 90.0; th <= 180.0; th += 10.0) {
            for (double rho = -180.0; rho < 180.0; rho += 30.0) {
                auto rec = grid->QueryNearest(
                    HBondDonorClass::AlphaHydrogen,
                    HBondAcceptorClass::BackboneCarbonyl,
                    LarsenHBondGeometry{r, th, rho});
                if (!rec.IsHit()) continue;
                ++n_total;
                if (rec.any_corner_imputed) ++n_imputed;
            }
        }
    }
    EXPECT_GT(n_total, 100) << "expected many hits in coarse scan";
    EXPECT_GT(n_imputed, 0)
        << "ALA donor archive should have some imputed cells; "
           "none reported across " << n_total << " queries — "
           "validity_mask wiring may be broken.";
}


// At a well-supported geometry, the mask should report no imputed
// corners. Picks the dead center of the well-populated theta=150°
// region (NMA donor, all rho values dense).
TEST_F(LarsenHBondGridTest, CentralGeometryHasNoImputedCorners) {
    auto rec = grid->QueryNearest(
        HBondDonorClass::AmideHydrogen,
        HBondAcceptorClass::BackboneCarbonyl,
        LarsenHBondGeometry{2.0, 150.0, 0.0});
    ASSERT_TRUE(rec.IsHit());
    EXPECT_FALSE(rec.any_corner_imputed)
        << "central NMA geometry should not hit imputed cells";
}


// Sanity: ALA donor at the "tight near-linear H-bond" geometry should
// produce a non-zero Hα contribution (Larsen Δσ_HαB is biggest in this
// region). The parser smoke showed isotropic ~-2 ppm at this geometry.
TEST_F(LarsenHBondGridTest, TightHBondGeometryGivesNonzeroHα) {
    auto rec = grid->QueryNearest(
        HBondDonorClass::AlphaHydrogen,
        HBondAcceptorClass::BackboneCarbonyl,
        LarsenHBondGeometry{2.0, 180.0, 0.0});
    ASSERT_TRUE(rec.IsHit());
    double iso_ha = rec.donor_HA.trace() / 3.0;
    EXPECT_GT(std::abs(iso_ha), 0.5)
        << "expected |Δσ_HA isotropic| > 0.5 ppm at tight linear H-bond; "
           "got " << iso_ha << " ppm";
}


// F3: ComputeLarsenDonorFrame returns an orthonormal rotation matrix
// with z aligned to (anchor → H) and x in the (anchor, third, H) plane,
// pointing in the (H − third) direction projected orthogonal to z.
// Tests x AND y AND z direction to lock in the sign convention (a
// y-sign-flipped frame would still be orthonormal with det +1, but
// would silently swap two coords on every transformed tensor).
TEST(LarsenHBondGridHelpers, DonorFrameIsOrthonormalAndAxesCorrect) {
    Vec3 H(1.0, 0.0, 0.0);
    Vec3 anchor(0.0, 0.0, 0.0);    // → z axis goes anchor→H = +x
    Vec3 third(0.0, 1.0, 0.0);     // in xy-plane, off z-axis
    Mat3 R = ComputeLarsenDonorFrame(H, anchor, third);

    // R should be orthonormal: R * R.T = I.
    Mat3 RRT = R * R.transpose();
    EXPECT_NEAR((RRT - Mat3::Identity()).norm(), 0.0, 1e-10);

    // Determinant should be +1 (proper rotation).
    EXPECT_NEAR(R.determinant(), 1.0, 1e-10);

    // z-axis (row 2 of R) should point along (H - anchor) = +x_lab.
    Vec3 z_lab = R.row(2);
    EXPECT_NEAR(z_lab.x(), 1.0, 1e-10);
    EXPECT_NEAR(z_lab.y(), 0.0, 1e-10);
    EXPECT_NEAR(z_lab.z(), 0.0, 1e-10);

    // x-axis is the projection of (H − third) = (1,−1,0) onto the
    // plane orthogonal to z (which is the +x lab axis), normalized.
    // The component perpendicular to z is (0,−1,0), so x_lab = (0,−1,0).
    Vec3 x_lab = R.row(0);
    EXPECT_NEAR(x_lab.x(), 0.0, 1e-10);
    EXPECT_NEAR(x_lab.y(), -1.0, 1e-10);
    EXPECT_NEAR(x_lab.z(), 0.0, 1e-10);

    // y = z × x = (1,0,0) × (0,−1,0) = (0·0 − 0·−1, 0·0 − 1·0,
    //                                    1·−1 − 0·0) = (0, 0, −1).
    Vec3 y_lab = R.row(1);
    EXPECT_NEAR(y_lab.x(), 0.0, 1e-10);
    EXPECT_NEAR(y_lab.y(), 0.0, 1e-10);
    EXPECT_NEAR(y_lab.z(), -1.0, 1e-10);
}


// Degenerate ComputeLarsenDonorFrame inputs return identity (with a
// logged warning) rather than NaN-poisoning the result.
TEST(LarsenHBondGridHelpers, DonorFrameDegenerateReturnsIdentity) {
    Vec3 H(1.0, 0.0, 0.0);
    Vec3 anchor_same(1.0, 0.0, 0.0);   // coincident with H
    Vec3 third(0.0, 1.0, 0.0);
    Mat3 R1 = ComputeLarsenDonorFrame(H, anchor_same, third);
    EXPECT_NEAR((R1 - Mat3::Identity()).norm(), 0.0, 1e-12);

    Vec3 anchor(0.0, 0.0, 0.0);
    Vec3 third_same(1.0, 0.0, 0.0);    // coincident with H
    Mat3 R2 = ComputeLarsenDonorFrame(H, anchor, third_same);
    EXPECT_NEAR((R2 - Mat3::Identity()).norm(), 0.0, 1e-12);

    Vec3 third_collinear(2.0, 0.0, 0.0);  // on the anchor→H line
    Mat3 R3 = ComputeLarsenDonorFrame(H, anchor, third_collinear);
    EXPECT_NEAR((R3 - Mat3::Identity()).norm(), 0.0, 1e-12);
}


// F3 + geometry contract: ComputeLarsenHBondGeometry computes the
// expected angle on a synthetic geometry, and θ-range is correct.
TEST(LarsenHBondGridHelpers, GeometryComputesCorrectAngleRange) {
    // Bent geometry: H at origin, O at (1,0,0), C at (1,1,0) → angle
    // at O (between O→H and O→C) is 90°. Third off-plane for a
    // well-defined dihedral.
    Vec3 H(0.0, 0.0, 0.0);
    Vec3 O(1.0, 0.0, 0.0);
    Vec3 C_bent(1.0, 1.0, 0.0);
    Vec3 third_bent(2.0, 1.0, 1.0);
    auto g = ComputeLarsenHBondGeometry(H, O, C_bent, third_bent);
    EXPECT_NEAR(g.r_angstrom, 1.0, 1e-10);
    EXPECT_NEAR(g.theta_deg, 90.0, 1e-6);
}


// IUPAC dihedral sign is locked in by a non-degenerate analytic case.
// A future refactor of the dihedral sign (e.g. swapping atan2(y,x) →
// atan2(-y,x)) would silently re-flip the grid lookup back to Larsen's
// filename convention, querying the wrong ρ lobe. This test catches that.
TEST(LarsenHBondGridHelpers, DihedralSignIsIupac) {
    // H=(0,0,0), O=(1,0,0), C=(1,1,0), third=(2,2,1).
    //   b1 = O − H = (1,0,0)
    //   b2 = C − O = (0,1,0)
    //   b3 = third − C = (1,1,1)
    //   n1 = b1 × b2 = (0,0,1)
    //   n2 = b2 × b3 = (1,0,-1)
    //   m1 = n1 × b2̂ = (0,0,1) × (0,1,0) = (-1,0,0)
    //   x = n1·n2 = -1, y = m1·n2 = -1
    //   ρ = atan2(-1, -1) = -135°
    Vec3 H(0.0, 0.0, 0.0);
    Vec3 O(1.0, 0.0, 0.0);
    Vec3 C(1.0, 1.0, 0.0);
    Vec3 third(2.0, 2.0, 1.0);
    auto g = ComputeLarsenHBondGeometry(H, O, C, third);
    EXPECT_NEAR(g.rho_deg, -135.0, 1e-6)
        << "IUPAC dihedral convention should give -135° at this geometry";

    // Mirror through z-plane: third = (2,2,−1). ρ should flip sign.
    Vec3 third_mirror(2.0, 2.0, -1.0);
    auto g2 = ComputeLarsenHBondGeometry(H, O, C, third_mirror);
    EXPECT_NEAR(g2.rho_deg, +135.0, 1e-6)
        << "z-mirror should flip dihedral sign";
}


// RotateTensorToProteinLabFrame applies σ_lab = Rᵀ · σ_canonical · R.
// Identity rotation preserves σ. A 90° z-rotation permutes diagonal
// entries as expected.
TEST(LarsenHBondGridHelpers, RotateTensorToProteinLabFrameMath) {
    Mat3 sigma_canonical;
    sigma_canonical << 100.0,   0.0,  0.0,
                         0.0,  50.0,  0.0,
                         0.0,   0.0, 10.0;

    // Identity: σ_lab == σ_canonical.
    Mat3 sigma_lab_id = RotateTensorToProteinLabFrame(sigma_canonical,
                                                      Mat3::Identity());
    EXPECT_NEAR((sigma_lab_id - sigma_canonical).norm(), 0.0, 1e-12);

    // 90° rotation around z. R has rows [x; y; z] where x = (0,1,0)
    // (e_x_canonical = e_y_lab), y = (-1,0,0), z = (0,0,1).
    //   R = [[0, 1, 0], [-1, 0, 0], [0, 0, 1]]
    //   σ_lab = R.T · σ_canonical · R
    //         = [[ 50,  0,  0],
    //            [  0,100,  0],
    //            [  0,  0, 10]]
    // (in the lab frame, the principal axis with value 100 is now along
    //  +y_lab; the 50 is along +x_lab.)
    Mat3 R_z90;
    R_z90 <<  0.0,  1.0, 0.0,
             -1.0,  0.0, 0.0,
              0.0,  0.0, 1.0;
    Mat3 sigma_lab_z90 = RotateTensorToProteinLabFrame(sigma_canonical, R_z90);
    Mat3 expected;
    expected <<  50.0,   0.0,  0.0,
                  0.0, 100.0,  0.0,
                  0.0,   0.0, 10.0;
    EXPECT_NEAR((sigma_lab_z90 - expected).norm(), 0.0, 1e-10);
}

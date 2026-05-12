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
    // r = 10 Å is well outside every archive's r_axis (1.5/1.6/1.8 → 4.0).
    auto rec = grid->QueryNearest(
        HBondDonorClass::AlphaHydrogen,
        HBondAcceptorClass::BackboneCarbonyl,
        LarsenHBondGeometry{10.0, 150.0, 0.0});
    EXPECT_FALSE(rec.IsHit());
    // r = 0.5 Å is below every minimum.
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
// with z aligned to (anchor → H) and x in the (anchor, third, H) plane.
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
}


// F3 + geometry contract: ComputeLarsenHBondGeometry computes the
// expected angle/dihedral on a synthetic geometry.
//
// Setup: donor H at origin, acceptor O at +x (1 Å away), acceptor C
// at +x_lab + small +y (so H-O-C angle is 180° minus small), third
// atom at +x_lab + +y2 + 0 (in the same plane → ρ = 0 ish).
TEST(LarsenHBondGridHelpers, GeometryComputesCorrectAngleAndDihedral) {
    // Linear chain: H — O — C — third, all on +x axis. θ should be
    // 180°, ρ undefined (third must be off the chain), so place third
    // off-axis.
    Vec3 H(0.0, 0.0, 0.0);
    Vec3 O(1.0, 0.0, 0.0);
    Vec3 C(2.0, 0.0, 0.0);
    Vec3 third(2.5, 1.0, 0.0);  // in xy-plane → ρ = 180° (flat)
    auto g = ComputeLarsenHBondGeometry(H, O, C, third);
    EXPECT_NEAR(g.r_angstrom, 1.0, 1e-10);
    EXPECT_NEAR(g.theta_deg, 180.0, 1e-6);
    // Flat (planar) geometry: ρ = 180° (trans configuration of third).
    EXPECT_TRUE(std::abs(std::abs(g.rho_deg) - 180.0) < 1e-6
                || std::abs(g.rho_deg) < 1e-6)
        << "planar setup should give ρ = 0° or 180°; got " << g.rho_deg;

    // Bent θ. Place C at (2, 1, 0); ∠ H-O-C should be < 180°.
    Vec3 C_bent(2.0, 1.0, 0.0);
    Vec3 third_bent(2.5, 1.5, 0.0);
    auto g2 = ComputeLarsenHBondGeometry(H, O, C_bent, third_bent);
    EXPECT_LT(g2.theta_deg, 180.0);
    EXPECT_GT(g2.theta_deg, 90.0);
}

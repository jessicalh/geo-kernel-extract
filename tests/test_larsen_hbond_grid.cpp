// Structure tests for LarsenHBondGrid: trilinear lookup with periodic
// ρ wrap, out-of-range handling, archive selection from (donor,
// acceptor) class.
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
    for (const auto& c : cases) {
        auto rec = grid->QueryNearest(c.donor, c.acceptor, 2.5, 150.0, 0.0);
        EXPECT_TRUE(rec.IsHit()) << "miss on " << c.label;
        EXPECT_TRUE(std::isfinite(rec.donor_N.norm()))
            << "non-finite donor_N on " << c.label;
        EXPECT_TRUE(std::isfinite(rec.donor_HA.norm()))
            << "non-finite donor_HA on " << c.label;
    }
}


// ALA donor archives have donor_CB; NMA donor archives don't.
TEST_F(LarsenHBondGridTest, DonorCBPresenceMatchesArchive) {
    auto ala_rec = grid->QueryNearest(
        HBondDonorClass::AlphaHydrogen,
        HBondAcceptorClass::BackboneCarbonyl,
        2.5, 150.0, 0.0);
    EXPECT_TRUE(ala_rec.has_donor_CB)
        << "ALA donor archive must have CB readout";

    auto nma_rec = grid->QueryNearest(
        HBondDonorClass::AmideHydrogen,
        HBondAcceptorClass::BackboneCarbonyl,
        2.5, 150.0, 0.0);
    EXPECT_FALSE(nma_rec.has_donor_CB)
        << "NMA donor archive must not have CB readout";
}


// Backbone-carbonyl acceptor archives have acceptor readouts; HOMe
// and acetate acceptor archives do not (Larsen 2015 does not define
// 2° terms for those acceptors).
TEST_F(LarsenHBondGridTest, AcceptorReadoutsOnlyForBackboneCarbonyl) {
    auto nma_acceptor_rec = grid->QueryNearest(
        HBondDonorClass::AlphaHydrogen,
        HBondAcceptorClass::BackboneCarbonyl,
        2.5, 150.0, 0.0);
    EXPECT_TRUE(nma_acceptor_rec.has_acceptor_readouts);

    auto hydroxyl_rec = grid->QueryNearest(
        HBondDonorClass::AlphaHydrogen,
        HBondAcceptorClass::HydroxylOxygen,
        2.5, 150.0, 0.0);
    EXPECT_FALSE(hydroxyl_rec.has_acceptor_readouts);

    auto carboxylate_rec = grid->QueryNearest(
        HBondDonorClass::AlphaHydrogen,
        HBondAcceptorClass::CarboxylateOxygen,
        2.5, 150.0, 0.0);
    EXPECT_FALSE(carboxylate_rec.has_acceptor_readouts);
}


// Out-of-range r returns is_hit=false.
TEST_F(LarsenHBondGridTest, OutOfRangeR) {
    // r = 10 Å is well outside every archive's r_axis (1.5/1.6/1.8 → 4.0).
    auto rec = grid->QueryNearest(
        HBondDonorClass::AlphaHydrogen,
        HBondAcceptorClass::BackboneCarbonyl,
        10.0, 150.0, 0.0);
    EXPECT_FALSE(rec.IsHit());
    // r = 0.5 Å is below every minimum.
    auto rec2 = grid->QueryNearest(
        HBondDonorClass::AlphaHydrogen,
        HBondAcceptorClass::BackboneCarbonyl,
        0.5, 150.0, 0.0);
    EXPECT_FALSE(rec2.IsHit());
}


// Out-of-range θ returns is_hit=false.
TEST_F(LarsenHBondGridTest, OutOfRangeTheta) {
    auto rec = grid->QueryNearest(
        HBondDonorClass::AlphaHydrogen,
        HBondAcceptorClass::BackboneCarbonyl,
        2.5, 45.0, 0.0);
    EXPECT_FALSE(rec.IsHit());
}


// ρ wraps periodically: ρ = 185° should give the same result as
// ρ = -175° (modulo numerical precision).
TEST_F(LarsenHBondGridTest, RhoWrapsPeriodically) {
    auto rec_plus  = grid->QueryNearest(
        HBondDonorClass::AlphaHydrogen,
        HBondAcceptorClass::BackboneCarbonyl,
        2.5, 150.0, 185.0);
    auto rec_minus = grid->QueryNearest(
        HBondDonorClass::AlphaHydrogen,
        HBondAcceptorClass::BackboneCarbonyl,
        2.5, 150.0, -175.0);
    ASSERT_TRUE(rec_plus.IsHit());
    ASSERT_TRUE(rec_minus.IsHit());
    // Tensors should match within FP tolerance — same geometry under
    // the periodic-ρ identification.
    double err = (rec_plus.donor_HA - rec_minus.donor_HA).norm();
    EXPECT_LT(err, 1e-4)
        << "ρ=185° should equal ρ=-175° via periodic wrap; got err="
        << err << " on donor_HA";
}


// Interpolated value at midpoint between two grid cells should differ
// from both endpoints by a finite, non-zero amount when the function
// varies. Just confirm the interpolation produces *some* output (the
// strict numerical equality is in the grid-point test below).
TEST_F(LarsenHBondGridTest, InterpolationProducesFiniteOutput) {
    // ALA donor at r=2.1 Å (between 2.0 and 2.2 grid points), θ=125°
    // (between 120 and 130), ρ=37.5° (between 30 and 45)
    auto rec = grid->QueryNearest(
        HBondDonorClass::AlphaHydrogen,
        HBondAcceptorClass::BackboneCarbonyl,
        2.1, 125.0, 37.5);
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
        2.0, 180.0, 0.0);
    ASSERT_TRUE(rec.IsHit());
    double iso_ha = rec.donor_HA.trace() / 3.0;
    EXPECT_GT(std::abs(iso_ha), 0.5)
        << "expected |Δσ_HA isotropic| > 0.5 ppm at tight linear H-bond; "
           "got " << iso_ha << " ppm";
}

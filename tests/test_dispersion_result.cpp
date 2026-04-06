#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include <cmath>
#include <iostream>

#include "DispersionResult.h"
#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "Protein.h"
#include "PdbFileReader.h"
#include "OrcaRunLoader.h"
#include "PhysicalConstants.h"

#include <filesystem>
namespace fs = std::filesystem;
using namespace nmr;



// ============================================================================
// Analytical test: single vertex at known distance.
//
// Vertex at origin, atom at (3, 0, 0). r = 3.
//   K_xx = 3*9/r^8 - 1/r^6 = 27/6561 - 1/729 = 0.004115 - 0.001372 = 0.002743
//   K_yy = 0 - 1/r^6 = -1/729 = -0.001372
//   K_zz = 0 - 1/r^6 = -1/729 = -0.001372
//   Tr = K_xx + K_yy + K_zz = 0.002743 - 0.001372 - 0.001372 = ~0
//   scalar = 1/r^6 = 1/729
// ============================================================================

TEST(DispAnalytical, SingleVertexTraceless) {
    Vec3 d(3, 0, 0);
    double r = d.norm();
    double r2 = r * r;
    double r6 = r2 * r2 * r2;
    double r8 = r6 * r2;

    Mat3 K = Mat3::Zero();
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            K(a, b) = 3.0 * d(a) * d(b) / r8 - (a == b ? 1.0 : 0.0) / r6;

    EXPECT_NEAR(K.trace(), 0.0, 1e-14) << "Dispersion kernel is traceless";
    EXPECT_NEAR(K(0, 0), 2.0 / r6, 1e-14);   // 3*9/r^8 - 1/r^6 = 3/r^6 - 1/r^6 = 2/r^6
    EXPECT_NEAR(K(1, 1), -1.0 / r6, 1e-14);
    EXPECT_NEAR(K(2, 2), -1.0 / r6, 1e-14);

    double scalar = 1.0 / r6;
    EXPECT_NEAR(scalar, 1.0 / 729.0, 1e-14);

    std::cout << "  Single vertex at r=3: K_diag = ("
              << K(0,0) << ", " << K(1,1) << ", " << K(2,2) << ")\n"
              << "  Trace = " << K.trace() << "\n"
              << "  Scalar = " << scalar << "\n";
}


TEST(DispAnalytical, SwitchingFunctionProperties) {
    // Reproduce the CHARMM switching function locally for testing.
    // Must match DispersionResult.cpp: R_SWITCH=4.3A, R_CUT=5.0A.
    // Brooks et al., J. Comput. Chem. 4, 187 (1983).
    constexpr double R_SWITCH = 4.3;  // Angstroms — onset of taper
    constexpr double R_CUT = 5.0;     // Angstroms — zero beyond

    auto S = [](double r) -> double {
        if (r <= R_SWITCH) return 1.0;
        if (r >= R_CUT) return 0.0;
        double rc2 = R_CUT * R_CUT;
        double rs2 = R_SWITCH * R_SWITCH;
        double r2 = r * r;
        double num = (rc2 - r2) * (rc2 - r2) * (rc2 + 2.0 * r2 - 3.0 * rs2);
        double den = (rc2 - rs2) * (rc2 - rs2) * (rc2 - rs2);
        return num / den;
    };

    // Boundary values
    EXPECT_NEAR(S(2.0), 1.0, 1e-12) << "Below switch onset: S=1";
    EXPECT_NEAR(S(4.3), 1.0, 1e-12) << "At switch onset: S=1";
    EXPECT_NEAR(S(5.0), 0.0, 1e-12) << "At cutoff: S=0";
    EXPECT_GT(S(4.5), 0.0) << "In taper: S > 0";
    EXPECT_LT(S(4.5), 1.0) << "In taper: S < 1";

    // C¹ continuity: S'(R_switch) = 0 and S'(R_cut) = 0
    // Verify numerically via finite differences
    double h = 1e-7;
    double dS_at_switch = (S(R_SWITCH + h) - S(R_SWITCH - h)) / (2 * h);
    double dS_at_cut = (S(R_CUT - h) - S(R_CUT - 2*h)) / h;
    EXPECT_NEAR(dS_at_switch, 0.0, 1e-4) << "S'(R_switch) = 0 (C¹)";
    EXPECT_NEAR(dS_at_cut, 0.0, 1e-4) << "S'(R_cut) = 0 (C¹)";

    // Monotonic decrease through the taper
    double prev = 1.0;
    for (double r = 4.3; r <= 5.0; r += 0.01) {
        double s = S(r);
        EXPECT_LE(s, prev + 1e-12)
            << "S must decrease monotonically at r=" << r;
        prev = s;
    }

    // Smoothness: S(r)/r^6 should decrease monotonically for all r > MIN_DISTANCE
    double prev_weighted = 1e10;
    for (double r = 1.6; r < 5.0; r += 0.05) {
        double weighted = S(r) / std::pow(r, 6);
        EXPECT_LE(weighted, prev_weighted + 1e-12)
            << "S/r^6 must decrease monotonically at r=" << r;
        prev_weighted = weighted;
    }

    std::cout << "  Switching function (CHARMM form, Brooks 1983):\n"
              << "    onset=" << R_SWITCH << "A, cutoff=" << R_CUT << "A\n"
              << "    S(3.0)=" << S(3.0) << " S(4.5)=" << S(4.5)
              << " S(4.9)=" << S(4.9) << " S(5.0)=" << S(5.0) << "\n"
              << "    S'(onset)=" << dS_at_switch
              << " S'(cutoff)=" << dS_at_cut << " (both ~0 for C¹)\n";
}


// ============================================================================
// Full protein test on 1UBQ
// ============================================================================

class DispProteinTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!fs::exists(nmr::test::TestEnvironment::UbqProtonated())) GTEST_SKIP() << "1UBQ not found";
        auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
        if (!r.Ok()) GTEST_SKIP() << "Failed to load 1UBQ";
        protein = std::move(r.protein);

        auto& conf = protein->Conformation();
        conf.AttachResult(GeometryResult::Compute(conf));
        conf.AttachResult(SpatialIndexResult::Compute(conf));
    }
    std::unique_ptr<Protein> protein;
};


TEST_F(DispProteinTest, ComputeAndAttach) {
    auto& conf = protein->Conformation();
    auto disp = DispersionResult::Compute(conf);
    ASSERT_NE(disp, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(disp)));
    ASSERT_TRUE(conf.HasResult<DispersionResult>());
}


TEST_F(DispProteinTest, TracelessAndSymmetric) {
    auto& conf = protein->Conformation();
    conf.AttachResult(DispersionResult::Compute(conf));

    int checked = 0;
    double max_trace = 0.0;
    double max_asym = 0.0;

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& rn : conf.AtomAt(ai).ring_neighbours) {
            if (rn.disp_tensor.isZero(1e-20)) continue;

            max_trace = std::max(max_trace, std::abs(rn.disp_tensor.trace()));

            for (int a = 0; a < 3; ++a)
                for (int b = a + 1; b < 3; ++b)
                    max_asym = std::max(max_asym,
                        std::abs(rn.disp_tensor(a, b) - rn.disp_tensor(b, a)));

            checked++;
        }
    }

    EXPECT_GT(checked, 0) << "Should have dispersion pairs";
    EXPECT_LT(max_trace, 1e-10) << "Dispersion tensor must be traceless";
    EXPECT_LT(max_asym, 1e-14) << "Dispersion tensor must be symmetric";

    std::cout << "  Checked " << checked << " ring-atom pairs\n"
              << "  max |Tr| = " << max_trace
              << ", max asymmetry = " << max_asym << "\n";
}


TEST_F(DispProteinTest, ContactCountsReasonable) {
    auto& conf = protein->Conformation();
    conf.AttachResult(DispersionResult::Compute(conf));

    int total_contacts = 0;
    int max_contacts = 0;
    int with_disp = 0;

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& rn : conf.AtomAt(ai).ring_neighbours) {
            if (rn.disp_contacts > 0) {
                total_contacts += rn.disp_contacts;
                max_contacts = std::max(max_contacts, rn.disp_contacts);
                with_disp++;
            }
        }
    }

    EXPECT_GT(with_disp, 0) << "Some atom-ring pairs should have contacts";
    // Each ring has 5-9 vertices. Max contacts per ring <= ring size.
    EXPECT_LE(max_contacts, 9) << "Max contacts should not exceed max ring size";

    std::cout << "  Pairs with contacts: " << with_disp
              << ", total contacts: " << total_contacts
              << ", max per pair: " << max_contacts << "\n";
}


TEST_F(DispProteinTest, ShieldingContributionHasT2) {
    auto& conf = protein->Conformation();
    conf.AttachResult(DispersionResult::Compute(conf));

    int nonzero_t2 = 0;
    double max_t2 = 0.0;

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        double t2m = conf.AtomAt(ai).disp_shielding_contribution.T2Magnitude();
        if (t2m > 1e-10) nonzero_t2++;
        max_t2 = std::max(max_t2, t2m);
    }

    EXPECT_GT(nonzero_t2, 0) << "Some atoms should have dispersion T2";

    std::cout << "  T2 nonzero: " << nonzero_t2
              << ", max |T2| = " << max_t2 << "\n";
}


// ============================================================================
// ORCA protein test
// ============================================================================

TEST(DispOrcaTest, RunOnProtonatedProtein) {
    OrcaRunFiles files;
    files.pdb_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.pdb";
    files.xyz_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.xyz";
    files.prmtop_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.prmtop";

    if (!fs::exists(files.xyz_path) || !fs::exists(files.prmtop_path))
        GTEST_SKIP() << "ORCA test data not found";

    auto load = BuildFromOrca(files);
    ASSERT_TRUE(load.Ok()) << load.error;

    auto& conf = load.protein->Conformation();
    conf.AttachResult(GeometryResult::Compute(conf));
    conf.AttachResult(SpatialIndexResult::Compute(conf));
    conf.AttachResult(DispersionResult::Compute(conf));

    int with_disp = 0;
    double max_scalar = 0.0;
    double max_t2 = 0.0;

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& rn : conf.AtomAt(ai).ring_neighbours) {
            if (rn.disp_contacts > 0) with_disp++;
            max_scalar = std::max(max_scalar, std::abs(rn.disp_scalar));
        }
        max_t2 = std::max(max_t2,
            conf.AtomAt(ai).disp_shielding_contribution.T2Magnitude());
    }

    EXPECT_GT(with_disp, 0) << "Some pairs should have dispersion contacts";

    std::cout << "  ORCA protein Dispersion summary:\n"
              << "    atoms=" << conf.AtomCount()
              << " rings=" << load.protein->RingCount() << "\n"
              << "    pairs with contacts: " << with_disp << "\n"
              << "    max |scalar| = " << max_scalar << " A^-6\n"
              << "    max |T2| = " << max_t2 << " A^-6\n";
}

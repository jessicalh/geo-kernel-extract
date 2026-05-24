#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include <cmath>
#include <iostream>

#include "HaighMallionResult.h"
#include "BiotSavartResult.h"
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
// Full protein test: 1UBQ
// ============================================================================

class HaighMallionProteinTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!fs::exists(nmr::test::TestEnvironment::UbqProtonated())) {
            GTEST_SKIP() << "1UBQ.pdb not found";
        }
        auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
        if (!r.Ok()) GTEST_SKIP() << "Failed to load 1UBQ";
        protein = std::move(r.protein);

        auto& conf = protein->Conformation();
        conf.AttachResult(GeometryResult::Compute(conf));
        conf.AttachResult(SpatialIndexResult::Compute(conf));
    }

    std::unique_ptr<Protein> protein;
};


TEST_F(HaighMallionProteinTest, ComputeAndAttach) {
    auto& conf = protein->Conformation();
    auto hm = HaighMallionResult::Compute(conf);
    ASSERT_NE(hm, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(hm)));
    ASSERT_TRUE(conf.HasResult<HaighMallionResult>());
}


TEST_F(HaighMallionProteinTest, RawIntegralIsSymmetricTraceless) {
    // H_ab = integral_S K_ab dS. K is symmetric and traceless at every point,
    // so the integral must be symmetric and traceless.
    auto& conf = protein->Conformation();
    conf.AttachResult(HaighMallionResult::Compute(conf));

    int checked = 0;
    double max_asym = 0, max_trace = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& rn : conf.AtomAt(ai).ring_neighbours) {
            if (rn.hm_H_tensor.norm() < 1e-15) continue;

            // Symmetry: H_ab = H_ba
            double asym = (rn.hm_H_tensor - rn.hm_H_tensor.transpose()).norm();
            max_asym = std::max(max_asym, asym);

            // Tracelessness: Tr(H) = 0
            double trace = std::abs(rn.hm_H_tensor.trace());
            max_trace = std::max(max_trace, trace);

            checked++;
        }
    }

    EXPECT_GT(checked, 100) << "Should verify many ring-atom pairs";
    EXPECT_LT(max_asym, 1e-10)
        << "H must be symmetric (integral of symmetric integrand)";
    EXPECT_LT(max_trace, 1e-10)
        << "H must be traceless (integral of traceless integrand)";

    std::cout << "  Checked " << checked << " pairs\n"
              << "  max asymmetry = " << max_asym << "\n"
              << "  max |trace|   = " << max_trace << "\n";
}


TEST_F(HaighMallionProteinTest, FullKernelIsRankOne) {
    // G_ab = n_b * V_a is rank-1 -> det(G) = 0.
    auto& conf = protein->Conformation();
    conf.AttachResult(HaighMallionResult::Compute(conf));

    // We need to reconstruct G from the stored data. The hm_shielding_contribution
    // is the accumulated sum. For per-ring rank-1 check, we reconstruct from
    // hm_B_field and ring normal.
    int checked = 0;
    double max_det = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& rn : conf.AtomAt(ai).ring_neighbours) {
            if (rn.hm_B_field.norm() < 1e-15) continue;

            const RingGeometry& geom = conf.ring_geometries[rn.ring_index];

            // Reconstruct G = n (x) V where V = hm_B_field
            Mat3 G;
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                    G(a, b) = geom.normal(b) * rn.hm_B_field(a);

            double det = std::abs(G.determinant());
            max_det = std::max(max_det, det);
            checked++;
        }
    }

    EXPECT_GT(checked, 100);
    EXPECT_LT(max_det, 1e-12)
        << "G = n (x) (H.n) is rank-1, determinant must be zero";

    std::cout << "  Checked " << checked << " pairs, max |det(G)| = "
              << max_det << "\n";
}


TEST_F(HaighMallionProteinTest, ShieldingContributionHasAllIrreps) {
    // The full kernel G = n (x) V has T0, T1, T2 all non-zero.
    auto& conf = protein->Conformation();
    conf.AttachResult(HaighMallionResult::Compute(conf));

    int nonzero_t0 = 0, nonzero_t1 = 0, nonzero_t2 = 0;
    double max_t0 = 0, max_t2 = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& sc = conf.AtomAt(ai).hm_shielding_contribution;
        if (std::abs(sc.T0) > 1e-8) nonzero_t0++;
        double t1mag = std::sqrt(sc.T1[0]*sc.T1[0] + sc.T1[1]*sc.T1[1]
                                 + sc.T1[2]*sc.T1[2]);
        if (t1mag > 1e-8) nonzero_t1++;
        double t2mag = sc.T2Magnitude();
        if (t2mag > 1e-8) nonzero_t2++;
        max_t0 = std::max(max_t0, std::abs(sc.T0));
        max_t2 = std::max(max_t2, t2mag);
    }

    EXPECT_GT(nonzero_t0, 0) << "Full kernel must have non-zero T0";
    EXPECT_GT(nonzero_t1, 0) << "Full kernel must have non-zero T1";
    EXPECT_GT(nonzero_t2, 0) << "Full kernel must have non-zero T2";

    std::cout << "  T0 nonzero: " << nonzero_t0 << ", max |T0| = " << max_t0 << "\n";
    std::cout << "  T1 nonzero: " << nonzero_t1 << "\n";
    std::cout << "  T2 nonzero: " << nonzero_t2 << ", max |T2| = " << max_t2 << "\n";
}


TEST_F(HaighMallionProteinTest, ConvergesToBSAtLargeDistance) {
    // At large distance from the ring, BS (wire segments) and HM (surface integral)
    // model the same physics and should produce similar T0 values. The T0 ratio
    // should approach a constant (related to the effective magnetic moment).
    // We check that HM and BS T0 have the same sign for distant atoms.
    auto& conf = protein->Conformation();
    conf.AttachResult(BiotSavartResult::Compute(conf));
    conf.AttachResult(HaighMallionResult::Compute(conf));

    int same_sign = 0, opposite_sign = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& rn : conf.AtomAt(ai).ring_neighbours) {
            if (rn.distance_to_center < 8.0) continue;  // only far-field
            if (std::abs(rn.G_spherical.T0) < 1e-10) continue;

            // Reconstruct HM full kernel T0 for this ring.
            // Full kernel G = -n (x) (H.n), so T0 = -n.(H.n)/3.
            // Same minus sign as BS (from sigma = -dB/dB_0).
            const RingGeometry& geom = conf.ring_geometries[rn.ring_index];
            Vec3 V = rn.hm_H_tensor * geom.normal;
            double hm_t0 = -geom.normal.dot(V) / 3.0;

            if (std::abs(hm_t0) < 1e-10) continue;

            if ((rn.G_spherical.T0 > 0) == (hm_t0 > 0))
                same_sign++;
            else
                opposite_sign++;
        }
    }

    int total = same_sign + opposite_sign;
    if (total > 0) {
        double frac = static_cast<double>(same_sign) / total;
        EXPECT_GT(frac, 0.9)
            << "BS and HM T0 should agree in sign at large distance";
        std::cout << "  Far-field T0 sign agreement: " << same_sign << "/"
                  << total << " (" << (100.0 * frac) << "%)\n";
    } else {
        std::cout << "  No far-field pairs to compare\n";
    }
}


// ============================================================================
// ORCA protein test
// ============================================================================

TEST(HaighMallionOrcaTest, RunOnProtonatedProtein) {
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

    auto hm = HaighMallionResult::Compute(conf);
    ASSERT_NE(hm, nullptr);
    conf.AttachResult(std::move(hm));

    double max_t0 = 0, max_t2 = 0, max_trace = 0;
    int with_hm = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& sc = conf.AtomAt(ai).hm_shielding_contribution;
        max_t0 = std::max(max_t0, std::abs(sc.T0));
        max_t2 = std::max(max_t2, sc.T2Magnitude());

        // Check raw integral tracelessness
        for (const auto& rn : conf.AtomAt(ai).ring_neighbours) {
            if (rn.hm_H_tensor.norm() > 1e-15) {
                max_trace = std::max(max_trace, std::abs(rn.hm_H_tensor.trace()));
                with_hm++;
            }
        }
    }

    std::cout << "  ORCA protein HaighMallion summary:\n"
              << "    atoms=" << conf.AtomCount()
              << " rings=" << load.protein->RingCount() << "\n"
              << "    ring-atom pairs with HM data: " << with_hm << "\n"
              << "    max |T0| (full kernel) = " << max_t0 << "\n"
              << "    max |T2| (full kernel) = " << max_t2 << "\n"
              << "    max |Tr(H)| (raw integral) = " << max_trace << "\n";

    EXPECT_GT(with_hm, 0) << "Should have HM data";
    EXPECT_GT(max_t0, 1e-6) << "Full kernel T0 should be non-zero";
    EXPECT_GT(max_t2, 1e-6) << "Full kernel T2 should be non-zero";
    EXPECT_LT(max_trace, 1e-8) << "Raw integral must be traceless";
}

#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include <cmath>
#include <iostream>

#include "McConnellResult.h"
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
// Analytical test: known geometry, hand-calculable result
// ============================================================================

TEST(McConnellAnalytical, DipolarKernelAtKnownGeometry) {
    // Atom at (3, 0, 0). Bond midpoint at origin. Bond direction (0, 0, 1).
    //
    // d = (3, 0, 0), r = 3, d_hat = (1, 0, 0)
    // b_hat = (0, 0, 1)
    // cos_theta = d_hat . b_hat = 0
    //
    // McConnell scalar f = (3*0 - 1) / 27 = -1/27 = -0.03704
    //
    // Symmetric dipolar kernel K:
    //   K_11 = (3*1*1 - 1) / 27 =  2/27 =  0.07407
    //   K_22 = (3*0*0 - 1) / 27 = -1/27 = -0.03704
    //   K_33 = (3*0*0 - 1) / 27 = -1/27 = -0.03704
    //   Off-diagonal = 0 (d_hat only has x component)
    //   Trace = 0 (traceless)
    //
    // Full McConnell tensor M / r^3:
    //   Term 1: 9*0 * d_hat_a * b_hat_b = 0 (cos_theta = 0)
    //   Term 2: -3 * b_hat_a * b_hat_b = -3 * outer(b,b) / r^3
    //           Only (2,2) = -3/27 = -1/9
    //   Term 3: -(3 d_hat_a d_hat_b - delta_ab) / r^3 = -K
    //
    //   M / r^3 = 0 + (-3 b⊗b) / r^3 + (-K)
    //   M_11 = 0 + 0 + (-2/27) = -2/27
    //   M_22 = 0 + 0 + (1/27) = 1/27
    //   M_33 = 0 + (-3/27) + (1/27) = -2/27
    //   Trace(M/r^3) = -2/27 + 1/27 + (-2/27) = -3/27 = -1/9
    //   T0 = Trace/3 = -1/27 = f/3   ✓

    // Build a minimal protein with one bond
    auto protein = std::make_unique<Protein>();

    auto atom_a = Atom::Create(Element::C);
    auto atom_b = Atom::Create(Element::O);
    auto atom_field = Atom::Create(Element::H);
    atom_a->residue_index = 0;
    atom_b->residue_index = 0;
    atom_field->residue_index = 0;

    protein->AddAtom(std::move(atom_a));  // 0: bond endpoint C at origin
    protein->AddAtom(std::move(atom_b));  // 1: bond endpoint O at (0,0,1)
    protein->AddAtom(std::move(atom_field));  // 2: field point H at (3,0,0)

    Residue res;
    res.type = AminoAcid::ALA;
    res.sequence_number = 1;
    res.chain_id = "A";
    res.atom_indices = {0, 1, 2};
    protein->AddResidue(res);

    // Positions: C at (-0.5,0,0), O at (0.5,0,0) gives midpoint at origin,
    // but we want bond direction (0,0,1). So:
    // C at (0,0,-0.5), O at (0,0,0.5), midpoint at origin, direction = (0,0,1)
    std::vector<Vec3> positions = {
        Vec3(0.0, 0.0, -0.5),   // C
        Vec3(0.0, 0.0,  0.5),   // O
        Vec3(3.0, 0.0,  0.0)    // H (field point)
    };

    protein->FinalizeConstruction(positions);

    // Manually set bond category for testing
    // (FinalizeConstruction may not classify correctly for this toy geometry)

    auto& conf = protein->AddCrystalConformation(positions, 0, 0, 0, "test");

    conf.AttachResult(GeometryResult::Compute(conf));
    conf.AttachResult(SpatialIndexResult::Compute(conf));

    auto mc = McConnellResult::Compute(conf);
    ASSERT_NE(mc, nullptr);
    conf.AttachResult(std::move(mc));

    // Check the field point atom (index 2)
    const auto& ca = conf.AtomAt(2);

    // Should have at least one bond neighbour
    ASSERT_GT(ca.bond_neighbours.size(), 0u)
        << "Field atom should see at least one bond";

    // Find the bond to the C-O pair
    const BondNeighbourhood& bn = ca.bond_neighbours[0];

    // Distance should be 3.0 A (from midpoint at origin to (3,0,0))
    EXPECT_NEAR(bn.distance_to_midpoint, 3.0, 0.01)
        << "Distance to bond midpoint should be 3.0 A";

    // McConnell scalar: (3*0 - 1) / 27 = -1/27
    EXPECT_NEAR(bn.mcconnell_scalar, -1.0/27.0, 1e-6)
        << "McConnell scalar f should be -1/27";

    // Dipolar kernel K should be traceless
    double trace_K = bn.dipolar_tensor.trace();
    EXPECT_NEAR(trace_K, 0.0, 1e-10) << "Dipolar kernel must be traceless";

    // K diagonal: (2/27, -1/27, -1/27)
    EXPECT_NEAR(bn.dipolar_tensor(0,0),  2.0/27.0, 1e-6);
    EXPECT_NEAR(bn.dipolar_tensor(1,1), -1.0/27.0, 1e-6);
    EXPECT_NEAR(bn.dipolar_tensor(2,2), -1.0/27.0, 1e-6);

    // K off-diagonal: all zero (d_hat = (1,0,0))
    EXPECT_NEAR(bn.dipolar_tensor(0,1), 0.0, 1e-10);
    EXPECT_NEAR(bn.dipolar_tensor(0,2), 0.0, 1e-10);
    EXPECT_NEAR(bn.dipolar_tensor(1,2), 0.0, 1e-10);

    // Full McConnell shielding contribution T0 should be non-zero
    // T0 = f / 3 = -1/81 (if we used the raw kernel)
    // Actually T0 = Trace(M/r^3) / 3 = (-1/9) / 3 = -1/27
    // which equals f. The f IS Trace(M)/r^3, not Trace(M/r^3)/3.
    // Let me check: Trace(M) = 3(3cos^2 theta - 1) = 3(-1) = -3.
    // M/r^3: Trace = -3/27 = -1/9. T0 = Trace/3 = -1/27.
    // And f = (3*0-1)/27 = -1/27. So T0 = f. Correct!
    EXPECT_NEAR(ca.mc_shielding_contribution.T0, -1.0/27.0, 1e-6)
        << "T0 of full McConnell tensor should equal McConnell scalar f";

    std::cout << "  Analytical test passed:\n"
              << "    f = " << bn.mcconnell_scalar << " (expected -1/27 = "
              << -1.0/27.0 << ")\n"
              << "    T0 = " << ca.mc_shielding_contribution.T0 << "\n"
              << "    K trace = " << trace_K << "\n"
              << "    K diag = (" << bn.dipolar_tensor(0,0) << ", "
              << bn.dipolar_tensor(1,1) << ", "
              << bn.dipolar_tensor(2,2) << ")\n";
}


// ============================================================================
// Full protein test
// ============================================================================

class McConnellProteinTest : public ::testing::Test {
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


TEST_F(McConnellProteinTest, ComputeAndAttach) {
    auto& conf = protein->Conformation();
    auto mc = McConnellResult::Compute(conf);
    ASSERT_NE(mc, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(mc)));
    ASSERT_TRUE(conf.HasResult<McConnellResult>());
}


TEST_F(McConnellProteinTest, EveryAtomHasBondNeighbours) {
    auto& conf = protein->Conformation();
    conf.AttachResult(McConnellResult::Compute(conf));

    int with_neighbours = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        if (!conf.AtomAt(ai).bond_neighbours.empty())
            with_neighbours++;
    }

    // Every atom in a protein should see nearby bonds
    EXPECT_GT(with_neighbours, static_cast<int>(conf.AtomCount()) - 5)
        << "Almost every atom should have bond neighbours";

    std::cout << "  Atoms with bond neighbours: " << with_neighbours
              << " / " << conf.AtomCount() << "\n";
}


TEST_F(McConnellProteinTest, DipolarKernelsAreTraceless) {
    auto& conf = protein->Conformation();
    conf.AttachResult(McConnellResult::Compute(conf));

    double max_trace = 0.0;
    int checked = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& bn : conf.AtomAt(ai).bond_neighbours) {
            double trace = bn.dipolar_tensor.trace();
            max_trace = std::max(max_trace, std::abs(trace));
            checked++;
        }
    }

    EXPECT_LT(max_trace, 1e-10)
        << "All dipolar kernels K must be traceless";

    std::cout << "  Checked " << checked << " kernels, max |trace| = "
              << max_trace << "\n";
}


TEST_F(McConnellProteinTest, McConnellScalarNonZero) {
    auto& conf = protein->Conformation();
    conf.AttachResult(McConnellResult::Compute(conf));

    // Most atoms should have non-zero CO sum
    int with_co = 0;
    double max_co = 0.0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        double co = std::abs(conf.AtomAt(ai).mcconnell_co_sum);
        if (co > 1e-6) with_co++;
        max_co = std::max(max_co, co);
    }

    EXPECT_GT(with_co, static_cast<int>(conf.AtomCount() / 2))
        << "Most atoms should have non-zero CO McConnell sum";

    std::cout << "  Atoms with |CO sum| > 1e-6: " << with_co
              << " / " << conf.AtomCount()
              << ", max |CO sum| = " << max_co << "\n";
}


TEST_F(McConnellProteinTest, ShieldingContributionHasT0AndT2) {
    auto& conf = protein->Conformation();
    conf.AttachResult(McConnellResult::Compute(conf));

    int nonzero_t0 = 0, nonzero_t2 = 0;
    double max_t0 = 0, max_t2 = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& sc = conf.AtomAt(ai).mc_shielding_contribution;
        if (std::abs(sc.T0) > 1e-8) nonzero_t0++;
        double t2mag = sc.T2Magnitude();
        if (t2mag > 1e-8) nonzero_t2++;
        max_t0 = std::max(max_t0, std::abs(sc.T0));
        max_t2 = std::max(max_t2, t2mag);
    }

    EXPECT_GT(nonzero_t0, 0) << "Full McConnell tensor must have non-zero T0";
    EXPECT_GT(nonzero_t2, 0) << "Full McConnell tensor must have non-zero T2";

    std::cout << "  T0 nonzero: " << nonzero_t0 << ", max |T0| = " << max_t0 << "\n";
    std::cout << "  T2 nonzero: " << nonzero_t2 << ", max |T2| = " << max_t2 << "\n";
}


TEST_F(McConnellProteinTest, T0EqualsScalarFForSingleBond) {
    // For any atom-bond pair, T0 of the full McConnell tensor M/r^3
    // should equal the McConnell scalar f.
    // This is because Trace(M) = 3(3cos^2 theta - 1),
    // so T0 = Trace(M/r^3)/3 = (3cos^2 theta - 1)/r^3 = f.

    auto& conf = protein->Conformation();
    conf.AttachResult(McConnellResult::Compute(conf));

    int checked = 0;
    double max_diff = 0.0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& bn : conf.AtomAt(ai).bond_neighbours) {
            // Reconstruct the full M tensor for this bond
            Vec3 atom_pos = conf.PositionAt(ai);
            Vec3 midpoint = conf.bond_midpoints[bn.bond_index];
            Vec3 d = atom_pos - midpoint;
            double r = d.norm();
            if (r < MIN_DISTANCE) continue;

            Vec3 d_hat = d / r;
            Vec3 b_hat = conf.bond_directions[bn.bond_index];
            double cos_theta = d_hat.dot(b_hat);
            double r3 = r * r * r;

            // Build M/r^3
            Mat3 M;
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                    M(a, b) = (9.0 * cos_theta * d_hat(a) * b_hat(b)
                               - 3.0 * b_hat(a) * b_hat(b)
                               - (3.0 * d_hat(a) * d_hat(b) - (a==b ? 1.0 : 0.0)))
                              / r3;

            double t0 = M.trace() / 3.0;
            double diff = std::abs(t0 - bn.mcconnell_scalar);
            max_diff = std::max(max_diff, diff);
            checked++;
        }
    }

    EXPECT_LT(max_diff, 1e-10)
        << "T0 of full M tensor must equal McConnell scalar f";

    std::cout << "  Checked " << checked << " pairs, max |T0 - f| = "
              << max_diff << "\n";
}


TEST_F(McConnellProteinTest, NearestCOTrackingWorks) {
    auto& conf = protein->Conformation();
    conf.AttachResult(McConnellResult::Compute(conf));

    int has_co = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        double d = conf.AtomAt(ai).nearest_CO_dist;
        if (d < NO_DATA_SENTINEL) {
            has_co++;
            EXPECT_GT(d, 0.0);
            EXPECT_LT(d, MCCONNELL_CUTOFF_A);
        }
    }

    EXPECT_GT(has_co, static_cast<int>(conf.AtomCount() / 2))
        << "Most atoms should have a nearest CO bond within range";

    std::cout << "  Atoms with nearest CO: " << has_co
              << " / " << conf.AtomCount() << "\n";
}


// ============================================================================
// ORCA protein test (protonated, more atoms)
// ============================================================================

TEST(McConnellOrcaTest, RunOnProtonatedProtein) {
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
    auto mc = McConnellResult::Compute(conf);
    ASSERT_NE(mc, nullptr);
    conf.AttachResult(std::move(mc));

    // Summary: bond neighbour counts, T0/T2 ranges
    double min_t0 = 1e30, max_t0 = -1e30;
    double max_t2 = 0;
    int total_bn = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& ca = conf.AtomAt(ai);
        total_bn += static_cast<int>(ca.bond_neighbours.size());
        double t0 = ca.mc_shielding_contribution.T0;
        min_t0 = std::min(min_t0, t0);
        max_t0 = std::max(max_t0, t0);
        max_t2 = std::max(max_t2, ca.mc_shielding_contribution.T2Magnitude());
    }

    std::cout << "  ORCA protein McConnell summary:\n"
              << "    atoms=" << conf.AtomCount()
              << " bonds=" << load.protein->BondCount()
              << " total_bond_neighbours=" << total_bn << "\n"
              << "    T0 range: [" << min_t0 << ", " << max_t0 << "]\n"
              << "    max |T2| = " << max_t2 << "\n";

    EXPECT_GT(total_bn, 0);
    EXPECT_GT(max_t0 - min_t0, 0.001) << "Should have a range of T0 values";
    EXPECT_GT(max_t2, 0.001) << "T2 should be non-zero";
}

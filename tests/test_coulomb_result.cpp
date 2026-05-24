#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include <cmath>
#include <iostream>

#include "CoulombResult.h"
#include "ChargeAssignmentResult.h"
#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "ApbsFieldResult.h"
#include "Protein.h"
#include "PdbFileReader.h"
#include "OrcaRunLoader.h"
#include "ChargeSource.h"
#include "PhysicalConstants.h"

#include <filesystem>
namespace fs = std::filesystem;
using namespace nmr;



// ============================================================================
// Analytical test: two point charges, hand-calculable E-field and EFG
// ============================================================================

TEST(CoulombAnalytical, TwoChargesKnownGeometry) {
    // Atom 0: charge +0.5 e at (0, 0, 0)
    // Atom 1: charge -0.3 e at (0, 0, 0)  -- dummy, both at same position
    // Actually: set up atoms with known positions and charges.
    //
    // Field point (atom 2) at (3, 0, 0), charge 0 (observer).
    // Source atom 0 at origin, charge +0.5 e.
    // Source atom 1 at (0, 3, 0), charge -0.3 e.
    //
    // E at atom 2 from atom 0:
    //   r = (3,0,0) - (0,0,0) = (3,0,0), |r| = 3
    //   E_0 = ke * 0.5 * (3,0,0) / 27 = ke * (0.5/27) * (3,0,0)
    //       = ke * (1/54, 0, 0) * 3 = ke * (1/18, 0, 0)
    //   Wait: E = ke * q * r / r^3 = ke * 0.5 * (3,0,0) / 27
    //       = ke * (1.5/27, 0, 0) = ke * (1/18, 0, 0)
    //
    // E at atom 2 from atom 1:
    //   r = (3,0,0) - (0,3,0) = (3,-3,0), |r| = 3*sqrt(2)
    //   r^3 = 27 * 2*sqrt(2) = 54*sqrt(2)
    //   E_1 = ke * (-0.3) * (3,-3,0) / (54*sqrt(2))
    //       = ke * (-0.3/54/sqrt(2)) * (3, -3, 0)
    //
    // Let's just verify the properties rather than exact values:
    // 1. EFG is traceless
    // 2. E is non-zero
    // 3. The total has the right order of magnitude

    auto protein = std::make_unique<Protein>();

    auto a0 = Atom::Create(Element::C);
    auto a1 = Atom::Create(Element::N);
    auto a2 = Atom::Create(Element::H);
    a0->residue_index = 0;
    a1->residue_index = 0;
    a2->residue_index = 0;

    protein->AddAtom(std::move(a0));
    protein->AddAtom(std::move(a1));
    protein->AddAtom(std::move(a2));

    Residue res;
    res.type = AminoAcid::ALA;
    res.sequence_number = 1;
    res.chain_id = "A";
    res.atom_indices = {0, 1, 2};
    protein->AddResidue(res);

    std::vector<Vec3> positions = {
        Vec3(0.0, 0.0, 0.0),   // atom 0: source, q = +0.5
        Vec3(0.0, 3.0, 0.0),   // atom 1: source, q = -0.3
        Vec3(3.0, 0.0, 0.0)    // atom 2: observer, q = 0
    };

    protein->FinalizeConstruction(positions);
    auto& conf = protein->AddCrystalConformation(positions, 0, 0, 0, "test");

    // Attach dependencies
    conf.AttachResult(GeometryResult::Compute(conf));

    // ChargeAssignmentResult: stub first, then overwrite with known values.
    // Stub assigns uniform 0.1; we overwrite AFTER attachment.
    conf.AttachResult(ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams().c_str()));

    conf.MutableAtomAt(0).partial_charge = 0.5;
    conf.MutableAtomAt(1).partial_charge = -0.3;
    conf.MutableAtomAt(2).partial_charge = 0.0;

    conf.AttachResult(SpatialIndexResult::Compute(conf));

    auto coulomb = CoulombResult::Compute(conf);
    ASSERT_NE(coulomb, nullptr);
    conf.AttachResult(std::move(coulomb));

    // Check observer atom (index 2)
    const auto& ca = conf.AtomAt(2);

    // E-field should be non-zero
    EXPECT_GT(ca.coulomb_E_magnitude, 0.0)
        << "Observer should see non-zero E-field";

    // EFG should be traceless
    double trace = ca.coulomb_EFG_total.trace();
    EXPECT_NEAR(trace, 0.0, 1e-10)
        << "EFG must be traceless (Gauss's law)";

    // Verify analytically: E at atom 2 from atom 0 only
    // r = (3,0,0), |r| = 3, r^3 = 27
    // E = ke * 0.5 * (3,0,0) / 27 = ke * (1/18, 0, 0)
    double E_x_from_0 = COULOMB_KE * 0.5 * 3.0 / 27.0;  // ke/18
    double E_y_from_0 = 0.0;

    // E from atom 1: r = (3,-3,0), |r| = 3*sqrt(2), r^3 = 54*sqrt(2)
    double r1 = 3.0 * std::sqrt(2.0);
    double r1_3 = r1 * r1 * r1;
    double E_x_from_1 = COULOMB_KE * (-0.3) * 3.0 / r1_3;
    double E_y_from_1 = COULOMB_KE * (-0.3) * (-3.0) / r1_3;

    double E_x_expected = E_x_from_0 + E_x_from_1;
    double E_y_expected = E_y_from_0 + E_y_from_1;

    EXPECT_NEAR(ca.coulomb_E_total(0), E_x_expected, 1e-8)
        << "E_x should match analytical value";
    EXPECT_NEAR(ca.coulomb_E_total(1), E_y_expected, 1e-8)
        << "E_y should match analytical value";
    EXPECT_NEAR(ca.coulomb_E_total(2), 0.0, 1e-10)
        << "E_z should be zero (all atoms in xy plane)";

    // Verify EFG diagonal from atom 0 alone (d along x):
    // V_xx = ke * q * (3*1 - 1)/r^3 = ke * 0.5 * 2/27
    // V_yy = ke * q * (3*0 - 1)/r^3 = ke * 0.5 * (-1)/27
    // V_zz = ke * q * (3*0 - 1)/r^3 = ke * 0.5 * (-1)/27
    // Sum(V_diag) = ke * 0.5 * (2 - 1 - 1)/27 = 0  ✓

    std::cout << "  Analytical Coulomb test:\n"
              << "    E = (" << ca.coulomb_E_total(0) << ", "
              << ca.coulomb_E_total(1) << ", "
              << ca.coulomb_E_total(2) << ") V/A\n"
              << "    |E| = " << ca.coulomb_E_magnitude << " V/A\n"
              << "    EFG trace = " << trace << "\n"
              << "    EFG T2 magnitude = "
              << ca.coulomb_EFG_total_spherical.T2Magnitude() << "\n";
}


// ============================================================================
// Full protein test
// ============================================================================

class CoulombProteinTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!fs::exists(nmr::test::TestEnvironment::UbqProtonated()) || !fs::exists(nmr::test::TestEnvironment::Ff14sbParams())) {
            GTEST_SKIP() << "1UBQ.pdb or ff14sb_params.dat not found";
        }
        auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
        if (!r.Ok()) GTEST_SKIP() << "Failed to load 1UBQ";
        protein = std::move(r.protein);

        auto& conf = protein->Conformation();
        conf.AttachResult(GeometryResult::Compute(conf));
        conf.AttachResult(ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams()));
        conf.AttachResult(SpatialIndexResult::Compute(conf));
    }

    std::unique_ptr<Protein> protein;
};


TEST_F(CoulombProteinTest, ComputeAndAttach) {
    auto& conf = protein->Conformation();
    auto coulomb = CoulombResult::Compute(conf);
    ASSERT_NE(coulomb, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(coulomb)));
    ASSERT_TRUE(conf.HasResult<CoulombResult>());
}


TEST_F(CoulombProteinTest, EFieldNonZeroForAllAtoms) {
    auto& conf = protein->Conformation();
    conf.AttachResult(CoulombResult::Compute(conf));

    int nonzero = 0;
    double sum_mag = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        double mag = conf.AtomAt(ai).coulomb_E_magnitude;
        if (mag > 1e-10) nonzero++;
        sum_mag += mag;
    }

    double mean = sum_mag / static_cast<double>(conf.AtomCount());
    EXPECT_EQ(nonzero, static_cast<int>(conf.AtomCount()))
        << "Every atom should have non-zero E-field";

    // Mean E should be on the order of 0.1-10 V/A for a protein
    EXPECT_GT(mean, 0.01) << "Mean |E| too small";
    EXPECT_LT(mean, 100.0) << "Mean |E| too large";

    std::cout << "  Non-zero E: " << nonzero << " / " << conf.AtomCount()
              << ", mean |E| = " << mean << " V/A\n";
}


TEST_F(CoulombProteinTest, EFGIsTraceless) {
    auto& conf = protein->Conformation();
    conf.AttachResult(CoulombResult::Compute(conf));

    double max_trace = 0.0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        double trace = std::abs(conf.AtomAt(ai).coulomb_EFG_total.trace());
        max_trace = std::max(max_trace, trace);
    }

    EXPECT_LT(max_trace, 1e-8)
        << "All EFG tensors must be traceless (Gauss's law)";

    std::cout << "  Max |EFG trace| = " << max_trace << "\n";
}


TEST_F(CoulombProteinTest, DecompositionSumsToTotal) {
    auto& conf = protein->Conformation();
    conf.AttachResult(CoulombResult::Compute(conf));

    double max_E_diff = 0.0;
    double max_EFG_diff = 0.0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& ca = conf.AtomAt(ai);

        // E_total = E_backbone + E_sidechain + E_aromatic
        Vec3 E_sum = ca.coulomb_E_backbone + ca.coulomb_E_sidechain
                   + ca.coulomb_E_aromatic;
        double E_diff = (ca.coulomb_E_total - E_sum).norm();
        max_E_diff = std::max(max_E_diff, E_diff);

        // EFG_total = EFG_backbone + EFG_sidechain + EFG_aromatic
        // (sidechain EFG not stored separately, but total - backbone - aromatic
        //  should give it)
        Mat3 EFG_sum = ca.coulomb_EFG_backbone + ca.coulomb_EFG_aromatic;
        // We don't have EFG_sidechain stored as Mat3, but the total should
        // equal backbone + sidechain + aromatic. Check E decomposition only.
    }

    EXPECT_LT(max_E_diff, 1e-8)
        << "E decomposition must sum to total";

    std::cout << "  Max |E_total - (bb + sc + arom)| = " << max_E_diff << "\n";
}


TEST_F(CoulombProteinTest, BackboneFractionIsReasonable) {
    auto& conf = protein->Conformation();
    conf.AttachResult(CoulombResult::Compute(conf));

    double sum_proj = 0;
    int count = 0;
    int positive = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        double proj = conf.AtomAt(ai).coulomb_E_backbone_frac;
        // Projection of E_backbone along E_total direction.
        // Positive = backbone aligned with total. Negative = opposed.
        // Bounded by |E_backbone|.
        sum_proj += proj;
        if (proj > 0) positive++;
        count++;
    }

    double mean_proj = sum_proj / count;
    // Most atoms should have positive projection (backbone dominates)
    EXPECT_GT(positive, count / 2)
        << "Majority of atoms should have backbone aligned with total E";
    EXPECT_GT(mean_proj, 0.0)
        << "Mean backbone projection should be positive";

    std::cout << "  Mean backbone E projection = " << mean_proj << " V/A\n"
              << "  Positive projection: " << positive << " / " << count << "\n";
}


TEST_F(CoulombProteinTest, T2IsNonZero) {
    auto& conf = protein->Conformation();
    conf.AttachResult(CoulombResult::Compute(conf));

    int nonzero_t2 = 0;
    double max_t2 = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        double t2 = conf.AtomAt(ai).coulomb_EFG_total_spherical.T2Magnitude();
        if (t2 > 1e-8) nonzero_t2++;
        max_t2 = std::max(max_t2, t2);
    }

    EXPECT_EQ(nonzero_t2, static_cast<int>(conf.AtomCount()))
        << "Every atom should have non-zero EFG T2";
    EXPECT_GT(max_t2, 0.01) << "Max T2 should be appreciable";

    std::cout << "  T2 nonzero: " << nonzero_t2
              << ", max |T2| = " << max_t2 << "\n";
}


// ============================================================================
// ORCA protein test (protonated, with prmtop charges)
// ============================================================================

TEST(CoulombOrcaTest, RunOnProtonatedProtein) {
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

    PrmtopChargeSource charge_source(files.prmtop_path);
    conf.AttachResult(ChargeAssignmentResult::Compute(conf, charge_source));
    conf.AttachResult(SpatialIndexResult::Compute(conf));

    auto coulomb = CoulombResult::Compute(conf);
    ASSERT_NE(coulomb, nullptr);
    conf.AttachResult(std::move(coulomb));

    // Summary
    double min_E = 1e30, max_E = 0;
    double max_t2 = 0;
    double max_trace = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& ca = conf.AtomAt(ai);
        double E_mag = ca.coulomb_E_magnitude;
        min_E = std::min(min_E, E_mag);
        max_E = std::max(max_E, E_mag);
        max_t2 = std::max(max_t2, ca.coulomb_EFG_total_spherical.T2Magnitude());
        max_trace = std::max(max_trace, std::abs(ca.coulomb_EFG_total.trace()));
    }

    std::cout << "  ORCA protein Coulomb summary:\n"
              << "    atoms=" << conf.AtomCount() << "\n"
              << "    |E| range: [" << min_E << ", " << max_E << "] V/A\n"
              << "    max |T2| = " << max_t2 << " V/A^2\n"
              << "    max |EFG trace| = " << max_trace << "\n";

    // Sanity checks
    EXPECT_GT(max_E, 0.1) << "Should have appreciable E-fields";
    EXPECT_LT(max_trace, 1e-8) << "EFG must be traceless";
    EXPECT_GT(max_t2, 0.01) << "T2 should be non-zero";
}


// ============================================================================
// APBS comparison: if both Coulomb and APBS are present, verify that the
// solvent contribution (APBS - vacuum) is physically reasonable.
// ============================================================================

TEST(CoulombApbsComparison, SolventContributionIsReasonable) {
    if (!fs::exists(nmr::test::TestEnvironment::UbqProtonated()) || !fs::exists(nmr::test::TestEnvironment::Ff14sbParams()))
        GTEST_SKIP() << "Test data not found";

    auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
    if (!r.Ok()) GTEST_SKIP() << r.error;

    auto& conf = r.protein->Conformation();
    conf.AttachResult(GeometryResult::Compute(conf));
    conf.AttachResult(ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams()));
    conf.AttachResult(SpatialIndexResult::Compute(conf));

    // APBS first (needs charges)
    auto apbs = ApbsFieldResult::Compute(conf);
    if (!apbs) GTEST_SKIP() << "APBS failed";
    conf.AttachResult(std::move(apbs));

    // Then Coulomb (will compute solvent = APBS - vacuum)
    conf.AttachResult(CoulombResult::Compute(conf));

    int has_solvent = 0;
    double mean_ratio = 0;
    int count = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& ca = conf.AtomAt(ai);
        double solv_mag = ca.coulomb_E_solvent.norm();
        double vac_mag = ca.coulomb_E_magnitude;
        if (solv_mag > 1e-10) has_solvent++;
        if (vac_mag > 1e-6) {
            mean_ratio += solv_mag / vac_mag;
            count++;
        }
    }
    if (count > 0) mean_ratio /= count;

    std::cout << "  APBS vs vacuum Coulomb:\n"
              << "    Atoms with non-zero solvent E: " << has_solvent
              << " / " << conf.AtomCount() << "\n"
              << "    Mean |E_solvent|/|E_vacuum| = " << mean_ratio << "\n";

    // Solvation should modify the E-field (screening), so solvent
    // contribution should be non-zero for most atoms
    EXPECT_GT(has_solvent, static_cast<int>(conf.AtomCount() / 2))
        << "Most atoms should have non-zero solvent contribution";
}

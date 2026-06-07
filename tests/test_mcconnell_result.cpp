#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iterator>
#include <cstdio>

#include "McConnellResult.h"
#include "MopacResult.h"
#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "Protein.h"
#include "PdbFileReader.h"
#include "OrcaRunLoader.h"
#include "PhysicalConstants.h"

#include <filesystem>
#include <Eigen/Geometry>
#include <unistd.h>
namespace fs = std::filesystem;
using namespace nmr;

namespace {

double MaxAbs(const Mat3& m) {
    double v = 0.0;
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c)
            v = std::max(v, std::abs(m(r, c)));
    return v;
}

bool SphericalAllZero(const SphericalTensor& st) {
    if (st.T0 != 0.0) return false;
    for (double v : st.T1) if (v != 0.0) return false;
    for (double v : st.T2) if (v != 0.0) return false;
    return true;
}

std::string ReadText(const fs::path& path) {
    std::ifstream in(path);
    return std::string(std::istreambuf_iterator<char>(in),
                       std::istreambuf_iterator<char>());
}

std::unique_ptr<Protein> BuildSyntheticPeptideCOProtein(
        bool degenerate = false,
        const Mat3& rotation = Mat3::Identity()) {
    auto protein = std::make_unique<Protein>();

    auto atom_n = Atom::Create(Element::N);
    auto atom_c = Atom::Create(Element::C);
    auto atom_o = Atom::Create(Element::O);
    auto atom_h = Atom::Create(Element::H);
    atom_n->pdb_atom_name = "N";
    atom_c->pdb_atom_name = "C";
    atom_o->pdb_atom_name = "O";
    atom_h->pdb_atom_name = "H";
    atom_n->residue_index = 0;
    atom_c->residue_index = 0;
    atom_o->residue_index = 0;
    atom_h->residue_index = 0;

    protein->AddAtom(std::move(atom_n));  // 0: amide N
    protein->AddAtom(std::move(atom_c));  // 1: carbonyl C
    protein->AddAtom(std::move(atom_o));  // 2: carbonyl O
    protein->AddAtom(std::move(atom_h));  // 3: target probe

    Residue res;
    res.type = AminoAcid::ALA;
    res.sequence_number = 1;
    res.chain_id = "A";
    res.atom_indices = {0, 1, 2, 3};
    res.N = 0;
    res.C = 1;
    res.O = 2;
    res.H = 3;
    protein->AddResidue(res);

    std::vector<Vec3> positions = {
        degenerate ? Vec3(-1.45, 0.0, 0.0) : Vec3(0.0, 1.45, 0.0),
        Vec3(0.0, 0.0, 0.0),
        Vec3(1.2, 0.0, 0.0),
        Vec3(0.6, 3.0, 0.0)
    };
    for (Vec3& pos : positions) pos = rotation * pos;

    protein->FinalizeConstruction(positions);
    protein->AddCrystalConformation(positions, 0, 0, 0, "synthetic_peptide_co");
    return protein;
}

size_t FirstBondWithCategory(const Protein& protein, BondCategory category) {
    for (size_t bi = 0; bi < protein.BondCount(); ++bi) {
        if (protein.BondAt(bi).category == category)
            return bi;
    }
    return SIZE_MAX;
}

size_t FirstBondBetween(const Protein& protein, size_t a, size_t b) {
    for (size_t bi = 0; bi < protein.BondCount(); ++bi) {
        const Bond& bond = protein.BondAt(bi);
        if ((bond.atom_index_a == a && bond.atom_index_b == b) ||
            (bond.atom_index_a == b && bond.atom_index_b == a))
            return bi;
    }
    return SIZE_MAX;
}

}  // namespace



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
    // Unit axial source shape:
    //   Qhat = u u^T - I/3 = diag(-1/3, -1/3, 2/3)
    // Clean McConnell response:
    //   A = D(r) Qhat
    // PCS scalar branch:
    //   T0(A) = trace(A)/3 = n^T Qhat n / r^3 = (-1/3)/27 = -1/81
    //
    // Symmetric dipolar kernel K:
    //   K_11 = (3*1*1 - 1) / 27 =  2/27 =  0.07407
    //   K_22 = (3*0*0 - 1) / 27 = -1/27 = -0.03704
    //   K_33 = (3*0*0 - 1) / 27 = -1/27 = -0.03704
    //   Off-diagonal = 0 (d_hat only has x component)
    //   Trace = 0 (traceless)
    //
    // Response A = D Qhat:
    //   diag(-2/81, 1/81, -2/81), so trace(A)/3 = -1/81.

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

    // PCS scalar: n^T Qhat n / r^3 = -1/81.
    EXPECT_NEAR(bn.mcconnell_scalar, -1.0/81.0, 1e-6)
        << "McConnell PCS scalar should be -1/81";

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

    EXPECT_NEAR(ca.mc_shielding_contribution.T0, -1.0/81.0, 1e-6)
        << "T0 of D(r)Qhat should equal the PCS scalar";

    std::cout << "  Analytical test passed:\n"
              << "    PCS scalar = " << bn.mcconnell_scalar
              << " (expected -1/81 = " << -1.0/81.0 << ")\n"
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


TEST_F(McConnellProteinTest, PCSScalarMatchesTraceBranchForSingleBond) {
    auto& conf = protein->Conformation();
    conf.AttachResult(McConnellResult::Compute(conf));

    int checked = 0;
    double max_diff = 0.0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& bn : conf.AtomAt(ai).bond_neighbours) {
            Vec3 atom_pos = conf.PositionAt(ai);
            Vec3 midpoint = conf.bond_midpoints[bn.bond_index];
            Vec3 axis = conf.bond_directions[bn.bond_index];
            auto kernel = McConnellResult::ComputePairKernel(
                atom_pos, midpoint, axis);
            double r = kernel.distance;
            if (r < MIN_DISTANCE) continue;

            double r3 = r * r * r;
            double pcs = kernel.direction.dot(
                kernel.source_shape * kernel.direction) / r3;
            double t0 = SphericalTensor::Decompose(kernel.response).T0;
            max_diff = std::max(max_diff, std::abs(t0 - bn.mcconnell_scalar));
            max_diff = std::max(max_diff, std::abs(t0 - pcs));
            checked++;
        }
    }

    EXPECT_LT(max_diff, 1e-10)
        << "T0(DQhat) must equal the PCS scalar n^T Qhat n/r^3";

    std::cout << "  Checked " << checked
              << " pairs, max PCS scalar branch diff = "
              << max_diff << "\n";
}


TEST(McConnellImplementationChecks, RotationEquivariance) {
    Vec3 source_center(0.3, -0.7, 1.1);
    Vec3 source_axis(0.4, 0.8, -0.2);
    Vec3 target(3.2, 1.4, -0.6);

    const auto k = McConnellResult::ComputePairKernel(
        target, source_center, source_axis);
    const SphericalTensor st = SphericalTensor::Decompose(k.response);

    const Mat3 R = Eigen::AngleAxisd(
        0.73, Vec3(0.2, -0.4, 0.9).normalized()).toRotationMatrix();
    const auto kr = McConnellResult::ComputePairKernel(
        R * target, R * source_center, R * source_axis);
    const SphericalTensor sr = SphericalTensor::Decompose(kr.response);

    EXPECT_NEAR(sr.T0, st.T0, 1e-12);

    Vec3 t1(st.T1[0], st.T1[1], st.T1[2]);
    Vec3 t1_rot(sr.T1[0], sr.T1[1], sr.T1[2]);
    EXPECT_LT((t1_rot - R * t1).norm(), 1e-12);

    SphericalTensor st_t2 = st;
    st_t2.T0 = 0.0;
    st_t2.T1 = {0.0, 0.0, 0.0};
    SphericalTensor sr_t2 = sr;
    sr_t2.T0 = 0.0;
    sr_t2.T1 = {0.0, 0.0, 0.0};
    const Mat3 expected_t2 = R * st_t2.Reconstruct() * R.transpose();
    EXPECT_LT(MaxAbs(sr_t2.Reconstruct() - expected_t2), 1e-12);

    EXPECT_LT(MaxAbs(kr.response - R * k.response * R.transpose()), 1e-12);
}


TEST(McConnellImplementationChecks, PairPCSScalarIdentity) {
    const auto k = McConnellResult::ComputePairKernel(
        Vec3(2.0, -1.0, 3.0),
        Vec3(-0.5, 0.25, 0.75),
        Vec3(0.3, 0.7, 0.2));
    ASSERT_GT(k.distance, 0.0);
    const double r3 = k.distance * k.distance * k.distance;
    const double pcs = k.direction.dot(k.source_shape * k.direction) / r3;
    const double t0 = SphericalTensor::Decompose(k.response).T0;
    EXPECT_NEAR(t0, k.response.trace() / 3.0, 1e-15);
    EXPECT_NEAR(t0, pcs, 1e-15);
}


TEST(McConnellImplementationChecks, PeptideCORhombicSourceShapePinnedAnalytic) {
    auto protein = BuildSyntheticPeptideCOProtein();
    auto& conf = protein->Conformation();
    conf.AttachResult(GeometryResult::Compute(conf));
    const size_t co_bond = FirstBondWithCategory(*protein, BondCategory::PeptideCO);
    ASSERT_NE(co_bond, SIZE_MAX);

    const auto k = McConnellResult::ComputePeptideCORhombicPairKernel(
        conf, co_bond, conf.PositionAt(3));

    // C->O = +x, plane normal = +z, e_in = +y.  From
    // (-5.4, +4.0, -14) the mean-subtracted principal values are
    // (-4/15, 137/15, -133/15); dividing by axial scale 137/10 gives
    // (out, para, in) = (-8/411, 2/3, -266/411).
    Mat3 expected = Mat3::Zero();
    expected(0, 0) = 2.0 / 3.0;
    expected(1, 1) = -266.0 / 411.0;
    expected(2, 2) = -8.0 / 411.0;

    EXPECT_LT(MaxAbs(k.source_shape - expected), 1e-12);
    EXPECT_NEAR(k.source_shape.trace(), 0.0, 1e-15);
}


TEST(McConnellImplementationChecks, PeptideCORhombicSourceShapeRotationEquivariance) {
    const Mat3 R = Eigen::AngleAxisd(
        0.91, Vec3(0.3, -0.7, 0.2).normalized()).toRotationMatrix();

    auto protein = BuildSyntheticPeptideCOProtein();
    auto rotated_protein = BuildSyntheticPeptideCOProtein(false, R);
    auto& conf = protein->Conformation();
    auto& rotated_conf = rotated_protein->Conformation();
    conf.AttachResult(GeometryResult::Compute(conf));
    rotated_conf.AttachResult(GeometryResult::Compute(rotated_conf));
    const size_t co_bond = FirstBondWithCategory(*protein, BondCategory::PeptideCO);
    const size_t rotated_co_bond =
        FirstBondWithCategory(*rotated_protein, BondCategory::PeptideCO);
    ASSERT_NE(co_bond, SIZE_MAX);
    ASSERT_NE(rotated_co_bond, SIZE_MAX);

    const auto k = McConnellResult::ComputePeptideCORhombicPairKernel(
        conf, co_bond, conf.PositionAt(3));
    const auto kr = McConnellResult::ComputePeptideCORhombicPairKernel(
        rotated_conf, rotated_co_bond, rotated_conf.PositionAt(3));

    EXPECT_LT(MaxAbs(kr.source_shape -
                     R * k.source_shape * R.transpose()), 1e-12);
    EXPECT_LT(MaxAbs(kr.response -
                     R * k.response * R.transpose()), 1e-12);
}


TEST(McConnellImplementationChecks, PeptideCORhombicContributionNonZeroAndTraceless) {
    auto protein = BuildSyntheticPeptideCOProtein();
    auto& conf = protein->Conformation();
    conf.AttachResult(GeometryResult::Compute(conf));
    conf.AttachResult(SpatialIndexResult::Compute(conf));
    conf.AttachResult(McConnellResult::Compute(conf));

    const size_t co_bond = FirstBondWithCategory(*protein, BondCategory::PeptideCO);
    ASSERT_NE(co_bond, SIZE_MAX);
    const auto k = McConnellResult::ComputePeptideCORhombicPairKernel(
        conf, co_bond, conf.PositionAt(3));
    EXPECT_NEAR(k.source_shape.trace(), 0.0, 1e-15);

    // The emitted rhombic array is the additive delta relative to the
    // legacy axial PeptideCO shape.  Here n = e_in = +y, r = 3 A, so
    // D = diag(-1/27, 2/27, -1/27) and
    // Q_delta = diag(0, -43/137, +43/137).
    Mat3 expected_delta = Mat3::Zero();
    expected_delta(1, 1) = -86.0 / 3699.0;
    expected_delta(2, 2) = -43.0 / 3699.0;

    const Mat3 actual =
        conf.AtomAt(3).mcconnell_peptide_co_rhombic.Reconstruct();
    EXPECT_GT(MaxAbs(actual), 1e-6);
    EXPECT_LT(MaxAbs(actual - expected_delta), 1e-12);
}


TEST(McConnellImplementationChecks, PeptideCORhombicFallsBackToAxial) {
    {
        auto protein = BuildSyntheticPeptideCOProtein();
        auto& conf = protein->Conformation();
        conf.AttachResult(GeometryResult::Compute(conf));
        const size_t nc_bond = FirstBondBetween(*protein, 0, 1);
        ASSERT_NE(nc_bond, SIZE_MAX);
        ASSERT_NE(protein->BondAt(nc_bond).category, BondCategory::PeptideCO);

        const auto axial = McConnellResult::ComputePairKernel(
            conf.PositionAt(3),
            conf.bond_midpoints[nc_bond],
            conf.bond_directions[nc_bond]);
        const auto rhombic = McConnellResult::ComputePeptideCORhombicPairKernel(
            conf, nc_bond, conf.PositionAt(3));

        EXPECT_LT(MaxAbs(rhombic.source_shape - axial.source_shape), 1e-12);
        EXPECT_LT(MaxAbs(rhombic.response - axial.response), 1e-12);
    }

    {
        auto protein = BuildSyntheticPeptideCOProtein(true);
        auto& conf = protein->Conformation();
        conf.AttachResult(GeometryResult::Compute(conf));
        const size_t co_bond = FirstBondWithCategory(*protein, BondCategory::PeptideCO);
        ASSERT_NE(co_bond, SIZE_MAX);

        const auto axial = McConnellResult::ComputePairKernel(
            conf.PositionAt(3),
            conf.bond_midpoints[co_bond],
            conf.bond_directions[co_bond]);
        const auto rhombic = McConnellResult::ComputePeptideCORhombicPairKernel(
            conf, co_bond, conf.PositionAt(3));

        EXPECT_LT(MaxAbs(rhombic.source_shape - axial.source_shape), 1e-12);
        EXPECT_LT(MaxAbs(rhombic.response - axial.response), 1e-12);

        conf.AttachResult(SpatialIndexResult::Compute(conf));
        conf.AttachResult(McConnellResult::Compute(conf));
        EXPECT_LT(MaxAbs(conf.AtomAt(3)
                             .mcconnell_peptide_co_rhombic.Reconstruct()),
                  1e-12);
    }
}


TEST_F(McConnellProteinTest, AromaticCategoryIsExactlyZero) {
    const Protein& p = protein->Conformation().ProteinRef();
    int aromatic_bonds = 0;
    for (size_t bi = 0; bi < p.BondCount(); ++bi)
        if (p.BondAt(bi).category == BondCategory::Aromatic)
            ++aromatic_bonds;
    ASSERT_GT(aromatic_bonds, 0);

    auto& conf = protein->Conformation();
    conf.AttachResult(McConnellResult::Compute(conf));
    const size_t cat = static_cast<size_t>(
        McConnellSourceCategory::AromaticZeroed);
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        EXPECT_TRUE(SphericalAllZero(
            conf.AtomAt(ai).mcconnell_source_tensors[cat][0]));
        EXPECT_TRUE(SphericalAllZero(
            conf.AtomAt(ai).mcconnell_source_tensors[cat][1]));
    }
}


TEST_F(McConnellProteinTest, MissingMopacZerosOnlyBoChannel) {
    auto& conf = protein->Conformation();
    ASSERT_FALSE(conf.HasResult<MopacResult>());
    conf.AttachResult(McConnellResult::Compute(conf));

    double max_fixed_t2 = 0.0;
    double max_bo_abs = 0.0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& ca = conf.AtomAt(ai);
        for (size_t c = 0; c < kMcConnellSourceCategoryCount; ++c) {
            max_fixed_t2 = std::max(
                max_fixed_t2,
                ca.mcconnell_source_tensors[c][0].T2Magnitude());
            const auto& bo = ca.mcconnell_source_tensors[c][1];
            max_bo_abs = std::max(max_bo_abs, std::abs(bo.T0));
            for (double v : bo.T1) max_bo_abs = std::max(max_bo_abs, std::abs(v));
            for (double v : bo.T2) max_bo_abs = std::max(max_bo_abs, std::abs(v));
        }
    }
    EXPECT_GT(max_fixed_t2, 0.0);
    EXPECT_EQ(max_bo_abs, 0.0);
}


TEST_F(McConnellProteinTest, NearFieldAuditCountsAreReported) {
    auto& conf = protein->Conformation();
    conf.AttachResult(McConnellResult::Compute(conf));

    int accepted = 0;
    int rejected = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        accepted += conf.AtomAt(ai).mcconnell_near_field_accepted_lt3A;
        rejected += conf.AtomAt(ai).mcconnell_near_field_rejected_lt3A;
    }
    EXPECT_GE(accepted, 0);
    EXPECT_GE(rejected, 0);
    std::cout << "  McConnell near-field audit (<3 A): accepted="
              << accepted << " rejected=" << rejected << "\n";
}


TEST_F(McConnellProteinTest, WriteFeaturesEmitsTwentyOneArraysAndManifest) {
    auto& conf = protein->Conformation();
    conf.AttachResult(McConnellResult::Compute(conf));

    const fs::path out_dir = fs::temp_directory_path() /
        ("mcconnell_features_" + std::to_string(::getpid()));
    fs::create_directories(out_dir);
    const auto& mc = conf.Result<McConnellResult>();
    EXPECT_EQ(mc.WriteFeatures(conf, out_dir.string()), 21);

    for (size_t c = 0; c < kMcConnellSourceCategoryCount; ++c) {
        const auto cat = static_cast<McConnellSourceCategory>(c);
        for (size_t h = 0; h < kMcConnellChannelCount; ++h) {
            const auto channel = static_cast<McConnellChannel>(h);
            const fs::path p = out_dir / (
                std::string("mc_") + McConnellSourceCategoryStem(cat) +
                "_" + McConnellChannelStem(channel) + ".npy");
            EXPECT_TRUE(fs::exists(p)) << p;
        }
    }
    EXPECT_TRUE(fs::exists(out_dir / "mc_peptide_co_rhombic.npy"));
    EXPECT_TRUE(fs::exists(out_dir / "mc_nearest_co_dir.npy"));
    EXPECT_TRUE(fs::exists(out_dir / "mc_nearest_co_midpoint.npy"));
    EXPECT_TRUE(fs::exists(out_dir / "mc_nearest_co_T2.npy"));
    EXPECT_TRUE(fs::exists(out_dir / "mc_nearest_cn_T2.npy"));
    EXPECT_TRUE(fs::exists(out_dir / "mc_bond_neighbors.npy"));

    const fs::path manifest = out_dir / "extraction_manifest.json";
    ASSERT_TRUE(fs::exists(manifest));
    const std::string text = ReadText(manifest);
    EXPECT_NE(text.find("\"source_model\": \"unit susceptibility shape; axial scale learned; peptide C=O rhombic scale pinned\""),
              std::string::npos);
    EXPECT_NE(text.find("\"bo_source\": \"MOPAC Wiberg bond order\""),
              std::string::npos);
    EXPECT_NE(text.find("\"aromatic_zeroed_when_ring_active\": true"),
              std::string::npos);
    EXPECT_NE(text.find("\"irrep_layout\": \"0e,1e_x,1e_y,1e_z,2e_m-2..+2\""),
              std::string::npos);
    EXPECT_NE(text.find("\"units\": \"Angstrom^-3\""),
              std::string::npos);
    EXPECT_NE(text.find("\"rhombic_status\": \"peptide_co_pinned_present\""),
              std::string::npos);
    EXPECT_NE(text.find("\"rhombic_array\": \"mc_peptide_co_rhombic.npy\""),
              std::string::npos);
    EXPECT_NE(text.find("\"source\": \"Hooper & Kaiser 1965 Table III, EF-corrected acetamide A, Abraham-anchored sign\""),
              std::string::npos);
    EXPECT_TRUE(fs::exists(out_dir / "mc_nearfield_counts.npy"));

    for (size_t c = 0; c < kMcConnellSourceCategoryCount; ++c) {
        const auto cat = static_cast<McConnellSourceCategory>(c);
        for (size_t h = 0; h < kMcConnellChannelCount; ++h) {
            const auto channel = static_cast<McConnellChannel>(h);
            const fs::path p = out_dir / (
                std::string("mc_") + McConnellSourceCategoryStem(cat) +
                "_" + McConnellChannelStem(channel) + ".npy");
            std::remove(p.string().c_str());
        }
    }
    std::remove((out_dir / "mc_peptide_co_rhombic.npy").string().c_str());
    std::remove((out_dir / "mc_nearfield_counts.npy").string().c_str());
    std::remove(manifest.string().c_str());
    ::rmdir(out_dir.string().c_str());
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

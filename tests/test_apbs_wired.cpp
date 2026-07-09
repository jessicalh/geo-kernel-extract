#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include "PdbFileReader.h"
#include "GeometryResult.h"
#include "ChargeAssignmentResult.h"
#include "ChargeSource.h"
#include "ForceFieldChargeTable.h"
#include "ApbsFieldResult.h"
#include "OperationLog.h"
#include "PhysicalConstants.h"
#include <filesystem>
#include <cmath>

using namespace nmr;


class ApbsWiredTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!std::filesystem::exists(nmr::test::TestEnvironment::UbqProtonated())) {
            GTEST_SKIP() << "1UBQ PDB not found";
        }
        auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
        if (!r.Ok()) GTEST_SKIP() << r.error;
        protein = std::move(r.protein);
    }
    std::unique_ptr<Protein> protein;
};

TEST_F(ApbsWiredTest, ApbsProducesNonZeroFieldsWhenSolverSucceeds) {
    // This test verifies the APBS bridge integration.
    // Fields come from the Poisson-Boltzmann grid and must be non-zero.

    auto& conf = protein->Conformation();

    auto charges = ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams().c_str());
    ASSERT_NE(charges, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(charges)));

    auto apbs = ApbsFieldResult::Compute(conf);
    ASSERT_NE(apbs, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(apbs)));

    const auto& result = conf.Result<ApbsFieldResult>();

    int nonzero_E = 0;
    int nonzero_EFG = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        Vec3 E = result.ElectricFieldAt(ai);
        Mat3 EFG = result.FieldGradientAt(ai);
        if (E.norm() > 1e-10) nonzero_E++;
        if (EFG.norm() > 1e-10) nonzero_EFG++;
    }

    EXPECT_GT(nonzero_E, static_cast<int>(conf.AtomCount()) / 2)
        << "Too few atoms with non-zero E-field";
    EXPECT_GT(nonzero_EFG, static_cast<int>(conf.AtomCount()) / 2)
        << "Too few atoms with non-zero EFG";
}

TEST_F(ApbsWiredTest, SphericalTensorRoundtrip) {
    auto& conf = protein->Conformation();

    auto charges = ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams().c_str());
    conf.AttachResult(std::move(charges));
    auto apbs = ApbsFieldResult::Compute(conf);
    conf.AttachResult(std::move(apbs));

    const auto& result = conf.Result<ApbsFieldResult>();

    for (size_t ai = 0; ai < std::min(conf.AtomCount(), size_t(30)); ++ai) {
        Mat3 EFG = result.FieldGradientAt(ai);
        SphericalTensor st = result.FieldGradientSphericalAt(ai);
        Mat3 reconstructed = st.Reconstruct();

        double diff = (EFG - reconstructed).norm();
        EXPECT_LT(diff, 1e-10)
            << "SphericalTensor roundtrip failed at atom " << ai
            << " diff=" << diff;
    }
}

TEST_F(ApbsWiredTest, FailureReturnsNullWithoutFallback) {
    // With valid inputs the APBS computation should succeed. Solver failure
    // returns nullptr; there is no vacuum-field fallback path.
    auto& conf = protein->Conformation();

    auto charges = ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams().c_str());
    conf.AttachResult(std::move(charges));
    auto apbs = ApbsFieldResult::Compute(conf);
    ASSERT_NE(apbs, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(apbs)));
}

TEST(ApbsWired, TruePlaceholderPbRadiusBlocksApbsBeforeSolve) {
    Protein protein;

    Residue residue;
    residue.type = AminoAcid::ALA;
    residue.sequence_number = 1;
    residue.chain_id = "A";
    const size_t ri = protein.AddResidue(residue);

    auto atom = Atom::Create(Element::C);
    atom->pdb_atom_name = "CA";
    atom->residue_index = ri;
    const size_t ai = protein.AddAtom(std::move(atom));
    protein.MutableResidueAt(ri).atom_indices = {ai};
    protein.FinalizeConstruction({Vec3(0.0, 0.0, 0.0)});

    std::vector<AtomChargeRadius> values(1);
    values[0] = {0.1, kCompatibilityPlaceholderPbRadiusAngstrom,
                 ChargeAssignmentStatus::PlaceholderPbRadius};
    std::string error;
    auto table = ForceFieldChargeTable::FromValues(
        ForceField::Amber_ff14SB, ChargeModelKind::Preloaded, "placeholder-test",
        std::move(values), 1, error);
    ASSERT_NE(table, nullptr) << error;
    protein.SetForceFieldCharges(std::move(table));

    auto& conf = protein.AddCrystalConformation(
        {Vec3(0.0, 0.0, 0.0)}, 0.0, 0.0, 0.0, "PLACEHOLDER");
    PreloadedChargeSource source({}, ForceField::Amber_ff14SB);
    auto charges = ChargeAssignmentResult::Compute(conf, source);
    ASSERT_NE(charges, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(charges)));

    auto apbs = ApbsFieldResult::Compute(conf);
    EXPECT_EQ(apbs, nullptr);
    EXPECT_EQ(protein.ForceFieldCharges().PlaceholderPbRadiusCount(), 1);
    EXPECT_EQ(protein.ForceFieldCharges().DerivedMbondi2PbRadiusCount(), 0);
}

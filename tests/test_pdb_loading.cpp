#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include "PdbFileReader.h"
#include <filesystem>

using namespace nmr;

// Path to 1UBQ test PDB

class PdbLoadingTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!std::filesystem::exists(nmr::test::TestEnvironment::UbqProtonated())) {
            GTEST_SKIP() << "1UBQ PDB not found at " << nmr::test::TestEnvironment::UbqProtonated();
        }
    }
};

TEST_F(PdbLoadingTest, LoadsSuccessfully) {
    auto result = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(result.Ok()) << result.error;
}

TEST_F(PdbLoadingTest, HasResidues) {
    auto result = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(result.Ok());
    // 1UBQ has 76 residues
    EXPECT_GE(result.protein->ResidueCount(), 70);
    EXPECT_LE(result.protein->ResidueCount(), 80);
}

TEST_F(PdbLoadingTest, HasAtoms) {
    auto result = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(result.Ok());
    // 1UBQ has ~600 heavy atoms, more with H
    EXPECT_GT(result.protein->AtomCount(), 500);
}

TEST_F(PdbLoadingTest, HasAromaticRings) {
    auto result = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(result.Ok());
    // 1UBQ: F4(PHE), F45(PHE), Y59(TYR), H68(HIS) = at least 4 residues
    // PHE has 1 ring, TYR has 1, HIS has 1
    EXPECT_GE(result.protein->RingCount(), 3);
}

TEST_F(PdbLoadingTest, HasCovalentBonds) {
    auto result = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(result.Ok());
    EXPECT_GT(result.protein->BondCount(), 100);
}

TEST_F(PdbLoadingTest, HasConformation) {
    auto result = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(result.Ok());
    EXPECT_GE(result.protein->ConformationCount(), 1u);

    auto& conf = result.protein->Conformation();
    EXPECT_EQ(conf.AtomCount(), result.protein->AtomCount());
}

TEST_F(PdbLoadingTest, BackboneIndicesCached) {
    auto result = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(result.Ok());

    const Residue& first = result.protein->ResidueAt(0);
    EXPECT_NE(first.N, Residue::NONE);
    EXPECT_NE(first.CA, Residue::NONE);
    EXPECT_NE(first.C, Residue::NONE);
}

TEST_F(PdbLoadingTest, BuildContextSet) {
    auto result = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(result.Ok());
    EXPECT_FALSE(result.protein->BuildContext().pdb_source.empty());
}

TEST_F(PdbLoadingTest, AtomElementsCorrect) {
    auto result = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(result.Ok());

    // Check that backbone N of first residue is nitrogen
    const Residue& first = result.protein->ResidueAt(0);
    if (first.N != Residue::NONE) {
        EXPECT_EQ(result.protein->AtomAt(first.N).element, Element::N);
    }
}

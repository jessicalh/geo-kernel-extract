#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include "Atom.h"
#include "PdbFileReader.h"
#include "GeometryResult.h"
#include <filesystem>

using namespace nmr;

// ============================================================================
// Flat Atom hierarchy tests
// ============================================================================

TEST(AtomFlatTest, HydrogenHasSettableParentIndex) {
    auto h = Atom::Create(Element::H);
    // Default parent_atom_index is SIZE_MAX (no parent assigned yet)
    EXPECT_EQ(h->parent_atom_index, SIZE_MAX);

    // Set the parent
    h->parent_atom_index = 42;
    EXPECT_EQ(h->parent_atom_index, 42u);
}

TEST(AtomFlatTest, CarbonHasDefaultParentIndex) {
    auto c = Atom::Create(Element::C);
    // Non-hydrogen atoms should have parent_atom_index == SIZE_MAX
    EXPECT_EQ(c->parent_atom_index, SIZE_MAX);
}

TEST(AtomFlatTest, NitrogenHasDefaultParentIndex) {
    auto n = Atom::Create(Element::N);
    EXPECT_EQ(n->parent_atom_index, SIZE_MAX);
}

TEST(AtomFlatTest, CovalentRadiusCorrectForAllElements) {
    EXPECT_DOUBLE_EQ(Atom::Create(Element::H)->CovalentRadius(), 0.31);
    EXPECT_DOUBLE_EQ(Atom::Create(Element::C)->CovalentRadius(), 0.76);
    EXPECT_DOUBLE_EQ(Atom::Create(Element::N)->CovalentRadius(), 0.71);
    EXPECT_DOUBLE_EQ(Atom::Create(Element::O)->CovalentRadius(), 0.66);
    EXPECT_DOUBLE_EQ(Atom::Create(Element::S)->CovalentRadius(), 1.05);
}

TEST(AtomFlatTest, ElectronegativityCorrectForAllElements) {
    EXPECT_DOUBLE_EQ(Atom::Create(Element::H)->Electronegativity(), 2.20);
    EXPECT_DOUBLE_EQ(Atom::Create(Element::C)->Electronegativity(), 2.55);
    EXPECT_DOUBLE_EQ(Atom::Create(Element::N)->Electronegativity(), 3.04);
    EXPECT_DOUBLE_EQ(Atom::Create(Element::O)->Electronegativity(), 3.44);
    EXPECT_DOUBLE_EQ(Atom::Create(Element::S)->Electronegativity(), 2.58);
}

TEST(AtomFlatTest, HBondDonorForNAndO) {
    EXPECT_TRUE(Atom::Create(Element::N)->IsHBondDonorElement());
    EXPECT_TRUE(Atom::Create(Element::O)->IsHBondDonorElement());
    EXPECT_FALSE(Atom::Create(Element::C)->IsHBondDonorElement());
    EXPECT_FALSE(Atom::Create(Element::H)->IsHBondDonorElement());
    EXPECT_FALSE(Atom::Create(Element::S)->IsHBondDonorElement());
}

TEST(AtomFlatTest, HBondAcceptorForNAndO) {
    EXPECT_TRUE(Atom::Create(Element::N)->IsHBondAcceptorElement());
    EXPECT_TRUE(Atom::Create(Element::O)->IsHBondAcceptorElement());
    EXPECT_FALSE(Atom::Create(Element::C)->IsHBondAcceptorElement());
    EXPECT_FALSE(Atom::Create(Element::H)->IsHBondAcceptorElement());
    EXPECT_FALSE(Atom::Create(Element::S)->IsHBondAcceptorElement());
}

TEST(AtomFlatTest, CreateFromString) {
    auto c = Atom::Create("C");
    EXPECT_EQ(c->element, Element::C);

    auto h = Atom::Create("H");
    EXPECT_EQ(h->element, Element::H);

    auto unk = Atom::Create("X");
    EXPECT_EQ(unk->element, Element::Unknown);
}

TEST(AtomFlatTest, NoDynamicCastNeeded) {
    // All atoms are the same concrete Atom class.
    // parent_atom_index is directly accessible without downcasting.
    auto h = Atom::Create(Element::H);
    h->parent_atom_index = 99;
    EXPECT_EQ(h->parent_atom_index, 99u);

    // For non-H atoms, parent_atom_index is SIZE_MAX (unused)
    auto c = Atom::Create(Element::C);
    EXPECT_EQ(c->parent_atom_index, SIZE_MAX);

    // Both are the same type (no virtual dispatch needed)
    EXPECT_EQ(h->CovalentRadius(), CovalentRadiusForElement(Element::H));
    EXPECT_EQ(c->CovalentRadius(), CovalentRadiusForElement(Element::C));
}

// ============================================================================
// Integration: PDB loading still works with flat atoms
// ============================================================================


TEST(AtomFlatTest, PdbLoadingStillWorks) {
    if (!std::filesystem::exists(nmr::test::TestEnvironment::UbqProtonated())) {
        GTEST_SKIP() << "1UBQ PDB not found";
    }

    auto result = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(result.Ok()) << result.error;

    // Basic sanity checks
    EXPECT_GT(result.protein->AtomCount(), 500u);
    EXPECT_GT(result.protein->BondCount(), 100u);
    EXPECT_GT(result.protein->RingCount(), 0u);

    // Hydrogen parent assignment: if the PDB has H atoms, at least some
    // should have parent_atom_index set. Heavy-atom-only crystal structures
    // (like 1UBQ) have no hydrogens, so the test passes trivially.
    int h_count = 0, h_with_parent = 0;
    for (size_t ai = 0; ai < result.protein->AtomCount(); ++ai) {
        const auto& atom = result.protein->AtomAt(ai);
        if (atom.element == Element::H) {
            h_count++;
            if (atom.parent_atom_index != SIZE_MAX)
                h_with_parent++;
        }
    }
    if (h_count > 0) {
        EXPECT_GT(h_with_parent, 0) << "PDB has " << h_count
            << " H atoms but none have parent_atom_index set";
    }
}

TEST(AtomFlatTest, RingDetectionStillWorks) {
    if (!std::filesystem::exists(nmr::test::TestEnvironment::UbqProtonated())) {
        GTEST_SKIP() << "1UBQ PDB not found";
    }

    auto result = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(result.Ok()) << result.error;

    // Ring detection should still produce typed ring objects
    EXPECT_GE(result.protein->RingCount(), 3u);

    // Check that ring atoms are valid indices
    for (size_t ri = 0; ri < result.protein->RingCount(); ++ri) {
        const auto& ring = result.protein->RingAt(ri);
        for (size_t ai : ring.atom_indices) {
            EXPECT_LT(ai, result.protein->AtomCount())
                << "Ring " << ri << " has invalid atom index " << ai;
        }
    }
}

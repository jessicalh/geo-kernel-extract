// tests/test_atom_semantic_predicates.cpp
//
// Unit tests for the compositional predicates on AtomSemanticTable.
// Specifically covers the three predicates added to support the
// LarsenHBondGrid HBondHα design (see
// spec/plan/larsen-hbond-shielding-design-2026-05-11.md §"Substrate
// predicate additions"):
//
//   - IsAnyAlphaHydrogen        — non-GLY HA + GLY HA2/HA3
//   - IsSidechainCarboxylateOxygen — ASP OD1/OD2 + GLU OE1/OE2 + C-term
//     carboxylate O's (anything with PlanarGroupKind::Carboxylate and
//     Element::O, including C-term deprotonated/protonated overrides).
//   - IsSidechainAmideOxygen     — ASN OD1 + GLN OE1
//
// Strategy: load 1ubq_protonated.pdb (the canonical AMBER-prepared
// fixture); for each predicate scan every atom and count matches,
// then compare against counts derived from the 1UBQ residue table.

#include <gtest/gtest.h>

#include <filesystem>

#include "TestEnvironment.h"
#include "Protein.h"
#include "PdbFileReader.h"
#include "AminoAcidType.h"
#include "Atom.h"
#include "Residue.h"
#include "LegacyAmberTopology.h"
#include "SemanticEnums.h"

using nmr::AminoAcid;
using nmr::Atom;
using nmr::AtomSemanticTable;
using nmr::Protein;
using nmr::Residue;
using nmr::ResidueTerminalState;

namespace {

struct PredicateCounts {
    int any_alpha_h = 0;
    int sidechain_carboxylate_o = 0;
    int sidechain_amide_o = 0;

    // Residue-type tallies — used to derive expected counts directly
    // from the 1UBQ residue composition rather than hard-coding numbers
    // a maintainer would have to remember to update if the fixture
    // changed.
    int gly_count = 0;
    int non_gly_count = 0;
    int asp_count = 0;
    int glu_count = 0;
    int asn_count = 0;
    int gln_count = 0;
    int c_term_count = 0;
};

PredicateCounts ScanProtein(const Protein& protein) {
    PredicateCounts c;
    const auto& topology = protein.LegacyAmber();
    EXPECT_TRUE(topology.HasAtomSemantic())
        << "LegacyAmberTopology must populate semantic store on "
           "FinalizeConstruction";

    for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        if (res.type == AminoAcid::Unknown) continue;
        if (res.type == AminoAcid::GLY) ++c.gly_count;
        else                            ++c.non_gly_count;
        if (res.type == AminoAcid::ASP) ++c.asp_count;
        if (res.type == AminoAcid::GLU) ++c.glu_count;
        if (res.type == AminoAcid::ASN) ++c.asn_count;
        if (res.type == AminoAcid::GLN) ++c.gln_count;
        if (res.terminal_state == ResidueTerminalState::CTerminus ||
            res.terminal_state == ResidueTerminalState::NAndCTerminus) {
            ++c.c_term_count;
        }

        for (size_t ai : res.atom_indices) {
            const AtomSemanticTable& row = topology.SemanticAt(ai);
            if (row.IsAnyAlphaHydrogen())            ++c.any_alpha_h;
            if (row.IsSidechainCarboxylateOxygen())  ++c.sidechain_carboxylate_o;
            if (row.IsSidechainAmideOxygen())        ++c.sidechain_amide_o;
        }
    }
    return c;
}

class AtomSemanticPredicatesTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!std::filesystem::exists(nmr::test::TestEnvironment::UbqProtonated())) {
            GTEST_SKIP() << "1UBQ PDB not found at "
                         << nmr::test::TestEnvironment::UbqProtonated();
        }
        auto r = nmr::BuildFromProtonatedPdb(
            nmr::test::TestEnvironment::UbqProtonated());
        ASSERT_TRUE(r.Ok()) << r.error;
        protein = std::move(r.protein);
    }

    std::unique_ptr<Protein> protein;
};

}  // namespace


// ============================================================================
// IsAnyAlphaHydrogen — every residue contributes 1 HA, GLY contributes 2.
// ============================================================================

TEST_F(AtomSemanticPredicatesTest, IsAnyAlphaHydrogenCountsHaPerResidue) {
    auto c = ScanProtein(*protein);
    const int expected = c.non_gly_count + 2 * c.gly_count;
    EXPECT_EQ(c.any_alpha_h, expected)
        << "Each non-GLY residue should contribute 1 alpha H (HA); "
           "each GLY should contribute 2 (HA2 + HA3). "
           "non-GLY=" << c.non_gly_count
        << " GLY=" << c.gly_count;
    // Sanity: 1UBQ has 76 residues w/ 4 or 6 GLY depending on labeling.
    EXPECT_GE(c.any_alpha_h, 76)
        << "At least one alpha H per residue expected";
}


// ============================================================================
// IsSidechainCarboxylateOxygen — ASP OD1/OD2, GLU OE1/OE2, plus C-term
// carboxylate O atoms (which the generated cap tables stamp with
// PlanarGroupKind::Carboxylate on BOTH O and OXT).
// ============================================================================

TEST_F(AtomSemanticPredicatesTest, IsSidechainCarboxylateOxygenCountsAspGluAndCterm) {
    auto c = ScanProtein(*protein);
    const int expected =
        2 * c.asp_count  // OD1 + OD2
      + 2 * c.glu_count  // OE1 + OE2
      + 2 * c.c_term_count;  // C-term carboxylate O + OXT
    EXPECT_EQ(c.sidechain_carboxylate_o, expected)
        << "ASP " << c.asp_count << " * 2 + GLU " << c.glu_count
        << " * 2 + C-term " << c.c_term_count << " * 2 = " << expected;
    // 1UBQ has 5 ASP + 6 GLU + 1 C-term → 10+12+2 = 24 expected.
    EXPECT_GE(c.sidechain_carboxylate_o, 10)
        << "1UBQ must show at least ASP+GLU contribution";
}


// ============================================================================
// IsSidechainAmideOxygen — ASN OD1 + GLN OE1 only (the corresponding
// ND2 / NE2 are nitrogens, filtered out by the element check).
// ============================================================================

TEST_F(AtomSemanticPredicatesTest, IsSidechainAmideOxygenCountsAsnGln) {
    auto c = ScanProtein(*protein);
    const int expected = c.asn_count + c.gln_count;
    EXPECT_EQ(c.sidechain_amide_o, expected)
        << "ASN " << c.asn_count << " + GLN " << c.gln_count
        << " sidechain amide O's expected";
}

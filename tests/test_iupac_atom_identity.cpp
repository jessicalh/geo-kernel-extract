#include <gtest/gtest.h>
#include "IupacAtomIdentity.h"

using namespace nmr;

TEST(IupacAtomIdentityTest, LookupAlaCA) {
    const auto* id = IupacAtomIdentity::Lookup("ALA", "CA");
    ASSERT_NE(id, nullptr);
    EXPECT_EQ(id->element, Element::C);
    EXPECT_TRUE(id->is_backbone);
    EXPECT_FALSE(id->is_ring_member);
}

TEST(IupacAtomIdentityTest, LookupPheCG) {
    const auto* id = IupacAtomIdentity::Lookup("PHE", "CG");
    ASSERT_NE(id, nullptr);
    EXPECT_EQ(id->element, Element::C);
    EXPECT_TRUE(id->is_ring_member);
}

TEST(IupacAtomIdentityTest, LookupUnknownAtomReturnsNull) {
    const auto* id = IupacAtomIdentity::Lookup("ALA", "ZZZ");
    EXPECT_EQ(id, nullptr);
}

TEST(IupacAtomIdentityTest, LookupUnknownResidueReturnsNull) {
    const auto* id = IupacAtomIdentity::Lookup("ZZZ", "CA");
    EXPECT_EQ(id, nullptr);
}

TEST(IupacAtomIdentityTest, AtomsForResidueGlyHasCorrectCount) {
    // GLY: N, CA, C, O, H, HA2, HA3 = 7 atoms
    const auto& atoms = IupacAtomIdentity::AtomsForResidue("GLY");
    EXPECT_EQ(atoms.size(), 7u);
}

TEST(IupacAtomIdentityTest, AtomsForResidueAla) {
    // ALA: N, CA, C, O, H, HA, CB, HB1, HB2, HB3 = 10 atoms
    const auto& atoms = IupacAtomIdentity::AtomsForResidue("ALA");
    EXPECT_EQ(atoms.size(), 10u);
}

TEST(IupacAtomIdentityTest, AtomsForUnknownResidueReturnsEmpty) {
    const auto& atoms = IupacAtomIdentity::AtomsForResidue("ZZZ");
    EXPECT_TRUE(atoms.empty());
}

TEST(IupacAtomIdentityTest, ValidateCorrect) {
    std::string err = IupacAtomIdentity::Validate("ALA", "CA", Element::C, true);
    EXPECT_TRUE(err.empty()) << err;
}

TEST(IupacAtomIdentityTest, ValidateWrongElement) {
    std::string err = IupacAtomIdentity::Validate("ALA", "CA", Element::N, true);
    EXPECT_FALSE(err.empty());
    EXPECT_NE(err.find("element mismatch"), std::string::npos);
}

TEST(IupacAtomIdentityTest, ValidateWrongBackbone) {
    // CB in ALA is sidechain (is_backbone = false)
    std::string err = IupacAtomIdentity::Validate("ALA", "CB", Element::C, true);
    EXPECT_FALSE(err.empty());
    EXPECT_NE(err.find("backbone flag mismatch"), std::string::npos);
}

TEST(IupacAtomIdentityTest, ValidateNonexistentAtom) {
    std::string err = IupacAtomIdentity::Validate("ALA", "ZZZ", Element::C, false);
    EXPECT_FALSE(err.empty());
    EXPECT_NE(err.find("not part of"), std::string::npos);
}

TEST(IupacAtomIdentityTest, PheRingMembersIdentified) {
    // PHE ring members: CG, CD1, CD2, CE1, CE2, CZ (6 carbons)
    const auto& atoms = IupacAtomIdentity::AtomsForResidue("PHE");
    int ring_count = 0;
    for (const auto& a : atoms) {
        if (a.is_ring_member) ring_count++;
    }
    // PHE has 6 ring carbons
    EXPECT_EQ(ring_count, 6);
}

TEST(IupacAtomIdentityTest, HisRingMembersIdentified) {
    // HIS ring members: CG, ND1, CE1, NE2, CD2 (5 atoms)
    const auto& atoms = IupacAtomIdentity::AtomsForResidue("HIS");
    int ring_count = 0;
    for (const auto& a : atoms) {
        if (a.is_ring_member) ring_count++;
    }
    EXPECT_EQ(ring_count, 5);
}

TEST(IupacAtomIdentityTest, BackboneAtomIdentity) {
    // Test backbone atoms across several amino acids
    for (const char* res : {"ALA", "PHE", "GLY", "TRP", "HIS"}) {
        const auto* n = IupacAtomIdentity::Lookup(res, "N");
        ASSERT_NE(n, nullptr) << res << " N";
        EXPECT_EQ(n->element, Element::N) << res << " N element";
        EXPECT_TRUE(n->is_backbone) << res << " N backbone";

        const auto* ca = IupacAtomIdentity::Lookup(res, "CA");
        ASSERT_NE(ca, nullptr) << res << " CA";
        EXPECT_EQ(ca->element, Element::C) << res << " CA element";
        EXPECT_TRUE(ca->is_backbone) << res << " CA backbone";

        const auto* c = IupacAtomIdentity::Lookup(res, "C");
        ASSERT_NE(c, nullptr) << res << " C";
        EXPECT_EQ(c->element, Element::C) << res << " C element";
        EXPECT_TRUE(c->is_backbone) << res << " C backbone";

        const auto* o = IupacAtomIdentity::Lookup(res, "O");
        ASSERT_NE(o, nullptr) << res << " O";
        EXPECT_EQ(o->element, Element::O) << res << " O element";
        EXPECT_TRUE(o->is_backbone) << res << " O backbone";
    }
}

TEST(IupacAtomIdentityTest, AllTwentyResiduesHaveAtoms) {
    const char* all20[] = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
        "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
        "THR", "TRP", "TYR", "VAL"
    };
    for (const char* res : all20) {
        const auto& atoms = IupacAtomIdentity::AtomsForResidue(res);
        EXPECT_GT(atoms.size(), 0u) << res << " has no atoms";
    }
}

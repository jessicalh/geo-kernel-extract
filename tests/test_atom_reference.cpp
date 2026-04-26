// Tests for AtomReference — typed cross-protein atom identity.
//
// Verifies: equality semantics; hash usability for unordered_map; that
// MakeAtomReference rebuilds the same reference for the same atom across
// two independent loads of the same PDB; that the BuildAtomReferenceMap
// pattern reliably matches atoms across two protein instances.

#include "TestEnvironment.h"
#include <gtest/gtest.h>

#include "AtomReference.h"
#include "Atom.h"
#include "Residue.h"
#include "Protein.h"
#include "PdbFileReader.h"

#include <unordered_map>

using namespace nmr;

namespace {
size_t FindResidue(const Protein& p, int seq_num, AminoAcid aa) {
    for (size_t ri = 0; ri < p.ResidueCount(); ++ri) {
        const Residue& r = p.ResidueAt(ri);
        if (r.sequence_number == seq_num && r.type == aa) return ri;
    }
    return SIZE_MAX;
}
size_t FindAtomInResidue(const Protein& p, size_t ri, const std::string& name) {
    const Residue& r = p.ResidueAt(ri);
    for (size_t ai : r.atom_indices) {
        if (p.AtomAt(ai).iupac_name == name) return ai;
    }
    return SIZE_MAX;
}
}  // namespace


TEST(AtomReferenceTest, EqualityRequiresAllFields) {
    AtomReference a{AminoAcid::PHE, 4, "A", "HB3"};
    AtomReference b{AminoAcid::PHE, 4, "A", "HB3"};
    AtomReference c{AminoAcid::PHE, 4, "A", "HB2"};   // different name
    AtomReference d{AminoAcid::PHE, 5, "A", "HB3"};   // different position
    AtomReference e{AminoAcid::TYR, 4, "A", "HB3"};   // different type
    AtomReference f{AminoAcid::PHE, 4, "B", "HB3"};   // different chain

    EXPECT_EQ(a, b);
    EXPECT_NE(a, c);
    EXPECT_NE(a, d);
    EXPECT_NE(a, e);
    EXPECT_NE(a, f);
}


TEST(AtomReferenceTest, HashAndUnorderedMap) {
    std::unordered_map<AtomReference, size_t> map;
    map[AtomReference{AminoAcid::PHE, 4, "", "HB3"}] = 42;
    map[AtomReference{AminoAcid::PHE, 4, "", "HB2"}] = 43;
    map[AtomReference{AminoAcid::TYR, 4, "", "HB3"}] = 44;

    auto it = map.find(AtomReference{AminoAcid::PHE, 4, "", "HB3"});
    ASSERT_NE(it, map.end());
    EXPECT_EQ(it->second, 42u);

    EXPECT_EQ(map.size(), 3u);
    EXPECT_EQ(map.count(AtomReference{AminoAcid::PHE, 4, "", "HB1"}), 0u);
}


TEST(AtomReferenceTest, MakeFromProteinAtomMatchesExpected) {
    auto result = BuildFromProtonatedPdb(test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(result.ok) << result.error;
    const Protein& p = *result.protein;

    // PHE-4 HB3 in 1ubq
    size_t phe_idx = FindResidue(p, 4, AminoAcid::PHE);
    ASSERT_NE(phe_idx, SIZE_MAX);
    size_t hb3_idx = FindAtomInResidue(p, phe_idx, "HB3");
    ASSERT_NE(hb3_idx, SIZE_MAX);

    AtomReference ref = MakeAtomReference(p, hb3_idx);
    EXPECT_EQ(ref.residue_type, AminoAcid::PHE);
    EXPECT_EQ(ref.residue_position, 4);
    EXPECT_EQ(ref.atom_name, IupacAtomName("HB3"));
}


TEST(AtomReferenceTest, MapAcrossTwoLoadsOfSameStructure) {
    // Loading the same PDB twice must produce AtomReferences that match
    // every atom one-to-one — the use case for cross-protein DFT compare.
    auto a = BuildFromProtonatedPdb(test::TestEnvironment::UbqProtonated());
    auto b = BuildFromProtonatedPdb(test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(a.ok); ASSERT_TRUE(b.ok);

    auto a_map = BuildAtomReferenceMap(*a.protein);
    auto b_map = BuildAtomReferenceMap(*b.protein);

    // Both maps should have one entry per atom — no collisions.
    EXPECT_EQ(a_map.size(), a.protein->AtomCount());
    EXPECT_EQ(b_map.size(), b.protein->AtomCount());

    // Every reference in a should map to the same reference in b.
    for (size_t ai = 0; ai < a.protein->AtomCount(); ++ai) {
        AtomReference ref = MakeAtomReference(*a.protein, ai);
        auto it = b_map.find(ref);
        ASSERT_NE(it, b_map.end()) << "missing in b: " << ref;
        const Atom& aa = a.protein->AtomAt(ai);
        const Atom& bb = b.protein->AtomAt(it->second);
        EXPECT_EQ(aa.element, bb.element);
        EXPECT_EQ(aa.iupac_name, bb.iupac_name);
    }
}


TEST(AtomReferenceTest, OrderingTotalAndConsistent) {
    // std::less should be a strict total order — same fields → equal,
    // different fields → ordered consistently.
    AtomReference r1{AminoAcid::PHE, 4, "", "HB2"};
    AtomReference r2{AminoAcid::PHE, 4, "", "HB3"};
    AtomReference r3{AminoAcid::PHE, 5, "", "HB2"};

    // Within same residue, atom_name orders.
    EXPECT_TRUE(r1 < r2);
    EXPECT_FALSE(r2 < r1);

    // Across residue positions, position orders.
    EXPECT_TRUE(r1 < r3);
    EXPECT_FALSE(r3 < r1);

    // Strict (no equal-comparing-less).
    EXPECT_FALSE(r1 < r1);
}

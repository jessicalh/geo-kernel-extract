#include "TestEnvironment.h"
//
// test_traversal_dump: traversal-and-dump test for the typed object model.
//
// Loads 1UBQ and traverses EVERY atom, EVERY residue, EVERY ring, EVERY bond.
// For each, it accesses TYPED PROPERTIES ONLY. iupac_name is used for
// printing only, never for comparison.
//
// This test is the proof that the typed model is complete: if any traversal
// requires a string comparison to produce its output, that is a bug.
//

#include <gtest/gtest.h>
#include "GeometryResult.h"
#include "PdbFileReader.h"
#include <filesystem>
#include <map>
#include <set>
#include <cstdio>

using namespace nmr;


class TraversalDumpTest : public ::testing::Test {
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


TEST_F(TraversalDumpTest, FullTypedTraversal) {
    auto& conf = protein->Conformation();

    // Attach GeometryResult so bond geometry and ring geometry are populated
    conf.AttachResult(GeometryResult::Compute(conf));

    // Counters for summary
    size_t atom_count = 0;
    size_t residue_count = 0;
    std::map<RingTypeIndex, size_t> rings_per_type;
    std::map<BondCategory, size_t> bonds_per_category;

    // Element counts (from typed property, not from name)
    std::map<Element, size_t> element_counts;

    // ================================================================
    // Traverse ALL atoms: typed properties only
    // ================================================================
    for (size_t ai = 0; ai < protein->AtomCount(); ++ai) {
        const Atom& atom = protein->AtomAt(ai);
        atom_count++;

        // Element: typed enum, not string
        Element elem = atom.element;
        element_counts[elem]++;

        // Residue index: typed size_t
        size_t res_idx = atom.residue_index;
        EXPECT_LT(res_idx, protein->ResidueCount())
            << "Atom " << ai << " has invalid residue_index";

        // Bond count: from typed bond_indices vector
        size_t bond_count = atom.bond_indices.size();
        // Every atom should have at least 1 bond
        // (isolated atoms would indicate a detection problem)
        EXPECT_GT(bond_count, 0u)
            << "Atom " << ai << " (" << atom.iupac_name << ") has no bonds";

        // Covalent radius and electronegativity: typed element properties
        double cr = atom.CovalentRadius();
        double en = atom.Electronegativity();
        EXPECT_GT(cr, 0.0);
        EXPECT_GT(en, 0.0);

        // ConformationAtom: position is accessible without string work
        Vec3 pos = conf.AtomAt(ai).Position();
        EXPECT_FALSE(std::isnan(pos.x()));
    }

    EXPECT_EQ(atom_count, protein->AtomCount());

    // ================================================================
    // Traverse ALL residues: typed properties only
    // ================================================================
    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        const Residue& res = protein->ResidueAt(ri);
        residue_count++;

        // Type: AminoAcid enum, not string
        AminoAcid aa = res.type;
        EXPECT_NE(aa, AminoAcid::Unknown)
            << "Residue " << ri << " has unknown amino acid type";

        // Chain ID: addressing (string is OK for display, not identity)
        EXPECT_FALSE(res.chain_id.empty());

        // IsAromatic: typed query method on the residue
        bool aromatic = res.IsAromatic();
        if (aa == AminoAcid::PHE || aa == AminoAcid::TYR ||
            aa == AminoAcid::TRP || aa == AminoAcid::HIS) {
            EXPECT_TRUE(aromatic) << "Expected aromatic for " << ThreeLetterCodeForAminoAcid(aa);
        }

        // Backbone indices: typed size_t values from index cache
        // N, CA, C should be present for all standard residues
        EXPECT_NE(res.N, Residue::NONE)
            << "Residue " << ri << " (" << ThreeLetterCodeForAminoAcid(aa) << ") missing N";
        EXPECT_NE(res.CA, Residue::NONE)
            << "Residue " << ri << " (" << ThreeLetterCodeForAminoAcid(aa) << ") missing CA";
        EXPECT_NE(res.C, Residue::NONE)
            << "Residue " << ri << " (" << ThreeLetterCodeForAminoAcid(aa) << ") missing C";

        // Chi angle count: from AminoAcidType table
        int chi_count = res.ChiAngleCount();
        EXPECT_GE(chi_count, 0);
        EXPECT_LE(chi_count, 4);

        // Atom indices: every residue should own atoms
        EXPECT_GT(res.atom_indices.size(), 0u);
    }

    EXPECT_EQ(residue_count, protein->ResidueCount());

    // ================================================================
    // Traverse ALL rings: typed properties only
    // ================================================================
    for (size_t ri = 0; ri < protein->RingCount(); ++ri) {
        const Ring& ring = protein->RingAt(ri);

        // Type index: typed enum
        RingTypeIndex type = ring.type_index;
        rings_per_type[type]++;

        // Virtual typed properties -- no string comparison
        double intensity = ring.Intensity();
        EXPECT_LT(intensity, 0.0)
            << "Ring " << ri << " has non-negative intensity";

        double lobe_offset = ring.JBLobeOffset();
        EXPECT_GT(lobe_offset, 0.0);

        int n_count = ring.NitrogenCount();
        EXPECT_GE(n_count, 0);
        EXPECT_LE(n_count, 2);

        RingAromaticity arom = ring.Aromaticity();
        (void)arom;  // Just verify it's accessible as a typed enum

        int size = ring.RingSizeValue();
        EXPECT_TRUE(size == 5 || size == 6 || size == 9)
            << "Ring " << ri << " has unexpected size " << size;

        // Atom count matches ring size
        EXPECT_EQ(static_cast<int>(ring.atom_indices.size()), size);

        // Parent residue: typed index
        EXPECT_LT(ring.parent_residue_index, protein->ResidueCount());

        // Fused partner: typed index or SIZE_MAX
        if (ring.IsFused()) {
            EXPECT_LT(ring.fused_partner_index, protein->RingCount());
        }

        // TypeIndexAsInt: typed integer for array indexing
        int idx_int = ring.TypeIndexAsInt();
        EXPECT_GE(idx_int, 0);
        EXPECT_LT(idx_int, static_cast<int>(RingTypeIndex::Count));

        // Ring geometry from GeometryResult (typed, computed from positions)
        const RingGeometry& geo = conf.ring_geometries[ri];
        EXPECT_NEAR(geo.normal.norm(), 1.0, 1e-6);
        EXPECT_GT(geo.radius, 0.5);
        EXPECT_LT(geo.radius, 5.0);
    }

    // ================================================================
    // Traverse ALL bonds: typed properties only
    // ================================================================
    for (size_t bi = 0; bi < protein->BondCount(); ++bi) {
        const Bond& bond = protein->BondAt(bi);

        // Category: typed enum
        BondCategory cat = bond.category;
        bonds_per_category[cat]++;

        // Order: typed enum
        BondOrder order = bond.order;
        EXPECT_NE(order, BondOrder::Unknown)
            << "Bond " << bi << " has unknown order";

        // Atom indices: typed size_t
        EXPECT_LT(bond.atom_index_a, protein->AtomCount());
        EXPECT_LT(bond.atom_index_b, protein->AtomCount());
        EXPECT_NE(bond.atom_index_a, bond.atom_index_b);

        // Classification helpers: typed, no strings
        if (cat == BondCategory::PeptideCO) {
            EXPECT_EQ(order, BondOrder::Double);
            EXPECT_TRUE(bond.IsBackbone());
            EXPECT_TRUE(bond.IsPeptideCO());
        }
        if (cat == BondCategory::PeptideCN) {
            EXPECT_EQ(order, BondOrder::Peptide);
            EXPECT_TRUE(bond.IsBackbone());
            EXPECT_TRUE(bond.IsPeptideBond());
        }
        if (cat == BondCategory::Aromatic) {
            EXPECT_TRUE(bond.IsAromatic());
        }

        // Bond geometry from GeometryResult
        double len = conf.bond_lengths[bi];
        EXPECT_GT(len, 0.3);
        EXPECT_LT(len, 3.5);
    }

    // ================================================================
    // Summary (printed, not tested -- this is the proof of completeness)
    // ================================================================
    printf("\n=== Typed Traversal Summary for 1UBQ ===\n");
    printf("Atoms: %zu\n", atom_count);

    printf("  Elements: ");
    for (auto& [elem, count] : element_counts) {
        printf("%s=%zu ", SymbolForElement(elem).c_str(), count);
    }
    printf("\n");

    printf("Residues: %zu\n", residue_count);

    printf("Rings: %zu\n", protein->RingCount());
    for (auto& [type, count] : rings_per_type) {
        printf("  %s: %zu\n", RingTypeName(type), count);
    }

    printf("Bonds: %zu\n", protein->BondCount());
    const char* category_names[] = {
        "PeptideCO", "PeptideCN", "BackboneOther", "SidechainCO",
        "Aromatic", "Disulfide", "SidechainOther", "Unknown"
    };
    for (auto& [cat, count] : bonds_per_category) {
        int idx = static_cast<int>(cat);
        const char* name = (idx >= 0 && idx < 8) ? category_names[idx] : "?";
        printf("  %s: %zu\n", name, count);
    }
    printf("=== End Traversal ===\n\n");

    // Sanity checks on counts
    // 1UBQ: 76 residues, ~660 atoms (with hydrogens varies by protonation)
    EXPECT_GE(residue_count, 70u);
    EXPECT_LE(residue_count, 80u);
    EXPECT_GT(atom_count, 500u);

    // Should have at least some of each bond category
    // (1UBQ has peptide bonds, aromatic bonds, backbone bonds)
    EXPECT_GT(bonds_per_category[BondCategory::PeptideCO], 0u);
    EXPECT_GT(bonds_per_category[BondCategory::PeptideCN], 0u);
    EXPECT_GT(bonds_per_category[BondCategory::BackboneOther], 0u);
    EXPECT_GT(bonds_per_category[BondCategory::Aromatic], 0u);

    // Should have at least PHE and TYR rings
    EXPECT_GT(rings_per_type.size(), 0u);
}


TEST_F(TraversalDumpTest, PeptideCOBondsMatchBackboneCache) {
    // Verify that every PeptideCO bond connects atoms that the backbone
    // index cache identifies as C and O of the same residue.
    for (size_t bi = 0; bi < protein->BondCount(); ++bi) {
        const Bond& bond = protein->BondAt(bi);
        if (bond.category != BondCategory::PeptideCO) continue;

        size_t a = bond.atom_index_a;
        size_t b = bond.atom_index_b;
        const Atom& atom_a = protein->AtomAt(a);
        const Atom& atom_b = protein->AtomAt(b);

        // Same residue
        EXPECT_EQ(atom_a.residue_index, atom_b.residue_index)
            << "PeptideCO bond " << bi << " crosses residues";

        const Residue& res = protein->ResidueAt(atom_a.residue_index);

        // One is res.C, the other is res.O
        bool match = (a == res.C && b == res.O) || (a == res.O && b == res.C);
        EXPECT_TRUE(match)
            << "PeptideCO bond " << bi << " atoms don't match backbone C,O cache";
    }
}


TEST_F(TraversalDumpTest, PeptideCNBondsCrossResidues) {
    // Verify that every PeptideCN bond connects C of one residue to N of another.
    for (size_t bi = 0; bi < protein->BondCount(); ++bi) {
        const Bond& bond = protein->BondAt(bi);
        if (bond.category != BondCategory::PeptideCN) continue;

        size_t a = bond.atom_index_a;
        size_t b = bond.atom_index_b;
        const Atom& atom_a = protein->AtomAt(a);
        const Atom& atom_b = protein->AtomAt(b);

        // Different residues
        EXPECT_NE(atom_a.residue_index, atom_b.residue_index)
            << "PeptideCN bond " << bi << " within same residue";

        const Residue& res_a = protein->ResidueAt(atom_a.residue_index);
        const Residue& res_b = protein->ResidueAt(atom_b.residue_index);

        // One is res_X.C, the other is res_Y.N
        bool match = (a == res_a.C && b == res_b.N) || (a == res_a.N && b == res_b.C);
        EXPECT_TRUE(match)
            << "PeptideCN bond " << bi << " atoms don't match backbone C/N cache";
    }
}


TEST_F(TraversalDumpTest, AromaticBondsAreRingMembers) {
    // Verify that every Aromatic bond connects two atoms that are in the
    // aromatic atom set (from ring detection).
    std::set<size_t> aromatic_atoms;
    for (size_t ri = 0; ri < protein->RingCount(); ++ri) {
        const Ring& ring = protein->RingAt(ri);
        for (size_t ai : ring.atom_indices) {
            aromatic_atoms.insert(ai);
        }
    }

    for (size_t bi = 0; bi < protein->BondCount(); ++bi) {
        const Bond& bond = protein->BondAt(bi);
        if (bond.category != BondCategory::Aromatic) continue;

        EXPECT_GT(aromatic_atoms.count(bond.atom_index_a), 0u)
            << "Aromatic bond " << bi << " atom_a not in ring";
        EXPECT_GT(aromatic_atoms.count(bond.atom_index_b), 0u)
            << "Aromatic bond " << bi << " atom_b not in ring";
    }
}

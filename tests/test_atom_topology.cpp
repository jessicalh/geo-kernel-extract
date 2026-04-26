// Smoke tests for the IUPAC AtomTopology stamping layer (Markley 1998).
//
// Loads 1ubq_protonated.pdb, walks specific residues, and verifies that
// Protein::StampAtomTopology produced the expected typed payload on Atom,
// that ring memberships are populated, partner atoms are linked, and
// residue context is correct.

#include "TestEnvironment.h"
#include <gtest/gtest.h>

#include "AtomTopology.h"
#include "Atom.h"
#include "Residue.h"
#include "Protein.h"
#include "PdbFileReader.h"
#include "AminoAcidType.h"

using namespace nmr;

namespace {

// Locate a residue by sequence number + amino-acid type. Returns SIZE_MAX
// if not found. (1ubq has unique sequence numbers within its single chain.)
size_t FindResidue(const Protein& p, int seq_num, AminoAcid aa) {
    for (size_t ri = 0; ri < p.ResidueCount(); ++ri) {
        const Residue& r = p.ResidueAt(ri);
        if (r.sequence_number == seq_num && r.type == aa) return ri;
    }
    return SIZE_MAX;
}

// Locate an atom in a residue by PDB name. Returns SIZE_MAX if not found.
size_t FindAtomInResidue(const Protein& p, size_t ri, const std::string& name) {
    const Residue& r = p.ResidueAt(ri);
    for (size_t ai : r.atom_indices) {
        if (p.AtomAt(ai).iupac_name == name) return ai;
    }
    return SIZE_MAX;
}

}  // namespace


TEST(AtomTopologyStampingTest, BackboneStampedAndAlpha) {
    auto result = BuildFromProtonatedPdb(test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(result.ok) << result.error;
    const Protein& p = *result.protein;

    // Check backbone N / C / O on the first residue (Met-1).
    size_t met_idx = FindResidue(p, 1, AminoAcid::MET);
    ASSERT_NE(met_idx, SIZE_MAX);

    size_t n_idx  = FindAtomInResidue(p, met_idx, "N");
    size_t ca_idx = FindAtomInResidue(p, met_idx, "CA");
    size_t c_idx  = FindAtomInResidue(p, met_idx, "C");
    size_t o_idx  = FindAtomInResidue(p, met_idx, "O");
    ASSERT_NE(n_idx, SIZE_MAX);
    ASSERT_NE(ca_idx, SIZE_MAX);
    ASSERT_NE(c_idx, SIZE_MAX);
    ASSERT_NE(o_idx, SIZE_MAX);

    EXPECT_TRUE(p.AtomAt(n_idx).topology.stamped);
    EXPECT_EQ(p.AtomAt(n_idx).topology.locant, Locant::Backbone);
    EXPECT_EQ(p.AtomAt(c_idx).topology.locant, Locant::Backbone);
    EXPECT_EQ(p.AtomAt(o_idx).topology.locant, Locant::Backbone);

    // CA / HA carry Locant::Alpha (start of the IUPAC Greek-letter chain).
    EXPECT_EQ(p.AtomAt(ca_idx).topology.locant, Locant::Alpha);
    size_t ha_idx = FindAtomInResidue(p, met_idx, "HA");
    if (ha_idx != SIZE_MAX) {
        EXPECT_EQ(p.AtomAt(ha_idx).topology.locant, Locant::Alpha);
    }
}


TEST(AtomTopologyStampingTest, PheBetaProchiralAndRing) {
    auto result = BuildFromProtonatedPdb(test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(result.ok) << result.error;
    const Protein& p = *result.protein;

    // Phe-4 in 1ubq.
    size_t phe_idx = FindResidue(p, 4, AminoAcid::PHE);
    ASSERT_NE(phe_idx, SIZE_MAX);

    size_t hb2_idx = FindAtomInResidue(p, phe_idx, "HB2");
    size_t hb3_idx = FindAtomInResidue(p, phe_idx, "HB3");
    size_t cd1_idx = FindAtomInResidue(p, phe_idx, "CD1");
    size_t cz_idx  = FindAtomInResidue(p, phe_idx, "CZ");
    ASSERT_NE(hb2_idx, SIZE_MAX);
    ASSERT_NE(hb3_idx, SIZE_MAX);
    ASSERT_NE(cd1_idx, SIZE_MAX);
    ASSERT_NE(cz_idx, SIZE_MAX);

    const Atom& hb2 = p.AtomAt(hb2_idx);
    const Atom& hb3 = p.AtomAt(hb3_idx);
    const Atom& cd1 = p.AtomAt(cd1_idx);
    const Atom& cz  = p.AtomAt(cz_idx);

    // β methylene: prochiral with QB pseudoatom, partner-linked.
    EXPECT_TRUE(hb2.topology.stamped);
    EXPECT_EQ(hb2.topology.locant, Locant::Beta);
    EXPECT_EQ(hb2.topology.diastereotopic_index, DiastereotopicIndex::Two);
    EXPECT_EQ(hb2.topology.prochiral_stereo, ProchiralStereo::ProS);
    EXPECT_EQ(hb2.topology.pseudoatom_class, PseudoatomClass::QB);
    EXPECT_EQ(hb2.partner_atom_index, hb3_idx);

    EXPECT_EQ(hb3.topology.diastereotopic_index, DiastereotopicIndex::Three);
    EXPECT_EQ(hb3.topology.prochiral_stereo, ProchiralStereo::ProR);
    EXPECT_EQ(hb3.partner_atom_index, hb2_idx);

    // Ring atoms: CD1 is δ branch 1, CZ is ζ. Both members of the Phe ring.
    EXPECT_EQ(cd1.topology.locant, Locant::Delta);
    EXPECT_EQ(cd1.topology.branch_index, BranchIndex::One);
    EXPECT_EQ(cd1.topology.ring_position, RingPosition::Member);
    EXPECT_FALSE(cd1.ring_indices.empty());

    EXPECT_EQ(cz.topology.locant, Locant::Zeta);
    EXPECT_EQ(cz.topology.ring_position, RingPosition::Member);
    EXPECT_FALSE(cz.ring_indices.empty());
}


TEST(AtomTopologyStampingTest, AlaMethylEquivalent) {
    auto result = BuildFromProtonatedPdb(test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(result.ok) << result.error;
    const Protein& p = *result.protein;

    // Ala-28 in 1ubq.
    size_t ala_idx = FindResidue(p, 28, AminoAcid::ALA);
    ASSERT_NE(ala_idx, SIZE_MAX);

    size_t hb1_idx = FindAtomInResidue(p, ala_idx, "HB1");
    size_t hb2_idx = FindAtomInResidue(p, ala_idx, "HB2");
    size_t hb3_idx = FindAtomInResidue(p, ala_idx, "HB3");
    ASSERT_NE(hb1_idx, SIZE_MAX);
    ASSERT_NE(hb2_idx, SIZE_MAX);
    ASSERT_NE(hb3_idx, SIZE_MAX);

    for (size_t ai : {hb1_idx, hb2_idx, hb3_idx}) {
        const Atom& a = p.AtomAt(ai);
        EXPECT_TRUE(a.topology.stamped);
        EXPECT_EQ(a.topology.locant, Locant::Beta);
        EXPECT_EQ(a.topology.prochiral_stereo, ProchiralStereo::Equivalent);
        EXPECT_EQ(a.topology.pseudoatom_class, PseudoatomClass::MB);
    }
}


TEST(AtomTopologyStampingTest, LysFourProchiralPairs) {
    auto result = BuildFromProtonatedPdb(test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(result.ok) << result.error;
    const Protein& p = *result.protein;

    // Lys-6 in 1ubq.
    size_t lys_idx = FindResidue(p, 6, AminoAcid::LYS);
    ASSERT_NE(lys_idx, SIZE_MAX);

    struct Pair { const char* h2; const char* h3; Locant loc; PseudoatomClass q;
                  ProchiralStereo s2; ProchiralStereo s3; };
    Pair pairs[] = {
        {"HB2", "HB3", Locant::Beta,    PseudoatomClass::QB, ProchiralStereo::ProS, ProchiralStereo::ProR},
        {"HG2", "HG3", Locant::Gamma,   PseudoatomClass::QG, ProchiralStereo::ProR, ProchiralStereo::ProS},
        {"HD2", "HD3", Locant::Delta,   PseudoatomClass::QD, ProchiralStereo::ProS, ProchiralStereo::ProR},
        {"HE2", "HE3", Locant::Epsilon, PseudoatomClass::QE, ProchiralStereo::ProR, ProchiralStereo::ProS},
    };

    for (const Pair& pair : pairs) {
        size_t i2 = FindAtomInResidue(p, lys_idx, pair.h2);
        size_t i3 = FindAtomInResidue(p, lys_idx, pair.h3);
        ASSERT_NE(i2, SIZE_MAX) << pair.h2;
        ASSERT_NE(i3, SIZE_MAX) << pair.h3;

        EXPECT_EQ(p.AtomAt(i2).topology.locant, pair.loc) << pair.h2;
        EXPECT_EQ(p.AtomAt(i2).topology.diastereotopic_index, DiastereotopicIndex::Two);
        EXPECT_EQ(p.AtomAt(i2).topology.prochiral_stereo, pair.s2) << pair.h2;
        EXPECT_EQ(p.AtomAt(i2).topology.pseudoatom_class, pair.q);
        EXPECT_EQ(p.AtomAt(i2).partner_atom_index, i3);

        EXPECT_EQ(p.AtomAt(i3).topology.diastereotopic_index, DiastereotopicIndex::Three);
        EXPECT_EQ(p.AtomAt(i3).topology.prochiral_stereo, pair.s3) << pair.h3;
        EXPECT_EQ(p.AtomAt(i3).partner_atom_index, i2);
    }
}


TEST(AtomTopologyStampingTest, ResidueContextChainEnds) {
    auto result = BuildFromProtonatedPdb(test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(result.ok) << result.error;
    const Protein& p = *result.protein;
    ASSERT_GT(p.ResidueCount(), 0u);

    // First residue: N-terminal, no prev.
    const Residue& first = p.ResidueAt(0);
    EXPECT_TRUE(first.is_n_terminal);
    EXPECT_EQ(first.prev_residue_type, AminoAcid::Unknown);
    if (p.ResidueCount() > 1) {
        EXPECT_EQ(first.next_residue_type, p.ResidueAt(1).type);
    }

    // Last residue: C-terminal, no next.
    size_t last = p.ResidueCount() - 1;
    const Residue& last_r = p.ResidueAt(last);
    EXPECT_TRUE(last_r.is_c_terminal);
    EXPECT_EQ(last_r.next_residue_type, AminoAcid::Unknown);
    if (last > 0) {
        EXPECT_EQ(last_r.prev_residue_type, p.ResidueAt(last - 1).type);
    }
}


// Coverage test: every entry in every standard amino acid template has
// `topology.stamped == true`. Catches a missing residue in the table or a
// missed atom inside a residue. Independent of any protein load — language-
// level over the static AMINO_ACID_TYPES.
TEST(AtomTopologyStampingTest, AllStandardAminoAcidsFullyPopulated) {
    int total_atoms = 0;
    int unstamped = 0;
    std::map<std::string, std::vector<std::string>> gaps;

    for (const AminoAcidType& aatype : AllAminoAcidTypes()) {
        if (aatype.index == AminoAcid::Unknown) continue;
        for (const AminoAcidAtom& a : aatype.atoms) {
            ++total_atoms;
            if (!a.topology.stamped) {
                ++unstamped;
                gaps[aatype.three_letter_code].push_back(a.name.AsString());
            }
        }
    }

    if (unstamped > 0) {
        for (const auto& kv : gaps) {
            std::string atoms;
            for (const auto& nm : kv.second) {
                if (!atoms.empty()) atoms += ", ";
                atoms += nm;
            }
            ADD_FAILURE() << kv.first << " missing topology for: " << atoms;
        }
    }

    EXPECT_EQ(unstamped, 0);
    EXPECT_GT(total_atoms, 200);  // sanity floor — sum across 20 AAs
}


// Coverage test: every atom of every residue in 1ubq stamps cleanly. Catches
// real-world misses from reduce-protonated input that wouldn't show up in the
// per-AA table walk above (e.g., a hydrogen produced by reduce that we forgot
// to put in the template).
TEST(AtomTopologyStampingTest, EveryAtomInUbqStamped) {
    auto result = BuildFromProtonatedPdb(test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(result.ok) << result.error;
    const Protein& p = *result.protein;

    int unstamped = 0;
    std::map<std::string, std::set<std::string>> first_gaps;  // residue → atom names
    for (size_t ai = 0; ai < p.AtomCount(); ++ai) {
        const Atom& a = p.AtomAt(ai);
        if (!a.topology.stamped) {
            ++unstamped;
            const Residue& r = p.ResidueAt(a.residue_index);
            const auto& aatype = r.AminoAcidInfo();
            first_gaps[aatype.three_letter_code].insert(a.iupac_name.AsString());
        }
    }

    if (unstamped > 0) {
        for (const auto& kv : first_gaps) {
            std::string atoms;
            for (const auto& nm : kv.second) {
                if (!atoms.empty()) atoms += ", ";
                atoms += nm;
            }
            ADD_FAILURE() << kv.first << " in 1ubq has unstamped atoms: " << atoms;
        }
    }

    EXPECT_EQ(unstamped, 0)
        << "Total unstamped atoms in 1ubq: " << unstamped << " of " << p.AtomCount();
}


TEST(AtomTopologyStampingTest, ChiPositionsForLysSidechain) {
    auto result = BuildFromProtonatedPdb(test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(result.ok) << result.error;
    const Protein& p = *result.protein;

    // Lys χ angles: χ1 N-CA-CB-CG, χ2 CA-CB-CG-CD, χ3 CB-CG-CD-CE, χ4 CG-CD-CE-NZ.
    size_t lys_idx = FindResidue(p, 6, AminoAcid::LYS);
    ASSERT_NE(lys_idx, SIZE_MAX);

    size_t cb_idx = FindAtomInResidue(p, lys_idx, "CB");
    ASSERT_NE(cb_idx, SIZE_MAX);

    // CB participates in χ1 (position 2) and χ2 (position 1).
    const Atom& cb = p.AtomAt(cb_idx);
    EXPECT_EQ(cb.topology.chi_position[0], 2);   // χ1 position 2
    EXPECT_EQ(cb.topology.chi_position[1], 1);   // χ2 position 1
    EXPECT_EQ(cb.topology.chi_position[2], 0);   // χ3 position 0
    EXPECT_EQ(cb.topology.chi_position[3], -1);  // not in χ4
}

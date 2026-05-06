// tests/topology/test_legacy_amber_semantic_integration.cpp
//
// Corpus / audit test for the LegacyAmberTopology semantic-substrate
// integration. Loads representative AMBER-prepared protein PDBs through
// the canonical PdbFileReader path, runs Protein::FinalizeConstruction,
// and asserts:
//
//   1. Coverage: every residue with canonical AMBER atom naming has a
//      non-default-constructed semantic record on every atom.
//   2. No residue triggers the fail-loudly path (the std::abort in
//      ComposeAtomSemantic). Tested implicitly: if any standard atom
//      goes unresolved by LookupBy, the process aborts before this
//      test reads a result.
//   3. Spot-checks against the 2026-05-05 substrate's known chemistry:
//        * NTERM_CHARGED N has formal_charge = +1.
//        * NTERM_CHARGED H1/H2/H3 carry PolarHKind::AmmoniumNH.
//        * CTERM_DEPROTONATED C and O sit in PlanarGroupKind::Carboxylate.
//        * GLY HA2 has ProchiralStereo::ProR.
//   4. Across the standard 20, every backbone amide H carries
//      PolarHKind::BackboneAmide (where present; PRO has no backbone H).
//   5. 1Z9B fleet input.pdb — the LYS-labelled-LYN-chemistry case —
//      now PASSES via Bundle B's post-protonation re-canonicalisation.
//      Spot-checks confirm every LYS residue with HZ1+HZ2 has variant
//      0 (LYN) and the canonical HZ2/HZ3 names; substrate composition
//      completes without FATAL.
//   6. Typed CacheResidueBackboneIndices produces the expected res.HA
//      slot for Glycine (HA2), the absent-by-design res.H for Pro,
//      and proper res.CB / res.HA on a representative non-Gly residue.
//   7. HIS variant coverage (HID/HIE/HIP) — substrate row chemistry
//      checked where the fixture has them.
//
// Per spec/plan/topology-encoding-dependencies-2026-05-05.md §H.5.

#include <gtest/gtest.h>

#include <filesystem>
#include <set>
#include <string>

#include "TestEnvironment.h"
#include "Protein.h"
#include "PdbFileReader.h"
#include "AminoAcidType.h"
#include "Atom.h"
#include "Residue.h"
#include "LegacyAmberTopology.h"
#include "SemanticEnums.h"
#include "generated/LegacyAmberSemanticTables.h"

namespace {

using nmr::AminoAcid;
using nmr::AminoAcidType;
using nmr::Atom;
using nmr::AtomMechanicalIdentity;
using nmr::AtomSemanticTable;
using nmr::BackboneRole;
using nmr::BranchAddress;
using nmr::DiastereotopicIndex;
using nmr::Element;
using nmr::Locant;
using nmr::PlanarGroupKind;
using nmr::PolarHKind;
using nmr::ProchiralStereo;
using nmr::Protein;
using nmr::Residue;
using nmr::ResidueTerminalState;
using nmr::RingSystemKind;
using nmr::RingPositionLabel;

namespace gen = nmr::topology_generated;

struct CoverageStats {
    int residues_checked = 0;
    int atoms_checked = 0;
    int backbone_amide_h_seen = 0;
    int backbone_amide_h_with_correct_polar_h = 0;
    int nterm_residues_seen = 0;
    int cterm_residues_seen = 0;
    int gly_ha2_seen = 0;
    int gly_ha2_with_pro_r = 0;
    int his_residues_seen = 0;
    int lyn_residues_seen = 0;
};

void AuditProtein(const Protein& protein, CoverageStats& stats) {
    ASSERT_TRUE(protein.HasTopology()) << "Protein must have a topology";
    const auto& topology = protein.LegacyAmber();
    ASSERT_TRUE(topology.HasAtomSemantic())
        << "LegacyAmberTopology must have a populated semantic store after "
           "FinalizeConstruction on a real AMBER-prepared protein";
    ASSERT_EQ(topology.AtomSemantic().size(), protein.AtomCount())
        << "Semantic store size must equal atom count";

    for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        if (res.type == AminoAcid::Unknown) continue;

        ++stats.residues_checked;
        if (res.terminal_state == ResidueTerminalState::NTerminus ||
            res.terminal_state == ResidueTerminalState::NAndCTerminus) {
            ++stats.nterm_residues_seen;
        }
        if (res.terminal_state == ResidueTerminalState::CTerminus ||
            res.terminal_state == ResidueTerminalState::NAndCTerminus) {
            ++stats.cterm_residues_seen;
        }
        if (res.type == AminoAcid::HIS) ++stats.his_residues_seen;
        if (res.type == AminoAcid::LYS &&
            res.protonation_variant_index == 0 &&
            res.protonation_state_resolved) {
            ++stats.lyn_residues_seen;
        }

        for (size_t ai : res.atom_indices) {
            ++stats.atoms_checked;
            const AtomSemanticTable& row = topology.SemanticAt(ai);
            const Atom& atom = protein.AtomAt(ai);

            // Element must be populated (substrate never emits Unknown
            // element; the FATAL+abort path catches misses upstream).
            EXPECT_NE(Element::Unknown, row.element)
                << "Default-constructed semantic row on atom '"
                << atom.pdb_atom_name << "' in "
                << res.AminoAcidInfo().three_letter_code
                << " seq " << res.sequence_number;

            // The element field on the row must agree with the runtime
            // atom's element. AtomMechanicalIdentity uses element as
            // the first lookup-key field; mismatch means a wrong row
            // was selected.
            EXPECT_EQ(atom.element, row.element)
                << "Element mismatch on atom '" << atom.pdb_atom_name
                << "' in " << res.AminoAcidInfo().three_letter_code
                << " seq " << res.sequence_number;

            if (atom.pdb_atom_name == "H" || atom.pdb_atom_name == "HN") {
                ++stats.backbone_amide_h_seen;
                if (row.polar_h == PolarHKind::BackboneAmide) {
                    ++stats.backbone_amide_h_with_correct_polar_h;
                }
            }

            if (res.type == AminoAcid::GLY && atom.pdb_atom_name == "HA2") {
                ++stats.gly_ha2_seen;
                if (row.prochiral == ProchiralStereo::ProR) {
                    ++stats.gly_ha2_with_pro_r;
                }
            }
        }
    }
}

size_t FindAtomIndex(const Protein& protein, size_t residue_index,
                     const std::string& target) {
    const Residue& res = protein.ResidueAt(residue_index);
    for (size_t ai : res.atom_indices) {
        if (protein.AtomAt(ai).pdb_atom_name == target) return ai;
    }
    return SIZE_MAX;
}

bool ResidueHasAtom(const Protein& protein, size_t ri, const std::string& name) {
    return FindAtomIndex(protein, ri, name) != SIZE_MAX;
}

// ============================================================================
// 1UBQ — full-coverage audit
// ============================================================================

TEST(LegacyAmberSemanticIntegration, UbqProtonatedFullCoverage) {
    if (!std::filesystem::exists(nmr::test::TestEnvironment::UbqProtonated())) {
        GTEST_SKIP() << "1UBQ PDB not found at "
                     << nmr::test::TestEnvironment::UbqProtonated();
    }
    auto result = nmr::BuildFromProtonatedPdb(
        nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(result.Ok()) << result.error;
    const Protein& protein = *result.protein;

    CoverageStats stats;
    AuditProtein(protein, stats);

    EXPECT_GT(stats.residues_checked, 70)
        << "1UBQ should have ~76 substrate-covered residues";
    EXPECT_GT(stats.atoms_checked, 1000)
        << "1UBQ has ~1231 atoms; expect substantial coverage";

    EXPECT_GT(stats.backbone_amide_h_seen, 0);
    EXPECT_EQ(stats.backbone_amide_h_seen,
              stats.backbone_amide_h_with_correct_polar_h)
        << "Every backbone amide H seen on 1UBQ must carry "
           "PolarHKind::BackboneAmide";
}


// ============================================================================
// 1UBQ — NTERM / CTERM cap composition spot-checks
// ============================================================================

TEST(LegacyAmberSemanticIntegration, UbqNTermAndCTermSpotChecks) {
    if (!std::filesystem::exists(nmr::test::TestEnvironment::UbqProtonated())) {
        GTEST_SKIP() << "1UBQ PDB not found";
    }
    auto result = nmr::BuildFromProtonatedPdb(
        nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(result.Ok()) << result.error;
    const Protein& protein = *result.protein;
    const auto& topology = protein.LegacyAmber();

    size_t n_term_idx = SIZE_MAX;
    size_t c_term_idx = SIZE_MAX;
    for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        if (n_term_idx == SIZE_MAX &&
            (res.terminal_state == ResidueTerminalState::NTerminus ||
             res.terminal_state == ResidueTerminalState::NAndCTerminus)) {
            n_term_idx = ri;
        }
        if (res.terminal_state == ResidueTerminalState::CTerminus ||
            res.terminal_state == ResidueTerminalState::NAndCTerminus) {
            c_term_idx = ri;
        }
    }
    ASSERT_NE(SIZE_MAX, n_term_idx);
    ASSERT_NE(SIZE_MAX, c_term_idx);

    // NTERM N: substrate composition flips formal_charge -> +1 via
    // LookupCap(NtermCharged) -> ApplyCapDelta on the chain N's row.
    const size_t n_atom = FindAtomIndex(protein, n_term_idx, "N");
    ASSERT_NE(SIZE_MAX, n_atom);
    EXPECT_EQ(1, topology.SemanticAt(n_atom).formal_charge)
        << "NTERM_CHARGED N must have formal_charge = +1 after cap composition";

    // NTERM H1/H2/H3: cap-only atoms whose substrate row carries
    // PolarHKind::AmmoniumNH. NTERM_CHARGED is the AMBER ff14SB default.
    int nterm_h_count = 0;
    int nterm_h_with_ammonium_polar_h = 0;
    for (const std::string& name : {"H1", "H2", "H3"}) {
        const size_t ai = FindAtomIndex(protein, n_term_idx, name);
        if (ai == SIZE_MAX) continue;
        ++nterm_h_count;
        if (topology.SemanticAt(ai).polar_h == PolarHKind::AmmoniumNH) {
            ++nterm_h_with_ammonium_polar_h;
        }
    }
    EXPECT_GE(nterm_h_count, 2);
    EXPECT_EQ(nterm_h_count, nterm_h_with_ammonium_polar_h);

    // CTERM C and O: chain entries' planar_group is PeptideAmide;
    // ApplyCapDelta(CtermDeprotonated) overlays Carboxylate.
    const size_t c_atom = FindAtomIndex(protein, c_term_idx, "C");
    const size_t o_atom = FindAtomIndex(protein, c_term_idx, "O");
    ASSERT_NE(SIZE_MAX, c_atom);
    ASSERT_NE(SIZE_MAX, o_atom);
    EXPECT_EQ(PlanarGroupKind::Carboxylate,
              topology.SemanticAt(c_atom).planar_group);
    EXPECT_EQ(PlanarGroupKind::Carboxylate,
              topology.SemanticAt(o_atom).planar_group);

    // CTERM OXT: cap-only; formal_charge = -1 in CtermDeprotonated.
    const size_t oxt = FindAtomIndex(protein, c_term_idx, "OXT");
    if (oxt != SIZE_MAX) {
        EXPECT_EQ(-1, topology.SemanticAt(oxt).formal_charge);
        EXPECT_EQ(PlanarGroupKind::Carboxylate,
                  topology.SemanticAt(oxt).planar_group);
    }
}


// ============================================================================
// 1Z9B fleet input — LYN re-canonicalisation un-skip
// ============================================================================

TEST(LegacyAmberSemanticIntegration, FleetAmberZ9bLynRecanonicalisation) {
    const std::string fixture =
        nmr::test::TestEnvironment::FleetAmberData() +
        "/1Z9B_6577/prep_run_20260501T143103Z/input.pdb";
    if (!std::filesystem::exists(fixture)) {
        GTEST_SKIP() << "fleet_amber 1Z9B input.pdb fixture not found at "
                     << fixture;
    }
    auto result = nmr::BuildFromProtonatedPdb(fixture);
    ASSERT_TRUE(result.Ok())
        << "Bundle B post-protonation re-canonicalisation should make "
           "1Z9B input.pdb load. Error: " << result.error;
    const Protein& protein = *result.protein;

    CoverageStats stats;
    AuditProtein(protein, stats);

    EXPECT_GT(stats.residues_checked, 50);
    EXPECT_GT(stats.atoms_checked, 800);
    EXPECT_GT(stats.backbone_amide_h_seen, 30);
    EXPECT_EQ(stats.backbone_amide_h_seen,
              stats.backbone_amide_h_with_correct_polar_h);

    // 1Z9B has 8 LYS residues, all of which carry HZ1+HZ2 only (LYN
    // chemistry mislabelled as LYS in the input). Bundle B's post-
    // protonation re-canonicalisation must rewrite HZ1->HZ2 and
    // HZ2->HZ3 after variant_idx=0 (LYN) is set.
    int lys_seen = 0;
    int lys_with_canonical_hz23 = 0;
    int lys_lyn_variant = 0;
    for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        if (res.type != AminoAcid::LYS) continue;
        ++lys_seen;
        if (res.protonation_variant_index == 0 &&
            res.protonation_state_resolved) {
            ++lys_lyn_variant;
        }
        const bool has_hz2 = ResidueHasAtom(protein, ri, "HZ2");
        const bool has_hz3 = ResidueHasAtom(protein, ri, "HZ3");
        const bool has_hz1 = ResidueHasAtom(protein, ri, "HZ1");
        if (has_hz2 && has_hz3 && !has_hz1) ++lys_with_canonical_hz23;
    }
    EXPECT_GT(lys_seen, 0);
    EXPECT_EQ(lys_seen, lys_lyn_variant)
        << "Every 1Z9B LYS residue with HZ1+HZ2 only should resolve to "
           "variant 0 (LYN)";
    EXPECT_EQ(lys_seen, lys_with_canonical_hz23)
        << "Every LYN-resolved residue must have canonical HZ2/HZ3 atom "
           "names after Bundle B's post-protonation re-canonicalisation";

    // Spot-check the LYN substrate row chemistry: HZ2/HZ3 should carry
    // PolarHKind::AmineNH (not AmmoniumNH) per the 2026-05-05 substrate
    // dependencies §A.1; NZ formal_charge = 0.
    int hz_count = 0;
    int hz_with_amine_nh = 0;
    int nz_zero_charge = 0;
    int nz_total = 0;
    for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        if (res.type != AminoAcid::LYS) continue;
        if (res.protonation_variant_index != 0) continue;
        for (const std::string& name : {"HZ2", "HZ3"}) {
            const size_t ai = FindAtomIndex(protein, ri, name);
            if (ai == SIZE_MAX) continue;
            ++hz_count;
            if (protein.LegacyAmber().SemanticAt(ai).polar_h ==
                PolarHKind::AmineNH) {
                ++hz_with_amine_nh;
            }
        }
        const size_t nz_ai = FindAtomIndex(protein, ri, "NZ");
        if (nz_ai != SIZE_MAX) {
            ++nz_total;
            if (protein.LegacyAmber().SemanticAt(nz_ai).formal_charge == 0) {
                ++nz_zero_charge;
            }
        }
    }
    EXPECT_GT(hz_count, 0);
    EXPECT_EQ(hz_count, hz_with_amine_nh)
        << "LYN HZ2/HZ3 must carry PolarHKind::AmineNH";
    EXPECT_GT(nz_total, 0);
    EXPECT_EQ(nz_total, nz_zero_charge)
        << "LYN NZ must have formal_charge = 0";
}


// ============================================================================
// 1P9J fleet input — full coverage + GLY HA2 check
// ============================================================================

TEST(LegacyAmberSemanticIntegration, FleetAmberP9jCoverage) {
    const std::string fixture =
        nmr::test::TestEnvironment::FleetAmberData() +
        "/1P9J_5801/prep_run_20260501T141627Z/input.pdb";
    if (!std::filesystem::exists(fixture)) {
        GTEST_SKIP() << "fleet_amber 1P9J input.pdb fixture not found at "
                     << fixture;
    }
    auto result = nmr::BuildFromProtonatedPdb(fixture);
    ASSERT_TRUE(result.Ok()) << result.error;
    const Protein& protein = *result.protein;

    CoverageStats stats;
    AuditProtein(protein, stats);
    EXPECT_GT(stats.residues_checked, 30);
    EXPECT_GT(stats.atoms_checked, 500);
    EXPECT_EQ(stats.backbone_amide_h_seen,
              stats.backbone_amide_h_with_correct_polar_h);

    if (stats.gly_ha2_seen > 0) {
        EXPECT_EQ(stats.gly_ha2_seen, stats.gly_ha2_with_pro_r)
            << "Every GLY HA2 must carry ProchiralStereo::ProR per "
               "Markley Fig 1 alternation inversion for glycine";
    }
}


// ============================================================================
// HIS variant chemistry — exercise HID/HIE/HIP rows on whichever
// fixture has them
// ============================================================================

TEST(LegacyAmberSemanticIntegration, HisVariantChemistry) {
    if (!std::filesystem::exists(nmr::test::TestEnvironment::UbqProtonated())) {
        GTEST_SKIP() << "1UBQ PDB not found";
    }
    auto result = nmr::BuildFromProtonatedPdb(
        nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(result.Ok()) << result.error;
    const Protein& protein = *result.protein;
    const auto& topology = protein.LegacyAmber();

    int his_seen = 0;
    int his_hid_seen = 0;
    int his_hie_seen = 0;
    int his_hip_seen = 0;
    for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        if (res.type != AminoAcid::HIS) continue;
        ++his_seen;
        ASSERT_TRUE(res.protonation_state_resolved);
        ASSERT_GE(res.protonation_variant_index, 0);
        const int variant = res.protonation_variant_index;
        if (variant == 0) ++his_hid_seen;
        if (variant == 1) ++his_hie_seen;
        if (variant == 2) ++his_hip_seen;

        const bool has_hd1 = ResidueHasAtom(protein, ri, "HD1");
        const bool has_he2 = ResidueHasAtom(protein, ri, "HE2");
        if (variant == 0) {  // HID: Nδ1 protonated, Nε2 not.
            EXPECT_TRUE(has_hd1) << "HID must have HD1";
            EXPECT_FALSE(has_he2) << "HID must NOT have HE2";
        } else if (variant == 1) {  // HIE: Nε2 protonated, Nδ1 not.
            EXPECT_FALSE(has_hd1) << "HIE must NOT have HD1";
            EXPECT_TRUE(has_he2) << "HIE must have HE2";
        } else if (variant == 2) {  // HIP: both protonated.
            EXPECT_TRUE(has_hd1);
            EXPECT_TRUE(has_he2);
        }

        // For HIP, Nε2 carries +1 formal charge.
        if (variant == 2) {
            const size_t ne2 = FindAtomIndex(protein, ri, "NE2");
            if (ne2 != SIZE_MAX) {
                EXPECT_EQ(1, topology.SemanticAt(ne2).formal_charge)
                    << "HIP Nε2 must have formal_charge = +1";
            }
        }
    }
    if (his_seen == 0) GTEST_SKIP() << "1UBQ has no HIS residues";
    EXPECT_GT(his_seen, 0);
}


// ============================================================================
// Typed CacheResidueBackboneIndices — Pro res.H absent, Gly res.HA = HA2
// ============================================================================

TEST(LegacyAmberSemanticIntegration, TypedBackboneCachePro) {
    if (!std::filesystem::exists(nmr::test::TestEnvironment::UbqProtonated())) {
        GTEST_SKIP() << "1UBQ PDB not found";
    }
    auto result = nmr::BuildFromProtonatedPdb(
        nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(result.Ok()) << result.error;
    const Protein& protein = *result.protein;

    int pro_seen = 0;
    int pro_h_none = 0;
    int pro_with_ca = 0;
    for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        if (res.type != AminoAcid::PRO) continue;
        ++pro_seen;
        if (res.H == Residue::NONE) ++pro_h_none;
        if (res.CA != Residue::NONE) ++pro_with_ca;
    }
    if (pro_seen == 0) GTEST_SKIP() << "No Pro in fixture";
    EXPECT_EQ(pro_seen, pro_h_none)
        << "Every Pro must have res.H == NONE (PRO chain table drops "
           "the backbone amide H per substrate dependencies §H.10)";
    EXPECT_EQ(pro_seen, pro_with_ca)
        << "Every Pro must still cache res.CA";
}


TEST(LegacyAmberSemanticIntegration, TypedBackboneCacheNonPro) {
    if (!std::filesystem::exists(nmr::test::TestEnvironment::UbqProtonated())) {
        GTEST_SKIP() << "1UBQ PDB not found";
    }
    auto result = nmr::BuildFromProtonatedPdb(
        nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(result.Ok()) << result.error;
    const Protein& protein = *result.protein;
    const auto& topology = protein.LegacyAmber();

    // Pick the first non-Pro internal residue and verify backbone cache
    // matches typed substrate.
    bool checked_one = false;
    for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        if (res.type == AminoAcid::PRO) continue;
        if (res.type == AminoAcid::GLY) continue;
        if (res.type == AminoAcid::Unknown) continue;
        if (res.terminal_state != ResidueTerminalState::Internal) continue;
        EXPECT_NE(Residue::NONE, res.N);
        EXPECT_NE(Residue::NONE, res.CA);
        EXPECT_NE(Residue::NONE, res.C);
        EXPECT_NE(Residue::NONE, res.O);
        EXPECT_NE(Residue::NONE, res.H);
        EXPECT_NE(Residue::NONE, res.HA);
        EXPECT_NE(Residue::NONE, res.CB);

        EXPECT_EQ(BackboneRole::AlphaHydrogen,
                  topology.SemanticAt(res.HA).backbone_role);
        EXPECT_EQ(BackboneRole::AmideHydrogen,
                  topology.SemanticAt(res.H).backbone_role);
        // CB carries Locant::Beta + Element::C.
        const AtomSemanticTable& cb = topology.SemanticAt(res.CB);
        EXPECT_EQ(Element::C, cb.element);
        EXPECT_EQ(Locant::Beta, cb.locant);
        checked_one = true;
        break;
    }
    EXPECT_TRUE(checked_one) << "Should have found a non-Gly/non-Pro internal residue";
}


TEST(LegacyAmberSemanticIntegration, TypedBackboneCacheGlyHa) {
    // Look at a fixture that has Gly. 1P9J and 1Z9B both have GLY.
    const std::string fixture =
        nmr::test::TestEnvironment::FleetAmberData() +
        "/1P9J_5801/prep_run_20260501T141627Z/input.pdb";
    if (!std::filesystem::exists(fixture)) {
        GTEST_SKIP() << "fleet_amber 1P9J input.pdb not found";
    }
    auto result = nmr::BuildFromProtonatedPdb(fixture);
    ASSERT_TRUE(result.Ok()) << result.error;
    const Protein& protein = *result.protein;
    const auto& topology = protein.LegacyAmber();

    int gly_seen = 0;
    int gly_with_ha2 = 0;
    int gly_cb_none = 0;
    for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        if (res.type != AminoAcid::GLY) continue;
        ++gly_seen;
        if (res.CB == Residue::NONE) ++gly_cb_none;
        if (res.HA != Residue::NONE) {
            // Gly's HA should be the HA2 atom — Locant::Alpha +
            // DiastereotopicIndex::Position2.
            const AtomSemanticTable& ha = topology.SemanticAt(res.HA);
            if (ha.locant == Locant::Alpha &&
                ha.di_index == DiastereotopicIndex::Position2) {
                ++gly_with_ha2;
            }
        }
    }
    if (gly_seen == 0) GTEST_SKIP() << "No Gly in fixture";
    EXPECT_EQ(gly_seen, gly_with_ha2)
        << "Every Gly's res.HA must point at HA2 (Locant::Alpha + Position2)";
    EXPECT_EQ(gly_seen, gly_cb_none)
        << "Every Gly's res.CB must be NONE (no CB in Gly substrate table)";
}

}  // namespace

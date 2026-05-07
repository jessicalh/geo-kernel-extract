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


// ============================================================================
// PRO N-terminus pyrrolidine ring preservation — codex Finding 1 regression
//
// Cap-table N entries carry RingSystemKind::NotInRing as a placeholder
// (cap atoms are not in any ring). For PRO at the N-terminus, the chain
// N row carries RingSystemKind::Pyrrolidine_Pro — the chain row's ring
// membership is the truth. ApplyCapDelta must NOT clobber ring_position
// during the backbone-cap overlay; it overlays only the chemistry-level
// fields (planar_group, polar_h, formal_charge, pseudoatom).
//
// The fleet_amber fixtures (1Z9B, 1P9J) do not start with PRO, so this
// regression is exercised on a synthetic fixture built in-test.
//
// Per spec/plan/topology-encoding-dependencies-2026-05-05.md §H.5
// (post-2026-05-06 codex Finding 1 amendment).
// ============================================================================

namespace {

// Build a minimal PRO-only protein. PRO at residue 1 makes it
// NTerminus + CTerminus (single chain, single residue → NAndCTerminus,
// per Protein::ResolveResidueTerminalStates). This drives both NTERM
// and CTERM cap composition; we assert chain N's ring membership stays.
//
// Atom set: PRO chain table (14 atoms: N, CA, C, O, HA, CB, HB2, HB3,
// CG, HG2, HG3, CD, HD2, HD3) plus N-terminal cap H1/H2/H3.
//
// Positions are zeros — the substrate composition path needs only typed
// atom-name plus residue type plus terminal_state to look up the right
// rows; bond detection is incidental.
std::unique_ptr<nmr::Protein> MakeNTermProProtein() {
    auto pp = std::make_unique<nmr::Protein>();
    nmr::Protein& p = *pp;

    nmr::Residue pro;
    pro.type = nmr::AminoAcid::PRO;
    pro.sequence_number = 1;
    pro.chain_id = "A";
    const size_t ri = p.AddResidue(pro);

    // PRO chain atoms + NTERM cap Hs.
    const std::pair<const char*, nmr::Element> atom_set[] = {
        {"N",   nmr::Element::N},
        {"CA",  nmr::Element::C},
        {"C",   nmr::Element::C},
        {"O",   nmr::Element::O},
        {"HA",  nmr::Element::H},
        {"CB",  nmr::Element::C},
        {"HB2", nmr::Element::H},
        {"HB3", nmr::Element::H},
        {"CG",  nmr::Element::C},
        {"HG2", nmr::Element::H},
        {"HG3", nmr::Element::H},
        {"CD",  nmr::Element::C},
        {"HD2", nmr::Element::H},
        {"HD3", nmr::Element::H},
        {"H1",  nmr::Element::H},
        {"H2",  nmr::Element::H},
        {"H3",  nmr::Element::H},
    };

    std::vector<size_t> atom_indices;
    atom_indices.reserve(sizeof(atom_set) / sizeof(atom_set[0]));
    for (const auto& kv : atom_set) {
        auto a = nmr::Atom::Create(kv.second);
        a->pdb_atom_name = kv.first;
        a->residue_index = ri;
        atom_indices.push_back(p.AddAtom(std::move(a)));
    }
    p.MutableResidueAt(ri).atom_indices = atom_indices;

    return pp;
}

}  // namespace

// ============================================================================
// ParseAtomName cap-only flag population — codex Finding 3 single-parse
//
// `ComposeAtomSemantic` parses each atom name ONCE per atom and dispatches
// from typed flags. This unit test verifies that `ParseAtomName` populates
// the cap-only flags correctly so the dispatch path is sound.
// ============================================================================

TEST(LegacyAmberSemanticIntegration, ParseAtomNamePopulatesCapOnlyFlags) {
    namespace gen = nmr::topology_generated;

    // NTERM cap-only set: H1, H2, H3.
    for (const std::string& nm : {"H1", "H2", "H3"}) {
        const gen::ParsedAtomName p = gen::ParseAtomName(nm, "");
        EXPECT_TRUE(p.is_cap_only_n) << nm << " must set is_cap_only_n";
        EXPECT_FALSE(p.is_cap_only_c) << nm << " must NOT set is_cap_only_c";
        EXPECT_TRUE(p.is_n_terminus) << nm << " must set is_n_terminus";
        EXPECT_FALSE(p.is_c_terminus) << nm << " must NOT set is_c_terminus";
    }

    // CTERM cap-only set: OXT, HXT.
    for (const std::string& nm : {"OXT", "HXT"}) {
        const gen::ParsedAtomName p = gen::ParseAtomName(nm, "");
        EXPECT_TRUE(p.is_cap_only_c) << nm << " must set is_cap_only_c";
        EXPECT_FALSE(p.is_cap_only_n) << nm << " must NOT set is_cap_only_n";
        EXPECT_TRUE(p.is_c_terminus) << nm << " must set is_c_terminus";
        EXPECT_FALSE(p.is_n_terminus) << nm << " must NOT set is_n_terminus";
    }

    // Backbone amide H: NOT cap-only (chain table carries this).
    {
        const gen::ParsedAtomName p = gen::ParseAtomName("H", "");
        EXPECT_FALSE(p.is_cap_only_n);
        EXPECT_FALSE(p.is_cap_only_c);
        EXPECT_TRUE(p.is_backbone);
        EXPECT_EQ(BackboneRole::AmideHydrogen, p.backbone_role);
    }

    // HID-variant ring atom HD1: NOT cap-only.
    {
        const gen::ParsedAtomName p = gen::ParseAtomName("HD1", "ND1");
        EXPECT_FALSE(p.is_cap_only_n);
        EXPECT_FALSE(p.is_cap_only_c);
        EXPECT_FALSE(p.is_backbone);
    }

    // H2N / OT1 / OT2 are NOT recognised as cap-only (CharmmLegacy
    // cleanup, codex Finding 2). They parse as element-only sidechain
    // atoms — neither cap-only-n nor cap-only-c.
    for (const std::string& nm : {"H2N", "OT1", "OT2"}) {
        const gen::ParsedAtomName p = gen::ParseAtomName(nm, "");
        EXPECT_FALSE(p.is_cap_only_n)
            << nm << " must NOT set is_cap_only_n (CharmmLegacy retired)";
        EXPECT_FALSE(p.is_cap_only_c)
            << nm << " must NOT set is_cap_only_c (CharmmLegacy retired)";
        EXPECT_FALSE(p.is_n_terminus)
            << nm << " must NOT set is_n_terminus (CharmmLegacy retired)";
        EXPECT_FALSE(p.is_c_terminus)
            << nm << " must NOT set is_c_terminus (CharmmLegacy retired)";
    }
}

TEST(LegacyAmberSemanticIntegration, IsCapOnlyAtomNameDelegatesToParser) {
    namespace gen = nmr::topology_generated;

    // The string-taking wrapper agrees with the typed parser flags.
    EXPECT_TRUE(gen::IsCapOnlyAtomName("H1"));
    EXPECT_TRUE(gen::IsCapOnlyAtomName("H2"));
    EXPECT_TRUE(gen::IsCapOnlyAtomName("H3"));
    EXPECT_TRUE(gen::IsCapOnlyAtomName("OXT"));
    EXPECT_TRUE(gen::IsCapOnlyAtomName("HXT"));

    // Backbone, chain, and sidechain names: NOT cap-only.
    EXPECT_FALSE(gen::IsCapOnlyAtomName("N"));
    EXPECT_FALSE(gen::IsCapOnlyAtomName("CA"));
    EXPECT_FALSE(gen::IsCapOnlyAtomName("H"));
    EXPECT_FALSE(gen::IsCapOnlyAtomName("HA"));
    EXPECT_FALSE(gen::IsCapOnlyAtomName("CB"));
    EXPECT_FALSE(gen::IsCapOnlyAtomName("HD1"));

    // CharmmLegacy retired (codex Finding 2): H2N / OT1 / OT2 are NOT cap-only.
    EXPECT_FALSE(gen::IsCapOnlyAtomName("H2N"));
    EXPECT_FALSE(gen::IsCapOnlyAtomName("OT1"));
    EXPECT_FALSE(gen::IsCapOnlyAtomName("OT2"));
}


// ============================================================================
// AminoAcid::Unknown with named atoms — codex Finding 4 fail-loud
//
// `ComposeAtomSemantic` previously skipped Unknown residues silently;
// their atom slots in the result vector kept default-constructed rows
// (Element::Unknown etc.). HasAtomSemantic() then reported true for
// every atom, including the unknown-residue atoms, and calculators
// silently consumed default rows. Fix: abort if an Unknown residue
// carries any named atom. Toy-fixture stub-proteins (atoms without
// pdb_atom_name) keep working because the stub-fixture guard fires
// first.
// ============================================================================

namespace {

std::unique_ptr<nmr::Protein> MakeMixedKnownAndUnknownProtein() {
    auto pp = std::make_unique<nmr::Protein>();
    nmr::Protein& p = *pp;

    // Residue 1: ALA (known, with named atoms).
    nmr::Residue ala;
    ala.type = nmr::AminoAcid::ALA;
    ala.sequence_number = 1;
    ala.chain_id = "A";
    const size_t ri_ala = p.AddResidue(ala);
    {
        const std::pair<const char*, nmr::Element> ala_set[] = {
            {"N",   nmr::Element::N},
            {"CA",  nmr::Element::C},
            {"C",   nmr::Element::C},
            {"O",   nmr::Element::O},
            {"H",   nmr::Element::H},
            {"HA",  nmr::Element::H},
            {"CB",  nmr::Element::C},
            {"HB1", nmr::Element::H},
            {"HB2", nmr::Element::H},
            {"HB3", nmr::Element::H},
        };
        std::vector<size_t> atom_indices;
        for (const auto& kv : ala_set) {
            auto a = nmr::Atom::Create(kv.second);
            a->pdb_atom_name = kv.first;
            a->residue_index = ri_ala;
            atom_indices.push_back(p.AddAtom(std::move(a)));
        }
        p.MutableResidueAt(ri_ala).atom_indices = atom_indices;
    }

    // Residue 2: AminoAcid::Unknown but with named atoms.
    nmr::Residue unk;
    unk.type = nmr::AminoAcid::Unknown;
    unk.sequence_number = 2;
    unk.chain_id = "A";
    const size_t ri_unk = p.AddResidue(unk);
    {
        auto a = nmr::Atom::Create(nmr::Element::C);
        a->pdb_atom_name = "X1";  // arbitrary non-standard atom name
        a->residue_index = ri_unk;
        p.MutableResidueAt(ri_unk).atom_indices.push_back(
            p.AddAtom(std::move(a)));
    }

    return pp;
}

// Codex Finding F4 (2026-05-06): the prior stub-fixture guard skipped
// AminoAcid::Unknown residues when checking "does this protein have
// real atom names?" — so an all-Unknown-but-named-atoms protein
// silently treated as a stub fixture, returning an empty substrate
// without firing the fail-loud Unknown-residue loop. The redefined
// guard checks "any atom anywhere has a non-empty pdb_atom_name?";
// the all-Unknown case now flows through to the fail-loud loop.
std::unique_ptr<nmr::Protein> MakeAllUnknownNamedProtein() {
    auto pp = std::make_unique<nmr::Protein>();
    nmr::Protein& p = *pp;

    // Single residue, Unknown type, but with a named atom.
    nmr::Residue unk;
    unk.type = nmr::AminoAcid::Unknown;
    unk.sequence_number = 1;
    unk.chain_id = "A";
    const size_t ri = p.AddResidue(unk);
    {
        auto a = nmr::Atom::Create(nmr::Element::C);
        a->pdb_atom_name = "Q1";  // non-standard but non-empty
        a->residue_index = ri;
        p.MutableResidueAt(ri).atom_indices.push_back(
            p.AddAtom(std::move(a)));
    }

    return pp;
}

}  // namespace

using LegacyAmberSemanticIntegrationDeathTest = ::testing::Test;

TEST_F(LegacyAmberSemanticIntegrationDeathTest, UnknownResidueWithNamedAtomsAborts) {
    auto pp = MakeMixedKnownAndUnknownProtein();
    nmr::Protein& p = *pp;

    std::vector<nmr::Vec3> positions(p.AtomCount(), nmr::Vec3::Zero());
    EXPECT_DEATH(
        p.FinalizeConstruction(positions),
        "AminoAcid::Unknown residue.*carries .* named atom");
}

// Codex Finding F4 (2026-05-06): all-Unknown-but-named-atoms protein
// now hits the fail-loud loop instead of the stub-fixture short-circuit.
TEST_F(LegacyAmberSemanticIntegrationDeathTest, AllUnknownNamedAtomsAborts) {
    auto pp = MakeAllUnknownNamedProtein();
    nmr::Protein& p = *pp;

    std::vector<nmr::Vec3> positions(p.AtomCount(), nmr::Vec3::Zero());
    EXPECT_DEATH(
        p.FinalizeConstruction(positions),
        "AminoAcid::Unknown residue.*carries .* named atom");
}


TEST(LegacyAmberSemanticIntegration, ProNTermPreservesPyrrolidineRing) {
    auto pp = MakeNTermProProtein();
    nmr::Protein& p = *pp;

    // FinalizeConstruction with zero positions: substrate composition
    // runs from typed atom names; geometric bond detection does not
    // affect the chain-N row's ring_position.
    std::vector<nmr::Vec3> positions(p.AtomCount(), nmr::Vec3::Zero());
    p.FinalizeConstruction(positions);

    ASSERT_TRUE(p.HasTopology());
    const auto& topology = p.LegacyAmber();
    ASSERT_TRUE(topology.HasAtomSemantic())
        << "Substrate must populate for a chain-named PRO residue";

    // Find the chain N atom of PRO residue 1.
    const size_t n_atom = FindAtomIndex(p, /*residue_index=*/0, "N");
    ASSERT_NE(SIZE_MAX, n_atom) << "PRO residue must carry chain N";

    const AtomSemanticTable& sem_n = topology.SemanticAt(n_atom);

    // Locked: chain N's ring_position survives the backbone cap overlay.
    EXPECT_EQ(RingSystemKind::Pyrrolidine_Pro, sem_n.ring_position.primary.ring)
        << "PRO N at NTERM must retain Pyrrolidine_Pro after ApplyCapDelta. "
           "If this fails, ApplyCapDelta is overwriting ring_position from "
           "the cap entry's NotInRing placeholder — see "
           "spec/plan/topology-encoding-dependencies-2026-05-05.md §H.5 "
           "and codex-review Finding 1.";
    EXPECT_EQ(RingPositionLabel::ProRingNitrogen,
              sem_n.ring_position.primary.position)
        << "PRO N retains the Pro-specific ProRingNitrogen ring position "
           "label (chemistry-distinct from generic Heteroatom_NoH; surfaces "
           "the in-ring secondary-amine character of the Pro backbone N "
           "per Vega & Boyer 1979, Schubert et al. 2002).";

    // Locked: backbone_role survives (chain-controlled identity field).
    EXPECT_EQ(BackboneRole::Nitrogen, sem_n.backbone_role);

    // Locked: cap-controlled chemistry fields ARE applied. NTERM_CHARGED
    // is the AMBER ff14SB default; chain N becomes formal_charge = +1.
    // This verifies the rest of the cap overlay still fired (the bug was
    // in over-applying ring_position; the chemistry overlay must still
    // work).
    EXPECT_EQ(1, sem_n.formal_charge)
        << "PRO N at NTERM_CHARGED must have formal_charge = +1 after "
           "cap composition; cap-delta of formal_charge survives.";
}


// ============================================================================
// TRP indole 9-atom perimeter — Slice A added Indole_Trp_9 in the tertiary
// RingMembership slot for the 9 perimeter heavy atoms (CG, CD1, NE1, CE2,
// CZ2, CH2, CZ3, CE3, CD2) plus their attached H atoms (HD1, HE1, HZ2,
// HH2, HZ3, HE3). Encodes the conjugated π current circuit (Case 1995,
// J. Biomol. NMR 6, 341-346) as typed substrate alongside the chemical
// 5-ring + 6-ring decomposition.
//
// The TRP bridgeheads CE2 and CD2 are the canonical three-membership
// case: primary=Indole_Trp_5 + secondary=Indole_Trp_6 + tertiary=
// Indole_Trp_9. Non-bridgehead perimeter atoms have primary + tertiary
// populated with secondary EMPTY.
// ============================================================================
TEST(LegacyAmberSemanticIntegration, TrpPerimeterAtomsCarryIndoleTrp9Tertiary) {
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

    const std::vector<std::string> perimeter_heavy = {
        "CG", "CD1", "NE1", "CE2", "CZ2", "CH2", "CZ3", "CE3", "CD2"};
    const std::vector<std::string> bridgeheads = {"CE2", "CD2"};
    const std::vector<std::string> non_bridge_perimeter = {
        "CG", "CD1", "NE1", "CZ2", "CH2", "CZ3", "CE3"};

    int trp_seen = 0;
    for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        if (res.type != AminoAcid::TRP) continue;
        ++trp_seen;

        // Every perimeter heavy atom MUST carry the tertiary slot.
        for (const auto& nm : perimeter_heavy) {
            const size_t ai = FindAtomIndex(protein, ri, nm);
            ASSERT_NE(SIZE_MAX, ai) << "TRP must carry " << nm;
            const auto& sem = topology.SemanticAt(ai);
            EXPECT_TRUE(sem.ring_position.HasTertiaryRing())
                << "TRP perimeter atom " << nm << " must have tertiary slot populated";
            EXPECT_EQ(RingSystemKind::Indole_Trp_9,
                      sem.ring_position.tertiary.ring)
                << "TRP perimeter atom " << nm << " tertiary ring must be Indole_Trp_9";
            EXPECT_EQ(RingPositionLabel::PerimeterMember,
                      sem.ring_position.tertiary.position)
                << "TRP perimeter atom " << nm << " tertiary position must be PerimeterMember";
            EXPECT_EQ(9, static_cast<int>(sem.ring_position.tertiary.ring_size))
                << "TRP perimeter atom " << nm << " tertiary ring_size must be 9";
        }

        // Bridgeheads CE2 and CD2 occupy primary + secondary + tertiary.
        // MembershipCount() must be 3.
        for (const auto& nm : bridgeheads) {
            const size_t ai = FindAtomIndex(protein, ri, nm);
            ASSERT_NE(SIZE_MAX, ai);
            const auto& rp = topology.SemanticAt(ai).ring_position;
            EXPECT_EQ(3, rp.MembershipCount())
                << "TRP bridgehead " << nm << " must be in 3 rings (5+6+9)";
            EXPECT_EQ(RingSystemKind::Indole_Trp_5, rp.primary.ring);
            EXPECT_EQ(RingSystemKind::Indole_Trp_6, rp.secondary.ring);
            EXPECT_EQ(RingSystemKind::Indole_Trp_9, rp.tertiary.ring);
        }

        // Non-bridgehead perimeter atoms have primary + tertiary populated;
        // secondary is EMPTY. MembershipCount() must be 2 (the case the
        // old InTwoRings() predicate got wrong).
        for (const auto& nm : non_bridge_perimeter) {
            const size_t ai = FindAtomIndex(protein, ri, nm);
            ASSERT_NE(SIZE_MAX, ai);
            const auto& rp = topology.SemanticAt(ai).ring_position;
            EXPECT_EQ(2, rp.MembershipCount())
                << "TRP non-bridgehead perimeter atom " << nm
                << " must be in 2 rings (primary 5/6-ring + tertiary 9-ring); "
                   "this is the case where the old InTwoRings() returned the "
                   "wrong answer";
            EXPECT_TRUE(rp.HasPrimaryRing());
            EXPECT_FALSE(rp.HasSecondaryRing())
                << "TRP non-bridgehead " << nm
                << " has empty secondary slot; HasSecondaryRing() must be false";
            EXPECT_TRUE(rp.HasTertiaryRing());
        }
    }
    EXPECT_GT(trp_seen, 0) << "1P9J fixture must contain at least one TRP";
}

}  // namespace

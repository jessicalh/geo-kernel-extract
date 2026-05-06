#include <gtest/gtest.h>
#include "NamingRegistry.h"
#include "AminoAcidType.h"

using namespace nmr;

// ============================================================================
// NamingRegistry — residue-name translation tests
//
// The registry is now ONLY responsible for residue-name translation
// (HIS <-> HSD/HSE/HSP, CYS <-> CYX/CYM/CYS2 etc.). Atom-name
// canonicalisation lives on NamingApplicator below.
// ============================================================================

TEST(NamingRegistryTest, IsKnownForAllStandardAminoAcids) {
    auto& reg = GlobalNamingRegistry();
    const char* standard[] = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
        "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
        "THR", "TRP", "TYR", "VAL"
    };
    for (const char* name : standard) {
        EXPECT_TRUE(reg.IsKnownResidueName(name)) << name << " should be known";
    }
}

TEST(NamingRegistryTest, UnknownNamesRejected) {
    auto& reg = GlobalNamingRegistry();
    EXPECT_FALSE(reg.IsKnownResidueName("ZZZ"));
    EXPECT_FALSE(reg.IsKnownResidueName("FOO"));
    EXPECT_FALSE(reg.IsKnownResidueName(""));
}

TEST(NamingRegistryTest, HisToAmberProducesVariants) {
    auto& reg = GlobalNamingRegistry();
    EXPECT_EQ(reg.ResolveForTool("HIS", ToolContext::Amber, "delta"),   "HID");
    EXPECT_EQ(reg.ResolveForTool("HIS", ToolContext::Amber, "epsilon"), "HIE");
    EXPECT_EQ(reg.ResolveForTool("HIS", ToolContext::Amber, "doubly"),  "HIP");
}

TEST(NamingRegistryTest, HisToCharmmProducesVariants) {
    auto& reg = GlobalNamingRegistry();
    EXPECT_EQ(reg.ResolveForTool("HIS", ToolContext::Charmm, "delta"),   "HSD");
    EXPECT_EQ(reg.ResolveForTool("HIS", ToolContext::Charmm, "epsilon"), "HSE");
    EXPECT_EQ(reg.ResolveForTool("HIS", ToolContext::Charmm, "doubly"),  "HSP");
}

TEST(NamingRegistryTest, ToCanonicalMapsVariants) {
    auto& reg = GlobalNamingRegistry();
    EXPECT_EQ(reg.ToCanonical("HSD"), "HIS");
    EXPECT_EQ(reg.ToCanonical("HSE"), "HIS");
    EXPECT_EQ(reg.ToCanonical("HSP"), "HIS");
    EXPECT_EQ(reg.ToCanonical("HID"), "HIS");
    EXPECT_EQ(reg.ToCanonical("HIE"), "HIS");
    EXPECT_EQ(reg.ToCanonical("HIP"), "HIS");
    EXPECT_EQ(reg.ToCanonical("CYX"), "CYS");
    EXPECT_EQ(reg.ToCanonical("CYM"), "CYS");
    EXPECT_EQ(reg.ToCanonical("ASH"), "ASP");
    EXPECT_EQ(reg.ToCanonical("GLH"), "GLU");
    EXPECT_EQ(reg.ToCanonical("LYN"), "LYS");
    EXPECT_EQ(reg.ToCanonical("TYM"), "TYR");
    EXPECT_EQ(reg.ToCanonical("MSE"), "MET");
}

TEST(NamingRegistryTest, ToCanonicalUnknownReturnsEmpty) {
    auto& reg = GlobalNamingRegistry();
    EXPECT_EQ(reg.ToCanonical("ZZZ"), "");
    EXPECT_EQ(reg.ToCanonical("FOO"), "");
}

TEST(NamingRegistryTest, StandardNamesMappedToThemselves) {
    auto& reg = GlobalNamingRegistry();
    EXPECT_EQ(reg.ToCanonical("ALA"), "ALA");
    EXPECT_EQ(reg.ToCanonical("PHE"), "PHE");
    EXPECT_EQ(reg.ToCanonical("HIS"), "HIS");
}

TEST(NamingRegistryTest, CaseInsensitive) {
    auto& reg = GlobalNamingRegistry();
    EXPECT_TRUE(reg.IsKnownResidueName("ala"));
    EXPECT_TRUE(reg.IsKnownResidueName("Ala"));
    EXPECT_EQ(reg.ToCanonical("hsd"), "HIS");
}

TEST(NamingRegistryTest, NonTitratableReturnsCanonical) {
    auto& reg = GlobalNamingRegistry();
    // ALA has no variants -- ResolveForTool returns the canonical name
    EXPECT_EQ(reg.ResolveForTool("ALA", ToolContext::Amber), "ALA");
    EXPECT_EQ(reg.ResolveForTool("GLY", ToolContext::Charmm), "GLY");
}

TEST(NamingRegistryTest, KnownVariantNames) {
    auto& reg = GlobalNamingRegistry();
    EXPECT_TRUE(reg.IsKnownResidueName("HID"));
    EXPECT_TRUE(reg.IsKnownResidueName("HIE"));
    EXPECT_TRUE(reg.IsKnownResidueName("HIP"));
    EXPECT_TRUE(reg.IsKnownResidueName("ASH"));
    EXPECT_TRUE(reg.IsKnownResidueName("CYX"));
    EXPECT_TRUE(reg.IsKnownResidueName("MSE"));
}


// ============================================================================
// Helpers for NamingApplicator tests
// ============================================================================

namespace {

// Build a NamingContext with a single atom + a sibling-set snapshot.
NamingContext MakeContext(const std::string& input_name,
                          AminoAcid residue_type,
                          NamingSource source,
                          std::set<std::string> siblings,
                          int variant_index = -1,
                          TerminalState ts = TerminalState::Internal) {
    NamingContext ctx;
    ctx.source = source;
    ctx.input_name = input_name;
    ctx.residue_type = residue_type;
    ctx.variant_index = variant_index;
    ctx.terminal_state = ts;
    ctx.sibling_input_names = std::move(siblings);
    ctx.residue_sequence_number = 1;
    ctx.chain_id = "A";
    return ctx;
}

// Produce a sibling set from AminoAcidType::atoms — i.e. the canonical
// chain inventory for the residue type. Used by idempotency tests.
std::set<std::string> CanonicalSiblingSet(AminoAcid residue_type) {
    const AminoAcidType& aatype = GetAminoAcidType(residue_type);
    std::set<std::string> siblings;
    for (const auto& a : aatype.atoms) siblings.insert(a.name);
    return siblings;
}

}  // namespace


// ============================================================================
// NamingApplicator — Bundle B atom-name canonicalisation tests
//
// Verifies the rule-application object model introduced 2026-05-06
// (spec/plan/naming-applicator-architecture-sketch-2026-05-06.md).
// Each test exercises a specific predicate path: pass-through on
// canonical inputs, sibling-aware shifts on pre-Markley/pdb2gmx-RTP
// patterns, and the resolution-method branches.
// ============================================================================

TEST(NamingApplicatorTest, RulesAreLoaded) {
    const auto& app = GlobalNamingApplicator();
    // ~20 rules expected after the 2026-05-06 install across the
    // Amber/pdb2gmx-RTP + Markley1998 vocabularies. (The 3 historic
    // CharmmLegacy rules were removed 2026-05-06 with codex Finding 2
    // — no active load path tags inputs CharmmLegacy, fleet_amber
    // fixtures don't contain HN/OT1/OT2.) The exact count is brittle;
    // just assert a substantive load.
    EXPECT_GE(app.RuleCount(), 15u);
}


// ----------------------------------------------------------------------------
// LYN HZ shift — sibling-aware (HZ1+HZ2, no HZ3 ⇒ pre-Markley LYN)
// ----------------------------------------------------------------------------

TEST(NamingApplicatorLyn, PreMarkleyLynHz1ToHz2) {
    const auto& app = GlobalNamingApplicator();
    // Siblings {HZ1, HZ2, no HZ3} signal pre-Markley LYN: the residue
    // has LYN chemistry but non-canonical naming. HZ1 shifts to HZ2.
    auto ctx = MakeContext("HZ1", AminoAcid::LYS,
                           NamingSource::CifppPdbInput,
                           {"NZ", "HZ1", "HZ2"});
    EXPECT_EQ(app.Apply(ctx), "HZ2");
}

TEST(NamingApplicatorLyn, PreMarkleyLynHz2ToHz3) {
    const auto& app = GlobalNamingApplicator();
    auto ctx = MakeContext("HZ2", AminoAcid::LYS,
                           NamingSource::CifppPdbInput,
                           {"NZ", "HZ1", "HZ2"});
    EXPECT_EQ(app.Apply(ctx), "HZ3");
}

TEST(NamingApplicatorLyn, CanonicalLynPassThrough) {
    const auto& app = GlobalNamingApplicator();
    // Siblings {HZ2, HZ3, no HZ1} ⇒ canonical LYN: no shift.
    auto ctx_hz2 = MakeContext("HZ2", AminoAcid::LYS,
                               NamingSource::AmberFf14SBCanonical,
                               {"NZ", "HZ2", "HZ3"},
                               /*variant_index=*/0);
    EXPECT_EQ(app.Apply(ctx_hz2), "HZ2");

    auto ctx_hz3 = MakeContext("HZ3", AminoAcid::LYS,
                               NamingSource::AmberFf14SBCanonical,
                               {"NZ", "HZ2", "HZ3"},
                               /*variant_index=*/0);
    EXPECT_EQ(app.Apply(ctx_hz3), "HZ3");
}

TEST(NamingApplicatorLyn, CanonicalChargedLysPassThrough) {
    const auto& app = GlobalNamingApplicator();
    // Siblings {HZ1, HZ2, HZ3} ⇒ canonical charged LYS: all three Hs
    // pass through unchanged.
    for (const std::string& nm : {"HZ1", "HZ2", "HZ3"}) {
        auto ctx = MakeContext(nm, AminoAcid::LYS,
                               NamingSource::AmberFf14SBCanonical,
                               {"NZ", "HZ1", "HZ2", "HZ3"});
        EXPECT_EQ(app.Apply(ctx), nm)
            << "Canonical LYS HZ atom " << nm << " should pass through";
    }
}


// ----------------------------------------------------------------------------
// GLY HA collapse — Markley 1998 §2.1.2
// ----------------------------------------------------------------------------

TEST(NamingApplicatorGly, PreMarkleyHaToHa2) {
    const auto& app = GlobalNamingApplicator();
    // Siblings have collapsed "HA" only — pre-Markley fixture.
    auto ctx = MakeContext("HA", AminoAcid::GLY,
                           NamingSource::CifppPdbInput,
                           {"N", "H", "CA", "HA", "C", "O"});
    EXPECT_EQ(app.Apply(ctx), "HA2");
}

TEST(NamingApplicatorGly, CanonicalGlyHa2Ha3PassThrough) {
    const auto& app = GlobalNamingApplicator();
    // Canonical Gly siblings include both HA2 and HA3.
    auto ctx_ha2 = MakeContext("HA2", AminoAcid::GLY,
                               NamingSource::AmberFf14SBCanonical,
                               {"N", "H", "CA", "HA2", "HA3", "C", "O"});
    EXPECT_EQ(app.Apply(ctx_ha2), "HA2");

    auto ctx_ha3 = MakeContext("HA3", AminoAcid::GLY,
                               NamingSource::AmberFf14SBCanonical,
                               {"N", "H", "CA", "HA2", "HA3", "C", "O"});
    EXPECT_EQ(app.Apply(ctx_ha3), "HA3");
}


// ----------------------------------------------------------------------------
// pdb2gmx-AMBER-RTP shifts — sibling-aware
// ----------------------------------------------------------------------------

TEST(NamingApplicatorPdb2gmxRtp, ProHd1Hd2ShiftsToHd2Hd3) {
    const auto& app = GlobalNamingApplicator();
    // PRO with siblings {HD1, HD2, no HD3} — pdb2gmx-AMBER-RTP shape.
    std::set<std::string> siblings = {"N", "CA", "HA", "C", "O",
                                       "CB", "HB2", "HB3",
                                       "CG", "HG2", "HG3",
                                       "CD", "HD1", "HD2"};
    auto ctx_hd1 = MakeContext("HD1", AminoAcid::PRO,
                               NamingSource::Pdb2gmxAmberRtpDeviation,
                               siblings);
    EXPECT_EQ(app.Apply(ctx_hd1), "HD2");

    auto ctx_hd2 = MakeContext("HD2", AminoAcid::PRO,
                               NamingSource::Pdb2gmxAmberRtpDeviation,
                               siblings);
    EXPECT_EQ(app.Apply(ctx_hd2), "HD3");
}

TEST(NamingApplicatorPdb2gmxRtp, ProCanonicalHd2Hd3PassThrough) {
    const auto& app = GlobalNamingApplicator();
    // Canonical PRO siblings include HD2+HD3 (no HD1) — predicate
    // should not fire.
    std::set<std::string> siblings = {"N", "CA", "HA", "C", "O",
                                       "CB", "HB2", "HB3",
                                       "CG", "HG2", "HG3",
                                       "CD", "HD2", "HD3"};
    auto ctx = MakeContext("HD2", AminoAcid::PRO,
                           NamingSource::CifppPdbInput,
                           siblings);
    EXPECT_EQ(app.Apply(ctx), "HD2");
}

TEST(NamingApplicatorPdb2gmxRtp, IleHdMethylShifts) {
    const auto& app = GlobalNamingApplicator();
    // ILE with siblings {HD1, HD2, HD3} but NOT HD11/HD12/HD13.
    std::set<std::string> siblings = {"N", "CA", "HA", "C", "O",
                                       "CB", "HB",
                                       "CG1", "HG12", "HG13",
                                       "CG2", "HG21", "HG22", "HG23",
                                       "CD", "HD1", "HD2", "HD3"};
    EXPECT_EQ(app.Apply(MakeContext("HD1", AminoAcid::ILE,
                                     NamingSource::Pdb2gmxAmberRtpDeviation,
                                     siblings)), "HD11");
    EXPECT_EQ(app.Apply(MakeContext("HD2", AminoAcid::ILE,
                                     NamingSource::Pdb2gmxAmberRtpDeviation,
                                     siblings)), "HD12");
    EXPECT_EQ(app.Apply(MakeContext("HD3", AminoAcid::ILE,
                                     NamingSource::Pdb2gmxAmberRtpDeviation,
                                     siblings)), "HD13");
    EXPECT_EQ(app.Apply(MakeContext("CD", AminoAcid::ILE,
                                     NamingSource::Pdb2gmxAmberRtpDeviation,
                                     siblings)), "CD1");
}

TEST(NamingApplicatorPdb2gmxRtp, IleHg11Hg12Shift) {
    const auto& app = GlobalNamingApplicator();
    // ILE γ1-methylene shape: HG11+HG12 instead of canonical HG12+HG13.
    std::set<std::string> siblings = {"N", "CA", "HA", "C", "O",
                                       "CB", "HB",
                                       "CG1", "HG11", "HG12",
                                       "CG2", "HG21", "HG22", "HG23",
                                       "CD1", "HD11", "HD12", "HD13"};
    EXPECT_EQ(app.Apply(MakeContext("HG11", AminoAcid::ILE,
                                     NamingSource::Pdb2gmxAmberRtpDeviation,
                                     siblings)), "HG12");
    EXPECT_EQ(app.Apply(MakeContext("HG12", AminoAcid::ILE,
                                     NamingSource::Pdb2gmxAmberRtpDeviation,
                                     siblings)), "HG13");
}

TEST(NamingApplicatorPdb2gmxRtp, LysDeltaEpsilonMethyleneShifts) {
    const auto& app = GlobalNamingApplicator();
    // LYS pdb2gmx siblings: HD1+HD2 (no HD3) and HE1+HE2 (no HE3).
    std::set<std::string> siblings = {"N", "CA", "HA", "C", "O",
                                       "CB", "HB2", "HB3",
                                       "CG", "HG2", "HG3",
                                       "CD", "HD1", "HD2",
                                       "CE", "HE1", "HE2",
                                       "NZ", "HZ1", "HZ2", "HZ3"};
    EXPECT_EQ(app.Apply(MakeContext("HD1", AminoAcid::LYS,
                                     NamingSource::Pdb2gmxAmberRtpDeviation,
                                     siblings)), "HD2");
    EXPECT_EQ(app.Apply(MakeContext("HD2", AminoAcid::LYS,
                                     NamingSource::Pdb2gmxAmberRtpDeviation,
                                     siblings)), "HD3");
    EXPECT_EQ(app.Apply(MakeContext("HE1", AminoAcid::LYS,
                                     NamingSource::Pdb2gmxAmberRtpDeviation,
                                     siblings)), "HE2");
    EXPECT_EQ(app.Apply(MakeContext("HE2", AminoAcid::LYS,
                                     NamingSource::Pdb2gmxAmberRtpDeviation,
                                     siblings)), "HE3");
}

TEST(NamingApplicatorPdb2gmxRtp, ArgDeltaMethyleneShift) {
    const auto& app = GlobalNamingApplicator();
    std::set<std::string> siblings = {"N", "CA", "HA", "C", "O",
                                       "CB", "HB2", "HB3",
                                       "CG", "HG2", "HG3",
                                       "CD", "HD1", "HD2",
                                       "NE", "HE",
                                       "CZ", "NH1", "HH11", "HH12",
                                       "NH2", "HH21", "HH22"};
    EXPECT_EQ(app.Apply(MakeContext("HD1", AminoAcid::ARG,
                                     NamingSource::Pdb2gmxAmberRtpDeviation,
                                     siblings)), "HD2");
    EXPECT_EQ(app.Apply(MakeContext("HD2", AminoAcid::ARG,
                                     NamingSource::Pdb2gmxAmberRtpDeviation,
                                     siblings)), "HD3");
}


// ----------------------------------------------------------------------------
// Source-aware correctness — same input, different source ⇒ different fate
// ----------------------------------------------------------------------------

TEST(NamingApplicatorSourceAware, IleHd1UnderDifferentSources) {
    const auto& app = GlobalNamingApplicator();
    std::set<std::string> pdb2gmx_siblings = {"CB", "HB",
                                                "CG1", "HG12", "HG13",
                                                "CG2", "HG21", "HG22", "HG23",
                                                "CD", "HD1", "HD2", "HD3"};
    std::set<std::string> canonical_siblings = {"CB", "HB",
                                                  "CG1", "HG12", "HG13",
                                                  "CG2", "HG21", "HG22", "HG23",
                                                  "CD1", "HD11", "HD12", "HD13"};

    // Under Pdb2gmxAmberRtpDeviation, HD1 with pdb2gmx siblings shifts.
    auto ctx_pdb2gmx = MakeContext("HD1", AminoAcid::ILE,
                                    NamingSource::Pdb2gmxAmberRtpDeviation,
                                    pdb2gmx_siblings);
    EXPECT_EQ(app.Apply(ctx_pdb2gmx), "HD11");

    // Under CifppPdbInput with canonical siblings, HD11 (canonical) is
    // pass-through. The shift rule's predicate requires HD1+HD2+HD3 in
    // siblings — canonical siblings break that condition.
    auto ctx_cifpp = MakeContext("HD11", AminoAcid::ILE,
                                  NamingSource::CifppPdbInput,
                                  canonical_siblings);
    EXPECT_EQ(app.Apply(ctx_cifpp), "HD11");
}


// ----------------------------------------------------------------------------
// Idempotency — every canonical-state input passes through unchanged
// ----------------------------------------------------------------------------

TEST(NamingApplicatorIdempotency, CanonicalChainAtomsPassThrough) {
    const auto& app = GlobalNamingApplicator();
    // Iterate all 20 standard residues; for each, test every canonical
    // chain atom (from AminoAcidType::atoms) under each NamingSource.
    const NamingSource sources[] = {
        NamingSource::AmberFf14SBCanonical,
        NamingSource::Pdb2gmxAmberRtpDeviation,
        NamingSource::CifppPdbInput,
        NamingSource::OrcaEcho,
    };

    for (const AminoAcidType& aatype : AllAminoAcidTypes()) {
        if (aatype.index == AminoAcid::Unknown) continue;
        std::set<std::string> siblings = CanonicalSiblingSet(aatype.index);
        for (const auto& a : aatype.atoms) {
            for (NamingSource src : sources) {
                auto ctx = MakeContext(a.name, aatype.index, src, siblings);
                const std::string out = app.Apply(ctx);
                EXPECT_EQ(out, a.name)
                    << "canonical " << aatype.three_letter_code << "/"
                    << a.name << " under source " << NamingSourceName(src)
                    << " should pass through unchanged";
            }
        }
    }
}

TEST(NamingApplicatorIdempotency, LoadTimeToleranceAcceptsVariantExtensionAtoms) {
    const auto& app = GlobalNamingApplicator();
    // SCOPE: load-time tolerance window — variant_index == -1
    // (unresolved). Loaders cannot know the resolved tautomer until
    // ResolveProtonationStates fires, so the oracle accepts the union
    // of variant-extension atoms permissively. The variant-strict
    // checks (codex Finding 5) are exercised in the per-variant tests
    // below (HieRejectsHd1, HipAcceptsHd1AndHe2, AshAcceptsHd2,
    // GlhAcceptsHe2).
    struct VariantCase {
        AminoAcid res;
        const char* atom;
        std::set<std::string> siblings;
        const char* note;
    };
    const VariantCase cases[] = {
        // HID-form HIS: chain + HD1 present, no HE2.
        {AminoAcid::HIS, "HD1",
         {"N", "H", "CA", "HA", "CB", "HB2", "HB3", "CG",
          "ND1", "HD1", "CD2", "HD2", "CE1", "HE1", "NE2"},
         "HID HD1"},
        // HIP-form HIS: chain + HD1 + HE2.
        {AminoAcid::HIS, "HD1",
         {"N", "H", "CA", "HA", "CB", "HB2", "HB3", "CG",
          "ND1", "HD1", "CD2", "HD2", "CE1", "HE1", "NE2", "HE2"},
         "HIP HD1"},
        // ASH-form ASP: chain + HD2.
        {AminoAcid::ASP, "HD2",
         {"N", "H", "CA", "HA", "CB", "HB2", "HB3",
          "CG", "OD1", "OD2", "HD2"},
         "ASH HD2"},
        // GLH-form GLU: chain + HE2.
        {AminoAcid::GLU, "HE2",
         {"N", "H", "CA", "HA", "CB", "HB2", "HB3",
          "CG", "HG2", "HG3", "CD", "OE1", "OE2", "HE2"},
         "GLH HE2"},
    };

    for (const VariantCase& c : cases) {
        // variant_index defaults to -1 — load-time tolerance window.
        auto ctx = MakeContext(c.atom, c.res,
                               NamingSource::AmberFf14SBCanonical,
                               c.siblings);
        EXPECT_EQ(app.Apply(ctx), c.atom)
            << "Variant-extension atom " << c.note
            << " should pass through under load-time tolerance "
               "(variant_index = -1)";
    }
}


// ----------------------------------------------------------------------------
// Variant-aware canonicality oracle — codex Finding 5
//
// IsCanonical is variant-aware AFTER protonation resolution. Under
// variant_index >= 0 (resolved), only the variant's true atoms pass;
// HD1 on HIE-resolved HIS is NOT canonical (HD1 belongs to HID/HIP).
// ----------------------------------------------------------------------------

TEST(NamingApplicatorVariantAware, HipAcceptsHd1AndHe2) {
    const auto& app = GlobalNamingApplicator();
    // HIP (variant_index = 2): both ND1 and NE2 are protonated. HD1 and
    // HE2 are both canonical.
    std::set<std::string> siblings = {"N", "H", "CA", "HA",
                                       "CB", "HB2", "HB3", "CG",
                                       "ND1", "HD1", "CD2", "HD2",
                                       "CE1", "HE1", "NE2", "HE2"};
    {
        auto ctx_hd1 = MakeContext("HD1", AminoAcid::HIS,
                                    NamingSource::AmberFf14SBCanonical,
                                    siblings, /*variant_index=*/2);
        EXPECT_EQ(app.Apply(ctx_hd1), "HD1");
    }
    {
        auto ctx_he2 = MakeContext("HE2", AminoAcid::HIS,
                                    NamingSource::AmberFf14SBCanonical,
                                    siblings, /*variant_index=*/2);
        EXPECT_EQ(app.Apply(ctx_he2), "HE2");
    }
}

TEST(NamingApplicatorVariantAware, AshAcceptsHd2) {
    const auto& app = GlobalNamingApplicator();
    // ASH (variant_index = 0 for ASP): OD2 protonated. HD2 is canonical.
    std::set<std::string> siblings = {"N", "H", "CA", "HA",
                                       "CB", "HB2", "HB3",
                                       "CG", "OD1", "OD2", "HD2"};
    auto ctx = MakeContext("HD2", AminoAcid::ASP,
                           NamingSource::AmberFf14SBCanonical,
                           siblings, /*variant_index=*/0);
    EXPECT_EQ(app.Apply(ctx), "HD2");
}

TEST(NamingApplicatorVariantAware, GlhAcceptsHe2) {
    const auto& app = GlobalNamingApplicator();
    // GLH (variant_index = 0 for GLU): OE2 protonated. HE2 is canonical.
    std::set<std::string> siblings = {"N", "H", "CA", "HA",
                                       "CB", "HB2", "HB3",
                                       "CG", "HG2", "HG3",
                                       "CD", "OE1", "OE2", "HE2"};
    auto ctx = MakeContext("HE2", AminoAcid::GLU,
                           NamingSource::AmberFf14SBCanonical,
                           siblings, /*variant_index=*/0);
    EXPECT_EQ(app.Apply(ctx), "HE2");
}

TEST(NamingApplicatorVariantAwareDeathTest, HieRejectsHd1) {
    const auto& app = GlobalNamingApplicator();
    // HIE (variant_index = 1): NE2 protonated, ND1 not. HD1 is NOT
    // canonical for HIE; no rule fires; the canonicality oracle returns
    // false; FailUnresolved aborts.
    std::set<std::string> siblings = {"N", "H", "CA", "HA",
                                       "CB", "HB2", "HB3", "CG",
                                       "ND1", "CD2", "HD2",
                                       "CE1", "HE1", "NE2", "HE2",
                                       "HD1"};  // illegal extra HD1 on HIE
    auto ctx = MakeContext("HD1", AminoAcid::HIS,
                           NamingSource::AmberFf14SBCanonical,
                           siblings, /*variant_index=*/1);
    EXPECT_DEATH(app.Apply(ctx),
                 "no rule applies and input is not canonical");
}


// ----------------------------------------------------------------------------
// Canonicality oracle — direct exercise
// ----------------------------------------------------------------------------

TEST(NamingApplicatorOracle, RecognisesChainAtoms) {
    const auto& app = GlobalNamingApplicator();
    auto ctx = MakeContext("CA", AminoAcid::ALA,
                           NamingSource::AmberFf14SBCanonical,
                           {});
    EXPECT_TRUE(app.IsCanonical(ctx));

    ctx = MakeContext("CB", AminoAcid::PHE,
                      NamingSource::AmberFf14SBCanonical, {});
    EXPECT_TRUE(app.IsCanonical(ctx));
}

TEST(NamingApplicatorOracle, RecognisesCapAtoms) {
    const auto& app = GlobalNamingApplicator();
    // H2N is intentionally NOT in this set: removed 2026-05-06 with
    // CharmmLegacy cleanup (codex Finding 2). Any reappearance of H2N
    // in a future load path needs a concrete source tag and rule.
    for (const std::string& nm : {"H1", "H2", "H3", "OXT", "HXT"}) {
        auto ctx = MakeContext(nm, AminoAcid::ALA,
                               NamingSource::AmberFf14SBCanonical, {});
        EXPECT_TRUE(app.IsCanonical(ctx))
            << "Cap atom name '" << nm << "' should be canonical";
    }
}

TEST(NamingApplicatorOracle, RejectsUnknownAtomName) {
    const auto& app = GlobalNamingApplicator();
    auto ctx = MakeContext("ZZZZ", AminoAcid::ALA,
                           NamingSource::AmberFf14SBCanonical, {});
    EXPECT_FALSE(app.IsCanonical(ctx));
}


// ----------------------------------------------------------------------------
// Resolve-method branches — explicit per-case correctness
// ----------------------------------------------------------------------------

TEST(NamingApplicatorResolve, EmptyMapCanonicalInputPassesThrough) {
    const auto& app = GlobalNamingApplicator();
    // Atom CA on ALA: no rule fires (no shift rule matches CA on ALA).
    // Map empty; IsCanonical returns true; resolver returns input.
    auto ctx = MakeContext("CA", AminoAcid::ALA,
                           NamingSource::AmberFf14SBCanonical,
                           {"N", "CA", "C", "O", "H", "HA",
                            "CB", "HB1", "HB2", "HB3"});
    EXPECT_EQ(app.Apply(ctx), "CA");
}

TEST(NamingApplicatorResolve, SingleRuleFiringReturnsItsOutput) {
    const auto& app = GlobalNamingApplicator();
    // GlyHaToHa2_PreMarkley is the only rule that fires for GLY/HA
    // input under CifppPdbInput when siblings have HA only (no HA2/HA3).
    // This exercises the "exactly one rule" branch in Resolve().
    auto ctx = MakeContext("HA", AminoAcid::GLY,
                           NamingSource::CifppPdbInput,
                           {"N", "H", "CA", "HA", "C", "O"});
    EXPECT_EQ(app.Apply(ctx), "HA2");
}

TEST(NamingApplicatorResolve, MultipleRulesAgreeingReturnsSharedOutput) {
    // Construct a NamingApplicator with two test-only rules from
    // different NamingSource tags that BOTH fire on the same context
    // and propose the SAME output. Asserts the resolver's branch 3a
    // path returns the shared output.
    //
    // The production rule set does not currently reach this branch
    // (see NamingRegistry.cpp::Resolve() codex Finding 6 comment); the
    // test-only `CustomRules` constructor is the canonical exercise.
    NamingApplicator::CustomRules custom;
    custom.rules.push_back(NamingRule{
        NamingSource::AmberFf14SBCanonical,
        "TestRuleA_AmberCanonicalAlaCa",
        "test fixture: AMBER ff14SB canonical pass-through for ALA/CA",
        [](const NamingContext& c) {
            return c.residue_type == AminoAcid::ALA && c.input_name == "CA";
        },
        [](const NamingContext& c) { return c.input_name; },
    });
    custom.rules.push_back(NamingRule{
        NamingSource::Markley1998,
        "TestRuleB_MarkleyAlaCa",
        "test fixture: Markley 1998 pass-through for ALA/CA — agrees",
        [](const NamingContext& c) {
            return c.residue_type == AminoAcid::ALA && c.input_name == "CA";
        },
        [](const NamingContext& c) { return c.input_name; },
    });
    const NamingApplicator test_app(std::move(custom));

    auto ctx = MakeContext("CA", AminoAcid::ALA,
                           NamingSource::AmberFf14SBCanonical,
                           {"N", "CA", "C", "O"});

    // Both rules are in the rule list.
    ASSERT_EQ(test_app.RuleCount(), 2u);

    // The resolver's branch 3a (all-agree) returns the shared output.
    EXPECT_EQ(test_app.Apply(ctx), "CA");
}


// ----------------------------------------------------------------------------
// Fail-on-unknown — death test for unresolvable input
// ----------------------------------------------------------------------------

using NamingApplicatorDeathTest = ::testing::Test;

TEST_F(NamingApplicatorDeathTest, UnknownAtomNameUnderUnknownSourceAborts) {
    const auto& app = GlobalNamingApplicator();
    // An atom name no rule matches AND not present in the canonicality
    // oracle. ZZZZ is fictional, residue ALA, no source rules match.
    auto ctx = MakeContext("ZZZZ", AminoAcid::ALA,
                           NamingSource::Unknown,
                           {"ZZZZ", "CA"});
    EXPECT_DEATH(app.Apply(ctx),
                 "no rule applies and input is not canonical");
}


// ----------------------------------------------------------------------------
// ApplyResidue — sibling snapshot semantics
// ----------------------------------------------------------------------------

TEST(NamingApplicatorApplyResidue, LynPreMarkleyShiftWholeResidue) {
    const auto& app = GlobalNamingApplicator();
    // Pre-Markley LYN siblings: HZ1+HZ2, no HZ3. Other side-chain atoms
    // in canonical AMBER ff14SB form (HD2/HD3, HE2/HE3) so the test
    // isolates the LYN HZ shift behaviour from any Pdb2gmx β/γ-methylene
    // shift behaviour. ApplyResidue snapshots siblings once; both HZ
    // atoms see the original {..., HZ1, HZ2} set.
    //
    // Source = CifppPdbInput: rules tagged Pdb2gmxAmberRtpDeviation do
    // not fire (their source-tag predicate returns false). The LYN HZ
    // shift rules are tagged AmberFf14SBCanonical (source-agnostic) so
    // they DO fire on this fixture, which is the test's purpose.
    //
    // Note on chemistry: AMBER ff14SB canonical LYS β/γ/δ/ε methylene
    // numbering starts at 2 (not 1) per Markley 1998 §2.1.2 — HD1 is
    // NOT canonical for LYS. The fixture uses HD2/HD3 + HE2/HE3 to
    // stay on canonical ground for everything except the LYN HZ pair.
    std::vector<std::string> input_names = {
        "N", "H", "CA", "HA", "CB", "HB2", "HB3",
        "CG", "HG2", "HG3", "CD", "HD2", "HD3",
        "CE", "HE2", "HE3", "NZ", "HZ1", "HZ2"};
    std::vector<std::string> parent_names(input_names.size());

    const auto outs = app.ApplyResidue(
        input_names, parent_names,
        AminoAcid::LYS,
        /*variant_index=*/-1,
        TerminalState::Internal,
        NamingSource::CifppPdbInput,
        /*sequence_number=*/28,
        "A");
    ASSERT_EQ(outs.size(), input_names.size());

    // LYN HZ shift fires (AmberFf14SBCanonical source-agnostic):
    // HZ1 -> HZ2, HZ2 -> HZ3.
    EXPECT_EQ(outs[input_names.size() - 2], "HZ2");
    EXPECT_EQ(outs[input_names.size() - 1], "HZ3");

    // All other (canonical) atoms pass through unchanged.
    auto find = [&](const std::string& n) -> std::string {
        for (size_t i = 0; i < input_names.size(); ++i)
            if (input_names[i] == n) return outs[i];
        return "<missing>";
    };
    EXPECT_EQ(find("HD2"), "HD2");
    EXPECT_EQ(find("HD3"), "HD3");
    EXPECT_EQ(find("HE2"), "HE2");
    EXPECT_EQ(find("HE3"), "HE3");
}

TEST(NamingApplicatorApplyResidue, IleSnapshotIsIndependentOfPerAtomOrder) {
    const auto& app = GlobalNamingApplicator();
    // ILE pdb2gmx-RTP shape. Snapshot at start of pass; per-atom
    // application uses snapshot. Outputs reorder relative to input
    // would NOT be visible in apply order: we confirm outputs match
    // the snapshot-driven expectations regardless of position.
    std::vector<std::string> input_names = {
        "N", "H", "CA", "HA", "C", "O",
        "CB", "HB",
        "CG1", "HG11", "HG12",
        "CG2", "HG21", "HG22", "HG23",
        "CD", "HD1", "HD2", "HD3"};
    std::vector<std::string> parent_names(input_names.size());

    const auto outs = app.ApplyResidue(
        input_names, parent_names,
        AminoAcid::ILE,
        /*variant_index=*/-1,
        TerminalState::Internal,
        NamingSource::Pdb2gmxAmberRtpDeviation,
        /*sequence_number=*/1,
        "A");
    ASSERT_EQ(outs.size(), input_names.size());

    auto find = [&](const std::string& n) -> std::string {
        for (size_t i = 0; i < input_names.size(); ++i)
            if (input_names[i] == n) return outs[i];
        return "<missing>";
    };
    EXPECT_EQ(find("CD"),  "CD1");
    EXPECT_EQ(find("HD1"), "HD11");
    EXPECT_EQ(find("HD2"), "HD12");
    EXPECT_EQ(find("HD3"), "HD13");
    EXPECT_EQ(find("HG11"), "HG12");
    EXPECT_EQ(find("HG12"), "HG13");
    // HG21..HG23 are canonical; pass-through.
    EXPECT_EQ(find("HG21"), "HG21");
    EXPECT_EQ(find("HG22"), "HG22");
    EXPECT_EQ(find("HG23"), "HG23");
}

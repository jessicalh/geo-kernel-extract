#include <gtest/gtest.h>
#include "NamingRegistry.h"

using namespace nmr;

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

TEST(NamingRegistryTest, TranslateHNFromCharmmToStandard) {
    auto& reg = GlobalNamingRegistry();
    EXPECT_EQ(reg.TranslateAtomName("HN", "ALA", ToolContext::Charmm, ToolContext::Standard), "H");
}

TEST(NamingRegistryTest, TranslateHToCharmmGivesHN) {
    auto& reg = GlobalNamingRegistry();
    EXPECT_EQ(reg.TranslateAtomName("H", "ALA", ToolContext::Standard, ToolContext::Charmm), "HN");
}

TEST(NamingRegistryTest, SameContextNoTranslation) {
    auto& reg = GlobalNamingRegistry();
    EXPECT_EQ(reg.TranslateAtomName("CA", "ALA",
        ToolContext::Standard, ToolContext::Standard), "CA");
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
// CanonicaliseAmberAtomName tests — Phase 1 substrate-runtime canonicalisation.
//
// Verifies the two-stage canonicalisation pipeline introduced for the
// chemistry-substrate runtime integration: Stage 1 collapses the wire-
// format universe (CHARMM-port HN -> H, OT1/OT2 -> O/OXT), Stage 2 fires
// residue-context-keyed rules (LYN HZ1/HZ2 -> HZ2/HZ3, GLY HA -> HA2)
// regardless of source context.
// ============================================================================

TEST(NamingRegistryCanonicaliseTest, StandardInputPassThrough) {
    auto& reg = GlobalNamingRegistry();
    // Standard / IUPAC inputs need no translation.
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("CA",  "ALA", ToolContext::Standard), "CA");
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("N",   "ALA", ToolContext::Standard), "N");
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("H",   "ALA", ToolContext::Standard), "H");
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("HA",  "ALA", ToolContext::Standard), "HA");
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("HG21","THR", ToolContext::Standard), "HG21");
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("HZ1", "LYS", ToolContext::Standard), "HZ1");
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("HZ2", "LYS", ToolContext::Standard), "HZ2");
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("HZ3", "LYS", ToolContext::Standard), "HZ3");
}

TEST(NamingRegistryCanonicaliseTest, CharmmHnCollapsesToH) {
    auto& reg = GlobalNamingRegistry();
    // Backbone amide H — CHARMM convention spells it HN; canonical AMBER is H.
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("HN", "ALA", ToolContext::Charmm), "H");
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("HN", "GLY", ToolContext::Charmm), "H");
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("HN", "PRO", ToolContext::Charmm), "H");
}

TEST(NamingRegistryCanonicaliseTest, CtermCharmmOxygensCollapseToAmber) {
    auto& reg = GlobalNamingRegistry();
    // CHARMM C-terminal carboxylate: OT1 (chain carbonyl O) and OT2
    // (carboxyl O). AMBER canonical is O / OXT.
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("OT1", "ALA", ToolContext::Charmm), "O");
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("OT2", "ALA", ToolContext::Charmm), "OXT");
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("OT1", "GLY", ToolContext::Charmm), "O");
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("OT2", "GLY", ToolContext::Charmm), "OXT");
}

TEST(NamingRegistryCanonicaliseTest, LynHzPreMarkleyToCanonical) {
    auto& reg = GlobalNamingRegistry();
    // LYN canonical is HZ2/HZ3 (NH2, two Hs, HZ1 absent at deprotonation).
    // Pre-Markley fixtures sometimes carry HZ1/HZ2 on residues recorded
    // as LYN. The rule fires in Stage 2 regardless of source context.
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("HZ1", "LYN", ToolContext::Standard), "HZ2");
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("HZ2", "LYN", ToolContext::Standard), "HZ3");
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("HZ1", "LYN", ToolContext::Amber), "HZ2");
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("HZ2", "LYN", ToolContext::Amber), "HZ3");
    // Critical chain-canonical: when the residue is LYS, the LYN rule
    // does NOT fire and HZ1/HZ2/HZ3 pass through (LYS NH3+ chain
    // chemistry retains all three H atoms).
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("HZ1", "LYS", ToolContext::Standard), "HZ1");
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("HZ2", "LYS", ToolContext::Standard), "HZ2");
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("HZ3", "LYS", ToolContext::Standard), "HZ3");
}

TEST(NamingRegistryCanonicaliseTest, GlyHaCollapsedMethyleneToCanonical) {
    auto& reg = GlobalNamingRegistry();
    // GLY canonical is HA2/HA3. Pre-Markley fixtures with a single
    // collapsed "HA" on Gly map to the IUPAC pro-R form HA2.
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("HA",  "GLY", ToolContext::Standard), "HA2");
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("HA1", "GLY", ToolContext::Standard), "HA3");
    // Non-Gly residues canonically use HA: rule does NOT fire.
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("HA",  "ALA", ToolContext::Standard), "HA");
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("HA",  "MET", ToolContext::Standard), "HA");
}

TEST(NamingRegistryCanonicaliseTest, NTermCapH1StaysH1NotBackboneH) {
    auto& reg = GlobalNamingRegistry();
    // The N-terminal cap atom H1 is a DISTINCT cap-only atom in
    // AMBER ff14SB (kCapNtermCharged: H1, H2, H3 on the +1 ammonium
    // nitrogen), NOT an alias of the backbone amide H. The earlier
    // OrcaRunLoader inline `H || HN || H1` cache populator confused
    // these; canonicalisation must NOT rewrite H1 -> H.
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("H1", "ALA", ToolContext::Charmm),  "H1");
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("H1", "ALA", ToolContext::Standard), "H1");
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("H2", "ALA", ToolContext::Standard), "H2");
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("H3", "ALA", ToolContext::Standard), "H3");
}

TEST(NamingRegistryCanonicaliseTest, LynRulesDoNotChainRewrite) {
    auto& reg = GlobalNamingRegistry();
    // The LYN HZ rules live only in the (Standard, Amber) keyspace,
    // so a CHARMM source flowing HZ1 through Stage 1 (Charmm,Standard)
    // and Stage 2 (Standard,Amber) should NOT chain-rewrite:
    // HZ1 -> HZ1 (Stage 1 no rule) -> HZ2 (Stage 2 LYN rule). Final HZ2.
    // Without this discipline, an additional (Charmm,Standard) rule
    // would chain HZ1 -> HZ2 -> HZ3 in one pass.
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("HZ1", "LYN", ToolContext::Charmm), "HZ2");
    EXPECT_EQ(reg.CanonicaliseAmberAtomName("HZ2", "LYN", ToolContext::Charmm), "HZ3");
}

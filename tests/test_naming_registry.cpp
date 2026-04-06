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

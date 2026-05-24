#include <gtest/gtest.h>
#include "Types.h"
#include "AminoAcidType.h"

using namespace nmr;

TEST(AminoAcid, ThreeLetterCodeRoundtrip) {
    for (int i = 0; i < 20; ++i) {
        AminoAcid aa = static_cast<AminoAcid>(i);
        std::string code = ThreeLetterCodeForAminoAcid(aa);
        AminoAcid back = AminoAcidFromThreeLetterCode(code);
        EXPECT_EQ(aa, back) << "Failed for " << code;
    }
}

TEST(AminoAcid, AromaticResiduePHE) {
    EXPECT_TRUE(IsAromaticAminoAcid(AminoAcid::PHE));
    EXPECT_TRUE(IsAromaticAminoAcid(AminoAcid::TYR));
    EXPECT_TRUE(IsAromaticAminoAcid(AminoAcid::TRP));
    EXPECT_TRUE(IsAromaticAminoAcid(AminoAcid::HIS));
}

TEST(AminoAcid, NonAromaticResidue) {
    EXPECT_FALSE(IsAromaticAminoAcid(AminoAcid::ALA));
    EXPECT_FALSE(IsAromaticAminoAcid(AminoAcid::GLY));
    EXPECT_FALSE(IsAromaticAminoAcid(AminoAcid::LEU));
}

TEST(AminoAcid, ProtonationVariants) {
    AminoAcid his = AminoAcidFromThreeLetterCode("HID");
    EXPECT_EQ(his, AminoAcid::HIS);

    AminoAcid ash = AminoAcidFromThreeLetterCode("ASH");
    EXPECT_EQ(ash, AminoAcid::ASP);
}

// Bundle C / Slice B (2026-05-07): the previous PHERingDefined and
// TRPThreeRings tests asserted on AminoAcidType::rings (the
// AminoAcidRing struct + per-residue rings[] arrays). Both were
// deleted along with the string-based DetectAromaticRings. The
// equivalent post-Bundle-C assertions ("PHE produces one
// PheBenzeneRing", "TRP produces TrpBenzene + TrpPyrrole +
// TrpPerimeter") run end-to-end on real fixtures (1UBQ, 1P9J) via
// existing tests in test_pdb_loading and test_foundation_results,
// since rings now require positions + substrate. The static
// presence-of-table assertions have no equivalent because the
// table itself has been retired.

TEST(AminoAcidType, PHEIsAromatic) {
    const auto& phe = GetAminoAcidType(AminoAcid::PHE);
    EXPECT_TRUE(phe.is_aromatic);
}

TEST(AminoAcidType, TRPIsAromatic) {
    const auto& trp = GetAminoAcidType(AminoAcid::TRP);
    EXPECT_TRUE(trp.is_aromatic);
}

TEST(AminoAcidType, HISVariants) {
    const auto& his = GetAminoAcidType(AminoAcid::HIS);
    EXPECT_TRUE(his.is_titratable);
    EXPECT_EQ(his.variants.size(), 3);
    EXPECT_EQ(std::string(his.variants[0].name), "HID");
}

TEST(AminoAcidType, ProNoAmideH) {
    const auto& pro = GetAminoAcidType(AminoAcid::PRO);
    EXPECT_FALSE(pro.has_amide_H);
}

TEST(AminoAcidType, GlyNoChiAngles) {
    const auto& gly = GetAminoAcidType(AminoAcid::GLY);
    EXPECT_EQ(gly.chi_angle_count, 0);
}

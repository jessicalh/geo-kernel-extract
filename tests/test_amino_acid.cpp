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

TEST(AminoAcidType, PHERingDefined) {
    const auto& phe = GetAminoAcidType(AminoAcid::PHE);
    EXPECT_TRUE(phe.is_aromatic);
    EXPECT_EQ(phe.rings.size(), 1);
    EXPECT_EQ(phe.rings[0].type_index, RingTypeIndex::PheBenzene);
    EXPECT_EQ(phe.rings[0].atom_names.size(), 6);
}

TEST(AminoAcidType, TRPThreeRings) {
    const auto& trp = GetAminoAcidType(AminoAcid::TRP);
    EXPECT_TRUE(trp.is_aromatic);
    EXPECT_EQ(trp.rings.size(), 3);
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

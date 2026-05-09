#include <gtest/gtest.h>

#include "generated/LegacyAmberSemanticTables.h"

namespace {

using nmr::AminoAcid;
using nmr::AtomMechanicalIdentity;
using nmr::AtomSemanticTable;
using nmr::BackboneRole;
using nmr::BranchAddress;
using nmr::DiastereotopicIndex;
using nmr::Element;
using nmr::Locant;
using nmr::PlanarGroupKind;
using nmr::PlanarStereo;
using nmr::PolarHKind;
using nmr::ProchiralStereo;
using nmr::PseudoatomKind;
using nmr::RingPositionLabel;
using nmr::RingSystemKind;
using nmr::TerminalState;

namespace gen = nmr::topology_generated;

AtomMechanicalIdentity Id(Element element,
                          Locant locant,
                          BranchAddress branch,
                          DiastereotopicIndex di_index,
                          BackboneRole backbone_role) {
    return AtomMechanicalIdentity{element, locant, branch, di_index,
                                  backbone_role};
}

AtomMechanicalIdentity IdentityOf(const AtomSemanticTable& row) {
    return AtomMechanicalIdentity{row.element, row.locant, row.branch,
                                  row.di_index, row.backbone_role};
}

bool SamePseudoatom(const nmr::PseudoatomMembership& a,
                    const nmr::PseudoatomMembership& b) {
    return a.kind == b.kind && a.locant == b.locant &&
           a.branch == b.branch && a.in_super_group == b.in_super_group;
}

bool SameRingMembership(const nmr::RingMembership& a,
                        const nmr::RingMembership& b) {
    return a.ring == b.ring && a.position == b.position &&
           a.ring_size == b.ring_size && a.aromatic == b.aromatic &&
           a.planar == b.planar && a.n_heteroatoms == b.n_heteroatoms;
}

bool SameRingPosition(const nmr::RingPosition& a,
                      const nmr::RingPosition& b) {
    return SameRingMembership(a.primary, b.primary) &&
           SameRingMembership(a.secondary, b.secondary);
}

TEST(TopologySemanticApi, LookupByFailsFastAtResidueAndVariantLevel) {
    const auto ala_n = Id(Element::N, Locant::None, {0, 0},
                          DiastereotopicIndex::None,
                          BackboneRole::Nitrogen);

    ASSERT_NE(nullptr, gen::LookupBy(AminoAcid::ALA, gen::kBaseVariantIdx, ala_n));
    EXPECT_EQ(nullptr, gen::LookupBy(AminoAcid::Unknown,
                                     gen::kBaseVariantIdx, ala_n));
    EXPECT_EQ(nullptr, gen::LookupBy(AminoAcid::ALA, 0, ala_n))
        << "ALA has no variant 0; chain callers must use kBaseVariantIdx";
    EXPECT_EQ(nullptr, gen::LookupBy(AminoAcid::ALA, 254, ala_n));

    AtomMechanicalIdentity impossible = ala_n;
    impossible.element = Element::Unknown;
    EXPECT_EQ(nullptr, gen::LookupBy(AminoAcid::ALA,
                                     gen::kBaseVariantIdx, impossible));
}

TEST(TopologySemanticApi, LookupByVariantIndicesSelectVariantChemistry) {
    const auto asp_od2 = Id(Element::O, Locant::Delta, {2, 0},
                            DiastereotopicIndex::None, BackboneRole::None);
    const auto ash_hd2 = Id(Element::H, Locant::Delta, {2, 0},
                            DiastereotopicIndex::None, BackboneRole::None);
    const auto glu_oe2 = Id(Element::O, Locant::Epsilon, {2, 0},
                            DiastereotopicIndex::None, BackboneRole::None);
    const auto glh_he2 = Id(Element::H, Locant::Epsilon, {2, 0},
                            DiastereotopicIndex::None, BackboneRole::None);

    const auto* asp_base = gen::LookupBy(AminoAcid::ASP,
                                         gen::kBaseVariantIdx, asp_od2);
    const auto* ash = gen::LookupBy(AminoAcid::ASP, 0, asp_od2);
    ASSERT_NE(nullptr, asp_base);
    ASSERT_NE(nullptr, ash);
    EXPECT_EQ(-1, asp_base->formal_charge);
    EXPECT_EQ(0, ash->formal_charge);
    EXPECT_EQ(nullptr, gen::LookupBy(AminoAcid::ASP,
                                     gen::kBaseVariantIdx, ash_hd2));
    EXPECT_NE(nullptr, gen::LookupBy(AminoAcid::ASP, 0, ash_hd2));

    const auto* glu_base = gen::LookupBy(AminoAcid::GLU,
                                         gen::kBaseVariantIdx, glu_oe2);
    const auto* glh = gen::LookupBy(AminoAcid::GLU, 0, glu_oe2);
    ASSERT_NE(nullptr, glu_base);
    ASSERT_NE(nullptr, glh);
    EXPECT_EQ(-1, glu_base->formal_charge);
    EXPECT_EQ(0, glh->formal_charge);
    EXPECT_EQ(nullptr, gen::LookupBy(AminoAcid::GLU,
                                     gen::kBaseVariantIdx, glh_he2));
    EXPECT_NE(nullptr, gen::LookupBy(AminoAcid::GLU, 0, glh_he2));
}

TEST(TopologySemanticApi, LookupByHistidineVariantsKeepTautomerAndChargeSemantics) {
    const auto nd1 = Id(Element::N, Locant::Delta, {1, 0},
                        DiastereotopicIndex::None, BackboneRole::None);
    const auto ne2 = Id(Element::N, Locant::Epsilon, {2, 0},
                        DiastereotopicIndex::None, BackboneRole::None);
    const auto hd1 = Id(Element::H, Locant::Delta, {1, 0},
                        DiastereotopicIndex::None, BackboneRole::None);
    const auto he2 = Id(Element::H, Locant::Epsilon, {2, 0},
                        DiastereotopicIndex::None, BackboneRole::None);

    ASSERT_NE(nullptr, gen::LookupBy(AminoAcid::HIS, 0, hd1));
    EXPECT_EQ(nullptr, gen::LookupBy(AminoAcid::HIS, 0, he2));
    EXPECT_EQ(nullptr, gen::LookupBy(AminoAcid::HIS, 1, hd1));
    ASSERT_NE(nullptr, gen::LookupBy(AminoAcid::HIS, 1, he2));
    ASSERT_NE(nullptr, gen::LookupBy(AminoAcid::HIS, 2, hd1));
    ASSERT_NE(nullptr, gen::LookupBy(AminoAcid::HIS, 2, he2));

    const auto* hid_nd1 = gen::LookupBy(AminoAcid::HIS, 0, nd1);
    const auto* hid_ne2 = gen::LookupBy(AminoAcid::HIS, 0, ne2);
    const auto* hip_nd1 = gen::LookupBy(AminoAcid::HIS, 2, nd1);
    const auto* hip_ne2 = gen::LookupBy(AminoAcid::HIS, 2, ne2);
    ASSERT_NE(nullptr, hid_nd1);
    ASSERT_NE(nullptr, hid_ne2);
    ASSERT_NE(nullptr, hip_nd1);
    ASSERT_NE(nullptr, hip_ne2);
    EXPECT_EQ(0, hid_nd1->formal_charge);
    EXPECT_EQ(0, hid_ne2->formal_charge);
    EXPECT_EQ(0, hip_nd1->formal_charge);
    EXPECT_EQ(1, hip_ne2->formal_charge);
}

TEST(TopologySemanticApi, LookupCapCoversTerminalStatesAndMissesCleanly) {
    const auto nterm_n = Id(Element::N, Locant::None, {0, 0},
                            DiastereotopicIndex::None, BackboneRole::Nitrogen);
    const auto nterm_h = Id(Element::H, Locant::None, {0, 0},
                            DiastereotopicIndex::None, BackboneRole::None);
    const auto cterm_c = Id(Element::C, Locant::None, {0, 0},
                            DiastereotopicIndex::None,
                            BackboneRole::CarbonylCarbon);
    const auto cterm_o = Id(Element::O, Locant::None, {0, 0},
                            DiastereotopicIndex::None,
                            BackboneRole::CarbonylOxygen);
    const auto cterm_oxt = Id(Element::O, Locant::None, {0, 0},
                              DiastereotopicIndex::None, BackboneRole::None);
    const auto cterm_hxt = Id(Element::H, Locant::None, {0, 0},
                              DiastereotopicIndex::None, BackboneRole::None);

    ASSERT_NE(nullptr, gen::LookupCap(TerminalState::NtermCharged, nterm_n));
    ASSERT_NE(nullptr, gen::LookupCap(TerminalState::NtermCharged, nterm_h));
    EXPECT_EQ(PolarHKind::AmmoniumNH,
              gen::LookupCap(TerminalState::NtermCharged, nterm_h)->polar_h);
    ASSERT_NE(nullptr, gen::LookupCap(TerminalState::NtermNeutral, nterm_h));
    EXPECT_EQ(PolarHKind::AmineNH,
              gen::LookupCap(TerminalState::NtermNeutral, nterm_h)->polar_h);

    ASSERT_NE(nullptr, gen::LookupCap(TerminalState::CtermDeprotonated, cterm_c));
    ASSERT_NE(nullptr, gen::LookupCap(TerminalState::CtermDeprotonated, cterm_o));
    ASSERT_NE(nullptr, gen::LookupCap(TerminalState::CtermDeprotonated, cterm_oxt));
    EXPECT_EQ(-1, gen::LookupCap(TerminalState::CtermDeprotonated,
                                 cterm_oxt)->formal_charge);

    ASSERT_NE(nullptr, gen::LookupCap(TerminalState::CtermProtonated, cterm_oxt));
    ASSERT_NE(nullptr, gen::LookupCap(TerminalState::CtermProtonated, cterm_hxt));
    EXPECT_EQ(0, gen::LookupCap(TerminalState::CtermProtonated,
                                cterm_oxt)->formal_charge);
    EXPECT_EQ(PolarHKind::CarboxylOH,
              gen::LookupCap(TerminalState::CtermProtonated, cterm_hxt)->polar_h);

    EXPECT_EQ(nullptr, gen::LookupCap(TerminalState::Internal, nterm_n));
    EXPECT_EQ(nullptr, gen::LookupCap(TerminalState::NtermNeutral, cterm_oxt));
    EXPECT_EQ(nullptr, gen::LookupCap(TerminalState::CtermDeprotonated, nterm_h));
}

TEST(TopologySemanticApi, ApplyCapDeltaChangesOnlyDeclaredDeltaFields) {
    AtomSemanticTable chain;
    chain.element = Element::N;
    chain.locant = Locant::Gamma;
    chain.branch = {1, 2};
    chain.di_index = DiastereotopicIndex::Position3;
    chain.backbone_role = BackboneRole::Nitrogen;
    chain.prochiral = ProchiralStereo::ProR;
    chain.planar_group = PlanarGroupKind::PeptideAmide;
    chain.planar_stereo = PlanarStereo::E;
    chain.pseudoatom = {PseudoatomKind::M, 3, 1, true};
    chain.polar_h = PolarHKind::NotPolar;
    chain.ring_position.primary = {RingSystemKind::Pyrrolidine_Pro,
                                   RingPositionLabel::Saturated, 5,
                                   false, false, 1};
    chain.ring_position.secondary = {RingSystemKind::NotInRing,
                                     RingPositionLabel::NotInRing, 0,
                                     false, false, 0};
    chain.aromatic = true;
    chain.formal_charge = -1;
    chain.is_exchangeable = false;
    chain.equivalence_class = 42;

    AtomSemanticTable cap = chain;
    cap.element = Element::O;
    cap.locant = Locant::Zeta;
    cap.branch = {2, 1};
    cap.di_index = DiastereotopicIndex::Position2;
    cap.backbone_role = BackboneRole::CarbonylOxygen;
    cap.prochiral = ProchiralStereo::ProS;
    cap.planar_group = PlanarGroupKind::Carboxylate;
    cap.planar_stereo = PlanarStereo::Z;
    cap.pseudoatom = {PseudoatomKind::Q, 6, 2, false};
    cap.polar_h = PolarHKind::AmineNH;
    cap.ring_position.primary = {RingSystemKind::Imidazole_His,
                                 RingPositionLabel::Heteroatom_NH, 5,
                                 true, true, 2};
    cap.aromatic = false;
    cap.formal_charge = 1;
    cap.is_exchangeable = false;
    cap.equivalence_class = 7;

    const AtomMechanicalIdentity original_identity = IdentityOf(chain);
    const auto original_prochiral = chain.prochiral;
    const auto original_planar_stereo = chain.planar_stereo;
    const bool original_aromatic = chain.aromatic;
    const auto original_equivalence_class = chain.equivalence_class;

    gen::ApplyCapDelta(chain, cap);

    EXPECT_EQ(original_identity, IdentityOf(chain));
    EXPECT_EQ(original_prochiral, chain.prochiral);
    EXPECT_EQ(original_planar_stereo, chain.planar_stereo);
    EXPECT_EQ(original_aromatic, chain.aromatic);
    EXPECT_EQ(original_equivalence_class, chain.equivalence_class);

    EXPECT_EQ(cap.planar_group, chain.planar_group);
    EXPECT_TRUE(SamePseudoatom(cap.pseudoatom, chain.pseudoatom));
    EXPECT_EQ(cap.polar_h, chain.polar_h);
    EXPECT_TRUE(SameRingPosition(cap.ring_position, chain.ring_position));
    EXPECT_EQ(cap.formal_charge, chain.formal_charge);
    EXPECT_TRUE(chain.is_exchangeable)
        << "ApplyCapDelta recomputes exchangeability from polar_h instead "
        << "of copying cap_delta.is_exchangeable blindly";
}

TEST(TopologySemanticApi, ApplyCapDeltaRoundTripsRealTerminalOverrides) {
    const auto chain_n_id = Id(Element::N, Locant::None, {0, 0},
                               DiastereotopicIndex::None,
                               BackboneRole::Nitrogen);
    const auto chain_o_id = Id(Element::O, Locant::None, {0, 0},
                               DiastereotopicIndex::None,
                               BackboneRole::CarbonylOxygen);

    const auto* ala_n =
        gen::LookupBy(AminoAcid::ALA, gen::kBaseVariantIdx, chain_n_id);
    const auto* nterm_charged_n =
        gen::LookupCap(TerminalState::NtermCharged, chain_n_id);
    ASSERT_NE(nullptr, ala_n);
    ASSERT_NE(nullptr, nterm_charged_n);

    AtomSemanticTable composed_n = *ala_n;
    const auto original_n_identity = IdentityOf(composed_n);
    const auto original_n_equivalence_class = composed_n.equivalence_class;
    gen::ApplyCapDelta(composed_n, *nterm_charged_n);
    EXPECT_EQ(original_n_identity, IdentityOf(composed_n));
    EXPECT_EQ(original_n_equivalence_class, composed_n.equivalence_class);
    EXPECT_EQ(1, composed_n.formal_charge);
    EXPECT_EQ(nterm_charged_n->planar_group, composed_n.planar_group);

    const auto* ala_o =
        gen::LookupBy(AminoAcid::ALA, gen::kBaseVariantIdx, chain_o_id);
    const auto* cterm_o =
        gen::LookupCap(TerminalState::CtermDeprotonated, chain_o_id);
    ASSERT_NE(nullptr, ala_o);
    ASSERT_NE(nullptr, cterm_o);

    AtomSemanticTable composed_o = *ala_o;
    const auto original_o_identity = IdentityOf(composed_o);
    const auto original_o_equivalence_class = composed_o.equivalence_class;
    gen::ApplyCapDelta(composed_o, *cterm_o);
    EXPECT_EQ(original_o_identity, IdentityOf(composed_o));
    EXPECT_EQ(original_o_equivalence_class, composed_o.equivalence_class);
    EXPECT_EQ(PlanarGroupKind::Carboxylate, composed_o.planar_group);
    EXPECT_EQ(0, composed_o.formal_charge);
}

}  // namespace

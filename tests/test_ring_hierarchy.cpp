#include <gtest/gtest.h>
#include "Ring.h"
#include "SemanticEnums.h"

using namespace nmr;

TEST(RingHierarchy, PheBenzeneProperties) {
    auto ring = CreateRing(RingTypeIndex::PheBenzene);
    EXPECT_DOUBLE_EQ(ring->Intensity(), -12.0);
    EXPECT_DOUBLE_EQ(ring->JBLobeOffset(), 0.64);
    EXPECT_EQ(ring->NitrogenCount(), 0);
    EXPECT_EQ(ring->Aromaticity(), RingAromaticity::Full);
    EXPECT_EQ(ring->RingSizeValue(), 6);
    EXPECT_EQ(std::string(ring->TypeName()), "PHE");
}

TEST(RingHierarchy, TyrPhenolProperties) {
    auto ring = CreateRing(RingTypeIndex::TyrPhenol);
    EXPECT_DOUBLE_EQ(ring->Intensity(), -11.28);
    EXPECT_DOUBLE_EQ(ring->JBLobeOffset(), 0.64);
    EXPECT_EQ(ring->NitrogenCount(), 0);
    EXPECT_EQ(ring->Aromaticity(), RingAromaticity::Full);
    EXPECT_EQ(ring->RingSizeValue(), 6);
}

TEST(RingHierarchy, TrpPyrroleProperties) {
    auto ring = CreateRing(RingTypeIndex::TrpPyrrole);
    EXPECT_DOUBLE_EQ(ring->Intensity(), -6.72);
    EXPECT_DOUBLE_EQ(ring->JBLobeOffset(), 0.52);
    EXPECT_EQ(ring->NitrogenCount(), 1);
    EXPECT_EQ(ring->Aromaticity(), RingAromaticity::Reduced);
    EXPECT_EQ(ring->RingSizeValue(), 5);
}

TEST(RingHierarchy, HisImidazoleProperties) {
    auto ring = CreateRing(RingTypeIndex::HisImidazole);
    EXPECT_DOUBLE_EQ(ring->Intensity(), -5.16);
    EXPECT_DOUBLE_EQ(ring->JBLobeOffset(), 0.50);
    EXPECT_EQ(ring->NitrogenCount(), 2);
    EXPECT_EQ(ring->Aromaticity(), RingAromaticity::Weak);
    EXPECT_EQ(ring->RingSizeValue(), 5);
}

TEST(RingHierarchy, HidImidazoleProperties) {
    auto ring = CreateRing(RingTypeIndex::HidImidazole);
    EXPECT_DOUBLE_EQ(ring->Intensity(), -5.16);
    EXPECT_EQ(ring->NitrogenCount(), 2);
    EXPECT_EQ(std::string(ring->TypeName()), "HID");
}

TEST(RingHierarchy, HieImidazoleProperties) {
    auto ring = CreateRing(RingTypeIndex::HieImidazole);
    EXPECT_DOUBLE_EQ(ring->Intensity(), -5.16);
    EXPECT_EQ(std::string(ring->TypeName()), "HIE");
}

TEST(RingHierarchy, IndolePerimeterProperties) {
    auto ring = CreateRing(RingTypeIndex::TrpPerimeter);
    EXPECT_DOUBLE_EQ(ring->Intensity(), -19.2);
    EXPECT_DOUBLE_EQ(ring->JBLobeOffset(), 0.60);
    EXPECT_EQ(ring->NitrogenCount(), 1);
    EXPECT_EQ(ring->Aromaticity(), RingAromaticity::Full);
    EXPECT_EQ(ring->RingSizeValue(), 9);
}

TEST(RingHierarchy, TrpBenzeneProperties) {
    auto ring = CreateRing(RingTypeIndex::TrpBenzene);
    EXPECT_DOUBLE_EQ(ring->Intensity(), -12.48);
    EXPECT_EQ(ring->RingSizeValue(), 6);
}

TEST(RingHierarchy, AllRingTypesCreatable) {
    for (int i = 0; i < static_cast<int>(RingTypeIndex::Count); ++i) {
        auto ring = CreateRing(static_cast<RingTypeIndex>(i));
        EXPECT_NE(ring, nullptr);
        EXPECT_EQ(ring->TypeIndexAsInt(), i);
    }
}

TEST(RingHierarchy, AromaticTypeCountBoundary) {
    // The named constant codifies the boundary between aromatic ring
    // chemistries (indices 0..kAromaticRingTypeCount-1) and saturated
    // ring chemistries. Pro pyrrolidine is the canonical saturated
    // case; calculator per-aromatic-type accumulator arrays gate Pro
    // out via `if (ti < kAromaticRingTypeCount)`.
    EXPECT_EQ(kAromaticRingTypeCount, 8);
    EXPECT_EQ(static_cast<int>(RingTypeIndex::ProPyrrolidine), kAromaticRingTypeCount);
    EXPECT_GE(static_cast<int>(RingTypeIndex::Count), kAromaticRingTypeCount + 1);
}

TEST(RingHierarchy, ProPyrrolidineProperties) {
    // Pro pyrrolidine: saturated 5-ring; ring current is identically
    // zero (Joule & Mills 2010 ch. 7). Literal 0.0 values, not
    // CalculatorConfig — physics, not calibration.
    auto ring = CreateRing(RingTypeIndex::ProPyrrolidine);
    EXPECT_DOUBLE_EQ(ring->Intensity(), 0.0);
    EXPECT_DOUBLE_EQ(ring->LiteratureIntensity(), 0.0);
    EXPECT_DOUBLE_EQ(ring->JBLobeOffset(), 0.0);
    EXPECT_EQ(ring->NitrogenCount(), 1);
    EXPECT_EQ(ring->Aromaticity(), RingAromaticity::None);
    EXPECT_EQ(ring->RingSizeValue(), 5);
}

TEST(RingHierarchy, PolymorphicCall) {
    // Verify polymorphism: call through base pointer
    std::unique_ptr<Ring> ring = CreateRing(RingTypeIndex::PheBenzene);
    EXPECT_DOUBLE_EQ(ring->Intensity(), -12.0);

    ring = CreateRing(RingTypeIndex::HisImidazole);
    EXPECT_DOUBLE_EQ(ring->Intensity(), -5.16);
}

TEST(RingHierarchy, ComputeGeometryRegularHexagon) {
    auto ring = CreateRing(RingTypeIndex::PheBenzene);
    // Create a regular hexagon in the xy-plane
    double R = 1.4;  // Benzene C-C distance
    ring->atom_indices = {0, 1, 2, 3, 4, 5};
    std::vector<Vec3> positions(6);
    for (int i = 0; i < 6; ++i) {
        double angle = i * 2.0 * 3.14159265358979 / 6.0;
        positions[i] = Vec3(R * std::cos(angle), R * std::sin(angle), 0.0);
    }

    auto geo = ring->ComputeGeometry(positions);

    EXPECT_NEAR(geo.center.norm(), 0.0, 1e-10);
    EXPECT_NEAR(std::abs(geo.normal.z()), 1.0, 1e-10);
    EXPECT_NEAR(geo.radius, R, 1e-10);
    EXPECT_EQ(geo.vertices.size(), 6);
}


// ============================================================================
// Substrate predicate semantics — Slice A added the tertiary RingMembership
// slot for atoms in three rings (TRP indole bridgeheads CE2/CD2 carry
// primary=Indole_Trp_5 + secondary=Indole_Trp_6 + tertiary=Indole_Trp_9;
// non-bridgehead perimeter atoms have primary + tertiary populated with
// secondary EMPTY). The new predicates HasPrimaryRing/HasSecondaryRing/
// HasTertiaryRing/MembershipCount/IsInAnyRing are slot-based; they
// replace the old InAnyRing/InTwoRings/InThreeRings predicates that
// returned the wrong answer for the non-bridgehead perimeter case
// (InTwoRings would return false even though the atom is in two rings,
// because secondary is the empty slot).
// ============================================================================

namespace {
constexpr RingMembership MakeRing(RingSystemKind kind,
                                  RingPositionLabel pos,
                                  uint8_t size) {
    return RingMembership{kind, pos, size, /*aromatic*/ true,
                          /*planar*/ true, /*n_het*/ 1};
}
}  // namespace

TEST(RingPositionPredicates, IsPopulatedDistinguishesDefaultFromPopulated) {
    RingMembership empty{};
    RingMembership populated = MakeRing(RingSystemKind::Indole_Trp_5,
                                        RingPositionLabel::Ipso, 5);
    EXPECT_FALSE(empty.IsPopulated());
    EXPECT_TRUE(populated.IsPopulated());
}

TEST(RingPositionPredicates, MembershipCountAndSlotPredicatesNoRing) {
    RingPosition rp{};
    EXPECT_FALSE(rp.HasPrimaryRing());
    EXPECT_FALSE(rp.HasSecondaryRing());
    EXPECT_FALSE(rp.HasTertiaryRing());
    EXPECT_EQ(rp.MembershipCount(), 0);
    EXPECT_FALSE(rp.IsInAnyRing());
}

TEST(RingPositionPredicates, MembershipCountAndSlotPredicatesPheRingAtom) {
    // Phe ring atom: primary populated, secondary + tertiary empty.
    // Models any aromatic-ring atom that's not in a fused or perimeter
    // construction.
    RingPosition rp{};
    rp.primary = MakeRing(RingSystemKind::Benzene_Phe,
                          RingPositionLabel::Ipso, 6);
    EXPECT_TRUE(rp.HasPrimaryRing());
    EXPECT_FALSE(rp.HasSecondaryRing());
    EXPECT_FALSE(rp.HasTertiaryRing());
    EXPECT_EQ(rp.MembershipCount(), 1);
    EXPECT_TRUE(rp.IsInAnyRing());
}

TEST(RingPositionPredicates, MembershipCountAndSlotPredicatesTrpBridgehead) {
    // TRP CE2 / CD2: in 5-ring (primary, smaller-ring convention),
    // 6-ring (secondary), AND the 9-atom indole perimeter (tertiary).
    // MembershipCount must report 3.
    RingPosition rp{};
    rp.primary   = MakeRing(RingSystemKind::Indole_Trp_5,
                            RingPositionLabel::BridgeFusion, 5);
    rp.secondary = MakeRing(RingSystemKind::Indole_Trp_6,
                            RingPositionLabel::BridgeFusion, 6);
    rp.tertiary  = MakeRing(RingSystemKind::Indole_Trp_9,
                            RingPositionLabel::PerimeterMember, 9);
    EXPECT_TRUE(rp.HasPrimaryRing());
    EXPECT_TRUE(rp.HasSecondaryRing());
    EXPECT_TRUE(rp.HasTertiaryRing());
    EXPECT_EQ(rp.MembershipCount(), 3);
    EXPECT_TRUE(rp.IsInAnyRing());
}

TEST(RingPositionPredicates, MembershipCountAndSlotPredicatesTrpNonBridgePerimeter) {
    // TRP CG / CD1 / NE1 / CE3 / CZ2 / CH2 / CZ3: in primary (5- or
    // 6-ring) + tertiary (9-perimeter); secondary is EMPTY. This is
    // the case where the old InTwoRings() predicate returned the
    // wrong answer (secondary empty → "false", but the atom IS in
    // two rings via primary + tertiary). MembershipCount must report
    // 2; HasSecondaryRing must report false; HasTertiaryRing true.
    RingPosition rp{};
    rp.primary  = MakeRing(RingSystemKind::Indole_Trp_5,
                           RingPositionLabel::Ipso, 5);
    rp.tertiary = MakeRing(RingSystemKind::Indole_Trp_9,
                           RingPositionLabel::PerimeterMember, 9);
    EXPECT_TRUE(rp.HasPrimaryRing());
    EXPECT_FALSE(rp.HasSecondaryRing());
    EXPECT_TRUE(rp.HasTertiaryRing());
    EXPECT_EQ(rp.MembershipCount(), 2);
    EXPECT_TRUE(rp.IsInAnyRing());
}

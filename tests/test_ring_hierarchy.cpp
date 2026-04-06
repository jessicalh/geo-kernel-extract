#include <gtest/gtest.h>
#include "Ring.h"

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

TEST(RingHierarchy, AllEightTypesCreatable) {
    for (int i = 0; i < static_cast<int>(RingTypeIndex::Count); ++i) {
        auto ring = CreateRing(static_cast<RingTypeIndex>(i));
        EXPECT_NE(ring, nullptr);
        EXPECT_EQ(ring->TypeIndexAsInt(), i);
    }
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

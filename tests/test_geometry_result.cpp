#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include "GeometryResult.h"
#include "PdbFileReader.h"
#include <filesystem>

using namespace nmr;


class GeometryResultTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!std::filesystem::exists(nmr::test::TestEnvironment::UbqProtonated())) {
            GTEST_SKIP() << "1UBQ PDB not found";
        }
        auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
        if (!r.Ok()) GTEST_SKIP() << r.error;
        protein = std::move(r.protein);
    }
    std::unique_ptr<Protein> protein;
};

TEST_F(GeometryResultTest, ComputeAndAttach) {
    auto& conf = protein->Conformation();
    auto geo = GeometryResult::Compute(conf);
    EXPECT_TRUE(conf.AttachResult(std::move(geo)));
    EXPECT_TRUE(conf.HasResult<GeometryResult>());
}

TEST_F(GeometryResultTest, RingGeometryComputed) {
    auto& conf = protein->Conformation();
    auto geo = GeometryResult::Compute(conf);
    conf.AttachResult(std::move(geo));

    for (size_t ri = 0; ri < protein->RingCount(); ++ri) {
        const auto& rg = conf.ring_geometries[ri];
        // Center should not be zero (rings have non-zero positions)
        EXPECT_GT(rg.center.norm(), 0.1);
        // Normal should be unit length
        EXPECT_NEAR(rg.normal.norm(), 1.0, 1e-6);
        // Radius should be positive and reasonable (1-3 Angstroms)
        EXPECT_GT(rg.radius, 0.5);
        EXPECT_LT(rg.radius, 5.0);
    }
}

TEST_F(GeometryResultTest, BondGeometryComputed) {
    auto& conf = protein->Conformation();
    auto geo = GeometryResult::Compute(conf);
    conf.AttachResult(std::move(geo));

    // Bond lengths should be in [0.5, 3.0] Angstroms
    for (size_t bi = 0; bi < protein->BondCount(); ++bi) {
        EXPECT_GT(conf.bond_lengths[bi], 0.3);
        EXPECT_LT(conf.bond_lengths[bi], 3.5);
    }
}

TEST_F(GeometryResultTest, GlobalGeometry) {
    auto& conf = protein->Conformation();
    auto geo = GeometryResult::Compute(conf);
    conf.AttachResult(std::move(geo));

    EXPECT_GT(conf.radius_of_gyration, 1.0);
    EXPECT_GT(conf.center_of_geometry.norm(), 0.0);
}

TEST_F(GeometryResultTest, PrebuiltCollections) {
    auto& conf = protein->Conformation();
    auto geo = GeometryResult::Compute(conf);
    conf.AttachResult(std::move(geo));

    // Should have at least some rings by type
    size_t total_rings = 0;
    for (auto& [type, indices] : conf.rings_by_type)
        total_rings += indices.size();
    EXPECT_EQ(total_rings, protein->RingCount());

    // Should have residues by type
    EXPECT_FALSE(conf.residues_by_type.empty());
}

TEST_F(GeometryResultTest, RingPairs) {
    auto& conf = protein->Conformation();
    auto geo = GeometryResult::Compute(conf);
    conf.AttachResult(std::move(geo));

    size_t n = protein->RingCount();
    EXPECT_EQ(conf.ring_pairs.size(), n * (n - 1) / 2);
}

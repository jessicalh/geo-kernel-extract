#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include "GeometryResult.h"
#include "DemoResult.h"
#include "PdbFileReader.h"
#include <filesystem>
#include <random>

using namespace nmr;


class TwoConformationsTest : public ::testing::Test {
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

TEST_F(TwoConformationsTest, TwoConformationsIndependent) {
    // Crystal conformation already exists from loading
    auto& conf1 = protein->Conformation();

    // Create a second conformation by jittering positions
    std::vector<Vec3> jittered = conf1.Positions();
    std::mt19937 rng(42);
    std::normal_distribution<double> noise(0.0, 0.1);
    for (auto& p : jittered) {
        p.x() += noise(rng);
        p.y() += noise(rng);
        p.z() += noise(rng);
    }

    auto& conf2 = protein->AddDerived(conf1, "jittered_0.1A", std::move(jittered));

    EXPECT_EQ(protein->ConformationCount(), 2);

    // Attach GeometryResult and DemoResult to BOTH
    conf1.AttachResult(GeometryResult::Compute(conf1));
    conf1.AttachResult(DemoResult::Compute(conf1));

    conf2.AttachResult(GeometryResult::Compute(conf2));
    conf2.AttachResult(DemoResult::Compute(conf2));

    // Both should have results
    EXPECT_TRUE(conf1.HasResult<GeometryResult>());
    EXPECT_TRUE(conf1.HasResult<DemoResult>());
    EXPECT_TRUE(conf2.HasResult<GeometryResult>());
    EXPECT_TRUE(conf2.HasResult<DemoResult>());

    // Results should be DIFFERENT (positions are different)
    auto& demo1 = conf1.Result<DemoResult>();
    auto& demo2 = conf2.Result<DemoResult>();

    // At least some atoms should have different nearest ring distances
    int diff_count = 0;
    for (size_t ai = 0; ai < conf1.AtomCount(); ++ai) {
        double d1 = demo1.NearestRingDistance(ai);
        double d2 = demo2.NearestRingDistance(ai);
        if (std::abs(d1 - d2) > 0.001) diff_count++;
    }
    EXPECT_GT(diff_count, 0) << "Jittered conformation should differ";

    // Ring geometry should differ
    if (protein->RingCount() > 0) {
        Vec3 c1 = conf1.ring_geometries[0].center;
        Vec3 c2 = conf2.ring_geometries[0].center;
        EXPECT_GT((c1 - c2).norm(), 0.001);
    }

    // Global geometry should differ
    double rg1 = conf1.radius_of_gyration;
    double rg2 = conf2.radius_of_gyration;
    EXPECT_GT(std::abs(rg1 - rg2), 0.0);
}

TEST_F(TwoConformationsTest, ResultsAccumulateIndependently) {
    auto& conf1 = protein->Conformation();

    std::vector<Vec3> jittered = conf1.Positions();
    for (auto& p : jittered) p.x() += 0.5;

    auto& conf2 = protein->AddDerived(conf1, "shifted_x", std::move(jittered));

    // Attach only GeometryResult to conf1
    conf1.AttachResult(GeometryResult::Compute(conf1));
    EXPECT_TRUE(conf1.HasResult<GeometryResult>());
    EXPECT_FALSE(conf1.HasResult<DemoResult>());

    // Attach both to conf2
    conf2.AttachResult(GeometryResult::Compute(conf2));
    conf2.AttachResult(DemoResult::Compute(conf2));
    EXPECT_TRUE(conf2.HasResult<GeometryResult>());
    EXPECT_TRUE(conf2.HasResult<DemoResult>());

    // conf1 should still not have DemoResult
    EXPECT_FALSE(conf1.HasResult<DemoResult>());
}

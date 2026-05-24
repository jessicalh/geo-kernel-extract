#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include "DemoResult.h"
#include "GeometryResult.h"
#include "PdbFileReader.h"
#include <filesystem>

using namespace nmr;


class DemoResultTest : public ::testing::Test {
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

TEST_F(DemoResultTest, RequiresGeometry) {
    auto& conf = protein->Conformation();

    // DemoResult requires GeometryResult -- should fail without it
    auto demo = DemoResult::Compute(conf);
    EXPECT_FALSE(conf.AttachResult(std::move(demo)));
}

TEST_F(DemoResultTest, ComputeAndAttach) {
    auto& conf = protein->Conformation();

    conf.AttachResult(GeometryResult::Compute(conf));
    auto demo = DemoResult::Compute(conf);
    EXPECT_TRUE(conf.AttachResult(std::move(demo)));
    EXPECT_TRUE(conf.HasResult<DemoResult>());
}

TEST_F(DemoResultTest, NearestRingDistances) {
    auto& conf = protein->Conformation();

    conf.AttachResult(GeometryResult::Compute(conf));
    conf.AttachResult(DemoResult::Compute(conf));

    auto& demo = conf.Result<DemoResult>();

    // All atoms should have some distance value
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        double d = demo.NearestRingDistance(ai);
        EXPECT_GE(d, 0.0);
    }

    // Atoms in an aromatic residue should be near a ring
    for (size_t ri = 0; ri < protein->RingCount(); ++ri) {
        const auto& ring = protein->RingAt(ri);
        for (size_t ai : ring.atom_indices) {
            // Ring atoms should be very close to their own ring center
            double d = demo.NearestRingDistance(ai);
            EXPECT_LT(d, 5.0) << "Ring atom " << ai << " too far from ring center";
        }
    }
}

TEST_F(DemoResultTest, SphericalTensorRoundtrip) {
    auto& conf = protein->Conformation();

    conf.AttachResult(GeometryResult::Compute(conf));
    conf.AttachResult(DemoResult::Compute(conf));

    auto& demo = conf.Result<DemoResult>();

    // The test decomposition should roundtrip
    Mat3 original;
    original << 1.0, 0.5, 0.3,
                0.2, 2.0, 0.7,
                0.1, 0.4, 3.0;

    Mat3 reconstructed = demo.TestReconstructed();

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(reconstructed(i,j), original(i,j), 1e-12)
                << "Mismatch at (" << i << "," << j << ")";
}

TEST_F(DemoResultTest, SphericalTensorComponents) {
    auto& conf = protein->Conformation();

    conf.AttachResult(GeometryResult::Compute(conf));
    conf.AttachResult(DemoResult::Compute(conf));

    auto& demo = conf.Result<DemoResult>();
    const auto& st = demo.TestDecomposition();

    // T0 = trace/3 = (1+2+3)/3 = 2.0
    EXPECT_NEAR(st.T0, 2.0, 1e-12);

    // T1 should be non-zero (asymmetric matrix)
    double T1_mag = 0;
    for (double v : st.T1) T1_mag += v * v;
    EXPECT_GT(T1_mag, 0.0);
}

TEST_F(DemoResultTest, ConformationAtomFieldsPopulated) {
    auto& conf = protein->Conformation();

    conf.AttachResult(GeometryResult::Compute(conf));
    conf.AttachResult(DemoResult::Compute(conf));

    // Check that DemoResult wrote to ConformationAtom fields
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& ca = conf.AtomAt(ai);
        // demo_nearest_ring_distance should have been set
        // (either a real distance or sentinel)
        if (protein->RingCount() > 0) {
            EXPECT_GT(ca.demo_nearest_ring_distance, 0.0);
        }
    }
}

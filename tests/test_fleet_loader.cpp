#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include "GromacsEnsembleLoader.h"
#include "ConformationResult.h"
#include "ChargeAssignmentResult.h"
#include "OperationRunner.h"
#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "EnrichmentResult.h"
#include "ApbsFieldResult.h"
#include "BiotSavartResult.h"

#include <filesystem>

namespace fs = std::filesystem;
using namespace nmr;

// Fleet test data lives in the repo

class FleetLoaderTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!fs::exists(nmr::test::TestEnvironment::FleetData())) {
            GTEST_SKIP() << "Fleet test data not found at " << nmr::test::TestEnvironment::FleetData();
        }
    }

    FleetPaths PathsFor(const std::string& protein_id) {
        FleetPaths p;
        p.sampled_poses_dir = std::string(nmr::test::TestEnvironment::FleetData()) + "/" + protein_id + "/poses";
        p.tpr_path = std::string(nmr::test::TestEnvironment::FleetData()) + "/" + protein_id + "/params/prod.tpr";
        p.force_field = ForceField::CHARMM36m;
        return p;
    }
};


TEST_F(FleetLoaderTest, LoadsSuccessfully) {
    auto result = BuildFromGromacs(PathsFor("1A6J_5789"));
    ASSERT_TRUE(result.ok) << result.error;
    EXPECT_GT(result.protein->AtomCount(), 100u);
    EXPECT_GT(result.protein->ResidueCount(), 10u);
}

TEST_F(FleetLoaderTest, HasTenFrames) {
    auto result = BuildFromGromacs(PathsFor("1A6J_5789"));
    ASSERT_TRUE(result.ok) << result.error;
    EXPECT_EQ(result.protein->MDFrameCount(), 10u);
    EXPECT_EQ(result.protein->ConformationCount(), 10u);
}

TEST_F(FleetLoaderTest, FrameMetadataPresent) {
    auto result = BuildFromGromacs(PathsFor("1A6J_5789"));
    ASSERT_TRUE(result.ok) << result.error;

    auto& frame0 = result.protein->MDFrameAt(0);
    EXPECT_GE(frame0.Walker(), 0);
    EXPECT_GT(frame0.TimePicoseconds(), 0.0);
    EXPECT_GT(frame0.BoltzmannWeight(), 0.0);
    EXPECT_GT(frame0.RmsdNanometres(), 0.0);
    EXPECT_GT(frame0.RadiusOfGyrationNm(), 0.0);
}

TEST_F(FleetLoaderTest, AtomCountConsistentAcrossFrames) {
    auto result = BuildFromGromacs(PathsFor("1A6J_5789"));
    ASSERT_TRUE(result.ok) << result.error;

    for (size_t i = 0; i < result.protein->MDFrameCount(); ++i) {
        EXPECT_EQ(result.protein->MDFrameAt(i).AtomCount(),
                  result.protein->AtomCount())
            << "Frame " << i << " atom count mismatch";
    }
}

TEST_F(FleetLoaderTest, PositionsDifferBetweenFrames) {
    auto result = BuildFromGromacs(PathsFor("1A6J_5789"));
    ASSERT_TRUE(result.ok) << result.error;
    ASSERT_GE(result.protein->MDFrameCount(), 2u);

    auto& f0 = result.protein->MDFrameAt(0);
    auto& f1 = result.protein->MDFrameAt(1);

    bool any_differ = false;
    for (size_t i = 0; i < f0.AtomCount(); ++i) {
        if ((f0.PositionAt(i) - f1.PositionAt(i)).norm() > 0.01) {
            any_differ = true;
            break;
        }
    }
    EXPECT_TRUE(any_differ) << "Frames 0 and 1 have identical positions";
}

TEST_F(FleetLoaderTest, TopologySharedAcrossFrames) {
    auto result = BuildFromGromacs(PathsFor("1A6J_5789"));
    ASSERT_TRUE(result.ok) << result.error;

    EXPECT_GT(result.protein->BondCount(), 100u);
    EXPECT_GT(result.protein->RingCount(), 0u);
}

TEST_F(FleetLoaderTest, ChargesReturned) {
    auto result = BuildFromGromacs(PathsFor("1A6J_5789"));
    ASSERT_TRUE(result.ok) << result.error;
    ASSERT_NE(result.charges, nullptr);

    // Verify charges load through ChargeAssignmentResult
    auto& conf = result.protein->Conformation();
    conf.AttachResult(GeometryResult::Compute(conf));
    auto ca = ChargeAssignmentResult::Compute(conf, *result.charges);
    ASSERT_NE(ca, nullptr);
    conf.AttachResult(std::move(ca));

    double sum = 0.0;
    for (size_t i = 0; i < conf.AtomCount(); ++i)
        sum += std::abs(conf.AtomAt(i).partial_charge);
    EXPECT_GT(sum, 0.0);
}

TEST_F(FleetLoaderTest, ConformationAccessOrthogonal) {
    auto result = BuildFromGromacs(PathsFor("1A6J_5789"));
    ASSERT_TRUE(result.ok) << result.error;

    auto& primary = result.protein->Conformation();
    auto& frame0 = result.protein->MDFrameAt(0);
    EXPECT_EQ(&primary, static_cast<ProteinConformation*>(&frame0));
}

TEST_F(FleetLoaderTest, FoundationResultsAttach) {
    auto result = BuildFromGromacs(PathsFor("1A6J_5789"));
    ASSERT_TRUE(result.ok) << result.error;

    auto& frame = result.protein->Conformation();
    ASSERT_TRUE(frame.AttachResult(GeometryResult::Compute(frame)));
    ASSERT_TRUE(frame.AttachResult(SpatialIndexResult::Compute(frame)));
    ASSERT_TRUE(frame.AttachResult(EnrichmentResult::Compute(frame)));

    EXPECT_TRUE(frame.HasResult<GeometryResult>());
    EXPECT_TRUE(frame.HasResult<SpatialIndexResult>());
    EXPECT_TRUE(frame.HasResult<EnrichmentResult>());
}

TEST_F(FleetLoaderTest, HisVariantDetected) {
    // 1AEP has HSP (doubly protonated) and plain "HIS" (epsilon)
    auto paths = PathsFor("1AEP_4814");
    if (!fs::exists(paths.sampled_poses_dir)) GTEST_SKIP();

    auto result = BuildFromGromacs(paths);
    ASSERT_TRUE(result.ok) << result.error;

    bool found_doubly = false;
    bool found_epsilon = false;
    for (size_t ri = 0; ri < result.protein->ResidueCount(); ++ri) {
        const auto& res = result.protein->ResidueAt(ri);
        if (res.type == AminoAcid::HIS) {
            if (res.protonation_variant_index == 2) found_doubly = true;
            if (res.protonation_variant_index == 1) found_epsilon = true;
        }
    }
    EXPECT_TRUE(found_doubly) << "Expected HSP (doubly protonated HIS)";
    EXPECT_TRUE(found_epsilon) << "Expected HSE (epsilon protonated HIS)";
}


// ============================================================================
// Full pipeline: load, run all 8 calculators on all 10 frames, write features
// ============================================================================

TEST_F(FleetLoaderTest, FullPipelineAllFrames) {
    auto result = BuildFromGromacs(PathsFor("1A6J_5789"));
    ASSERT_TRUE(result.ok) << result.error;
    ASSERT_EQ(result.protein->MDFrameCount(), 10u);

    // Full pipeline options: charges, APBS, MOPAC, all calculators
    RunOptions opts;
    opts.charge_source = result.charges.get();
    opts.net_charge = result.net_charge;

    // Run full pipeline on all 10 frames
    auto run_results = RunAllFrames(*result.protein, opts);
    ASSERT_EQ(run_results.size(), 10u);

    // Every frame should have: foundation (4) + charges + MOPAC + APBS + 8 calculators = 15
    for (size_t i = 0; i < run_results.size(); ++i) {
        EXPECT_GE(run_results[i].attached.size(), 13u)
            << "Frame " << i << ": expected 13+ results, got " << run_results[i].attached.size();
    }

    // Verify specific results on frame 0
    auto& f0 = result.protein->MDFrameAt(0);
    EXPECT_TRUE(f0.HasResult<BiotSavartResult>());
    EXPECT_TRUE(f0.HasResult<ApbsFieldResult>())
        << "APBS must run — it's the solvated E-field the model needs";

    // Write features for each frame
    std::string output_base = "/tmp/fleet_test_output/1A6J_5789";
    fs::create_directories(output_base);

    int total_arrays = 0;
    for (size_t i = 0; i < result.protein->MDFrameCount(); ++i) {
        auto& frame = result.protein->MDFrameAt(i);
        char frame_dir[256];
        std::snprintf(frame_dir, sizeof(frame_dir), "%s/frame_%03zu",
                      output_base.c_str(), i + 1);
        fs::create_directories(frame_dir);

        int arrays = ConformationResult::WriteAllFeatures(frame, frame_dir);
        EXPECT_GT(arrays, 25) << "Frame " << i << ": too few arrays written";
        total_arrays += arrays;
    }

    std::cout << "\n  Fleet pipeline complete: 1A6J_5789\n"
              << "  Frames: " << result.protein->MDFrameCount() << "\n"
              << "  Atoms: " << result.protein->AtomCount() << "\n"
              << "  Total NPY arrays: " << total_arrays << "\n"
              << "  Results per frame: " << run_results[0].attached.size() << "\n"
              << "  Attached: ";
    for (const auto& name : run_results[0].attached) std::cout << name << " ";
    std::cout << "\n";

    // Clean up
    fs::remove_all("/tmp/fleet_test_output");
}

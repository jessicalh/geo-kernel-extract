//
// test_smoke_fes_fleet.cpp
//
// Smoke test for fleet extraction with fes-sampler naming convention.
//
// Loads 2 small proteins (1I8X_4351: 422 atoms, 1HD6_4820: 526 atoms)
// with descriptive PDB names (boltzmann_minimum, max_dev_phi_*, etc.)
// via BuildFromGromacs, runs the full pipeline on all 6 poses each,
// and verifies that output directories carry the source PDB names.
//
// This exercises the production extraction path: descriptive pose names
// from fes-sampler → CHARMM36m charges from TPR → all extractors
// (including MOPAC) → named output directories with NPY files.
//
// Runtime: ~15-20 minutes (12 MOPAC runs on small proteins).
// Run alone — do not combine with other MOPAC-heavy tests.
//
//   ./fes_fleet_smoke_tests
//
// Separate executable to avoid inflating the main test suite.
//

#include "TestEnvironment.h"
#include <gtest/gtest.h>

#include "GromacsEnsembleLoader.h"
#include "ConformationResult.h"
#include "OperationRunner.h"
#include "OperationLog.h"
#include "RuntimeEnvironment.h"

#include <filesystem>
#include <fstream>
#include <set>
#include <algorithm>

namespace fs = std::filesystem;
using namespace nmr;


static std::set<std::string> NpyFiles(const std::string& dir) {
    std::set<std::string> result;
    if (!fs::exists(dir)) return result;
    for (const auto& entry : fs::directory_iterator(dir)) {
        if (entry.path().extension() == ".npy")
            result.insert(entry.path().filename().string());
    }
    return result;
}

static bool VerifyNpyMagic(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    char buf[6];
    if (!in.read(buf, 6)) return false;
    return buf[0] == '\x93' && buf[1] == 'N' && buf[2] == 'U'
        && buf[3] == 'M' && buf[4] == 'P' && buf[5] == 'Y';
}


// Run one protein through the fleet pipeline and validate output naming.
static void RunFesProtein(const std::string& protein_id,
                          const std::string& fes_base,
                          const std::string& output_base,
                          size_t expected_poses) {
    std::string poses_dir = fes_base + "/" + protein_id + "/poses";
    std::string tpr_path  = fes_base + "/" + protein_id + "/params/prod.tpr";

    if (!fs::exists(tpr_path)) {
        GTEST_SKIP() << protein_id << " TPR not found at " << tpr_path;
    }
    if (!fs::exists(poses_dir)) {
        GTEST_SKIP() << protein_id << " poses not found at " << poses_dir;
    }

    // ---- Load ----
    FleetPaths paths;
    paths.sampled_poses_dir = poses_dir;
    paths.tpr_path = tpr_path;
    paths.force_field = ForceField::CHARMM36m;

    auto build = BuildFromGromacs(paths);
    ASSERT_TRUE(build.Ok()) << protein_id << ": " << build.error;
    ASSERT_EQ(build.protein->MDFrameCount(), expected_poses)
        << protein_id << ": expected " << expected_poses << " poses";

    // Pose names should be populated (descriptive naming, not pose_NNN)
    ASSERT_EQ(build.pose_names.size(), expected_poses)
        << protein_id << ": pose_names not populated";

    // Every pose name should contain the protein_id
    for (const auto& name : build.pose_names) {
        EXPECT_NE(name.find(protein_id), std::string::npos)
            << "Pose name missing protein_id: " << name;
    }

    std::cout << "\n  " << protein_id << ": "
              << build.protein->AtomCount() << " atoms, "
              << build.protein->RingCount() << " rings, "
              << expected_poses << " poses\n";
    for (const auto& name : build.pose_names) {
        std::cout << "    " << name << "\n";
    }

    // ---- Run pipeline ----
    RunOptions opts;
    opts.charge_source = build.charges.get();
    opts.net_charge = build.net_charge;

    auto run_results = RunAllFrames(*build.protein, opts);
    ASSERT_EQ(run_results.size(), expected_poses);
    for (size_t i = 0; i < run_results.size(); ++i) {
        EXPECT_TRUE(run_results[i].Ok())
            << protein_id << " pose " << i << ": " << run_results[i].error;
        EXPECT_GE(run_results[i].attached.size(), 13u)
            << protein_id << " pose " << i << ": too few results";
    }

    // ---- Write features with named output dirs ----
    std::string protein_out = output_base + "/" + protein_id;
    fs::create_directories(protein_out);

    int total_arrays = 0;
    for (size_t i = 0; i < build.protein->MDFrameCount(); ++i) {
        auto& frame = build.protein->MDFrameAt(i);

        // Build output dir: {pose_name}_frame{NNN}
        std::string subdir;
        if (i < build.pose_names.size()) {
            char tag[32];
            std::snprintf(tag, sizeof(tag), "_frame%03zu", i + 1);
            subdir = protein_out + "/" + build.pose_names[i] + tag;
        } else {
            char frame_dir[512];
            std::snprintf(frame_dir, sizeof(frame_dir), "%s/frame_%03zu",
                          protein_out.c_str(), i + 1);
            subdir = frame_dir;
        }
        fs::create_directories(subdir);

        int arrays = ConformationResult::WriteAllFeatures(frame, subdir);
        EXPECT_GE(arrays, 35) << protein_id << " pose " << i << ": too few arrays";
        total_arrays += arrays;

        // Validate NPY files
        auto files = NpyFiles(subdir);
        for (const auto& name : files) {
            std::string path = subdir + "/" + name;
            EXPECT_GT(fs::file_size(path), 0u) << name << " is empty";
            EXPECT_TRUE(VerifyNpyMagic(path)) << name << " bad NPY magic";
        }

        // Verify ring_contributions has 59 columns
        std::string rc_path = subdir + "/ring_contributions.npy";
        if (fs::exists(rc_path)) {
            // Quick header check: read shape from NPY header
            std::ifstream rc(rc_path, std::ios::binary);
            char hdr[128];
            rc.read(hdr, 128);
            std::string hdr_str(hdr, 128);
            EXPECT_NE(hdr_str.find("59"), std::string::npos)
                << "ring_contributions missing 59 columns in " << subdir;
        }
    }

    // Verify output directory names are descriptive (not frame_NNN)
    for (const auto& entry : fs::directory_iterator(protein_out)) {
        if (!entry.is_directory()) continue;
        std::string dirname = entry.path().filename().string();
        EXPECT_NE(dirname.find(protein_id), std::string::npos)
            << "Output dir should contain protein_id: " << dirname;
        EXPECT_NE(dirname.find("_frame"), std::string::npos)
            << "Output dir should contain _frame tag: " << dirname;
    }

    std::cout << "  " << protein_id << ": " << total_arrays
              << " arrays across " << expected_poses << " poses\n"
              << "  Output: " << protein_out << "\n";
}


TEST(SmokeFesFleet, NamedPoses) {
    RuntimeEnvironment::Load();

    std::string fes_base;
#ifdef NMR_TEST_DATA_DIR
    fes_base = std::string(NMR_TEST_DATA_DIR) + "/fes_fleet";
#endif
    if (fes_base.empty() || !fs::exists(fes_base)) {
        GTEST_SKIP() << "FES fleet test data not found at " << fes_base;
    }

    std::string output_base;
#ifdef NMR_TEST_DATA_DIR
    output_base = std::string(NMR_TEST_DATA_DIR) + "/../golden/smoke/fes_fleet";
#endif
    fs::create_directories(output_base);

    std::string log_path = output_base + "/log.jsonl";
    OperationLog::ConfigureFile(log_path);

    RunFesProtein("1I8X_4351", fes_base, output_base, 6);
    RunFesProtein("1HD6_4820", fes_base, output_base, 6);

    OperationLog::CloseFile();
}

//
// test_gromacs_streaming: tests for GromacsProtein fleet path and
// GromacsFrameHandler trajectory path.
//

#include "GromacsProtein.h"
#include "GromacsFrameHandler.h"
#include "GromacsFinalResult.h"
#include "OperationRunner.h"
#include "OperationLog.h"

#include <gtest/gtest.h>
#include <filesystem>
#include <iostream>
#include <set>

namespace fs = std::filesystem;

#ifndef NMR_TEST_DATA_DIR
#error "NMR_TEST_DATA_DIR must be defined"
#endif

static const std::string FLEET_DIR =
    std::string(NMR_TEST_DATA_DIR) + "/fleet_test_large";

static const std::string FULLSYS_DIR =
    std::string(NMR_TEST_DATA_DIR) + "/fleet_test_fullsys/1ZR7_6721";


// ============================================================================
// Fleet path: Build from pre-extracted PDB poses
// ============================================================================

TEST(GromacsStreaming, FleetBuildAndCompute) {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);

    nmr::FleetPaths paths;
    paths.sampled_poses_dir = FLEET_DIR + "/initial-samples";
    paths.tpr_path = FLEET_DIR + "/run_params/prod.tpr";
    paths.force_field = nmr::ForceField::CHARMM36m;

    nmr::GromacsProtein gp;
    ASSERT_TRUE(gp.Build(paths)) << gp.error();

    ASSERT_GT(gp.protein().AtomCount(), 0u);
    ASSERT_EQ(gp.AtomCount(), gp.protein().AtomCount());

    // Run calculators on conformation 0
    nmr::RunOptions opts;
    opts.skip_mopac = true;
    opts.skip_coulomb = true;
    opts.charge_source = gp.charges();

    auto& conf0 = gp.protein().ConformationAt(0);
    nmr::RunResult rr = nmr::OperationRunner::Run(conf0, opts);
    ASSERT_TRUE(rr.Ok()) << rr.error;

    std::cout << "Fleet: " << gp.protein().AtomCount() << " atoms, "
              << gp.protein().ConformationCount() << " poses, "
              << rr.attached.size() << " results\n";
}


// ============================================================================
// Trajectory path: full-system XTC with water + ions
// ============================================================================

TEST(GromacsStreaming, TrajectoryBuildAndScan) {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);

    std::string tpr = FULLSYS_DIR + "/run_params/prod_fullsys.tpr";
    std::string xtc = FULLSYS_DIR + "/walker_0/md.xtc";

    if (!fs::exists(tpr) || !fs::exists(xtc)) {
        GTEST_SKIP() << "Full-system test data not found";
    }

    // --- 1. Build adapter from TPR ---
    nmr::GromacsProtein gp;
    ASSERT_TRUE(gp.BuildFromTrajectory(tpr)) << gp.error();

    // Protein not finalized yet — no bonds, no rings
    EXPECT_EQ(gp.protein().AtomCount(), 479u);
    EXPECT_EQ(gp.protein().BondCount(), 0u);

    // --- 2. Open trajectory (finalizes protein from first frame) ---
    nmr::GromacsFrameHandler handler(gp);
    ASSERT_TRUE(handler.Open(xtc, tpr)) << handler.error();

    // Protein finalized: bonds, rings, conformation 0
    EXPECT_GT(gp.protein().BondCount(), 0u);
    EXPECT_GT(gp.protein().RingCount(), 0u);
    EXPECT_EQ(gp.protein().ConformationCount(), 1u);
    EXPECT_EQ(gp.AtomCount(), 479u);

    std::cout << "Trajectory: " << gp.protein().AtomCount() << " atoms, "
              << gp.protein().BondCount() << " bonds, "
              << gp.protein().RingCount() << " rings\n";

    // --- 3. Scan 10 more frames (lightweight, no NPY) ---
    nmr::RunOptions scan_opts;
    scan_opts.charge_source = gp.charges();
    scan_opts.skip_mopac = true;
    scan_opts.skip_apbs = true;
    scan_opts.skip_coulomb = true;
    scan_opts.skip_dssp = true;

    size_t scanned = 0;
    for (int i = 0; i < 10; ++i) {
        if (!handler.Next(scan_opts)) break;
        ++scanned;
    }
    ASSERT_GT(scanned, 0u);

    std::cout << "Scanned " << scanned << " frames after frame 0\n";

    // --- 4. Check accumulators have data ---
    // Frame 0 + scanned frames should all be accumulated
    size_t expected_frames = 1 + scanned;
    EXPECT_EQ(gp.AtomAt(0).sasa.count, static_cast<int>(expected_frames));
    EXPECT_EQ(gp.AtomAt(0).water_n_first.count,
              static_cast<int>(expected_frames));

    // RMSF should be non-negative
    EXPECT_GE(gp.AtomAt(0).RMSF(), 0.0);

    // --- 5. Write catalog ---
    std::string temp_dir = "/tmp/nmr_trajectory_test_" +
                           std::to_string(getpid());
    fs::create_directories(temp_dir);

    std::string catalog = temp_dir + "/atom_catalog.csv";
    gp.WriteCatalog(catalog);
    ASSERT_TRUE(fs::exists(catalog));
    EXPECT_GT(fs::file_size(catalog), 1000u);

    std::cout << "Catalog: " << fs::file_size(catalog) << " bytes\n";

    // --- 6. Clean up ---
    std::string rm_cmd = "rm -rf " + temp_dir;
    (void)system(rm_cmd.c_str());
}


// ============================================================================
// Trajectory pass 2: reopen and extract selected frames with NPY
// ============================================================================

TEST(GromacsStreaming, TrajectoryExtractSelectedFrames) {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);

    std::string tpr = FULLSYS_DIR + "/run_params/prod_fullsys.tpr";
    std::string xtc = FULLSYS_DIR + "/walker_0/md.xtc";

    if (!fs::exists(tpr) || !fs::exists(xtc)) {
        GTEST_SKIP() << "Full-system test data not found";
    }

    // Build + open
    nmr::GromacsProtein gp;
    ASSERT_TRUE(gp.BuildFromTrajectory(tpr)) << gp.error();

    nmr::GromacsFrameHandler handler(gp);
    ASSERT_TRUE(handler.Open(xtc, tpr)) << handler.error();

    // Scan 20 frames
    nmr::RunOptions scan_opts;
    scan_opts.charge_source = gp.charges();
    scan_opts.skip_mopac = true;
    scan_opts.skip_apbs = true;
    scan_opts.skip_coulomb = true;
    scan_opts.skip_dssp = true;

    size_t total = 1;  // frame 0
    while (total < 20 && handler.Next(scan_opts)) ++total;

    std::cout << "Scanned " << total << " frames\n";

    // Reopen for pass 2
    ASSERT_TRUE(handler.Reopen()) << handler.error();

    // Extract frames 3 and 7 with full calculators
    std::set<size_t> selected = {3, 7};
    std::string temp_dir = "/tmp/nmr_extract_test_" +
                           std::to_string(getpid());
    fs::create_directories(temp_dir);

    nmr::RunOptions extract_opts;
    extract_opts.charge_source = gp.charges();
    extract_opts.skip_mopac = true;
    extract_opts.skip_coulomb = true;

    size_t extracted = 0;
    for (size_t fi = 1; fi < total; ++fi) {
        if (selected.count(fi)) {
            ASSERT_TRUE(handler.Next(extract_opts, temp_dir))
                << "failed to extract frame " << fi;
            ++extracted;
        } else {
            handler.Skip();
        }
    }

    EXPECT_EQ(extracted, 2u);

    // Verify NPY directories exist
    EXPECT_TRUE(fs::exists(temp_dir + "/frame_0003"));
    EXPECT_TRUE(fs::exists(temp_dir + "/frame_0007"));

    // Verify NPY files were written
    size_t npy_count = 0;
    for (auto& entry : fs::recursive_directory_iterator(temp_dir)) {
        if (entry.path().extension() == ".npy") ++npy_count;
    }
    EXPECT_GT(npy_count, 0u);

    std::cout << "Extracted " << extracted << " frames, "
              << npy_count << " NPY files\n";

    // Clean up
    std::string rm_cmd = "rm -rf " + temp_dir;
    (void)system(rm_cmd.c_str());
}

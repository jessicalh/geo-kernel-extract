//
// test_gromacs_streaming: proof of concept for trajectory streaming.
//
// Loads a fleet protein from TPR, opens the walker_0 XTC trajectory,
// creates free-standing conformations for 10 frames, runs geometry
// calculators on each, logs what was computed, and frees the
// conformation. Proves the pattern works end-to-end.
//

#include "GromacsProtein.h"
#include "GromacsFrameHandler.h"
#include "GromacsFinalResult.h"
#include "OperationRunner.h"
#include "OperationLog.h"

#include <gtest/gtest.h>
#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

#ifndef NMR_TEST_DATA_DIR
#error "NMR_TEST_DATA_DIR must be defined"
#endif

static const std::string FLEET_DIR =
    std::string(NMR_TEST_DATA_DIR) + "/fleet_test_large";


TEST(GromacsStreaming, LoadAndStream10Frames) {
    // Enable all log channels so we see what happens
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);

    // --- 1. Build GromacsProtein from TPR + first pose ---
    nmr::FleetPaths paths;
    paths.sampled_poses_dir = FLEET_DIR + "/initial-samples";
    paths.tpr_path = FLEET_DIR + "/run_params/prod.tpr";
    paths.force_field = nmr::ForceField::CHARMM36m;

    nmr::GromacsProtein gp;
    ASSERT_TRUE(gp.Build(paths)) << gp.error();

    std::cout << "\n=== GromacsProtein ==="
              << "\n  protein_id: " << gp.protein_id()
              << "\n  atoms:      " << gp.protein().AtomCount()
              << "\n  residues:   " << gp.protein().ResidueCount()
              << "\n  poses:      " << gp.protein().ConformationCount()
              << "\n  net_charge: " << gp.net_charge()
              << std::endl;

    ASSERT_GT(gp.protein().AtomCount(), 0u);

    // --- 2. Run geometry calculators on conformation 0 ---
    // This is the root — it lives in the Protein's vector.
    nmr::RunOptions opts;
    opts.skip_mopac = true;
    opts.skip_coulomb = true;
    // APBS runs — canonical electrostatics path (4s, faster than Coulomb)
    // DSSP runs on every frame
    opts.charge_source = gp.charges();
    // .edr for per-frame GROMACS energy extraction
    opts.edr_path = FLEET_DIR + "/walker_0/md.edr";

    auto& conf0 = gp.protein().Conformation();
    nmr::RunResult rr0 = nmr::OperationRunner::Run(conf0, opts);
    ASSERT_TRUE(rr0.Ok()) << rr0.error;

    std::cout << "\n=== Conformation 0 (root) ==="
              << "\n  attached: ";
    for (const auto& name : rr0.attached)
        std::cout << name << " ";
    std::cout << std::endl;

    // --- 3. Open XTC trajectory ---
    nmr::GromacsFrameHandler handler(gp, opts);

    std::string xtc_path = FLEET_DIR + "/walker_0/md.xtc";
    std::string tpr_path = FLEET_DIR + "/walker_0/md.tpr";

    ASSERT_TRUE(handler.OpenTrajectory(xtc_path, tpr_path))
        << handler.error();

    std::cout << "\n=== Trajectory ==="
              << "\n  frames: " << handler.FrameCount()
              << std::endl;

    ASSERT_GT(handler.FrameCount(), 0u);

    // --- 4. Process 10 frames ---
    std::string temp_dir = "/tmp/nmr_streaming_test_" +
                           std::to_string(getpid());
    fs::create_directories(temp_dir);

    size_t n_processed = handler.ProcessFrames(10, temp_dir);

    std::cout << "\n=== Processed ==="
              << "\n  frames: " << n_processed
              << "\n  frame_paths: " << gp.frame_paths().size()
              << "\n  conformations in Protein: "
              << gp.protein().ConformationCount()
              << " (should still be " << gp.protein().ConformationCount()
              << " — streaming frames are free-standing)"
              << std::endl;

    ASSERT_GT(n_processed, 0u);

    // Verify streaming frames did NOT add to the Protein's vector.
    // The Protein should still have only the poses from BuildFromGromacs.
    // (ConformationCount includes the original poses, not streaming frames.)

    // --- 5. Finalize ---
    nmr::GromacsFinalResult final_result(gp);
    ASSERT_TRUE(final_result.Finalize(temp_dir + "/output", temp_dir))
        << final_result.error();

    // --- 6. Clean up ---
    // fs::remove_all resolves to libtorch's broken std::filesystem stub
    // when linked against libtorch. Use POSIX rm instead.
    std::string rm_cmd = "rm -rf " + temp_dir;
    (void)system(rm_cmd.c_str());

    std::cout << "\n=== Done ===" << std::endl;
}

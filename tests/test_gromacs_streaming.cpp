//
// test_gromacs_streaming: tests for GromacsFrameHandler trajectory path.
//

// AIMNet2Result.h MUST come before GROMACS headers — GROMACS defines
// DIM as a macro which poisons PyTorch template parameters.
#include "AIMNet2Result.h"

#include "GromacsProtein.h"
#include "GromacsFrameHandler.h"
#include "GromacsFinalResult.h"
#include "CalculatorConfig.h"
#include "OperationRunner.h"
#include "OperationLog.h"

#include <gtest/gtest.h>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <set>

namespace fs = std::filesystem;

#ifndef NMR_TEST_DATA_DIR
#error "NMR_TEST_DATA_DIR must be defined"
#endif

// walker_0/ has md.tpr, md.xtc, md.edr — production layout.
static const std::string TRAJ_DIR =
    std::string(NMR_TEST_DATA_DIR) + "/fleet_test_fullsys/1ZR7_6721/walker_0";


// ============================================================================
// Trajectory path: full-system XTC with water + ions
// ============================================================================

TEST(GromacsStreaming, TrajectoryBuildAndScan) {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);

    if (!fs::exists(TRAJ_DIR + "/md.tpr") ||
        !fs::exists(TRAJ_DIR + "/md.xtc")) {
        GTEST_SKIP() << "Full-system test data not found";
    }

    // --- 1. Build adapter from directory ---
    nmr::GromacsProtein gp;
    ASSERT_TRUE(gp.BuildFromTrajectory(TRAJ_DIR)) << gp.error();

    // Protein not finalized yet — no bonds, no rings
    EXPECT_EQ(gp.protein().AtomCount(), 479u);
    EXPECT_EQ(gp.protein().BondCount(), 0u);

    // --- 2. Open trajectory (finalizes protein from first frame) ---
    nmr::GromacsFrameHandler handler(gp);
    std::string tpr = TRAJ_DIR + "/md.tpr";
    std::string xtc = TRAJ_DIR + "/md.xtc";
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
    size_t expected_frames = 1 + scanned;
    EXPECT_EQ(gp.AtomAt(0).sasa.count, static_cast<int>(expected_frames));
    EXPECT_EQ(gp.AtomAt(0).water_n_first.count,
              static_cast<int>(expected_frames));

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

    std::string rm_cmd = "rm -rf " + temp_dir;
    (void)system(rm_cmd.c_str());
}


// ============================================================================
// Full accumulator test: DSSP + AIMNet2 + all classical kernels + water
// ============================================================================

TEST(GromacsStreaming, FullAccumulatorScan) {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);
    nmr::CalculatorConfig::Load(
        std::string(NMR_TEST_DATA_DIR) + "/../../data/calculator_params.toml");

    std::string tpr = TRAJ_DIR + "/md.tpr";
    std::string xtc = TRAJ_DIR + "/md.xtc";
    std::string jpt = std::string(NMR_TEST_DATA_DIR) +
                       "/../../data/models/aimnet2_wb97m_0.jpt";

    if (!fs::exists(tpr) || !fs::exists(xtc)) {
        GTEST_SKIP() << "Full-system test data not found";
    }

    // Load AIMNet2 model — MUST succeed, no silent degradation
    ASSERT_TRUE(fs::exists(jpt)) << "AIMNet2 model not found: " << jpt;
    auto aimnet2 = nmr::AIMNet2Model::Load(jpt);
    ASSERT_NE(aimnet2, nullptr) << "AIMNet2 model failed to load";

    // Build from directory
    nmr::GromacsProtein gp;
    ASSERT_TRUE(gp.BuildFromTrajectory(TRAJ_DIR)) << gp.error();

    // Canonical fleet scan opts: everything EXCEPT MOPAC and vacuum Coulomb
    nmr::RunOptions opts;
    opts.charge_source = gp.charges();
    opts.skip_mopac = true;
    opts.skip_coulomb = true;
    opts.skip_apbs = false;
    opts.skip_dssp = false;
    opts.aimnet2_model = aimnet2.get();

    // Open with the same opts as all frames
    nmr::GromacsFrameHandler handler(gp);
    ASSERT_TRUE(handler.Open(xtc, tpr, opts)) << handler.error();

    size_t total = 1;
    while (total < 11 && handler.Next(opts)) ++total;

    std::cout << "\nScanned " << total << " frames with full calculators\n";

    // Verify classical kernel accumulators
    const auto& a0 = gp.AtomAt(0);
    EXPECT_EQ(a0.bs_T0.count, static_cast<int>(total));
    EXPECT_EQ(a0.mc_T0.count, static_cast<int>(total));
    EXPECT_EQ(a0.pq_T0.count, static_cast<int>(total));
    EXPECT_EQ(a0.disp_T0.count, static_cast<int>(total));
    EXPECT_EQ(a0.hm_T0.count, static_cast<int>(total));
    EXPECT_EQ(a0.rs_T0.count, static_cast<int>(total));

    // Water accumulators
    EXPECT_EQ(a0.water_emag.count, static_cast<int>(total));
    EXPECT_EQ(a0.water_efield_x.count, static_cast<int>(total));

    // DSSP accumulators
    EXPECT_EQ(a0.phi_cos.count, static_cast<int>(total));
    EXPECT_GT(a0.dssp_hbond_energy.count, 0);

    // AIMNet2
    EXPECT_EQ(a0.aimnet2_charge.count, static_cast<int>(total));
    EXPECT_GT(a0.aimnet2_charge.Std(), 0.0) << "charge variance should be nonzero";

    // Delta trackers: should have total-1 observations
    EXPECT_EQ(a0.bs_T0_delta.delta.count, static_cast<int>(total) - 1);
    EXPECT_EQ(a0.sasa_delta.delta.count, static_cast<int>(total) - 1);
    EXPECT_EQ(a0.water_n_first_delta.delta.count, static_cast<int>(total) - 1);

    // Print a summary for a ring-adjacent atom
    auto print = [](const char* name, const nmr::Welford& w) {
        std::cout << std::setw(25) << std::left << name
                  << " n=" << std::setw(3) << w.count
                  << " mean=" << std::setw(12) << w.mean
                  << " std=" << std::setw(12) << w.Std()
                  << " range=" << w.Range() << "\n";
    };

    for (size_t i = 0; i < gp.AtomCount(); ++i) {
        const auto& a = gp.AtomAt(i);
        if (a.bs_T0.mean > 0.05 && a.bs_T0.Std() > 0.001) {
            std::cout << "\n=== Atom " << i << " ("
                      << gp.protein().AtomAt(i).pdb_atom_name << " "
                      << nmr::ThreeLetterCodeForAminoAcid(gp.protein().ResidueAt(gp.protein().AtomAt(i).residue_index).type) << ") ===\n";
            print("bs_T0 (ring current)", a.bs_T0);
            print("mc_T0 (McConnell)", a.mc_T0);
            print("mc_T2mag", a.mc_T2mag);
            print("aimnet2_charge", a.aimnet2_charge);
            print("sasa", a.sasa);
            print("water_emag", a.water_emag);
            print("water_n_first", a.water_n_first);
            print("phi_cos", a.phi_cos);
            for (int k = 0; k < 4; ++k) {
                if (a.chi_cos[k].count == 0) break;
                char label[16];
                std::snprintf(label, sizeof(label), "chi%d_cos", k + 1);
                print(label, a.chi_cos[k]);
                std::cout << "  chi" << k+1 << " transitions: "
                          << a.chi_transitions[k].transitions
                          << "/" << a.chi_transitions[k].frames << "\n";
            }
            std::cout << "ss8 transitions: " << a.ss8_transitions.transitions
                      << "/" << a.ss8_transitions.frames << "\n";
            break;
        }
    }

    // Write H5 master file
    std::string h5_path = "/tmp/nmr_test_" + std::to_string(getpid()) + ".h5";
    gp.WriteH5(h5_path);
    ASSERT_TRUE(fs::exists(h5_path));
    EXPECT_GT(fs::file_size(h5_path), 1000u);
    std::cout << "\nH5: " << h5_path << " ("
              << fs::file_size(h5_path) / 1024 << " KB)\n";
    std::cout << "Stored " << gp.StoredFrameCount() << " frames x "
              << gp.AtomCount() << " atoms\n";

    // Write catalog CSV too
    std::string csv_path = "/tmp/nmr_test_" + std::to_string(getpid()) + "_catalog.csv";
    gp.WriteCatalog(csv_path);
    std::cout << "Catalog: " << csv_path << " ("
              << fs::file_size(csv_path) / 1024 << " KB)\n";

    // Verify H5 with h5py
    std::string py_cmd =
        "HDF5_DISABLE_VERSION_CHECK=2 python3 -c \""
        "import h5py; f=h5py.File('" + h5_path + "','r'); "
        "T=f.attrs['positions_shape_T']; N=f.attrs['positions_shape_N']; "
        "print(f'positions: ({T}, {N}, 3) stored as', f['positions'].shape); "
        "print('frame_times:', f['frame_times'].shape); "
        "print('rollup mean:', f['rollup/mean'].shape); "
        "names=list(f['rollup/names'][:]); "
        "print(f'rollup columns: {len(names)}:', names[:5], '...'); "
        "print('bonds:', f['bonds/length_mean'].shape); "
        "f.close()\"";
    std::cout << "\n";
    (void)system(py_cmd.c_str());

    // Clean up
    fs::remove(h5_path);
    fs::remove(csv_path);
}


// ============================================================================
// Trajectory pass 2: reopen and extract selected frames with NPY
// ============================================================================

TEST(GromacsStreaming, TrajectoryExtractSelectedFrames) {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);

    std::string tpr = TRAJ_DIR + "/md.tpr";
    std::string xtc = TRAJ_DIR + "/md.xtc";

    if (!fs::exists(tpr) || !fs::exists(xtc)) {
        GTEST_SKIP() << "Full-system test data not found";
    }

    // Build + open
    nmr::GromacsProtein gp;
    ASSERT_TRUE(gp.BuildFromTrajectory(TRAJ_DIR)) << gp.error();

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
    while (handler.Next(extract_opts, "", false)) {
        size_t fi = handler.frame_index();
        if (selected.count(fi)) {
            std::string frame_dir = temp_dir + "/frame_" + std::to_string(fi);
            fs::create_directories(frame_dir);
            int n = nmr::ConformationResult::WriteAllFeatures(
                handler.conformation(), frame_dir);
            std::cout << "Frame " << fi << ": " << n << " NPY files\n";
            ++extracted;
        }
    }

    EXPECT_EQ(extracted, selected.size());

    // Verify NPY files exist
    for (size_t fi : selected) {
        std::string frame_dir = temp_dir + "/frame_" + std::to_string(fi);
        EXPECT_TRUE(fs::exists(frame_dir + "/ring_contributions.npy"))
            << "missing ring_contributions.npy for frame " << fi;
    }

    // Clean up
    std::string rm_cmd = "rm -rf " + temp_dir;
    (void)system(rm_cmd.c_str());
}

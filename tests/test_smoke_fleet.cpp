//
// test_smoke_fleet.cpp
//
// End-to-end smoke test for the GROMACS fleet extraction pipeline.
//
// Loads 1A6J_5789 (2480 atoms, 10 poses) via BuildFromGromacs, runs the
// full extraction pipeline on every pose via RunAllFrames, writes NPY
// features per frame, validates logs and output.
//
// This is expensive (~40 min) because every pose runs MOPAC + APBS.
// It exercises the exact path that the 700-protein x 10-pose production
// extraction uses. Run it:
//   - After changing GromacsEnsembleLoader or TPR reading
//   - After changing OperationRunner::RunEnsemble
//   - After changing any extractor (if the quick smoke_tests pass first)
//   - Before launching a production extraction batch
//
// Separate executable: ./fleet_smoke_tests
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
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <set>
#include <algorithm>
#include <cstring>

namespace fs = std::filesystem;
using namespace nmr;


// ============================================================================
// Helpers (same as test_smoke.cpp — duplicated to keep executables independent)
// ============================================================================

static std::string TimestampDir() {
    auto now = std::chrono::system_clock::now();
    auto time = std::chrono::system_clock::to_time_t(now);
    std::tm tm_buf;
    localtime_r(&time, &tm_buf);
    std::ostringstream oss;
    oss << std::put_time(&tm_buf, "%Y-%m-%d_%H%M%S");
    return oss.str();
}

struct LogEntry {
    std::string level;
    std::string op;
    std::string detail;
};

static std::vector<LogEntry> ParseLog(const std::string& path) {
    std::vector<LogEntry> entries;
    std::ifstream in(path);
    std::string line;
    while (std::getline(in, line)) {
        LogEntry e;
        auto extract = [&](const std::string& key) -> std::string {
            std::string needle = "\"" + key + "\":\"";
            auto pos = line.find(needle);
            if (pos == std::string::npos) return "";
            pos += needle.size();
            auto end = line.find('"', pos);
            if (end == std::string::npos) return "";
            return line.substr(pos, end - pos);
        };
        e.level = extract("level");
        e.op = extract("op");
        e.detail = extract("detail");
        entries.push_back(std::move(e));
    }
    return entries;
}

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

static bool FilesIdentical(const std::string& a, const std::string& b) {
    std::ifstream fa(a, std::ios::binary | std::ios::ate);
    std::ifstream fb(b, std::ios::binary | std::ios::ate);
    if (!fa.is_open() || !fb.is_open()) return false;
    if (fa.tellg() != fb.tellg()) return false;
    fa.seekg(0); fb.seekg(0);
    const size_t BUF = 8192;
    char ba[BUF], bb[BUF];
    while (fa && fb) {
        fa.read(ba, BUF);
        fb.read(bb, BUF);
        auto n = fa.gcount();
        if (n != fb.gcount()) return false;
        if (std::memcmp(ba, bb, static_cast<size_t>(n)) != 0) return false;
    }
    return true;
}


// ============================================================================
// Fleet smoke test
// ============================================================================

TEST(SmokeFleet, AllPoses) {
    RuntimeEnvironment::Load();

    // Check fleet test data is available
    std::string fleet_base = test::TestEnvironment::FleetData();
    if (fleet_base.empty() || !fs::exists(fleet_base)) {
        GTEST_SKIP() << "Fleet test data not found";
    }

    FleetPaths paths;
    paths.sampled_poses_dir = fleet_base + "/1A6J_5789/poses";
    paths.tpr_path = fleet_base + "/1A6J_5789/params/prod.tpr";
    paths.force_field = ForceField::CHARMM36m;

    if (!fs::exists(paths.tpr_path)) {
        GTEST_SKIP() << "1A6J_5789 TPR not found at " << paths.tpr_path;
    }

    // ---- Load ----
    auto build = BuildFromGromacs(paths);
    ASSERT_TRUE(build.Ok()) << build.error;
    ASSERT_EQ(build.protein->MDFrameCount(), 10u)
        << "Expected 10 poses, got " << build.protein->MDFrameCount();

    std::cout << "\n  Fleet: 1A6J_5789\n"
              << "  Atoms: " << build.protein->AtomCount() << "\n"
              << "  Residues: " << build.protein->ResidueCount() << "\n"
              << "  Rings: " << build.protein->RingCount() << "\n"
              << "  Poses: " << build.protein->MDFrameCount() << "\n";

    // ---- Output directory ----
    std::string data_dir;
#ifdef NMR_TEST_DATA_DIR
    data_dir = NMR_TEST_DATA_DIR;
#endif
    std::string run_dir = data_dir + "/../golden/smoke/" + TimestampDir() + "/fleet";
    fs::create_directories(run_dir);

    std::string log_path = run_dir + "/log.jsonl";
    OperationLog::ConfigureFile(log_path);

    // ---- Run pipeline on all poses ----
    RunOptions opts;
    opts.charge_source = build.charges.get();
    opts.net_charge = build.net_charge;

    auto run_results = RunAllFrames(*build.protein, opts);
    ASSERT_EQ(run_results.size(), 10u);

    // Every frame should succeed with 13+ results
    for (size_t i = 0; i < run_results.size(); ++i) {
        EXPECT_TRUE(run_results[i].Ok())
            << "Frame " << i << " error: " << run_results[i].error;
        EXPECT_GE(run_results[i].attached.size(), 13u)
            << "Frame " << i << ": expected 13+ results, got "
            << run_results[i].attached.size();
    }

    // Log what frame 0 attached (all frames should be identical)
    std::cout << "  Results per frame (" << run_results[0].attached.size() << "):\n";
    for (const auto& name : run_results[0].attached) {
        std::cout << "    " << name << "\n";
    }

    // ---- Write features per frame ----
    int total_arrays = 0;
    size_t min_per_frame = 35;

    for (size_t i = 0; i < build.protein->MDFrameCount(); ++i) {
        auto& frame = build.protein->MDFrameAt(i);
        char frame_dir[512];
        std::snprintf(frame_dir, sizeof(frame_dir), "%s/frame_%03zu",
                      run_dir.c_str(), i + 1);
        fs::create_directories(frame_dir);

        int arrays = ConformationResult::WriteAllFeatures(frame, frame_dir);
        EXPECT_GE(arrays, static_cast<int>(min_per_frame))
            << "Frame " << i << ": too few arrays";
        total_arrays += arrays;

        // Validate NPY files in this frame
        auto files = NpyFiles(frame_dir);
        for (const auto& name : files) {
            std::string path = std::string(frame_dir) + "/" + name;
            EXPECT_GT(fs::file_size(path), 0u)
                << "Frame " << i << " " << name << " is empty";
            EXPECT_TRUE(VerifyNpyMagic(path))
                << "Frame " << i << " " << name << " bad NPY magic";
        }
    }

    OperationLog::CloseFile();

    // ---- Validate log ----
    auto entries = ParseLog(log_path);
    ASSERT_GT(entries.size(), 0u) << "Log file is empty";

    int error_count = 0;
    for (const auto& e : entries) {
        if (e.level == "ERROR") {
            error_count++;
            ADD_FAILURE() << "Log ERROR: [" << e.op << "] " << e.detail;
        }
    }

    // Scope matching
    std::vector<std::string> open_scopes;
    for (const auto& e : entries) {
        if (e.op.find("[BEGIN]") != std::string::npos) {
            std::string base = e.op.substr(0, e.op.find(" [BEGIN]"));
            open_scopes.push_back(base);
        } else if (e.op.find("[END]") != std::string::npos) {
            std::string base = e.op.substr(0, e.op.find(" [END]"));
            auto it = std::find(open_scopes.rbegin(), open_scopes.rend(), base);
            if (it != open_scopes.rend()) {
                open_scopes.erase((it + 1).base());
            }
        }
    }
    EXPECT_EQ(open_scopes.size(), 0u) << "Unclosed log scopes";

    // ---- Binary comparison to blessed baseline (frame 1 only) ----
    std::string blessed_frame1;
#ifdef NMR_TEST_DATA_DIR
    blessed_frame1 = std::string(NMR_TEST_DATA_DIR)
        + "/../golden/blessed/fleet/frame_001";
#endif

    if (!blessed_frame1.empty() && fs::exists(blessed_frame1)) {
        std::string run_frame1 = run_dir + "/frame_001";
        auto run_files = NpyFiles(run_frame1);
        auto blessed_files = NpyFiles(blessed_frame1);

        int identical = 0, different = 0;
        for (const auto& f : blessed_files) {
            if (!run_files.count(f)) {
                ADD_FAILURE() << "Missing from run: " << f;
                continue;
            }
            if (FilesIdentical(run_frame1 + "/" + f, blessed_frame1 + "/" + f)) {
                identical++;
            } else {
                different++;
                ADD_FAILURE() << "BINARY DIFF (frame_001): " << f;
            }
        }
        for (const auto& f : run_files) {
            if (!blessed_files.count(f))
                std::cout << "  NEW file: " << f << "\n";
        }
        std::cout << "  Binary comparison (frame_001): "
                  << identical << " identical, " << different << " different\n";
    } else {
        std::cout << "  No blessed fleet baseline — skipping binary comparison.\n"
                  << "  To bless:\n"
                  << "    mkdir -p " << (blessed_frame1.empty() ? "<blessed>/fleet" : blessed_frame1.substr(0, blessed_frame1.rfind('/'))) << "\n"
                  << "    cp -r " << run_dir << "/frame_001 " << blessed_frame1 << "\n";
    }

    // ---- Summary ----
    std::cout << "\n  Fleet smoke complete:\n"
              << "  Poses: " << build.protein->MDFrameCount() << "\n"
              << "  Total NPY arrays: " << total_arrays << "\n"
              << "  Log: " << entries.size() << " entries, "
              << error_count << " errors\n"
              << "  Output: " << run_dir << "\n";
}

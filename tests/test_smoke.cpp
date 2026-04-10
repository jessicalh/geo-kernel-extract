//
// test_smoke.cpp
//
// End-to-end smoke test for the extraction pipeline.
//
// Two configurations exercise both pipeline branches:
//   SmokeNoDft  — 1UBQ + ff14SB charges, no ORCA path.
//   SmokeWithDft — P84477 from consolidated, with ORCA NMR output.
//
// Each run:
//   1. Loads protein and runs the pipeline (all extractors, in order)
//   2. Calls WriteAllFeatures to a timestamped output directory
//   3. Captures the operation log to log.jsonl in the same directory
//   4. Validates the log: no errors, expected calculators present
//   5. Validates the NPY output: file count, non-zero sizes, header sanity
//   6. Binary-compares to blessed baseline (if available)
//
// This is the test you run when you change an extractor, add a parameter,
// or touch the pipeline. If it passes, the system still produces the same
// bytes it did before your change.
//

#include "TestEnvironment.h"
#include <gtest/gtest.h>

#include "OperationRunner.h"
#include "ConformationResult.h"
#include "PdbFileReader.h"
#include "OrcaRunLoader.h"
#include "ChargeSource.h"
#include "RuntimeEnvironment.h"
#include "OperationLog.h"

#include <filesystem>
#include <fstream>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <set>

namespace fs = std::filesystem;
using namespace nmr;


// ============================================================================
// Helpers
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


// Parse a log.jsonl file and extract levels and operation names.
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
        // Minimal JSON extraction — fields are always in fixed order from BuildJson.
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


// Collect .npy filenames from a directory.
static std::set<std::string> NpyFiles(const std::string& dir) {
    std::set<std::string> result;
    if (!fs::exists(dir)) return result;
    for (const auto& entry : fs::directory_iterator(dir)) {
        if (entry.path().extension() == ".npy")
            result.insert(entry.path().filename().string());
    }
    return result;
}


// Read first 10 bytes of a .npy file and verify magic.
static bool VerifyNpyMagic(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    char buf[6];
    if (!in.read(buf, 6)) return false;
    return buf[0] == '\x93' && buf[1] == 'N' && buf[2] == 'U'
        && buf[3] == 'M' && buf[4] == 'P' && buf[5] == 'Y';
}


// Binary-compare two files. Returns true if identical.
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
// Smoke test fixture
// ============================================================================

class SmokeTest : public ::testing::Test {
protected:
    // The root for all smoke output: tests/golden/smoke/
    std::string smoke_root_;
    // This run's timestamped directory
    std::string run_dir_;

    void SetUp() override {
        RuntimeEnvironment::Load();

        // Build output paths
        std::string data_dir;
#ifdef NMR_TEST_DATA_DIR
        data_dir = NMR_TEST_DATA_DIR;
#endif
        smoke_root_ = data_dir + "/../golden/smoke";
        run_dir_ = smoke_root_ + "/" + TimestampDir();
    }

    // Run the full pipeline, write output, validate everything.
    // Returns the output directory path.
    std::string RunSmoke(
        const std::string& label,
        ProteinConformation& conf,
        const RunOptions& opts,
        size_t min_results,
        size_t min_npy_files,
        const std::string& blessed_dir)
    {
        std::string out_dir = run_dir_ + "/" + label;
        fs::create_directories(out_dir);

        std::string log_path = out_dir + "/log.jsonl";
        OperationLog::ConfigureFile(log_path);

        // ---- Phase 1: Run the pipeline ----
        auto result = OperationRunner::Run(conf, opts);
        EXPECT_TRUE(result.Ok()) << "Pipeline error: " << result.error;
        EXPECT_GE(result.attached.size(), min_results)
            << "Expected " << min_results << "+ results, got "
            << result.attached.size();

        // Log what attached
        std::cout << "\n  [" << label << "] Pipeline results ("
                  << result.attached.size() << "):\n";
        for (const auto& name : result.attached) {
            std::cout << "    " << name << "\n";
        }

        // ---- Phase 2: Write all features ----
        int arrays = ConformationResult::WriteAllFeatures(conf, out_dir);
        EXPECT_GE(arrays, static_cast<int>(min_npy_files))
            << "Expected " << min_npy_files << "+ NPY arrays";

        OperationLog::CloseFile();

        // ---- Phase 3: Validate log ----
        ValidateLog(log_path, result.attached);

        // ---- Phase 4: Validate NPY files ----
        ValidateNpy(out_dir, min_npy_files);

        // ---- Phase 5: Binary comparison to blessed baseline ----
        if (!blessed_dir.empty() && fs::exists(blessed_dir)) {
            BinaryCompare(out_dir, blessed_dir);
        } else {
            std::cout << "  [" << label << "] No blessed baseline at "
                      << blessed_dir << " — skipping binary comparison.\n"
                      << "  To bless this run:\n"
                      << "    cp -r " << out_dir << " " << blessed_dir << "\n";
        }

        std::cout << "  [" << label << "] Output: " << out_dir << "\n";
        return out_dir;
    }

private:
    void ValidateLog(const std::string& log_path,
                     const std::vector<std::string>& expected_results) {
        auto entries = ParseLog(log_path);
        ASSERT_GT(entries.size(), 0u) << "Log file is empty";

        // No errors
        int error_count = 0;
        for (const auto& e : entries) {
            if (e.level == "ERROR") {
                error_count++;
                ADD_FAILURE() << "Log ERROR: [" << e.op << "] " << e.detail;
            }
        }
        std::cout << "  Log: " << entries.size() << " entries, "
                  << error_count << " errors\n";

        // Every expected result name appears in the log
        for (const auto& name : expected_results) {
            bool found = false;
            for (const auto& e : entries) {
                if (e.op.find(name) != std::string::npos ||
                    e.detail.find(name) != std::string::npos) {
                    found = true;
                    break;
                }
            }
            EXPECT_TRUE(found) << name << " not mentioned in log";
        }

        // Scope matching: every [BEGIN] has a [END]
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
        EXPECT_EQ(open_scopes.size(), 0u)
            << "Unclosed log scopes: " << open_scopes.size();
        for (const auto& s : open_scopes) {
            std::cout << "  UNCLOSED SCOPE: " << s << "\n";
        }
    }

    void ValidateNpy(const std::string& dir, size_t min_count) {
        auto files = NpyFiles(dir);
        EXPECT_GE(files.size(), min_count)
            << "Expected " << min_count << "+ .npy files, got " << files.size();

        int valid = 0;
        for (const auto& name : files) {
            std::string path = dir + "/" + name;
            auto sz = fs::file_size(path);
            EXPECT_GT(sz, 0u) << name << " is empty";
            EXPECT_TRUE(VerifyNpyMagic(path)) << name << " has bad NPY magic";
            if (sz > 0) valid++;
        }
        std::cout << "  NPY: " << files.size() << " files, "
                  << valid << " valid\n";
    }

    void BinaryCompare(const std::string& run_dir,
                       const std::string& blessed_dir) {
        auto run_files = NpyFiles(run_dir);
        auto blessed_files = NpyFiles(blessed_dir);

        // Check for missing files in either direction
        for (const auto& f : blessed_files) {
            EXPECT_TRUE(run_files.count(f))
                << "Missing from run: " << f << " (present in blessed)";
        }
        for (const auto& f : run_files) {
            if (!blessed_files.count(f)) {
                std::cout << "  NEW file (not in blessed): " << f << "\n";
            }
        }

        // Binary compare common files
        int identical = 0, different = 0;
        for (const auto& f : blessed_files) {
            if (!run_files.count(f)) continue;
            std::string a = run_dir + "/" + f;
            std::string b = blessed_dir + "/" + f;
            if (FilesIdentical(a, b)) {
                identical++;
            } else {
                different++;
                ADD_FAILURE() << "BINARY DIFF: " << f
                    << " (run=" << fs::file_size(a)
                    << " blessed=" << fs::file_size(b) << " bytes)";
            }
        }
        std::cout << "  Binary comparison: " << identical << " identical, "
                  << different << " different\n";
    }
};


// ============================================================================
// SmokeNoDft: 1UBQ + ff14SB file charges, no ORCA
// ============================================================================

TEST_F(SmokeTest, NoDft) {
    // Load 1UBQ
    if (test::TestEnvironment::UbqProtonated().empty() ||
        !fs::exists(test::TestEnvironment::UbqProtonated())) {
        GTEST_SKIP() << "1UBQ protonated PDB not available";
    }

    auto build = BuildFromProtonatedPdb(test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(build.Ok()) << build.error;
    auto& conf = build.protein->Conformation();

    // Set up charges from ff14SB parameter file
    RunOptions opts;
    std::unique_ptr<ParamFileChargeSource> charges;
    if (!test::TestEnvironment::Ff14sbParams().empty() &&
        fs::exists(test::TestEnvironment::Ff14sbParams())) {
        charges = std::make_unique<ParamFileChargeSource>(
            test::TestEnvironment::Ff14sbParams());
        opts.charge_source = charges.get();
    }

    // Blessed baseline
    std::string blessed;
#ifdef NMR_TEST_DATA_DIR
    blessed = std::string(NMR_TEST_DATA_DIR) + "/../golden/blessed/nodft";
#endif

    // Expect: foundation (4) + MOPAC + APBS + 6 ring calculators
    //         + Coulomb + MopacCoulomb + MopacMcConnell + HBond = ~17
    // NPY: 40+ arrays (no orca_*.npy)
    RunSmoke("nodft", conf, opts, 13, 35, blessed);
}


// ============================================================================
// SmokeWithDft: P84477 + prmtop charges + ORCA NMR
// ============================================================================

TEST_F(SmokeTest, WithDft) {
    std::string dir = std::string(test::TestEnvironment::Consolidated()) + "P84477/";
    if (!fs::exists(dir)) GTEST_SKIP() << "P84477 consolidated data not found";

    // Load protein from ORCA run
    OrcaRunFiles files;
    files.pdb_path = dir + "P84477_WT.pdb";
    files.xyz_path = dir + "P84477_WT.xyz";
    files.prmtop_path = dir + "P84477_WT.prmtop";

    auto build = BuildFromOrca(files);
    ASSERT_TRUE(build.Ok()) << build.error;
    auto& conf = build.protein->Conformation();

    RunOptions opts;
    opts.charge_source = build.charges.get();
    opts.net_charge = build.net_charge;

    // Find ORCA NMR output
    for (const auto& entry : fs::directory_iterator(dir)) {
        std::string name = entry.path().filename().string();
        if (name.find("P84477_WT") == 0 && name.find("_nmr.out") != std::string::npos) {
            opts.orca_nmr_path = entry.path().string();
            break;
        }
    }

    // Blessed baseline
    std::string blessed;
#ifdef NMR_TEST_DATA_DIR
    blessed = std::string(NMR_TEST_DATA_DIR) + "/../golden/blessed/withdft";
#endif

    // Expect: everything NoDft has + OrcaShieldingResult = ~18
    // NPY: 45+ arrays (including orca_*.npy)
    RunSmoke("withdft", conf, opts, 14, 40, blessed);
}

//
// test_job_spec: end-to-end tests for JobSpec parsing, validation,
// and execution through all five modes.
//
// Separate executable from nmr_tests so we can run these (~10 min
// with MOPAC) without the full 2.5-hour suite.
//
// Build:  cmake --build build --target job_spec_tests
// Run:    build/job_spec_tests
//

#include <gtest/gtest.h>

#include "JobSpec.h"
#include "BuildResult.h"
#include "PdbFileReader.h"
#include "OrcaRunLoader.h"
#include "GromacsEnsembleLoader.h"
#include "OperationRunner.h"
#include "ConformationResult.h"
#include "RuntimeEnvironment.h"
#include "CalculatorConfig.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"

#include <filesystem>
#include <cstring>
#include <vector>
#include <string>

namespace fs = std::filesystem;
using namespace nmr;

// ============================================================================
// Test data paths — resolved from NMR_TEST_DATA_DIR at compile time
// ============================================================================

#ifndef NMR_TEST_DATA_DIR
#error "NMR_TEST_DATA_DIR must be defined"
#endif

static const std::string DATA = NMR_TEST_DATA_DIR;

// Test data inventory:
//   tests/data/external/1UBQ.pdb                    — bare PDB (no H)
//   tests/data/1ubq_protonated.pdb                  — pre-protonated PDB
//   tests/data/orca/A0A7C5FAR6_WT.xyz/.prmtop/_nmr.out — ORCA WT
//   tests/data/orca/A0A7C5FAR6_ALA.xyz/.prmtop/_nmr.out — ORCA ALA
//   tests/data/fleet/1A6J_5789/params/prod.tpr      — GROMACS TPR
//   tests/data/fleet/1A6J_5789/poses/                — poses + ensemble.json


// ============================================================================
// Helper: simulate argv from a vector of strings
// ============================================================================

struct FakeArgv {
    std::vector<std::string> storage;
    std::vector<char*> ptrs;

    FakeArgv(std::initializer_list<std::string> args)
        : storage(args) {
        for (auto& s : storage)
            ptrs.push_back(s.data());
    }

    int argc() const { return static_cast<int>(ptrs.size()); }
    char** argv() { return ptrs.data(); }
};


// ============================================================================
// Parsing tests — no file I/O, just argv → JobSpec
// ============================================================================

TEST(JobSpecParse, PdbMode) {
    FakeArgv a{"nmr_extract", "--pdb", "/tmp/test.pdb", "--pH", "6.5",
               "--output", "/tmp/out", "--config", "/tmp/params.toml"};
    auto spec = ParseJobSpec(a.argc(), a.argv());

    EXPECT_TRUE(spec.Ok());
    EXPECT_EQ(spec.mode, JobMode::Pdb);
    EXPECT_EQ(spec.pdb_path, "/tmp/test.pdb");
    EXPECT_DOUBLE_EQ(spec.pH, 6.5);
    EXPECT_EQ(spec.output_dir, "/tmp/out");
    EXPECT_EQ(spec.config_path, "/tmp/params.toml");
}

TEST(JobSpecParse, ProtonatedPdbMode) {
    FakeArgv a{"nmr_extract", "--protonated-pdb", "/tmp/prot.pdb",
               "--output", "/tmp/out"};
    auto spec = ParseJobSpec(a.argc(), a.argv());

    EXPECT_TRUE(spec.Ok());
    EXPECT_EQ(spec.mode, JobMode::ProtonatedPdb);
    EXPECT_EQ(spec.pdb_path, "/tmp/prot.pdb");
}

TEST(JobSpecParse, OrcaMode) {
    FakeArgv a{"nmr_extract", "--orca", "--root", "/data/A0A7C5FAR6_WT",
               "--output", "/tmp/out"};
    auto spec = ParseJobSpec(a.argc(), a.argv());

    EXPECT_TRUE(spec.Ok());
    EXPECT_EQ(spec.mode, JobMode::Orca);
    EXPECT_EQ(spec.orca_files.xyz_path, "/data/A0A7C5FAR6_WT.xyz");
    EXPECT_EQ(spec.orca_files.prmtop_path, "/data/A0A7C5FAR6_WT.prmtop");
    EXPECT_EQ(spec.orca_files.nmr_out_path, "/data/A0A7C5FAR6_WT_nmr.out");
}

TEST(JobSpecParse, MutantMode) {
    FakeArgv a{"nmr_extract", "--mutant",
               "--wt", "/data/FOO_WT", "--ala", "/data/FOO_ALA",
               "--output", "/tmp/out"};
    auto spec = ParseJobSpec(a.argc(), a.argv());

    EXPECT_TRUE(spec.Ok());
    EXPECT_EQ(spec.mode, JobMode::Mutant);
    EXPECT_EQ(spec.wt_files.xyz_path, "/data/FOO_WT.xyz");
    EXPECT_EQ(spec.ala_files.prmtop_path, "/data/FOO_ALA.prmtop");
}

TEST(JobSpecParse, FleetModeReturnsError) {
    // --fleet was removed 2026-04-12. Passing it should give a clear error
    // pointing to --trajectory as the replacement.
    FakeArgv a{"nmr_extract", "--fleet",
               "--tpr", "/data/prod.tpr", "--poses", "/data/poses",
               "--output", "/tmp/out"};
    auto spec = ParseJobSpec(a.argc(), a.argv());

    EXPECT_FALSE(spec.Ok());
    EXPECT_NE(spec.error.find("removed"), std::string::npos);
    EXPECT_NE(spec.error.find("--trajectory"), std::string::npos);
}

TEST(JobSpecParse, MissingModeFails) {
    FakeArgv a{"nmr_extract", "--output", "/tmp/out"};
    auto spec = ParseJobSpec(a.argc(), a.argv());
    EXPECT_FALSE(spec.Ok());
}

TEST(JobSpecParse, OrcaMissingRootFails) {
    FakeArgv a{"nmr_extract", "--orca", "--output", "/tmp/out"};
    auto spec = ParseJobSpec(a.argc(), a.argv());
    EXPECT_FALSE(spec.Ok());
    EXPECT_NE(spec.error.find("--root"), std::string::npos);
}

TEST(JobSpecParse, MutantMissingAlaFails) {
    FakeArgv a{"nmr_extract", "--mutant", "--wt", "/data/FOO_WT",
               "--output", "/tmp/out"};
    auto spec = ParseJobSpec(a.argc(), a.argv());
    EXPECT_FALSE(spec.Ok());
}

TEST(JobSpecParse, HelpReturnsNone) {
    FakeArgv a{"nmr_extract", "--help"};
    auto spec = ParseJobSpec(a.argc(), a.argv());
    EXPECT_EQ(spec.mode, JobMode::None);
}

TEST(JobSpecParse, NoOutputIsOk) {
    // Viewer mode: no --output is valid
    FakeArgv a{"nmr_extract", "--pdb", "/tmp/test.pdb"};
    auto spec = ParseJobSpec(a.argc(), a.argv());
    EXPECT_TRUE(spec.Ok());
    EXPECT_TRUE(spec.output_dir.empty());
}


// ============================================================================
// Validation tests — check file existence against real test data
// ============================================================================

TEST(JobSpecValidate, PdbValid) {
    FakeArgv a{"test", "--pdb", (DATA + "/external/1UBQ.pdb").c_str(),
               "--output", "/tmp/jobspec_test_pdb"};
    auto spec = ParseJobSpec(a.argc(), a.argv());
    ASSERT_TRUE(spec.Ok());
    EXPECT_TRUE(ValidateJobSpec(spec));
    EXPECT_TRUE(spec.error.empty());
}

TEST(JobSpecValidate, PdbMissingFile) {
    FakeArgv a{"test", "--pdb", "/nonexistent/file.pdb",
               "--output", "/tmp/out"};
    auto spec = ParseJobSpec(a.argc(), a.argv());
    ASSERT_TRUE(spec.Ok());
    EXPECT_FALSE(ValidateJobSpec(spec));
    EXPECT_NE(spec.error.find("not found"), std::string::npos);
}

TEST(JobSpecValidate, OrcaValid) {
    std::string root = DATA + "/orca/A0A7C5FAR6_WT";
    FakeArgv a{"test", "--orca", "--root", root.c_str(),
               "--output", "/tmp/jobspec_test_orca"};
    auto spec = ParseJobSpec(a.argc(), a.argv());
    ASSERT_TRUE(spec.Ok());
    EXPECT_TRUE(ValidateJobSpec(spec));
    // NMR .out should be found
    EXPECT_FALSE(spec.orca_files.nmr_out_path.empty());
}

TEST(JobSpecValidate, OrcaMissingPrmtopReportsError) {
    // Use xyz root that exists but with a mangled prmtop name
    FakeArgv a{"test", "--orca", "--root", "/tmp/nonexistent_root",
               "--output", "/tmp/out"};
    auto spec = ParseJobSpec(a.argc(), a.argv());
    ASSERT_TRUE(spec.Ok());
    EXPECT_FALSE(ValidateJobSpec(spec));
    EXPECT_NE(spec.error.find("not found"), std::string::npos);
}

TEST(JobSpecValidate, MutantValid) {
    std::string wt  = DATA + "/orca/A0A7C5FAR6_WT";
    std::string ala = DATA + "/orca/A0A7C5FAR6_ALA";
    FakeArgv a{"test", "--mutant", "--wt", wt.c_str(), "--ala", ala.c_str(),
               "--output", "/tmp/jobspec_test_mutant"};
    auto spec = ParseJobSpec(a.argc(), a.argv());
    ASSERT_TRUE(spec.Ok());
    EXPECT_TRUE(ValidateJobSpec(spec));
}

// JobSpecValidate::FleetValid removed 2026-04-12 — --fleet mode removed.
// BuildFromGromacs library path is exercised by FleetLibraryDirect below.


// ============================================================================
// End-to-end tests: parse → validate → build → run → write features
//
// These call the real builders and OperationRunner. Each takes ~2 min
// because of MOPAC. This is the price of testing the real pipeline.
// ============================================================================

class JobSpecE2E : public ::testing::Test {
    // RuntimeEnvironment::Load() is called in test_main.cpp before all tests.
};


TEST_F(JobSpecE2E, PdbEndToEnd) {
    std::string out = "/tmp/jobspec_e2e_pdb";
    FakeArgv a{"test", "--pdb", (DATA + "/external/1UBQ.pdb").c_str(),
               "--output", out.c_str()};

    auto spec = ParseJobSpec(a.argc(), a.argv());
    ASSERT_TRUE(spec.Ok());
    ASSERT_TRUE(ValidateJobSpec(spec));

    // Build
    auto build = BuildFromPdb(spec.pdb_path, spec.pH);
    ASSERT_TRUE(build.Ok()) << build.error;

    auto& conf = build.protein->Conformation();

    // Run
    RunOptions opts;
    opts.charge_source = build.charges.get();
    opts.net_charge = build.net_charge;
    auto result = OperationRunner::Run(conf, opts);
    EXPECT_TRUE(result.Ok()) << result.error;
    EXPECT_GT(result.attached.size(), 5u);

    // Write features
    fs::create_directories(out);
    int arrays = ConformationResult::WriteAllFeatures(conf, out);
    EXPECT_GT(arrays, 0);
}


TEST_F(JobSpecE2E, ProtonatedPdbEndToEnd) {
    std::string out = "/tmp/jobspec_e2e_protonated";
    FakeArgv a{"test", "--protonated-pdb",
               (DATA + "/1ubq_protonated.pdb").c_str(),
               "--output", out.c_str()};

    auto spec = ParseJobSpec(a.argc(), a.argv());
    ASSERT_TRUE(spec.Ok());
    ASSERT_TRUE(ValidateJobSpec(spec));

    auto build = BuildFromProtonatedPdb(spec.pdb_path);
    ASSERT_TRUE(build.Ok()) << build.error;

    auto& conf = build.protein->Conformation();

    RunOptions opts;
    opts.charge_source = build.charges.get();
    opts.net_charge = build.net_charge;
    auto result = OperationRunner::Run(conf, opts);
    EXPECT_TRUE(result.Ok()) << result.error;

    fs::create_directories(out);
    int arrays = ConformationResult::WriteAllFeatures(conf, out);
    EXPECT_GT(arrays, 0);
}


TEST_F(JobSpecE2E, OrcaEndToEnd) {
    std::string root = DATA + "/orca/A0A7C5FAR6_WT";
    std::string out = "/tmp/jobspec_e2e_orca";
    FakeArgv a{"test", "--orca", "--root", root.c_str(),
               "--output", out.c_str()};

    auto spec = ParseJobSpec(a.argc(), a.argv());
    ASSERT_TRUE(spec.Ok());
    ASSERT_TRUE(ValidateJobSpec(spec));

    auto build = BuildFromOrca(spec.orca_files);
    ASSERT_TRUE(build.Ok()) << build.error;

    auto& conf = build.protein->Conformation();

    RunOptions opts;
    opts.charge_source = build.charges.get();
    opts.net_charge = build.net_charge;
    if (!spec.orca_files.nmr_out_path.empty())
        opts.orca_nmr_path = spec.orca_files.nmr_out_path;

    auto result = OperationRunner::Run(conf, opts);
    EXPECT_TRUE(result.Ok()) << result.error;

    fs::create_directories(out);
    int arrays = ConformationResult::WriteAllFeatures(conf, out);
    EXPECT_GT(arrays, 0);
}


TEST_F(JobSpecE2E, MutantEndToEnd) {
    std::string wt_root  = DATA + "/orca/A0A7C5FAR6_WT";
    std::string ala_root = DATA + "/orca/A0A7C5FAR6_ALA";
    std::string out = "/tmp/jobspec_e2e_mutant";
    FakeArgv a{"test", "--mutant", "--wt", wt_root.c_str(),
               "--ala", ala_root.c_str(), "--output", out.c_str()};

    auto spec = ParseJobSpec(a.argc(), a.argv());
    ASSERT_TRUE(spec.Ok());
    ASSERT_TRUE(ValidateJobSpec(spec));

    auto wt_build  = BuildFromOrca(spec.wt_files);
    auto ala_build = BuildFromOrca(spec.ala_files);
    ASSERT_TRUE(wt_build.Ok())  << wt_build.error;
    ASSERT_TRUE(ala_build.Ok()) << ala_build.error;

    auto& wt_conf  = wt_build.protein->Conformation();
    auto& ala_conf = ala_build.protein->Conformation();

    RunOptions wt_opts;
    wt_opts.charge_source = wt_build.charges.get();
    wt_opts.net_charge = wt_build.net_charge;
    if (!spec.wt_files.nmr_out_path.empty())
        wt_opts.orca_nmr_path = spec.wt_files.nmr_out_path;

    RunOptions ala_opts;
    ala_opts.charge_source = ala_build.charges.get();
    ala_opts.net_charge = ala_build.net_charge;
    if (!spec.ala_files.nmr_out_path.empty())
        ala_opts.orca_nmr_path = spec.ala_files.nmr_out_path;

    auto result = OperationRunner::RunMutantComparison(
        wt_conf, wt_opts, ala_conf, ala_opts);
    EXPECT_TRUE(result.Ok()) << result.error;

    fs::create_directories(out);
    int arrays = ConformationResult::WriteAllFeatures(wt_conf, out);
    EXPECT_GT(arrays, 0);
}


TEST_F(JobSpecE2E, FleetLibraryDirect) {
    // --fleet CLI mode removed 2026-04-12, but BuildFromGromacs stays in the
    // library for backward compatibility with pre-extracted PDB pose data.
    // This test calls it directly, bypassing JobSpec.
    std::string tpr   = DATA + "/fleet/1A6J_5789/params/prod.tpr";
    std::string poses = DATA + "/fleet/1A6J_5789/poses";
    std::string out   = "/tmp/jobspec_e2e_fleet";

    FleetPaths paths;
    paths.tpr_path          = tpr;
    paths.sampled_poses_dir = poses;

    auto build = BuildFromGromacs(paths);
    ASSERT_TRUE(build.Ok()) << build.error;

    EXPECT_GT(build.protein->ConformationCount(), 1u);

    RunOptions opts;
    opts.charge_source = build.charges.get();
    opts.net_charge = build.net_charge;

    auto results = OperationRunner::RunEnsemble(*build.protein, opts);
    EXPECT_FALSE(results.empty());

    fs::create_directories(out);
    int total = 0;
    for (size_t i = 0; i < build.protein->ConformationCount(); ++i) {
        auto& conf = build.protein->ConformationAt(i);
        std::string frame_dir = out + "/frame_" + std::to_string(i + 1);
        fs::create_directories(frame_dir);
        total += ConformationResult::WriteAllFeatures(conf, frame_dir);
    }
    EXPECT_GT(total, 0);
}

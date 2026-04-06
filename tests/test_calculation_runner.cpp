#include "TestEnvironment.h"
//
// test_operation_runner.cpp
//
// Tests the OperationRunner sequences against real data.
// Mirrors the use cases in spec/USE_CASES.md.
//

#include <gtest/gtest.h>
#include <filesystem>

#include "OperationRunner.h"
#include "PdbFileReader.h"
#include "OrcaRunLoader.h"
#include "ChargeSource.h"
#include "BiotSavartResult.h"
#include "CoulombResult.h"
#include "HBondResult.h"
#include "OrcaShieldingResult.h"
#include "MutationDeltaResult.h"

namespace fs = std::filesystem;
using namespace nmr;



// ============================================================================
// Use case A: protonated PDB with charges (no DFT)
// ============================================================================

TEST(OperationRunnerTest, RunWithCharges) {
    if (!fs::exists(nmr::test::TestEnvironment::UbqProtonated())) GTEST_SKIP() << "1UBQ not found";

    auto build = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(build.Ok()) << build.error;
    auto& conf = build.protein->Conformation();

    RunOptions opts;
    opts.charge_source = build.charges.get();
    opts.net_charge = build.net_charge;

    auto result = OperationRunner::Run(conf, opts);
    ASSERT_TRUE(result.Ok()) << result.error;

    EXPECT_TRUE(conf.HasResult<BiotSavartResult>());
    EXPECT_TRUE(conf.HasResult<CoulombResult>());
    EXPECT_TRUE(conf.HasResult<HBondResult>());
    EXPECT_GT(result.attached.size(), 10u);

    std::cout << "  RunWithCharges: " << result.attached.size() << " results\n";
}


// ============================================================================
// Use case B: single DFT with charges + ORCA shielding
// ============================================================================

TEST(OperationRunnerTest, RunSingleDft) {
    OrcaRunFiles files;
    files.pdb_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.pdb";
    files.xyz_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.xyz";
    files.prmtop_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.prmtop";
    std::string nmr_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT_nmr.out";

    if (!fs::exists(files.prmtop_path)) GTEST_SKIP();
    if (!fs::exists(files.xyz_path)) GTEST_SKIP();

    auto build = BuildFromOrca(files);
    ASSERT_TRUE(build.Ok()) << build.error;
    auto& conf = build.protein->Conformation();

    RunOptions opts;
    opts.charge_source = build.charges.get();
    opts.net_charge = build.net_charge;
    if (fs::exists(nmr_path)) opts.orca_nmr_path = nmr_path;

    auto result = OperationRunner::Run(conf, opts);
    ASSERT_TRUE(result.Ok()) << result.error;

    EXPECT_TRUE(conf.HasResult<BiotSavartResult>());
    EXPECT_TRUE(conf.HasResult<CoulombResult>());

    if (!opts.orca_nmr_path.empty()) {
        EXPECT_TRUE(conf.HasResult<OrcaShieldingResult>());
    }

    std::cout << "  RunSingleDft: " << result.attached.size() << " results\n";
}


// ============================================================================
// Use case C: WT + ALA mutant comparison
// ============================================================================

TEST(OperationRunnerTest, RunMutantComparison) {
    std::string dir = std::string(nmr::test::TestEnvironment::Consolidated()) + "P84477/";
    if (!fs::exists(dir)) GTEST_SKIP() << "P84477 not found";

    // WT
    OrcaRunFiles wt_files;
    wt_files.pdb_path = dir + "P84477_WT.pdb";
    wt_files.xyz_path = dir + "P84477_WT.xyz";
    wt_files.prmtop_path = dir + "P84477_WT.prmtop";
    std::string wt_nmr = dir + "P84477_WT_20260311_024106_nmr.out";

    // ALA
    OrcaRunFiles ala_files;
    ala_files.pdb_path = dir + "P84477_ALA.pdb";
    ala_files.xyz_path = dir + "P84477_ALA.xyz";
    ala_files.prmtop_path = dir + "P84477_ALA.prmtop";
    std::string ala_nmr = dir + "P84477_ALA_20260311_025630_nmr.out";

    if (!fs::exists(wt_files.prmtop_path) || !fs::exists(ala_files.prmtop_path))
        GTEST_SKIP();

    auto wt_build = BuildFromOrca(wt_files);
    auto ala_build = BuildFromOrca(ala_files);
    ASSERT_TRUE(wt_build.Ok()) << wt_build.error;
    ASSERT_TRUE(ala_build.Ok()) << ala_build.error;

    auto& wt_conf = wt_build.protein->Conformation();
    auto& ala_conf = ala_build.protein->Conformation();

    RunOptions wt_opts;
    wt_opts.charge_source = wt_build.charges.get();
    wt_opts.net_charge = wt_build.net_charge;
    if (fs::exists(wt_nmr)) wt_opts.orca_nmr_path = wt_nmr;

    RunOptions ala_opts;
    ala_opts.charge_source = ala_build.charges.get();
    ala_opts.net_charge = ala_build.net_charge;
    if (fs::exists(ala_nmr)) ala_opts.orca_nmr_path = ala_nmr;

    auto result = OperationRunner::RunMutantComparison(
        wt_conf, wt_opts, ala_conf, ala_opts);
    ASSERT_TRUE(result.Ok()) << result.error;

    EXPECT_TRUE(wt_conf.HasResult<BiotSavartResult>());
    EXPECT_TRUE(wt_conf.HasResult<CoulombResult>());

    if (fs::exists(wt_nmr) && fs::exists(ala_nmr)) {
        EXPECT_TRUE(wt_conf.HasResult<MutationDeltaResult>());
    }

    std::cout << "  MutantComparison: " << result.attached.size()
              << " results on WT\n";
}

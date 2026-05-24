#include "TestEnvironment.h"
//
// test_write_features.cpp
//
// Runs the full pipeline on P84477 and writes all features to NPY files.
// This is the baseline for the topology refactor: after refactoring,
// the output directory should be byte-identical.
//

#include <gtest/gtest.h>
#include <filesystem>

#include "ConformationResult.h"
#include "OperationRunner.h"
#include "OrcaRunLoader.h"
#include "ChargeSource.h"
#include "OperationLog.h"

namespace fs = std::filesystem;
using namespace nmr;



TEST(WriteFeatures, P84477Baseline) {
    std::string dir = std::string(nmr::test::TestEnvironment::Consolidated()) + "P84477/";
    if (!fs::exists(dir)) GTEST_SKIP() << "P84477 not found";

    // Load WT with full pipeline
    OrcaRunFiles files;
    files.pdb_path = dir + "P84477_WT.pdb";
    files.xyz_path = dir + "P84477_WT.xyz";
    files.prmtop_path = dir + "P84477_WT.prmtop";

    auto load = BuildFromOrca(files);
    ASSERT_TRUE(load.Ok()) << load.error;
    auto& conf = load.protein->Conformation();

    RunOptions opts;
    opts.charge_source = load.charges.get();
    opts.net_charge = load.net_charge;

    // Find NMR output
    for (const auto& entry : fs::directory_iterator(dir)) {
        std::string name = entry.path().filename().string();
        if (name.find("P84477_WT") == 0 && name.find("_nmr.out") != std::string::npos) {
            opts.orca_nmr_path = entry.path().string();
            break;
        }
    }

    auto result = OperationRunner::Run(conf, opts);
    ASSERT_TRUE(result.Ok()) << result.error;

    // Write all features
    int arrays = ConformationResult::WriteAllFeatures(conf, nmr::test::TestEnvironment::BaselineFeatures());
    EXPECT_GT(arrays, 25) << "Expected 25+ arrays from 8 calculators + identity";

    std::cout << "\n  WriteFeatures baseline written to " << nmr::test::TestEnvironment::BaselineFeatures() << "\n"
              << "  Arrays: " << arrays << "\n"
              << "  Atoms: " << conf.AtomCount() << "\n"
              << "  Results attached: " << result.attached.size() << "\n";

    // List what we wrote
    for (const auto& entry : fs::directory_iterator(nmr::test::TestEnvironment::BaselineFeatures())) {
        auto sz = fs::file_size(entry);
        std::cout << "    " << entry.path().filename().string()
                  << " (" << sz << " bytes)\n";
    }
}

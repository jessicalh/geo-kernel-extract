//
// Test main: loads RuntimeEnvironment, TestEnvironment, and
// CalculatorConfig before any tests. CalculatorConfig is loaded so
// that JobSpecE2E and any other test exercising ValidateJobSpec can
// resolve the AIMNet2 TOML default (aimnet2_model_path) — required
// since the 2026-05-04 contract scope generalisation per
// feedback_aimnet2_required_no_weasel.
//

#include "RuntimeEnvironment.h"
#include "TestEnvironment.h"
#include "CalculatorConfig.h"
#include <gtest/gtest.h>
#include <string>

#ifndef NMR_TEST_DATA_DIR
#error "NMR_TEST_DATA_DIR must be defined"
#endif

int main(int argc, char** argv) {
    nmr::RuntimeEnvironment::Load();
    nmr::test::TestEnvironment::Load();
    nmr::CalculatorConfig::Load(
        std::string(NMR_TEST_DATA_DIR) + "/../../data/calculator_params.toml");
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

//
// Test main: loads RuntimeEnvironment and TestEnvironment before any tests.
//

#include "RuntimeEnvironment.h"
#include "TestEnvironment.h"
#include <gtest/gtest.h>

int main(int argc, char** argv) {
    nmr::RuntimeEnvironment::Load();
    nmr::test::TestEnvironment::Load();
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

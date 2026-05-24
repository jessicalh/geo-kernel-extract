#include <gtest/gtest.h>
#include "RuntimeEnvironment.h"
#include <filesystem>

using namespace nmr;
namespace fs = std::filesystem;

// Load() is called by test_main.cpp before any tests run.
// These tests verify the live surface (ff14sb_params, tmpdir).
// MOPAC is linked (libmopac.so), not a binary path — no accessor needed.
// Removed: xtb (replaced by libmopac), propka, pdb2pqr, tleap, kaml, gmx.

TEST(RuntimeEnvironment, LoadSucceeds) {
    // Load was already called by test_main. RequireLoaded should pass.
    EXPECT_TRUE(RuntimeEnvironment::RequireLoaded());
}

TEST(RuntimeEnvironment, Ff14sbParamsExists) {
    const auto& path = RuntimeEnvironment::Ff14sbParams();
    EXPECT_FALSE(path.empty()) << "ff14sb_params not configured";
    EXPECT_TRUE(fs::exists(path)) << "ff14sb_params not found at: " << path;
}

TEST(RuntimeEnvironment, VerifyFindsNoMissing) {
    auto missing = RuntimeEnvironment::Verify();
    for (const auto& m : missing) {
        ADD_FAILURE() << "Missing: " << m;
    }
}

TEST(RuntimeEnvironment, TempFilePathIncludesGuidAndProtein) {
    std::string path = RuntimeEnvironment::TempFilePath("1UBQ", "protonated.pdb");
    EXPECT_NE(path.find("1UBQ"), std::string::npos);
    EXPECT_NE(path.find("protonated.pdb"), std::string::npos);
    EXPECT_NE(path.find(RuntimeEnvironment::TmpDir()), std::string::npos);
}

TEST(RuntimeEnvironment, TempFilePathIsUniquePerProcess) {
    std::string p1 = RuntimeEnvironment::TempFilePath("1UBQ", "a.pdb");
    std::string p2 = RuntimeEnvironment::TempFilePath("2GB1", "a.pdb");
    EXPECT_NE(p1, p2);
    auto extract_guid = [](const std::string& path) {
        auto pos = path.rfind('/');
        auto under = path.find('_', pos + 1);
        return path.substr(pos + 1, under - pos - 1);
    };
    EXPECT_EQ(extract_guid(p1), extract_guid(p2));
}

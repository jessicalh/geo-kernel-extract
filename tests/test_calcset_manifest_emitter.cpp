#include <gtest/gtest.h>
#include <nlohmann/json.hpp>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <regex>
#include <string>

#include "CalcsetManifestEmitter.h"
#include "CategoryInfoProjection.h"
#include "PdbFileReader.h"
#include "TestEnvironment.h"
#include "TopologySidecar.h"

namespace fs = std::filesystem;

namespace {

fs::path MakeTempRoot() {
    const auto stamp = std::chrono::steady_clock::now()
        .time_since_epoch()
        .count();
    return fs::temp_directory_path() /
        ("calcset_manifest_emitter_test_" + std::to_string(stamp));
}

void RemoveFlatTempDir(const fs::path& root) {
    std::error_code ec;
    if (fs::is_directory(root, ec)) {
        for (const auto& entry : fs::directory_iterator(root, ec)) {
            if (ec) break;
            fs::remove(entry.path(), ec);
        }
    }
    fs::remove(root, ec);
}

TEST(CalcsetManifestEmitterTest, SinglePoseManifestParsesAndPointsAtSidecar) {
    if (!fs::exists(nmr::test::TestEnvironment::UbqProtonated())) {
        GTEST_SKIP() << "1UBQ PDB not found";
    }

    auto build = nmr::BuildFromProtonatedPdb(
        nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(build.Ok()) << build.error;

    const fs::path root = MakeTempRoot();
    fs::create_directories(root);

    ASSERT_EQ(1, nmr::CategoryInfoProjection::WriteFeatures(
        *build.protein, root.string()));
    ASSERT_EQ(5, nmr::TopologySidecar::WriteFeatures(
        *build.protein, root.string(), "1UBQ"));

    nmr::CalcsetManifestEmitter::SinglePoseArtifacts artifacts;
    artifacts.pose_kind = "protonated_pdb";
    artifacts.pose_dir = root;
    artifacts.extraction_manifest = root / "extraction_manifest.json";

    const auto result = nmr::CalcsetManifestEmitter::WriteSinglePose(
        root,
        nmr::CalcsetManifestEmitter::IdentityFromProteinId("1UBQ"),
        artifacts);
    ASSERT_TRUE(result.ok) << result.error;

    EXPECT_TRUE(std::regex_match(
        result.path.filename().string(),
        std::regex(R"(^1UBQ_[0-9]{8}T[0-9]{6}Z\.lgs$)")));

    std::ifstream in(result.path);
    ASSERT_TRUE(in.good());
    const auto j = nlohmann::ordered_json::parse(in);

    ASSERT_TRUE(j.is_object());
    EXPECT_EQ(1, j.at("schema_version").get<int>());
    EXPECT_EQ("single_pose", j.at("kind").get<std::string>());
    EXPECT_EQ("1UBQ", j.at("dataset_id").get<std::string>());
    EXPECT_EQ("1UBQ", j.at("protein_id").get<std::string>());
    EXPECT_EQ("1UBQ calcset", j.at("human_name").get<std::string>());

    ASSERT_TRUE(j.contains("single_pose"));
    const auto& single = j.at("single_pose");
    EXPECT_EQ("protonated_pdb", single.at("pose_kind").get<std::string>());

    const fs::path manifest_rel =
        single.at("extraction_manifest").get<std::string>();
    EXPECT_FALSE(manifest_rel.is_absolute());
    EXPECT_TRUE(fs::is_regular_file(root / manifest_rel))
        << "unresolved artifact pointer: " << manifest_rel;

    ASSERT_TRUE(j.contains("metadata"));
    EXPECT_TRUE(j.at("metadata").contains("generated_at_utc"));
    EXPECT_TRUE(j.at("metadata").contains("lgs_writer"));

    RemoveFlatTempDir(root);
}

}  // namespace

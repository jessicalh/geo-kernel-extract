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

void RemoveTempTree(const fs::path& root) {
    std::error_code ec;
    if (!fs::is_directory(root, ec)) {
        fs::remove(root, ec);
        return;
    }
    for (const auto& entry : fs::directory_iterator(root, ec)) {
        if (ec) break;
        RemoveTempTree(entry.path());
    }
    fs::remove(root, ec);
}

void WriteTextFile(const fs::path& path, const std::string& text) {
    fs::create_directories(path.parent_path());
    std::ofstream out(path);
    out << text;
}

nlohmann::ordered_json ReadJson(const fs::path& path) {
    std::ifstream in(path);
    EXPECT_TRUE(in.good()) << path;
    return nlohmann::ordered_json::parse(in);
}

void ExpectCurrentAppMetadata(const nlohmann::ordered_json& j) {
    ASSERT_TRUE(j.contains("metadata"));
    const auto& metadata = j.at("metadata");
    EXPECT_TRUE(metadata.contains("generated_at_utc"));
    EXPECT_TRUE(metadata.contains("lgs_writer"));
    EXPECT_EQ("0.0.1", metadata.at("app_version").get<std::string>());
    EXPECT_EQ("whole_app", metadata.at("app_version_scope").get<std::string>());
    EXPECT_FALSE(metadata.contains("producer_extractor_version"));
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

    const auto j = ReadJson(result.path);

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

    ExpectCurrentAppMetadata(j);

    const auto sidecar = ReadJson(root / "extraction_manifest.json");
    EXPECT_FALSE(sidecar.contains("extractor_version"));
    EXPECT_FALSE(sidecar.contains("app_version"));

    RemoveFlatTempDir(root);
}

TEST(CalcsetManifestEmitterTest, TrajectoryManifestUsesWholeAppMetadata) {
    const fs::path root = MakeTempRoot();
    fs::create_directories(root / "md");
    fs::create_directories(root / "extract");
    WriteTextFile(root / "topol.top", "; topology\n");
    WriteTextFile(root / "extract" / "trajectory.h5", "not parsed\n");
    WriteTextFile(root / "extract" / "extraction_manifest.json",
                  R"({"schema_version":"1.0","extractor":"nmr_extract","protein_id":"P"})");

    nmr::CalcsetManifestEmitter::TrajectoryArtifacts artifacts;
    artifacts.md_dir = root / "md";
    artifacts.topology_top = root / "topol.top";
    artifacts.extraction_dir = root / "extract";
    artifacts.trajectory_h5 = root / "extract" / "trajectory.h5";
    artifacts.extraction_manifest = root / "extract" / "extraction_manifest.json";
    artifacts.frame_dt_ps = 10.0;

    const auto result = nmr::CalcsetManifestEmitter::WriteTrajectory(
        root,
        nmr::CalcsetManifestEmitter::IdentityFromProteinId("P"),
        artifacts);
    ASSERT_TRUE(result.ok) << result.error;

    const auto j = ReadJson(result.path);
    EXPECT_EQ("trajectory", j.at("kind").get<std::string>());
    ExpectCurrentAppMetadata(j);
    EXPECT_FALSE(j.at("metadata").contains("producer_extractor_version"));

    RemoveTempTree(root);
}

TEST(CalcsetManifestEmitterTest, MutantPairManifestUsesWholeAppMetadata) {
    const fs::path root = MakeTempRoot();
    fs::create_directories(root);
    WriteTextFile(root / "WT.lgs", "{}\n");
    WriteTextFile(root / "ALA.lgs", "{}\n");

    nmr::CalcsetManifestEmitter::MutantPairArtifacts artifacts;
    artifacts.wt_lgs = root / "WT.lgs";
    artifacts.ala_lgs = root / "ALA.lgs";

    const auto result = nmr::CalcsetManifestEmitter::WriteMutantPair(
        root,
        nmr::CalcsetManifestEmitter::IdentityFromProteinId("WT_ALA"),
        artifacts);
    ASSERT_TRUE(result.ok) << result.error;

    const auto j = ReadJson(result.path);
    EXPECT_EQ("mutant_pair", j.at("kind").get<std::string>());
    ExpectCurrentAppMetadata(j);
    EXPECT_FALSE(j.at("metadata").contains("producer_extractor_version"));

    RemoveTempTree(root);
}

}  // namespace

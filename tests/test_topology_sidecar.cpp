// tests/test_topology_sidecar.cpp
//
// Tests for src/TopologySidecar.{h,cpp}.
//
// Coverage:
//   1. WriteFeatures emits four files: bonds.npy, rings.npy,
//      ring_membership.npy, extraction_manifest.json.
//   2. bonds.npy row count == LegacyAmber::BondCount().
//   3. rings.npy row count == AromaticRingCount() + SaturatedRingCount().
//   4. ring_membership.npy row count == sum of per-ring atom_indices.size().
//   5. extraction_manifest.json contains the expected top-level keys and
//      axis-size fields.
//   6. Each bond's atom_index_a / atom_index_b < AtomCount() (boundary
//      sanity at the byte level — codex's first-pass validation gate).
//

#include <gtest/gtest.h>

#include <array>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "TestEnvironment.h"
#include "TopologySidecar.h"
#include "Protein.h"
#include "PdbFileReader.h"
#include "LegacyAmberTopology.h"
#include "RingTopology.h"
#include "Ring.h"

namespace fs = std::filesystem;

namespace {

using nmr::Protein;
using nmr::TopologySidecar;

// Read whole file into a vector of bytes.
std::vector<char> ReadFileBytes(const fs::path& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return {};
    return std::vector<char>(std::istreambuf_iterator<char>(f),
                              std::istreambuf_iterator<char>());
}

// Read whole file into a string (text).
std::string ReadFileText(const fs::path& path) {
    std::ifstream f(path);
    if (!f) return {};
    return std::string(std::istreambuf_iterator<char>(f),
                        std::istreambuf_iterator<char>());
}

// Parse NPY 1.0 header to get the declared shape's first dimension.
// Returns SIZE_MAX on parse failure.
std::size_t ReadNpyRowCount(const std::vector<char>& bytes) {
    if (bytes.size() < 10) return SIZE_MAX;
    if (static_cast<unsigned char>(bytes[0]) != 0x93) return SIZE_MAX;
    const uint16_t header_len = static_cast<uint8_t>(bytes[8]) |
                                 (static_cast<uint8_t>(bytes[9]) << 8);
    if (bytes.size() < static_cast<std::size_t>(10 + header_len)) return SIZE_MAX;
    const std::string header(bytes.data() + 10, header_len);
    const auto pos = header.find("'shape': (");
    if (pos == std::string::npos) return SIZE_MAX;
    const auto start = pos + 10;
    const auto comma = header.find(',', start);
    if (comma == std::string::npos) return SIZE_MAX;
    return std::stoul(header.substr(start, comma - start));
}

class TopologySidecarTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!fs::exists(nmr::test::TestEnvironment::UbqProtonated())) {
            GTEST_SKIP() << "1UBQ PDB not found";
        }
        auto build = nmr::BuildFromProtonatedPdb(
            nmr::test::TestEnvironment::UbqProtonated());
        ASSERT_TRUE(build.Ok()) << build.error;
        protein_ = std::move(build.protein);

        dir_ = fs::temp_directory_path() /
            ("topology_sidecar_test_" + std::to_string(::getpid()));
        fs::create_directories(dir_);

        const int written = TopologySidecar::WriteFeatures(
            *protein_, dir_.string(), "1UBQ");
        ASSERT_EQ(5, written);
    }
    void TearDown() override {
        // NOTE: cannot call fs::remove_all here -- libtorch ships a
        // broken std::filesystem::remove_all that segfaults via a
        // misresolved symbol when the test executable links torch
        // (see feedback_libtorch_broken_filesystem memory). Remove
        // individual files instead.
        std::error_code ec;
        fs::remove(dir_ / "bonds.npy", ec);
        fs::remove(dir_ / "rings.npy", ec);
        fs::remove(dir_ / "ring_membership.npy", ec);
        fs::remove(dir_ / "extraction_manifest.json", ec);
        fs::remove(dir_, ec);
    }

    std::unique_ptr<Protein> protein_;
    fs::path dir_;
};

// ============================================================================
// File emission
// ============================================================================

TEST_F(TopologySidecarTest, EmitsAllFiveFiles) {
    EXPECT_TRUE(fs::exists(dir_ / "residues.npy"));
    EXPECT_TRUE(fs::exists(dir_ / "bonds.npy"));
    EXPECT_TRUE(fs::exists(dir_ / "rings.npy"));
    EXPECT_TRUE(fs::exists(dir_ / "ring_membership.npy"));
    EXPECT_TRUE(fs::exists(dir_ / "extraction_manifest.json"));
}

// ============================================================================
// Axis row counts
// ============================================================================

TEST_F(TopologySidecarTest, BondsAxisCountMatchesTopology) {
    const auto bytes = ReadFileBytes(dir_ / "bonds.npy");
    const std::size_t rows = ReadNpyRowCount(bytes);
    EXPECT_EQ(protein_->LegacyAmber().BondCount(), rows);
}

TEST_F(TopologySidecarTest, RingsAxisCountMatchesTopology) {
    const auto bytes = ReadFileBytes(dir_ / "rings.npy");
    const std::size_t rows = ReadNpyRowCount(bytes);
    const auto& rt = protein_->LegacyAmber().Rings();
    EXPECT_EQ(rt.AromaticCount() + rt.SaturatedCount(), rows);
}

TEST_F(TopologySidecarTest, RingMembershipAxisCountMatchesTopology) {
    const auto bytes = ReadFileBytes(dir_ / "ring_membership.npy");
    const std::size_t rows = ReadNpyRowCount(bytes);

    std::size_t expected = 0;
    const auto& rt = protein_->LegacyAmber().Rings();
    for (std::size_t i = 0; i < rt.AromaticCount(); ++i) {
        expected += rt.AromaticAt(i).atom_indices.size();
    }
    for (std::size_t i = 0; i < rt.SaturatedCount(); ++i) {
        expected += rt.SaturatedAt(i).atom_indices.size();
    }
    EXPECT_EQ(expected, rows);
}

// ============================================================================
// Bond endpoint sanity (codex first-pass validation gate)
// ============================================================================

TEST_F(TopologySidecarTest, BondEndpointsAreWithinAtomAxis) {
    const auto bytes = ReadFileBytes(dir_ / "bonds.npy");
    ASSERT_FALSE(bytes.empty());
    const uint16_t header_len = static_cast<uint8_t>(bytes[8]) |
                                 (static_cast<uint8_t>(bytes[9]) << 8);
    const std::size_t data_start = 10 + header_len;
    const std::size_t row_size = 4 + 4 + 4 + 1 + 1 + 1 + 1 + 1 + 1;
    const std::size_t N = (bytes.size() - data_start) / row_size;
    ASSERT_EQ(N, protein_->LegacyAmber().BondCount());

    const std::size_t atom_count = protein_->AtomCount();
    for (std::size_t bi = 0; bi < N; ++bi) {
        const char* row = bytes.data() + data_start + bi * row_size;
        int32_t bond_index, atom_a, atom_b;
        std::memcpy(&bond_index, row + 0, 4);
        std::memcpy(&atom_a, row + 4, 4);
        std::memcpy(&atom_b, row + 8, 4);
        EXPECT_EQ(static_cast<int32_t>(bi), bond_index);
        EXPECT_GE(atom_a, 0);
        EXPECT_LT(static_cast<std::size_t>(atom_a), atom_count);
        EXPECT_GE(atom_b, 0);
        EXPECT_LT(static_cast<std::size_t>(atom_b), atom_count);
    }
}

// ============================================================================
// Manifest
// ============================================================================

TEST_F(TopologySidecarTest, ManifestHasExpectedKeys) {
    const std::string j = ReadFileText(dir_ / "extraction_manifest.json");
    ASSERT_FALSE(j.empty());

    EXPECT_NE(j.find("\"schema_version\""), std::string::npos);
    EXPECT_NE(j.find("\"extractor\""), std::string::npos);
    EXPECT_NE(j.find("\"generated_at_utc\""), std::string::npos);
    EXPECT_NE(j.find("\"protein_id\": \"1UBQ\""), std::string::npos);
    EXPECT_NE(j.find("\"topology\""), std::string::npos);
    EXPECT_NE(j.find("\"axis_sizes\""), std::string::npos);
    EXPECT_NE(j.find("\"axis_alignment\""), std::string::npos);
    EXPECT_NE(j.find("\"has_atom_semantic\""), std::string::npos);
    EXPECT_NE(j.find("\"atom\""), std::string::npos);
    EXPECT_NE(j.find("\"residue\""), std::string::npos);
    EXPECT_NE(j.find("\"bond\""), std::string::npos);
    EXPECT_NE(j.find("\"ring\""), std::string::npos);
    EXPECT_NE(j.find("\"aromatic_ring\""), std::string::npos);
    EXPECT_NE(j.find("\"saturated_ring\""), std::string::npos);
}

TEST_F(TopologySidecarTest, ManifestAxisSizesMatchProtein) {
    const std::string j = ReadFileText(dir_ / "extraction_manifest.json");
    ASSERT_FALSE(j.empty());

    // Lightweight string-find checks. We're not parsing JSON; we're
    // verifying that the declared axis sizes match the protein's actual
    // counts.
    const std::string atom_decl =
        "\"atom\": " + std::to_string(protein_->AtomCount());
    EXPECT_NE(j.find(atom_decl), std::string::npos)
        << "manifest must declare atom axis size " << protein_->AtomCount();

    const std::string res_decl =
        "\"residue\": " + std::to_string(protein_->ResidueCount());
    EXPECT_NE(j.find(res_decl), std::string::npos);

    const std::string bond_decl =
        "\"bond\": " + std::to_string(protein_->LegacyAmber().BondCount());
    EXPECT_NE(j.find(bond_decl), std::string::npos);
}

// ============================================================================
// Idempotency
// ============================================================================

TEST_F(TopologySidecarTest, IdempotentRewrite) {
    const auto bytes_before = ReadFileBytes(dir_ / "bonds.npy");
    ASSERT_FALSE(bytes_before.empty());

    // Run a second time. Manifest's generated_at_utc will differ; bonds,
    // rings, and ring_membership are topology-only and must be byte-
    // identical.
    const int written = TopologySidecar::WriteFeatures(
        *protein_, dir_.string(), "1UBQ");
    EXPECT_EQ(5, written);

    const auto bytes_after = ReadFileBytes(dir_ / "bonds.npy");
    EXPECT_EQ(bytes_before, bytes_after);
}

// ============================================================================
// Manifest content — value-bearing assertions
//
// The previous round of tests checked only key presence; that caught
// missing fields but missed the bad axis_alignment claim (residue_type.npy
// was wrongly called residue-axis). Pull integer values out by hand for
// axis_sizes and assert against the protein's actual surface; assert
// alignment claims by their specific anchoring fact.
// ============================================================================

// Hand-rolled "value-after-key" extractor for the flat axis_sizes block.
// Codex's JSON manifest is project-controlled; full JSON parsing is
// overkill for a self-emitted format.
int ExtractAxisSize(const std::string& j, const std::string& axis) {
    const std::string needle = "\"" + axis + "\":";
    auto pos = j.find(needle);
    if (pos == std::string::npos) return -1;
    pos += needle.size();
    while (pos < j.size() && std::isspace(static_cast<unsigned char>(j[pos]))) ++pos;
    int value = 0;
    bool any = false;
    while (pos < j.size() && std::isdigit(static_cast<unsigned char>(j[pos]))) {
        value = value * 10 + (j[pos] - '0');
        ++pos;
        any = true;
    }
    return any ? value : -1;
}

TEST_F(TopologySidecarTest, ManifestAxisSizesExactlyMatchProtein) {
    const std::string j = ReadFileText(dir_ / "extraction_manifest.json");
    ASSERT_FALSE(j.empty());

    EXPECT_EQ(static_cast<int>(protein_->AtomCount()),
              ExtractAxisSize(j, "atom"));
    EXPECT_EQ(static_cast<int>(protein_->ResidueCount()),
              ExtractAxisSize(j, "residue"));
    EXPECT_EQ(static_cast<int>(protein_->LegacyAmber().BondCount()),
              ExtractAxisSize(j, "bond"));
    EXPECT_EQ(static_cast<int>(protein_->LegacyAmber().AromaticRingCount()),
              ExtractAxisSize(j, "aromatic_ring"));
    EXPECT_EQ(static_cast<int>(protein_->LegacyAmber().SaturatedRingCount()),
              ExtractAxisSize(j, "saturated_ring"));
    EXPECT_EQ(static_cast<int>(protein_->LegacyAmber().AromaticRingCount() +
                                  protein_->LegacyAmber().SaturatedRingCount()),
              ExtractAxisSize(j, "ring"));
    const auto aromatic = protein_->LegacyAmber().AromaticRingCount();
    EXPECT_EQ(static_cast<int>(aromatic * (aromatic > 0 ? aromatic - 1 : 0) / 2),
              ExtractAxisSize(j, "ring_pair"));

    // ring_membership axis size must match the sum of per-ring vertex counts.
    std::size_t expected_membership = 0;
    const auto& rt = protein_->LegacyAmber().Rings();
    for (std::size_t i = 0; i < rt.AromaticCount(); ++i) {
        expected_membership += rt.AromaticAt(i).atom_indices.size();
    }
    for (std::size_t i = 0; i < rt.SaturatedCount(); ++i) {
        expected_membership += rt.SaturatedAt(i).atom_indices.size();
    }
    EXPECT_EQ(static_cast<int>(expected_membership),
              ExtractAxisSize(j, "ring_membership"));
}

TEST_F(TopologySidecarTest, ManifestEnumVocabPresent) {
    const std::string j = ReadFileText(dir_ / "extraction_manifest.json");
    ASSERT_FALSE(j.empty());

    // enum_vocab block resolves codex's enum_vocab_refs requirement.
    // Consumers reading bond_order == 4 need to know that's Peptide.
    EXPECT_NE(j.find("\"enum_vocab\""), std::string::npos);
    EXPECT_NE(j.find("\"terminal_state\""), std::string::npos);
    EXPECT_NE(j.find("\"bond_order\""), std::string::npos);
    EXPECT_NE(j.find("\"bond_category\""), std::string::npos);
    EXPECT_NE(j.find("\"ring_kind\""), std::string::npos);
    EXPECT_NE(j.find("\"ring_type_index\""), std::string::npos);
    EXPECT_NE(j.find("\"residue_type\""), std::string::npos);
    // Spot-check specific value bindings. (Pretty-printed by nlohmann
    // ordered_json: "key": "value" with a space after the colon.)
    EXPECT_NE(j.find("\"4\": \"Peptide\""), std::string::npos);
    EXPECT_NE(j.find("\"0\": \"PheBenzene\""), std::string::npos);
    EXPECT_NE(j.find("\"14\": \"PRO\""), std::string::npos);
}

TEST_F(TopologySidecarTest, ManifestAxisAlignmentClaimsAreCorrect) {
    const std::string j = ReadFileText(dir_ / "extraction_manifest.json");
    ASSERT_FALSE(j.empty());

    // residue alignment: residues.npy is the canonical residue axis;
    // residue_type.npy in the identity block is atom-axis. The previous
    // wording said "residue_type.npy is residue-axis" which is wrong;
    // catch any regression.
    EXPECT_NE(j.find("residues.npy is the canonical residue axis"),
              std::string::npos);
    EXPECT_NE(j.find("residue_type.npy / residue_index.npy in the identity block are atom-axis"),
              std::string::npos);

    // bond alignment: bonds.npy is the canonical bond axis.
    EXPECT_NE(j.find("bonds.npy is the canonical bond axis"),
              std::string::npos);

    // ring alignment: aromatic first then saturated, ring_id absolute.
    EXPECT_NE(j.find("aromatic rings first"),
              std::string::npos);
    EXPECT_NE(j.find("ring_id is the absolute row index"),
              std::string::npos);

    // aromatic_ring alignment: ring_geometry.npy ↔ rings.npy aromatic prefix.
    EXPECT_NE(j.find("ring_geometry.npy is aromatic-only"),
              std::string::npos);
}

// Residues axis count
TEST_F(TopologySidecarTest, ResiduesAxisCountMatchesTopology) {
    const auto bytes = ReadFileBytes(dir_ / "residues.npy");
    const std::size_t rows = ReadNpyRowCount(bytes);
    EXPECT_EQ(protein_->ResidueCount(), rows);
}

// Durable, fixture-independent freeze gate for the producer-side calibration
// and provisional-constant contract (M10/M12).  This exercises the actual
// TopologySidecar writer on a valid zero-atom topology rather than rebuilding
// its JSON formula in the test.
TEST(TopologySidecarManifestContract, EmitsResolvedRingScalingAndStatuses) {
    Protein protein;
    protein.FinalizeConstruction({});

    const fs::path dir = fs::temp_directory_path() /
        ("topology_manifest_contract_" + std::to_string(::getpid()));
    fs::create_directories(dir);
    ASSERT_EQ(5, TopologySidecar::WriteFeatures(
        protein, dir.string(), "manifest-contract"));

    std::ifstream in(dir / "extraction_manifest.json");
    ASSERT_TRUE(in.is_open());
    const nlohmann::ordered_json manifest =
        nlohmann::ordered_json::parse(in);

    const auto& ring = manifest.at("feature_metadata").at("ring_current");
    EXPECT_EQ(ring.at("ring_type_order"),
              nlohmann::ordered_json({"PheBenzene", "TyrPhenol",
                                      "TrpBenzene", "TrpPyrrole",
                                      "TrpPerimeter", "HisImidazole",
                                      "HidImidazole", "HieImidazole"}));
    EXPECT_EQ(ring.at("resolved_intensity_nA_per_T"),
              nlohmann::ordered_json({-12.0, -11.28, -12.48, -6.72,
                                      -19.2, -5.16, -5.16, -5.16}));
    EXPECT_EQ(ring.at("resolved_jb_lobe_offset_A"),
              nlohmann::ordered_json({0.64, 0.64, 0.64, 0.52,
                                      0.60, 0.50, 0.50, 0.50}));
    ASSERT_EQ(ring.at("resolved_intensity_source").size(), 8u);
    for (const auto& source : ring.at("resolved_intensity_source")) {
        EXPECT_TRUE(source == "default" || source == "toml");
    }

    const auto& intensity_constants =
        ring.at("constants").at("ring_current_intensity");
    ASSERT_EQ(intensity_constants.size(), 8u);
    EXPECT_EQ(intensity_constants.at(4).at("status"), "estimated");
    EXPECT_EQ(intensity_constants.at(4).at("literature_status"),
              "good_enough");
    EXPECT_EQ(intensity_constants.at(6).at("status"), "good_enough");
    EXPECT_EQ(intensity_constants.at(7).at("status"), "good_enough");
    EXPECT_EQ(intensity_constants.at(4).at("ring_type"), "TrpPerimeter");

    const auto& pi = manifest.at("feature_metadata")
                         .at("pi_quadrupole").at("physical_prefactor");
    EXPECT_DOUBLE_EQ(pi.at("value").get<double>(), 1.0);
    EXPECT_EQ(pi.at("status"), "deferred_learnable");
    EXPECT_EQ(pi.at("affected_ring_types").size(), 8u);

    const auto& larsen =
        manifest.at("feature_metadata").at("larsen_hbond");
    EXPECT_EQ(larsen.at("imputed_policy"), "emitted_unmasked");
    EXPECT_EQ(larsen.at("acceptor_archive_approximations")
                  .at("SidechainCarbonyl"),
              "BackboneCarbonyl/NMA archive approximation");
    const std::string larsen_pair_axis =
        manifest.at("axis_alignment").at("larsen_hbond_pair")
            .get<std::string>();
    EXPECT_NE(larsen_pair_axis.find("larsen_hbond_pairs_index.npy"),
              std::string::npos);
    for (const char* disposition : std::array{
             "0=missing_frame_atoms",
             "1=theta_out_of_range",
             "2=grid_miss",
             "3=success",
             "4=invalid_frame",
             "5=carboxylate_symmetry_filtered"}) {
        EXPECT_NE(larsen_pair_axis.find(disposition), std::string::npos)
            << disposition;
    }
    EXPECT_NE(manifest.at("axis_alignment").at("atom")
                  .get<std::string>().find("separate larsen_hbond_pair axis"),
              std::string::npos);

    std::error_code ec;
    fs::remove(dir / "residues.npy", ec);
    fs::remove(dir / "bonds.npy", ec);
    fs::remove(dir / "rings.npy", ec);
    fs::remove(dir / "ring_membership.npy", ec);
    fs::remove(dir / "extraction_manifest.json", ec);
    fs::remove(dir, ec);
}

}  // namespace

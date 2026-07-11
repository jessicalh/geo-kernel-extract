#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>

#include "RingSusceptibilityResult.h"
#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "CalculatorConfig.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "PdbFileReader.h"
#include "OrcaRunLoader.h"
#include "PhysicalConstants.h"

#include <filesystem>
namespace fs = std::filesystem;
using namespace nmr;

namespace {

struct NpyDoubleArray {
    std::vector<size_t> shape;
    std::vector<double> values;
};

std::string Trim(std::string value) {
    auto is_space = [](unsigned char c) { return std::isspace(c) != 0; };
    value.erase(value.begin(),
                std::find_if(value.begin(), value.end(),
                             [&](char c) { return !is_space(c); }));
    value.erase(std::find_if(value.rbegin(), value.rend(),
                             [&](char c) { return !is_space(c); }).base(),
                value.end());
    return value;
}

NpyDoubleArray ReadFloat64Npy(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    EXPECT_TRUE(in.is_open()) << path;
    NpyDoubleArray out;
    if (!in.is_open()) return out;
    char magic[6] = {};
    in.read(magic, 6);
    EXPECT_EQ(std::string(magic, 6), std::string("\x93NUMPY", 6));
    char version[2] = {};
    in.read(version, 2);
    std::uint16_t header_length = 0;
    in.read(reinterpret_cast<char*>(&header_length), sizeof(header_length));
    std::string header(header_length, '\0');
    in.read(header.data(), header_length);
    const auto begin = header.find('(');
    const auto end = header.find(')', begin);
    EXPECT_NE(begin, std::string::npos);
    EXPECT_NE(end, std::string::npos);
    if (begin == std::string::npos || end == std::string::npos) return out;
    std::stringstream stream(header.substr(begin + 1, end - begin - 1));
    std::string token;
    size_t count = 1;
    while (std::getline(stream, token, ',')) {
        token = Trim(token);
        if (token.empty()) continue;
        const size_t extent = static_cast<size_t>(std::stoull(token));
        out.shape.push_back(extent);
        count *= extent;
    }
    if (out.shape.empty()) count = 0;
    out.values.resize(count);
    if (count > 0) {
        in.read(reinterpret_cast<char*>(out.values.data()),
                static_cast<std::streamsize>(count * sizeof(double)));
    }
    return out;
}

}  // namespace



// ============================================================================
// Analytical test: verify the kept scalar rescue.
// ============================================================================

TEST(RingChiAnalytical, PointCenterKernelAcceptsInsideFiniteRingRadius) {
    const auto kernel = ring_susceptibility_detail::ComputeKernel(
        Vec3(0.0, 0.0, 0.5), Vec3::Zero(), Vec3(0.0, 0.0, 1.0));
    EXPECT_TRUE(kernel.valid);
    EXPECT_DOUBLE_EQ(kernel.distance, 0.5);
    EXPECT_TRUE(std::isfinite(kernel.scalar));
    EXPECT_NE(kernel.scalar, 0.0);

    const auto singular = ring_susceptibility_detail::ComputeKernel(
        Vec3(0.0, 0.0,
             CalculatorConfig::Get("singularity_guard_distance")),
        Vec3::Zero(), Vec3(0.0, 0.0, 1.0));
    EXPECT_FALSE(singular.valid);
}

TEST(RingChiAnalytical, ScalarMatchesRingContributionColumn7Formula) {
    if (!fs::exists(nmr::test::TestEnvironment::UbqProtonated())) GTEST_SKIP() << "1UBQ not found";
    auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
    if (!r.Ok()) GTEST_SKIP() << r.error;

    auto& conf = r.protein->Conformation();
    conf.AttachResult(GeometryResult::Compute(conf));
    conf.AttachResult(SpatialIndexResult::Compute(conf));
    conf.AttachResult(RingSusceptibilityResult::Compute(conf));

    ASSERT_GT(r.protein->RingCount(), 0u) << "1UBQ should have rings";

    int checked = 0;
    double max_diff = 0.0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& rn : conf.AtomAt(ai).ring_neighbours) {
            if (std::abs(rn.chi_scalar) < 1e-15) continue;

            const RingGeometry& geom = conf.ring_geometries[rn.ring_index];
            Vec3 d = conf.PositionAt(ai) - geom.center;
            double rm = d.norm();
            if (rm < MIN_DISTANCE) continue;

            Vec3 d_hat = d / rm;
            double cos_theta = d_hat.dot(geom.normal);
            double r3 = rm * rm * rm;

            double f = (3.0 * cos_theta * cos_theta - 1.0) / r3;
            double col7_formula = 0.0;
            if (rn.distance_to_center > 1e-12) {
                double cos_th = rn.z / rn.distance_to_center;
                double rr3 = rn.distance_to_center * rn.distance_to_center
                           * rn.distance_to_center;
                if (rr3 > 1e-30)
                    col7_formula = (3.0 * cos_th * cos_th - 1.0) / rr3;
            }
            double diff = std::abs(col7_formula - f);
            max_diff = std::max(max_diff, diff);

            EXPECT_NEAR(rn.chi_scalar, f, 1e-10);
            EXPECT_NEAR(col7_formula, f, 1e-10);
            checked++;
        }
    }

    EXPECT_GT(checked, 100) << "Should verify many ring-atom pairs";
    EXPECT_LT(max_diff, 1e-10) << "Column-7 formula must equal f at machine precision";

    std::cout << "  Verified ring-chi scalar rescue on " << checked
              << " ring-atom pairs, max |col7-f| = " << max_diff << "\n";
}


TEST(RingChiEmission, SparseAndPerTypeArraysReadBackInFinalRowOrder) {
    Protein protein;
    ProteinConformation conf(
        &protein, {Vec3::Zero()}, "ringchi-write-readback");
    auto& rows = conf.MutableAtomAt(0).ring_neighbours;
    rows.emplace_back();
    rows.back().ring_index = 4;
    rows.back().ring_type = RingTypeIndex::HisImidazole;
    rows.back().chi_scalar = -0.125;
    rows.emplace_back();
    rows.back().ring_index = 1;
    rows.back().ring_type = RingTypeIndex::PheBenzene;
    rows.back().chi_scalar = 0.75;

    RingSusceptibilityResult result;
    const fs::path output_dir =
        nmr::test::TestEnvironment::TempPath("ringchi_emit_readback");
    fs::create_directories(output_dir);
    ASSERT_EQ(result.WriteFeatures(conf, output_dir.string()), 2);

    const auto sparse = ReadFloat64Npy(output_dir / "ringchi_scalar.npy");
    const auto dense =
        ReadFloat64Npy(output_dir / "ringchi_per_type_T0.npy");
    EXPECT_EQ(sparse.shape, (std::vector<size_t>{2}));
    ASSERT_EQ(sparse.values.size(), 2u);
    EXPECT_DOUBLE_EQ(sparse.values[0], -0.125);
    EXPECT_DOUBLE_EQ(sparse.values[1], 0.75);
    EXPECT_EQ(dense.shape, (std::vector<size_t>{1, 8}));
    ASSERT_EQ(dense.values.size(), 8u);
    EXPECT_DOUBLE_EQ(
        dense.values[static_cast<int>(RingTypeIndex::HisImidazole)],
        -0.125);
    EXPECT_DOUBLE_EQ(
        dense.values[static_cast<int>(RingTypeIndex::PheBenzene)],
        0.75);

    ProteinConformation empty_conf(
        &protein, {Vec3::Zero()}, "ringchi-zero-row-readback");
    const fs::path empty_dir =
        nmr::test::TestEnvironment::TempPath("ringchi_emit_zero_rows");
    fs::create_directories(empty_dir);
    ASSERT_EQ(result.WriteFeatures(empty_conf, empty_dir.string()), 2);
    EXPECT_EQ(ReadFloat64Npy(empty_dir / "ringchi_scalar.npy").shape,
              (std::vector<size_t>{0}));
    EXPECT_EQ(ReadFloat64Npy(
                  empty_dir / "ringchi_per_type_T0.npy").shape,
              (std::vector<size_t>{1, 8}));
}


// ============================================================================
// Full protein test
// ============================================================================

class RingChiProteinTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!fs::exists(nmr::test::TestEnvironment::UbqProtonated())) {
            GTEST_SKIP() << "1UBQ.pdb not found";
        }
        auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
        if (!r.Ok()) GTEST_SKIP() << "Failed to load 1UBQ";
        protein = std::move(r.protein);

        auto& conf = protein->Conformation();
        conf.AttachResult(GeometryResult::Compute(conf));
        conf.AttachResult(SpatialIndexResult::Compute(conf));
    }

    std::unique_ptr<Protein> protein;
};


TEST_F(RingChiProteinTest, ComputeAndAttach) {
    auto& conf = protein->Conformation();
    auto rchi = RingSusceptibilityResult::Compute(conf);
    ASSERT_NE(rchi, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(rchi)));
    ASSERT_TRUE(conf.HasResult<RingSusceptibilityResult>());
}


TEST_F(RingChiProteinTest, ScalarMatchesColumn7Formula) {
    auto& conf = protein->Conformation();
    conf.AttachResult(RingSusceptibilityResult::Compute(conf));

    int checked = 0;
    double max_diff = 0.0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& rn : conf.AtomAt(ai).ring_neighbours) {
            if (std::abs(rn.chi_scalar) < 1e-15) continue;
            double col7_formula = 0.0;
            if (rn.distance_to_center > 1e-12) {
                double cos_th = rn.z / rn.distance_to_center;
                double r3 = rn.distance_to_center * rn.distance_to_center
                          * rn.distance_to_center;
                if (r3 > 1e-30)
                    col7_formula = (3.0 * cos_th * cos_th - 1.0) / r3;
            }
            double diff = std::abs(col7_formula - rn.chi_scalar);
            max_diff = std::max(max_diff, diff);
            checked++;
        }
    }

    EXPECT_GT(checked, 0) << "Should have checked some ring-atom pairs";
    EXPECT_LT(max_diff, 1e-10)
        << "Ring contribution column-7 formula must equal scalar f";

    std::cout << "  Checked " << checked << " pairs, max |col7 - f| = "
              << max_diff << "\n";
}


TEST_F(RingChiProteinTest, ScalarIsPopulated) {
    auto& conf = protein->Conformation();
    conf.AttachResult(RingSusceptibilityResult::Compute(conf));

    int nonzero_scalar = 0;
    double max_scalar = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& rn : conf.AtomAt(ai).ring_neighbours) {
            const double s = std::abs(rn.chi_scalar);
            if (s > 1e-8) nonzero_scalar++;
            max_scalar = std::max(max_scalar, s);
        }
    }

    EXPECT_GT(nonzero_scalar, 0) << "Ring chi scalar should be non-zero";

    std::cout << "  Scalar nonzero: " << nonzero_scalar
              << ", max |scalar| = " << max_scalar << "\n";
}


TEST_F(RingChiProteinTest, AzimuthIsPopulatedWithoutBiotSavart) {
    auto& conf = protein->Conformation();
    ASSERT_TRUE(conf.AttachResult(RingSusceptibilityResult::Compute(conf)));

    size_t checked = 0;
    size_t non_default = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const RingNeighbourhood& row : conf.AtomAt(ai).ring_neighbours) {
            if (row.rho <= CalculatorConfig::Get(
                    "near_zero_vector_norm_threshold")) continue;
            EXPECT_TRUE(std::isfinite(row.cos_phi));
            EXPECT_TRUE(std::isfinite(row.sin_phi));
            EXPECT_NEAR(row.cos_phi * row.cos_phi +
                            row.sin_phi * row.sin_phi,
                        1.0, 1e-12);
            if (std::abs(row.cos_phi - 1.0) > 1e-8 ||
                std::abs(row.sin_phi) > 1e-8) {
                ++non_default;
            }
            ++checked;
        }
    }
    EXPECT_GT(checked, 0u);
    EXPECT_GT(non_default, 0u);
}


// ============================================================================
// ORCA protein test
// ============================================================================

TEST(RingChiOrcaTest, RunOnProtonatedProtein) {
    OrcaRunFiles files;
    files.pdb_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.pdb";
    files.xyz_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.xyz";
    files.prmtop_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.prmtop";

    if (!fs::exists(files.xyz_path) || !fs::exists(files.prmtop_path))
        GTEST_SKIP() << "ORCA test data not found";

    auto load = BuildFromOrca(files);
    ASSERT_TRUE(load.Ok()) << load.error;

    auto& conf = load.protein->Conformation();
    conf.AttachResult(GeometryResult::Compute(conf));
    conf.AttachResult(SpatialIndexResult::Compute(conf));

    auto rchi = RingSusceptibilityResult::Compute(conf);
    ASSERT_NE(rchi, nullptr);
    conf.AttachResult(std::move(rchi));

    // Summary
    double max_scalar = 0;
    int with_rings = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& atom = conf.AtomAt(ai);
        if (!atom.ring_neighbours.empty()) with_rings++;
        for (const auto& rn : atom.ring_neighbours)
            max_scalar = std::max(max_scalar, std::abs(rn.chi_scalar));
    }

    std::cout << "  ORCA protein RingSusceptibility summary:\n"
              << "    atoms=" << conf.AtomCount()
              << " rings=" << load.protein->RingCount() << "\n"
              << "    atoms with ring neighbours: " << with_rings << "\n"
              << "    max |scalar| = " << max_scalar << " A^-3\n";

    EXPECT_GT(with_rings, 0) << "Some atoms should see rings";
    EXPECT_GT(max_scalar, 0.001) << "Scalar should be non-zero";
}

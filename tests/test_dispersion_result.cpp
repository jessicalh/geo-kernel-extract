#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <array>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include "DispersionResult.h"
#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "Protein.h"
#include "PdbFileReader.h"
#include "OrcaRunLoader.h"
#include "PhysicalConstants.h"

#include <filesystem>
namespace fs = std::filesystem;
using namespace nmr;

namespace {

struct NpyPayload {
    std::vector<size_t> shape;
    std::vector<char> payload;
};

std::string Trim(std::string s) {
    auto is_space = [](unsigned char c) { return std::isspace(c) != 0; };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                                    [&](char c) { return !is_space(c); }));
    s.erase(std::find_if(s.rbegin(), s.rend(),
                         [&](char c) { return !is_space(c); }).base(), s.end());
    return s;
}

NpyPayload ReadNpyPayload(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    EXPECT_TRUE(in.is_open()) << path;
    NpyPayload out;
    if (!in.is_open()) return out;

    char magic[6] = {};
    in.read(magic, 6);
    EXPECT_EQ(std::string(magic, 6), std::string("\x93NUMPY", 6));
    char version[2] = {};
    in.read(version, 2);
    EXPECT_EQ(version[0], 1);
    EXPECT_EQ(version[1], 0);
    uint16_t header_len = 0;
    in.read(reinterpret_cast<char*>(&header_len), sizeof(header_len));
    std::string header(header_len, '\0');
    in.read(header.data(), header_len);

    auto shape_begin = header.find('(');
    auto shape_end = header.find(')', shape_begin);
    EXPECT_NE(shape_begin, std::string::npos);
    EXPECT_NE(shape_end, std::string::npos);
    if (shape_begin == std::string::npos || shape_end == std::string::npos)
        return out;
    std::stringstream ss(header.substr(shape_begin + 1,
                                      shape_end - shape_begin - 1));
    std::string token;
    while (std::getline(ss, token, ',')) {
        token = Trim(token);
        if (!token.empty()) out.shape.push_back(static_cast<size_t>(std::stoull(token)));
    }

    out.payload.assign(std::istreambuf_iterator<char>(in),
                       std::istreambuf_iterator<char>());
    return out;
}

}  // namespace



// ============================================================================
// Analytical test: single vertex at known distance.
//
// Vertex at origin, atom at (3, 0, 0). r = 3.
//   K_xx = 3*9/r^8 - 1/r^6 = 27/6561 - 1/729 = 0.004115 - 0.001372 = 0.002743
//   K_yy = 0 - 1/r^6 = -1/729 = -0.001372
//   K_zz = 0 - 1/r^6 = -1/729 = -0.001372
//   Tr = K_xx + K_yy + K_zz = 0.002743 - 0.001372 - 0.001372 = ~0
//   scalar = 1/r^6 = 1/729
// ============================================================================

TEST(DispAnalytical, SingleVertexTraceless) {
    Vec3 d(3, 0, 0);
    double r = d.norm();
    double r2 = r * r;
    double r6 = r2 * r2 * r2;
    double r8 = r6 * r2;

    Mat3 K = Mat3::Zero();
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            K(a, b) = 3.0 * d(a) * d(b) / r8 - (a == b ? 1.0 : 0.0) / r6;

    EXPECT_NEAR(K.trace(), 0.0, 1e-14) << "Dispersion kernel is traceless";
    EXPECT_NEAR(K(0, 0), 2.0 / r6, 1e-14);   // 3*9/r^8 - 1/r^6 = 3/r^6 - 1/r^6 = 2/r^6
    EXPECT_NEAR(K(1, 1), -1.0 / r6, 1e-14);
    EXPECT_NEAR(K(2, 2), -1.0 / r6, 1e-14);

    double scalar = 1.0 / r6;
    EXPECT_NEAR(scalar, 1.0 / 729.0, 1e-14);

    std::cout << "  Single vertex at r=3: K_diag = ("
              << K(0,0) << ", " << K(1,1) << ", " << K(2,2) << ")\n"
              << "  Trace = " << K.trace() << "\n"
              << "  Scalar = " << scalar << "\n";
}


TEST(DispAnalytical, SwitchingFunctionProperties) {
    // Reproduce the CHARMM switching function locally for testing.
    // Must match DispersionResult.cpp: R_SWITCH=4.3A, R_CUT=5.0A.
    // Brooks et al., J. Comput. Chem. 4, 187 (1983).
    constexpr double R_SWITCH = 4.3;  // Angstroms — onset of taper
    constexpr double R_CUT = 5.0;     // Angstroms — zero beyond

    auto S = [](double r) -> double {
        if (r <= R_SWITCH) return 1.0;
        if (r >= R_CUT) return 0.0;
        double rc2 = R_CUT * R_CUT;
        double rs2 = R_SWITCH * R_SWITCH;
        double r2 = r * r;
        double num = (rc2 - r2) * (rc2 - r2) * (rc2 + 2.0 * r2 - 3.0 * rs2);
        double den = (rc2 - rs2) * (rc2 - rs2) * (rc2 - rs2);
        return num / den;
    };

    // Boundary values
    EXPECT_NEAR(S(2.0), 1.0, 1e-12) << "Below switch onset: S=1";
    EXPECT_NEAR(S(4.3), 1.0, 1e-12) << "At switch onset: S=1";
    EXPECT_NEAR(S(5.0), 0.0, 1e-12) << "At cutoff: S=0";
    EXPECT_GT(S(4.5), 0.0) << "In taper: S > 0";
    EXPECT_LT(S(4.5), 1.0) << "In taper: S < 1";

    // C¹ continuity: S'(R_switch) = 0 and S'(R_cut) = 0
    // Verify numerically via finite differences
    double h = 1e-7;
    double dS_at_switch = (S(R_SWITCH + h) - S(R_SWITCH - h)) / (2 * h);
    double dS_at_cut = (S(R_CUT - h) - S(R_CUT - 2*h)) / h;
    EXPECT_NEAR(dS_at_switch, 0.0, 1e-4) << "S'(R_switch) = 0 (C¹)";
    EXPECT_NEAR(dS_at_cut, 0.0, 1e-4) << "S'(R_cut) = 0 (C¹)";

    // Monotonic decrease through the taper
    double prev = 1.0;
    for (double r = 4.3; r <= 5.0; r += 0.01) {
        double s = S(r);
        EXPECT_LE(s, prev + 1e-12)
            << "S must decrease monotonically at r=" << r;
        prev = s;
    }

    // Smoothness: S(r)/r^6 should decrease monotonically for all r > MIN_DISTANCE
    double prev_weighted = 1e10;
    for (double r = 1.6; r < 5.0; r += 0.05) {
        double weighted = S(r) / std::pow(r, 6);
        EXPECT_LE(weighted, prev_weighted + 1e-12)
            << "S/r^6 must decrease monotonically at r=" << r;
        prev_weighted = weighted;
    }

    std::cout << "  Switching function (CHARMM form, Brooks 1983):\n"
              << "    onset=" << R_SWITCH << "A, cutoff=" << R_CUT << "A\n"
              << "    S(3.0)=" << S(3.0) << " S(4.5)=" << S(4.5)
              << " S(4.9)=" << S(4.9) << " S(5.0)=" << S(5.0) << "\n"
              << "    S'(onset)=" << dS_at_switch
              << " S'(cutoff)=" << dS_at_cut << " (both ~0 for C¹)\n";
}


TEST(DispAnalytical, VertexKernelAcceptsPointInsideRingCenterRadius) {
    RingGeometry geom;
    geom.center = Vec3::Zero();
    geom.normal = Vec3(0.0, 0.0, 1.0);
    geom.radius = 1.4;
    const double pi = std::acos(-1.0);
    for (int i = 0; i < 6; ++i) {
        const double angle = 2.0 * pi * static_cast<double>(i) / 6.0;
        geom.vertices.emplace_back(geom.radius * std::cos(angle),
                                   geom.radius * std::sin(angle), 0.0);
    }

    const Vec3 point(0.0, 0.0, 0.5);
    ASSERT_LT((point - geom.center).norm(), geom.radius);
    const auto vertices =
        dispersion_detail::ComputeRingVertices(point, geom);
    ASSERT_EQ(vertices.size(), geom.vertices.size());

    double scalar_sum = 0.0;
    for (const auto& vertex : vertices) {
        EXPECT_TRUE(vertex.valid)
            << "finite vertex distance, not center distance, owns validity";
        EXPECT_TRUE(std::isfinite(vertex.scalar));
        scalar_sum += vertex.scalar;
    }
    EXPECT_GT(scalar_sum, 0.0);
}


// ============================================================================
// Full protein test on 1UBQ
// ============================================================================

class DispProteinTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!fs::exists(nmr::test::TestEnvironment::UbqProtonated())) GTEST_SKIP() << "1UBQ not found";
        auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
        if (!r.Ok()) GTEST_SKIP() << "Failed to load 1UBQ";
        protein = std::move(r.protein);

        auto& conf = protein->Conformation();
        conf.AttachResult(GeometryResult::Compute(conf));
        conf.AttachResult(SpatialIndexResult::Compute(conf));
    }
    std::unique_ptr<Protein> protein;
};


TEST_F(DispProteinTest, ComputeAndAttach) {
    auto& conf = protein->Conformation();
    auto disp = DispersionResult::Compute(conf);
    ASSERT_NE(disp, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(disp)));
    ASSERT_TRUE(conf.HasResult<DispersionResult>());
}


TEST_F(DispProteinTest, ScalarsAndContactsPopulated) {
    auto& conf = protein->Conformation();
    conf.AttachResult(DispersionResult::Compute(conf));

    int with_scalar = 0;
    double max_scalar = 0.0;

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& rn : conf.AtomAt(ai).ring_neighbours) {
            if (rn.disp_scalar > 0.0) {
                with_scalar++;
                max_scalar = std::max(max_scalar, rn.disp_scalar);
            }
        }
    }

    EXPECT_GT(with_scalar, 0) << "Should have dispersion scalar pairs";
    EXPECT_GT(max_scalar, 0.0) << "Dispersion scalar should be positive";

    std::cout << "  Pairs with scalar: " << with_scalar
              << ", max scalar = " << max_scalar << "\n";
}


TEST_F(DispProteinTest, ContactCountsReasonable) {
    auto& conf = protein->Conformation();
    conf.AttachResult(DispersionResult::Compute(conf));

    int total_contacts = 0;
    int max_contacts = 0;
    int with_disp = 0;

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& rn : conf.AtomAt(ai).ring_neighbours) {
            if (rn.disp_contacts > 0) {
                total_contacts += rn.disp_contacts;
                max_contacts = std::max(max_contacts, rn.disp_contacts);
                with_disp++;
            }
        }
    }

    EXPECT_GT(with_disp, 0) << "Some atom-ring pairs should have contacts";
    // Each ring has 5-9 vertices. Max contacts per ring <= ring size.
    EXPECT_LE(max_contacts, 9) << "Max contacts should not exceed max ring size";

    std::cout << "  Pairs with contacts: " << with_disp
              << ", total contacts: " << total_contacts
              << ", max per pair: " << max_contacts << "\n";
}


TEST_F(DispProteinTest, PerTypeT0AccumulatesScalars) {
    auto& conf = protein->Conformation();
    conf.AttachResult(DispersionResult::Compute(conf));

    int checked_atoms = 0;

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& ca = conf.AtomAt(ai);
        std::array<double, 8> sparse_sum = {};
        for (const auto& rn : ca.ring_neighbours) {
            int ti = static_cast<int>(rn.ring_type);
            if (ti >= 0 && ti < 8)
                sparse_sum[ti] += rn.disp_scalar;
        }
        for (int ti = 0; ti < 8; ++ti)
            EXPECT_NEAR(ca.per_type_disp_scalar_sum[ti], sparse_sum[ti], 1e-14);
        checked_atoms++;
    }

    EXPECT_GT(checked_atoms, 0) << "Should check atoms";

    std::cout << "  Checked per-type dispersion scalar sums on "
              << checked_atoms << " atoms\n";
}


TEST_F(DispProteinTest, AromaticR6ProximityAliasIsPayloadIdentical) {
    auto& conf = protein->Conformation();
    auto disp = DispersionResult::Compute(conf);
    ASSERT_NE(disp, nullptr);
    const std::string out_dir =
        nmr::test::TestEnvironment::TempPath("dispersion_alias_out");
    fs::create_directories(out_dir);
    ASSERT_EQ(disp->WriteFeatures(conf, out_dir), 2);

    const fs::path legacy_path =
        fs::path(out_dir) / "disp_per_type_T0.npy";
    const fs::path alias_path =
        fs::path(out_dir) / "aromatic_r6_proximity_per_type_T0.npy";
    ASSERT_TRUE(fs::exists(legacy_path));
    ASSERT_TRUE(fs::exists(alias_path));

    const auto legacy = ReadNpyPayload(legacy_path);
    const auto alias = ReadNpyPayload(alias_path);
    const size_t N = conf.AtomCount();
    EXPECT_EQ(legacy.shape, (std::vector<size_t>{N, 8}));
    EXPECT_EQ(alias.shape, (std::vector<size_t>{N, 8}));
    ASSERT_EQ(legacy.payload.size(), alias.payload.size());

    if (legacy.payload != alias.payload) {
        const size_t n_double = legacy.payload.size() / sizeof(double);
        for (size_t i = 0; i < n_double; ++i) {
            double a = 0.0, b = 0.0;
            std::memcpy(&a, legacy.payload.data() + i * sizeof(double),
                        sizeof(double));
            std::memcpy(&b, alias.payload.data() + i * sizeof(double),
                        sizeof(double));
            if (a != b) {
                ADD_FAILURE()
                    << "dispersion alias payload differs at element " << i
                    << ": legacy=" << a << " alias=" << b;
                break;
            }
        }
    }
    EXPECT_EQ(legacy.payload, alias.payload);
}


// ============================================================================
// ORCA protein test
// ============================================================================

TEST(DispOrcaTest, RunOnProtonatedProtein) {
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
    conf.AttachResult(DispersionResult::Compute(conf));

    int with_disp = 0;
    double max_scalar = 0.0;

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& rn : conf.AtomAt(ai).ring_neighbours) {
            if (rn.disp_contacts > 0) with_disp++;
            max_scalar = std::max(max_scalar, std::abs(rn.disp_scalar));
        }
    }

    EXPECT_GT(with_disp, 0) << "Some pairs should have dispersion contacts";

    std::cout << "  ORCA protein Dispersion summary:\n"
              << "    atoms=" << conf.AtomCount()
              << " rings=" << load.protein->RingCount() << "\n"
              << "    pairs with contacts: " << with_disp << "\n"
              << "    max |scalar| = " << max_scalar << " A^-6\n";
}

#include "CategoryInfoProjection.h"
#include "ConformationResult.h"
#include "DsspResult.h"
#include "FrameNpyEmitter.h"
#include "GeometryResult.h"
#include "GromacsFrameHandler.h"
#include "GromacsFramePullResult.h"
#include "HBondResult.h"
#include "HydrationGeometryResult.h"
#include "HydrationShellResult.h"
#include "LocalBackboneGeometryResult.h"
#include "LocalBackboneGeometryTrajectoryResult.h"
#include "OperationRunner.h"
#include "PhysicalConstants.h"
#include "PdbFileReader.h"
#include "PlanarGeometryResult.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "RunConfiguration.h"
#include "SasaResult.h"
#include "Session.h"
#include "SpatialIndexResult.h"
#include "TestEnvironment.h"
#include "TopologySidecar.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"
#include "WaterFieldResult.h"
#include "WaterHBondGeometryResult.h"

#include <gtest/gtest.h>

#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <Eigen/Geometry>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>

// Keep GROMACS last: its legacy headers define DIM as a macro.
#include "gromacs/fileio/trrio.h"
#include "gromacs/utility/vectypes.h"

#ifndef NMR_TEST_DATA_DIR
#error "NMR_TEST_DATA_DIR must be defined"
#endif
#ifndef NMR_TEST_PYTHON_EXECUTABLE
#error "NMR_TEST_PYTHON_EXECUTABLE must be defined"
#endif
#ifndef NMR_LOCAL_GEOMETRY_ACCEPTANCE_SCRIPT
#error "NMR_LOCAL_GEOMETRY_ACCEPTANCE_SCRIPT must be defined"
#endif
#ifndef NMR_LOCAL_GEOMETRY_TRAJECTORY_FIXTURE
#error "NMR_LOCAL_GEOMETRY_TRAJECTORY_FIXTURE must be defined"
#endif

namespace fs = std::filesystem;

namespace {

struct NpyArray {
    std::string descr;
    std::vector<std::size_t> shape;
    std::vector<unsigned char> bytes;
};

std::string Trim(std::string text) {
    const auto first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, last - first + 1);
}

NpyArray ReadNpy(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    EXPECT_TRUE(in.is_open()) << path;
    NpyArray out;
    if (!in) return out;

    char magic[6]{};
    in.read(magic, 6);
    EXPECT_EQ(std::string(magic, 6), std::string("\x93NUMPY", 6)) << path;
    std::uint8_t major = 0, minor = 0;
    in.read(reinterpret_cast<char*>(&major), 1);
    in.read(reinterpret_cast<char*>(&minor), 1);
    EXPECT_EQ(major, 1u) << path;
    EXPECT_EQ(minor, 0u) << path;
    std::uint16_t header_length = 0;
    in.read(reinterpret_cast<char*>(&header_length), sizeof(header_length));
    std::string header(header_length, '\0');
    in.read(header.data(), static_cast<std::streamsize>(header.size()));

    const std::string descr_key = "'descr': '";
    auto descr_start = header.find(descr_key);
    EXPECT_NE(descr_start, std::string::npos) << path;
    if (descr_start == std::string::npos) return out;
    descr_start += descr_key.size();
    const auto descr_end = header.find('\'', descr_start);
    out.descr = header.substr(descr_start, descr_end - descr_start);

    const auto shape_key = header.find("'shape':");
    const auto shape_open = header.find('(', shape_key);
    const auto shape_close = header.find(')', shape_open);
    EXPECT_NE(shape_key, std::string::npos) << path;
    EXPECT_NE(shape_open, std::string::npos) << path;
    EXPECT_NE(shape_close, std::string::npos) << path;
    std::stringstream shape_stream(
        header.substr(shape_open + 1, shape_close - shape_open - 1));
    std::string token;
    while (std::getline(shape_stream, token, ',')) {
        token = Trim(token);
        if (!token.empty()) out.shape.push_back(std::stoull(token));
    }

    in.seekg(0, std::ios::end);
    const auto end = in.tellg();
    const auto data_start = static_cast<std::streamoff>(10 + header_length);
    EXPECT_GE(end, data_start) << path;
    if (end < data_start) return out;
    out.bytes.resize(static_cast<std::size_t>(end - data_start));
    in.seekg(data_start);
    in.read(reinterpret_cast<char*>(out.bytes.data()),
            static_cast<std::streamsize>(out.bytes.size()));
    EXPECT_TRUE(in.good() || in.eof()) << path;
    return out;
}

template <typename T>
std::vector<T> Values(const NpyArray& array) {
    EXPECT_EQ(array.bytes.size() % sizeof(T), 0u);
    std::vector<T> values(array.bytes.size() / sizeof(T));
    if (!values.empty())
        std::memcpy(values.data(), array.bytes.data(), array.bytes.size());
    return values;
}

void ExpectShape(const NpyArray& array,
                 std::initializer_list<std::size_t> expected,
                 const char* descr) {
    EXPECT_EQ(array.shape, std::vector<std::size_t>(expected));
    EXPECT_EQ(array.descr, descr);
}

void CheckHBondTypedRows(const nmr::Protein& protein,
                         const fs::path& output_dir,
                         std::size_t expected_count = SIZE_MAX) {
    const NpyArray index_file = ReadNpy(output_dir / "hbond_pairs_index.npy");
    EXPECT_EQ(index_file.descr, "<i4");
    EXPECT_EQ(index_file.shape.size(), 2u);
    if (index_file.shape.size() != 2u) return;
    EXPECT_EQ(index_file.shape[1], 6u);
    if (index_file.shape[1] != 6u) return;
    if (expected_count != SIZE_MAX)
        EXPECT_EQ(index_file.shape[0], expected_count);

    const auto rows = Values<std::int32_t>(index_file);
    EXPECT_EQ(rows.size(), index_file.shape[0] * 6);
    if (rows.size() != index_file.shape[0] * 6) return;
    for (std::size_t hi = 0; hi < index_file.shape[0]; ++hi) {
        const std::int32_t donor_ri = rows[hi * 6 + 0];
        const std::int32_t donor_n = rows[hi * 6 + 1];
        const std::int32_t donor_h = rows[hi * 6 + 2];
        const std::int32_t acceptor_ri = rows[hi * 6 + 3];
        const std::int32_t acceptor_o = rows[hi * 6 + 4];
        EXPECT_GE(donor_ri, 0) << "hbond_row=" << hi;
        EXPECT_GE(acceptor_ri, 0) << "hbond_row=" << hi;
        if (donor_ri < 0 || acceptor_ri < 0 ||
            static_cast<std::size_t>(donor_ri) >= protein.ResidueCount() ||
            static_cast<std::size_t>(acceptor_ri) >= protein.ResidueCount()) {
            continue;
        }
        const nmr::Residue& donor = protein.ResidueAt(
            static_cast<std::size_t>(donor_ri));
        const nmr::Residue& acceptor = protein.ResidueAt(
            static_cast<std::size_t>(acceptor_ri));
        EXPECT_EQ(donor_n, static_cast<std::int32_t>(donor.N))
            << "hbond_row=" << hi;
        EXPECT_EQ(donor_h, static_cast<std::int32_t>(donor.H))
            << "hbond_row=" << hi;
        EXPECT_EQ(acceptor_o, static_cast<std::int32_t>(acceptor.O))
            << "hbond_row=" << hi;
    }
}

double IndependentDihedral(const nmr::Vec3& p0, const nmr::Vec3& p1,
                           const nmr::Vec3& p2, const nmr::Vec3& p3) {
    const nmr::Vec3 central = p2 - p1;
    const double central_norm = central.norm();
    if (!p0.allFinite() || !p1.allFinite() ||
        !p2.allFinite() || !p3.allFinite() || central_norm <= 1e-12) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const nmr::Vec3 axis = central / central_norm;
    const nmr::Vec3 first = p0 - p1;
    const nmr::Vec3 last = p3 - p2;
    const nmr::Vec3 v = first - first.dot(axis) * axis;
    const nmr::Vec3 w = last - last.dot(axis) * axis;
    if (v.norm() <= 1e-12 || w.norm() <= 1e-12)
        return std::numeric_limits<double>::quiet_NaN();
    return std::atan2(v.cross(axis).dot(w), v.dot(w));
}

void ExpectNearOrBothNan(double actual, double expected, double tolerance) {
    if (std::isnan(expected)) {
        EXPECT_TRUE(std::isnan(actual));
    } else {
        EXPECT_NEAR(actual, expected, tolerance);
    }
}

std::string ShellQuote(const fs::path& path) {
    std::string quoted = "'";
    for (char ch : path.string()) {
        if (ch == '\'') quoted += "'\\''";
        else quoted += ch;
    }
    quoted += "'";
    return quoted;
}

fs::path FreshRoot() {
    const auto stamp = std::chrono::steady_clock::now()
                           .time_since_epoch().count();
    const fs::path root = fs::temp_directory_path() /
        ("npy_local_geometry_acceptance_" + std::to_string(::getpid()) +
         "_" + std::to_string(stamp));
    fs::create_directories(root);
    return root;
}

nmr::Mat3 DirectionalAcceptanceRotation() {
    // Fixed independently chosen proper transform. The Python serialized-
    // boundary oracle constructs the same Rodrigues rotation from these
    // recorded constants.
    return Eigen::AngleAxisd(
        0.731,
        nmr::Vec3(0.37, -0.61, 0.704).normalized()).toRotationMatrix();
}

std::vector<nmr::Vec3> DirectionalAcceptanceTransform(
        const nmr::ProteinConformation& conf) {
    const nmr::Mat3 rotation = DirectionalAcceptanceRotation();
    const nmr::Vec3 translation(3.25, -1.75, 2.5);
    std::vector<nmr::Vec3> transformed;
    transformed.reserve(conf.AtomCount());
    for (const nmr::Vec3& position : conf.Positions())
        transformed.push_back(rotation * position + translation);
    return transformed;
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

void ShiftTriclinicImage(std::vector<float>& xyz,
                         std::size_t atom_index,
                         const rvec box[3],
                         int na, int nb, int nc) {
    for (int d = 0; d < 3; ++d) {
        xyz[atom_index * 3 + static_cast<std::size_t>(d)] +=
            static_cast<float>(na) * box[0][d] +
            static_cast<float>(nb) * box[1][d] +
            static_cast<float>(nc) * box[2][d];
    }
}

void AttachExplicitSolventOwners(
        nmr::ProteinConformation& conf,
        const nmr::SolventEnvironment& solvent) {
    ASSERT_TRUE(conf.AttachResult(nmr::GeometryResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(nmr::SpatialIndexResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(nmr::SasaResult::Compute(conf)));

    auto field = nmr::WaterFieldResult::Compute(conf, solvent);
    ASSERT_NE(field, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(field)));

    auto hbond = nmr::WaterHBondGeometryResult::Compute(conf, solvent);
    ASSERT_NE(hbond, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(hbond)));

    auto shell = nmr::HydrationShellResult::Compute(conf, solvent);
    ASSERT_NE(shell, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(shell)));

    auto hydration = nmr::HydrationGeometryResult::Compute(conf, solvent);
    ASSERT_NE(hydration, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(hydration)));
}

void WriteExplicitSolventOwners(const nmr::ProteinConformation& conf,
                                const fs::path& output_dir) {
    fs::create_directories(output_dir);
    ASSERT_EQ(conf.Result<nmr::WaterFieldResult>().WriteFeatures(
                  conf, output_dir.string()), 9);
    ASSERT_EQ(conf.Result<nmr::WaterHBondGeometryResult>().WriteFeatures(
                  conf, output_dir.string()), 3);
    ASSERT_EQ(conf.Result<nmr::HydrationShellResult>().WriteFeatures(
                  conf, output_dir.string()), 1);
    ASSERT_EQ(conf.Result<nmr::HydrationGeometryResult>().WriteFeatures(
                  conf, output_dir.string()), 1);
}

void ExpectPeriodicNpyEqual(const fs::path& baseline_dir,
                            const fs::path& shifted_dir,
                            const char* filename,
                            double absolute_tolerance,
                            double relative_tolerance) {
    const NpyArray baseline = ReadNpy(baseline_dir / filename);
    const NpyArray shifted = ReadNpy(shifted_dir / filename);
    ASSERT_EQ(shifted.descr, baseline.descr) << filename;
    ASSERT_EQ(shifted.shape, baseline.shape) << filename;
    ASSERT_EQ(shifted.bytes.size(), baseline.bytes.size()) << filename;

    if (baseline.descr != "<f8") {
        EXPECT_EQ(shifted.bytes, baseline.bytes) << filename;
        std::cout << "PBC_INVARIANCE " << filename << " exact_bytes\n";
        return;
    }

    const auto a = Values<double>(baseline);
    const auto b = Values<double>(shifted);
    ASSERT_EQ(b.size(), a.size()) << filename;
    double max_absolute_delta = 0.0;
    double max_scaled_delta = 0.0;
    std::size_t max_absolute_index = 0;
    std::size_t mismatch_count = 0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (std::isnan(a[i]) || std::isnan(b[i])) {
            if (!(std::isnan(a[i]) && std::isnan(b[i]))) ++mismatch_count;
            continue;
        }
        if (std::isinf(a[i]) || std::isinf(b[i])) {
            if (a[i] != b[i]) ++mismatch_count;
            continue;
        }
        const double delta = std::abs(a[i] - b[i]);
        const double scale = std::max(std::abs(a[i]), std::abs(b[i]));
        if (delta > max_absolute_delta) {
            max_absolute_delta = delta;
            max_absolute_index = i;
        }
        max_scaled_delta = std::max(
            max_scaled_delta,
            delta / std::max(1.0, scale));
        if (delta > absolute_tolerance + relative_tolerance * scale)
            ++mismatch_count;
    }
    EXPECT_EQ(mismatch_count, 0u)
        << filename << " max_abs=" << max_absolute_delta
        << " max_scaled=" << max_scaled_delta
        << " max_index=" << max_absolute_index
        << " baseline=" << a[max_absolute_index]
        << " shifted=" << b[max_absolute_index];
    std::cout << "PBC_INVARIANCE " << filename
              << " max_abs=" << max_absolute_delta
              << " mismatches=" << mismatch_count << '\n';
}

void ExpectExactFloatColumns(const fs::path& baseline_dir,
                             const fs::path& shifted_dir,
                             const char* filename,
                             std::initializer_list<std::size_t> columns) {
    const NpyArray baseline = ReadNpy(baseline_dir / filename);
    const NpyArray shifted = ReadNpy(shifted_dir / filename);
    ASSERT_EQ(baseline.descr, "<f8") << filename;
    ASSERT_EQ(shifted.descr, baseline.descr) << filename;
    ASSERT_EQ(shifted.shape, baseline.shape) << filename;
    ASSERT_EQ(baseline.shape.size(), 2u) << filename;
    const auto a = Values<double>(baseline);
    const auto b = Values<double>(shifted);
    const std::size_t width = baseline.shape[1];
    std::size_t mismatches = 0;
    for (std::size_t row = 0; row < baseline.shape[0]; ++row) {
        for (const std::size_t column : columns) {
            ASSERT_LT(column, width) << filename;
            mismatches += (a[row * width + column] !=
                           b[row * width + column]);
        }
    }
    EXPECT_EQ(mismatches, 0u) << filename;
    std::cout << "PBC_INVARIANCE " << filename
              << " identity_columns_exact mismatches=" << mismatches << '\n';
}

void CheckCoordinateTorsionOracles(
        const nmr::Protein& protein,
        const fs::path& output_dir,
        const std::array<std::size_t, 7>* exact_counts = nullptr) {
    const std::size_t R = protein.ResidueCount();
    const std::size_t N = protein.AtomCount();
    const NpyArray pos_file = ReadNpy(output_dir / "pos.npy");
    const NpyArray observed_file = ReadNpy(output_dir / "dssp_observed.npy");
    const NpyArray angle_file = ReadNpy(output_dir / "dssp_torsion_angle.npy");
    const NpyArray sin_file = ReadNpy(output_dir / "dssp_torsion_sin.npy");
    const NpyArray cos_file = ReadNpy(output_dir / "dssp_torsion_cos.npy");
    const NpyArray valid_file = ReadNpy(output_dir / "dssp_torsion_valid.npy");
    ExpectShape(pos_file, {N, 3}, "<f8");
    ExpectShape(observed_file, {N}, "|i1");
    ExpectShape(angle_file, {R, 6}, "<f8");
    ExpectShape(sin_file, {R, 6}, "<f8");
    ExpectShape(cos_file, {R, 6}, "<f8");
    ExpectShape(valid_file, {R, 6}, "|u1");
    const auto pos = Values<double>(pos_file);
    const auto observed = Values<std::int8_t>(observed_file);
    const auto angle = Values<double>(angle_file);
    const auto sine = Values<double>(sin_file);
    const auto cosine = Values<double>(cos_file);
    const auto valid = Values<std::uint8_t>(valid_file);
    if (pos.size() != N * 3 || observed.size() != N ||
        angle.size() != R * 6 || sine.size() != R * 6 ||
        cosine.size() != R * 6 || valid.size() != R * 6) return;

    auto position = [&](std::size_t atom_index) {
        return nmr::Vec3(pos[atom_index * 3 + 0],
                         pos[atom_index * 3 + 1],
                         pos[atom_index * 3 + 2]);
    };
    auto residue_observed = [&](const nmr::Residue& residue) {
        return !residue.atom_indices.empty() &&
            observed[residue.atom_indices.front()] != 0;
    };

    std::array<std::size_t, 6> checked{0, 0, 0, 0, 0, 0};
    for (std::size_t ri = 0; ri < R; ++ri) {
        const nmr::Residue& residue = protein.ResidueAt(ri);
        const auto predecessor = protein.BackbonePredecessor(ri);
        const auto successor = protein.BackboneSuccessor(ri);
        std::array<double, 6> expected;
        expected.fill(std::numeric_limits<double>::quiet_NaN());
        if (predecessor && residue_observed(residue)) {
            const nmr::Residue& prev = protein.ResidueAt(*predecessor);
            expected[0] = IndependentDihedral(
                position(prev.C), position(residue.N),
                position(residue.CA), position(residue.C));
        }
        if (successor && residue_observed(residue)) {
            const nmr::Residue& next = protein.ResidueAt(*successor);
            expected[1] = IndependentDihedral(
                position(residue.N), position(residue.CA),
                position(residue.C), position(next.N));
        }
        for (std::size_t chi = 0; chi < 4; ++chi) {
            if (!residue.chi[chi].Valid()) continue;
            const auto& a = residue.chi[chi].a;
            expected[2 + chi] = IndependentDihedral(
                position(a[0]), position(a[1]), position(a[2]), position(a[3]));
        }

        for (std::size_t column = 0; column < 6; ++column) {
            const std::size_t offset = ri * 6 + column;
            const bool expected_valid = std::isfinite(expected[column]);
            EXPECT_EQ(valid[offset], expected_valid ? 1u : 0u)
                << "residue=" << ri << " torsion_column=" << column;
            if (!expected_valid) {
                EXPECT_TRUE(std::isnan(angle[offset]));
                EXPECT_TRUE(std::isnan(sine[offset]));
                EXPECT_TRUE(std::isnan(cosine[offset]));
                continue;
            }
            const double tolerance = column < 2 ? 0.05 : 1e-12;
            EXPECT_NEAR(angle[offset], expected[column], tolerance)
                << "residue=" << ri << " torsion_column=" << column;
            EXPECT_NEAR(sine[offset], std::sin(expected[column]), tolerance);
            EXPECT_NEAR(cosine[offset], std::cos(expected[column]), tolerance);
            ++checked[column];
        }
    }

    const NpyArray omega_file = ReadNpy(output_dir / "omega_actual.npy");
    const NpyArray omega_sin_file = ReadNpy(output_dir / "omega_sin.npy");
    const NpyArray omega_cos_file = ReadNpy(output_dir / "omega_cos.npy");
    const NpyArray omega_valid_file = ReadNpy(output_dir / "omega_valid.npy");
    ExpectShape(omega_file, {R}, "<f8");
    ExpectShape(omega_sin_file, {R}, "<f8");
    ExpectShape(omega_cos_file, {R}, "<f8");
    ExpectShape(omega_valid_file, {R}, "|u1");
    const auto omega = Values<double>(omega_file);
    const auto omega_sin = Values<double>(omega_sin_file);
    const auto omega_cos = Values<double>(omega_cos_file);
    const auto omega_valid = Values<std::uint8_t>(omega_valid_file);
    if (omega.size() != R || omega_sin.size() != R ||
        omega_cos.size() != R || omega_valid.size() != R) return;
    std::size_t omega_checked = 0;
    for (std::size_t ri = 0; ri < R; ++ri) {
        const auto successor = protein.BackboneSuccessor(ri);
        double expected = std::numeric_limits<double>::quiet_NaN();
        if (successor) {
            const nmr::Residue& residue = protein.ResidueAt(ri);
            const nmr::Residue& next = protein.ResidueAt(*successor);
            if (residue.CA != nmr::Residue::NONE &&
                next.CA != nmr::Residue::NONE) {
                expected = IndependentDihedral(
                    position(residue.CA), position(residue.C),
                    position(next.N), position(next.CA));
            }
        }
        const bool expected_valid = std::isfinite(expected);
        EXPECT_EQ(omega_valid[ri], expected_valid ? 1u : 0u);
        ExpectNearOrBothNan(omega[ri], expected, 1e-12);
        if (expected_valid) {
            EXPECT_NEAR(omega_sin[ri], std::sin(expected), 1e-12);
            EXPECT_NEAR(omega_cos[ri], std::cos(expected), 1e-12);
            ++omega_checked;
        } else {
            EXPECT_TRUE(std::isnan(omega_sin[ri]));
            EXPECT_TRUE(std::isnan(omega_cos[ri]));
        }
    }

    if (exact_counts) {
        for (std::size_t column = 0; column < 6; ++column)
            EXPECT_EQ(checked[column], (*exact_counts)[column]);
        EXPECT_EQ(omega_checked, (*exact_counts)[6]);
    } else {
        EXPECT_GT(checked[0], 0u);
        EXPECT_GT(checked[1], 0u);
        EXPECT_GT(omega_checked, 0u);
    }
}

void CheckStaticTorsionAndPartnerOracles(
        const nmr::ProteinConformation& conf,
        const fs::path& output_dir) {
    const nmr::Protein& protein = conf.ProteinRef();
    const std::size_t R = protein.ResidueCount();
    const std::size_t N = conf.AtomCount();
    const auto& dssp = conf.Result<nmr::DsspResult>();

    const std::array<std::size_t, 7> expected_counts{
        75, 75, 68, 54, 24, 11, 75
    };
    CheckCoordinateTorsionOracles(protein, output_dir, &expected_counts);

    const NpyArray partner_file = ReadNpy(
        output_dir / "dssp_hbond_partner_residue_index.npy");
    const NpyArray energy_file = ReadNpy(output_dir / "dssp_hbond_energy.npy");
    ExpectShape(partner_file, {N, 4}, "<i4");
    ExpectShape(energy_file, {N, 4}, "<f8");
    const auto partner = Values<std::int32_t>(partner_file);
    const auto energy = Values<double>(energy_file);
    std::size_t present[4] = {0, 0, 0, 0};
    std::size_t absent[4] = {0, 0, 0, 0};
    for (std::size_t ai = 0; ai < N; ++ai) {
        const std::size_t ri = protein.AtomAt(ai).residue_index;
        const nmr::DsspResidue& dr = dssp.AllResidues()[ri];
        const nmr::DsspResidue::HBondPartner* slots[4] = {
            &dr.acceptors[0], &dr.acceptors[1], &dr.donors[0], &dr.donors[1]
        };
        for (std::size_t slot = 0; slot < 4; ++slot) {
            const std::size_t offset = ai * 4 + slot;
            const bool has_partner = dr.observed &&
                slots[slot]->residue_index != SIZE_MAX &&
                slots[slot]->residue_index < R;
            EXPECT_EQ(partner[offset], has_partner
                ? static_cast<std::int32_t>(slots[slot]->residue_index) : -1);
            if (dr.observed) {
                EXPECT_DOUBLE_EQ(energy[offset], slots[slot]->energy);
            } else {
                EXPECT_TRUE(std::isnan(energy[offset]));
            }
            if (has_partner) ++present[slot];
            else ++absent[slot];
        }
    }
    for (std::size_t slot = 0; slot < 4; ++slot) {
        EXPECT_GT(present[slot], 0u) << "slot=" << slot;
    }
    EXPECT_GT(absent[0] + absent[1] + absent[2] + absent[3], 0u);

}

void CheckLocalBackboneAgainstTrajectoryH5(
        const nmr::TrajectoryProtein& tp,
        const nmr::LocalBackboneGeometryTrajectoryResult& result,
        const fs::path& frame_dir,
        const fs::path& h5_path) {
    {
        HighFive::File file(h5_path.string(), HighFive::File::Truncate);
        result.WriteH5Group(tp, file);
    }
    HighFive::File file(h5_path.string(), HighFive::File::ReadOnly);
    auto group = file.getGroup("/trajectory/local_backbone_geometry");
    const std::size_t R = tp.ProteinRef().ResidueCount();
    ASSERT_EQ(result.NumFrames(), 1u);

    for (const char* stem : {"tau_N_CA_C", "angle_N_CA_CB", "angle_CB_CA_C",
                              "angle_Cprev_N_CA", "angle_CA_C_Nnext",
                              "cb_deviation"}) {
        std::vector<double> h5(R);
        group.getDataSet(stem).read(h5.data());
        const NpyArray npy_file = ReadNpy(frame_dir / (std::string(stem) + ".npy"));
        const NpyArray valid_file = ReadNpy(
            frame_dir / (std::string(stem) + "_valid.npy"));
        ExpectShape(npy_file, {R}, "<f8");
        ExpectShape(valid_file, {R}, "|u1");
        const auto npy = Values<double>(npy_file);
        const auto valid = Values<std::uint8_t>(valid_file);
        for (std::size_t ri = 0; ri < R; ++ri) {
            ExpectNearOrBothNan(npy[ri], h5[ri], 0.0);
            EXPECT_EQ(valid[ri], std::isfinite(h5[ri]) ? 1u : 0u);
        }
    }

    std::vector<double> h5_vector(R * 3);
    group.getDataSet("cb_local_vector").read(h5_vector.data());
    const NpyArray vector_file = ReadNpy(frame_dir / "cb_residual_vector.npy");
    const NpyArray vector_valid_file = ReadNpy(
        frame_dir / "cb_residual_vector_valid.npy");
    ExpectShape(vector_file, {R, 3}, "<f8");
    ExpectShape(vector_valid_file, {R}, "|u1");
    const auto vector = Values<double>(vector_file);
    const auto valid = Values<std::uint8_t>(vector_valid_file);
    for (std::size_t ri = 0; ri < R; ++ri) {
        bool all_finite = true;
        for (std::size_t component = 0; component < 3; ++component) {
            const std::size_t offset = ri * 3 + component;
            ExpectNearOrBothNan(vector[offset], h5_vector[offset], 0.0);
            all_finite = all_finite && std::isfinite(h5_vector[offset]);
        }
        EXPECT_EQ(valid[ri], all_finite ? 1u : 0u);
    }
}

}  // namespace


TEST(NpyLocalGeometryAcceptance, StaticAndTrajectoryProductionPaths) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::FrameNpyEmitter::Reset();
    const fs::path root = FreshRoot();

    // Ordinary static production loader -> OperationRunner -> owner emitters.
    auto build = nmr::BuildFromProtonatedPdb(
        (fs::path(NMR_TEST_DATA_DIR) / "1ubq_protonated.pdb").string());
    ASSERT_TRUE(build.Ok()) << build.error;
    auto& static_conf = build.protein->Conformation();
    nmr::RunOptions static_options;
    static_options.skip_mopac = true;
    static_options.skip_apbs = true;
    static_options.skip_coulomb = true;
    static_options.skip_dssp = false;
    const auto static_run = nmr::OperationRunner::Run(static_conf, static_options);
    ASSERT_TRUE(static_run.Ok()) << static_run.error;
    ASSERT_TRUE(static_conf.HasResult<nmr::GeometryResult>());
    ASSERT_TRUE(static_conf.HasResult<nmr::LocalBackboneGeometryResult>());
    ASSERT_TRUE(static_conf.HasResult<nmr::DsspResult>());
    ASSERT_TRUE(static_conf.HasResult<nmr::PlanarGeometryResult>());
    ASSERT_TRUE(static_conf.HasResult<nmr::HBondResult>());

    const fs::path static_dir = root / "static";
    fs::create_directories(static_dir);
    ASSERT_EQ(nmr::CategoryInfoProjection::WriteFeatures(
                  static_conf.ProteinRef(), static_dir.string()), 1);
    ASSERT_EQ(nmr::TopologySidecar::WriteFeatures(
                  static_conf.ProteinRef(), static_dir.string()), 5);
    ASSERT_GT(nmr::ConformationResult::WriteAllFeatures(
                  static_conf, static_dir.string()), 0);
    CheckHBondTypedRows(
        static_conf.ProteinRef(), static_dir,
        static_conf.Result<nmr::HBondResult>().HBondCount());
    CheckStaticTorsionAndPartnerOracles(static_conf, static_dir);

    // Rerun the complete reduced production stack on transformed INPUT.
    // This is deliberately not an algebra-only rotation of emitted values:
    // every owning calculator reconstructs its output from Qx+t.
    nmr::ProteinConformation transformed_conf(
        &static_conf.ProteinRef(),
        DirectionalAcceptanceTransform(static_conf),
        "directional acceptance proper rotation + translation");
    const auto transformed_run = nmr::OperationRunner::Run(
        transformed_conf, static_options);
    ASSERT_TRUE(transformed_run.Ok()) << transformed_run.error;
    const fs::path transformed_dir = root / "static_transformed";
    fs::create_directories(transformed_dir);
    ASSERT_EQ(nmr::CategoryInfoProjection::WriteFeatures(
                  transformed_conf.ProteinRef(), transformed_dir.string()), 1);
    ASSERT_EQ(nmr::TopologySidecar::WriteFeatures(
                  transformed_conf.ProteinRef(), transformed_dir.string()), 5);
    ASSERT_GT(nmr::ConformationResult::WriteAllFeatures(
                  transformed_conf, transformed_dir.string()), 0);

    // Ordinary trajectory production loader -> Trajectory::Run ->
    // OperationRunner -> FrameNpyEmitter.  The committed fixture is one
    // complete frame, so there is no skip or cadence shortcut.
    const fs::path fixture_root = NMR_LOCAL_GEOMETRY_TRAJECTORY_FIXTURE;
    const fs::path production_dir = fixture_root / "production";
    const fs::path tpr_path = production_dir / "production.tpr";
    const fs::path trr_path = production_dir / "production.trr";
    ASSERT_TRUE(fs::is_regular_file(tpr_path));
    ASSERT_TRUE(fs::is_regular_file(trr_path));
    ASSERT_TRUE(fs::is_regular_file(fixture_root / "topol.top"));

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(production_dir.string())) << tp.Error();
    const fs::path trajectory_dir = root / "trajectory";
    nmr::FrameNpyEmitter::Config emitter_config;
    emitter_config.output_dir = trajectory_dir;
    nmr::FrameNpyEmitter::Configure(tp.ProteinRef(), emitter_config);

    nmr::RunConfiguration trajectory_config;
    trajectory_config.SetName("NpyLocalGeometryAcceptance");
    auto& frame_options = trajectory_config.MutablePerFrameRunOptions();
    frame_options.skip_mopac = true;
    frame_options.skip_apbs = true;
    frame_options.skip_coulomb = true;
    frame_options.skip_dssp = false;
    trajectory_config.RequireConformationResult(typeid(nmr::GeometryResult));
    trajectory_config.RequireConformationResult(
        typeid(nmr::LocalBackboneGeometryResult));
    trajectory_config.RequireConformationResult(typeid(nmr::DsspResult));
    trajectory_config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& protein)
            -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::LocalBackboneGeometryTrajectoryResult::Create(protein);
        });

    nmr::Trajectory trajectory(trr_path, tpr_path, {});
    nmr::Session session;
    ASSERT_EQ(trajectory.Run(tp, trajectory_config, session, {}, trajectory_dir),
              nmr::kOk);
    ASSERT_EQ(trajectory.FrameCount(), 1u);
    ASSERT_EQ(trajectory.FrameIndices().size(), 1u);
    const std::size_t frame_index = trajectory.FrameIndices().front();
    std::ostringstream frame_name;
    frame_name << "frame_";
    frame_name.width(6);
    frame_name.fill('0');
    frame_name << frame_index;
    const fs::path frame_dir = trajectory_dir / frame_name.str();
    ASSERT_TRUE(fs::is_directory(frame_dir));

    const auto& local_trajectory =
        tp.Result<nmr::LocalBackboneGeometryTrajectoryResult>();
    CheckLocalBackboneAgainstTrajectoryH5(
        tp, local_trajectory, frame_dir, root / "local_backbone.h5");
    CheckHBondTypedRows(tp.ProteinRef(), frame_dir);
    CheckCoordinateTorsionOracles(tp.ProteinRef(), frame_dir);

    // NumPy is the final reader: every new static and per-frame array is
    // opened with allow_pickle=False, then direct coordinate checks run.
    const std::string command =
        ShellQuote(NMR_TEST_PYTHON_EXECUTABLE) + " " +
        ShellQuote(NMR_LOCAL_GEOMETRY_ACCEPTANCE_SCRIPT) + " " +
        ShellQuote(static_dir) + " " + ShellQuote(transformed_dir) + " " +
        ShellQuote(frame_dir);
    EXPECT_EQ(std::system(command.c_str()), 0) << command;

    nmr::FrameNpyEmitter::Reset();
    RemoveTempTree(root);
}


TEST(SolventPbcInvariance, AllAffectedArraysAndProteinPath) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    const fs::path root = FreshRoot();
    const fs::path fixture_root = NMR_LOCAL_GEOMETRY_TRAJECTORY_FIXTURE;
    const fs::path production_dir = fixture_root / "production";
    const fs::path tpr_path = production_dir / "production.tpr";
    const fs::path source_trr_path = production_dir / "production.trr";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(production_dir.string())) << tp.Error();
    const auto& topology = tp.SysReader().Topology();
    ASSERT_GT(topology.protein_count, 0u);
    ASSERT_GT(topology.water_count, 0u);
    ASSERT_GT(topology.ion_count, 0u);

    int64_t step = 0;
    real time = 0;
    real lambda = 0;
    rvec box[3];
    int natoms = static_cast<int>(topology.total_atoms);
    std::vector<float> baseline_xyz(topology.total_atoms * 3);
    gmx_trr_read_single_frame(
        source_trr_path, &step, &time, &lambda, box, &natoms,
        reinterpret_cast<rvec*>(baseline_xyz.data()), nullptr, nullptr);
    ASSERT_EQ(natoms, static_cast<int>(topology.total_atoms));
    ASSERT_NE(box[2][0], 0.0f)
        << "fixture must exercise a triclinic, not diagonal-only, cell";

    // Pin the separately retained full topology to the validated wholer:
    // one water H is moved across a cell face and must return beside its O.
    std::vector<float> wholer_probe = baseline_xyz;
    const std::size_t probe_oxygen = topology.water_O_start;
    ShiftTriclinicImage(
        wholer_probe, probe_oxygen + 1, box, 1, 0, 0);
    float wholer_box[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) wholer_box[i][j] = box[i][j];
    }
    ASSERT_TRUE(tp.SysReader().MakeSystemWhole(wholer_probe, wholer_box));
    const nmr::Vec3 probe_oh(
        static_cast<double>(wholer_probe[(probe_oxygen + 1) * 3]) -
            static_cast<double>(wholer_probe[probe_oxygen * 3]),
        static_cast<double>(wholer_probe[(probe_oxygen + 1) * 3 + 1]) -
            static_cast<double>(wholer_probe[probe_oxygen * 3 + 1]),
        static_cast<double>(wholer_probe[(probe_oxygen + 1) * 3 + 2]) -
            static_cast<double>(wholer_probe[probe_oxygen * 3 + 2]));
    EXPECT_LT(probe_oh.norm() * nmr::ANGSTROM_PER_NANOMETRE, 2.0);
    std::cout << "PBC_INVARIANCE full_topology_wholer repaired_probe_water\n";

    std::vector<float> shifted_xyz = baseline_xyz;
    const std::array<std::array<int, 3>, 6> molecule_shifts{{
        {{ 1,  0, -1}}, {{-1,  1,  0}}, {{ 0, -1,  1}},
        {{ 1, -1,  1}}, {{-1,  0, -1}}, {{ 0,  1, -1}},
    }};
    for (std::size_t wi = 0; wi < topology.water_count; ++wi) {
        const std::size_t oxygen = topology.water_O_start + wi * 3;
        const auto& shift = molecule_shifts[wi % molecule_shifts.size()];
        for (std::size_t site = 0; site < 3; ++site) {
            ShiftTriclinicImage(shifted_xyz, oxygen + site, box,
                                shift[0], shift[1], shift[2]);
        }
        // The whole water was moved to an equivalent image; now put each H
        // in another equivalent image too. This is the torn-water case that
        // corrupts WaterMolecule::Dipole when solvent bypasses do_pbc_mtop.
        ShiftTriclinicImage(shifted_xyz, oxygen + 1, box, 1, -1, 0);
        ShiftTriclinicImage(shifted_xyz, oxygen + 2, box, -1, 0, 1);
    }
    for (std::size_t ii = 0; ii < topology.ion_count; ++ii) {
        const auto& shift = molecule_shifts[ii % molecule_shifts.size()];
        ShiftTriclinicImage(shifted_xyz, topology.ion_start + ii, box,
                            shift[0], shift[1], shift[2]);
    }

    const fs::path baseline_trr = root / "baseline.trr";
    const fs::path shifted_trr = root / "shifted.trr";
    gmx_trr_write_single_frame(
        baseline_trr, step, time, lambda, box, natoms,
        reinterpret_cast<const rvec*>(baseline_xyz.data()), nullptr, nullptr);
    gmx_trr_write_single_frame(
        shifted_trr, step, time, lambda, box, natoms,
        reinterpret_cast<const rvec*>(shifted_xyz.data()), nullptr, nullptr);

    nmr::GromacsFrameHandler baseline_handler(tp);
    nmr::GromacsFrameHandler shifted_handler(tp);
    ASSERT_TRUE(baseline_handler.Open(baseline_trr.string(), tpr_path.string()))
        << baseline_handler.error();
    ASSERT_TRUE(shifted_handler.Open(shifted_trr.string(), tpr_path.string()))
        << shifted_handler.error();
    ASSERT_TRUE(baseline_handler.ReadNextFrame()) << baseline_handler.error();
    ASSERT_TRUE(shifted_handler.ReadNextFrame()) << shifted_handler.error();

    const auto& baseline_solvent = baseline_handler.Solvent();
    const auto& shifted_solvent = shifted_handler.Solvent();
    ASSERT_EQ(shifted_solvent.waters.size(), baseline_solvent.waters.size());
    double max_dipole_delta = 0.0;
    for (std::size_t wi = 0; wi < baseline_solvent.waters.size(); ++wi) {
        const auto& baseline_water = baseline_solvent.waters[wi];
        const auto& shifted_water = shifted_solvent.waters[wi];
        EXPECT_LT((baseline_water.H1_pos - baseline_water.O_pos).norm(), 2.0);
        EXPECT_LT((baseline_water.H2_pos - baseline_water.O_pos).norm(), 2.0);
        EXPECT_LT((shifted_water.H1_pos - shifted_water.O_pos).norm(), 2.0);
        EXPECT_LT((shifted_water.H2_pos - shifted_water.O_pos).norm(), 2.0);
        max_dipole_delta = std::max(
            max_dipole_delta,
            (baseline_water.Dipole() - shifted_water.Dipole()).norm());
    }
    EXPECT_LE(max_dipole_delta, 5.0e-5);
    std::cout << "PBC_INVARIANCE stored_water_dipoles max_abs="
              << max_dipole_delta << '\n';

    const auto& baseline_protein = baseline_handler.ProteinPositions();
    const auto& shifted_protein = shifted_handler.ProteinPositions();
    ASSERT_EQ(shifted_protein.size(), baseline_protein.size());
    for (std::size_t ai = 0; ai < baseline_protein.size(); ++ai) {
        EXPECT_EQ(shifted_protein[ai].x(), baseline_protein[ai].x()) << ai;
        EXPECT_EQ(shifted_protein[ai].y(), baseline_protein[ai].y()) << ai;
        EXPECT_EQ(shifted_protein[ai].z(), baseline_protein[ai].z()) << ai;
    }

    // The new full-system pass must not become an alternative protein path.
    // Re-run the previously validated protein-only wholer directly and require
    // the handler's emitted protein coordinates to be exactly its output.
    std::vector<float> validated_protein(
        baseline_xyz.begin() + topology.protein_start * 3,
        baseline_xyz.begin() +
            (topology.protein_start + topology.protein_count) * 3);
    float protein_box[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) protein_box[i][j] = box[i][j];
    }
    ASSERT_TRUE(tp.SysReader().MakeProteinWhole(
        validated_protein, protein_box));
    for (std::size_t ai = 0; ai < baseline_protein.size(); ++ai) {
        for (std::size_t d = 0; d < 3; ++d) {
            EXPECT_EQ(baseline_protein[ai][static_cast<Eigen::Index>(d)],
                      static_cast<double>(validated_protein[ai * 3 + d]) *
                          nmr::ANGSTROM_PER_NANOMETRE)
                << ai << "," << d;
        }
    }
    std::cout << "PBC_INVARIANCE protein_positions exact_validated_path\n";

    tp.Seed(baseline_protein, baseline_handler.Time());
    auto baseline_conf = tp.TickConformation(baseline_protein);
    auto shifted_conf = tp.TickConformation(shifted_protein);
    AttachExplicitSolventOwners(*baseline_conf, baseline_handler.Solvent());
    AttachExplicitSolventOwners(*shifted_conf, shifted_handler.Solvent());

    ASSERT_TRUE(baseline_handler.Solvent().HasPeriodicCell());
    ASSERT_TRUE(shifted_handler.Solvent().HasPeriodicCell());
    EXPECT_TRUE((baseline_handler.Solvent().BoxMatrix().array() ==
                 baseline_handler.BoxMatrix().array()).all());
    EXPECT_TRUE((shifted_handler.Solvent().BoxMatrix().array() ==
                 shifted_handler.BoxMatrix().array()).all());

    const fs::path baseline_dir = root / "baseline_arrays";
    const fs::path shifted_dir = root / "shifted_arrays";
    WriteExplicitSolventOwners(*baseline_conf, baseline_dir);
    WriteExplicitSolventOwners(*shifted_conf, shifted_dir);
    auto baseline_frame = nmr::GromacsFramePullResult::Compute(
        *baseline_conf, nullptr, &baseline_handler.BoxMatrix());
    auto shifted_frame = nmr::GromacsFramePullResult::Compute(
        *shifted_conf, nullptr, &shifted_handler.BoxMatrix());
    ASSERT_NE(baseline_frame, nullptr);
    ASSERT_NE(shifted_frame, nullptr);
    ASSERT_EQ(baseline_frame->WriteFeatures(
                  *baseline_conf, baseline_dir.string()), 1);
    ASSERT_EQ(shifted_frame->WriteFeatures(
                  *shifted_conf, shifted_dir.string()), 1);

    const NpyArray emitted_box = ReadNpy(baseline_dir / "gromacs_box.npy");
    ASSERT_EQ(emitted_box.shape, (std::vector<std::size_t>{1, 9}));
    const auto emitted_box_values = Values<double>(emitted_box);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            const std::size_t column = static_cast<std::size_t>(i * 3 + j);
            EXPECT_EQ(emitted_box_values[column],
                      baseline_handler.BoxMatrix()(i, j));
            EXPECT_EQ(emitted_box_values[column],
                      static_cast<double>(box[j][i]) *
                          nmr::ANGSTROM_PER_NANOMETRE);
        }
    }
    std::cout << "PBC_INVARIANCE gromacs_box serialized_layout_units_exact\n";

    const NpyArray baseline_shell = ReadNpy(
        baseline_dir / "water_shell_counts.npy");
    const auto shell_values = Values<double>(baseline_shell);
    ASSERT_TRUE(std::any_of(shell_values.begin(), shell_values.end(),
                            [](double value) { return value > 0.0; }));
    const NpyArray baseline_candidates = ReadNpy(
        baseline_dir / "water_hbond_candidates.npy");
    ASSERT_EQ(baseline_candidates.shape.size(), 2u);
    ASSERT_GT(baseline_candidates.shape[0], 0u);

    constexpr double kAbs = 2.0e-4;
    constexpr double kRel = 2.0e-6;
    for (const char* filename : {
             "water_efield.npy",
             "water_efg.npy",
             "water_efg_first.npy",
             "water_efield_first.npy",
             "water_shell_counts.npy",
             "water_efield_clamp_mask.npy",
             "water_efield_clamp_scale.npy",
             "water_efield_first_clamp_mask.npy",
             "water_efield_first_clamp_scale.npy",
             "hydration_shell.npy",
             "water_polarization.npy",
             "water_hbond_counts.npy",
             "gromacs_box.npy",
         }) {
        ExpectPeriodicNpyEqual(
            baseline_dir, shifted_dir, filename, kAbs, kRel);
    }
    // Candidate identities, modes, roles, pass decisions, and summary indices
    // are exact. Continuous H-bond geometry includes angles in degrees; its
    // 0.002 tolerance reflects arithmetic on the single-precision TRR inputs.
    ExpectExactFloatColumns(
        baseline_dir, shifted_dir, "water_hbond_candidates.npy",
        {0, 1, 2, 3, 4, 15});
    ExpectPeriodicNpyEqual(
        baseline_dir, shifted_dir, "water_hbond_candidates.npy", 2.0e-3, kRel);
    ExpectExactFloatColumns(
        baseline_dir, shifted_dir, "water_hbond_nearest.npy",
        {3, 4, 5, 6, 7});
    ExpectPeriodicNpyEqual(
        baseline_dir, shifted_dir, "water_hbond_nearest.npy", 2.0e-3, kRel);

    RemoveTempTree(root);
}

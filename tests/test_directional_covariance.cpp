#include "DirectionalTestHelpers.h"
#include "TestEnvironment.h"

#include <gtest/gtest.h>

#include "BiotSavartResult.h"
#include "ChargeAssignmentResult.h"
#include "CoulombResult.h"
#include "DsspResult.h"
#include "EnrichmentResult.h"
#include "EeqCoulombResult.h"
#include "EeqResult.h"
#include "GeometryResult.h"
#include "HBondResult.h"
#include "HaighMallionResult.h"
#include "HydrationGeometryResult.h"
#include "LocalBackboneGeometryResult.h"
#include "McConnellResult.h"
#include "PdbFileReader.h"
#include "PiQuadrupoleLocalTensorResult.h"
#include "PiQuadrupoleResult.h"
#include "PlanarGeometryResult.h"
#include "Protein.h"
#include "DispersionResult.h"
#include "RingSusceptibilityResult.h"
#include "SasaResult.h"
#include "SidechainCarbonylAnisotropyResult.h"
#include "SpatialIndexResult.h"
#include "WaterHBondGeometryResult.h"
#include "WaterFieldResult.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <unistd.h>

namespace {

using nmr::Mat3;
using nmr::Protein;
using nmr::ProteinConformation;
using nmr::SphericalTensor;
using nmr::Vec3;
using nmr::test::directional::Axial;
using nmr::test::directional::CircularDifference;
using nmr::test::directional::EvenRank2;
using nmr::test::directional::Near;
using nmr::test::directional::NearMatrix;
using nmr::test::directional::NearVector;
using nmr::test::directional::OrthogonalTransform;
using nmr::test::directional::Polar;
using nmr::test::directional::Position;
using nmr::test::directional::RotateNativeT2;
using nmr::test::directional::SeededTransform;
using nmr::test::directional::T1Vector;
using nmr::test::directional::T2Matrix;

// Recorded tolerances.  Analytic distance/vector/tensor calculators are
// expected to differ only by floating-point reassociation after a rigid O(3)
// transform.  HM additionally uses adaptive surface quadrature, hence its
// slightly wider (still sub-ppb relative) bound.  SASA has a separate,
// explicitly discretization-derived bound at its test site below.
constexpr double kGeometryAbs = 2.0e-11;
constexpr double kGeometryRel = 2.0e-12;
constexpr double kKernelAbs = 2.0e-12;
constexpr double kKernelRel = 2.0e-10;
constexpr double kHmAbs = 2.0e-11;
constexpr double kHmRel = 8.0e-9;
constexpr double kElectrostaticAbs = 2.0e-9;
constexpr double kElectrostaticRel = 2.0e-10;
constexpr double kStructuralZeroAbs = 2.0e-9;

struct SerializedNpy {
    std::string descr;
    std::vector<std::size_t> shape;
    std::vector<char> payload;
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

SerializedNpy ReadSerializedNpy(const std::filesystem::path& path) {
    SerializedNpy out;
    std::ifstream input(path, std::ios::binary);
    EXPECT_TRUE(input.is_open()) << path;
    if (!input.is_open()) return out;

    char magic[6] = {};
    input.read(magic, sizeof(magic));
    EXPECT_EQ(std::string(magic, sizeof(magic)),
              std::string("\x93NUMPY", sizeof(magic)));
    unsigned char version[2] = {};
    input.read(reinterpret_cast<char*>(version), sizeof(version));
    EXPECT_TRUE((version[0] == 1 || version[0] == 2) && version[1] == 0);

    std::uint32_t header_length = 0;
    if (version[0] == 1) {
        unsigned char bytes[2] = {};
        input.read(reinterpret_cast<char*>(bytes), sizeof(bytes));
        header_length = static_cast<std::uint32_t>(bytes[0]) |
                        (static_cast<std::uint32_t>(bytes[1]) << 8);
    } else {
        unsigned char bytes[4] = {};
        input.read(reinterpret_cast<char*>(bytes), sizeof(bytes));
        header_length = static_cast<std::uint32_t>(bytes[0]) |
                        (static_cast<std::uint32_t>(bytes[1]) << 8) |
                        (static_cast<std::uint32_t>(bytes[2]) << 16) |
                        (static_cast<std::uint32_t>(bytes[3]) << 24);
    }

    std::string header(header_length, '\0');
    input.read(header.data(), static_cast<std::streamsize>(header.size()));
    const std::string descr_key = "'descr': '";
    const std::size_t descr_begin = header.find(descr_key);
    EXPECT_NE(descr_begin, std::string::npos) << header;
    if (descr_begin != std::string::npos) {
        const std::size_t value_begin = descr_begin + descr_key.size();
        const std::size_t value_end = header.find('\'', value_begin);
        EXPECT_NE(value_end, std::string::npos) << header;
        if (value_end != std::string::npos)
            out.descr = header.substr(value_begin, value_end - value_begin);
    }

    const std::size_t shape_begin = header.find('(');
    const std::size_t shape_end = header.find(')', shape_begin);
    EXPECT_NE(shape_begin, std::string::npos) << header;
    EXPECT_NE(shape_end, std::string::npos) << header;
    if (shape_begin != std::string::npos && shape_end != std::string::npos) {
        std::stringstream shape_stream(
            header.substr(shape_begin + 1, shape_end - shape_begin - 1));
        std::string token;
        while (std::getline(shape_stream, token, ',')) {
            token = Trim(token);
            if (!token.empty())
                out.shape.push_back(static_cast<std::size_t>(
                    std::stoull(token)));
        }
    }

    out.payload.assign(std::istreambuf_iterator<char>(input),
                       std::istreambuf_iterator<char>());
    return out;
}

template <typename T>
std::vector<T> NpyValues(const SerializedNpy& array,
                         const char* expected_descr) {
    EXPECT_EQ(array.descr, expected_descr);
    EXPECT_EQ(array.payload.size() % sizeof(T), 0u);
    std::vector<T> values(array.payload.size() / sizeof(T));
    if (!values.empty())
        std::memcpy(values.data(), array.payload.data(), array.payload.size());
    return values;
}

SphericalTensor UnpackFull9(const double* packed) {
    SphericalTensor out;
    out.T0 = packed[0];
    for (int component = 0; component < 3; ++component)
        out.T1[component] = packed[1 + component];
    for (int component = 0; component < 5; ++component)
        out.T2[component] = packed[4 + component];
    return out;
}

void ExpectVector(const Vec3& actual, const Vec3& expected,
                  double abs_tolerance, double rel_tolerance,
                  const char* quantity);
void ExpectEvenSphericalTensor(const SphericalTensor& transformed,
                               const SphericalTensor& original,
                               const OrthogonalTransform& x,
                               double abs_tolerance,
                               double rel_tolerance,
                               const char* quantity);

void CheckSerializedPolarNpy(const std::filesystem::path& original_dir,
                             const std::filesystem::path& moved_dir,
                             const char* filename,
                             std::size_t row_count,
                             const OrthogonalTransform& x,
                             double abs_tolerance,
                             double rel_tolerance) {
    const SerializedNpy source_npy = ReadSerializedNpy(original_dir / filename);
    const SerializedNpy moved_npy = ReadSerializedNpy(moved_dir / filename);
    ASSERT_EQ(source_npy.shape,
              (std::vector<std::size_t>{row_count, 3u})) << filename;
    ASSERT_EQ(moved_npy.shape, source_npy.shape) << filename;
    const auto source = NpyValues<double>(source_npy, "<f8");
    const auto moved = NpyValues<double>(moved_npy, "<f8");
    for (std::size_t row = 0; row < row_count; ++row) {
        const Vec3 a(source[row * 3], source[row * 3 + 1],
                     source[row * 3 + 2]);
        const Vec3 b(moved[row * 3], moved[row * 3 + 1],
                     moved[row * 3 + 2]);
        ExpectVector(b, Polar(x, a), abs_tolerance, rel_tolerance, filename);
    }
}

void CheckSerializedFull9Npy(const std::filesystem::path& original_dir,
                             const std::filesystem::path& moved_dir,
                             const char* filename,
                             std::size_t row_count,
                             const OrthogonalTransform& x,
                             double abs_tolerance,
                             double rel_tolerance,
                             bool structural_t0_t1_zero) {
    const SerializedNpy source_npy = ReadSerializedNpy(original_dir / filename);
    const SerializedNpy moved_npy = ReadSerializedNpy(moved_dir / filename);
    ASSERT_EQ(source_npy.shape,
              (std::vector<std::size_t>{row_count, 9u})) << filename;
    ASSERT_EQ(moved_npy.shape, source_npy.shape) << filename;
    const auto source = NpyValues<double>(source_npy, "<f8");
    const auto moved = NpyValues<double>(moved_npy, "<f8");
    for (std::size_t row = 0; row < row_count; ++row) {
        const SphericalTensor a = UnpackFull9(source.data() + row * 9);
        const SphericalTensor b = UnpackFull9(moved.data() + row * 9);
        ExpectEvenSphericalTensor(b, a, x, abs_tolerance, rel_tolerance,
                                  filename);
        if (structural_t0_t1_zero) {
            EXPECT_LE(std::abs(b.T0), kStructuralZeroAbs) << filename;
            EXPECT_LE(T1Vector(b).norm(), kStructuralZeroAbs) << filename;
        }
    }
}

void CheckSerializedNativeT2Npy(const std::filesystem::path& original_dir,
                                const std::filesystem::path& moved_dir,
                                const char* filename,
                                std::size_t row_count,
                                const OrthogonalTransform& x,
                                double abs_tolerance,
                                double rel_tolerance) {
    const SerializedNpy source_npy = ReadSerializedNpy(original_dir / filename);
    const SerializedNpy moved_npy = ReadSerializedNpy(moved_dir / filename);
    ASSERT_EQ(source_npy.shape,
              (std::vector<std::size_t>{row_count, 5u})) << filename;
    ASSERT_EQ(moved_npy.shape, source_npy.shape) << filename;
    const auto source = NpyValues<double>(source_npy, "<f8");
    const auto moved = NpyValues<double>(moved_npy, "<f8");
    for (std::size_t row = 0; row < row_count; ++row) {
        SphericalTensor a;
        SphericalTensor b;
        for (std::size_t component = 0; component < 5; ++component) {
            a.T2[component] = source[row * 5 + component];
            b.T2[component] = moved[row * 5 + component];
        }
        const SphericalTensor expected = RotateNativeT2(x, a);
        for (std::size_t component = 0; component < 5; ++component) {
            EXPECT_TRUE(Near(b.T2[component], expected.T2[component],
                             abs_tolerance, rel_tolerance))
                << filename << " row=" << row
                << " component=" << component;
        }
    }
}

std::filesystem::path FreshOutputDirectory(const std::string& stem) {
    const std::filesystem::path path =
        nmr::test::TestEnvironment::TempPath(stem);
    std::error_code error;
    std::filesystem::create_directories(path, error);
    EXPECT_FALSE(error) << path << ": " << error.message();
    return path;
}

void RemoveOutputDirectory(const std::filesystem::path& path) {
    // Do not use std::filesystem::remove_all in this torch-linked binary:
    // libtorch's bundled filesystem symbols have caused teardown crashes in
    // this test target.  Every producer here writes only regular files into
    // its dedicated flat directory, so remove those files explicitly.
    std::error_code error;
    for (const auto& entry : std::filesystem::directory_iterator(path, error)) {
        EXPECT_TRUE(entry.is_regular_file()) << entry.path();
        EXPECT_EQ(std::remove(entry.path().string().c_str()), 0)
            << entry.path();
    }
    EXPECT_FALSE(error) << path << ": " << error.message();
    EXPECT_EQ(::rmdir(path.string().c_str()), 0) << path;
}

void ExpectScalarOrSameNaN(double actual, double expected,
                           double abs_tolerance, double rel_tolerance,
                           const char* quantity) {
    if (std::isnan(expected)) {
        EXPECT_TRUE(std::isnan(actual)) << quantity;
        return;
    }
    ASSERT_TRUE(std::isfinite(expected)) << quantity;
    EXPECT_TRUE(std::isfinite(actual)) << quantity;
    EXPECT_TRUE(Near(actual, expected, abs_tolerance, rel_tolerance))
        << quantity << " actual=" << actual << " expected=" << expected;
}

void ExpectVector(const Vec3& actual, const Vec3& expected,
                  double abs_tolerance, double rel_tolerance,
                  const char* quantity) {
    EXPECT_TRUE(NearVector(actual, expected, abs_tolerance, rel_tolerance))
        << quantity << " actual=" << actual.transpose()
        << " expected=" << expected.transpose()
        << " error=" << (actual - expected).norm()
        << " tolerance="
        << abs_tolerance + rel_tolerance * expected.norm();
}

void ExpectMatrix(const Mat3& actual, const Mat3& expected,
                  double abs_tolerance, double rel_tolerance,
                  const char* quantity) {
    EXPECT_TRUE(NearMatrix(actual, expected, abs_tolerance, rel_tolerance))
        << quantity << " error=" << (actual - expected).norm()
        << " tolerance="
        << abs_tolerance + rel_tolerance * expected.norm();
}

void ExpectEvenSphericalTensor(const SphericalTensor& transformed,
                               const SphericalTensor& original,
                               const OrthogonalTransform& x,
                               double abs_tolerance,
                               double rel_tolerance,
                               const char* quantity) {
    const Mat3 expected = EvenRank2(x, original.Reconstruct());
    ExpectMatrix(transformed.Reconstruct(), expected,
                 abs_tolerance, rel_tolerance, quantity);
    EXPECT_TRUE(Near(transformed.T0, original.T0,
                     abs_tolerance, rel_tolerance))
        << quantity << " T0 must be invariant";
    ExpectVector(T1Vector(transformed), Axial(x, T1Vector(original)),
                 abs_tolerance, rel_tolerance, quantity);

    // Independent project-native five-component path: reconstruct only T2,
    // rotate Cartesian Q*T*Q^T, decompose, and compare with the T2 produced by
    // rerunning the calculator on transformed coordinates.
    const SphericalTensor expected_t2 = RotateNativeT2(x, original);
    for (int component = 0; component < 5; ++component) {
        EXPECT_TRUE(Near(transformed.T2[component],
                         expected_t2.T2[component],
                         abs_tolerance, rel_tolerance))
            << quantity << " native T2 component " << component;
    }
}

void ExpectEfg(const SphericalTensor& transformed,
               const SphericalTensor& original,
               const OrthogonalTransform& x,
               const char* quantity) {
    ExpectEvenSphericalTensor(transformed, original, x,
                              kElectrostaticAbs, kElectrostaticRel,
                              quantity);
    EXPECT_LE(std::abs(original.T0), kStructuralZeroAbs) << quantity;
    EXPECT_LE(T1Vector(original).norm(), kStructuralZeroAbs) << quantity;
    EXPECT_LE(std::abs(transformed.T0), kStructuralZeroAbs) << quantity;
    EXPECT_LE(T1Vector(transformed).norm(), kStructuralZeroAbs) << quantity;
}

bool AttachGeometryAndSpatial(ProteinConformation& conf) {
    auto geometry = nmr::GeometryResult::Compute(conf);
    if (!geometry || !conf.AttachResult(std::move(geometry))) return false;
    auto spatial = nmr::SpatialIndexResult::Compute(conf);
    return spatial && conf.AttachResult(std::move(spatial));
}

bool AttachRingStack(ProteinConformation& conf) {
    if (!AttachGeometryAndSpatial(conf)) return false;
    auto bs = nmr::BiotSavartResult::Compute(conf);
    if (!bs || !conf.AttachResult(std::move(bs))) return false;
    auto hm = nmr::HaighMallionResult::Compute(conf);
    if (!hm || !conf.AttachResult(std::move(hm))) return false;
    auto mc = nmr::McConnellResult::Compute(conf);
    if (!mc || !conf.AttachResult(std::move(mc))) return false;
    auto pq = nmr::PiQuadrupoleResult::Compute(conf);
    return pq && conf.AttachResult(std::move(pq));
}

std::unique_ptr<Protein> LoadUbqOrSkip() {
    const std::string path = nmr::test::TestEnvironment::UbqProtonated();
    if (path.empty() || !std::filesystem::exists(path)) return nullptr;
    auto loaded = nmr::BuildFromProtonatedPdb(path);
    if (!loaded.Ok()) return nullptr;
    return std::move(loaded.protein);
}

const nmr::AtomNeighbour* FindSpatialNeighbour(
        const ProteinConformation& conf, std::size_t atom,
        std::size_t neighbour) {
    for (const nmr::AtomNeighbour& row :
         conf.AtomAt(atom).spatial_neighbours) {
        if (row.atom_index == neighbour) return &row;
    }
    return nullptr;
}

void CheckGeometrySpatialAndRingStack(const ProteinConformation& original,
                                      const ProteinConformation& transformed,
                                      const OrthogonalTransform& x) {
    ASSERT_EQ(original.AtomCount(), transformed.AtomCount());
    const Protein& protein = original.ProteinRef();

    for (std::size_t atom = 0; atom < original.AtomCount(); ++atom) {
        SCOPED_TRACE("atom=" + std::to_string(atom));
        ExpectVector(transformed.PositionAt(atom),
                     Position(x, original.PositionAt(atom)),
                     kGeometryAbs, kGeometryRel, "pos.npy");

        const auto& original_rows = original.AtomAt(atom).spatial_neighbours;
        const auto& transformed_rows =
            transformed.AtomAt(atom).spatial_neighbours;
        ASSERT_EQ(transformed_rows.size(), original_rows.size())
            << "spatial_neighbors.npy sparse identity changed";
        for (const nmr::AtomNeighbour& row : original_rows) {
            const nmr::AtomNeighbour* moved =
                FindSpatialNeighbour(transformed, atom, row.atom_index);
            ASSERT_NE(moved, nullptr)
                << "spatial_neighbors.npy lost atom pair";
            EXPECT_TRUE(Near(moved->distance, row.distance,
                             kGeometryAbs, kGeometryRel));
            ExpectVector(moved->direction, Polar(x, row.direction),
                         kGeometryAbs, kGeometryRel,
                         "spatial_neighbors[:,2:5]");
        }

        const auto& original_atom = original.AtomAt(atom);
        const auto& transformed_atom = transformed.AtomAt(atom);
        ExpectVector(transformed_atom.total_B_field,
                     Axial(x, original_atom.total_B_field),
                     kKernelAbs, kKernelRel, "bs_total_B.npy");
        ExpectEvenSphericalTensor(transformed_atom.bs_shielding_contribution,
                                  original_atom.bs_shielding_contribution,
                                  x, kKernelAbs, kKernelRel,
                                  "bs_shielding.npy");
        ExpectEvenSphericalTensor(transformed_atom.hm_shielding_contribution,
                                  original_atom.hm_shielding_contribution,
                                  x, kHmAbs, kHmRel,
                                  "hm_shielding.npy");

        ASSERT_EQ(transformed_atom.ring_neighbours.size(),
                  original_atom.ring_neighbours.size());
        for (std::size_t row_index = 0;
             row_index < original_atom.ring_neighbours.size(); ++row_index) {
            const auto& row = original_atom.ring_neighbours[row_index];
            const auto& moved =
                transformed_atom.ring_neighbours[row_index];
            ASSERT_EQ(moved.ring_index, row.ring_index);
            SCOPED_TRACE("ring=" + std::to_string(row.ring_index));
            EXPECT_TRUE(Near(moved.distance_to_center,
                             row.distance_to_center,
                             kGeometryAbs, kGeometryRel));
            ExpectVector(moved.direction_to_center,
                         Polar(x, row.direction_to_center),
                         kGeometryAbs, kGeometryRel,
                         "ring_direction_to_center.npy");
            EXPECT_TRUE(Near(moved.z,
                             x.Determinant() * row.z,
                             kGeometryAbs, kGeometryRel))
                << "ring_contributions z is a pseudoscalar";
            EXPECT_TRUE(Near(moved.theta, row.theta,
                             kGeometryAbs, kGeometryRel))
                << "ring_contributions theta=atan2(rho,abs(z)) is invariant";
            ExpectVector(moved.B_field, Axial(x, row.B_field),
                         kKernelAbs, kKernelRel,
                         "bs_ring_B_field.npy");
            EXPECT_TRUE(Near(moved.B_cylindrical.x(),
                             x.Determinant() * row.B_cylindrical.x(),
                             kKernelAbs, kKernelRel))
                << "bs_ring_B_cylindrical.npy col0 is a pseudoscalar";
            EXPECT_DOUBLE_EQ(row.B_cylindrical.y(), 0.0)
                << "bs_ring_B_cylindrical.npy col1 structural zero";
            EXPECT_DOUBLE_EQ(moved.B_cylindrical.y(), 0.0)
                << "bs_ring_B_cylindrical.npy col1 structural zero";
            EXPECT_TRUE(Near(moved.B_cylindrical.z(),
                             row.B_cylindrical.z(),
                             kKernelAbs, kKernelRel))
                << "bs_ring_B_cylindrical.npy col2 is invariant";
            ExpectMatrix(moved.G_tensor, EvenRank2(x, row.G_tensor),
                         kKernelAbs, kKernelRel,
                         "ring_contributions BS G");
            ExpectMatrix(moved.hm_H_tensor,
                         EvenRank2(x, row.hm_H_tensor),
                         kHmAbs, kHmRel,
                         "ring_contributions HM H");
            ExpectVector(moved.hm_B_field, Axial(x, row.hm_B_field),
                         kHmAbs, kHmRel,
                         "hm_ring_B_field.npy");
            ExpectMatrix(moved.hm_G_tensor,
                         EvenRank2(x, row.hm_G_tensor),
                         kHmAbs, kHmRel,
                         "ring_contributions HM G");
            EXPECT_LE(std::abs(moved.hm_H_spherical.T0),
                      kStructuralZeroAbs);
            EXPECT_LE(T1Vector(moved.hm_H_spherical).norm(),
                      kStructuralZeroAbs);
            EXPECT_TRUE(Near(moved.quad_scalar, row.quad_scalar,
                             kKernelAbs, kKernelRel))
                << "piquad scalar must be O(3)-invariant";
        }

        for (std::size_t category = 0;
             category < nmr::kMcConnellSourceCategoryCount; ++category) {
            for (std::size_t channel = 0;
                 channel < nmr::kMcConnellChannelCount; ++channel) {
                ExpectEvenSphericalTensor(
                    transformed_atom.mcconnell_source_tensors[category][channel],
                    original_atom.mcconnell_source_tensors[category][channel],
                    x, kKernelAbs, kKernelRel,
                    "mc_<category>_<channel>.npy");
            }
        }
        ExpectEvenSphericalTensor(
            transformed_atom.mcconnell_peptide_co_rhombic,
            original_atom.mcconnell_peptide_co_rhombic,
            x, kKernelAbs, kKernelRel,
            "mc_peptide_co_rhombic.npy");
        if (original_atom.nearest_CO_dist > 0.0) {
            ExpectVector(transformed_atom.nearest_CO_midpoint,
                         Position(x, original_atom.nearest_CO_midpoint),
                         kGeometryAbs, kGeometryRel,
                         "mc_nearest_co_midpoint.npy");
            ExpectVector(transformed_atom.dir_nearest_CO,
                         Polar(x, original_atom.dir_nearest_CO),
                         kGeometryAbs, kGeometryRel,
                         "mc_nearest_co_dir.npy");
            ExpectEvenSphericalTensor(transformed_atom.T2_CO_nearest,
                                      original_atom.T2_CO_nearest,
                                      x, kKernelAbs, kKernelRel,
                                      "mc_nearest_co_T2.npy");
        }
        if (original_atom.nearest_CN_dist > 0.0) {
            ExpectEvenSphericalTensor(transformed_atom.T2_CN_nearest,
                                      original_atom.T2_CN_nearest,
                                      x, kKernelAbs, kKernelRel,
                                      "mc_nearest_cn_T2.npy");
        }
    }

    ASSERT_EQ(original.bond_directions.size(),
              transformed.bond_directions.size());
    for (std::size_t bond = 0; bond < original.bond_directions.size(); ++bond) {
        SCOPED_TRACE("bond=" + std::to_string(bond));
        EXPECT_EQ(transformed.bond_geometry_valid[bond],
                  original.bond_geometry_valid[bond]);
        EXPECT_TRUE(Near(transformed.bond_lengths[bond],
                         original.bond_lengths[bond],
                         kGeometryAbs, kGeometryRel));
        ExpectVector(transformed.bond_directions[bond],
                     Polar(x, original.bond_directions[bond]),
                     kGeometryAbs, kGeometryRel, "bond_direction.npy");
        ExpectVector(transformed.bond_midpoints[bond],
                     Position(x, original.bond_midpoints[bond]),
                     kGeometryAbs, kGeometryRel, "bond midpoint");
    }

    ASSERT_EQ(original.ring_geometries.size(),
              transformed.ring_geometries.size());
    for (std::size_t ring = 0; ring < protein.RingCount(); ++ring) {
        const auto& geometry = original.ring_geometries[ring];
        const auto& moved = transformed.ring_geometries[ring];
        SCOPED_TRACE("ring_geometry=" + std::to_string(ring));
        ExpectVector(moved.center, Position(x, geometry.center),
                     kGeometryAbs, kGeometryRel,
                     "ring_geometry[:,3:6]");
        ExpectVector(moved.normal, Axial(x, geometry.normal),
                     kGeometryAbs, kGeometryRel,
                     "ring_geometry[:,6:9]");
        EXPECT_TRUE(Near(moved.radius, geometry.radius,
                         kGeometryAbs, kGeometryRel));
    }

    // ring_pair_geometry signed-offset columns are produced directly from
    // these centers/normals by WriteAllFeatures.  Pin their odd parity here.
    for (std::size_t a = 0; a < protein.RingCount(); ++a) {
        for (std::size_t b = a + 1; b < protein.RingCount(); ++b) {
            const Vec3 delta = original.ring_geometries[b].center -
                               original.ring_geometries[a].center;
            const Vec3 moved_delta = transformed.ring_geometries[b].center -
                                     transformed.ring_geometries[a].center;
            const double za = delta.dot(original.ring_geometries[a].normal);
            const double zb = delta.dot(original.ring_geometries[b].normal);
            EXPECT_TRUE(Near(moved_delta.dot(
                                 transformed.ring_geometries[a].normal),
                             x.Determinant() * za,
                             kGeometryAbs, kGeometryRel));
            EXPECT_TRUE(Near(moved_delta.dot(
                                 transformed.ring_geometries[b].normal),
                             x.Determinant() * zb,
                             kGeometryAbs, kGeometryRel));
        }
    }
}

bool AttachPlanar(ProteinConformation& conf) {
    auto geometry = nmr::GeometryResult::Compute(conf);
    if (!geometry || !conf.AttachResult(std::move(geometry))) return false;
    auto enrichment = nmr::EnrichmentResult::Compute(conf);
    if (!enrichment || !conf.AttachResult(std::move(enrichment))) return false;
    auto planar = nmr::PlanarGeometryResult::Compute(conf);
    return planar && conf.AttachResult(std::move(planar));
}

void CheckPlanar(const nmr::PlanarGeometryResult& original,
                 const nmr::PlanarGeometryResult& transformed,
                 const OrthogonalTransform& x) {
    ASSERT_EQ(transformed.OmegaActual().size(), original.OmegaActual().size());
    for (std::size_t i = 0; i < original.OmegaActual().size(); ++i) {
        const double a = original.OmegaActual()[i];
        const double b = transformed.OmegaActual()[i];
        if (!std::isfinite(a)) {
            EXPECT_FALSE(std::isfinite(b));
            continue;
        }
        const double expected = x.IsImproper() ? -a : a;
        EXPECT_NEAR(CircularDifference(b, expected), 0.0, 2.0e-11)
            << "omega_actual.npy residue " << i;
        const double da = original.OmegaDeviation()[i];
        const double db = transformed.OmegaDeviation()[i];
        const double expected_d = x.IsImproper() ? -da : da;
        EXPECT_NEAR(CircularDifference(db, expected_d), 0.0, 2.0e-11)
            << "omega_deviation.npy residue " << i;
    }

    ASSERT_EQ(transformed.AromaticChi2().size(),
              original.AromaticChi2().size());
    for (std::size_t i = 0; i < original.AromaticChi2().size(); ++i) {
        const double a = original.AromaticChi2()[i];
        const double b = transformed.AromaticChi2()[i];
        if (!std::isfinite(a)) {
            EXPECT_FALSE(std::isfinite(b));
            continue;
        }
        EXPECT_NEAR(CircularDifference(b, x.IsImproper() ? -a : a),
                    0.0, 2.0e-11)
            << "aromatic_chi2.npy ring " << i;
    }

    ASSERT_EQ(transformed.PuckerQ().size(), original.PuckerQ().size());
    for (std::size_t i = 0; i < original.PuckerQ().size(); ++i) {
        const double q = original.PuckerQ()[i];
        const double moved_q = transformed.PuckerQ()[i];
        if (!std::isfinite(q)) {
            EXPECT_FALSE(std::isfinite(moved_q));
            continue;
        }
        EXPECT_NEAR(moved_q, q, 2.0e-11) << "pucker_Q.npy ring " << i;
        const double theta = original.PuckerTheta()[i];
        const double moved_theta = transformed.PuckerTheta()[i];
        double expected = theta + (x.IsImproper() ? 180.0 : 0.0);
        expected = std::fmod(expected, 360.0);
        double difference = std::fmod(moved_theta - expected + 540.0, 360.0)
                          - 180.0;
        EXPECT_NEAR(difference, 0.0, 2.0e-9)
            << "pucker_theta.npy ring " << i;
    }
}

void AssignDeterministicCharges(ProteinConformation& conf) {
    conf.ForceAttachResultForTesting(
        std::make_unique<nmr::ChargeAssignmentResult>());
    for (std::size_t i = 0; i < conf.AtomCount(); ++i) {
        const int centered = static_cast<int>(i % 9) - 4;
        conf.MutableAtomAt(i).partial_charge = 0.0375 * centered;
    }
}

bool AttachElectrostatics(ProteinConformation& conf) {
    auto geometry = nmr::GeometryResult::Compute(conf);
    if (!geometry || !conf.AttachResult(std::move(geometry))) return false;
    AssignDeterministicCharges(conf);
    auto spatial = nmr::SpatialIndexResult::Compute(conf);
    if (!spatial || !conf.AttachResult(std::move(spatial))) return false;
    auto coulomb = nmr::CoulombResult::Compute(conf);
    if (!coulomb || !conf.AttachResult(std::move(coulomb))) return false;
    auto eeq = nmr::EeqResult::Compute(conf, 0);
    if (!eeq || !conf.AttachResult(std::move(eeq))) return false;
    auto eeq_coulomb = nmr::EeqCoulombResult::Compute(conf);
    return eeq_coulomb && conf.AttachResult(std::move(eeq_coulomb));
}

void CheckElectrostatics(const ProteinConformation& original,
                         const ProteinConformation& transformed,
                         const OrthogonalTransform& x) {
    for (std::size_t atom = 0; atom < original.AtomCount(); ++atom) {
        SCOPED_TRACE("electrostatic atom=" + std::to_string(atom));
        const auto& a = original.AtomAt(atom);
        const auto& b = transformed.AtomAt(atom);
        ExpectVector(b.coulomb_E_total, Polar(x, a.coulomb_E_total),
                     kElectrostaticAbs, kElectrostaticRel,
                     "coulomb_E.npy");
        ExpectVector(b.coulomb_E_backbone,
                     Polar(x, a.coulomb_E_backbone),
                     kElectrostaticAbs, kElectrostaticRel,
                     "coulomb_E_backbone.npy");
        ExpectVector(b.coulomb_E_sidechain,
                     Polar(x, a.coulomb_E_sidechain),
                     kElectrostaticAbs, kElectrostaticRel,
                     "coulomb_E_sidechain.npy");
        ExpectVector(b.coulomb_E_aromatic,
                     Polar(x, a.coulomb_E_aromatic),
                     kElectrostaticAbs, kElectrostaticRel,
                     "coulomb_E_aromatic.npy");
        ExpectEfg(b.coulomb_EFG_total_spherical,
                  a.coulomb_EFG_total_spherical, x,
                  "coulomb_efg.npy / coulomb_efg_t2.npy");
        ExpectEfg(b.coulomb_EFG_backbone_spherical,
                  a.coulomb_EFG_backbone_spherical, x,
                  "coulomb_efg_backbone.npy");
        ExpectEfg(b.coulomb_EFG_sidechain_spherical,
                  a.coulomb_EFG_sidechain_spherical, x,
                  "coulomb_efg_sidechain.npy");
        ExpectEfg(b.coulomb_EFG_aromatic_spherical,
                  a.coulomb_EFG_aromatic_spherical, x,
                  "coulomb_efg_aromatic.npy");

        EXPECT_NEAR(b.eeq_charge, a.eeq_charge, 3.0e-12);
        ExpectVector(b.eeq_coulomb_E_total,
                     Polar(x, a.eeq_coulomb_E_total),
                     kElectrostaticAbs, kElectrostaticRel,
                     "eeq_coulomb_E.npy");
        ExpectVector(b.eeq_coulomb_E_backbone,
                     Polar(x, a.eeq_coulomb_E_backbone),
                     kElectrostaticAbs, kElectrostaticRel,
                     "eeq_coulomb_E_backbone.npy");
        ExpectVector(b.eeq_coulomb_E_sidechain,
                     Polar(x, a.eeq_coulomb_E_sidechain),
                     kElectrostaticAbs, kElectrostaticRel,
                     "eeq_coulomb_E_sidechain.npy");
        ExpectVector(b.eeq_coulomb_E_aromatic,
                     Polar(x, a.eeq_coulomb_E_aromatic),
                     kElectrostaticAbs, kElectrostaticRel,
                     "eeq_coulomb_E_aromatic.npy");
        ExpectEfg(b.eeq_coulomb_EFG_total_spherical,
                  a.eeq_coulomb_EFG_total_spherical, x,
                  "eeq_coulomb_efg.npy");
        ExpectEfg(b.eeq_coulomb_EFG_backbone_spherical,
                  a.eeq_coulomb_EFG_backbone_spherical, x,
                  "eeq_coulomb_efg_backbone.npy");
        ExpectEfg(b.eeq_coulomb_EFG_sidechain_spherical,
                  a.eeq_coulomb_EFG_sidechain_spherical, x,
                  "eeq_coulomb_efg_sidechain.npy");
        ExpectEfg(b.eeq_coulomb_EFG_aromatic_spherical,
                  a.eeq_coulomb_EFG_aromatic_spherical, x,
                  "eeq_coulomb_efg_aromatic.npy");
    }
}

std::unique_ptr<Protein> BuildOneAtomProtein(const Vec3& position) {
    auto protein = std::make_unique<Protein>();
    nmr::Residue residue;
    residue.type = nmr::AminoAcid::ALA;
    residue.sequence_number = 1;
    residue.chain_id = "A";
    const std::size_t residue_index = protein->AddResidue(std::move(residue));
    auto atom = nmr::Atom::Create(nmr::Element::C);
    atom->residue_index = residue_index;
    const std::size_t atom_index = protein->AddAtom(std::move(atom));
    protein->MutableResidueAt(residue_index).atom_indices.push_back(atom_index);
    const std::vector<Vec3> positions{position};
    protein->FinalizeConstruction(positions);
    protein->AddConformation(positions, "directional one-atom source");
    return protein;
}

nmr::SolventEnvironment BuildAsymmetricSolvent() {
    nmr::SolventEnvironment solvent;
    auto add_water = [&](const Vec3& oxygen,
                         const Vec3& h1_offset,
                         const Vec3& h2_offset) {
        nmr::WaterMolecule water;
        water.O_pos = oxygen;
        water.H1_pos = oxygen + h1_offset;
        water.H2_pos = oxygen + h2_offset;
        water.O_charge = -0.834;
        water.H_charge = 0.417;
        solvent.waters.push_back(water);
        solvent.water_O_positions.push_back(oxygen);
    };
    add_water(Vec3(2.2, 0.4, -0.3),
              Vec3(0.72, 0.42, 0.11), Vec3(-0.31, 0.81, -0.09));
    add_water(Vec3(-1.4, 2.3, 0.7),
              Vec3(0.34, -0.74, 0.29), Vec3(-0.66, -0.21, 0.48));
    add_water(Vec3(3.8, -1.2, 0.9),
              Vec3(-0.44, 0.69, -0.22), Vec3(0.58, 0.51, 0.31));
    return solvent;
}

struct HBondFixture {
    std::unique_ptr<Protein> protein;
    std::size_t remote_target = SIZE_MAX;
};

HBondFixture BuildHBondFixture() {
    HBondFixture fixture;
    fixture.protein = std::make_unique<Protein>();
    for (int residue_index = 0; residue_index < 4; ++residue_index) {
        nmr::Residue residue;
        residue.type = nmr::AminoAcid::ALA;
        residue.sequence_number = residue_index + 1;
        residue.chain_id = residue_index < 2 ? "A" :
                           residue_index == 2 ? "B" : "C";
        fixture.protein->AddResidue(std::move(residue));
    }

    std::vector<Vec3> positions;
    auto add_atom = [&](std::size_t residue_index, const char* name,
                        nmr::Element element, const Vec3& position) {
        auto atom = nmr::Atom::Create(element);
        atom->pdb_atom_name = name;
        atom->residue_index = residue_index;
        const std::size_t atom_index =
            fixture.protein->AddAtom(std::move(atom));
        fixture.protein->MutableResidueAt(residue_index)
            .atom_indices.push_back(atom_index);
        positions.push_back(position);
        return atom_index;
    };
    add_atom(0, "C", nmr::Element::C, Vec3(-1.33, 0.0, 0.0));
    const std::size_t donor_n =
        add_atom(1, "N", nmr::Element::N, Vec3(0.0, 0.0, 0.0));
    const std::size_t donor_h =
        add_atom(1, "H", nmr::Element::H, Vec3(1.0, 0.0, 0.0));
    const std::size_t acceptor_o =
        add_atom(2, "O", nmr::Element::O, Vec3(1.0, 2.0, 0.0));
    fixture.remote_target =
        add_atom(3, "C", nmr::Element::C, Vec3(2.0, 2.0, 0.0));

    fixture.protein->MutableResidueAt(1).N = donor_n;
    fixture.protein->MutableResidueAt(1).H = donor_h;
    fixture.protein->MutableResidueAt(2).O = acceptor_o;
    fixture.protein->FinalizeConstruction(positions);
    fixture.protein->AddConformation(positions, "directional H-bond source");
    return fixture;
}

bool AttachHBond(ProteinConformation& conf) {
    if (!AttachGeometryAndSpatial(conf)) return false;
    std::vector<nmr::DsspResidue> residues(conf.ProteinRef().ResidueCount());
    residues[1].observed = true;
    residues[1].acceptors[0].residue_index = 2;
    auto dssp = nmr::DsspResult::CreateForTesting(std::move(residues));
    if (!dssp || !conf.AttachResult(std::move(dssp))) return false;
    auto hbond = nmr::HBondResult::Compute(conf);
    return hbond && conf.AttachResult(std::move(hbond));
}

std::unique_ptr<Protein> BuildSasaFixture() {
    auto protein = std::make_unique<Protein>();
    nmr::Residue residue;
    residue.type = nmr::AminoAcid::CYS;
    residue.sequence_number = 1;
    residue.chain_id = "A";
    const std::size_t ri = protein->AddResidue(std::move(residue));
    std::vector<Vec3> positions = {
        Vec3(0.0, 0.0, 0.0),
        Vec3(2.65, 0.37, -0.21),
        Vec3(-0.42, 2.83, 0.64),
    };
    const std::array<nmr::Element, 3> elements = {
        nmr::Element::C, nmr::Element::O, nmr::Element::S};
    for (std::size_t i = 0; i < positions.size(); ++i) {
        auto atom = nmr::Atom::Create(elements[i]);
        atom->residue_index = ri;
        const std::size_t ai = protein->AddAtom(std::move(atom));
        protein->MutableResidueAt(ri).atom_indices.push_back(ai);
    }
    protein->FinalizeConstruction(positions);
    protein->AddConformation(positions, "directional SASA source");
    return protein;
}

bool AttachSasa(ProteinConformation& conf) {
    if (!AttachGeometryAndSpatial(conf)) return false;
    auto sasa = nmr::SasaResult::Compute(conf);
    return sasa && conf.AttachResult(std::move(sasa));
}

bool AttachLocalBackbone(ProteinConformation& conf) {
    auto result = nmr::LocalBackboneGeometryResult::Compute(conf);
    return result && conf.AttachResult(std::move(result));
}

bool AttachSidechainCarbonylStack(ProteinConformation& conf) {
    if (!AttachGeometryAndSpatial(conf)) return false;
    auto mcconnell = nmr::McConnellResult::Compute(conf);
    if (!mcconnell || !conf.AttachResult(std::move(mcconnell))) return false;
    auto sidechain = nmr::SidechainCarbonylAnisotropyResult::Compute(conf);
    return sidechain && conf.AttachResult(std::move(sidechain));
}

bool AttachPiLocalStack(ProteinConformation& conf) {
    if (!AttachGeometryAndSpatial(conf)) return false;
    auto bs = nmr::BiotSavartResult::Compute(conf);
    if (!bs || !conf.AttachResult(std::move(bs))) return false;
    auto hm = nmr::HaighMallionResult::Compute(conf);
    if (!hm || !conf.AttachResult(std::move(hm))) return false;
    auto susceptibility = nmr::RingSusceptibilityResult::Compute(conf);
    if (!susceptibility ||
        !conf.AttachResult(std::move(susceptibility))) return false;
    auto quadrupole = nmr::PiQuadrupoleResult::Compute(conf);
    if (!quadrupole || !conf.AttachResult(std::move(quadrupole))) return false;
    auto dispersion = nmr::DispersionResult::Compute(conf);
    if (!dispersion || !conf.AttachResult(std::move(dispersion))) return false;
    auto local = nmr::PiQuadrupoleLocalTensorResult::Compute(conf);
    return local && conf.AttachResult(std::move(local));
}

nmr::SolventEnvironment BuildTypedWaterHBondSolvent(
        const ProteinConformation& conf,
        std::size_t* acceptor_atom,
        std::size_t* polar_h_atom) {
    nmr::SolventEnvironment solvent;
    *acceptor_atom = SIZE_MAX;
    *polar_h_atom = SIZE_MAX;
    const Protein& protein = conf.ProteinRef();
    const auto& topology = protein.LegacyAmber();
    if (!topology.HasAtomSemantic()) return solvent;

    for (std::size_t atom = 0; atom < protein.AtomCount(); ++atom) {
        const auto& semantic = topology.SemanticAt(atom);
        if (*acceptor_atom == SIZE_MAX &&
            nmr::water_hbond_geometry_detail::ProteinAcceptorClass(
                semantic) != 0) {
            *acceptor_atom = atom;
        }
        if (*polar_h_atom == SIZE_MAX && semantic.IsPolarH() &&
            protein.AtomAt(atom).parent_atom_index != SIZE_MAX &&
            protein.AtomAt(atom).parent_atom_index < protein.AtomCount()) {
            *polar_h_atom = atom;
        }
    }
    if (*acceptor_atom == SIZE_MAX || *polar_h_atom == SIZE_MAX)
        return solvent;

    auto append_water = [&](const Vec3& oxygen,
                            const Vec3& hydrogen1,
                            const Vec3& hydrogen2) {
        nmr::WaterMolecule water;
        water.O_pos = oxygen;
        water.H1_pos = hydrogen1;
        water.H2_pos = hydrogen2;
        water.O_charge = -0.834;
        water.H_charge = 0.417;
        solvent.waters.push_back(water);
        solvent.water_O_positions.push_back(oxygen);
    };

    // Water-donor row: H1 lies between O and the typed acceptor, with H2
    // deliberately asymmetric so production's nearest-H choice is stable.
    const Vec3 acceptor = conf.PositionAt(*acceptor_atom);
    const Vec3 u = Vec3(0.73, -0.41, 0.55).normalized();
    const Vec3 perpendicular =
        (Vec3(-0.21, 0.64, 0.39) -
         Vec3(-0.21, 0.64, 0.39).dot(u) * u).normalized();
    const Vec3 donor_oxygen = acceptor + 2.80 * u;
    append_water(donor_oxygen,
                 acceptor + 1.82 * u,
                 donor_oxygen + 0.96 * perpendicular);

    // Protein-donor row: water O extends the typed parent->H direction.
    // The emitted water-H columns for every such mode-2 row are intentional
    // NaNs and are checked after the production writer runs.
    const std::size_t parent =
        protein.AtomAt(*polar_h_atom).parent_atom_index;
    const Vec3 donor_h = conf.PositionAt(*polar_h_atom);
    const Vec3 donor_heavy = conf.PositionAt(parent);
    const Vec3 donor_axis = (donor_h - donor_heavy).normalized();
    const Vec3 acceptor_oxygen = donor_h + 1.75 * donor_axis;
    append_water(acceptor_oxygen,
                 acceptor_oxygen + 0.83 * perpendicular,
                 acceptor_oxygen - 0.37 * perpendicular + 0.72 * u);
    return solvent;
}

}  // namespace


TEST(DirectionalCovariance,
     GeometrySpatialRingBsHmPiQuadrupoleAndMcConnellRerunO3) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto protein = LoadUbqOrSkip();
    if (!protein) GTEST_SKIP() << "1UBQ protonated fixture unavailable";

    ProteinConformation& original = protein->Conformation();
    ASSERT_TRUE(AttachRingStack(original));

    for (const bool improper : {false, true}) {
        const OrthogonalTransform x =
            SeededTransform(0xD1AEC710ULL, improper);
        ASSERT_LE((x.Q.transpose() * x.Q - Mat3::Identity()).norm(),
                  2.0e-15);
        ASSERT_NEAR(std::abs(x.Determinant()), 1.0, 2.0e-15);
        ProteinConformation& moved = protein->AddConformation(
            nmr::test::directional::Positions(x, original.Positions()),
            improper ? "directional improper rerun" :
                       "directional proper rerun");
        ASSERT_TRUE(AttachRingStack(moved));
        CheckGeometrySpatialAndRingStack(original, moved, x);
    }
}


TEST(DirectionalCovariance, PlanarSignedGeometryRerunO3) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto protein = LoadUbqOrSkip();
    if (!protein) GTEST_SKIP() << "1UBQ protonated fixture unavailable";

    ProteinConformation& original = protein->Conformation();
    ASSERT_TRUE(AttachPlanar(original));
    const auto& original_result =
        original.Result<nmr::PlanarGeometryResult>();

    for (const bool improper : {false, true}) {
        const OrthogonalTransform x =
            SeededTransform(0x91A4A2ULL, improper);
        ProteinConformation& moved = protein->AddConformation(
            nmr::test::directional::Positions(x, original.Positions()),
            improper ? "planar improper rerun" : "planar proper rerun");
        ASSERT_TRUE(AttachPlanar(moved));
        CheckPlanar(original_result,
                    moved.Result<nmr::PlanarGeometryResult>(), x);
    }
}


TEST(DirectionalCovariance, CoulombAndEeqCoulombRerunO3WithNativeT2) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto protein = LoadUbqOrSkip();
    if (!protein) GTEST_SKIP() << "1UBQ protonated fixture unavailable";

    ProteinConformation& original = protein->Conformation();
    ASSERT_TRUE(AttachElectrostatics(original));
    const std::filesystem::path original_dir =
        FreshOutputDirectory("directional_coulomb_original");
    ASSERT_EQ(original.Result<nmr::CoulombResult>().WriteFeatures(
                  original, original_dir.string()), 12);

    for (const bool improper : {false, true}) {
        const OrthogonalTransform x =
            SeededTransform(0xE1EC7A05ULL, improper);
        ProteinConformation& moved = protein->AddConformation(
            nmr::test::directional::Positions(x, original.Positions()),
            improper ? "electrostatic improper rerun" :
                       "electrostatic proper rerun");
        ASSERT_TRUE(AttachElectrostatics(moved));
        CheckElectrostatics(original, moved, x);

        const std::filesystem::path moved_dir = FreshOutputDirectory(
            improper ? "directional_coulomb_improper" :
                       "directional_coulomb_proper");
        ASSERT_EQ(moved.Result<nmr::CoulombResult>().WriteFeatures(
                      moved, moved_dir.string()), 12);
        for (const char* filename : {
                 "coulomb_E.npy", "coulomb_E_backbone.npy",
                 "coulomb_E_sidechain.npy", "coulomb_E_aromatic.npy"}) {
            CheckSerializedPolarNpy(
                original_dir, moved_dir, filename, original.AtomCount(), x,
                kElectrostaticAbs, kElectrostaticRel);
        }
        CheckSerializedFull9Npy(
            original_dir, moved_dir, "coulomb_efg.npy",
            original.AtomCount(), x, kElectrostaticAbs,
            kElectrostaticRel, true);
        for (const char* filename : {
                 "coulomb_efg_t2.npy", "coulomb_efg_backbone.npy",
                 "coulomb_efg_sidechain.npy",
                 "coulomb_efg_aromatic.npy"}) {
            CheckSerializedNativeT2Npy(
                original_dir, moved_dir, filename, original.AtomCount(), x,
                kElectrostaticAbs, kElectrostaticRel);
        }
        RemoveOutputDirectory(moved_dir);
    }
    RemoveOutputDirectory(original_dir);
}


TEST(DirectionalCovariance,
     WaterFieldAndHydrationGeometryRerunO3WithTransformedSolvent) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto protein = BuildOneAtomProtein(Vec3(0.35, -0.27, 0.19));
    ProteinConformation& original = protein->Conformation();
    ASSERT_TRUE(AttachGeometryAndSpatial(original));
    const nmr::SolventEnvironment solvent = BuildAsymmetricSolvent();
    auto water = nmr::WaterFieldResult::Compute(original, solvent);
    ASSERT_NE(water, nullptr);
    ASSERT_TRUE(original.AttachResult(std::move(water)));
    auto sasa = nmr::SasaResult::Compute(original);
    ASSERT_NE(sasa, nullptr);
    ASSERT_TRUE(original.AttachResult(std::move(sasa)));
    // HydrationGeometryResult owns the production SasaResult dependency.
    // Replace only its sampled normal with a deterministic non-axis vector so
    // this test isolates the hydration algebra from SASA's separate finite
    // Fibonacci-lattice covariance bound.
    original.MutableAtomAt(0).sasa_normal =
        Vec3(0.31, -0.44, 0.842).normalized();
    auto hydration =
        nmr::HydrationGeometryResult::Compute(original, solvent);
    ASSERT_NE(hydration, nullptr);
    ASSERT_TRUE(original.AttachResult(std::move(hydration)));
    const std::filesystem::path original_dir =
        FreshOutputDirectory("directional_water_original");
    ASSERT_EQ(original.Result<nmr::WaterFieldResult>().WriteFeatures(
                  original, original_dir.string()), 9);
    ASSERT_EQ(original.Result<nmr::HydrationGeometryResult>().WriteFeatures(
                  original, original_dir.string()), 1);
    const std::array<std::string, 5> water_companion_names = {
        "water_shell_counts.npy",
        "water_efield_clamp_mask.npy",
        "water_efield_clamp_scale.npy",
        "water_efield_first_clamp_mask.npy",
        "water_efield_first_clamp_scale.npy",
    };
    std::vector<SerializedNpy> water_companions_a;
    water_companions_a.reserve(water_companion_names.size());
    for (const auto& name : water_companion_names) {
        water_companions_a.push_back(ReadSerializedNpy(original_dir / name));
    }
    const auto shell_counts =
        NpyValues<double>(water_companions_a.front(), "<f8");
    ASSERT_TRUE(std::any_of(shell_counts.begin(), shell_counts.end(),
                            [](double value) { return value > 0.0; }))
        << "water_shell_counts.npy proof must include a nonempty shell";

    for (const bool improper : {false, true}) {
        const OrthogonalTransform x =
            SeededTransform(0x50A7E17ULL, improper);
        ProteinConformation& moved = protein->AddConformation(
            nmr::test::directional::Positions(x, original.Positions()),
            improper ? "water improper rerun" : "water proper rerun");
        ASSERT_TRUE(AttachGeometryAndSpatial(moved));
        const nmr::SolventEnvironment moved_solvent =
            nmr::test::directional::Solvent(x, solvent);
        auto moved_water =
            nmr::WaterFieldResult::Compute(moved, moved_solvent);
        ASSERT_NE(moved_water, nullptr);
        ASSERT_TRUE(moved.AttachResult(std::move(moved_water)));
        auto moved_sasa = nmr::SasaResult::Compute(moved);
        ASSERT_NE(moved_sasa, nullptr);
        ASSERT_TRUE(moved.AttachResult(std::move(moved_sasa)));
        moved.MutableAtomAt(0).sasa_normal =
            Polar(x, original.AtomAt(0).sasa_normal);
        auto moved_hydration =
            nmr::HydrationGeometryResult::Compute(moved, moved_solvent);
        ASSERT_NE(moved_hydration, nullptr);
        ASSERT_TRUE(moved.AttachResult(std::move(moved_hydration)));

        const auto& a = original.AtomAt(0);
        const auto& b = moved.AtomAt(0);
        ExpectVector(b.water_efield, Polar(x, a.water_efield),
                     kElectrostaticAbs, kElectrostaticRel,
                     "water_efield.npy");
        ExpectVector(b.water_efield_first,
                     Polar(x, a.water_efield_first),
                     kElectrostaticAbs, kElectrostaticRel,
                     "water_efield_first.npy");
        ExpectEfg(b.water_efg_spherical, a.water_efg_spherical, x,
                  "water_efg.npy");
        ExpectEfg(b.water_efg_first_spherical,
                  a.water_efg_first_spherical, x,
                  "water_efg_first.npy");
        ExpectVector(b.water_dipole_vector,
                     Polar(x, a.water_dipole_vector),
                     kGeometryAbs, kGeometryRel,
                     "water_polarization[:,0:3]");
        ExpectVector(b.water_surface_normal,
                     Polar(x, a.water_surface_normal),
                     kGeometryAbs, kGeometryRel,
                     "water_polarization[:,3:6]");
        EXPECT_NEAR(b.sasa_half_shell_asymmetry,
                    a.sasa_half_shell_asymmetry, 1.0e-12);
        EXPECT_NEAR(b.sasa_dipole_alignment,
                    a.sasa_dipole_alignment, 1.0e-12);

        const std::filesystem::path moved_dir = FreshOutputDirectory(
            improper ? "directional_water_improper" :
                       "directional_water_proper");
        ASSERT_EQ(moved.Result<nmr::WaterFieldResult>().WriteFeatures(
                      moved, moved_dir.string()), 9);
        ASSERT_EQ(moved.Result<nmr::HydrationGeometryResult>().WriteFeatures(
                      moved, moved_dir.string()), 1);
        for (std::size_t companion = 0;
             companion < water_companion_names.size(); ++companion) {
            const SerializedNpy moved_companion = ReadSerializedNpy(
                moved_dir / water_companion_names[companion]);
            EXPECT_EQ(moved_companion.descr,
                      water_companions_a[companion].descr)
                << water_companion_names[companion];
            EXPECT_EQ(moved_companion.shape,
                      water_companions_a[companion].shape)
                << water_companion_names[companion];
            EXPECT_EQ(moved_companion.payload,
                      water_companions_a[companion].payload)
                << water_companion_names[companion];
        }
        for (const char* filename : {
                 "water_efield.npy", "water_efield_first.npy"}) {
            CheckSerializedPolarNpy(
                original_dir, moved_dir, filename, original.AtomCount(), x,
                kElectrostaticAbs, kElectrostaticRel);
        }
        for (const char* filename : {
                 "water_efg.npy", "water_efg_first.npy"}) {
            CheckSerializedNativeT2Npy(
                original_dir, moved_dir, filename, original.AtomCount(), x,
                kElectrostaticAbs, kElectrostaticRel);
        }
        const SerializedNpy polarization_a_npy = ReadSerializedNpy(
            original_dir / "water_polarization.npy");
        const SerializedNpy polarization_b_npy = ReadSerializedNpy(
            moved_dir / "water_polarization.npy");
        ASSERT_EQ(polarization_a_npy.shape,
                  (std::vector<std::size_t>{original.AtomCount(), 10u}));
        ASSERT_EQ(polarization_b_npy.shape, polarization_a_npy.shape);
        const auto polarization_a =
            NpyValues<double>(polarization_a_npy, "<f8");
        const auto polarization_b =
            NpyValues<double>(polarization_b_npy, "<f8");
        for (std::size_t atom = 0; atom < original.AtomCount(); ++atom) {
            for (std::size_t offset : {0u, 3u}) {
                const Vec3 source(
                    polarization_a[atom * 10 + offset],
                    polarization_a[atom * 10 + offset + 1],
                    polarization_a[atom * 10 + offset + 2]);
                const Vec3 emitted(
                    polarization_b[atom * 10 + offset],
                    polarization_b[atom * 10 + offset + 1],
                    polarization_b[atom * 10 + offset + 2]);
                ExpectVector(emitted, Polar(x, source),
                             kGeometryAbs, kGeometryRel,
                             "water_polarization.npy");
            }
            for (std::size_t column = 6; column < 10; ++column) {
                EXPECT_TRUE(Near(polarization_b[atom * 10 + column],
                                 polarization_a[atom * 10 + column],
                                 kGeometryAbs, kGeometryRel))
                    << "water_polarization.npy scalar column=" << column;
            }
        }
        RemoveOutputDirectory(moved_dir);
    }
    RemoveOutputDirectory(original_dir);
}


TEST(DirectionalCovariance, HBondDirectionsRerunO3) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    HBondFixture fixture = BuildHBondFixture();
    ProteinConformation& original = fixture.protein->Conformation();
    ASSERT_TRUE(AttachHBond(original));
    ASSERT_GT(original.AtomAt(fixture.remote_target).hbond_nearest_dist, 0.0);

    for (const bool improper : {false, true}) {
        const OrthogonalTransform x =
            SeededTransform(0xAB01D5ULL, improper);
        ProteinConformation& moved = fixture.protein->AddConformation(
            nmr::test::directional::Positions(x, original.Positions()),
            improper ? "hbond improper rerun" : "hbond proper rerun");
        ASSERT_TRUE(AttachHBond(moved));
        for (std::size_t atom = 0; atom < original.AtomCount(); ++atom) {
            const auto& a = original.AtomAt(atom);
            const auto& b = moved.AtomAt(atom);
            EXPECT_TRUE(Near(b.hbond_nearest_dist, a.hbond_nearest_dist,
                             kGeometryAbs, kGeometryRel));
            EXPECT_TRUE(Near(b.hbond_mcconnell_scalar,
                             a.hbond_mcconnell_scalar,
                             kKernelAbs, kKernelRel));
            ExpectVector(b.hbond_nearest_dir,
                         Polar(x, a.hbond_nearest_dir),
                         kGeometryAbs, kGeometryRel,
                         "hbond_nearest_dir.npy");
        }
    }
}


TEST(DirectionalCovariance,
     SasaRerunO3BoundedByRecordedFibonacciDiscretizationError) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto protein = BuildSasaFixture();
    ProteinConformation& original = protein->Conformation();
    ASSERT_TRUE(AttachSasa(original));

    // SasaResult samples a finite, lab-fixed Fibonacci lattice.  It is not
    // algebraically O(3)-closed at finite n_points; these bounds record that
    // named quadrature cause rather than pretending machine precision.  The
    // area bound is 4.0 A^2: on this 92-point fixture the observed worst case
    // is exactly three exposed-sample cells (3.4938881219436553 A^2 for O),
    // while one cell is about 1.16 A^2.  The 0.15 unit-normal bound
    // corresponds to about 8.6 degrees.  Both numbers therefore name the
    // finite quadrature cause; neither is a loosened analytic-kernel bound.
    constexpr double kSasaAreaAbsA2 = 4.0;
    constexpr double kSasaAreaRel = 0.0;
    constexpr double kSasaNormalVectorError = 0.15;

    for (const bool improper : {false, true}) {
        const OrthogonalTransform x =
            SeededTransform(0x5A5A2026ULL, improper);
        ProteinConformation& moved = protein->AddConformation(
            nmr::test::directional::Positions(x, original.Positions()),
            improper ? "sasa improper rerun" : "sasa proper rerun");
        ASSERT_TRUE(AttachSasa(moved));
        for (std::size_t atom = 0; atom < original.AtomCount(); ++atom) {
            const auto& a = original.AtomAt(atom);
            const auto& b = moved.AtomAt(atom);
            EXPECT_TRUE(Near(b.atom_sasa, a.atom_sasa,
                             kSasaAreaAbsA2, kSasaAreaRel))
                << "atom_sasa.npy fixed-grid discrepancy at atom " << atom
                << " actual=" << b.atom_sasa
                << " expected=" << a.atom_sasa
                << " abs_error=" << std::abs(b.atom_sasa - a.atom_sasa)
                << " tolerance="
                << kSasaAreaAbsA2 + kSasaAreaRel * std::abs(a.atom_sasa);
            if (a.sasa_normal.norm() > 0.5 && b.sasa_normal.norm() > 0.5) {
                EXPECT_LE((b.sasa_normal - Polar(x, a.sasa_normal)).norm(),
                          kSasaNormalVectorError)
                    << "sasa_normal.npy is polar by its radial-sum Compute "
                       "body; finite Fibonacci sampling sets this tolerance";
            }
        }
    }
}


TEST(DirectionalCovariance,
     LocalBackboneCbResidualVectorRerunProperWithSerializedMask) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto protein = LoadUbqOrSkip();
    if (!protein) GTEST_SKIP() << "1UBQ protonated fixture unavailable";

    ProteinConformation& original = protein->Conformation();
    ASSERT_TRUE(AttachLocalBackbone(original));
    const std::filesystem::path original_dir =
        FreshOutputDirectory("directional_local_backbone_original");
    ASSERT_EQ(original.Result<nmr::LocalBackboneGeometryResult>()
                  .WriteFeatures(original, original_dir.string()),
              14);
    const SerializedNpy original_vector = ReadSerializedNpy(
        original_dir / "cb_residual_vector.npy");
    const SerializedNpy original_valid = ReadSerializedNpy(
        original_dir / "cb_residual_vector_valid.npy");
    const std::vector<double> a = NpyValues<double>(original_vector, "<f8");
    const std::vector<std::uint8_t> a_valid =
        NpyValues<std::uint8_t>(original_valid, "|u1");
    ASSERT_EQ(original_vector.shape,
              (std::vector<std::size_t>{protein->ResidueCount(), 3}));
    ASSERT_EQ(original_valid.shape,
              (std::vector<std::size_t>{protein->ResidueCount()}));

    // The ideal-CB kernel contains c cross n and therefore encodes the
    // protein's handedness.  Its ordinary polar law is valid for proper
    // rigid transforms; an improper transform is deliberately not claimed.
    const OrthogonalTransform x =
        SeededTransform(0xCBA11BACULL, false);
    ProteinConformation& moved = protein->AddConformation(
        nmr::test::directional::Positions(x, original.Positions()),
        "local-backbone proper rerun");
    ASSERT_TRUE(AttachLocalBackbone(moved));
    const std::filesystem::path moved_dir =
        FreshOutputDirectory("directional_local_backbone_moved");
    ASSERT_EQ(moved.Result<nmr::LocalBackboneGeometryResult>()
                  .WriteFeatures(moved, moved_dir.string()),
              14);
    const SerializedNpy moved_vector = ReadSerializedNpy(
        moved_dir / "cb_residual_vector.npy");
    const SerializedNpy moved_valid = ReadSerializedNpy(
        moved_dir / "cb_residual_vector_valid.npy");
    const std::vector<double> b = NpyValues<double>(moved_vector, "<f8");
    const std::vector<std::uint8_t> b_valid =
        NpyValues<std::uint8_t>(moved_valid, "|u1");
    ASSERT_EQ(moved_vector.shape, original_vector.shape);
    ASSERT_EQ(moved_valid.shape, original_valid.shape);
    ASSERT_EQ(b_valid, a_valid)
        << "cb_residual_vector_valid.npy row mask changed";

    std::size_t finite_rows = 0;
    for (std::size_t residue = 0; residue < a_valid.size(); ++residue) {
        SCOPED_TRACE("residue=" + std::to_string(residue));
        const Vec3 source(a[residue * 3], a[residue * 3 + 1],
                          a[residue * 3 + 2]);
        const Vec3 actual(b[residue * 3], b[residue * 3 + 1],
                          b[residue * 3 + 2]);
        if (a_valid[residue] == 0) {
            EXPECT_FALSE(source.allFinite());
            EXPECT_FALSE(actual.allFinite());
            for (int component = 0; component < 3; ++component)
                EXPECT_EQ(std::isnan(actual[component]),
                          std::isnan(source[component]));
            continue;
        }
        ++finite_rows;
        ASSERT_TRUE(source.allFinite());
        ASSERT_TRUE(actual.allFinite());
        ExpectVector(actual, Polar(x, source),
                     kGeometryAbs, kGeometryRel,
                     "cb_residual_vector.npy");
    }
    EXPECT_GT(finite_rows, 0u);
    RemoveOutputDirectory(original_dir);
    RemoveOutputDirectory(moved_dir);
}


TEST(DirectionalCovariance,
     SidechainCarbonylFrameAndFixedTensorRerunO3Serialized) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto protein = LoadUbqOrSkip();
    if (!protein) GTEST_SKIP() << "1UBQ protonated fixture unavailable";

    ProteinConformation& original = protein->Conformation();
    ASSERT_TRUE(AttachSidechainCarbonylStack(original));
    const std::filesystem::path original_dir =
        FreshOutputDirectory("directional_sidechain_co_original");
    ASSERT_EQ(original.Result<nmr::SidechainCarbonylAnisotropyResult>()
                  .WriteFeatures(original, original_dir.string()),
              6);
    const SerializedNpy source_rows_a = ReadSerializedNpy(
        original_dir / "sidechain_co_source_bonds.npy");
    const SerializedNpy frames_a_npy = ReadSerializedNpy(
        original_dir / "sidechain_co_frame.npy");
    const SerializedNpy quality_a_npy = ReadSerializedNpy(
        original_dir / "sidechain_co_frame_quality.npy");
    const SerializedNpy fixed_a_npy = ReadSerializedNpy(
        original_dir / "sidechain_co_fixed_T2.npy");
    const SerializedNpy bo_a_npy = ReadSerializedNpy(
        original_dir / "sidechain_co_bo_T2.npy");
    const auto source_a = NpyValues<std::int32_t>(source_rows_a, "<i4");
    const auto frames_a = NpyValues<double>(frames_a_npy, "<f8");
    const auto quality_a = NpyValues<double>(quality_a_npy, "<f8");
    const auto fixed_a = NpyValues<double>(fixed_a_npy, "<f8");
    const auto bo_a = NpyValues<double>(bo_a_npy, "<f8");
    ASSERT_EQ(source_rows_a.shape.size(), 2u);
    ASSERT_EQ(source_rows_a.shape[1], 8u);
    const std::size_t source_count = source_rows_a.shape[0];
    ASSERT_GT(source_count, 0u);
    ASSERT_EQ(frames_a_npy.shape,
              (std::vector<std::size_t>{source_count, 12}));
    ASSERT_EQ(quality_a_npy.shape,
              (std::vector<std::size_t>{source_count, 4}));
    ASSERT_EQ(fixed_a_npy.shape,
              (std::vector<std::size_t>{original.AtomCount(), 9}));
    ASSERT_EQ(bo_a_npy.shape, fixed_a_npy.shape);

    for (const bool improper : {false, true}) {
        const OrthogonalTransform x =
            SeededTransform(0x51DEC0A2ULL, improper);
        ProteinConformation& moved = protein->AddConformation(
            nmr::test::directional::Positions(x, original.Positions()),
            improper ? "sidechain CO improper rerun" :
                       "sidechain CO proper rerun");
        ASSERT_TRUE(AttachSidechainCarbonylStack(moved));
        const std::filesystem::path moved_dir = FreshOutputDirectory(
            improper ? "directional_sidechain_co_improper" :
                       "directional_sidechain_co_proper");
        ASSERT_EQ(moved.Result<nmr::SidechainCarbonylAnisotropyResult>()
                      .WriteFeatures(moved, moved_dir.string()),
                  6);

        const SerializedNpy source_rows_b = ReadSerializedNpy(
            moved_dir / "sidechain_co_source_bonds.npy");
        const SerializedNpy frames_b_npy = ReadSerializedNpy(
            moved_dir / "sidechain_co_frame.npy");
        const SerializedNpy quality_b_npy = ReadSerializedNpy(
            moved_dir / "sidechain_co_frame_quality.npy");
        const SerializedNpy fixed_b_npy = ReadSerializedNpy(
            moved_dir / "sidechain_co_fixed_T2.npy");
        const SerializedNpy bo_b_npy = ReadSerializedNpy(
            moved_dir / "sidechain_co_bo_T2.npy");
        const auto source_b =
            NpyValues<std::int32_t>(source_rows_b, "<i4");
        const auto frames_b = NpyValues<double>(frames_b_npy, "<f8");
        const auto quality_b = NpyValues<double>(quality_b_npy, "<f8");
        const auto fixed_b = NpyValues<double>(fixed_b_npy, "<f8");
        const auto bo_b = NpyValues<double>(bo_b_npy, "<f8");
        ASSERT_EQ(source_rows_b.shape, source_rows_a.shape);
        ASSERT_EQ(frames_b_npy.shape, frames_a_npy.shape);
        ASSERT_EQ(quality_b_npy.shape, quality_a_npy.shape);
        ASSERT_EQ(fixed_b_npy.shape, fixed_a_npy.shape);
        ASSERT_EQ(bo_b_npy.shape, bo_a_npy.shape);
        EXPECT_EQ(source_b, source_a)
            << "sidechain_co_source_bonds.npy row identity changed";

        std::size_t valid_frames = 0;
        for (std::size_t row = 0; row < source_count; ++row) {
            SCOPED_TRACE("sidechain source row=" + std::to_string(row));
            const bool valid = quality_a[row * 4 + 3] == 1.0;
            EXPECT_EQ(quality_b[row * 4 + 3], quality_a[row * 4 + 3]);
            for (int component = 0; component < 3; ++component) {
                ExpectScalarOrSameNaN(
                    quality_b[row * 4 + component],
                    quality_a[row * 4 + component],
                    kGeometryAbs, kGeometryRel,
                    "sidechain_co_frame_quality.npy");
            }

            const Vec3 origin_a(frames_a[row * 12],
                                frames_a[row * 12 + 1],
                                frames_a[row * 12 + 2]);
            const Vec3 origin_b(frames_b[row * 12],
                                frames_b[row * 12 + 1],
                                frames_b[row * 12 + 2]);
            if (origin_a.allFinite()) {
                ExpectVector(origin_b, Position(x, origin_a),
                             kGeometryAbs, kGeometryRel,
                             "sidechain_co_frame.npy origin");
            } else {
                EXPECT_FALSE(origin_b.allFinite());
            }
            if (!valid) {
                for (std::size_t component = 3; component < 12;
                     ++component) {
                    EXPECT_EQ(std::isnan(frames_b[row * 12 + component]),
                              std::isnan(frames_a[row * 12 + component]));
                }
                continue;
            }
            ++valid_frames;
            const Vec3 x_a(frames_a[row * 12 + 3],
                           frames_a[row * 12 + 4],
                           frames_a[row * 12 + 5]);
            const Vec3 y_a(frames_a[row * 12 + 6],
                           frames_a[row * 12 + 7],
                           frames_a[row * 12 + 8]);
            const Vec3 z_a(frames_a[row * 12 + 9],
                           frames_a[row * 12 + 10],
                           frames_a[row * 12 + 11]);
            const Vec3 x_b(frames_b[row * 12 + 3],
                           frames_b[row * 12 + 4],
                           frames_b[row * 12 + 5]);
            const Vec3 y_b(frames_b[row * 12 + 6],
                           frames_b[row * 12 + 7],
                           frames_b[row * 12 + 8]);
            const Vec3 z_b(frames_b[row * 12 + 9],
                           frames_b[row * 12 + 10],
                           frames_b[row * 12 + 11]);
            ExpectVector(x_b, Polar(x, x_a),
                         kGeometryAbs, kGeometryRel,
                         "sidechain_co_frame.npy x axis");
            ExpectVector(y_b, Polar(x, y_a),
                         kGeometryAbs, kGeometryRel,
                         "sidechain_co_frame.npy y axis");
            ExpectVector(z_b, Axial(x, z_a),
                         kGeometryAbs, kGeometryRel,
                         "sidechain_co_frame.npy z axis");
        }
        EXPECT_GT(valid_frames, 0u);

        for (std::size_t atom = 0; atom < original.AtomCount(); ++atom) {
            SCOPED_TRACE("sidechain tensor atom=" + std::to_string(atom));
            const SphericalTensor fixed_source =
                UnpackFull9(&fixed_a[atom * 9]);
            const SphericalTensor fixed_actual =
                UnpackFull9(&fixed_b[atom * 9]);
            ExpectEvenSphericalTensor(
                fixed_actual, fixed_source, x,
                kKernelAbs, kKernelRel,
                "sidechain_co_fixed_T2.npy (PackFull9)");
            for (std::size_t component = 0; component < 9; ++component) {
                EXPECT_TRUE(std::isnan(bo_a[atom * 9 + component]));
                EXPECT_TRUE(std::isnan(bo_b[atom * 9 + component]))
                    << "sidechain_co_bo_T2.npy must retain explicit "
                       "MOPAC-unavailable NaN";
            }
        }
        RemoveOutputDirectory(moved_dir);
    }
    RemoveOutputDirectory(original_dir);
}


TEST(DirectionalCovariance,
     PiQuadrupoleLocalMixedParityRerunO3Serialized) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto protein = LoadUbqOrSkip();
    if (!protein) GTEST_SKIP() << "1UBQ protonated fixture unavailable";

    ProteinConformation& original = protein->Conformation();
    ASSERT_TRUE(AttachPiLocalStack(original));
    const std::filesystem::path original_dir =
        FreshOutputDirectory("directional_piquad_local_original");
    ASSERT_EQ(original.Result<nmr::PiQuadrupoleLocalTensorResult>()
                  .WriteFeatures(original, original_dir.string()),
              4);
    const SerializedNpy full_a_npy = ReadSerializedNpy(
        original_dir / "piquad_local_tensor.npy");
    const SerializedNpy t2_a_npy = ReadSerializedNpy(
        original_dir / "piquad_local_T2.npy");
    const SerializedNpy frame_a_npy = ReadSerializedNpy(
        original_dir / "piquad_local_frame.npy");
    const SerializedNpy geometry_a_npy = ReadSerializedNpy(
        original_dir / "piquad_local_geometry.npy");
    const auto full_a = NpyValues<double>(full_a_npy, "<f8");
    const auto t2_a = NpyValues<double>(t2_a_npy, "<f8");
    const auto frame_a = NpyValues<double>(frame_a_npy, "<f8");
    const auto geometry_a = NpyValues<double>(geometry_a_npy, "<f8");
    ASSERT_EQ(full_a_npy.shape.size(), 2u);
    ASSERT_EQ(full_a_npy.shape[1], 9u);
    const std::size_t pair_count = full_a_npy.shape[0];
    ASSERT_GT(pair_count, 0u);
    ASSERT_EQ(t2_a_npy.shape,
              (std::vector<std::size_t>{pair_count, 5}));
    ASSERT_EQ(frame_a_npy.shape,
              (std::vector<std::size_t>{pair_count, 9}));
    ASSERT_EQ(geometry_a_npy.shape,
              (std::vector<std::size_t>{pair_count, 8}));

    for (const bool improper : {false, true}) {
        const OrthogonalTransform x =
            SeededTransform(0x91A0CA11ULL, improper);
        ProteinConformation& moved = protein->AddConformation(
            nmr::test::directional::Positions(x, original.Positions()),
            improper ? "pi local improper rerun" :
                       "pi local proper rerun");
        ASSERT_TRUE(AttachPiLocalStack(moved));
        const std::filesystem::path moved_dir = FreshOutputDirectory(
            improper ? "directional_piquad_local_improper" :
                       "directional_piquad_local_proper");
        ASSERT_EQ(moved.Result<nmr::PiQuadrupoleLocalTensorResult>()
                      .WriteFeatures(moved, moved_dir.string()),
                  4);
        const SerializedNpy full_b_npy = ReadSerializedNpy(
            moved_dir / "piquad_local_tensor.npy");
        const SerializedNpy t2_b_npy = ReadSerializedNpy(
            moved_dir / "piquad_local_T2.npy");
        const SerializedNpy frame_b_npy = ReadSerializedNpy(
            moved_dir / "piquad_local_frame.npy");
        const SerializedNpy geometry_b_npy = ReadSerializedNpy(
            moved_dir / "piquad_local_geometry.npy");
        const auto full_b = NpyValues<double>(full_b_npy, "<f8");
        const auto t2_b = NpyValues<double>(t2_b_npy, "<f8");
        const auto frame_b = NpyValues<double>(frame_b_npy, "<f8");
        const auto geometry_b = NpyValues<double>(geometry_b_npy, "<f8");
        ASSERT_EQ(full_b_npy.shape, full_a_npy.shape);
        ASSERT_EQ(t2_b_npy.shape, t2_a_npy.shape);
        ASSERT_EQ(frame_b_npy.shape, frame_a_npy.shape);
        ASSERT_EQ(geometry_b_npy.shape, geometry_a_npy.shape);

        Mat3 local_parity = Mat3::Identity();
        local_parity(2, 2) = x.Determinant();
        std::size_t valid_rows = 0;
        for (std::size_t row = 0; row < pair_count; ++row) {
            SCOPED_TRACE("pi local row=" + std::to_string(row));
            for (int identity = 0; identity < 3; ++identity)
                EXPECT_DOUBLE_EQ(geometry_b[row * 8 + identity],
                                 geometry_a[row * 8 + identity]);
            EXPECT_DOUBLE_EQ(geometry_b[row * 8 + 6],
                             geometry_a[row * 8 + 6]);
            EXPECT_DOUBLE_EQ(geometry_b[row * 8 + 7],
                             geometry_a[row * 8 + 7]);
            ExpectScalarOrSameNaN(geometry_b[row * 8 + 3],
                                  geometry_a[row * 8 + 3],
                                  kGeometryAbs, kGeometryRel,
                                  "piquad_local_geometry.npy distance");
            if (std::isnan(geometry_a[row * 8 + 4])) {
                EXPECT_TRUE(std::isnan(geometry_b[row * 8 + 4]));
            } else {
                EXPECT_TRUE(Near(geometry_b[row * 8 + 4],
                                 x.Determinant() *
                                     geometry_a[row * 8 + 4],
                                 kGeometryAbs, kGeometryRel))
                    << "piquad_local_geometry.npy cos_theta pseudoscalar";
            }
            ExpectScalarOrSameNaN(geometry_b[row * 8 + 5],
                                  geometry_a[row * 8 + 5],
                                  kKernelAbs, kKernelRel,
                                  "piquad_local_geometry.npy quad scalar");

            const bool valid = geometry_a[row * 8 + 6] == 1.0;
            if (!valid) {
                for (std::size_t component = 0; component < 9; ++component) {
                    EXPECT_EQ(std::isnan(full_b[row * 9 + component]),
                              std::isnan(full_a[row * 9 + component]));
                    EXPECT_EQ(std::isnan(frame_b[row * 9 + component]),
                              std::isnan(frame_a[row * 9 + component]));
                }
                for (std::size_t component = 0; component < 5; ++component)
                    EXPECT_EQ(std::isnan(t2_b[row * 5 + component]),
                              std::isnan(t2_a[row * 5 + component]));
                continue;
            }
            ++valid_rows;
            const Vec3 local_x_a(frame_a[row * 9],
                                 frame_a[row * 9 + 1],
                                 frame_a[row * 9 + 2]);
            const Vec3 local_y_a(frame_a[row * 9 + 3],
                                 frame_a[row * 9 + 4],
                                 frame_a[row * 9 + 5]);
            const Vec3 local_z_a(frame_a[row * 9 + 6],
                                 frame_a[row * 9 + 7],
                                 frame_a[row * 9 + 8]);
            const Vec3 local_x_b(frame_b[row * 9],
                                 frame_b[row * 9 + 1],
                                 frame_b[row * 9 + 2]);
            const Vec3 local_y_b(frame_b[row * 9 + 3],
                                 frame_b[row * 9 + 4],
                                 frame_b[row * 9 + 5]);
            const Vec3 local_z_b(frame_b[row * 9 + 6],
                                 frame_b[row * 9 + 7],
                                 frame_b[row * 9 + 8]);
            ExpectVector(local_x_b, Polar(x, local_x_a),
                         kGeometryAbs, kGeometryRel,
                         "piquad_local_frame.npy x axis");
            ExpectVector(local_y_b, Polar(x, local_y_a),
                         kGeometryAbs, kGeometryRel,
                         "piquad_local_frame.npy y axis");
            ExpectVector(local_z_b, Axial(x, local_z_a),
                         kGeometryAbs, kGeometryRel,
                         "piquad_local_frame.npy z axis");

            const SphericalTensor source =
                UnpackFull9(&full_a[row * 9]);
            const SphericalTensor actual =
                UnpackFull9(&full_b[row * 9]);
            const Mat3 expected_local =
                local_parity * source.Reconstruct() * local_parity;
            ExpectMatrix(actual.Reconstruct(), expected_local,
                         kKernelAbs, kKernelRel,
                         "piquad_local_tensor.npy mixed local parity");
            EXPECT_LE(std::abs(source.T0), kKernelAbs);
            EXPECT_LE(T1Vector(source).norm(), kKernelAbs);
            EXPECT_LE(std::abs(actual.T0), kKernelAbs);
            EXPECT_LE(T1Vector(actual).norm(), kKernelAbs);
            const SphericalTensor expected_native =
                SphericalTensor::Decompose(expected_local);
            for (int component = 0; component < 5; ++component) {
                EXPECT_TRUE(Near(t2_a[row * 5 + component],
                                 source.T2[component],
                                 kKernelAbs, kKernelRel));
                EXPECT_TRUE(Near(t2_b[row * 5 + component],
                                 expected_native.T2[component],
                                 kKernelAbs, kKernelRel))
                    << "piquad_local_T2.npy component " << component;
            }
        }
        EXPECT_GT(valid_rows, 0u);
        RemoveOutputDirectory(moved_dir);
    }
    RemoveOutputDirectory(original_dir);
}


TEST(DirectionalCovariance,
     WaterHBondCandidateAbsolutePositionsRerunO3Serialized) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    // EvaluateGeometry reports acos(dot) in degrees.  The deliberately
    // collinear 180-degree forcing row is ill-conditioned at |cos|=1:
    // O(epsilon) dot-product reassociation becomes O(sqrt(epsilon)) in acos.
    // The observed transformed discrepancy is 8.5377363e-7 degree; retain a
    // named 2e-6 degree endpoint bound while distances/positions stay on the
    // analytic geometry tolerance above.
    constexpr double kWaterAngleAbsDegree = 2.0e-6;
    auto protein = LoadUbqOrSkip();
    if (!protein) GTEST_SKIP() << "1UBQ protonated fixture unavailable";

    ProteinConformation& original = protein->Conformation();
    ASSERT_TRUE(AttachGeometryAndSpatial(original));
    std::size_t acceptor_atom = SIZE_MAX;
    std::size_t polar_h_atom = SIZE_MAX;
    const nmr::SolventEnvironment solvent = BuildTypedWaterHBondSolvent(
        original, &acceptor_atom, &polar_h_atom);
    ASSERT_NE(acceptor_atom, SIZE_MAX);
    ASSERT_NE(polar_h_atom, SIZE_MAX);
    ASSERT_EQ(solvent.waters.size(), 2u);
    auto result = nmr::WaterHBondGeometryResult::Compute(original, solvent);
    ASSERT_NE(result, nullptr);
    const std::filesystem::path original_dir =
        FreshOutputDirectory("directional_water_hbond_original");
    ASSERT_EQ(result->WriteFeatures(original, original_dir.string()), 3);
    const SerializedNpy candidates_a_npy = ReadSerializedNpy(
        original_dir / "water_hbond_candidates.npy");
    const SerializedNpy counts_a_npy = ReadSerializedNpy(
        original_dir / "water_hbond_counts.npy");
    const SerializedNpy nearest_a_npy = ReadSerializedNpy(
        original_dir / "water_hbond_nearest.npy");
    const auto candidates_a = NpyValues<double>(candidates_a_npy, "<f8");
    const auto counts_a = NpyValues<std::int32_t>(counts_a_npy, "<i4");
    const auto nearest_a = NpyValues<double>(nearest_a_npy, "<f8");
    ASSERT_EQ(candidates_a_npy.shape.size(), 2u);
    ASSERT_EQ(candidates_a_npy.shape[1], 16u);
    const std::size_t candidate_count = candidates_a_npy.shape[0];
    ASSERT_GT(candidate_count, 0u);

    for (const bool improper : {false, true}) {
        const OrthogonalTransform x =
            SeededTransform(0xA27B0D11ULL, improper);
        ProteinConformation& moved = protein->AddConformation(
            nmr::test::directional::Positions(x, original.Positions()),
            improper ? "water H-bond improper rerun" :
                       "water H-bond proper rerun");
        ASSERT_TRUE(AttachGeometryAndSpatial(moved));
        const nmr::SolventEnvironment moved_solvent =
            nmr::test::directional::Solvent(x, solvent);
        auto moved_result =
            nmr::WaterHBondGeometryResult::Compute(moved, moved_solvent);
        ASSERT_NE(moved_result, nullptr);
        const std::filesystem::path moved_dir = FreshOutputDirectory(
            improper ? "directional_water_hbond_improper" :
                       "directional_water_hbond_proper");
        ASSERT_EQ(moved_result->WriteFeatures(moved, moved_dir.string()), 3);
        const SerializedNpy candidates_b_npy = ReadSerializedNpy(
            moved_dir / "water_hbond_candidates.npy");
        const SerializedNpy counts_b_npy = ReadSerializedNpy(
            moved_dir / "water_hbond_counts.npy");
        const SerializedNpy nearest_b_npy = ReadSerializedNpy(
            moved_dir / "water_hbond_nearest.npy");
        const auto candidates_b =
            NpyValues<double>(candidates_b_npy, "<f8");
        const auto counts_b =
            NpyValues<std::int32_t>(counts_b_npy, "<i4");
        const auto nearest_b = NpyValues<double>(nearest_b_npy, "<f8");
        ASSERT_EQ(candidates_b_npy.shape, candidates_a_npy.shape);
        ASSERT_EQ(counts_b_npy.shape, counts_a_npy.shape);
        ASSERT_EQ(nearest_b_npy.shape, nearest_a_npy.shape);
        EXPECT_EQ(counts_b, counts_a)
            << "water_hbond_counts.npy mask/count rows changed";

        std::size_t mode1_rows = 0;
        std::size_t mode2_rows = 0;
        for (std::size_t row = 0; row < candidate_count; ++row) {
            SCOPED_TRACE("water H-bond row=" + std::to_string(row));
            for (int identity = 0; identity < 5; ++identity)
                EXPECT_DOUBLE_EQ(candidates_b[row * 16 + identity],
                                 candidates_a[row * 16 + identity]);
            for (int scalar = 5; scalar <= 8; ++scalar) {
                ExpectScalarOrSameNaN(
                    candidates_b[row * 16 + scalar],
                    candidates_a[row * 16 + scalar],
                    scalar == 8 ? kWaterAngleAbsDegree : kGeometryAbs,
                    scalar == 8 ? 0.0 : kGeometryRel,
                    "water_hbond_candidates.npy scalar geometry");
            }
            EXPECT_DOUBLE_EQ(candidates_b[row * 16 + 15],
                             candidates_a[row * 16 + 15]);

            const Vec3 oxygen_a(candidates_a[row * 16 + 9],
                                candidates_a[row * 16 + 10],
                                candidates_a[row * 16 + 11]);
            const Vec3 oxygen_b(candidates_b[row * 16 + 9],
                                candidates_b[row * 16 + 10],
                                candidates_b[row * 16 + 11]);
            ASSERT_TRUE(oxygen_a.allFinite());
            ASSERT_TRUE(oxygen_b.allFinite());
            ExpectVector(oxygen_b, Position(x, oxygen_a),
                         kGeometryAbs, kGeometryRel,
                         "water_hbond_candidates.npy water O position");

            const int mode =
                static_cast<int>(candidates_a[row * 16 + 3]);
            const Vec3 hydrogen_a(candidates_a[row * 16 + 12],
                                  candidates_a[row * 16 + 13],
                                  candidates_a[row * 16 + 14]);
            const Vec3 hydrogen_b(candidates_b[row * 16 + 12],
                                  candidates_b[row * 16 + 13],
                                  candidates_b[row * 16 + 14]);
            if (mode == 1) {
                ++mode1_rows;
                ASSERT_TRUE(hydrogen_a.allFinite());
                ASSERT_TRUE(hydrogen_b.allFinite());
                ExpectVector(hydrogen_b, Position(x, hydrogen_a),
                             kGeometryAbs, kGeometryRel,
                             "water_hbond_candidates.npy chosen H position");
            } else {
                ASSERT_EQ(mode, 2);
                ++mode2_rows;
                for (int component = 0; component < 3; ++component) {
                    EXPECT_TRUE(std::isnan(hydrogen_a[component]));
                    EXPECT_TRUE(std::isnan(hydrogen_b[component]))
                        << "mode-2 unavailable water-H position must stay NaN";
                }
            }
        }
        EXPECT_GT(mode1_rows, 0u);
        EXPECT_GT(mode2_rows, 0u);
        ASSERT_EQ(nearest_b.size(), nearest_a.size());
        for (std::size_t value = 0; value < nearest_a.size(); ++value) {
            ExpectScalarOrSameNaN(nearest_b[value], nearest_a[value],
                                  value % 8 == 2 ?
                                      kWaterAngleAbsDegree : kGeometryAbs,
                                  value % 8 == 2 ? 0.0 : kGeometryRel,
                                  "water_hbond_nearest.npy");
        }
        RemoveOutputDirectory(moved_dir);
    }
    RemoveOutputDirectory(original_dir);
}


TEST(DirectionalCovariance, NativeT2BasisRoundTripUsesCartesianQtQt) {
    SphericalTensor source;
    source.T2 = {0.37, -1.2, 0.84, 2.1, -0.43};
    for (const bool improper : {false, true}) {
        const OrthogonalTransform x =
            SeededTransform(0x725B4515ULL, improper);
        const SphericalTensor rotated = RotateNativeT2(x, source);
        EXPECT_NEAR(rotated.T0, 0.0, 2.0e-15);
        EXPECT_NEAR(T1Vector(rotated).norm(), 0.0, 2.0e-15);
        ExpectMatrix(T2Matrix(rotated), EvenRank2(x, T2Matrix(source)),
                     3.0e-15, 3.0e-15,
                     "project-native T2 reconstruct/rotate/decompose");
    }
}

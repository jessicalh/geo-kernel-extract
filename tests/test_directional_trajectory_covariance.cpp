#include "DirectionalTestHelpers.h"
#include "TestEnvironment.h"

#include <gtest/gtest.h>

#include "BiotSavartResult.h"
#include "BsShieldingTimeSeriesTrajectoryResult.h"
#include "BsWelfordTrajectoryResult.h"
#include "GeometryResult.h"
#include "HaighMallionResult.h"
#include "HmShieldingTimeSeriesTrajectoryResult.h"
#include "HmWelfordTrajectoryResult.h"
#include "HydrationGeometryResult.h"
#include "HydrationGeometryTimeSeriesTrajectoryResult.h"
#include "HydrationGeometryWelfordTrajectoryResult.h"
#include "LocalBackboneGeometryTrajectoryResult.h"
#include "McConnellResult.h"
#include "McConnellShieldingTimeSeriesTrajectoryResult.h"
#include "McConnellWelfordTrajectoryResult.h"
#include "PdbFileReader.h"
#include "PositionsTimeSeriesTrajectoryResult.h"
#include "ReorientationalDynamicsTrajectoryResult.h"
#include "SasaResult.h"
#include "SasaTimeSeriesTrajectoryResult.h"
#include "SasaWelfordTrajectoryResult.h"
#include "SolventEnvironment.h"
#include "SpatialIndexResult.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"
#include "WaterFieldResult.h"
#include "WaterFieldTimeSeriesTrajectoryResult.h"
#include "WaterFieldWelfordTrajectoryResult.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5File.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <typeindex>
#include <utility>
#include <vector>
#include <unistd.h>

namespace {

using nmr::Mat3;
using nmr::Protein;
using nmr::ProteinConformation;
using nmr::SphericalTensor;
using nmr::Trajectory;
using nmr::TrajectoryProtein;
using nmr::Vec3;
using nmr::test::directional::EvenRank2;
using nmr::test::directional::Axial;
using nmr::test::directional::Near;
using nmr::test::directional::NearMatrix;
using nmr::test::directional::NearVector;
using nmr::test::directional::OrthogonalTransform;
using nmr::test::directional::Polar;
using nmr::test::directional::Position;
using nmr::test::directional::RotateNativeT2;
using nmr::test::directional::SeededTransform;
using nmr::test::directional::T1Vector;

std::unique_ptr<Protein> LoadUbqOrSkip() {
    const std::string path = nmr::test::TestEnvironment::UbqProtonated();
    if (path.empty() || !std::filesystem::exists(path)) return nullptr;
    auto loaded = nmr::BuildFromProtonatedPdb(path);
    if (!loaded.Ok()) return nullptr;
    return std::move(loaded.protein);
}

template <typename T>
std::vector<T> ReadH5Flat(const std::string& path,
                          const std::string& dataset,
                          std::vector<std::size_t>* dimensions = nullptr) {
    HighFive::File file(path, HighFive::File::ReadOnly);
    auto data_set = file.getDataSet(dataset);
    const std::vector<std::size_t> dims =
        data_set.getSpace().getDimensions();
    if (dimensions) *dimensions = dims;
    const std::size_t count = std::accumulate(
        dims.begin(), dims.end(), std::size_t{1},
        std::multiplies<std::size_t>());
    std::vector<T> values(count);
    if (!values.empty()) data_set.read(values.data());
    return values;
}

std::string ReadH5StringAttribute(const std::string& path,
                                  const std::string& dataset,
                                  const std::string& attribute) {
    HighFive::File file(path, HighFive::File::ReadOnly);
    std::string value;
    file.getDataSet(dataset).getAttribute(attribute).read(value);
    return value;
}

std::string ReadH5GroupStringAttribute(const std::string& path,
                                       const std::string& group,
                                       const std::string& attribute) {
    HighFive::File file(path, HighFive::File::ReadOnly);
    std::string value;
    file.getGroup(group).getAttribute(attribute).read(value);
    return value;
}

bool ReadH5BoolAttribute(const std::string& path,
                         const std::string& dataset,
                         const std::string& attribute) {
    HighFive::File file(path, HighFive::File::ReadOnly);
    bool value = false;
    file.getDataSet(dataset).getAttribute(attribute).read(value);
    return value;
}

bool ReadH5GroupBoolAttribute(const std::string& path,
                              const std::string& group,
                              const std::string& attribute) {
    HighFive::File file(path, HighFive::File::ReadOnly);
    bool value = false;
    file.getGroup(group).getAttribute(attribute).read(value);
    return value;
}

std::string FreshH5Path(const std::string& stem) {
    const std::string path = nmr::test::TestEnvironment::TempPath(stem);
    (void)std::remove(path.c_str());
    return path;
}

void RemoveH5(const std::string& path) {
    EXPECT_EQ(std::remove(path.c_str()), 0) << path;
}

std::string ShellQuote(const std::string& value) {
    std::string quoted = "'";
    for (const char c : value) {
        if (c == '\'') quoted += "'\\''";
        else quoted += c;
    }
    quoted += "'";
    return quoted;
}

std::vector<double> ReadFloat64Npy(
        const std::filesystem::path& path,
        std::vector<std::size_t>* dimensions = nullptr) {
    std::ifstream input(path, std::ios::binary);
    EXPECT_TRUE(input.is_open()) << path;
    if (!input.is_open()) return {};
    char magic[6] = {};
    input.read(magic, sizeof(magic));
    EXPECT_EQ(std::string(magic, sizeof(magic)),
              std::string("\x93NUMPY", sizeof(magic)));
    unsigned char version[2] = {};
    input.read(reinterpret_cast<char*>(version), sizeof(version));
    EXPECT_TRUE((version[0] == 1u || version[0] == 2u) && version[1] == 0u);
    std::uint32_t header_length = 0;
    if (version[0] == 1u) {
        unsigned char bytes[2] = {};
        input.read(reinterpret_cast<char*>(bytes), sizeof(bytes));
        header_length = static_cast<std::uint32_t>(bytes[0]) |
                        (static_cast<std::uint32_t>(bytes[1]) << 8u);
    } else {
        unsigned char bytes[4] = {};
        input.read(reinterpret_cast<char*>(bytes), sizeof(bytes));
        header_length = static_cast<std::uint32_t>(bytes[0]) |
                        (static_cast<std::uint32_t>(bytes[1]) << 8u) |
                        (static_cast<std::uint32_t>(bytes[2]) << 16u) |
                        (static_cast<std::uint32_t>(bytes[3]) << 24u);
    }
    std::string header(header_length, '\0');
    input.read(header.data(), static_cast<std::streamsize>(header.size()));
    EXPECT_NE(header.find("'descr': '<f8'"), std::string::npos) << header;
    const std::size_t shape_begin = header.find('(');
    const std::size_t shape_end = header.find(')', shape_begin);
    EXPECT_NE(shape_begin, std::string::npos) << header;
    EXPECT_NE(shape_end, std::string::npos) << header;
    std::vector<std::size_t> parsed_dimensions;
    if (shape_begin != std::string::npos && shape_end != std::string::npos) {
        std::stringstream shape_stream(
            header.substr(shape_begin + 1, shape_end - shape_begin - 1));
        std::string token;
        while (std::getline(shape_stream, token, ',')) {
            const std::size_t first = token.find_first_not_of(" \t");
            if (first != std::string::npos)
                parsed_dimensions.push_back(
                    static_cast<std::size_t>(std::stoull(token.substr(first))));
        }
    }
    if (dimensions) *dimensions = parsed_dimensions;
    const std::size_t count = std::accumulate(
        parsed_dimensions.begin(), parsed_dimensions.end(), std::size_t{1},
        std::multiplies<std::size_t>());
    std::vector<double> values(count);
    if (!values.empty()) {
        input.read(reinterpret_cast<char*>(values.data()),
                   static_cast<std::streamsize>(values.size() * sizeof(double)));
    }
    EXPECT_TRUE(input.good() || input.eof()) << path;
    return values;
}

void ExpectNumpyRejectsPickleRequirement(const std::filesystem::path& path) {
    const std::string script =
        "import numpy as np,sys; a=np.load(sys.argv[1],allow_pickle=False); "
        "assert a.dtype==np.float64 and a.ndim==2";
    const std::string command = ShellQuote(NMR_TEST_PYTHON_EXECUTABLE) +
        " -c " + ShellQuote(script) + " " + ShellQuote(path.string());
    EXPECT_EQ(std::system(command.c_str()), 0) << command;
}

std::vector<std::vector<Vec3>> BuildInternalMotionFrames(
        const std::vector<Vec3>& reference) {
    std::vector<std::vector<Vec3>> frames;
    constexpr std::size_t kFrameCount = 6;
    frames.reserve(kFrameCount);
    for (std::size_t frame = 0; frame < kFrameCount; ++frame) {
        std::vector<Vec3> positions = reference;
        if (frame != 0) {
            for (std::size_t atom = 0; atom < positions.size(); ++atom) {
                const double a = 0.071 * static_cast<double>(atom) +
                                 0.47 * static_cast<double>(frame);
                const double b = 0.053 * static_cast<double>(atom) -
                                 0.31 * static_cast<double>(frame);
                positions[atom] += Vec3(
                    0.018 * std::sin(a),
                    0.014 * std::cos(b),
                    0.011 * std::sin(a + b));
            }
        }
        frames.push_back(std::move(positions));
    }
    return frames;
}

void ExpectSameIncludingNaN(const std::vector<double>& actual,
                            const std::vector<double>& expected,
                            double abs_tolerance,
                            double rel_tolerance,
                            const char* quantity) {
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t index = 0; index < expected.size(); ++index) {
        SCOPED_TRACE(std::string(quantity) + " index=" +
                     std::to_string(index));
        if (std::isnan(expected[index])) {
            EXPECT_TRUE(std::isnan(actual[index]));
        } else {
            ASSERT_TRUE(std::isfinite(expected[index]));
            ASSERT_TRUE(std::isfinite(actual[index]));
            EXPECT_TRUE(Near(actual[index], expected[index],
                             abs_tolerance, rel_tolerance))
                << "actual=" << actual[index]
                << " expected=" << expected[index];
        }
    }
}

template <typename T>
void ExpectExactH5(const std::string& source_path,
                   const std::string& moved_path,
                   const std::string& dataset,
                   const std::vector<std::size_t>& expected_dimensions) {
    std::vector<std::size_t> source_dimensions;
    std::vector<std::size_t> moved_dimensions;
    const auto source = ReadH5Flat<T>(
        source_path, dataset, &source_dimensions);
    const auto moved = ReadH5Flat<T>(
        moved_path, dataset, &moved_dimensions);
    EXPECT_EQ(source_dimensions, expected_dimensions) << dataset;
    EXPECT_EQ(moved_dimensions, expected_dimensions) << dataset;
    EXPECT_EQ(moved, source) << dataset;
}

void ExpectPolarVectorTimeSeries(
        const std::string& source_path,
        const std::string& moved_path,
        const std::string& dataset,
        const OrthogonalTransform& transform,
        std::size_t row_count,
        std::size_t frame_count,
        double abs_tolerance,
        double rel_tolerance) {
    std::vector<std::size_t> source_dimensions;
    std::vector<std::size_t> moved_dimensions;
    const auto source = ReadH5Flat<double>(
        source_path, dataset, &source_dimensions);
    const auto moved = ReadH5Flat<double>(
        moved_path, dataset, &moved_dimensions);
    const std::vector<std::size_t> expected_dimensions = {
        row_count, frame_count, 3u};
    ASSERT_EQ(source_dimensions, expected_dimensions) << dataset;
    ASSERT_EQ(moved_dimensions, expected_dimensions) << dataset;
    ASSERT_EQ(source.size(), moved.size());

    double signal = 0.0;
    for (std::size_t row = 0; row < row_count * frame_count; ++row) {
        SCOPED_TRACE(dataset + " row=" + std::to_string(row));
        const Vec3 source_vector(source[row * 3], source[row * 3 + 1],
                                 source[row * 3 + 2]);
        const Vec3 moved_vector(moved[row * 3], moved[row * 3 + 1],
                                moved[row * 3 + 2]);
        if (!source_vector.allFinite()) {
            EXPECT_TRUE(std::isnan(source_vector.x()));
            EXPECT_TRUE(std::isnan(source_vector.y()));
            EXPECT_TRUE(std::isnan(source_vector.z()));
            EXPECT_TRUE(std::isnan(moved_vector.x()));
            EXPECT_TRUE(std::isnan(moved_vector.y()));
            EXPECT_TRUE(std::isnan(moved_vector.z()));
            continue;
        }
        ASSERT_TRUE(moved_vector.allFinite());
        signal += source_vector.norm();
        EXPECT_TRUE(NearVector(moved_vector, Polar(transform, source_vector),
                               abs_tolerance, rel_tolerance))
            << "actual=" << moved_vector.transpose()
            << " expected="
            << Polar(transform, source_vector).transpose();
    }
    EXPECT_GT(signal, 1.0e-18) << dataset << " covariance is vacuous";
}

void ExpectFullTensorTimeSeries(
        const std::string& source_path,
        const std::string& moved_path,
        const std::string& dataset,
        const OrthogonalTransform& transform,
        std::size_t atom_count,
        std::size_t frame_count,
        double abs_tolerance,
        double rel_tolerance) {
    std::vector<std::size_t> source_dimensions;
    std::vector<std::size_t> moved_dimensions;
    const auto source_values = ReadH5Flat<double>(
        source_path, dataset, &source_dimensions);
    const auto moved_values = ReadH5Flat<double>(
        moved_path, dataset, &moved_dimensions);
    const std::vector<std::size_t> expected_dimensions = {
        atom_count, frame_count, 9u};
    ASSERT_EQ(source_dimensions, expected_dimensions) << dataset;
    ASSERT_EQ(moved_dimensions, expected_dimensions) << dataset;
    ASSERT_EQ(source_values.size(), moved_values.size());

    double t1_signal = 0.0;
    double t2_signal = 0.0;
    for (std::size_t row = 0; row < atom_count * frame_count; ++row) {
        SCOPED_TRACE(dataset + " row=" + std::to_string(row));
        const std::size_t base = row * 9;
        bool source_has_nan = false;
        for (std::size_t component = 0; component < 9; ++component)
            source_has_nan = source_has_nan ||
                             std::isnan(source_values[base + component]);
        if (source_has_nan) {
            for (std::size_t component = 0; component < 9; ++component) {
                EXPECT_TRUE(std::isnan(source_values[base + component]));
                EXPECT_TRUE(std::isnan(moved_values[base + component]));
            }
            continue;
        }

        SphericalTensor source;
        SphericalTensor actual;
        source.T0 = source_values[base];
        actual.T0 = moved_values[base];
        for (int component = 0; component < 3; ++component) {
            source.T1[component] = source_values[base + 1 + component];
            actual.T1[component] = moved_values[base + 1 + component];
        }
        for (int component = 0; component < 5; ++component) {
            source.T2[component] = source_values[base + 4 + component];
            actual.T2[component] = moved_values[base + 4 + component];
        }
        ASSERT_TRUE(std::isfinite(actual.T0));
        EXPECT_TRUE(Near(actual.T0, source.T0,
                         abs_tolerance, rel_tolerance));
        EXPECT_TRUE(NearVector(T1Vector(actual),
                               Axial(transform, T1Vector(source)),
                               abs_tolerance, rel_tolerance));
        const SphericalTensor expected_t2 =
            RotateNativeT2(transform, source);
        for (int component = 0; component < 5; ++component) {
            EXPECT_TRUE(Near(actual.T2[component],
                             expected_t2.T2[component],
                             abs_tolerance, rel_tolerance))
                << "native T2 component=" << component;
        }
        EXPECT_TRUE(NearMatrix(actual.Reconstruct(),
                               EvenRank2(transform, source.Reconstruct()),
                               abs_tolerance, rel_tolerance));
        t1_signal += T1Vector(source).norm();
        t2_signal += nmr::test::directional::T2Matrix(source).norm();
    }
    EXPECT_GT(t1_signal, 1.0e-18) << dataset << " T1 is vacuous";
    EXPECT_GT(t2_signal, 1.0e-18) << dataset << " T2 is vacuous";
}

void ExpectNativeT2TimeSeries(
        const std::string& source_path,
        const std::string& moved_path,
        const std::string& dataset,
        const OrthogonalTransform& transform,
        std::size_t atom_count,
        std::size_t frame_count,
        double abs_tolerance,
        double rel_tolerance) {
    std::vector<std::size_t> source_dimensions;
    std::vector<std::size_t> moved_dimensions;
    const auto source_values = ReadH5Flat<double>(
        source_path, dataset, &source_dimensions);
    const auto moved_values = ReadH5Flat<double>(
        moved_path, dataset, &moved_dimensions);
    const std::vector<std::size_t> expected_dimensions = {
        atom_count, frame_count, 5u};
    ASSERT_EQ(source_dimensions, expected_dimensions) << dataset;
    ASSERT_EQ(moved_dimensions, expected_dimensions) << dataset;
    ASSERT_EQ(source_values.size(), moved_values.size());

    double signal = 0.0;
    for (std::size_t row = 0; row < atom_count * frame_count; ++row) {
        SCOPED_TRACE(dataset + " row=" + std::to_string(row));
        const std::size_t base = row * 5;
        bool source_has_nan = false;
        for (std::size_t component = 0; component < 5; ++component)
            source_has_nan = source_has_nan ||
                             std::isnan(source_values[base + component]);
        if (source_has_nan) {
            for (std::size_t component = 0; component < 5; ++component) {
                EXPECT_TRUE(std::isnan(source_values[base + component]));
                EXPECT_TRUE(std::isnan(moved_values[base + component]));
            }
            continue;
        }
        SphericalTensor source;
        SphericalTensor actual;
        for (int component = 0; component < 5; ++component) {
            source.T2[component] = source_values[base + component];
            actual.T2[component] = moved_values[base + component];
        }
        const SphericalTensor expected = RotateNativeT2(transform, source);
        for (int component = 0; component < 5; ++component) {
            EXPECT_TRUE(Near(actual.T2[component],
                             expected.T2[component],
                             abs_tolerance, rel_tolerance));
        }
        EXPECT_TRUE(NearMatrix(
            nmr::test::directional::T2Matrix(actual),
            EvenRank2(transform,
                      nmr::test::directional::T2Matrix(source)),
            abs_tolerance, rel_tolerance));
        signal += nmr::test::directional::T2Matrix(source).norm();
    }
    EXPECT_GT(signal, 1.0e-18) << dataset << " T2 is vacuous";
}

bool AttachKernelOwners(ProteinConformation& conf) {
    auto geometry = nmr::GeometryResult::Compute(conf);
    if (!geometry || !conf.AttachResult(std::move(geometry))) return false;
    auto spatial = nmr::SpatialIndexResult::Compute(conf);
    if (!spatial || !conf.AttachResult(std::move(spatial))) return false;
    auto bs = nmr::BiotSavartResult::Compute(conf);
    if (!bs || !conf.AttachResult(std::move(bs))) return false;
    auto hm = nmr::HaighMallionResult::Compute(conf);
    if (!hm || !conf.AttachResult(std::move(hm))) return false;
    auto mc = nmr::McConnellResult::Compute(conf);
    return mc && conf.AttachResult(std::move(mc));
}

void ExpectNonClosedComponentStatisticsMetadata(
        const std::string& path,
        const std::string& group_path,
        const std::vector<std::string>& directional_prefixes) {
    HighFive::File file(path, HighFive::File::ReadOnly);
    auto group = file.getGroup(group_path);
    std::string scope;
    std::string component_law;
    group.getAttribute("irrep_metadata_scope").read(scope);
    group.getAttribute("componentwise_statistic_transformation")
        .read(component_law);
    EXPECT_NE(scope.find("only assembled"), std::string::npos);
    EXPECT_NE(scope.find("means"), std::string::npos);
    EXPECT_NE(component_law.find("no closed irrep transformation law"),
              std::string::npos);

    // Presence is asserted, covariance is deliberately not: these are
    // nonlinear componentwise statistics, so applying an irrep matrix to
    // m2/std/min/max/frame extrema would be a false scientific claim.
    for (const std::string& prefix : directional_prefixes) {
        for (const char* suffix : {
                 "_m2", "_std", "_min", "_max",
                 "_min_frame", "_max_frame"}) {
            EXPECT_TRUE(group.exist(prefix + suffix))
                << group_path << "/" << prefix << suffix;
        }
    }
}

void ExpectKernelWelfordDirectionalMetadata(
        const std::string& path, const std::string& group) {
    const auto attr = [&](const char* name) {
        return ReadH5GroupStringAttribute(path, group, name);
    };
    EXPECT_EQ(attr("coordinate_frame"), "conformation_cartesian_xyz");
    EXPECT_EQ(attr("source_tensor_basis"),
              "project_native_full9_spherical_tensor_v1");
    EXPECT_EQ(attr("source_tensor_component_order"),
              "T0,T1_x,T1_y,T1_z,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2");
    EXPECT_EQ(attr("source_tensor_parity"), "even");
    EXPECT_EQ(attr("source_tensor_transformation"),
              "even_rank2: T'=R T R^T");
    EXPECT_EQ(attr("source_tensor_structural_zero_components"), "none");
    EXPECT_EQ(attr("t1_basis"),
              "project_native_cartesian_levi_civita_dual_xyz_v1");
    EXPECT_EQ(attr("t1_component_order"), "T1_x,T1_y,T1_z");
    EXPECT_EQ(attr("t1_frame"), "conformation_cartesian_xyz");
    EXPECT_EQ(attr("t1_parity"), "even");
    EXPECT_EQ(attr("t1_transformation"),
              "axial_vector: a'=det(R) R a");
    EXPECT_EQ(attr("t1_semantics"),
              "Cartesian Levi-Civita dual x,y,z (not real-Y1m): "
              "a=((T_yz-T_zy)/2,(T_zx-T_xz)/2,(T_xy-T_yx)/2); "
              "axial a'=det(R) R a; generically nonzero");
    EXPECT_FALSE(ReadH5GroupBoolAttribute(path, group,
                                          "t1_structural_zero"));
    EXPECT_EQ(attr("t2_basis"),
              "project_native_t2_isometric_real_tesseral_v1");
    EXPECT_EQ(attr("t2_component_order"),
              "T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2");
    EXPECT_EQ(attr("t2_frame"), "conformation_cartesian_xyz");
    EXPECT_EQ(attr("t2_parity"), "even");
    EXPECT_EQ(attr("t2_transformation"),
              "even_rank2: T'=R T R^T");
    EXPECT_FALSE(ReadH5GroupBoolAttribute(path, group,
                                          "t2_structural_zero"));
    EXPECT_EQ(attr("normalization"), "isometric_real_sph");
    EXPECT_EQ(attr("directional_mean_transformation"),
              "t0_mean is invariant; assembled t1_mean is axial: "
              "a'=det(R) R a; assembled t2_mean is even rank-2: "
              "T'=R T R^T");
    EXPECT_NE(attr("normalization_scope").find("Levi-Civita convention"),
              std::string::npos);
    EXPECT_NE(attr("e3nn_export").find("project_t2_to_e3nn"),
              std::string::npos);
    EXPECT_NE(attr("directional_metadata_scope").find("only t1_mean and t2_mean"),
              std::string::npos);
}

void ExpectWaterWelfordDirectionalMetadata(const std::string& path) {
    const std::string group = "/trajectory/water_field_welford";
    const auto attr = [&](const char* name) {
        return ReadH5GroupStringAttribute(path, group, name);
    };
    EXPECT_EQ(attr("efield_coordinate_frame"),
              "conformation_cartesian_xyz");
    EXPECT_EQ(attr("efield_normalization"), "cartesian");
    EXPECT_EQ(attr("efield_parity"), "1o");
    EXPECT_EQ(attr("efield_mean_transformation"),
              "assembled efield_{x,y,z}_mean and "
              "efield_first_{x,y,z}_mean are polar-vector means: v'=R v");
    EXPECT_EQ(attr("efg_t2_basis"),
              "project_native_t2_isometric_real_tesseral_v1");
    EXPECT_EQ(attr("efg_t2_component_order"),
              "T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2");
    EXPECT_EQ(attr("efg_t2_frame"), "conformation_cartesian_xyz");
    EXPECT_EQ(attr("efg_t2_normalization"), "isometric_real_sph");
    EXPECT_EQ(attr("efg_t2_parity"), "even");
    EXPECT_EQ(attr("efg_t2_mean_transformation"),
              "assembled efg_t2_mean and efg_first_t2_mean are native "
              "isometric real-tesseral T2 means of even-rank-2 Cartesian "
              "tensors: T'=R T R^T");
    EXPECT_NE(attr("efg_e3nn_export").find("project_t2_to_e3nn"),
              std::string::npos);
    EXPECT_NE(attr("directional_metadata_scope").find(
                  "56 directional-statistic dataset paths only"),
              std::string::npos);
    EXPECT_TRUE(ReadH5GroupBoolAttribute(path, group,
                                         "efg_t0_structural_zero"));
    EXPECT_TRUE(ReadH5GroupBoolAttribute(path, group,
                                         "efg_t1_structural_zero"));
}

void ExpectReorientationTensorDirectionalMetadata(const std::string& path) {
    const std::string dataset =
        "/trajectory/reorientational_dynamics/bond_orientation_tensor";
    const auto attr = [&](const char* name) {
        return ReadH5StringAttribute(path, dataset, name);
    };
    EXPECT_EQ(attr("tensor_basis"), "cartesian_matrix_row_major");
    EXPECT_EQ(attr("tensor_component_order"),
              "XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ");
    EXPECT_EQ(attr("tensor_frame"),
              "reference_conformation_cartesian_xyz");
    EXPECT_EQ(attr("coordinate_frame"),
              "reference_conformation_cartesian_xyz");
    EXPECT_EQ(attr("tensor_parity"), "even");
    EXPECT_EQ(attr("parity"), "even");
    EXPECT_EQ(attr("irrep_layout"), "0e+2e");
    EXPECT_NE(attr("tensor_transformation").find("T'=R T R^T"),
              std::string::npos);
    EXPECT_EQ(attr("tensor_symmetry"), "symmetric");
    EXPECT_EQ(attr("tensor_trace"),
              "trace=1 for every finite physical row because each sample "
              "is a unit-vector outer product");
    EXPECT_FALSE(ReadH5BoolAttribute(path, dataset,
                                     "tensor_t0_structural_zero"));
    EXPECT_EQ(attr("tensor_t0_semantics"),
              "project SphericalTensor T0=trace(T)/3=1/3; fixed nonzero signal");
    EXPECT_TRUE(ReadH5BoolAttribute(path, dataset,
                                    "tensor_t1_structural_zero"));
    EXPECT_EQ(attr("tensor_t1_semantics"),
              "project Cartesian Levi-Civita dual "
              "(T_yz-T_zy,T_zx-T_xz,T_xy-T_yx)/2; zero by symmetry");
    EXPECT_FALSE(ReadH5BoolAttribute(path, dataset,
                                     "tensor_t2_structural_zero"));
    EXPECT_EQ(attr("tensor_structural_zero_components"), "T1 only");
    EXPECT_NE(attr("e3nn_export").find("project-to-e3nn basis conversion"),
              std::string::npos);
}

void ExpectTensorWelfordMeans(
        const std::string& source_path,
        const std::string& moved_path,
        const std::string& group_path,
        const OrthogonalTransform& x,
        double abs_tolerance,
        double rel_tolerance) {
    std::vector<std::size_t> t0_dims;
    std::vector<std::size_t> t1_dims;
    std::vector<std::size_t> t2_dims;
    const auto t0_a = ReadH5Flat<double>(
        source_path, group_path + "/t0_mean", &t0_dims);
    const auto t1_a = ReadH5Flat<double>(
        source_path, group_path + "/t1_mean", &t1_dims);
    const auto t2_a = ReadH5Flat<double>(
        source_path, group_path + "/t2_mean", &t2_dims);
    const auto t0_b = ReadH5Flat<double>(
        moved_path, group_path + "/t0_mean");
    const auto t1_b = ReadH5Flat<double>(
        moved_path, group_path + "/t1_mean");
    const auto t2_b = ReadH5Flat<double>(
        moved_path, group_path + "/t2_mean");
    const auto n_a = ReadH5Flat<std::size_t>(
        source_path, group_path + "/n_frames_per_atom");
    const auto n_b = ReadH5Flat<std::size_t>(
        moved_path, group_path + "/n_frames_per_atom");
    ASSERT_EQ(t0_dims.size(), 1u);
    ASSERT_EQ(t1_dims,
              (std::vector<std::size_t>{t0_dims[0], 3}));
    ASSERT_EQ(t2_dims,
              (std::vector<std::size_t>{t0_dims[0], 5}));
    ASSERT_EQ(t0_b.size(), t0_a.size());
    ASSERT_EQ(t1_b.size(), t1_a.size());
    ASSERT_EQ(t2_b.size(), t2_a.size());
    ASSERT_EQ(n_a.size(), t0_a.size());
    EXPECT_EQ(n_b, n_a);
    const std::size_t max_samples =
        *std::max_element(n_a.begin(), n_a.end());
    ASSERT_GT(max_samples, 0u);

    double t1_signal = 0.0;
    double t2_signal = 0.0;
    std::size_t conditionally_omitted_rows = 0;
    for (std::size_t atom = 0; atom < t0_a.size(); ++atom) {
        SCOPED_TRACE(group_path + " atom=" + std::to_string(atom));
        SphericalTensor source;
        SphericalTensor actual;
        source.T0 = t0_a[atom];
        actual.T0 = t0_b[atom];
        for (int component = 0; component < 3; ++component) {
            source.T1[component] = t1_a[atom * 3 + component];
            actual.T1[component] = t1_b[atom * 3 + component];
        }
        for (int component = 0; component < 5; ++component) {
            source.T2[component] = t2_a[atom * 5 + component];
            actual.T2[component] = t2_b[atom * 5 + component];
        }

        const bool source_has_nan =
            std::isnan(source.T0) ||
            std::any_of(source.T1.begin(), source.T1.end(),
                        [](double value) { return std::isnan(value); }) ||
            std::any_of(source.T2.begin(), source.T2.end(),
                        [](double value) { return std::isnan(value); });
        const bool actual_has_nan =
            std::isnan(actual.T0) ||
            std::any_of(actual.T1.begin(), actual.T1.end(),
                        [](double value) { return std::isnan(value); }) ||
            std::any_of(actual.T2.begin(), actual.T2.end(),
                        [](double value) { return std::isnan(value); });
        if (source_has_nan || actual_has_nan) {
            EXPECT_EQ(n_a[atom], 0u);
            EXPECT_TRUE(source_has_nan);
            EXPECT_TRUE(actual_has_nan);
            EXPECT_TRUE(std::isnan(source.T0));
            EXPECT_TRUE(std::isnan(actual.T0));
            for (int component = 0; component < 3; ++component) {
                EXPECT_TRUE(std::isnan(source.T1[component]));
                EXPECT_TRUE(std::isnan(actual.T1[component]));
            }
            for (int component = 0; component < 5; ++component) {
                EXPECT_TRUE(std::isnan(source.T2[component]));
                EXPECT_TRUE(std::isnan(actual.T2[component]));
            }
            ++conditionally_omitted_rows;
            continue;
        }

        EXPECT_GT(n_a[atom], 0u);
        if (n_a[atom] < max_samples) ++conditionally_omitted_rows;

        t1_signal += T1Vector(source).norm();
        t2_signal += nmr::test::directional::T2Matrix(source).norm();

        EXPECT_TRUE(Near(actual.T0, source.T0,
                         abs_tolerance, rel_tolerance));
        EXPECT_TRUE(NearVector(T1Vector(actual),
                               Axial(x, T1Vector(source)),
                               abs_tolerance, rel_tolerance));
        const SphericalTensor expected_t2 = RotateNativeT2(x, source);
        for (int component = 0; component < 5; ++component) {
            EXPECT_TRUE(Near(actual.T2[component],
                             expected_t2.T2[component],
                             abs_tolerance, rel_tolerance))
                << "native T2 component=" << component;
        }
        EXPECT_TRUE(NearMatrix(actual.Reconstruct(),
                               EvenRank2(x, source.Reconstruct()),
                               abs_tolerance, rel_tolerance));
    }
    EXPECT_GT(t1_signal, 1.0e-18)
        << group_path << " T1 mean covariance must be non-vacuous";
    EXPECT_GT(t2_signal, 1.0e-18)
        << group_path << " T2 mean covariance must be non-vacuous";
    if (group_path == "/trajectory/mc_welford") {
        EXPECT_GT(conditionally_omitted_rows, 0u)
            << "the terminal PeptideCO fixture must exercise conditional moments";
    } else {
        EXPECT_EQ(conditionally_omitted_rows, 0u) << group_path;
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
    protein->AddConformation(positions, "directional Welford one atom");
    return protein;
}

std::unique_ptr<Protein> BuildSasaTrajectoryProtein() {
    auto protein = std::make_unique<Protein>();
    nmr::Residue residue;
    residue.type = nmr::AminoAcid::CYS;
    residue.sequence_number = 1;
    residue.chain_id = "A";
    const std::size_t residue_index = protein->AddResidue(std::move(residue));
    const std::vector<nmr::Element> elements = {
        nmr::Element::C, nmr::Element::O, nmr::Element::N, nmr::Element::S};
    const std::vector<Vec3> positions = {
        Vec3(0.0, 0.0, 0.0), Vec3(2.35, 0.25, -0.10),
        Vec3(-1.85, 1.25, 0.30), Vec3(0.40, -2.30, 0.90)};
    for (const nmr::Element element : elements) {
        auto atom = nmr::Atom::Create(element);
        atom->residue_index = residue_index;
        const std::size_t atom_index = protein->AddAtom(std::move(atom));
        protein->MutableResidueAt(residue_index).atom_indices.push_back(
            atom_index);
    }
    protein->FinalizeConstruction(positions);
    protein->AddConformation(positions, "directional SASA trajectory source");
    return protein;
}

std::vector<std::vector<Vec3>> BuildSasaMotionFrames(
        const std::vector<Vec3>& reference) {
    std::vector<std::vector<Vec3>> frames;
    constexpr std::size_t kFrameCount = 4;
    frames.reserve(kFrameCount);
    for (std::size_t frame = 0; frame < kFrameCount; ++frame) {
        std::vector<Vec3> positions = reference;
        const double f = static_cast<double>(frame);
        positions[1] += Vec3(0.08 * f, -0.035 * f, 0.025 * f);
        positions[2] += Vec3(-0.025 * f, 0.07 * f, -0.03 * f);
        positions[3] += Vec3(0.04 * f, 0.02 * f, 0.065 * f);
        frames.push_back(std::move(positions));
    }
    return frames;
}

bool AttachSasaSource(ProteinConformation& conf) {
    auto geometry = nmr::GeometryResult::Compute(conf);
    if (!geometry || !conf.AttachResult(std::move(geometry))) return false;
    auto spatial = nmr::SpatialIndexResult::Compute(conf);
    if (!spatial || !conf.AttachResult(std::move(spatial))) return false;
    auto sasa = nmr::SasaResult::Compute(conf);
    return sasa && conf.AttachResult(std::move(sasa));
}

struct ExpectedScalarMoments {
    double mean = 0.0;
    double m2 = 0.0;
    double std = 0.0;
    double min = 0.0;
    double max = 0.0;
    std::size_t min_frame = 0;
    std::size_t max_frame = 0;
};

ExpectedScalarMoments ComputeExpectedScalarMoments(
        const std::vector<double>& values,
        const std::vector<std::size_t>& frame_indices) {
    EXPECT_FALSE(values.empty());
    EXPECT_EQ(values.size(), frame_indices.size());
    ExpectedScalarMoments out;
    if (values.empty()) return out;
    out.mean = std::accumulate(values.begin(), values.end(), 0.0) /
               static_cast<double>(values.size());
    out.min = values.front();
    out.max = values.front();
    out.min_frame = frame_indices.front();
    out.max_frame = frame_indices.front();
    for (std::size_t i = 0; i < values.size(); ++i) {
        const double centered = values[i] - out.mean;
        out.m2 += centered * centered;
        if (values[i] < out.min) {
            out.min = values[i];
            out.min_frame = frame_indices[i];
        }
        if (values[i] > out.max) {
            out.max = values[i];
            out.max_frame = frame_indices[i];
        }
    }
    out.std = values.size() > 1
        ? std::sqrt(out.m2 / static_cast<double>(values.size() - 1))
        : 0.0;
    return out;
}

void CheckSasaWelfordSerialization(
        const std::string& h5_path,
        const std::filesystem::path& npy_directory,
        std::size_t atom_count,
        std::size_t frame_count) {
    const std::string ts_group = "/trajectory/sasa_time_series/";
    const std::string w_group = "/trajectory/sasa_welford/";
    std::vector<std::size_t> raw_dimensions;
    const auto raw = ReadH5Flat<double>(
        h5_path, ts_group + "sasa", &raw_dimensions);
    ASSERT_EQ(raw_dimensions,
              (std::vector<std::size_t>{atom_count, frame_count}));
    const auto frame_indices = ReadH5Flat<std::size_t>(
        h5_path, ts_group + "frame_indices");
    const auto frame_times = ReadH5Flat<double>(
        h5_path, ts_group + "frame_times");
    ASSERT_EQ(frame_indices.size(), frame_count);
    ASSERT_EQ(frame_times.size(), frame_count);
    const auto n_frames_per_atom = ReadH5Flat<std::size_t>(
        h5_path, w_group + "n_frames_per_atom");
    const auto delta_n_per_atom = ReadH5Flat<std::size_t>(
        h5_path, w_group + "delta_n_per_atom");
    const auto dxdt_n_per_atom = ReadH5Flat<std::size_t>(
        h5_path, w_group + "dxdt_n_per_atom");
    ASSERT_EQ(n_frames_per_atom,
              std::vector<std::size_t>(atom_count, frame_count));
    ASSERT_EQ(delta_n_per_atom,
              std::vector<std::size_t>(atom_count, frame_count - 1u));
    ASSERT_EQ(dxdt_n_per_atom,
              std::vector<std::size_t>(atom_count, frame_count - 1u));

    auto check_channel = [&](const std::string& prefix,
                             const std::vector<double>& samples,
                             std::size_t sample_count,
                             const std::vector<std::size_t>& sample_frames) {
        ASSERT_EQ(samples.size(), atom_count * sample_count);
        ASSERT_EQ(sample_frames.size(), sample_count);
        const auto mean = ReadH5Flat<double>(h5_path, w_group + prefix + "_mean");
        const auto m2 = ReadH5Flat<double>(h5_path, w_group + prefix + "_m2");
        const auto stddev = ReadH5Flat<double>(h5_path, w_group + prefix + "_std");
        const auto min = ReadH5Flat<double>(h5_path, w_group + prefix + "_min");
        const auto max = ReadH5Flat<double>(h5_path, w_group + prefix + "_max");
        const auto min_frame = ReadH5Flat<std::size_t>(
            h5_path, w_group + prefix + "_min_frame");
        const auto max_frame = ReadH5Flat<std::size_t>(
            h5_path, w_group + prefix + "_max_frame");
        ASSERT_EQ(mean.size(), atom_count);
        for (std::size_t atom = 0; atom < atom_count; ++atom) {
            std::vector<double> atom_samples(sample_count);
            for (std::size_t sample = 0; sample < sample_count; ++sample)
                atom_samples[sample] = samples[atom * sample_count + sample];
            const ExpectedScalarMoments expected =
                ComputeExpectedScalarMoments(atom_samples, sample_frames);
            SCOPED_TRACE(prefix + " atom=" + std::to_string(atom));
            EXPECT_TRUE(Near(mean[atom], expected.mean, 1.0e-10, 2.0e-12));
            EXPECT_TRUE(Near(m2[atom], expected.m2, 1.0e-10, 2.0e-12));
            EXPECT_TRUE(Near(stddev[atom], expected.std, 1.0e-10, 2.0e-12));
            EXPECT_TRUE(Near(min[atom], expected.min, 1.0e-10, 2.0e-12));
            EXPECT_TRUE(Near(max[atom], expected.max, 1.0e-10, 2.0e-12));
            EXPECT_EQ(min_frame[atom], expected.min_frame);
            EXPECT_EQ(max_frame[atom], expected.max_frame);
        }
    };

    check_channel("sasa", raw, frame_count, frame_indices);
    const std::size_t delta_count = frame_count - 1;
    std::vector<double> delta(atom_count * delta_count);
    std::vector<double> abs_delta(atom_count * delta_count);
    std::vector<double> delta_squared(atom_count * delta_count);
    std::vector<double> dxdt(atom_count * delta_count);
    for (std::size_t atom = 0; atom < atom_count; ++atom) {
        for (std::size_t sample = 0; sample < delta_count; ++sample) {
            const double d = raw[atom * frame_count + sample + 1] -
                             raw[atom * frame_count + sample];
            const std::size_t offset = atom * delta_count + sample;
            delta[offset] = d;
            abs_delta[offset] = std::abs(d);
            delta_squared[offset] = d * d;
            dxdt[offset] = d /
                (frame_times[sample + 1] - frame_times[sample]);
        }
    }
    const std::vector<std::size_t> delta_frames(
        frame_indices.begin() + 1, frame_indices.end());
    check_channel("sasa_delta", delta, delta_count, delta_frames);
    check_channel("sasa_abs_delta", abs_delta, delta_count, delta_frames);
    check_channel("sasa_delta_squared", delta_squared, delta_count,
                  delta_frames);
    check_channel("sasa_dxdt", dxdt, delta_count, delta_frames);

    const auto rms = ReadH5Flat<double>(h5_path, w_group + "rms_delta");
    for (std::size_t atom = 0; atom < atom_count; ++atom) {
        double mean_delta_squared = 0.0;
        for (std::size_t sample = 0; sample < delta_count; ++sample)
            mean_delta_squared +=
                delta_squared[atom * delta_count + sample];
        mean_delta_squared /= static_cast<double>(delta_count);
        EXPECT_TRUE(Near(rms[atom], std::sqrt(mean_delta_squared),
                         1.0e-10, 2.0e-12));
    }

    for (const std::pair<std::string, std::string>& alias : {
             std::pair<std::string, std::string>{"mean", "sasa_mean"},
             {"std", "sasa_std"}, {"min", "sasa_min"},
             {"max", "sasa_max"}, {"delta_mean", "sasa_delta_mean"},
             {"delta_std", "sasa_delta_std"}}) {
        EXPECT_EQ(ReadH5Flat<double>(h5_path, w_group + alias.first),
                  ReadH5Flat<double>(h5_path, w_group + alias.second));
    }

    std::vector<std::size_t> npy_dimensions;
    const auto npy = ReadFloat64Npy(
        npy_directory / "sasa_welford.npy", &npy_dimensions);
    ASSERT_EQ(npy_dimensions,
              (std::vector<std::size_t>{atom_count, 7u}));
    ExpectNumpyRejectsPickleRequirement(
        npy_directory / "sasa_welford.npy");
    const std::array<std::string, 6> h5_columns = {
        "sasa_mean", "sasa_std", "sasa_min", "sasa_max",
        "sasa_delta_mean", "sasa_delta_std"};
    for (std::size_t column = 0; column < h5_columns.size(); ++column) {
        const auto expected = ReadH5Flat<double>(
            h5_path, w_group + h5_columns[column]);
        for (std::size_t atom = 0; atom < atom_count; ++atom)
            EXPECT_DOUBLE_EQ(npy[atom * 7 + column], expected[atom]);
    }
    for (std::size_t atom = 0; atom < atom_count; ++atom)
        EXPECT_DOUBLE_EQ(npy[atom * 7 + 6],
                         static_cast<double>(frame_count));
}

nmr::SolventEnvironment BuildWaterSolventFrame(const Vec3& atom_position,
                                                std::size_t frame) {
    nmr::SolventEnvironment solvent;
    auto add = [&](const Vec3& oxygen_offset,
                   const Vec3& h1_offset,
                   const Vec3& h2_offset) {
        nmr::WaterMolecule water;
        water.O_pos = atom_position + oxygen_offset;
        water.H1_pos = water.O_pos + h1_offset;
        water.H2_pos = water.O_pos + h2_offset;
        water.O_charge = -0.834;
        water.H_charge = 0.417;
        solvent.waters.push_back(water);
        solvent.water_O_positions.push_back(water.O_pos);
    };
    const double f = static_cast<double>(frame);
    add(Vec3(2.20 + 0.08 * f, 0.42 - 0.03 * f, -0.31 + 0.02 * f),
        Vec3(0.72, 0.42, 0.11), Vec3(-0.31, 0.81, -0.09));
    add(Vec3(-1.45 + 0.04 * f, 2.28 + 0.06 * f, 0.73 - 0.01 * f),
        Vec3(0.34, -0.74, 0.29), Vec3(-0.66, -0.21, 0.48));
    add(Vec3(3.72 - 0.05 * f, -1.16 + 0.02 * f, 0.91 + 0.04 * f),
        Vec3(-0.44, 0.69, -0.22), Vec3(0.58, 0.51, 0.31));
    return solvent;
}

bool AttachWaterAndHydrationOwners(
        ProteinConformation& conf,
        const nmr::SolventEnvironment& solvent,
        const Vec3& sasa_normal_input) {
    auto geometry = nmr::GeometryResult::Compute(conf);
    if (!geometry || !conf.AttachResult(std::move(geometry))) return false;
    auto spatial = nmr::SpatialIndexResult::Compute(conf);
    if (!spatial || !conf.AttachResult(std::move(spatial))) return false;
    auto sasa = nmr::SasaResult::Compute(conf);
    if (!sasa || !conf.AttachResult(std::move(sasa))) return false;

    // HydrationGeometry consumes the upstream SASA normal.  This assembly
    // test supplies a deterministically transformed upstream input so it
    // isolates the real HydrationGeometry + Welford paths; finite-grid SASA
    // covariance and its recorded bound are tested separately (DC-SA).
    conf.MutableAtomAt(0).sasa_normal = sasa_normal_input;
    auto water = nmr::WaterFieldResult::Compute(conf, solvent);
    if (!water || !conf.AttachResult(std::move(water))) return false;
    auto hydration = nmr::HydrationGeometryResult::Compute(conf, solvent);
    return hydration && conf.AttachResult(std::move(hydration));
}

std::vector<double> ReadAssembledXyzMean(
        const std::string& path,
        const std::string& group,
        const std::string& prefix) {
    const auto x = ReadH5Flat<double>(path, group + "/" + prefix + "_x_mean");
    const auto y = ReadH5Flat<double>(path, group + "/" + prefix + "_y_mean");
    const auto z = ReadH5Flat<double>(path, group + "/" + prefix + "_z_mean");
    EXPECT_EQ(y.size(), x.size());
    EXPECT_EQ(z.size(), x.size());
    std::vector<double> out(x.size() * 3);
    for (std::size_t atom = 0; atom < x.size(); ++atom) {
        out[atom * 3] = x[atom];
        out[atom * 3 + 1] = y[atom];
        out[atom * 3 + 2] = z[atom];
    }
    return out;
}

}  // namespace


TEST(DirectionalTrajectoryCovariance,
     SasaFiniteGridTimeSeriesAndWelfordRerunO3SerializedH5AndNpy) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    constexpr std::uint64_t kTransformSeed = 0x5A5A2026ULL;
    // The finite 92-point lab grid changes at most five exposed O-sphere
    // samples on this four-frame forcing sequence. One O sample is
    // 4*pi*(1.52+1.4)^2/92 = 1.164629373981218 A^2, so the observed maximum
    // is 5.82314686990609 A^2. The recorded 6 A^2 bound is the next concise
    // envelope and names the quadrature cause; it is not an analytic-kernel
    // tolerance.
    constexpr double kRecordedSasaAreaAbsA2 = 6.0;

    auto run = [&](const std::string& suffix,
                   const OrthogonalTransform* transform) {
        auto tp = TrajectoryProtein::CreateForTesting(
            BuildSasaTrajectoryProtein());
        EXPECT_NE(tp, nullptr);
        if (!tp) return std::pair<std::string, std::filesystem::path>{};
        const auto source_frames = BuildSasaMotionFrames(
            tp->CanonicalConformation().Positions());
        auto time_series =
            nmr::SasaTimeSeriesTrajectoryResult::Create(*tp);
        auto welford = nmr::SasaWelfordTrajectoryResult::Create(*tp);
        EXPECT_NE(time_series, nullptr);
        EXPECT_NE(welford, nullptr);
        Trajectory dummy("", "", "");
        for (std::size_t frame = 0; frame < source_frames.size(); ++frame) {
            const std::vector<Vec3> positions = transform
                ? nmr::test::directional::Positions(
                      *transform, source_frames[frame])
                : source_frames[frame];
            auto conf = tp->TickConformation(positions);
            EXPECT_TRUE(AttachSasaSource(*conf));
            time_series->Compute(*conf, *tp, dummy,
                                 50u + frame, 0.5 * frame);
            welford->Compute(*conf, *tp, dummy,
                             50u + frame, 0.5 * frame);
        }
        time_series->Finalize(*tp, dummy);
        welford->Finalize(*tp, dummy);
        const std::string h5_path = FreshH5Path(
            "directional_sasa_" + suffix + ".h5");
        {
            HighFive::File file(h5_path, HighFive::File::Truncate);
            time_series->WriteH5Group(*tp, file);
            welford->WriteH5Group(*tp, file);
        }
        const std::filesystem::path npy_directory =
            nmr::test::TestEnvironment::TempPath(
                "directional_sasa_" + suffix + "_npy");
        std::error_code error;
        std::filesystem::create_directories(npy_directory, error);
        EXPECT_FALSE(error) << error.message();
        EXPECT_EQ(welford->WriteFeatures(*tp, npy_directory.string()), 1);
        CheckSasaWelfordSerialization(
            h5_path, npy_directory, tp->AtomCount(), source_frames.size());
        return std::make_pair(h5_path, npy_directory);
    };

    const auto source = run("source", nullptr);
    ASSERT_FALSE(source.first.empty());
    const auto source_raw = ReadH5Flat<double>(
        source.first, "/trajectory/sasa_time_series/sasa");
    const auto source_npy = ReadFloat64Npy(
        source.second / "sasa_welford.npy");
    std::size_t changed_samples = 0;
    std::size_t changed_welford_values = 0;
    for (const bool improper : {false, true}) {
        const OrthogonalTransform transform =
            SeededTransform(kTransformSeed, improper);
        const auto moved = run(improper ? "improper" : "proper", &transform);
        ASSERT_FALSE(moved.first.empty());
        const auto moved_raw = ReadH5Flat<double>(
            moved.first, "/trajectory/sasa_time_series/sasa");
        ASSERT_EQ(moved_raw.size(), source_raw.size());
        for (std::size_t i = 0; i < source_raw.size(); ++i) {
            ASSERT_TRUE(std::isfinite(source_raw[i]));
            ASSERT_TRUE(std::isfinite(moved_raw[i]));
            const double difference = std::abs(moved_raw[i] - source_raw[i]);
            EXPECT_LE(difference, kRecordedSasaAreaAbsA2)
                << "finite 92-point lab-grid SASA envelope index=" << i;
            if (difference > 1.0e-12) ++changed_samples;
        }
        const auto moved_npy = ReadFloat64Npy(
            moved.second / "sasa_welford.npy");
        ASSERT_EQ(moved_npy.size(), source_npy.size());
        for (std::size_t row = 0; row < source_npy.size() / 7u; ++row) {
            for (std::size_t column = 0; column < 6u; ++column) {
                if (std::abs(moved_npy[row * 7u + column] -
                             source_npy[row * 7u + column]) > 1.0e-12) {
                    ++changed_welford_values;
                }
            }
            EXPECT_DOUBLE_EQ(moved_npy[row * 7u + 6u],
                             source_npy[row * 7u + 6u]);
        }
        for (const std::string& counter : {
                 "n_frames_per_atom", "delta_n_per_atom",
                 "dxdt_n_per_atom"}) {
            ExpectExactH5<std::size_t>(
                source.first, moved.first,
                "/trajectory/sasa_welford/" + counter,
                {source_npy.size() / 7u});
        }
        RemoveH5(moved.first);
        EXPECT_EQ(std::remove(
                      (moved.second / "sasa_welford.npy").string().c_str()),
                  0);
        EXPECT_EQ(::rmdir(moved.second.string().c_str()), 0);
    }
    EXPECT_GT(changed_samples, 0u)
        << "finite lab-fixed SASA grid forcing must be non-vacuous";
    EXPECT_GT(changed_welford_values, 0u)
        << "finite-grid orientation dependence must reach serialized Welford values";
    RemoveH5(source.first);
    EXPECT_EQ(std::remove(
                  (source.second / "sasa_welford.npy").string().c_str()),
              0);
    EXPECT_EQ(::rmdir(source.second.string().c_str()), 0);
}


TEST(DirectionalTrajectoryCovariance,
     PositionsTimeSeriesRerunO3SerializedH5) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto source_protein = LoadUbqOrSkip();
    if (!source_protein)
        GTEST_SKIP() << "1UBQ protonated fixture unavailable";
    auto source_tp = TrajectoryProtein::CreateForTesting(
        std::move(source_protein));
    ASSERT_NE(source_tp, nullptr);
    auto source_result =
        nmr::PositionsTimeSeriesTrajectoryResult::Create(*source_tp);
    ASSERT_NE(source_result, nullptr);
    Trajectory dummy("", "", "");
    source_result->Compute(source_tp->CanonicalConformation(), *source_tp,
                           dummy, 17, 3.25);
    source_result->Finalize(*source_tp, dummy);
    const std::string source_path =
        FreshH5Path("directional_positions_source.h5");
    {
        HighFive::File file(source_path, HighFive::File::Truncate);
        source_result->WriteH5Group(*source_tp, file);
    }
    std::vector<std::size_t> source_dims;
    const std::vector<double> source_xyz = ReadH5Flat<double>(
        source_path, "/trajectory/positions/xyz", &source_dims);
    ASSERT_EQ(source_dims,
              (std::vector<std::size_t>{source_tp->AtomCount(), 1, 3}));
    const std::string positions_dataset = "/trajectory/positions/xyz";
    EXPECT_EQ(ReadH5StringAttribute(
                  source_path, positions_dataset, "physical_class"),
              "cartesian_position");
    EXPECT_EQ(ReadH5StringAttribute(
                  source_path, positions_dataset, "coordinate_frame"),
              "conformation_cartesian_xyz");
    EXPECT_EQ(ReadH5StringAttribute(
                  source_path, positions_dataset, "parity"),
              "odd");
    EXPECT_EQ(ReadH5StringAttribute(
                  source_path, positions_dataset, "transformation"),
              "affine_position: p'=R p+t");
    EXPECT_EQ(ReadH5StringAttribute(
                  source_path, positions_dataset, "units"),
              "Angstrom");
    const std::vector<std::size_t> source_indices =
        ReadH5Flat<std::size_t>(
            source_path, "/trajectory/positions/frame_indices");
    const std::vector<double> source_times = ReadH5Flat<double>(
        source_path, "/trajectory/positions/frame_times");
    ASSERT_EQ(source_indices, (std::vector<std::size_t>{17}));
    ASSERT_EQ(source_times, (std::vector<double>{3.25}));

    for (const bool improper : {false, true}) {
        auto moved_protein = LoadUbqOrSkip();
        ASSERT_NE(moved_protein, nullptr);
        auto moved_tp = TrajectoryProtein::CreateForTesting(
            std::move(moved_protein));
        ASSERT_NE(moved_tp, nullptr);
        ASSERT_EQ(moved_tp->AtomCount(), source_tp->AtomCount());
        const OrthogonalTransform x =
            SeededTransform(0xB051710ULL, improper);
        auto moved_conf = moved_tp->TickConformation(
            nmr::test::directional::Positions(
                x, source_tp->CanonicalConformation().Positions()));
        auto moved_result =
            nmr::PositionsTimeSeriesTrajectoryResult::Create(*moved_tp);
        ASSERT_NE(moved_result, nullptr);
        moved_result->Compute(*moved_conf, *moved_tp, dummy, 17, 3.25);
        moved_result->Finalize(*moved_tp, dummy);
        const std::string moved_path = FreshH5Path(
            improper ? "directional_positions_improper.h5" :
                       "directional_positions_proper.h5");
        {
            HighFive::File file(moved_path, HighFive::File::Truncate);
            moved_result->WriteH5Group(*moved_tp, file);
        }
        std::vector<std::size_t> moved_dims;
        const std::vector<double> moved_xyz = ReadH5Flat<double>(
            moved_path, "/trajectory/positions/xyz", &moved_dims);
        ASSERT_EQ(moved_dims, source_dims);
        for (const std::pair<std::string, std::string>& attribute : {
                 std::pair<std::string, std::string>{
                     "physical_class", "cartesian_position"},
                 {"coordinate_frame", "conformation_cartesian_xyz"},
                 {"parity", "odd"},
                 {"transformation", "affine_position: p'=R p+t"},
                 {"units", "Angstrom"}}) {
            EXPECT_EQ(ReadH5StringAttribute(
                          moved_path, positions_dataset, attribute.first),
                      attribute.second);
        }
        EXPECT_EQ(ReadH5Flat<std::size_t>(
                      moved_path,
                      "/trajectory/positions/frame_indices"),
                  source_indices);
        EXPECT_EQ(ReadH5Flat<double>(
                      moved_path,
                      "/trajectory/positions/frame_times"),
                  source_times);
        for (std::size_t atom = 0; atom < source_tp->AtomCount(); ++atom) {
            const Vec3 source(source_xyz[atom * 3],
                              source_xyz[atom * 3 + 1],
                              source_xyz[atom * 3 + 2]);
            const Vec3 actual(moved_xyz[atom * 3],
                              moved_xyz[atom * 3 + 1],
                              moved_xyz[atom * 3 + 2]);
            const Vec3 expected = Position(x, source);
            EXPECT_LE((actual - expected).norm(), 2.0e-12)
                << "/trajectory/positions/xyz atom=" << atom;
        }
        RemoveH5(moved_path);
    }
    RemoveH5(source_path);
}


TEST(DirectionalTrajectoryCovariance,
     ReorientationalBodyTensorRerunO3SerializedH5) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto source_protein = LoadUbqOrSkip();
    if (!source_protein)
        GTEST_SKIP() << "1UBQ protonated fixture unavailable";
    auto source_tp = TrajectoryProtein::CreateForTesting(
        std::move(source_protein));
    ASSERT_NE(source_tp, nullptr);
    const auto frames = BuildInternalMotionFrames(
        source_tp->CanonicalConformation().Positions());
    auto source_result =
        nmr::ReorientationalDynamicsTrajectoryResult::Create(*source_tp);
    ASSERT_NE(source_result, nullptr);
    ASSERT_GT(source_result->NumVectors(), 0u);
    Trajectory dummy("", "", "");
    for (std::size_t frame = 0; frame < frames.size(); ++frame) {
        auto conf = source_tp->TickConformation(frames[frame]);
        source_result->Compute(*conf, *source_tp, dummy,
                               101 + frame, 0.5 * frame);
    }
    source_result->Finalize(*source_tp, dummy);
    const std::string source_path =
        FreshH5Path("directional_reorientation_source.h5");
    {
        HighFive::File file(source_path, HighFive::File::Truncate);
        source_result->WriteH5Group(*source_tp, file);
    }
    const std::string group = "/trajectory/reorientational_dynamics/";
    std::vector<std::size_t> tensor_dims;
    const std::vector<double> source_tensor = ReadH5Flat<double>(
        source_path, group + "bond_orientation_tensor", &tensor_dims);
    ASSERT_EQ(tensor_dims.size(), 3u);
    ASSERT_EQ(tensor_dims[0], source_result->NumVectors());
    ASSERT_EQ(tensor_dims[1], 3u);
    ASSERT_EQ(tensor_dims[2], 3u);
    ExpectReorientationTensorDirectionalMetadata(source_path);
    const std::vector<double> source_s2 = ReadH5Flat<double>(
        source_path, group + "order_parameter_S2");
    const std::vector<double> source_internal = ReadH5Flat<double>(
        source_path, group + "bond_vector_autocorrelation");
    const std::vector<double> source_lab = ReadH5Flat<double>(
        source_path, group + "bond_vector_autocorrelation_lab");
    const std::vector<std::uint8_t> source_kind =
        ReadH5Flat<std::uint8_t>(source_path, group + "vector_kind");
    const std::vector<std::int32_t> source_tail =
        ReadH5Flat<std::int32_t>(source_path, group + "tail_atom");
    const std::vector<std::int32_t> source_head =
        ReadH5Flat<std::int32_t>(source_path, group + "head_atom");
    const std::vector<std::int32_t> source_owner =
        ReadH5Flat<std::int32_t>(source_path, group + "owning_atom");
    const std::vector<std::int32_t> source_residue =
        ReadH5Flat<std::int32_t>(source_path, group + "residue_index");

    // Recorded SVD/Kabsch bounds for this six-frame forcing sequence.  The
    // common O(3) transform is conjugated through the proper Kabsch rotation;
    // SVD roundoff, rather than a scientific approximation, sets these
    // sub-nanometric tensor/TCF tolerances.
    constexpr double kBodyTensorAbs = 8.0e-11;
    constexpr double kBodyTensorRel = 2.0e-10;
    constexpr double kTcfAbs = 8.0e-11;

    for (const bool improper : {false, true}) {
        auto moved_protein = LoadUbqOrSkip();
        ASSERT_NE(moved_protein, nullptr);
        auto moved_tp = TrajectoryProtein::CreateForTesting(
            std::move(moved_protein));
        ASSERT_NE(moved_tp, nullptr);
        const OrthogonalTransform x =
            SeededTransform(0xAE071E17ULL, improper);
        auto moved_result =
            nmr::ReorientationalDynamicsTrajectoryResult::Create(*moved_tp);
        ASSERT_NE(moved_result, nullptr);
        ASSERT_EQ(moved_result->NumVectors(), source_result->NumVectors());
        for (std::size_t frame = 0; frame < frames.size(); ++frame) {
            auto conf = moved_tp->TickConformation(
                nmr::test::directional::Positions(x, frames[frame]));
            moved_result->Compute(*conf, *moved_tp, dummy,
                                  101 + frame, 0.5 * frame);
        }
        moved_result->Finalize(*moved_tp, dummy);
        const std::string moved_path = FreshH5Path(
            improper ? "directional_reorientation_improper.h5" :
                       "directional_reorientation_proper.h5");
        {
            HighFive::File file(moved_path, HighFive::File::Truncate);
            moved_result->WriteH5Group(*moved_tp, file);
        }
        std::vector<std::size_t> moved_dims;
        const std::vector<double> moved_tensor = ReadH5Flat<double>(
            moved_path, group + "bond_orientation_tensor", &moved_dims);
        ASSERT_EQ(moved_dims, tensor_dims);
        ExpectReorientationTensorDirectionalMetadata(moved_path);
        EXPECT_EQ(ReadH5Flat<std::uint8_t>(
                      moved_path, group + "vector_kind"),
                  source_kind);
        EXPECT_EQ(ReadH5Flat<std::int32_t>(
                      moved_path, group + "tail_atom"),
                  source_tail);
        EXPECT_EQ(ReadH5Flat<std::int32_t>(
                      moved_path, group + "head_atom"),
                  source_head);
        EXPECT_EQ(ReadH5Flat<std::int32_t>(
                      moved_path, group + "owning_atom"),
                  source_owner);
        EXPECT_EQ(ReadH5Flat<std::int32_t>(
                      moved_path, group + "residue_index"),
                  source_residue);

        for (std::size_t row = 0; row < tensor_dims[0]; ++row) {
            Mat3 source = Mat3::Zero();
            Mat3 actual = Mat3::Zero();
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    source(i, j) = source_tensor[row * 9 + i * 3 + j];
                    actual(i, j) = moved_tensor[row * 9 + i * 3 + j];
                }
            }
            const Mat3 expected = EvenRank2(x, source);
            EXPECT_TRUE(NearMatrix(actual, expected,
                                   kBodyTensorAbs, kBodyTensorRel))
                << "/trajectory/reorientational_dynamics/"
                   "bond_orientation_tensor row=" << row
                << " error=" << (actual - expected).norm();
            EXPECT_LE((actual - actual.transpose()).norm(),
                      kBodyTensorAbs);
            EXPECT_NEAR(actual.trace(), 1.0, kBodyTensorAbs);
        }
        ExpectSameIncludingNaN(
            ReadH5Flat<double>(moved_path,
                               group + "order_parameter_S2"),
            source_s2, kTcfAbs, 2.0e-10,
            "order_parameter_S2");
        ExpectSameIncludingNaN(
            ReadH5Flat<double>(moved_path,
                               group + "bond_vector_autocorrelation"),
            source_internal, kTcfAbs, 2.0e-10,
            "bond_vector_autocorrelation");
        ExpectSameIncludingNaN(
            ReadH5Flat<double>(moved_path,
                               group + "bond_vector_autocorrelation_lab"),
            source_lab, kTcfAbs, 2.0e-10,
            "bond_vector_autocorrelation_lab");
        RemoveH5(moved_path);
    }
    RemoveH5(source_path);
}


TEST(DirectionalTrajectoryCovariance,
     KernelOwnerWelfordAssembledMeansRerunO3SerializedH5) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto source_protein = LoadUbqOrSkip();
    if (!source_protein)
        GTEST_SKIP() << "1UBQ protonated fixture unavailable";
    auto source_tp = TrajectoryProtein::CreateForTesting(
        std::move(source_protein));
    ASSERT_NE(source_tp, nullptr);
    const auto frames = BuildInternalMotionFrames(
        source_tp->CanonicalConformation().Positions());
    auto source_bs = nmr::BsWelfordTrajectoryResult::Create(*source_tp);
    auto source_hm = nmr::HmWelfordTrajectoryResult::Create(*source_tp);
    auto source_mc = nmr::McConnellWelfordTrajectoryResult::Create(*source_tp);
    auto source_bs_ts =
        nmr::BsShieldingTimeSeriesTrajectoryResult::Create(*source_tp);
    auto source_hm_ts =
        nmr::HmShieldingTimeSeriesTrajectoryResult::Create(*source_tp);
    auto source_mc_ts =
        nmr::McConnellShieldingTimeSeriesTrajectoryResult::Create(*source_tp);
    ASSERT_NE(source_bs, nullptr);
    ASSERT_NE(source_hm, nullptr);
    ASSERT_NE(source_mc, nullptr);
    ASSERT_NE(source_bs_ts, nullptr);
    ASSERT_NE(source_hm_ts, nullptr);
    ASSERT_NE(source_mc_ts, nullptr);
    Trajectory dummy("", "", "");
    constexpr std::size_t kFrames = 3;
    for (std::size_t frame = 0; frame < kFrames; ++frame) {
        auto conf = source_tp->TickConformation(frames[frame]);
        ASSERT_TRUE(AttachKernelOwners(*conf));
        source_bs->Compute(*conf, *source_tp, dummy,
                           211 + frame, 0.25 * frame);
        source_hm->Compute(*conf, *source_tp, dummy,
                           211 + frame, 0.25 * frame);
        source_mc->Compute(*conf, *source_tp, dummy,
                           211 + frame, 0.25 * frame);
        source_bs_ts->Compute(*conf, *source_tp, dummy,
                              211 + frame, 0.25 * frame);
        source_hm_ts->Compute(*conf, *source_tp, dummy,
                              211 + frame, 0.25 * frame);
        source_mc_ts->Compute(*conf, *source_tp, dummy,
                              211 + frame, 0.25 * frame);
    }
    source_bs->Finalize(*source_tp, dummy);
    source_hm->Finalize(*source_tp, dummy);
    source_mc->Finalize(*source_tp, dummy);
    source_bs_ts->Finalize(*source_tp, dummy);
    source_hm_ts->Finalize(*source_tp, dummy);
    source_mc_ts->Finalize(*source_tp, dummy);
    const std::string source_path =
        FreshH5Path("directional_kernel_welford_source.h5");
    {
        HighFive::File file(source_path, HighFive::File::Truncate);
        source_bs->WriteH5Group(*source_tp, file);
        source_hm->WriteH5Group(*source_tp, file);
        source_mc->WriteH5Group(*source_tp, file);
        source_bs_ts->WriteH5Group(*source_tp, file);
        source_hm_ts->WriteH5Group(*source_tp, file);
        source_mc_ts->WriteH5Group(*source_tp, file);
    }
    for (const std::string& group : {
             "/trajectory/bs_welford",
             "/trajectory/hm_welford",
             "/trajectory/mc_welford"}) {
        ExpectKernelWelfordDirectionalMetadata(source_path, group);
        ExpectNonClosedComponentStatisticsMetadata(
            source_path, group, {"t1", "t2"});
        const auto counts = ReadH5Flat<std::size_t>(
            source_path, group + "/n_frames_per_atom");
        ASSERT_EQ(counts.size(), source_tp->AtomCount());
        std::size_t conditionally_omitted = 0;
        for (std::size_t count : counts) {
            if (group == "/trajectory/mc_welford") {
                EXPECT_LE(count, kFrames);
                if (count < kFrames) ++conditionally_omitted;
            } else {
                EXPECT_EQ(count, kFrames);
            }
        }
        if (group == "/trajectory/mc_welford") {
            EXPECT_GT(conditionally_omitted, 0u);
        }
    }

    for (const bool improper : {false, true}) {
        auto moved_protein = LoadUbqOrSkip();
        ASSERT_NE(moved_protein, nullptr);
        auto moved_tp = TrajectoryProtein::CreateForTesting(
            std::move(moved_protein));
        ASSERT_NE(moved_tp, nullptr);
        const OrthogonalTransform x =
            SeededTransform(0x0E1F04DULL, improper);
        auto moved_bs = nmr::BsWelfordTrajectoryResult::Create(*moved_tp);
        auto moved_hm = nmr::HmWelfordTrajectoryResult::Create(*moved_tp);
        auto moved_mc =
            nmr::McConnellWelfordTrajectoryResult::Create(*moved_tp);
        auto moved_bs_ts =
            nmr::BsShieldingTimeSeriesTrajectoryResult::Create(*moved_tp);
        auto moved_hm_ts =
            nmr::HmShieldingTimeSeriesTrajectoryResult::Create(*moved_tp);
        auto moved_mc_ts =
            nmr::McConnellShieldingTimeSeriesTrajectoryResult::Create(
                *moved_tp);
        ASSERT_NE(moved_bs, nullptr);
        ASSERT_NE(moved_hm, nullptr);
        ASSERT_NE(moved_mc, nullptr);
        ASSERT_NE(moved_bs_ts, nullptr);
        ASSERT_NE(moved_hm_ts, nullptr);
        ASSERT_NE(moved_mc_ts, nullptr);
        for (std::size_t frame = 0; frame < kFrames; ++frame) {
            auto conf = moved_tp->TickConformation(
                nmr::test::directional::Positions(x, frames[frame]));
            ASSERT_TRUE(AttachKernelOwners(*conf));
            moved_bs->Compute(*conf, *moved_tp, dummy,
                              211 + frame, 0.25 * frame);
            moved_hm->Compute(*conf, *moved_tp, dummy,
                              211 + frame, 0.25 * frame);
            moved_mc->Compute(*conf, *moved_tp, dummy,
                              211 + frame, 0.25 * frame);
            moved_bs_ts->Compute(*conf, *moved_tp, dummy,
                                 211 + frame, 0.25 * frame);
            moved_hm_ts->Compute(*conf, *moved_tp, dummy,
                                 211 + frame, 0.25 * frame);
            moved_mc_ts->Compute(*conf, *moved_tp, dummy,
                                 211 + frame, 0.25 * frame);
        }
        moved_bs->Finalize(*moved_tp, dummy);
        moved_hm->Finalize(*moved_tp, dummy);
        moved_mc->Finalize(*moved_tp, dummy);
        moved_bs_ts->Finalize(*moved_tp, dummy);
        moved_hm_ts->Finalize(*moved_tp, dummy);
        moved_mc_ts->Finalize(*moved_tp, dummy);
        const std::string moved_path = FreshH5Path(
            improper ? "directional_kernel_welford_improper.h5" :
                       "directional_kernel_welford_proper.h5");
        {
            HighFive::File file(moved_path, HighFive::File::Truncate);
            moved_bs->WriteH5Group(*moved_tp, file);
            moved_hm->WriteH5Group(*moved_tp, file);
            moved_mc->WriteH5Group(*moved_tp, file);
            moved_bs_ts->WriteH5Group(*moved_tp, file);
            moved_hm_ts->WriteH5Group(*moved_tp, file);
            moved_mc_ts->WriteH5Group(*moved_tp, file);
        }

        for (const std::string& group : {
                 "/trajectory/bs_welford",
                 "/trajectory/hm_welford",
                 "/trajectory/mc_welford"}) {
            ExpectKernelWelfordDirectionalMetadata(moved_path, group);
        }

        // These are means assembled from per-component Welford channels
        // after each real owner was rerun on every transformed frame.
        ExpectTensorWelfordMeans(
            source_path, moved_path, "/trajectory/bs_welford", x,
            3.0e-12, 3.0e-10);
        ExpectTensorWelfordMeans(
            source_path, moved_path, "/trajectory/hm_welford", x,
            3.0e-11, 1.0e-8);
        ExpectTensorWelfordMeans(
            source_path, moved_path, "/trajectory/mc_welford", x,
            3.0e-12, 3.0e-10);

        // Raw owner timelines use exactly PackFull9 at the serialized H5
        // boundary. The static BS/HM/MC owners above were rerun on every
        // transformed frame; no already-produced tensor is re-used.
        ExpectFullTensorTimeSeries(
            source_path, moved_path,
            "/trajectory/bs_shielding_time_series/xyz", x,
            source_tp->AtomCount(), kFrames, 3.0e-12, 3.0e-10);
        ExpectFullTensorTimeSeries(
            source_path, moved_path,
            "/trajectory/hm_shielding_time_series/xyz", x,
            source_tp->AtomCount(), kFrames, 3.0e-11, 1.0e-8);
        ExpectFullTensorTimeSeries(
            source_path, moved_path,
            "/trajectory/mc_shielding_time_series/xyz", x,
            source_tp->AtomCount(), kFrames, 3.0e-12, 3.0e-10);
        for (const std::string& raw_group : {
                 "/trajectory/bs_shielding_time_series",
                 "/trajectory/hm_shielding_time_series",
                 "/trajectory/mc_shielding_time_series"}) {
            ExpectExactH5<std::size_t>(
                source_path, moved_path, raw_group + "/frame_indices",
                {kFrames});
            ExpectExactH5<double>(
                source_path, moved_path, raw_group + "/frame_times",
                {kFrames});
        }
        RemoveH5(moved_path);
    }
    RemoveH5(source_path);
}


TEST(DirectionalTrajectoryCovariance,
     WaterHydrationOwnerWelfordAssembledMeansRerunO3SerializedH5) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    const Vec3 base_position(0.35, -0.27, 0.19);
    auto source_tp = TrajectoryProtein::CreateForTesting(
        BuildOneAtomProtein(base_position));
    ASSERT_NE(source_tp, nullptr);
    auto source_water =
        nmr::WaterFieldWelfordTrajectoryResult::Create(*source_tp);
    auto source_hydration =
        nmr::HydrationGeometryWelfordTrajectoryResult::Create(*source_tp);
    auto source_water_ts =
        nmr::WaterFieldTimeSeriesTrajectoryResult::Create(*source_tp);
    auto source_hydration_ts =
        nmr::HydrationGeometryTimeSeriesTrajectoryResult::Create(*source_tp);
    ASSERT_NE(source_water, nullptr);
    ASSERT_NE(source_hydration, nullptr);
    ASSERT_NE(source_water_ts, nullptr);
    ASSERT_NE(source_hydration_ts, nullptr);
    Trajectory dummy("", "", "");
    constexpr std::size_t kFrames = 3;
    std::array<Vec3, kFrames> positions;
    std::array<Vec3, kFrames> normals;
    std::array<nmr::SolventEnvironment, kFrames> solvents;
    for (std::size_t frame = 0; frame < kFrames; ++frame) {
        const double f = static_cast<double>(frame);
        positions[frame] = base_position +
            Vec3(0.025 * f, -0.018 * f, 0.011 * f);
        normals[frame] = Vec3(0.31 + 0.03 * f,
                              -0.44 + 0.02 * f,
                              0.842 - 0.01 * f).normalized();
        solvents[frame] = BuildWaterSolventFrame(positions[frame], frame);
        auto conf = source_tp->TickConformation({positions[frame]});
        ASSERT_TRUE(AttachWaterAndHydrationOwners(
            *conf, solvents[frame], normals[frame]));
        source_water->Compute(*conf, *source_tp, dummy,
                              307 + frame, 0.4 * frame);
        source_hydration->Compute(*conf, *source_tp, dummy,
                                  307 + frame, 0.4 * frame);
        source_water_ts->Compute(*conf, *source_tp, dummy,
                                 307 + frame, 0.4 * frame);
        source_hydration_ts->Compute(*conf, *source_tp, dummy,
                                     307 + frame, 0.4 * frame);
    }
    source_water->Finalize(*source_tp, dummy);
    source_hydration->Finalize(*source_tp, dummy);
    source_water_ts->Finalize(*source_tp, dummy);
    source_hydration_ts->Finalize(*source_tp, dummy);
    const std::string source_path =
        FreshH5Path("directional_water_hydration_welford_source.h5");
    {
        HighFive::File file(source_path, HighFive::File::Truncate);
        source_water->WriteH5Group(*source_tp, file);
        source_hydration->WriteH5Group(*source_tp, file);
        source_water_ts->WriteH5Group(*source_tp, file);
        source_hydration_ts->WriteH5Group(*source_tp, file);
    }
    ExpectNonClosedComponentStatisticsMetadata(
        source_path, "/trajectory/water_field_welford",
        {"efield_x", "efield_y", "efield_z",
         "efield_first_x", "efield_first_y", "efield_first_z",
         "efg_t2", "efg_first_t2"});
    ExpectNonClosedComponentStatisticsMetadata(
        source_path, "/trajectory/hydration_geometry_welford",
        {"dipole_vector_x", "dipole_vector_y", "dipole_vector_z",
         "surface_normal_x", "surface_normal_y", "surface_normal_z"});
    {
        HighFive::File file(source_path, HighFive::File::ReadOnly);
        auto group = file.getGroup("/trajectory/water_field_welford");
        bool t0_zero = false;
        bool t1_zero = false;
        group.getAttribute("efg_t0_structural_zero").read(t0_zero);
        group.getAttribute("efg_t1_structural_zero").read(t1_zero);
        EXPECT_TRUE(t0_zero);
        EXPECT_TRUE(t1_zero);
        auto raw_group =
            file.getGroup("/trajectory/water_field_time_series");
        t0_zero = false;
        t1_zero = false;
        raw_group.getAttribute("efg_t0_structural_zero").read(t0_zero);
        raw_group.getAttribute("efg_t1_structural_zero").read(t1_zero);
        EXPECT_TRUE(t0_zero);
        EXPECT_TRUE(t1_zero);
    }
    ExpectWaterWelfordDirectionalMetadata(source_path);
    const std::vector<std::size_t> expected_frame_indices = {307u, 308u, 309u};
    const std::vector<double> expected_frame_times = {0.0, 0.4, 0.8};
    for (const std::string& group : {
             "/trajectory/water_field_time_series",
             "/trajectory/hydration_geometry_time_series"}) {
        EXPECT_EQ(ReadH5Flat<std::size_t>(
                      source_path, group + "/frame_indices"),
                  expected_frame_indices) << group;
        EXPECT_EQ(ReadH5Flat<double>(source_path, group + "/frame_times"),
                  expected_frame_times) << group;
    }
    EXPECT_EQ(ReadH5Flat<std::size_t>(
                  source_path,
                  "/trajectory/water_field_welford/n_frames_per_atom"),
              (std::vector<std::size_t>{kFrames}));
    EXPECT_EQ(ReadH5Flat<std::uint8_t>(
                  source_path,
                  "/trajectory/water_field_welford/"
                  "source_attached_per_frame"),
              (std::vector<std::uint8_t>{1u, 1u, 1u}));
    EXPECT_EQ(ReadH5Flat<std::size_t>(
                  source_path,
                  "/trajectory/hydration_geometry_welford/"
                  "n_frames_per_atom"),
              (std::vector<std::size_t>{kFrames}));
    EXPECT_EQ(ReadH5Flat<std::size_t>(
                  source_path,
                  "/trajectory/hydration_geometry_welford/"
                  "delta_n_per_atom"),
              (std::vector<std::size_t>{kFrames - 1u}));
    EXPECT_EQ(ReadH5Flat<std::size_t>(
                  source_path,
                  "/trajectory/hydration_geometry_welford/"
                  "dxdt_n_per_atom"),
              (std::vector<std::size_t>{kFrames - 1u}));
    EXPECT_EQ(ReadH5Flat<std::uint8_t>(
                  source_path,
                  "/trajectory/hydration_geometry_welford/"
                  "source_attached_per_frame"),
              (std::vector<std::uint8_t>{1u, 1u, 1u}));

    for (const bool improper : {false, true}) {
        auto moved_tp = TrajectoryProtein::CreateForTesting(
            BuildOneAtomProtein(base_position));
        ASSERT_NE(moved_tp, nullptr);
        const OrthogonalTransform x =
            SeededTransform(0x0A7E2F1EULL, improper);
        auto moved_water =
            nmr::WaterFieldWelfordTrajectoryResult::Create(*moved_tp);
        auto moved_hydration =
            nmr::HydrationGeometryWelfordTrajectoryResult::Create(*moved_tp);
        auto moved_water_ts =
            nmr::WaterFieldTimeSeriesTrajectoryResult::Create(*moved_tp);
        auto moved_hydration_ts =
            nmr::HydrationGeometryTimeSeriesTrajectoryResult::Create(
                *moved_tp);
        ASSERT_NE(moved_water, nullptr);
        ASSERT_NE(moved_hydration, nullptr);
        ASSERT_NE(moved_water_ts, nullptr);
        ASSERT_NE(moved_hydration_ts, nullptr);
        for (std::size_t frame = 0; frame < kFrames; ++frame) {
            auto conf = moved_tp->TickConformation(
                {Position(x, positions[frame])});
            const auto moved_solvent =
                nmr::test::directional::Solvent(x, solvents[frame]);
            ASSERT_TRUE(AttachWaterAndHydrationOwners(
                *conf, moved_solvent, Polar(x, normals[frame])));
            moved_water->Compute(*conf, *moved_tp, dummy,
                                 307 + frame, 0.4 * frame);
            moved_hydration->Compute(*conf, *moved_tp, dummy,
                                     307 + frame, 0.4 * frame);
            moved_water_ts->Compute(*conf, *moved_tp, dummy,
                                    307 + frame, 0.4 * frame);
            moved_hydration_ts->Compute(*conf, *moved_tp, dummy,
                                        307 + frame, 0.4 * frame);
        }
        moved_water->Finalize(*moved_tp, dummy);
        moved_hydration->Finalize(*moved_tp, dummy);
        moved_water_ts->Finalize(*moved_tp, dummy);
        moved_hydration_ts->Finalize(*moved_tp, dummy);
        const std::string moved_path = FreshH5Path(
            improper ? "directional_water_hydration_welford_improper.h5" :
                       "directional_water_hydration_welford_proper.h5");
        {
            HighFive::File file(moved_path, HighFive::File::Truncate);
            moved_water->WriteH5Group(*moved_tp, file);
            moved_hydration->WriteH5Group(*moved_tp, file);
            moved_water_ts->WriteH5Group(*moved_tp, file);
            moved_hydration_ts->WriteH5Group(*moved_tp, file);
        }
        ExpectWaterWelfordDirectionalMetadata(moved_path);

        const std::string water_group =
            "/trajectory/water_field_welford";
        for (const std::string& prefix : {"efield", "efield_first"}) {
            const auto a = ReadAssembledXyzMean(
                source_path, water_group, prefix);
            const auto b = ReadAssembledXyzMean(
                moved_path, water_group, prefix);
            ASSERT_EQ(a.size(), 3u);
            ASSERT_EQ(b.size(), 3u);
            const Vec3 source(a[0], a[1], a[2]);
            const Vec3 actual(b[0], b[1], b[2]);
            EXPECT_TRUE(NearVector(actual, Polar(x, source),
                                   2.0e-9, 2.0e-10))
                << water_group << "/" << prefix << "_*_mean";
        }
        for (const std::string& prefix : {"efg_t2", "efg_first_t2"}) {
            const auto a = ReadH5Flat<double>(
                source_path, water_group + "/" + prefix + "_mean");
            const auto b = ReadH5Flat<double>(
                moved_path, water_group + "/" + prefix + "_mean");
            ASSERT_EQ(a.size(), 5u);
            ASSERT_EQ(b.size(), 5u);
            SphericalTensor source;
            SphericalTensor actual;
            for (int component = 0; component < 5; ++component) {
                source.T2[component] = a[component];
                actual.T2[component] = b[component];
            }
            const SphericalTensor expected = RotateNativeT2(x, source);
            for (int component = 0; component < 5; ++component) {
                EXPECT_TRUE(Near(actual.T2[component],
                                 expected.T2[component],
                                 2.0e-9, 2.0e-10));
            }
        }

        const std::string hydration_group =
            "/trajectory/hydration_geometry_welford";
        for (const std::string& prefix : {
                 "dipole_vector", "surface_normal"}) {
            const auto a = ReadAssembledXyzMean(
                source_path, hydration_group, prefix);
            const auto b = ReadAssembledXyzMean(
                moved_path, hydration_group, prefix);
            ASSERT_EQ(a.size(), 3u);
            ASSERT_EQ(b.size(), 3u);
            const Vec3 source(a[0], a[1], a[2]);
            const Vec3 actual(b[0], b[1], b[2]);
            EXPECT_TRUE(NearVector(actual, Polar(x, source),
                                   2.0e-11, 2.0e-12))
                << hydration_group << "/" << prefix << "_*_mean";
        }

        // Raw serialized timelines from the same real WaterField and
        // HydrationGeometry owners. Water EFG is native-T2-only: reconstruct
        // to Cartesian, apply Q T Q^T, and decompose back before comparing
        // the five stored components.
        const std::string water_ts_group =
            "/trajectory/water_field_time_series";
        ExpectPolarVectorTimeSeries(
            source_path, moved_path, water_ts_group + "/efield", x,
            1u, kFrames, 2.0e-9, 2.0e-10);
        ExpectPolarVectorTimeSeries(
            source_path, moved_path, water_ts_group + "/efield_first", x,
            1u, kFrames, 2.0e-9, 2.0e-10);
        ExpectNativeT2TimeSeries(
            source_path, moved_path, water_ts_group + "/efg", x,
            1u, kFrames, 2.0e-9, 2.0e-10);
        ExpectNativeT2TimeSeries(
            source_path, moved_path, water_ts_group + "/efg_first", x,
            1u, kFrames, 2.0e-9, 2.0e-10);
        ExpectExactH5<std::uint8_t>(
            source_path, moved_path,
            water_ts_group + "/source_attached_per_frame", {kFrames});
        ExpectExactH5<std::uint32_t>(
            source_path, moved_path, water_ts_group + "/n_first",
            {1u, kFrames});
        ExpectExactH5<std::uint32_t>(
            source_path, moved_path, water_ts_group + "/n_second",
            {1u, kFrames});
        ExpectExactH5<std::uint8_t>(
            source_path, moved_path,
            water_ts_group + "/efield_clamp_mask", {1u, kFrames});
        ExpectExactH5<std::uint8_t>(
            source_path, moved_path,
            water_ts_group + "/efield_first_clamp_mask", {1u, kFrames});
        ExpectExactH5<double>(
            source_path, moved_path,
            water_ts_group + "/efield_clamp_scale", {1u, kFrames});
        ExpectExactH5<double>(
            source_path, moved_path,
            water_ts_group + "/efield_first_clamp_scale", {1u, kFrames});
        ExpectExactH5<std::size_t>(
            source_path, moved_path, water_ts_group + "/frame_indices",
            {kFrames});
        ExpectExactH5<double>(
            source_path, moved_path, water_ts_group + "/frame_times",
            {kFrames});
        ExpectExactH5<std::size_t>(
            source_path, moved_path,
            water_group + "/n_frames_per_atom", {1u});
        ExpectExactH5<std::uint8_t>(
            source_path, moved_path,
            water_group + "/source_attached_per_frame", {kFrames});

        const std::string hydration_ts_group =
            "/trajectory/hydration_geometry_time_series";
        ExpectPolarVectorTimeSeries(
            source_path, moved_path,
            hydration_ts_group + "/dipole_vector", x,
            1u, kFrames, 2.0e-11, 2.0e-12);
        ExpectPolarVectorTimeSeries(
            source_path, moved_path,
            hydration_ts_group + "/surface_normal", x,
            1u, kFrames, 2.0e-11, 2.0e-12);
        ExpectExactH5<std::uint8_t>(
            source_path, moved_path,
            hydration_ts_group + "/source_attached_per_frame", {kFrames});
        ExpectExactH5<std::uint32_t>(
            source_path, moved_path,
            hydration_ts_group + "/first_shell_count", {1u, kFrames});
        ExpectExactH5<std::size_t>(
            source_path, moved_path,
            hydration_ts_group + "/frame_indices", {kFrames});
        ExpectExactH5<double>(
            source_path, moved_path,
            hydration_ts_group + "/frame_times", {kFrames});
        for (const std::string& counter : {
                 "n_frames_per_atom", "delta_n_per_atom",
                 "dxdt_n_per_atom"}) {
            ExpectExactH5<std::size_t>(
                source_path, moved_path,
                hydration_group + "/" + counter, {1u});
        }
        ExpectExactH5<std::uint8_t>(
            source_path, moved_path,
            hydration_group + "/source_attached_per_frame", {kFrames});
        for (const std::string& scalar : {
                 "half_shell_asymmetry", "dipole_alignment",
                 "dipole_coherence", "mean_net_dipole_eA",
                 "dipole_order_parameter"}) {
            const auto source_scalar = ReadH5Flat<double>(
                source_path, hydration_ts_group + "/" + scalar);
            const auto moved_scalar = ReadH5Flat<double>(
                moved_path, hydration_ts_group + "/" + scalar);
            ExpectSameIncludingNaN(
                moved_scalar, source_scalar, 2.0e-11, 2.0e-12,
                (hydration_ts_group + "/" + scalar).c_str());
        }
        RemoveH5(moved_path);
    }
    RemoveH5(source_path);
}


TEST(DirectionalTrajectoryCovariance,
     LocalBackboneCbResidualRerunProperAndImproperSerializedH5) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto source_protein = LoadUbqOrSkip();
    if (!source_protein)
        GTEST_SKIP() << "1UBQ protonated fixture unavailable";
    auto source_tp = TrajectoryProtein::CreateForTesting(
        std::move(source_protein));
    ASSERT_NE(source_tp, nullptr);
    const auto source_positions =
        source_tp->CanonicalConformation().Positions();
    auto source_result =
        nmr::LocalBackboneGeometryTrajectoryResult::Create(*source_tp);
    ASSERT_NE(source_result, nullptr);
    Trajectory dummy("", "", "");
    source_result->Compute(source_tp->CanonicalConformation(), *source_tp,
                           dummy, 401, 8.0);
    source_result->Finalize(*source_tp, dummy);
    const std::string source_path =
        FreshH5Path("directional_local_backbone_source.h5");
    {
        HighFive::File file(source_path, HighFive::File::Truncate);
        source_result->WriteH5Group(*source_tp, file);
    }

    const std::string group = "/trajectory/local_backbone_geometry";
    const std::size_t residue_count = source_tp->ProteinRef().ResidueCount();
    std::vector<std::size_t> source_dimensions;
    const auto source_vector = ReadH5Flat<double>(
        source_path, group + "/cb_local_vector", &source_dimensions);
    ASSERT_EQ(source_dimensions,
              (std::vector<std::size_t>{residue_count, 1u, 3u}));
    std::vector<std::size_t> source_deviation_dimensions;
    const auto source_deviation = ReadH5Flat<double>(
        source_path, group + "/cb_deviation", &source_deviation_dimensions);
    ASSERT_EQ(source_deviation_dimensions,
              (std::vector<std::size_t>{residue_count, 1u}));
    EXPECT_EQ(ReadH5StringAttribute(
                  source_path, group + "/cb_deviation", "coordinate_frame"),
              "intrinsic_chiral_lookup");
    EXPECT_EQ(ReadH5StringAttribute(
                  source_path, group + "/cb_deviation", "transformation"),
              "rotation-invariant under proper rotations; ideal-L-CB construction "
              "is chirality-conditioned and has no improper-transform contract");
    EXPECT_EQ(ReadH5StringAttribute(
                  source_path, group + "/cb_deviation", "parity"),
              "mixed");
    EXPECT_EQ(ReadH5StringAttribute(
                  source_path, group + "/cb_local_vector", "coordinate_frame"),
              "conformation_cartesian_xyz");
    EXPECT_EQ(ReadH5StringAttribute(
                  source_path, group + "/cb_local_vector", "transformation"),
              "polar displacement under proper rotations: v'=R v; the ideal-L-CB "
              "construction mixes bond vectors with a cross product and therefore "
              "has no single improper-transform parity");
    EXPECT_EQ(ReadH5StringAttribute(
                  source_path, group + "/cb_local_vector", "parity"),
              "mixed");

    constexpr std::uint64_t kTransformSeed = 0xCB10CA1B20260713ULL;
    constexpr double kProperAbsTolerance = 2.0e-11;
    constexpr double kProperRelTolerance = 2.0e-12;
    std::size_t improper_nonpolar_rows = 0;
    std::size_t improper_changed_deviation_rows = 0;
    for (const bool improper : {false, true}) {
        auto moved_protein = LoadUbqOrSkip();
        ASSERT_NE(moved_protein, nullptr);
        auto moved_tp = TrajectoryProtein::CreateForTesting(
            std::move(moved_protein));
        ASSERT_NE(moved_tp, nullptr);
        const OrthogonalTransform transform =
            SeededTransform(kTransformSeed, improper);
        auto moved_conf = moved_tp->TickConformation(
            nmr::test::directional::Positions(transform, source_positions));
        auto moved_result =
            nmr::LocalBackboneGeometryTrajectoryResult::Create(*moved_tp);
        ASSERT_NE(moved_result, nullptr);
        moved_result->Compute(*moved_conf, *moved_tp, dummy, 401, 8.0);
        moved_result->Finalize(*moved_tp, dummy);
        const std::string moved_path = FreshH5Path(
            improper ? "directional_local_backbone_improper.h5" :
                       "directional_local_backbone_proper.h5");
        {
            HighFive::File file(moved_path, HighFive::File::Truncate);
            moved_result->WriteH5Group(*moved_tp, file);
        }

        if (!improper) {
            // Despite the legacy `local` name, the residual is emitted in
            // global xyz and is a polar vector for SO(3) transforms.
            ExpectPolarVectorTimeSeries(
                source_path, moved_path, group + "/cb_local_vector",
                transform, residue_count, 1u,
                kProperAbsTolerance, kProperRelTolerance);
            const auto moved_deviation = ReadH5Flat<double>(
                moved_path, group + "/cb_deviation");
            ExpectSameIncludingNaN(
                moved_deviation, source_deviation,
                kProperAbsTolerance, kProperRelTolerance,
                "proper cb_deviation");
        } else {
            std::vector<std::size_t> moved_dimensions;
            const auto moved_vector = ReadH5Flat<double>(
                moved_path, group + "/cb_local_vector", &moved_dimensions);
            ASSERT_EQ(moved_dimensions, source_dimensions);
            ASSERT_EQ(moved_vector.size(), source_vector.size());
            for (std::size_t ri = 0; ri < residue_count; ++ri) {
                const Vec3 source(source_vector[ri * 3],
                                  source_vector[ri * 3 + 1],
                                  source_vector[ri * 3 + 2]);
                const Vec3 actual(moved_vector[ri * 3],
                                  moved_vector[ri * 3 + 1],
                                  moved_vector[ri * 3 + 2]);
                if (!source.allFinite()) {
                    EXPECT_TRUE(std::isnan(source.x()));
                    EXPECT_TRUE(std::isnan(source.y()));
                    EXPECT_TRUE(std::isnan(source.z()));
                    EXPECT_TRUE(std::isnan(actual.x()));
                    EXPECT_TRUE(std::isnan(actual.y()));
                    EXPECT_TRUE(std::isnan(actual.z()));
                    continue;
                }
                ASSERT_TRUE(actual.allFinite());
                if ((actual - Polar(transform, source)).norm() > 1.0e-4)
                    ++improper_nonpolar_rows;
            }
            std::vector<std::size_t> moved_deviation_dimensions;
            const auto moved_deviation = ReadH5Flat<double>(
                moved_path, group + "/cb_deviation",
                &moved_deviation_dimensions);
            ASSERT_EQ(moved_deviation_dimensions,
                      source_deviation_dimensions);
            ASSERT_EQ(moved_deviation.size(), source_deviation.size());
            for (std::size_t ri = 0; ri < residue_count; ++ri) {
                if (std::isnan(source_deviation[ri])) {
                    EXPECT_TRUE(std::isnan(moved_deviation[ri]));
                    continue;
                }
                ASSERT_TRUE(std::isfinite(source_deviation[ri]));
                ASSERT_TRUE(std::isfinite(moved_deviation[ri]));
                if (std::abs(moved_deviation[ri] - source_deviation[ri]) >
                    1.0e-6) {
                    ++improper_changed_deviation_rows;
                }
            }
        }

        // The fixed-sign cross(c,n) term constructs an L-amino-acid ideal
        // C-beta. It is SO(3)-covariant but chiral under reflection, so the
        // improper rerun above is intentionally checked only for honest
        // finite/NaN serialization and for absence of a false polar law.
        for (const std::string& mask : {
                 "has_N_CA_C", "has_CB", "has_prev_C", "has_next_N",
                 "is_glycine", "is_proline"}) {
            ExpectExactH5<std::uint8_t>(
                source_path, moved_path, group + "/" + mask,
                {residue_count});
        }
        ExpectExactH5<std::uint8_t>(
            source_path, moved_path,
            group + "/source_attached_per_frame", {1u});
        ExpectExactH5<std::size_t>(
            source_path, moved_path, group + "/frame_indices", {1u});
        ExpectExactH5<double>(
            source_path, moved_path, group + "/frame_times", {1u});
        RemoveH5(moved_path);
    }
    EXPECT_GT(improper_nonpolar_rows, 0u)
        << "improper chirality forcing must refute a polar-vector claim";
    EXPECT_GT(improper_changed_deviation_rows, 0u)
        << "improper chirality forcing must refute an O(3)-invariant scalar claim";
    RemoveH5(source_path);
}

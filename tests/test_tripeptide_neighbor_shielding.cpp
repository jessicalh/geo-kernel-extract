// Smoke test for TripeptideNeighborShieldingResult.
//
// Per Larsen 2015 Eq 3: each residue i receives Δσ_BB^{i-1}(i) +
// Δσ_BB^{i+1}(i) — read at the flanking ALA cap atoms of the
// (i±1)-centered tripeptides, with the AAA reference at standard
// angles (φ_std=-120°, ψ_std=140°) subtracted.
//
// On 1UBQ_pm6dh3plus.pdb we expect:
//   - residues_any ≥ 60 (almost every internal residue receives at
//     least one direction's contribution)
//   - atoms_accumulated > 200 (~7 cap atoms per residue × 60 residues)
//   - residual_vec_prev / residual_vec_next populated on most BB
//     atoms whose contributing direction hit
//   - SER residues' i-1 / i+1 directions, if SER itself is i±1, would
//     hit orca_input_orientation rows; assert at least one
//     ResidueMatch has frame_type=orca_input_orientation

#include "TestEnvironment.h"
#include <gtest/gtest.h>

#include "AminoAcidType.h"
#include "ConformationAtom.h"
#include "EnrichmentResult.h"
#include "GeometryResult.h"
#include "PdbFileReader.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "RuntimeEnvironment.h"
#include "Session.h"
#include "SpatialIndexResult.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"
#include "TripeptideDftTable.h"
#include "TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult.h"
#include "TripeptideNeighborResidualVecPrevTimeSeriesTrajectoryResult.h"
#include "TripeptideNeighborShieldingResult.h"
#include "TripeptideNeighborShieldingTimeSeriesTrajectoryResult.h"
#include "DirectionalTestHelpers.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5File.hpp>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace nmr;
namespace fs = std::filesystem;

#ifndef NMR_TEST_PYTHON_EXECUTABLE
#error "NMR_TEST_PYTHON_EXECUTABLE must be defined"
#endif
#ifndef NMR_NPY_ALLOW_PICKLE_FALSE_SCRIPT
#error "NMR_NPY_ALLOW_PICKLE_FALSE_SCRIPT must be defined"
#endif

namespace {

struct NpyArray {
    std::vector<size_t> shape;
    std::vector<char> bytes;
};

std::string Trim(std::string s) {
    auto is_space = [](unsigned char c) { return std::isspace(c) != 0; };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [&](char c) { return !is_space(c); }));
    s.erase(std::find_if(s.rbegin(), s.rend(), [&](char c) { return !is_space(c); }).base(), s.end());
    return s;
}

NpyArray ReadNpy(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    EXPECT_TRUE(in.is_open()) << path;
    NpyArray arr;
    if (!in.is_open())
        return arr;
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
        return arr;
    std::stringstream ss(header.substr(shape_begin + 1, shape_end - shape_begin - 1));
    std::string token;
    while (std::getline(ss, token, ',')) {
        token = Trim(token);
        if (!token.empty())
            arr.shape.push_back(static_cast<size_t>(std::stoull(token)));
    }
    arr.bytes.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    return arr;
}

template <typename T>
const T* DataAs(const NpyArray& arr) {
    return reinterpret_cast<const T*>(arr.bytes.data());
}

template <typename T>
std::vector<T> ReadH5Flat(const fs::path& path,
                          const std::string& dataset,
                          std::vector<std::size_t>* dimensions = nullptr) {
    HighFive::File file(path.string(), HighFive::File::ReadOnly);
    auto data_set = file.getDataSet(dataset);
    const auto dims = data_set.getSpace().getDimensions();
    if (dimensions) *dimensions = dims;
    const std::size_t count = std::accumulate(
        dims.begin(), dims.end(), std::size_t{1},
        std::multiplies<std::size_t>());
    std::vector<T> values(count);
    if (!values.empty()) data_set.read(values.data());
    return values;
}

bool RowHasFinite(const double* data, size_t row, size_t cols) {
    for (size_t c = 0; c < cols; ++c) {
        if (std::isfinite(data[row * cols + c]))
            return true;
    }
    return false;
}

SphericalTensor UnpackFull9(const double* data) {
    SphericalTensor tensor;
    tensor.T0 = data[0];
    for (size_t component = 0; component < 3; ++component)
        tensor.T1[component] = data[1 + component];
    for (size_t component = 0; component < 5; ++component)
        tensor.T2[component] = data[4 + component];
    return tensor;
}

void ExpectRawFull9TrajectoryMetadata(
    const fs::path& path, const std::string& group_path,
    const std::string& expected_transformation) {
    HighFive::File file(path.string(), HighFive::File::ReadOnly);
    auto group = file.getGroup(group_path);
    std::string parity;
    std::string coordinate_frame;
    std::string tensor_basis;
    std::string tensor_component_order;
    std::string tensor_frame;
    std::string tensor_t1_semantics;
    std::string tensor_structural_zero_components;
    std::string e3nn_export;
    std::string normalization_scope;
    std::string transformation;
    bool tensor_t1_structural_zero = true;
    group.getAttribute("parity").read(parity);
    group.getAttribute("coordinate_frame").read(coordinate_frame);
    group.getAttribute("tensor_basis").read(tensor_basis);
    group.getAttribute("tensor_component_order").read(
        tensor_component_order);
    group.getAttribute("tensor_frame").read(tensor_frame);
    group.getAttribute("tensor_t1_semantics").read(tensor_t1_semantics);
    group.getAttribute("tensor_t1_structural_zero").read(
        tensor_t1_structural_zero);
    group.getAttribute("tensor_structural_zero_components").read(
        tensor_structural_zero_components);
    group.getAttribute("e3nn_export").read(e3nn_export);
    group.getAttribute("normalization_scope").read(normalization_scope);
    group.getAttribute("transformation").read(transformation);

    EXPECT_EQ(parity, "mixed") << group_path;
    EXPECT_EQ(coordinate_frame, "conformation_cartesian_xyz")
        << group_path;
    EXPECT_EQ(tensor_basis,
              "project_native_full9_spherical_tensor_v1")
        << group_path;
    EXPECT_EQ(tensor_component_order,
              "T0,T1_x,T1_y,T1_z,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2")
        << group_path;
    EXPECT_EQ(tensor_frame, "conformation_cartesian_xyz") << group_path;
    EXPECT_EQ(tensor_t1_semantics,
              "Cartesian Levi-Civita dual x,y,z (not real-Y1m): "
              "a=((T_yz-T_zy)/2,(T_zx-T_xz)/2,(T_xy-T_yx)/2); "
              "axial a'=det(R) R a; generically nonzero")
        << group_path;
    EXPECT_FALSE(tensor_t1_structural_zero) << group_path;
    EXPECT_EQ(tensor_structural_zero_components, "none") << group_path;
    EXPECT_EQ(e3nn_export,
              "explicit project-basis to e3nn conversion required before use")
        << group_path;
    EXPECT_EQ(normalization_scope,
              "xyz tensor payload: T2 uses isometric real-tesseral "
              "normalization; T1 uses the tensor_t1_semantics convention")
        << group_path;
    EXPECT_EQ(transformation, expected_transformation) << group_path;
}

void ExpectSameOrNaN(double actual, double expected, double tolerance,
                     const std::string& context) {
    if (std::isnan(expected)) {
        EXPECT_TRUE(std::isnan(actual)) << context;
    } else {
        EXPECT_TRUE(std::isfinite(actual)) << context;
        EXPECT_NEAR(actual, expected, tolerance) << context;
    }
}

}  // namespace

// Forward-declare the file-local PRODUCTION dihedral helper (named per-file
// namespace, external linkage) so the fixed-coordinate ±60° test below pins
// production DIRECTLY, not a copy of it (vet finding 2026-07). A production
// sign regression now fails this non-skippable test.
namespace nmr {
namespace tripeptide_neighbor_dihedral {
double DihedralDegrees(const Vec3&, const Vec3&, const Vec3&, const Vec3&);
}
}  // namespace nmr

namespace nmr {
namespace tripeptide_neighbor_diagnostics {
std::array<double, 59> PackDiagnosticRow(const TripeptideNeighborShieldingResult::ResidueMatch&, std::uint8_t);
}
}  // namespace nmr

namespace nmr {
namespace tripeptide_neighbor_status {
int HisVariantDiagnosticCode(AminoAcid, int);
TripeptideMatchStatus UnsupportedHisStatus(AminoAcid, int);
}  // namespace tripeptide_neighbor_status
}  // namespace nmr


TEST(TripeptideNeighborDihedral, CanonicalFixedCoordinatesPinSignAndDegenerateNaN) {
    // Calls the PRODUCTION helper directly (not a copy), so a production sign
    // regression fails this non-skippable test (vet finding 2026-07).
    using nmr::tripeptide_neighbor_dihedral::DihedralDegrees;
    const Vec3 a(1.0, 0.0, 0.0);
    const Vec3 b(0.0, 0.0, 0.0);
    const Vec3 c(0.0, 0.0, 1.0);
    const double root3_over_2 = std::sqrt(3.0) / 2.0;

    // Analytic fixtures distinguish the canonical trajectory convention
    // from the former sign-reversed tripeptide triple product.
    EXPECT_NEAR(DihedralDegrees(a, b, c, Vec3(0.5, -root3_over_2, 1.0)), 60.0, 1e-12);
    EXPECT_NEAR(DihedralDegrees(a, b, c, Vec3(0.5, root3_over_2, 1.0)), -60.0, 1e-12);

    EXPECT_TRUE(std::isnan(DihedralDegrees(a, b, b, Vec3(0.0, 1.0, 0.0))));
    EXPECT_TRUE(
        std::isnan(DihedralDegrees(Vec3(0.0, 0.0, 0.0), Vec3(1.0, 0.0, 0.0), Vec3(2.0, 0.0, 0.0), Vec3(2.0, 1.0, 0.0))));
    EXPECT_TRUE(
        std::isnan(DihedralDegrees(Vec3(0.0, 1.0, 0.0), Vec3(0.0, 0.0, 0.0), Vec3(1.0, 0.0, 0.0), Vec3(2.0, 0.0, 0.0))));
}


TEST(TripeptideNeighborDiagnostics, ProductionPackerPinsExact59ColumnsStatusesAndGlobals) {
    TripeptideNeighborShieldingResult::ResidueMatch rm;

    auto fill = [](TripeptideNeighborShieldingResult::SideMatch& sm, double offset) {
        sm.calc_id = static_cast<int>(100 + offset);
        sm.neighbor_residue_index = static_cast<int>(10 + offset);
        sm.neighbor_residue_type_code = static_cast<int>(AminoAcid::LYS);
        sm.backbone_rmsd = 0.1 + offset;
        sm.n_atoms_matched = static_cast<int>(7 + offset);
        sm.natural_chi_axes = 4;
        sm.n_chi_axes_used = 1;
        sm.dropped_chi_distance_deg = 30.0 + offset;
        sm.phi_actual = -61.0 + offset;
        sm.psi_actual = 139.0 + offset;
        sm.phi_db = -60;
        sm.psi_db = 140;
        for (int k = 0; k < 4; ++k) {
            sm.chi_actual[k] = offset + 0.5 + k;
            sm.target_chi_grid[k] = 20 * (k + 1);
            sm.chi_db[k] = -20 * (k + 1);
        }
    };

    fill(rm.prev, 0.0);
    rm.prev.status = TripeptideMatchStatus::UnsupportedHid;
    rm.prev.frame_type = "gaussian_standard_orientation";
    rm.prev.his_variant_hint = 2;

    fill(rm.next, 1.0);
    rm.next.status = TripeptideMatchStatus::Ok;
    rm.next.frame_type = "orca_input_orientation";
    rm.next.his_variant_hint = 1;

    const auto row = nmr::tripeptide_neighbor_diagnostics::PackDiagnosticRow(rm, /*AAA OPBE=*/1);

    // Previous side, columns 0..27.
    EXPECT_DOUBLE_EQ(row[0], 2.0);
    EXPECT_DOUBLE_EQ(row[1], 100.0);
    EXPECT_DOUBLE_EQ(row[2], 0.0);
    EXPECT_DOUBLE_EQ(row[3], 10.0);
    EXPECT_DOUBLE_EQ(row[4], static_cast<double>(AminoAcid::LYS));
    EXPECT_DOUBLE_EQ(row[5], 0.1);
    EXPECT_DOUBLE_EQ(row[6], 1.0);
    EXPECT_DOUBLE_EQ(row[7], 7.0);
    EXPECT_DOUBLE_EQ(row[8], 4.0);
    EXPECT_DOUBLE_EQ(row[9], 1.0);
    EXPECT_DOUBLE_EQ(row[10], 30.0);
    EXPECT_DOUBLE_EQ(row[11], -61.0);
    EXPECT_DOUBLE_EQ(row[12], 139.0);
    EXPECT_DOUBLE_EQ(row[13], -60.0);
    EXPECT_DOUBLE_EQ(row[14], 140.0);
    for (int k = 0; k < 4; ++k) {
        EXPECT_DOUBLE_EQ(row[15 + k], 0.5 + k);
        EXPECT_DOUBLE_EQ(row[19 + k], 20.0 * (k + 1));
        EXPECT_DOUBLE_EQ(row[23 + k], -20.0 * (k + 1));
    }
    EXPECT_DOUBLE_EQ(row[27], 2.0);

    // Next side is the same exact layout at offset 28.
    EXPECT_DOUBLE_EQ(row[28 + 0], 1.0);
    EXPECT_DOUBLE_EQ(row[28 + 1], 101.0);
    EXPECT_DOUBLE_EQ(row[28 + 2], 1.0);
    EXPECT_DOUBLE_EQ(row[28 + 3], 11.0);
    EXPECT_DOUBLE_EQ(row[28 + 5], 1.1);
    EXPECT_DOUBLE_EQ(row[28 + 6], 2.0);
    EXPECT_DOUBLE_EQ(row[28 + 7], 8.0);
    EXPECT_DOUBLE_EQ(row[28 + 27], 1.0);

    // AAA method, mixed-method, any-match globals.
    EXPECT_DOUBLE_EQ(row[56], 1.0);
    EXPECT_DOUBLE_EQ(row[57], 1.0);
    EXPECT_DOUBLE_EQ(row[58], 1.0);

    // Exact frozen status integer values on both per-side blocks.
    const std::array<TripeptideMatchStatus, 5> statuses = {
        TripeptideMatchStatus::Miss,
        TripeptideMatchStatus::Ok,
        TripeptideMatchStatus::UnsupportedHid,
        TripeptideMatchStatus::UnsupportedHie,
        TripeptideMatchStatus::PerceptionFailed,
    };
    for (std::size_t i = 0; i < statuses.size(); ++i) {
        rm.prev.status = statuses[i];
        rm.next.status = statuses[statuses.size() - 1 - i];
        const auto status_row = nmr::tripeptide_neighbor_diagnostics::PackDiagnosticRow(rm, 1);
        EXPECT_DOUBLE_EQ(status_row[0], static_cast<double>(i));
        EXPECT_DOUBLE_EQ(status_row[2], i == 1 ? 1.0 : 0.0);
        EXPECT_DOUBLE_EQ(status_row[28], static_cast<double>(statuses.size() - 1 - i));
        EXPECT_DOUBLE_EQ(status_row[30], statuses.size() - 1 - i == 1 ? 1.0 : 0.0);
    }

    rm.prev.status = TripeptideMatchStatus::Ok;
    rm.prev.natural_chi_axes = 1;
    const auto shallow = nmr::tripeptide_neighbor_diagnostics::PackDiagnosticRow(rm, 1);
    EXPECT_TRUE(std::isfinite(shallow[19]));
    EXPECT_TRUE(std::isfinite(shallow[23]));
    for (int k = 1; k < 4; ++k) {
        EXPECT_TRUE(std::isnan(shallow[19 + k]));
        EXPECT_TRUE(std::isnan(shallow[23 + k]));
    }

    EXPECT_EQ(nmr::tripeptide_neighbor_status::UnsupportedHisStatus(AminoAcid::HIS, 0),
              TripeptideMatchStatus::UnsupportedHid);
    EXPECT_EQ(nmr::tripeptide_neighbor_status::UnsupportedHisStatus(AminoAcid::HIS, 1),
              TripeptideMatchStatus::UnsupportedHie);
    EXPECT_EQ(nmr::tripeptide_neighbor_status::UnsupportedHisStatus(AminoAcid::HIS, 2),
              TripeptideMatchStatus::Miss);
    EXPECT_EQ(nmr::tripeptide_neighbor_status::HisVariantDiagnosticCode(AminoAcid::HIS, 0), 2);
    EXPECT_EQ(nmr::tripeptide_neighbor_status::HisVariantDiagnosticCode(AminoAcid::HIS, 1), 3);
    EXPECT_EQ(nmr::tripeptide_neighbor_status::HisVariantDiagnosticCode(AminoAcid::HIS, 2), 1);
}


TEST(TripeptideNeighborM14, HidHieDirectionsAreCensoredWhileHipDirectionsCanMatch) {
    if (RuntimeEnvironment::TensorCs15Dsn().empty()) {
        GTEST_SKIP() << "tensorcs15 DSN not configured";
    }

    Session session;
    ASSERT_EQ(session.LoadTripeptideDftTable(), kOk) << session.LastError();
    const std::string pdb = nmr::test::TestEnvironment::UbqProtonated();
    auto built = BuildFromProtonatedPdb(pdb);
    ASSERT_TRUE(built.Ok()) << built.error;
    std::unique_ptr<Protein> protein = std::move(built.protein);
    auto& conf = protein->Conformation();
    ASSERT_TRUE(conf.AttachResult(GeometryResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(SpatialIndexResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(EnrichmentResult::Compute(conf)));

    std::size_t his_index = SIZE_MAX;
    for (std::size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        if (protein->ResidueAt(ri).type == AminoAcid::HIS) {
            his_index = ri;
            break;
        }
    }
    ASSERT_NE(his_index, SIZE_MAX);
    const auto before_his = protein->BackbonePredecessor(his_index);
    const auto after_his = protein->BackboneSuccessor(his_index);
    ASSERT_TRUE(before_his.has_value());
    ASSERT_TRUE(after_his.has_value());

    auto compute_variant = [&](int variant) {
        protein->MutableResidueAt(his_index).protonation_variant_index = variant;
        auto result = TripeptideNeighborShieldingResult::Compute(conf, *session.TripeptideDftTablePtr());
        EXPECT_NE(result, nullptr);
        return result;
    };

    auto hid = compute_variant(0);
    ASSERT_NE(hid, nullptr);
    EXPECT_EQ(hid->ResidueMatches()[*before_his].next.status, TripeptideMatchStatus::UnsupportedHid);
    EXPECT_EQ(hid->ResidueMatches()[*after_his].prev.status, TripeptideMatchStatus::UnsupportedHid);

    auto hie = compute_variant(1);
    ASSERT_NE(hie, nullptr);
    EXPECT_EQ(hie->ResidueMatches()[*before_his].next.status, TripeptideMatchStatus::UnsupportedHie);
    EXPECT_EQ(hie->ResidueMatches()[*after_his].prev.status, TripeptideMatchStatus::UnsupportedHie);

    auto hip = compute_variant(2);
    ASSERT_NE(hip, nullptr);
    EXPECT_EQ(hip->ResidueMatches()[*before_his].next.status, TripeptideMatchStatus::Ok);
    EXPECT_EQ(hip->ResidueMatches()[*after_his].prev.status, TripeptideMatchStatus::Ok);
    EXPECT_GT(hip->ResidueMatches()[*before_his].next.n_atoms_matched, 0);
    EXPECT_GT(hip->ResidueMatches()[*after_his].prev.n_atoms_matched, 0);
}


TEST(TripeptideNeighborDirectionalCovariance,
     ProductionLookupRerunAndAllSerializedOutputsUnderProperRigidTransform) {
    if (RuntimeEnvironment::TensorCs15Dsn().empty()) {
        GTEST_SKIP() << "tensorcs15 DSN not configured";
    }

    Session session;
    ASSERT_EQ(session.LoadTripeptideDftTable(), kOk) << session.LastError();
    auto built = BuildFromProtonatedPdb(
        nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(built.Ok()) << built.error;
    auto protein = std::move(built.protein);

    ProteinConformation& original = protein->Conformation();
    ASSERT_TRUE(original.AttachResult(GeometryResult::Compute(original)));
    ASSERT_TRUE(original.AttachResult(EnrichmentResult::Compute(original)));
    auto original_result = TripeptideNeighborShieldingResult::Compute(
        original, *session.TripeptideDftTablePtr());
    ASSERT_NE(original_result, nullptr);
    ASSERT_GT(original_result->AtomsAccumulated(), 0);

    // tensorcs15 stores chiral tripeptide poses rather than a mirrored
    // companion database.  The source-preserving contract is therefore
    // SO(3).  Seed and tolerances are frozen here for the pre-freeze gate.
    constexpr std::uint64_t kTransformSeed = 0x5452494E45494748ULL;
    const auto transform = nmr::test::directional::SeededTransform(
        kTransformSeed, false);
    ProteinConformation& moved = protein->AddConformation(
        nmr::test::directional::Positions(transform,
                                          original.Positions()),
        "tripeptide neighbor proper covariance rerun");
    ASSERT_TRUE(moved.AttachResult(GeometryResult::Compute(moved)));
    ASSERT_TRUE(moved.AttachResult(EnrichmentResult::Compute(moved)));
    auto moved_result = TripeptideNeighborShieldingResult::Compute(
        moved, *session.TripeptideDftTablePtr());
    ASSERT_NE(moved_result, nullptr);
    ASSERT_EQ(moved_result->ResiduesWithAnyNeighbor(),
              original_result->ResiduesWithAnyNeighbor());
    ASSERT_EQ(moved_result->AtomsAccumulated(),
              original_result->AtomsAccumulated());

    constexpr double kTensorAbsTolerance = 8.0e-8;
    constexpr double kTensorRelTolerance = 2.0e-10;
    constexpr double kVectorAbsTolerance = 8.0e-8;
    constexpr double kVectorRelTolerance = 2.0e-10;
    constexpr double kScalarTolerance = 8.0e-8;

    const fs::path original_dir = nmr::test::TestEnvironment::TempPath(
        "tripeptide_neighbor_directional_original");
    const fs::path moved_dir = nmr::test::TestEnvironment::TempPath(
        "tripeptide_neighbor_directional_moved");
    fs::create_directories(original_dir);
    fs::create_directories(moved_dir);
    ASSERT_EQ(original_result->WriteFeatures(original,
                                             original_dir.string()), 7);
    ASSERT_EQ(moved_result->WriteFeatures(moved, moved_dir.string()), 7);
    auto all_neighbor_feature_paths = [](const fs::path& dir) {
        return std::vector<fs::path>{
            dir / "tripeptide_neighbor_shielding.npy",
            dir / "tripeptide_neighbor_shielding_prev.npy",
            dir / "tripeptide_neighbor_shielding_next.npy",
            dir / "tripeptide_neighbor_residual_vec_prev.npy",
            dir / "tripeptide_neighbor_residual_vec_next.npy",
            dir / "tripeptide_neighbor_reference.npy",
            dir / "tripeptide_neighbor_diagnostics.npy",
        };
    };
    ASSERT_EQ(nmr::test::directional::RunNumpyAllowPickleFalse(
                  NMR_TEST_PYTHON_EXECUTABLE,
                  NMR_NPY_ALLOW_PICKLE_FALSE_SCRIPT,
                  all_neighbor_feature_paths(original_dir)),
              0);
    ASSERT_EQ(nmr::test::directional::RunNumpyAllowPickleFalse(
                  NMR_TEST_PYTHON_EXECUTABLE,
                  NMR_NPY_ALLOW_PICKLE_FALSE_SCRIPT,
                  all_neighbor_feature_paths(moved_dir)),
              0);

    // These comparisons start from bytes reread from both production NPY
    // writers.  No already-produced output is transformed and fed back to a
    // calculator: the expected law is applied only after both real Compute()
    // reruns have completed.
    int finite_tensor_rows = 0;
    for (const char* name : {
             "tripeptide_neighbor_shielding.npy",
             "tripeptide_neighbor_shielding_prev.npy",
             "tripeptide_neighbor_shielding_next.npy"}) {
        const NpyArray source = ReadNpy(original_dir / name);
        const NpyArray actual = ReadNpy(moved_dir / name);
        ASSERT_EQ(source.shape,
                  (std::vector<size_t>{original.AtomCount(), 9}))
            << name;
        ASSERT_EQ(actual.shape, source.shape) << name;
        ASSERT_EQ(source.bytes.size(), original.AtomCount() * 9 * sizeof(double));
        ASSERT_EQ(actual.bytes.size(), source.bytes.size());
        const double* source_data = DataAs<double>(source);
        const double* actual_data = DataAs<double>(actual);
        for (size_t atom_index = 0; atom_index < original.AtomCount();
             ++atom_index) {
            bool row_finite = false;
            for (size_t component = 0; component < 9; ++component) {
                const double source_value =
                    source_data[atom_index * 9 + component];
                const double actual_value =
                    actual_data[atom_index * 9 + component];
                EXPECT_EQ(std::isnan(actual_value),
                          std::isnan(source_value))
                    << name << " atom=" << atom_index
                    << " component=" << component;
                row_finite = row_finite || std::isfinite(source_value);
            }
            if (!row_finite)
                continue;
            ++finite_tensor_rows;
            const SphericalTensor source_tensor =
                UnpackFull9(source_data + atom_index * 9);
            const SphericalTensor actual_tensor =
                UnpackFull9(actual_data + atom_index * 9);
            const Mat3 expected_matrix = nmr::test::directional::EvenRank2(
                transform, source_tensor.Reconstruct());
            EXPECT_TRUE(nmr::test::directional::NearMatrix(
                actual_tensor.Reconstruct(), expected_matrix,
                kTensorAbsTolerance, kTensorRelTolerance))
                << name << " atom=" << atom_index;

            // Explicit native-T2 round trip: serialized native five-vector
            // -> Cartesian symmetric/traceless -> Q T Q^T -> native
            // five-vector, compared with the rerun's serialized five-vector.
            const SphericalTensor expected_t2 =
                nmr::test::directional::RotateNativeT2(
                    transform, source_tensor);
            for (size_t component = 0; component < 5; ++component) {
                EXPECT_NEAR(actual_tensor.T2[component],
                            expected_t2.T2[component],
                            kTensorAbsTolerance)
                    << name << " atom=" << atom_index
                    << " native_T2_component=" << component;
            }
            EXPECT_NEAR(actual_tensor.T0, source_tensor.T0,
                        kScalarTolerance)
                << name << " atom=" << atom_index << " T0";
            const Vec3 expected_t1 = nmr::test::directional::Axial(
                transform,
                nmr::test::directional::T1Vector(source_tensor));
            EXPECT_TRUE(nmr::test::directional::NearVector(
                nmr::test::directional::T1Vector(actual_tensor),
                expected_t1, kTensorAbsTolerance,
                kTensorRelTolerance))
                << name << " atom=" << atom_index << " T1";
        }
    }
    EXPECT_GT(finite_tensor_rows, 300);

    int finite_vector_rows = 0;
    for (const char* name : {
             "tripeptide_neighbor_residual_vec_prev.npy",
             "tripeptide_neighbor_residual_vec_next.npy"}) {
        const NpyArray source = ReadNpy(original_dir / name);
        const NpyArray actual = ReadNpy(moved_dir / name);
        ASSERT_EQ(source.shape,
                  (std::vector<size_t>{original.AtomCount(), 3}))
            << name;
        ASSERT_EQ(actual.shape, source.shape) << name;
        ASSERT_EQ(source.bytes.size(), original.AtomCount() * 3 * sizeof(double));
        ASSERT_EQ(actual.bytes.size(), source.bytes.size());
        const double* source_data = DataAs<double>(source);
        const double* actual_data = DataAs<double>(actual);
        for (size_t atom_index = 0; atom_index < original.AtomCount();
             ++atom_index) {
            bool row_finite = false;
            for (size_t component = 0; component < 3; ++component) {
                const double source_value =
                    source_data[atom_index * 3 + component];
                const double actual_value =
                    actual_data[atom_index * 3 + component];
                EXPECT_EQ(std::isnan(actual_value),
                          std::isnan(source_value))
                    << name << " atom=" << atom_index
                    << " component=" << component;
                row_finite = row_finite || std::isfinite(source_value);
            }
            if (!row_finite)
                continue;
            ++finite_vector_rows;
            const Vec3 source_vector(
                source_data[atom_index * 3 + 0],
                source_data[atom_index * 3 + 1],
                source_data[atom_index * 3 + 2]);
            const Vec3 actual_vector(
                actual_data[atom_index * 3 + 0],
                actual_data[atom_index * 3 + 1],
                actual_data[atom_index * 3 + 2]);
            EXPECT_TRUE(nmr::test::directional::NearVector(
                actual_vector,
                nmr::test::directional::Polar(transform, source_vector),
                kVectorAbsTolerance, kVectorRelTolerance))
                << name << " atom=" << atom_index;
        }
    }
    EXPECT_GT(finite_vector_rows, 200);

    for (const char* name : {
             "tripeptide_neighbor_reference.npy",
             "tripeptide_neighbor_diagnostics.npy"}) {
        const NpyArray source = ReadNpy(original_dir / name);
        const NpyArray actual = ReadNpy(moved_dir / name);
        ASSERT_EQ(actual.shape, source.shape) << name;
        ASSERT_EQ(actual.bytes.size(), source.bytes.size()) << name;
        ASSERT_EQ(source.bytes.size() % sizeof(double), 0U) << name;
        const size_t count = source.bytes.size() / sizeof(double);
        const double* source_data = DataAs<double>(source);
        const double* actual_data = DataAs<double>(actual);
        for (size_t component = 0; component < count; ++component) {
            ExpectSameOrNaN(actual_data[component], source_data[component],
                            kScalarTolerance,
                            std::string(name) + " component=" +
                                std::to_string(component));
        }
    }

    // Exercise the exact trajectory owners with two real, attached owner
    // results: source frame then the independently rerun proper transform.
    ASSERT_TRUE(original.AttachResult(std::move(original_result)));
    ASSERT_TRUE(moved.AttachResult(std::move(moved_result)));
    auto trajectory_protein =
        TrajectoryProtein::CreateForTesting(std::move(protein));
    ASSERT_NE(trajectory_protein, nullptr);
    auto shielding_ts =
        TripeptideNeighborShieldingTimeSeriesTrajectoryResult::Create(
            *trajectory_protein);
    auto prev_ts =
        TripeptideNeighborResidualVecPrevTimeSeriesTrajectoryResult::Create(
            *trajectory_protein);
    auto next_ts =
        TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult::Create(
            *trajectory_protein);
    ASSERT_NE(shielding_ts, nullptr);
    ASSERT_NE(prev_ts, nullptr);
    ASSERT_NE(next_ts, nullptr);
    Trajectory dummy("", "", "");
    for (auto* result : std::array<TrajectoryResult*, 3>{
             shielding_ts.get(), prev_ts.get(), next_ts.get()}) {
        result->Compute(original, *trajectory_protein, dummy, 43, 3.5);
        result->Compute(moved, *trajectory_protein, dummy, 47, 4.5);
        result->Finalize(*trajectory_protein, dummy);
    }
    const fs::path h5_path = nmr::test::TestEnvironment::TempPath(
        "tripeptide_neighbor_directional_raw.h5");
    (void)std::remove(h5_path.string().c_str());
    {
        HighFive::File file(h5_path.string(), HighFive::File::Truncate);
        shielding_ts->WriteH5Group(*trajectory_protein, file);
        prev_ts->WriteH5Group(*trajectory_protein, file);
        next_ts->WriteH5Group(*trajectory_protein, file);
    }

    const std::string shielding_path =
        "/trajectory/tripeptide_neighbor_shielding_time_series/xyz";
    const std::string prev_path =
        "/trajectory/tripeptide_neighbor_residual_vec_prev_time_series/xyz";
    const std::string next_path =
        "/trajectory/tripeptide_neighbor_residual_vec_next_time_series/xyz";
    std::vector<std::size_t> shielding_dims;
    std::vector<std::size_t> prev_dims;
    std::vector<std::size_t> next_dims;
    const auto h5_shielding = ReadH5Flat<double>(
        h5_path, shielding_path, &shielding_dims);
    const auto h5_prev = ReadH5Flat<double>(
        h5_path, prev_path, &prev_dims);
    const auto h5_next = ReadH5Flat<double>(
        h5_path, next_path, &next_dims);
    const std::size_t atom_count = original.AtomCount();
    EXPECT_EQ(shielding_dims,
              (std::vector<std::size_t>{atom_count, 2u, 9u}));
    EXPECT_EQ(prev_dims,
              (std::vector<std::size_t>{atom_count, 2u, 3u}));
    EXPECT_EQ(next_dims, prev_dims);
    ASSERT_EQ(h5_shielding.size(), atom_count * 2 * 9);
    ASSERT_EQ(h5_prev.size(), atom_count * 2 * 3);
    ASSERT_EQ(h5_next.size(), atom_count * 2 * 3);

    // This one group describes both serialized frames: the source rerun and
    // the independently transformed production rerun.
    ExpectRawFull9TrajectoryMetadata(
        h5_path,
        "/trajectory/tripeptide_neighbor_shielding_time_series",
        "even_rank2 under proper rotations: T'=R T R^T; typed tripeptide "
        "lookup/Kabsch alignment is L-amino-acid chirality-conditioned and "
        "has no improper-transform contract");

    for (const std::string& group : {
             std::string(
                 "/trajectory/tripeptide_neighbor_shielding_time_series"),
             std::string(
                 "/trajectory/tripeptide_neighbor_residual_vec_prev_time_series"),
             std::string(
                 "/trajectory/tripeptide_neighbor_residual_vec_next_time_series")}) {
        EXPECT_EQ(ReadH5Flat<std::size_t>(
                      h5_path, group + "/frame_indices"),
                  (std::vector<std::size_t>{43u, 47u})) << group;
        EXPECT_EQ(ReadH5Flat<double>(h5_path, group + "/frame_times"),
                  (std::vector<double>{3.5, 4.5})) << group;
        EXPECT_EQ(ReadH5Flat<std::uint8_t>(
                      h5_path, group + "/source_attached_per_frame"),
                  (std::vector<std::uint8_t>{1u, 1u})) << group;
    }

    int finite_h5_tensor_rows = 0;
    int finite_h5_prev_rows = 0;
    int finite_h5_next_rows = 0;
    for (std::size_t atom_index = 0; atom_index < atom_count;
         ++atom_index) {
        const std::size_t tensor_source_base = (atom_index * 2) * 9;
        const std::size_t tensor_moved_base =
            (atom_index * 2 + 1) * 9;
        const bool matched = original.AtomAt(atom_index)
                                 .tripeptide_neighbor_has_match;
        EXPECT_EQ(moved.AtomAt(atom_index).tripeptide_neighbor_has_match,
                  matched) << atom_index;
        for (std::size_t component = 0; component < 9; ++component) {
            EXPECT_EQ(std::isnan(h5_shielding[
                          tensor_moved_base + component]),
                      std::isnan(h5_shielding[
                          tensor_source_base + component]))
                << shielding_path << " atom=" << atom_index
                << " component=" << component;
            EXPECT_EQ(std::isnan(h5_shielding[
                          tensor_source_base + component]),
                      !matched)
                << shielding_path << " atom=" << atom_index
                << " component=" << component << " applicability";
        }
        if (matched) {
            ++finite_h5_tensor_rows;
            const SphericalTensor source_tensor =
                UnpackFull9(h5_shielding.data() + tensor_source_base);
            const SphericalTensor actual_tensor =
                UnpackFull9(h5_shielding.data() + tensor_moved_base);
            EXPECT_TRUE(nmr::test::directional::NearMatrix(
                actual_tensor.Reconstruct(),
                nmr::test::directional::EvenRank2(
                    transform, source_tensor.Reconstruct()),
                kTensorAbsTolerance, kTensorRelTolerance))
                << shielding_path << " atom=" << atom_index;
            EXPECT_NEAR(actual_tensor.T0, source_tensor.T0,
                        kTensorAbsTolerance)
                << shielding_path << " atom=" << atom_index << " T0";
            EXPECT_TRUE(nmr::test::directional::NearVector(
                nmr::test::directional::T1Vector(actual_tensor),
                nmr::test::directional::Axial(
                    transform,
                    nmr::test::directional::T1Vector(source_tensor)),
                kTensorAbsTolerance, kTensorRelTolerance))
                << shielding_path << " atom=" << atom_index << " T1";
            const SphericalTensor expected_t2 =
                nmr::test::directional::RotateNativeT2(
                    transform, source_tensor);
            for (std::size_t component = 0; component < 5; ++component) {
                EXPECT_NEAR(actual_tensor.T2[component],
                            expected_t2.T2[component],
                            kTensorAbsTolerance)
                    << shielding_path << " atom=" << atom_index
                    << " native_T2_component=" << component;
            }
        }

        auto check_residual = [&](const std::vector<double>& values,
                                  const std::string& path,
                                  int& finite_rows) {
            const std::size_t source_base = (atom_index * 2) * 3;
            const std::size_t moved_base = (atom_index * 2 + 1) * 3;
            bool finite = false;
            for (std::size_t component = 0; component < 3; ++component) {
                EXPECT_EQ(std::isnan(values[moved_base + component]),
                          std::isnan(values[source_base + component]))
                    << path << " atom=" << atom_index
                    << " component=" << component;
                finite = finite ||
                    std::isfinite(values[source_base + component]);
            }
            if (!finite) {
                for (std::size_t component = 0; component < 3;
                     ++component) {
                    EXPECT_TRUE(std::isnan(values[source_base + component]))
                        << path << " atom=" << atom_index
                        << " component=" << component;
                    EXPECT_TRUE(std::isnan(values[moved_base + component]))
                        << path << " atom=" << atom_index
                        << " moved component=" << component;
                }
                return;
            }
            ++finite_rows;
            for (std::size_t component = 0; component < 3; ++component) {
                ASSERT_TRUE(std::isfinite(values[source_base + component]))
                    << path << " atom=" << atom_index;
                ASSERT_TRUE(std::isfinite(values[moved_base + component]))
                    << path << " atom=" << atom_index;
            }
            const Vec3 source(values[source_base],
                              values[source_base + 1],
                              values[source_base + 2]);
            const Vec3 actual(values[moved_base],
                              values[moved_base + 1],
                              values[moved_base + 2]);
            EXPECT_TRUE(nmr::test::directional::NearVector(
                actual, nmr::test::directional::Polar(transform, source),
                kVectorAbsTolerance, kVectorRelTolerance))
                << path << " atom=" << atom_index;
        };
        check_residual(h5_prev, prev_path, finite_h5_prev_rows);
        check_residual(h5_next, next_path, finite_h5_next_rows);
    }
    EXPECT_GT(finite_h5_tensor_rows, 0);
    EXPECT_GT(finite_h5_prev_rows, 0);
    EXPECT_GT(finite_h5_next_rows, 0);
    EXPECT_EQ(std::remove(h5_path.string().c_str()), 0) << h5_path;

    for (const fs::path& dir : {original_dir, moved_dir}) {
        for (const char* name : {
                 "tripeptide_neighbor_shielding.npy",
                 "tripeptide_neighbor_shielding_prev.npy",
                 "tripeptide_neighbor_shielding_next.npy",
                 "tripeptide_neighbor_residual_vec_prev.npy",
                 "tripeptide_neighbor_residual_vec_next.npy",
                 "tripeptide_neighbor_reference.npy",
                 "tripeptide_neighbor_diagnostics.npy"}) {
            EXPECT_EQ(std::remove((dir / name).string().c_str()), 0)
                << name;
        }
        EXPECT_EQ(std::remove(dir.string().c_str()), 0) << dir;
    }
}


class TripeptideNeighborShieldingTest : public ::testing::Test {
protected:
    void SetUp() override {
        const std::string kLarsen1UbqPm6 = nmr::test::TestEnvironment::Larsen1UbqPm6Pdb();
        if (RuntimeEnvironment::TensorCs15Dsn().empty()) {
            GTEST_SKIP() << "tensorcs15 DSN not configured";
        }
        if (!fs::exists(kLarsen1UbqPm6)) {
            GTEST_SKIP() << "Larsen 1UBQ PM6-D3H+ PDB not at " << kLarsen1UbqPm6;
        }
        ASSERT_EQ(session.LoadTripeptideDftTable(), kOk) << session.LastError();
        ASSERT_TRUE(session.HasTripeptideDftTable());

        auto r = BuildFromProtonatedPdb(kLarsen1UbqPm6);
        ASSERT_TRUE(r.Ok()) << r.error;
        protein = std::move(r.protein);
    }

    Session session;
    std::unique_ptr<Protein> protein;
};


TEST_F(TripeptideNeighborShieldingTest, RunsOn1UbqPm6) {
    auto& conf = protein->Conformation();

    auto geo = GeometryResult::Compute(conf);
    ASSERT_NE(geo, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(geo)));

    auto sp = SpatialIndexResult::Compute(conf);
    ASSERT_NE(sp, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(sp)));

    auto en = EnrichmentResult::Compute(conf);
    ASSERT_NE(en, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(en)));

    auto tn = TripeptideNeighborShieldingResult::Compute(conf, *session.TripeptideDftTablePtr());
    ASSERT_NE(tn, nullptr);

    EXPECT_GE(tn->ResiduesWithAnyNeighbor(), 60) << "expected most internal residues to get >=1 neighbor "
                                                    "contribution; got "
                                                 << tn->ResiduesWithAnyNeighbor();
    EXPECT_GT(tn->AtomsAccumulated(), 200) << "expected >200 per-atom Δσ accumulations; got " << tn->AtomsAccumulated();

    // Exercise the production file-local helper through the DB row it
    // selected. ALA has no chi axes, so independently computed canonical
    // phi/psi identify one unambiguous 2-degree-grid row. The former
    // sign-reversed helper selects the opposite-angle row and a different
    // calc_id.
    auto query_ala_from_geometry = [&](size_t neighbor_index) {
        TripeptideDftRecord miss;
        const auto prev_idx = protein->BackbonePredecessor(neighbor_index);
        const auto next_idx = protein->BackboneSuccessor(neighbor_index);
        if (!prev_idx || !next_idx)
            return miss;
        const Residue& prev = protein->ResidueAt(*prev_idx);
        const Residue& neighbor = protein->ResidueAt(neighbor_index);
        const Residue& next = protein->ResidueAt(*next_idx);
        if (prev.C == Residue::NONE || neighbor.N == Residue::NONE || neighbor.CA == Residue::NONE
            || neighbor.C == Residue::NONE || next.N == Residue::NONE) {
            return miss;
        }
        const double phi = nmr::tripeptide_neighbor_dihedral::DihedralDegrees(conf.PositionAt(prev.C),
                                                                              conf.PositionAt(neighbor.N),
                                                                              conf.PositionAt(neighbor.CA),
                                                                              conf.PositionAt(neighbor.C));
        const double psi = nmr::tripeptide_neighbor_dihedral::DihedralDegrees(conf.PositionAt(neighbor.N),
                                                                              conf.PositionAt(neighbor.CA),
                                                                              conf.PositionAt(neighbor.C),
                                                                              conf.PositionAt(next.N));
        return session.TripeptideDftTablePtr()->QueryNearest('A',
                                                             phi,
                                                             psi,
                                                             0.0,
                                                             0.0,
                                                             0.0,
                                                             0.0,
                                                             /*n_chi_axes=*/0,
                                                             /*his_variant_hint=*/-1);
    };

    int n_ala_dihedral_queries_pinned = 0;
    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        const auto& match = tn->ResidueMatches()[ri];
        if (const auto prev_idx = protein->BackbonePredecessor(ri);
            prev_idx && protein->ResidueAt(*prev_idx).type == AminoAcid::ALA && match.prev_calc_id != 0) {
            const auto expected = query_ala_from_geometry(*prev_idx);
            ASSERT_TRUE(expected.IsHit());
            ASSERT_TRUE(expected.larsen.has_value());
            EXPECT_EQ(match.prev_calc_id, expected.calc_id) << "i-1 ALA neighbor at residue " << ri;
            ++n_ala_dihedral_queries_pinned;
        }
        if (const auto next_idx = protein->BackboneSuccessor(ri);
            next_idx && protein->ResidueAt(*next_idx).type == AminoAcid::ALA && match.next_calc_id != 0) {
            const auto expected = query_ala_from_geometry(*next_idx);
            ASSERT_TRUE(expected.IsHit());
            ASSERT_TRUE(expected.larsen.has_value());
            EXPECT_EQ(match.next_calc_id, expected.calc_id) << "i+1 ALA neighbor at residue " << ri;
            ++n_ala_dihedral_queries_pinned;
        }
    }
    EXPECT_GT(n_ala_dihedral_queries_pinned, 0) << "1UBQ should exercise an ALA neighbor DB lookup";

    // Frame_type discriminator: SER residues as i-1 or i+1 neighbors
    // should produce ResidueMatch entries with frame_type ==
    // "orca_input_orientation" (the project SER PBE regen).
    int n_pbe_dir = 0;
    int n_opbe_dir = 0;
    for (const auto& m : tn->ResidueMatches()) {
        for (const auto& ft : {m.prev_frame_type, m.next_frame_type}) {
            if (ft == "orca_input_orientation")
                ++n_pbe_dir;
            else if (ft == "gaussian_standard_orientation")
                ++n_opbe_dir;
        }
    }
    EXPECT_GT(n_pbe_dir, 0) << "at least one SER-side neighbor lookup should hit ORCA PBE; " << "got " << n_pbe_dir
                            << " PBE vs " << n_opbe_dir << " OPBE";
    EXPECT_GT(n_opbe_dir, n_pbe_dir) << "OPBE should dominate (only SER is PBE); " << "got " << n_opbe_dir << " OPBE, "
                                     << n_pbe_dir << " PBE";

    // Residual sanity: per-direction residual_vecs populated on
    // matched atoms. The cap-side Kabsch aligns flanking-ALA N/CA/C
    // onto residue i's N/CA/C; the AAA reference is also at standard
    // angles, so AXA and AAA residuals should both be small.
    //
    // Per M4 contract, the per-direction residual is NaN whenever a
    // direction had no contribution at all (terminal residue's
    // missing i-1 / i+1, or chain break). The previous formulation
    // here read .norm() on every has_match atom and compared the
    // result against 3.0 — NaN comparisons evaluate false in C++ so
    // an entire direction being absent everywhere would silently
    // pass the extreme bound check. The fix is to count finite
    // residuals separately and only fold finite vectors into the
    // max/extreme tallies; the test then also asserts non-trivial
    // finite coverage in each direction.
    int n_atoms_with_neighbor = 0;
    int n_finite_prev = 0, n_finite_next = 0;
    int n_extreme_prev = 0, n_extreme_next = 0;
    double max_prev = 0.0, max_next = 0.0;
    for (size_t i = 0; i < conf.AtomCount(); ++i) {
        const auto& ca = conf.AtomAt(i);
        if (!ca.tripeptide_neighbor_has_match)
            continue;
        ++n_atoms_with_neighbor;
        const double pr = ca.tripeptide_neighbor_residual_vec_prev.norm();
        const double nr = ca.tripeptide_neighbor_residual_vec_next.norm();
        if (std::isfinite(pr)) {
            ++n_finite_prev;
            if (pr > max_prev)
                max_prev = pr;
            if (pr > 3.0)
                ++n_extreme_prev;
        }
        if (std::isfinite(nr)) {
            ++n_finite_next;
            if (nr > max_next)
                max_next = nr;
            if (nr > 3.0)
                ++n_extreme_next;
        }
    }
    EXPECT_GT(n_atoms_with_neighbor, 200);
    EXPECT_GT(n_finite_prev, 100) << "expected most has_match atoms to also have a finite "
                                     "i-1 residual; got "
                                  << n_finite_prev << " finite of " << n_atoms_with_neighbor << " has_match atoms";
    EXPECT_GT(n_finite_next, 100) << "expected most has_match atoms to also have a finite "
                                     "i+1 residual; got "
                                  << n_finite_next << " finite of " << n_atoms_with_neighbor << " has_match atoms";
    EXPECT_LT(n_extreme_prev, 20) << "more than 20 atoms have prev-direction residual > 3 Å " << "(max=" << max_prev
                                  << " Å) — likely Kabsch issue";
    EXPECT_LT(n_extreme_next, 20) << "more than 20 atoms have next-direction residual > 3 Å " << "(max=" << max_next
                                  << " Å) — likely Kabsch issue";

    // NPY emission.
    const std::string out_dir = nmr::test::TestEnvironment::TempPath("tripeptide_neighbor_smoke_out");
    fs::create_directories(out_dir);
    int n_npy = tn->WriteFeatures(conf, out_dir);
    EXPECT_EQ(n_npy, 7);
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_neighbor_shielding.npy"));
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_neighbor_shielding_prev.npy"));
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_neighbor_shielding_next.npy"));
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_neighbor_residual_vec_prev.npy"));
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_neighbor_residual_vec_next.npy"));
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_neighbor_reference.npy"));
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_neighbor_diagnostics.npy"));
    EXPECT_FALSE(fs::exists(out_dir + "/tripeptide_neighbor_prev_status.npy"));
    EXPECT_FALSE(fs::exists(out_dir + "/tripeptide_neighbor_next_status.npy"));

    auto sum_arr = ReadNpy(out_dir + "/tripeptide_neighbor_shielding.npy");
    auto prev_arr = ReadNpy(out_dir + "/tripeptide_neighbor_shielding_prev.npy");
    auto next_arr = ReadNpy(out_dir + "/tripeptide_neighbor_shielding_next.npy");
    auto ref_arr = ReadNpy(out_dir + "/tripeptide_neighbor_reference.npy");
    auto diag_arr = ReadNpy(out_dir + "/tripeptide_neighbor_diagnostics.npy");
    const std::vector<size_t> tensor_shape{conf.AtomCount(), 9};
    ASSERT_EQ(sum_arr.shape, tensor_shape);
    ASSERT_EQ(prev_arr.shape, tensor_shape);
    ASSERT_EQ(next_arr.shape, tensor_shape);
    ASSERT_EQ(ref_arr.shape, (std::vector<size_t>{1, 5}));
    ASSERT_EQ(diag_arr.shape, (std::vector<size_t>{protein->ResidueCount(), 59}));

    const double* sum = DataAs<double>(sum_arr);
    const double* prev = DataAs<double>(prev_arr);
    const double* next = DataAs<double>(next_arr);
    const double* ref = DataAs<double>(ref_arr);
    int n_prev_rows = 0, n_next_rows = 0, n_reconstructed = 0;
    for (size_t i = 0; i < conf.AtomCount(); ++i) {
        const bool has_prev = RowHasFinite(prev, i, 9);
        const bool has_next = RowHasFinite(next, i, 9);
        if (has_prev)
            ++n_prev_rows;
        if (has_next)
            ++n_next_rows;
        if (!has_prev && !has_next) {
            for (int c = 0; c < 9; ++c)
                EXPECT_TRUE(std::isnan(sum[i * 9 + c]));
            continue;
        }
        ++n_reconstructed;
        for (int c = 0; c < 9; ++c) {
            const double expected = (has_prev ? prev[i * 9 + c] : 0.0) + (has_next ? next[i * 9 + c] : 0.0);
            EXPECT_NEAR(sum[i * 9 + c], expected, 1e-9) << "atom=" << i << " component=" << c;
        }
    }
    EXPECT_GT(n_prev_rows, 100);
    EXPECT_GT(n_next_rows, 100);
    EXPECT_GT(n_reconstructed, 200);
    EXPECT_GT(ref[0], 0.0);
    EXPECT_DOUBLE_EQ(ref[1], 1.0);
    EXPECT_DOUBLE_EQ(ref[2], -120.0);
    EXPECT_DOUBLE_EQ(ref[3], 140.0);
    EXPECT_DOUBLE_EQ(ref[4], 1.0);

    const double* diag = DataAs<double>(diag_arr);
    int n_fallback_sides = 0;
    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        const auto& m = tn->ResidueMatches()[ri];
        const double* row = &diag[ri * 59];
        for (const auto side : {std::pair<const TripeptideNeighborShieldingResult::SideMatch*, size_t>{&m.prev, 0},
                                std::pair<const TripeptideNeighborShieldingResult::SideMatch*, size_t>{&m.next, 28}}) {
            const auto& sm = *side.first;
            const size_t b = side.second;
            EXPECT_DOUBLE_EQ(row[b + 0], static_cast<double>(sm.status));
            EXPECT_DOUBLE_EQ(row[b + 1], static_cast<double>(sm.calc_id));
            EXPECT_DOUBLE_EQ(row[b + 2], sm.status == TripeptideMatchStatus::Ok ? 1.0 : 0.0);
            EXPECT_DOUBLE_EQ(row[b + 3], static_cast<double>(sm.neighbor_residue_index));
            EXPECT_DOUBLE_EQ(row[b + 4], static_cast<double>(sm.neighbor_residue_type_code));
            EXPECT_DOUBLE_EQ(row[b + 7], static_cast<double>(sm.n_atoms_matched));
            EXPECT_DOUBLE_EQ(row[b + 8], static_cast<double>(sm.natural_chi_axes));
            EXPECT_DOUBLE_EQ(row[b + 9], static_cast<double>(sm.n_chi_axes_used));
            EXPECT_DOUBLE_EQ(row[b + 27], static_cast<double>(sm.his_variant_hint));
            if (sm.status == TripeptideMatchStatus::Ok && sm.n_chi_axes_used < sm.natural_chi_axes) {
                ++n_fallback_sides;
                EXPECT_TRUE(std::isfinite(row[b + 10]));
                EXPECT_GE(row[b + 10], 0.0);
                EXPECT_DOUBLE_EQ(row[b + 10], sm.dropped_chi_distance_deg);
            }
            if (sm.status == TripeptideMatchStatus::UnsupportedHid) {
                EXPECT_DOUBLE_EQ(row[b + 0], 2.0);
                EXPECT_DOUBLE_EQ(row[b + 2], 0.0);
                EXPECT_DOUBLE_EQ(row[b + 27], 2.0);
            } else if (sm.status == TripeptideMatchStatus::UnsupportedHie) {
                EXPECT_DOUBLE_EQ(row[b + 0], 3.0);
                EXPECT_DOUBLE_EQ(row[b + 2], 0.0);
                EXPECT_DOUBLE_EQ(row[b + 27], 3.0);
            }
        }
        EXPECT_DOUBLE_EQ(row[56], ref[1]);
        const bool any_match = m.prev.status == TripeptideMatchStatus::Ok || m.next.status == TripeptideMatchStatus::Ok;
        EXPECT_DOUBLE_EQ(row[58], any_match ? 1.0 : 0.0);
    }
    EXPECT_GT(n_fallback_sides, 0) << "1UBQ should exercise at least one neighbor dropped-chi fallback";
}

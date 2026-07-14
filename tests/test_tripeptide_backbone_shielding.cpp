// Smoke test for TripeptideBackboneShieldingResult.
//
// Loads 1UBQ_pm6dh3plus.pdb (Larsen's PM6-D3H+ optimised geometry —
// /mnt/expansion/larsen_archive/structures/, the published-RMSD
// validation target from Larsen 2015) via BuildFromProtonatedPdb,
// attaches the dependency chain (Geometry, SpatialIndex, Enrichment),
// then runs TripeptideBackboneShieldingResult against the local
// tensorcs15 replica and verifies:
//
//   - residues attempted == 76 (1UBQ length); matched is most of them
//     (terminal pair always skipped because phi/psi requires both
//     flanking residues; LYS/ARG/etc. with chi3/chi4 set should hit
//     the chi3+chi4 lookup path)
//   - per-atom tensors are finite where has_match is true
//   - frame_type discriminator: any SER residue should hit
//     orca_input_orientation rows (project SER PBE regen)
//   - WriteFeatures emits the central tensor, residual, method tag,
//     distance, and atom-match metadata NPYs
//
// Skipped if the tensorcs15 DSN is not configured (test machine
// doesn't have the local DB), or if the Larsen archive PDB is not at
// the expected path.

#include "TestEnvironment.h"
#include <gtest/gtest.h>

#include "AminoAcidType.h"
#include "ConformationAtom.h"
#include "EnrichmentResult.h"
#include "GeometryResult.h"
#include "OperationLog.h"
#include "PdbFileReader.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "RuntimeEnvironment.h"
#include "Session.h"
#include "SpatialIndexResult.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"
#include "TripeptideBackboneShieldingResult.h"
#include "TripeptideBackboneMethodTagTimeSeriesTrajectoryResult.h"
#include "TripeptideBackboneResidualVecTimeSeriesTrajectoryResult.h"
#include "TripeptideBackboneShieldingTimeSeriesTrajectoryResult.h"
#include "TripeptideDftTable.h"
#include "DirectionalTestHelpers.h"

#include <libpq-fe.h>
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

SphericalTensor UnpackFull9(const double* values) {
    SphericalTensor tensor;
    tensor.T0 = values[0];
    for (std::size_t component = 0; component < 3; ++component)
        tensor.T1[component] = values[1 + component];
    for (std::size_t component = 0; component < 5; ++component)
        tensor.T2[component] = values[4 + component];
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

}  // namespace

// Forward-declare the file-local PRODUCTION dihedral helper (named per-file
// namespace, external linkage) so the fixed-coordinate ±60° test below pins
// production DIRECTLY, not a copy of it (vet finding 2026-07). A production
// sign regression now fails this non-skippable test.
namespace nmr {
namespace tripeptide_backbone_dihedral {
double DihedralDegrees(const Vec3&, const Vec3&, const Vec3&, const Vec3&);
}
}  // namespace nmr

// Production C3/M14 seams. These functions are defined in named, per-file
// namespaces and are called by QueryNearest/WriteFeatures themselves.
namespace nmr {
namespace tripeptide_dft_query {
double CircularDeltaDegrees(double, double);
double DroppedChiSquaredDistance(const std::array<int, 4>&, const std::array<int, 4>&, int, int);
std::string OmittedChiScoreSql(int, int, const std::array<int, 4>&);
std::string NearestOrderBySql();
}  // namespace tripeptide_dft_query
}  // namespace nmr

namespace nmr {
namespace tripeptide_backbone_status {
int HisVariantDiagnosticCode(AminoAcid, int);
TripeptideMatchStatus UnsupportedHisStatus(AminoAcid, int);
}  // namespace tripeptide_backbone_status
}  // namespace nmr

namespace nmr {
namespace tripeptide_backbone_diagnostics {
std::array<double, 28> PackDiagnosticRow(const TripeptideBackboneShieldingResult::ResidueMatch&);
}
}  // namespace nmr


TEST(TripeptideBackboneDihedral, CanonicalFixedCoordinatesPinSignAndDegenerateNaN) {
    // Calls the PRODUCTION helper directly (not a copy), so a production sign
    // regression fails this non-skippable test (vet finding 2026-07).
    using nmr::tripeptide_backbone_dihedral::DihedralDegrees;
    const Vec3 a(1.0, 0.0, 0.0);
    const Vec3 b(0.0, 0.0, 0.0);
    const Vec3 c(0.0, 0.0, 1.0);
    const double root3_over_2 = std::sqrt(3.0) / 2.0;

    // Analytic ±60° fixtures. The former tripeptide triple product gives the
    // opposite signs, so a sign regression fails here.
    EXPECT_NEAR(DihedralDegrees(a, b, c, Vec3(0.5, -root3_over_2, 1.0)), 60.0, 1e-12);
    EXPECT_NEAR(DihedralDegrees(a, b, c, Vec3(0.5, root3_over_2, 1.0)), -60.0, 1e-12);

    // Zero central bond and collinear plane-defining triples are undefined,
    // matching the trajectory convention.
    EXPECT_TRUE(std::isnan(DihedralDegrees(a, b, b, Vec3(0.0, 1.0, 0.0))));
    EXPECT_TRUE(
        std::isnan(DihedralDegrees(Vec3(0.0, 0.0, 0.0), Vec3(1.0, 0.0, 0.0), Vec3(2.0, 0.0, 0.0), Vec3(2.0, 1.0, 0.0))));
    EXPECT_TRUE(
        std::isnan(DihedralDegrees(Vec3(0.0, 1.0, 0.0), Vec3(0.0, 0.0, 0.0), Vec3(1.0, 0.0, 0.0), Vec3(2.0, 0.0, 0.0))));
}


TEST(TripeptideDftFallback, ProductionSqlAndCircularScorePreferNearestDroppedChiBeforeCalcId) {
    using nmr::tripeptide_dft_query::CircularDeltaDegrees;
    using nmr::tripeptide_dft_query::DroppedChiSquaredDistance;
    using nmr::tripeptide_dft_query::OmittedChiScoreSql;

    // Analytic wrap-around oracle: -170 and +170 are 20°, not 340° apart.
    EXPECT_DOUBLE_EQ(CircularDeltaDegrees(-170.0, 170.0), 20.0);
    EXPECT_DOUBLE_EQ(CircularDeltaDegrees(-180.0, 180.0), 0.0);

    // Same retained chi1; the lower-id candidate is deliberately farther
    // on omitted chi2..4 than the higher-id candidate. This calls the
    // production score used to audit QueryNearest's selected SQL row.
    const std::array<int, 4> target{40, -170, 170, -160};
    const std::array<int, 4> lower_id_far{40, 120, 100, 80};
    const std::array<int, 4> higher_id_near{40, -180, -180, -140};
    const double far_score = DroppedChiSquaredDistance(lower_id_far, target, /*natural=*/4, /*used=*/1);
    const double near_score = DroppedChiSquaredDistance(higher_id_near, target, /*natural=*/4, /*used=*/1);
    EXPECT_LT(near_score, far_score);
    EXPECT_NEAR(std::sqrt(near_score), std::sqrt(600.0), 1e-12);

    // Pin the exact production SQL construction: every omitted natural
    // axis participates via circular LEAST distance, and calc_id appears
    // only after the score in QueryNearest's fixed ORDER BY.
    const std::string score_sql = OmittedChiScoreSql(4, 1, target);
    EXPECT_EQ(score_sql.find("chi1"), std::string::npos);
    EXPECT_NE(score_sql.find("chi2"), std::string::npos);
    EXPECT_NE(score_sql.find("chi3"), std::string::npos);
    EXPECT_NE(score_sql.find("chi4"), std::string::npos);
    EXPECT_NE(score_sql.find("LEAST"), std::string::npos);
    EXPECT_NE(score_sql.find("-170"), std::string::npos);
    EXPECT_NE(score_sql.find("170"), std::string::npos);
    EXPECT_NE(score_sql.find("-160"), std::string::npos);
    EXPECT_EQ(OmittedChiScoreSql(4, 4, target), "0.0");
    EXPECT_EQ(nmr::tripeptide_dft_query::NearestOrderBySql(),
              "ORDER BY dropped_chi_score ASC, calc_id ASC LIMIT 1");
}


TEST(TripeptideBackboneDiagnostics, ProductionPackerPinsExact28ColumnsAndFiveStatusCodes) {
    TripeptideBackboneShieldingResult::ResidueMatch rm;
    rm.status = TripeptideMatchStatus::Ok;
    rm.calc_id = 42;
    rm.backbone_rmsd = 0.125;
    rm.ca_match_dist = 0.25;
    rm.frame_type = "gaussian_standard_orientation";
    rm.n_atoms_matched = 17;
    rm.natural_chi_axes = 4;
    rm.n_chi_axes_used = 1;
    rm.dropped_chi_distance_deg = 30.0;
    rm.phi_actual = -61.5;
    rm.psi_actual = 139.25;
    rm.phi_db = -60;
    rm.psi_db = 140;
    rm.residue_type_code = static_cast<int>(AminoAcid::HIS);
    rm.his_variant_hint = 1;
    for (int k = 0; k < 4; ++k) {
        rm.chi_actual[k] = 10.5 + k;
        rm.target_chi_grid[k] = 20 * (k + 1);
        rm.chi_db[k] = -20 * (k + 1);
    }

    const auto row = nmr::tripeptide_backbone_diagnostics::PackDiagnosticRow(rm);
    EXPECT_DOUBLE_EQ(row[0], 42.0);
    EXPECT_DOUBLE_EQ(row[1], 1.0);
    EXPECT_DOUBLE_EQ(row[2], 1.0);
    EXPECT_DOUBLE_EQ(row[3], 0.125);
    EXPECT_DOUBLE_EQ(row[4], 0.25);
    EXPECT_DOUBLE_EQ(row[5], 1.0);
    EXPECT_DOUBLE_EQ(row[6], 17.0);
    EXPECT_DOUBLE_EQ(row[7], 4.0);
    EXPECT_DOUBLE_EQ(row[8], 1.0);
    EXPECT_DOUBLE_EQ(row[9], 30.0);
    EXPECT_DOUBLE_EQ(row[10], -61.5);
    EXPECT_DOUBLE_EQ(row[11], 139.25);
    EXPECT_DOUBLE_EQ(row[12], -60.0);
    EXPECT_DOUBLE_EQ(row[13], 140.0);
    for (int k = 0; k < 4; ++k) {
        EXPECT_DOUBLE_EQ(row[14 + k], 10.5 + k);
        EXPECT_DOUBLE_EQ(row[18 + k], 20.0 * (k + 1));
        EXPECT_DOUBLE_EQ(row[22 + k], -20.0 * (k + 1));
    }
    EXPECT_DOUBLE_EQ(row[26], static_cast<double>(AminoAcid::HIS));
    EXPECT_DOUBLE_EQ(row[27], 1.0);

    const std::array<TripeptideMatchStatus, 5> statuses = {
        TripeptideMatchStatus::Miss,
        TripeptideMatchStatus::Ok,
        TripeptideMatchStatus::UnsupportedHid,
        TripeptideMatchStatus::UnsupportedHie,
        TripeptideMatchStatus::PerceptionFailed,
    };
    for (std::size_t i = 0; i < statuses.size(); ++i) {
        rm.status = statuses[i];
        const auto status_row = nmr::tripeptide_backbone_diagnostics::PackDiagnosticRow(rm);
        EXPECT_DOUBLE_EQ(status_row[1], static_cast<double>(i));
        EXPECT_DOUBLE_EQ(status_row[2], i == 1 ? 1.0 : 0.0);
    }

    // Chi target/DB columns are floating diagnostics: axes beyond the
    // residue's natural depth are unavailable and therefore remain NaN.
    rm.status = TripeptideMatchStatus::Ok;
    rm.natural_chi_axes = 1;
    const auto shallow = nmr::tripeptide_backbone_diagnostics::PackDiagnosticRow(rm);
    EXPECT_TRUE(std::isfinite(shallow[18]));
    EXPECT_TRUE(std::isfinite(shallow[22]));
    for (int k = 1; k < 4; ++k) {
        EXPECT_TRUE(std::isnan(shallow[18 + k]));
        EXPECT_TRUE(std::isnan(shallow[22 + k]));
    }

    // Execute the exact production censor/diagnostic mapping. HIP is not
    // coerced from HID/HIE: only variant 2 is the supported HIP code.
    EXPECT_EQ(nmr::tripeptide_backbone_status::UnsupportedHisStatus(AminoAcid::HIS, 0),
              TripeptideMatchStatus::UnsupportedHid);
    EXPECT_EQ(nmr::tripeptide_backbone_status::UnsupportedHisStatus(AminoAcid::HIS, 1),
              TripeptideMatchStatus::UnsupportedHie);
    EXPECT_EQ(nmr::tripeptide_backbone_status::UnsupportedHisStatus(AminoAcid::HIS, 2),
              TripeptideMatchStatus::Miss);
    EXPECT_EQ(nmr::tripeptide_backbone_status::HisVariantDiagnosticCode(AminoAcid::HIS, 0), 2);
    EXPECT_EQ(nmr::tripeptide_backbone_status::HisVariantDiagnosticCode(AminoAcid::HIS, 1), 3);
    EXPECT_EQ(nmr::tripeptide_backbone_status::HisVariantDiagnosticCode(AminoAcid::HIS, 2), 1);
}


TEST(TripeptideBackboneM14, HidHieAreExplicitlyCensoredWhileHipUsesTheProductionLookup) {
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

    auto compute_variant = [&](int variant) {
        protein->MutableResidueAt(his_index).protonation_variant_index = variant;
        auto result = TripeptideBackboneShieldingResult::Compute(conf, *session.TripeptideDftTablePtr());
        EXPECT_NE(result, nullptr);
        return result;
    };

    auto hid = compute_variant(0);
    ASSERT_NE(hid, nullptr);
    EXPECT_EQ(hid->ResidueMatches()[his_index].status, TripeptideMatchStatus::UnsupportedHid);
    EXPECT_EQ(hid->ResidueMatches()[his_index].his_variant_hint, 2);
    for (std::size_t ai : protein->ResidueAt(his_index).atom_indices) {
        EXPECT_FALSE(conf.AtomAt(ai).tripeptide_bb_has_match);
    }

    auto hie = compute_variant(1);
    ASSERT_NE(hie, nullptr);
    EXPECT_EQ(hie->ResidueMatches()[his_index].status, TripeptideMatchStatus::UnsupportedHie);
    EXPECT_EQ(hie->ResidueMatches()[his_index].his_variant_hint, 3);
    for (std::size_t ai : protein->ResidueAt(his_index).atom_indices) {
        EXPECT_FALSE(conf.AtomAt(ai).tripeptide_bb_has_match);
    }

    auto hip = compute_variant(2);
    ASSERT_NE(hip, nullptr);
    EXPECT_EQ(hip->ResidueMatches()[his_index].status, TripeptideMatchStatus::Ok);
    EXPECT_EQ(hip->ResidueMatches()[his_index].his_variant_hint, 1);
    EXPECT_GT(hip->ResidueMatches()[his_index].n_atoms_matched, 0);
}


TEST(TripeptideBackboneDirectionalCovariance,
     ProductionLookupAndKabschRerunAfterProperRigidTransform) {
    if (RuntimeEnvironment::TensorCs15Dsn().empty()) {
        GTEST_SKIP() << "tensorcs15 DSN not configured";
    }
    Session session;
    ASSERT_EQ(session.LoadTripeptideDftTable(), kOk) << session.LastError();
    auto build = BuildFromProtonatedPdb(
        nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(build.Ok()) << build.error;
    auto protein = std::move(build.protein);
    ProteinConformation& original = protein->Conformation();
    ASSERT_TRUE(original.AttachResult(GeometryResult::Compute(original)));
    ASSERT_TRUE(original.AttachResult(EnrichmentResult::Compute(original)));
    auto original_result = TripeptideBackboneShieldingResult::Compute(
        original, *session.TripeptideDftTablePtr());
    ASSERT_NE(original_result, nullptr);
    ASSERT_GT(original_result->AtomsAssigned(), 0);

    // The tensorcs15 source contains chiral tripeptide poses, not their
    // mirrored database.  Its global contract is therefore SO(3): an
    // improper transform would change signed lookup angles and is not a
    // source-preserving rerun.
    const auto transform = nmr::test::directional::SeededTransform(
        0x5452495045505449ULL, false);
    ProteinConformation& moved = protein->AddConformation(
        nmr::test::directional::Positions(transform,
                                          original.Positions()),
        "tripeptide proper covariance rerun");
    ASSERT_TRUE(moved.AttachResult(GeometryResult::Compute(moved)));
    ASSERT_TRUE(moved.AttachResult(EnrichmentResult::Compute(moved)));
    auto moved_result = TripeptideBackboneShieldingResult::Compute(
        moved, *session.TripeptideDftTablePtr());
    ASSERT_NE(moved_result, nullptr);
    ASSERT_EQ(moved_result->ResiduesMatched(),
              original_result->ResiduesMatched());
    ASSERT_EQ(moved_result->AtomsAssigned(),
              original_result->AtomsAssigned());

    constexpr double kTensorAbsTolerance = 8.0e-8;
    constexpr double kTensorRelTolerance = 2.0e-10;
    constexpr double kVectorAbsTolerance = 8.0e-8;
    constexpr double kVectorRelTolerance = 2.0e-10;
    constexpr double kScalarTolerance = 8.0e-8;

    ASSERT_EQ(moved_result->ResidueMatches().size(),
              original_result->ResidueMatches().size());
    for (size_t residue_index = 0;
         residue_index < original_result->ResidueMatches().size();
         ++residue_index) {
        const auto& a = original_result->ResidueMatches()[residue_index];
        const auto& b = moved_result->ResidueMatches()[residue_index];
        EXPECT_EQ(b.status, a.status) << residue_index;
        EXPECT_EQ(b.calc_id, a.calc_id) << residue_index;
        EXPECT_EQ(b.n_atoms_matched, a.n_atoms_matched) << residue_index;
        auto expect_same_or_nan = [&](double actual, double expected,
                                      const char* label) {
            if (std::isnan(expected)) {
                EXPECT_TRUE(std::isnan(actual))
                    << label << " residue=" << residue_index;
            } else {
                EXPECT_NEAR(actual, expected, kScalarTolerance)
                    << label << " residue=" << residue_index;
            }
        };
        expect_same_or_nan(b.phi_actual, a.phi_actual, "phi");
        expect_same_or_nan(b.psi_actual, a.psi_actual, "psi");
        for (size_t chi = 0; chi < 4; ++chi) {
            expect_same_or_nan(b.chi_actual[chi], a.chi_actual[chi],
                               "chi");
        }
    }

    for (size_t atom_index = 0; atom_index < original.AtomCount();
         ++atom_index) {
        const auto& a = original.AtomAt(atom_index);
        const auto& b = moved.AtomAt(atom_index);
        ASSERT_EQ(b.tripeptide_bb_has_match,
                  a.tripeptide_bb_has_match) << atom_index;
        EXPECT_EQ(b.tripeptide_bb_method_tag,
                  a.tripeptide_bb_method_tag) << atom_index;
        if (!a.tripeptide_bb_has_match) continue;
        EXPECT_TRUE(nmr::test::directional::NearMatrix(
            b.tripeptide_bb_shielding_tensor,
            nmr::test::directional::EvenRank2(
                transform, a.tripeptide_bb_shielding_tensor),
            kTensorAbsTolerance, kTensorRelTolerance))
            << "tripeptide_bb_shielding.npy atom=" << atom_index;
        EXPECT_TRUE(nmr::test::directional::NearVector(
            b.tripeptide_bb_residual_vec,
            nmr::test::directional::Polar(
                transform, a.tripeptide_bb_residual_vec),
            kVectorAbsTolerance, kVectorRelTolerance))
            << "tripeptide_bb_residual_vec.npy atom=" << atom_index;
        EXPECT_NEAR(b.tripeptide_bb_match_distance,
                    a.tripeptide_bb_match_distance,
                    kScalarTolerance) << atom_index;
    }

    const fs::path out_dir = nmr::test::TestEnvironment::TempPath(
        "tripeptide_bb_directional_covariance");
    fs::create_directories(out_dir);
    ASSERT_EQ(moved_result->WriteFeatures(moved, out_dir.string()), 6);
    ASSERT_EQ(nmr::test::directional::RunNumpyAllowPickleFalse(
                  NMR_TEST_PYTHON_EXECUTABLE,
                  NMR_NPY_ALLOW_PICKLE_FALSE_SCRIPT,
                  {out_dir / "tripeptide_bb_shielding.npy",
                   out_dir / "tripeptide_bb_residual_vec.npy",
                   out_dir / "tripeptide_bb_match_distance.npy",
                   out_dir / "tripeptide_bb_method_tag.npy",
                   out_dir / "tripeptide_bb_match_atoms.npy",
                   out_dir / "tripeptide_bb_diagnostics.npy"}),
              0);
    const NpyArray shielding =
        ReadNpy(out_dir / "tripeptide_bb_shielding.npy");
    const NpyArray residual =
        ReadNpy(out_dir / "tripeptide_bb_residual_vec.npy");
    const NpyArray match_distance =
        ReadNpy(out_dir / "tripeptide_bb_match_distance.npy");
    const NpyArray method_tag =
        ReadNpy(out_dir / "tripeptide_bb_method_tag.npy");
    const NpyArray match_atoms =
        ReadNpy(out_dir / "tripeptide_bb_match_atoms.npy");
    const NpyArray diagnostics =
        ReadNpy(out_dir / "tripeptide_bb_diagnostics.npy");
    ASSERT_EQ(shielding.shape,
              (std::vector<size_t>{moved.AtomCount(), 9}));
    ASSERT_EQ(residual.shape,
              (std::vector<size_t>{moved.AtomCount(), 3}));
    ASSERT_EQ(match_distance.shape,
              (std::vector<size_t>{moved.AtomCount()}));
    ASSERT_EQ(method_tag.shape,
              (std::vector<size_t>{moved.AtomCount()}));
    ASSERT_EQ(match_atoms.shape,
              (std::vector<size_t>{moved.AtomCount(), 5}));
    ASSERT_EQ(diagnostics.shape,
              (std::vector<size_t>{protein->ResidueCount(), 28}));
    const double* shielding_values = DataAs<double>(shielding);
    const double* residual_values = DataAs<double>(residual);
    const double* distance_values = DataAs<double>(match_distance);
    const std::int8_t* method_values = DataAs<std::int8_t>(method_tag);
    const double* match_values = DataAs<double>(match_atoms);
    const double* diagnostic_values = DataAs<double>(diagnostics);
    for (size_t atom_index = 0; atom_index < moved.AtomCount();
         ++atom_index) {
        const auto& source = original.AtomAt(atom_index);
        EXPECT_EQ(method_values[atom_index],
                  static_cast<std::int8_t>(source.tripeptide_bb_method_tag))
            << "tripeptide_bb_method_tag.npy atom=" << atom_index;
        EXPECT_DOUBLE_EQ(
            match_values[atom_index * 5],
            static_cast<double>(protein->AtomAt(atom_index).residue_index));
        EXPECT_DOUBLE_EQ(match_values[atom_index * 5 + 1],
                         source.tripeptide_bb_has_match ? 1.0 : 0.0);
        EXPECT_DOUBLE_EQ(match_values[atom_index * 5 + 4],
                         static_cast<double>(source.tripeptide_bb_method_tag));
        if (!source.tripeptide_bb_has_match) {
            for (size_t component = 0; component < 9; ++component) {
                EXPECT_TRUE(std::isnan(
                    shielding_values[atom_index * 9 + component]));
            }
            for (size_t component = 0; component < 3; ++component) {
                EXPECT_TRUE(std::isnan(
                    residual_values[atom_index * 3 + component]));
            }
            EXPECT_TRUE(std::isnan(distance_values[atom_index]));
            EXPECT_TRUE(std::isnan(match_values[atom_index * 5 + 3]));
            continue;
        }
        double expected_tensor[9] = {};
        SphericalTensor::Decompose(
            nmr::test::directional::EvenRank2(
                transform, source.tripeptide_bb_shielding_tensor))
            .PackFull9(expected_tensor);
        for (size_t component = 0; component < 9; ++component) {
            EXPECT_NEAR(shielding_values[atom_index * 9 + component],
                        expected_tensor[component],
                        kTensorAbsTolerance)
                << "tripeptide_bb_shielding.npy atom=" << atom_index
                << " component=" << component;
        }
        const Vec3 expected_residual =
            nmr::test::directional::Polar(
                transform, source.tripeptide_bb_residual_vec);
        for (size_t component = 0; component < 3; ++component) {
            EXPECT_NEAR(residual_values[atom_index * 3 + component],
                        expected_residual(component),
                        kVectorAbsTolerance)
                << "tripeptide_bb_residual_vec.npy atom=" << atom_index
                << " component=" << component;
        }
        EXPECT_NEAR(distance_values[atom_index],
                    source.tripeptide_bb_match_distance,
                    kScalarTolerance)
            << "tripeptide_bb_match_distance.npy atom=" << atom_index;
        EXPECT_NEAR(match_values[atom_index * 5 + 3],
                    source.tripeptide_bb_match_distance,
                    kScalarTolerance)
            << "tripeptide_bb_match_atoms.npy atom=" << atom_index;
    }
    ASSERT_EQ(moved_result->ResidueMatches().size(), protein->ResidueCount());
    for (size_t residue_index = 0; residue_index < protein->ResidueCount();
         ++residue_index) {
        const auto expected =
            tripeptide_backbone_diagnostics::PackDiagnosticRow(
                moved_result->ResidueMatches()[residue_index]);
        for (size_t column = 0; column < expected.size(); ++column) {
            const double actual =
                diagnostic_values[residue_index * expected.size() + column];
            if (std::isnan(expected[column])) {
                EXPECT_TRUE(std::isnan(actual))
                    << "tripeptide_bb_diagnostics.npy residue="
                    << residue_index << " column=" << column;
            } else {
                EXPECT_NEAR(actual, expected[column], kScalarTolerance)
                    << "tripeptide_bb_diagnostics.npy residue="
                    << residue_index << " column=" << column;
            }
        }
    }

    // Raw trajectory boundary: feed the two attached results produced by the
    // real owner reruns into each owning trajectory result.  Frame 0 is the
    // source conformation and frame 1 is the proper rigid transform, so the
    // serialized H5 laws are checked without rotating an output back into a
    // calculator.
    ASSERT_TRUE(original.AttachResult(std::move(original_result)));
    ASSERT_TRUE(moved.AttachResult(std::move(moved_result)));
    auto trajectory_protein =
        TrajectoryProtein::CreateForTesting(std::move(protein));
    ASSERT_NE(trajectory_protein, nullptr);
    auto shielding_ts =
        TripeptideBackboneShieldingTimeSeriesTrajectoryResult::Create(
            *trajectory_protein);
    auto residual_ts =
        TripeptideBackboneResidualVecTimeSeriesTrajectoryResult::Create(
            *trajectory_protein);
    auto method_ts =
        TripeptideBackboneMethodTagTimeSeriesTrajectoryResult::Create(
            *trajectory_protein);
    ASSERT_NE(shielding_ts, nullptr);
    ASSERT_NE(residual_ts, nullptr);
    ASSERT_NE(method_ts, nullptr);
    Trajectory dummy("", "", "");
    for (auto* result : std::array<TrajectoryResult*, 3>{
             shielding_ts.get(), residual_ts.get(), method_ts.get()}) {
        result->Compute(original, *trajectory_protein, dummy, 31, 1.25);
        result->Compute(moved, *trajectory_protein, dummy, 37, 2.75);
        result->Finalize(*trajectory_protein, dummy);
    }
    const fs::path h5_path = nmr::test::TestEnvironment::TempPath(
        "tripeptide_bb_directional_raw.h5");
    (void)std::remove(h5_path.string().c_str());
    {
        HighFive::File file(h5_path.string(), HighFive::File::Truncate);
        shielding_ts->WriteH5Group(*trajectory_protein, file);
        residual_ts->WriteH5Group(*trajectory_protein, file);
        method_ts->WriteH5Group(*trajectory_protein, file);
    }

    const std::string shielding_path =
        "/trajectory/tripeptide_bb_shielding_time_series/xyz";
    const std::string residual_path =
        "/trajectory/tripeptide_bb_residual_vec_time_series/xyz";
    const std::string method_path =
        "/trajectory/tripeptide_bb_method_tag_time_series/method_tag";
    std::vector<std::size_t> shielding_dims;
    std::vector<std::size_t> residual_dims;
    std::vector<std::size_t> method_dims;
    const auto h5_shielding = ReadH5Flat<double>(
        h5_path, shielding_path, &shielding_dims);
    const auto h5_residual = ReadH5Flat<double>(
        h5_path, residual_path, &residual_dims);
    const auto h5_method = ReadH5Flat<std::uint8_t>(
        h5_path, method_path, &method_dims);
    const std::size_t atom_count = original.AtomCount();
    EXPECT_EQ(shielding_dims,
              (std::vector<std::size_t>{atom_count, 2u, 9u}));
    EXPECT_EQ(residual_dims,
              (std::vector<std::size_t>{atom_count, 2u, 3u}));
    EXPECT_EQ(method_dims,
              (std::vector<std::size_t>{atom_count, 2u}));
    ASSERT_EQ(h5_shielding.size(), atom_count * 2 * 9);
    ASSERT_EQ(h5_residual.size(), atom_count * 2 * 3);
    ASSERT_EQ(h5_method.size(), atom_count * 2);

    // This one group describes both serialized frames: the source rerun and
    // the independently transformed production rerun.
    ExpectRawFull9TrajectoryMetadata(
        h5_path,
        "/trajectory/tripeptide_bb_shielding_time_series",
        "even_rank2 under proper rotations: T'=R T R^T; typed tripeptide "
        "lookup/Kabsch alignment is L-amino-acid chirality-conditioned and "
        "has no improper-transform contract");

    for (const std::string& group : {
             std::string("/trajectory/tripeptide_bb_shielding_time_series"),
             std::string("/trajectory/tripeptide_bb_residual_vec_time_series"),
             std::string("/trajectory/tripeptide_bb_method_tag_time_series")}) {
        EXPECT_EQ(ReadH5Flat<std::size_t>(
                      h5_path, group + "/frame_indices"),
                  (std::vector<std::size_t>{31u, 37u})) << group;
        EXPECT_EQ(ReadH5Flat<double>(h5_path, group + "/frame_times"),
                  (std::vector<double>{1.25, 2.75})) << group;
        EXPECT_EQ(ReadH5Flat<std::uint8_t>(
                      h5_path, group + "/source_attached_per_frame"),
                  (std::vector<std::uint8_t>{1u, 1u})) << group;
    }

    int finite_h5_tensor_rows = 0;
    int finite_h5_vector_rows = 0;
    for (std::size_t atom_index = 0; atom_index < atom_count;
         ++atom_index) {
        const std::size_t tensor_source_base = (atom_index * 2) * 9;
        const std::size_t tensor_moved_base =
            (atom_index * 2 + 1) * 9;
        const std::size_t vector_source_base = (atom_index * 2) * 3;
        const std::size_t vector_moved_base =
            (atom_index * 2 + 1) * 3;
        const bool matched = original.AtomAt(atom_index)
                                 .tripeptide_bb_has_match;
        EXPECT_EQ(moved.AtomAt(atom_index).tripeptide_bb_has_match,
                  matched) << atom_index;
        EXPECT_EQ(h5_method[atom_index * 2],
                  original.AtomAt(atom_index).tripeptide_bb_method_tag)
            << method_path << " atom=" << atom_index << " source";
        EXPECT_EQ(h5_method[atom_index * 2 + 1],
                  moved.AtomAt(atom_index).tripeptide_bb_method_tag)
            << method_path << " atom=" << atom_index << " moved";
        EXPECT_EQ(h5_method[atom_index * 2 + 1],
                  h5_method[atom_index * 2])
            << method_path << " atom=" << atom_index << " identity";

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
        for (std::size_t component = 0; component < 3; ++component) {
            EXPECT_EQ(std::isnan(h5_residual[
                          vector_moved_base + component]),
                      std::isnan(h5_residual[
                          vector_source_base + component]))
                << residual_path << " atom=" << atom_index
                << " component=" << component;
            EXPECT_EQ(std::isnan(h5_residual[
                          vector_source_base + component]),
                      !matched)
                << residual_path << " atom=" << atom_index
                << " component=" << component << " applicability";
        }
        if (!matched) continue;
        ++finite_h5_tensor_rows;
        ++finite_h5_vector_rows;
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
        const Vec3 source_vector(
            h5_residual[vector_source_base],
            h5_residual[vector_source_base + 1],
            h5_residual[vector_source_base + 2]);
        const Vec3 actual_vector(
            h5_residual[vector_moved_base],
            h5_residual[vector_moved_base + 1],
            h5_residual[vector_moved_base + 2]);
        EXPECT_TRUE(nmr::test::directional::NearVector(
            actual_vector,
            nmr::test::directional::Polar(transform, source_vector),
            kVectorAbsTolerance, kVectorRelTolerance))
            << residual_path << " atom=" << atom_index;
    }
    EXPECT_GT(finite_h5_tensor_rows, 0);
    EXPECT_GT(finite_h5_vector_rows, 0);
    EXPECT_EQ(std::remove(h5_path.string().c_str()), 0) << h5_path;

    for (const char* name : {
             "tripeptide_bb_shielding.npy",
             "tripeptide_bb_residual_vec.npy",
             "tripeptide_bb_match_distance.npy",
             "tripeptide_bb_method_tag.npy",
             "tripeptide_bb_match_atoms.npy",
             "tripeptide_bb_diagnostics.npy"}) {
        EXPECT_EQ(std::remove((out_dir / name).string().c_str()), 0)
            << name;
    }
    EXPECT_EQ(std::remove(out_dir.string().c_str()), 0) << out_dir;
}


class TripeptideBackboneShieldingTest : public ::testing::Test {
protected:
    void SetUp() override {
        const std::string kLarsen1UbqPm6 = nmr::test::TestEnvironment::Larsen1UbqPm6Pdb();
        if (RuntimeEnvironment::TensorCs15Dsn().empty()) {
            GTEST_SKIP() << "tensorcs15 DSN not configured "
                            "(set [databases].tensorcs15 in "
                            "~/.nmr_tools.toml)";
        }
        if (!fs::exists(kLarsen1UbqPm6)) {
            GTEST_SKIP() << "Larsen 1UBQ PM6-D3H+ PDB not found at " << kLarsen1UbqPm6
                         << " (download via larsen ERDA archive — see "
                            "/mnt/expansion/larsen_archive/README.md)";
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


TEST_F(TripeptideBackboneShieldingTest, QueryNearestSelectsNearestRoundedDroppedChiPoseInsteadOfLowestCalcId) {
    // Supplemental live-DB forcing function: independently locate an AKA
    // phi/psi group whose lowest calc_id has a different chi pose, then ask
    // production QueryNearest to drop all four chi equality constraints. A
    // second read-only query ranks the candidates with an independently
    // expressed modular circular-distance oracle.
    std::unique_ptr<PGconn, decltype(&PQfinish)> conn(PQconnectdb(RuntimeEnvironment::TensorCs15Dsn().c_str()), &PQfinish);
    ASSERT_NE(conn.get(), nullptr);
    ASSERT_EQ(PQstatus(conn.get()), CONNECTION_OK) << PQerrorMessage(conn.get());

    const char* oracle_sql = R"SQL(
        SELECT r.phi, r.psi, r.chi1, r.chi2, r.chi3, r.chi4,
               CASE WHEN (ROUND(r.chi1 / 20.0) * 20)::int = 180
                    THEN -180 ELSE (ROUND(r.chi1 / 20.0) * 20)::int END AS target_chi1,
               CASE WHEN (ROUND(r.chi2 / 20.0) * 20)::int = 180
                    THEN -180 ELSE (ROUND(r.chi2 / 20.0) * 20)::int END AS target_chi2,
               CASE WHEN (ROUND(r.chi3 / 20.0) * 20)::int = 180
                    THEN -180 ELSE (ROUND(r.chi3 / 20.0) * 20)::int END AS target_chi3,
               CASE WHEN (ROUND(r.chi4 / 20.0) * 20)::int = 180
                    THEN -180 ELSE (ROUND(r.chi4 / 20.0) * 20)::int END AS target_chi4,
               (SELECT MIN(g.calc_id)
                  FROM raw_dft_calculations g
                 WHERE g.tripeptide = r.tripeptide
                   AND g.phi = r.phi AND g.psi = r.psi) AS group_min
          FROM raw_dft_calculations r
         WHERE r.tripeptide = 'AKA'
           AND r.chi1 IS NOT NULL AND r.chi2 IS NOT NULL
           AND r.chi3 IS NOT NULL AND r.chi4 IS NOT NULL
           AND r.calc_id > (
               SELECT MIN(g.calc_id)
                 FROM raw_dft_calculations g
                WHERE g.tripeptide = r.tripeptide
                  AND g.phi = r.phi AND g.psi = r.psi)
         ORDER BY r.calc_id ASC
         LIMIT 1
    )SQL";
    std::unique_ptr<PGresult, decltype(&PQclear)> oracle(PQexec(conn.get(), oracle_sql), &PQclear);
    ASSERT_NE(oracle.get(), nullptr);
    ASSERT_EQ(PQresultStatus(oracle.get()), PGRES_TUPLES_OK) << PQerrorMessage(conn.get());
    ASSERT_EQ(PQntuples(oracle.get()), 1) << "tensorcs15 AKA should contain multiple chi poses per phi/psi";

    auto value = [&](int col) {
        EXPECT_FALSE(PQgetisnull(oracle.get(), 0, col));
        if (PQgetisnull(oracle.get(), 0, col))
            return 0;
        return std::stoi(PQgetvalue(oracle.get(), 0, col));
    };
    const int phi = value(0);
    const int psi = value(1);
    const std::array<int, 4> chi = {value(2), value(3), value(4), value(5)};
    const std::array<int, 4> target_chi = {value(6), value(7), value(8), value(9)};
    const int group_min = value(10);

    const TripeptideDftRecord selected = session.TripeptideDftTablePtr()->QueryNearest('K',
                                                                                       phi,
                                                                                       psi,
                                                                                       chi[0],
                                                                                       chi[1],
                                                                                       chi[2],
                                                                                       chi[3],
                                                                                       /*n_chi_axes=*/0,
                                                                                       /*his_variant_hint=*/-1);
    ASSERT_TRUE(selected.IsHit());
    EXPECT_EQ(selected.natural_chi_axes, 4);
    EXPECT_EQ(selected.n_chi_axes_used, 0);
    EXPECT_EQ(selected.target_phi_grid_deg, phi);
    EXPECT_EQ(selected.target_psi_grid_deg, psi);
    EXPECT_EQ(selected.target_chi_grid_deg, target_chi);

    // Independent circular-distance oracle: modulo-to-[-180,180) is a
    // different formulation from production's LEAST(abs, 360-abs).
    auto modular_delta_sq = [](int axis, int target) {
        const std::string column = "chi" + std::to_string(axis);
        const std::string delta = "ABS(MOD((" + column + " - (" + std::to_string(target)
                                  + ") + 540)::numeric, 360) - 180)";
        return "(" + delta + " * " + delta + ")";
    };
    std::string oracle_score;
    for (int k = 0; k < 4; ++k) {
        if (!oracle_score.empty())
            oracle_score += " + ";
        oracle_score += modular_delta_sq(k + 1, target_chi[k]);
    }
    std::ostringstream nearest_sql;
    nearest_sql << "SELECT calc_id, (" << oracle_score << ")::double precision "
                << "FROM raw_dft_calculations WHERE tripeptide='AKA' "
                << "AND phi=" << phi << " AND psi=" << psi << " "
                << "AND chi1 IS NOT NULL AND chi2 IS NOT NULL "
                << "AND chi3 IS NOT NULL AND chi4 IS NOT NULL "
                << "ORDER BY 2 ASC, calc_id ASC LIMIT 1";
    std::unique_ptr<PGresult, decltype(&PQclear)> nearest(
        PQexec(conn.get(), nearest_sql.str().c_str()), &PQclear);
    ASSERT_NE(nearest.get(), nullptr);
    ASSERT_EQ(PQresultStatus(nearest.get()), PGRES_TUPLES_OK) << PQerrorMessage(conn.get());
    ASSERT_EQ(PQntuples(nearest.get()), 1);
    const int expected_calc_id = std::stoi(PQgetvalue(nearest.get(), 0, 0));
    const double expected_score = std::stod(PQgetvalue(nearest.get(), 0, 1));
    ASSERT_NE(expected_calc_id, group_min)
        << "fixture must distinguish nearest omitted-chi selection from the old calc_id-only order";
    EXPECT_EQ(selected.calc_id, expected_calc_id);
    EXPECT_NEAR(selected.dropped_chi_distance_deg, std::sqrt(expected_score), 1e-12);
}


TEST_F(TripeptideBackboneShieldingTest, RunsOn1UbqPm6) {
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

    auto tbb = TripeptideBackboneShieldingResult::Compute(conf, *session.TripeptideDftTablePtr());
    ASSERT_NE(tbb, nullptr);

    // 1UBQ has 76 residues; ResiduesAttempted counts standard residues
    // (Unknown skipped). Termini fail phi/psi.
    EXPECT_EQ(tbb->ResiduesAttempted(), static_cast<int>(protein->ResidueCount()));
    EXPECT_GE(tbb->ResiduesMatched(), 60);  // most internals match
    EXPECT_GT(tbb->AtomsAssigned(), 0);

    // Backbone Kabsch RMSD should be small — the canonical tripeptide
    // and the protein backbone share the same N/CA/C bond geometry by
    // construction. ~0.5 Å is the empirical rough upper bound; if mean
    // exceeds 1 Å something is structurally wrong with the alignment.
    EXPECT_LT(tbb->MeanBackboneRmsd(), 1.0);

    // Exercise the production file-local helper through Compute. The
    // result keeps the actual protein phi/psi before its DB lookup; compare
    // those retained values with the independent coordinate oracle above.
    // The old (n1 x n2).b2hat formula flips every non-zero sign.
    int n_dihedrals_pinned = 0;
    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        const auto prev_idx = protein->BackbonePredecessor(ri);
        const auto next_idx = protein->BackboneSuccessor(ri);
        if (!prev_idx || !next_idx)
            continue;
        const Residue& prev = protein->ResidueAt(*prev_idx);
        const Residue& res = protein->ResidueAt(ri);
        const Residue& next = protein->ResidueAt(*next_idx);
        if (prev.C == Residue::NONE || res.N == Residue::NONE || res.CA == Residue::NONE || res.C == Residue::NONE
            || next.N == Residue::NONE) {
            continue;
        }

        const double phi = nmr::tripeptide_backbone_dihedral::DihedralDegrees(conf.PositionAt(prev.C),
                                                                              conf.PositionAt(res.N),
                                                                              conf.PositionAt(res.CA),
                                                                              conf.PositionAt(res.C));
        const double psi = nmr::tripeptide_backbone_dihedral::DihedralDegrees(conf.PositionAt(res.N),
                                                                              conf.PositionAt(res.CA),
                                                                              conf.PositionAt(res.C),
                                                                              conf.PositionAt(next.N));
        ASSERT_TRUE(std::isfinite(phi));
        ASSERT_TRUE(std::isfinite(psi));
        const auto& match = tbb->ResidueMatches()[ri];
        EXPECT_NEAR(match.phi_actual, phi, 1e-10) << "residue " << ri;
        EXPECT_NEAR(match.psi_actual, psi, 1e-10) << "residue " << ri;
        n_dihedrals_pinned += 2;
    }
    EXPECT_GT(n_dihedrals_pinned, 100) << "1UBQ should exercise both file-local backbone dihedral paths";

    // Spot check: every matched atom has finite tensor.
    int n_finite = 0, n_total = 0;
    for (size_t i = 0; i < conf.AtomCount(); ++i) {
        const auto& ca = conf.AtomAt(i);
        if (!ca.tripeptide_bb_has_match)
            continue;
        ++n_total;
        const auto& m = ca.tripeptide_bb_shielding_tensor;
        bool finite = true;
        for (int r = 0; r < 3 && finite; ++r)
            for (int c = 0; c < 3 && finite; ++c)
                if (!std::isfinite(m(r, c)))
                    finite = false;
        if (finite)
            ++n_finite;
    }
    EXPECT_EQ(n_finite, n_total) << "non-finite tensor at " << (n_total - n_finite) << " atoms";

    // Spot check: SER residues that hit the DB should carry the
    // orca_input_orientation tag (project SER PBE regen). 1UBQ has
    // 4 SER residues. Not all may hit (e.g. terminal SER skipped).
    int n_ser_attempted = 0, n_ser_matched_pbe = 0;
    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        if (protein->ResidueAt(ri).type != AminoAcid::SER)
            continue;
        ++n_ser_attempted;
        const auto& m = tbb->ResidueMatches()[ri];
        if (m.calc_id == 0)
            continue;
        if (m.frame_type == "orca_input_orientation") {
            ++n_ser_matched_pbe;
        }
    }
    EXPECT_GE(n_ser_attempted, 1) << "1UBQ should have at least 1 SER";
    EXPECT_GT(n_ser_matched_pbe, 0) << "SER residues that match should hit ORCA PBE rows "
                                    << "(frame_type=orca_input_orientation), " << "got " << n_ser_matched_pbe << "/"
                                    << n_ser_attempted;

    // Test NPY emission.
    const std::string out_dir = nmr::test::TestEnvironment::TempPath("tripeptide_bb_smoke_out");
    fs::create_directories(out_dir);
    int n_npy = tbb->WriteFeatures(conf, out_dir);
    EXPECT_EQ(n_npy, 6);
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_bb_shielding.npy"));
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_bb_residual_vec.npy"));
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_bb_match_distance.npy"));
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_bb_method_tag.npy"));
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_bb_match_atoms.npy"));
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_bb_diagnostics.npy"));
    EXPECT_FALSE(fs::exists(out_dir + "/tripeptide_bb_residue_status.npy"));

    auto match_atoms = ReadNpy(out_dir + "/tripeptide_bb_match_atoms.npy");
    ASSERT_EQ(match_atoms.shape, (std::vector<size_t>{conf.AtomCount(), 5}));
    const double* ma = DataAs<double>(match_atoms);
    int n_match_atom_rows = 0;
    for (size_t i = 0; i < conf.AtomCount(); ++i) {
        const auto& ca = conf.AtomAt(i);
        EXPECT_DOUBLE_EQ(ma[i * 5 + 0], static_cast<double>(protein->AtomAt(i).residue_index));
        if (ca.tripeptide_bb_method_tag == 0)
            continue;
        ++n_match_atom_rows;
        EXPECT_DOUBLE_EQ(ma[i * 5 + 1], 1.0);
        EXPECT_GT(ma[i * 5 + 2], 0.0);
        EXPECT_TRUE(std::isfinite(ma[i * 5 + 3]));
        EXPECT_DOUBLE_EQ(ma[i * 5 + 4], static_cast<double>(ca.tripeptide_bb_method_tag));
    }
    EXPECT_GT(n_match_atom_rows, 0);

    // Read back the single frozen residue-diagnostics schema. This pins
    // both C3 fallback fields and M14 censor codes at the emitted boundary.
    auto diagnostics = ReadNpy(out_dir + "/tripeptide_bb_diagnostics.npy");
    ASSERT_EQ(diagnostics.shape, (std::vector<size_t>{protein->ResidueCount(), 28}));
    const double* diag = DataAs<double>(diagnostics);
    int n_fallback_rows = 0;
    int n_his_censored = 0;
    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        const auto& m = tbb->ResidueMatches()[ri];
        const double* row = &diag[ri * 28];
        EXPECT_DOUBLE_EQ(row[0], static_cast<double>(m.calc_id));
        EXPECT_DOUBLE_EQ(row[1], static_cast<double>(m.status));
        EXPECT_DOUBLE_EQ(row[2], m.status == TripeptideMatchStatus::Ok ? 1.0 : 0.0);
        EXPECT_DOUBLE_EQ(row[6], static_cast<double>(m.n_atoms_matched));
        EXPECT_DOUBLE_EQ(row[7], static_cast<double>(m.natural_chi_axes));
        EXPECT_DOUBLE_EQ(row[8], static_cast<double>(m.n_chi_axes_used));
        EXPECT_DOUBLE_EQ(row[26], static_cast<double>(m.residue_type_code));
        EXPECT_DOUBLE_EQ(row[27], static_cast<double>(m.his_variant_hint));

        if (m.status == TripeptideMatchStatus::Ok && m.n_chi_axes_used < m.natural_chi_axes) {
            ++n_fallback_rows;
            EXPECT_TRUE(std::isfinite(row[9]));
            EXPECT_GE(row[9], 0.0);
            EXPECT_DOUBLE_EQ(row[9], m.dropped_chi_distance_deg);
            for (int k = 0; k < 4; ++k) {
                if (k < m.natural_chi_axes) {
                    EXPECT_DOUBLE_EQ(row[18 + k], static_cast<double>(m.target_chi_grid[k]));
                    EXPECT_DOUBLE_EQ(row[22 + k], static_cast<double>(m.chi_db[k]));
                } else {
                    EXPECT_TRUE(std::isnan(row[18 + k]));
                    EXPECT_TRUE(std::isnan(row[22 + k]));
                }
            }
        }

        const Residue& residue = protein->ResidueAt(ri);
        if (residue.type == AminoAcid::HIS && residue.protonation_variant_index == 0) {
            ++n_his_censored;
            EXPECT_DOUBLE_EQ(row[1], 2.0);  // unsupported_hid
            EXPECT_DOUBLE_EQ(row[2], 0.0);
            EXPECT_DOUBLE_EQ(row[27], 2.0);  // HID hint
        } else if (residue.type == AminoAcid::HIS && residue.protonation_variant_index == 1) {
            ++n_his_censored;
            EXPECT_DOUBLE_EQ(row[1], 3.0);  // unsupported_hie
            EXPECT_DOUBLE_EQ(row[2], 0.0);
            EXPECT_DOUBLE_EQ(row[27], 3.0);  // HIE hint
        }
    }
    EXPECT_GT(n_fallback_rows, 0) << "1UBQ should exercise at least one dropped-chi fallback";
    EXPECT_GE(n_his_censored, 0);  // exact codes are pinned non-skipping above

    // Residual sanity (post-Kabsch). The 3-point N/CA/C Kabsch has a
    // mean ~0.5 Å fit residual (chi-grid coarseness on flanking
    // residues' contribution to backbone C/O positions); individual
    // BB atoms can be 1-2 Å off. Above 3 Å on a backbone atom
    // indicates a Kabsch failure or mismapping — fail loud.
    //
    // The residual itself IS the ML feature (residual_vec on
    // ConformationAtom) — we don't gate ML emission on it. This test
    // just sanity-checks the distribution.
    int n_bb_extreme = 0;
    double max_bb_residual = 0.0;
    int n_bb_total = 0;
    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        const auto& res = protein->ResidueAt(ri);
        for (size_t slot : {res.N, res.CA, res.C, res.O}) {
            if (slot == Residue::NONE)
                continue;
            const auto& ca = conf.AtomAt(slot);
            if (!ca.tripeptide_bb_has_match)
                continue;
            ++n_bb_total;
            const double d = ca.tripeptide_bb_match_distance;
            if (d > max_bb_residual)
                max_bb_residual = d;
            if (d > 3.0)
                ++n_bb_extreme;
        }
    }
    EXPECT_GT(n_bb_total, 100) << "expected backbone atoms to mostly get tensors";
    EXPECT_LT(n_bb_extreme, 5) << "more than 5 backbone atoms have post-Kabsch residual > 3 Å "
                               << "— Kabsch failure or substrate disagreement (max=" << max_bb_residual << " Å)";
}

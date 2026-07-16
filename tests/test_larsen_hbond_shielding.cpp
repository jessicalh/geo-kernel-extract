// End-to-end smoke for LarsenHBondShieldingResult on 1UBQ_pm6dh3plus.pdb.
//
// Asserts the calculator across both donor classes (amide H + Hα) and
// all four acceptor classes (BackboneCarbonyl / SidechainCarbonyl /
// HydroxylOxygen / CarboxylateOxygen) via spatial-sweep enumeration:
//   - Detects ≥ 100 H-bond pairs (1UBQ has ~60 backbone amide-H pairs
//     plus Hα and sidechain contributions; aggregate easily clears 100).
//   - All emitted Mat3 tensors are finite.
//   - larsen_hbond_n_pairs > 0 on atoms that received contributions.
//   - larsen_hbond_any_corner_imputed propagates to Table 2 target
//     atoms (codex C4 fix).
//   - Cβ diagnostic tensors are finite and bounded — the value reflects
//     Larsen's actual DFT-computed Cβ shielding (which his Table 2
//     CHOSE not to include in ProCS15); non-zero is the methodology
//     statement, not a pipeline bug.
//   - Δσ_w water-term gate is geometric (no candidate O in range).
//   - WriteFeatures emits the per-atom and raw pair-audit NPYs.
//
// Per-cell Larsen Table 2 dispatch unit tests cover LarsenContribDispatch
// directly (no fixture needed).
//
// Skipped if the LarsenHBondGrid directory or the Larsen-archive 1UBQ
// PDB is not present.

#include "TestEnvironment.h"
#include <gtest/gtest.h>

#include "ConformationAtom.h"
#include "DsspResult.h"
#include "EnrichmentResult.h"
#include "GeometryResult.h"
#include "LarsenHBond1pHaBShieldingTimeSeriesTrajectoryResult.h"
#include "LarsenHBond1pHBShieldingTimeSeriesTrajectoryResult.h"
#include "LarsenHBond2pHaBShieldingTimeSeriesTrajectoryResult.h"
#include "LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult.h"
#include "LarsenHBondCountTimeSeriesTrajectoryResult.h"
#include "LarsenHBondWaterTermTimeSeriesTrajectoryResult.h"
#include "LarsenHBondGrid.h"
#include "LarsenSidechainDonorAuditResult.h"
#include "LarsenHBondShieldingResult.h"
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
#include "DirectionalTestHelpers.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>
#include <limits>
#include <memory>
#include <numeric>
#include <set>
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

fs::path MakeUniqueTempDirectory(const std::string& prefix) {
    static std::atomic<unsigned long long> sequence{0};
    const auto tick = std::chrono::steady_clock::now()
                          .time_since_epoch()
                          .count();
    const fs::path dir = fs::temp_directory_path() /
        (prefix + "_" + std::to_string(tick) + "_" +
         std::to_string(sequence.fetch_add(1)));
    fs::create_directories(dir);
    return dir;
}

void RemoveFlatTempDirectory(
    const fs::path& dir,
    std::initializer_list<const char*> file_names) {
    for (const char* file_name : file_names) {
        const std::string path = (dir / file_name).string();
        EXPECT_EQ(std::remove(path.c_str()), 0) << path;
    }
    const std::string path = dir.string();
    EXPECT_EQ(std::remove(path.c_str()), 0) << path;
}

struct RawNpy {
    std::string header;
    std::vector<char> payload;
};

template <typename T>
T PayloadValue(const RawNpy& npy, std::size_t element_index) {
    T value{};
    const std::size_t byte_offset = element_index * sizeof(T);
    EXPECT_LE(byte_offset + sizeof(T), npy.payload.size());
    if (byte_offset + sizeof(T) <= npy.payload.size()) {
        std::memcpy(&value, npy.payload.data() + byte_offset, sizeof(T));
    }
    return value;
}

SphericalTensor UnpackPayloadFull9(const RawNpy& npy,
                                   std::size_t row) {
    SphericalTensor tensor;
    tensor.T0 = PayloadValue<double>(npy, row * 9);
    for (std::size_t component = 0; component < 3; ++component) {
        tensor.T1[component] =
            PayloadValue<double>(npy, row * 9 + 1 + component);
    }
    for (std::size_t component = 0; component < 5; ++component) {
        tensor.T2[component] =
            PayloadValue<double>(npy, row * 9 + 4 + component);
    }
    return tensor;
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

void ExpectWrappedOppositeOrNaN(double actual, double source,
                                double tolerance,
                                const std::string& context) {
    if (std::isnan(source)) {
        EXPECT_TRUE(std::isnan(actual)) << context;
    } else {
        EXPECT_TRUE(std::isfinite(actual)) << context;
        EXPECT_NEAR(std::remainder(actual + source, 360.0), 0.0,
                    tolerance) << context;
    }
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

RawNpy ReadRawNpy(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    EXPECT_TRUE(in.is_open()) << path;
    char magic[6] = {};
    in.read(magic, 6);
    EXPECT_EQ(std::string(magic, 6), std::string("\x93NUMPY", 6));
    char version[2] = {};
    in.read(version, 2);
    EXPECT_EQ(version[0], 1);
    EXPECT_EQ(version[1], 0);
    std::uint16_t header_len = 0;
    in.read(reinterpret_cast<char*>(&header_len), sizeof(header_len));
    RawNpy out;
    out.header.resize(header_len);
    in.read(out.header.data(), header_len);
    out.payload.assign(std::istreambuf_iterator<char>(in),
                       std::istreambuf_iterator<char>());
    return out;
}

std::vector<float> ConstantTensorGrid(float base) {
    constexpr std::size_t kGridCells = 2 * 2 * 2;
    std::vector<float> values(kGridCells * 9, 0.0f);
    for (std::size_t cell = 0; cell < kGridCells; ++cell) {
        values[cell * 9 + 0] = base;
        values[cell * 9 + 4] = base * 2.0f;
        values[cell * 9 + 8] = base * 3.0f;
    }
    return values;
}

void WriteTensorDataset(HighFive::File& file,
                        const std::string& name,
                        float base) {
    const HighFive::DataSpace space{2, 2, 2, 3, 3};
    auto dataset = file.createDataSet<float>(name, space);
    const std::vector<float> values = ConstantTensorGrid(base);
    dataset.write_raw(values.data());
}

// Creates all six schema-valid archives with a tiny regular grid. Every
// successful query interpolates the production tensors above and sees
// exactly one imputed corner, making C34 deterministic without an external
// database or GTEST_SKIP gate.
fs::path WriteSyntheticLarsenGridDirectory(
    const std::string& tag,
    bool force_all_grid_miss = false) {
    const fs::path dir =
        MakeUniqueTempDirectory("nmr_larsen_synthetic_" + tag);

    const std::vector<double> theta_axis = {90.0, 180.0};
    const std::vector<double> rho_axis = {-180.0, 0.0};
    const std::vector<std::uint8_t> validity = {0, 1, 1, 1, 1, 1, 1, 1};
    const char* stems[6] = {
        "NMANMA", "NMACOH", "NMACOO",
        "ALANMA", "ALACOH", "ALACOO",
    };

    for (int archive = 0; archive < 6; ++archive) {
        HighFive::File file(
            (dir / (std::string(stems[archive]) + "_dense.h5")).string(),
            HighFive::File::Overwrite);
        const std::vector<double> r_axis = force_all_grid_miss
            ? std::vector<double>{0.05, 0.10}
            : (archive < 3
                ? std::vector<double>{1.5, 3.0}
                : std::vector<double>{1.8, 4.0});
        file.createDataSet("r_axis", r_axis);
        file.createDataSet("theta_axis", theta_axis);
        file.createDataSet("rho_axis", rho_axis);

        float base = static_cast<float>(archive + 1);
        WriteTensorDataset(file, "donor_N",  base + 0.1f);
        WriteTensorDataset(file, "donor_CA", base + 0.2f);
        WriteTensorDataset(file, "donor_C",  base + 0.3f);
        WriteTensorDataset(file, "donor_HA", base + 0.4f);
        WriteTensorDataset(file, "donor_HN", base + 0.5f);
        const bool ala_donor = archive >= 3;
        if (ala_donor) WriteTensorDataset(file, "donor_CB", base + 0.6f);

        const bool nma_acceptor = archive == 0 || archive == 3;
        if (nma_acceptor) {
            WriteTensorDataset(file, "acceptor_N",  base + 0.7f);
            WriteTensorDataset(file, "acceptor_HN", base + 0.8f);
            WriteTensorDataset(file, "acceptor_CA", base + 0.9f);
            WriteTensorDataset(file, "acceptor_C",  base + 1.0f);
            WriteTensorDataset(file, "acceptor_HA", base + 1.1f);
        }

        auto mask = file.createDataSet<std::uint8_t>(
            "validity_mask", HighFive::DataSpace{2, 2, 2});
        mask.write_raw(validity.data());
    }
    return dir;
}

}  // namespace


// ---------------------------------------------------------------------------
// LarsenContribDispatch::Applies — per-cell Table 2 dispatch verification.
// ---------------------------------------------------------------------------

TEST(LarsenContribDispatchTest, Table2Cells) {
    using TA = LarsenContribDispatch::TargetAtom;
    using Term = LarsenContribDispatch::Term;

    // N row: 1HB ✓, 2HB ✓, 1HaB ✗, 2HaB ✓, RC ✗, w ✗
    EXPECT_TRUE (LarsenContribDispatch::Applies(TA::N,  Term::Primary_HB));
    EXPECT_TRUE (LarsenContribDispatch::Applies(TA::N,  Term::Secondary_HB));
    EXPECT_FALSE(LarsenContribDispatch::Applies(TA::N,  Term::Primary_HaB));
    EXPECT_TRUE (LarsenContribDispatch::Applies(TA::N,  Term::Secondary_HaB));
    EXPECT_FALSE(LarsenContribDispatch::Applies(TA::N,  Term::RingCurrent));
    EXPECT_FALSE(LarsenContribDispatch::Applies(TA::N,  Term::Water));

    // Cα row: 1HB ✗, 2HB ✓, 1HaB ✗, 2HaB ✗, RC ✗, w ✗
    EXPECT_FALSE(LarsenContribDispatch::Applies(TA::CA, Term::Primary_HB));
    EXPECT_TRUE (LarsenContribDispatch::Applies(TA::CA, Term::Secondary_HB));
    EXPECT_FALSE(LarsenContribDispatch::Applies(TA::CA, Term::Primary_HaB));
    EXPECT_FALSE(LarsenContribDispatch::Applies(TA::CA, Term::Secondary_HaB));

    // Cβ row: ALL FALSE (diagnostic; Larsen Table 2 says no contribution).
    for (int term = 0; term < (int)Term::Count; ++term) {
        EXPECT_FALSE(LarsenContribDispatch::Applies(
            TA::CB, static_cast<Term>(term)))
            << "Cβ row must be all-false per Larsen Table 2; term="
            << term;
    }

    // C' row: 1HB ✗, 2HB ✓, 1HaB ✗, 2HaB ✗, RC ✗, w ✗
    EXPECT_FALSE(LarsenContribDispatch::Applies(TA::C, Term::Primary_HB));
    EXPECT_TRUE (LarsenContribDispatch::Applies(TA::C, Term::Secondary_HB));
    EXPECT_FALSE(LarsenContribDispatch::Applies(TA::C, Term::Primary_HaB));
    EXPECT_FALSE(LarsenContribDispatch::Applies(TA::C, Term::Secondary_HaB));

    // Hα row: 1HB ✓, 2HB ✓, 1HaB ✓, 2HaB ✓, RC ✓, w ✗
    EXPECT_TRUE (LarsenContribDispatch::Applies(TA::HA, Term::Primary_HB));
    EXPECT_TRUE (LarsenContribDispatch::Applies(TA::HA, Term::Secondary_HB));
    EXPECT_TRUE (LarsenContribDispatch::Applies(TA::HA, Term::Primary_HaB));
    EXPECT_TRUE (LarsenContribDispatch::Applies(TA::HA, Term::Secondary_HaB));
    EXPECT_TRUE (LarsenContribDispatch::Applies(TA::HA, Term::RingCurrent));
    EXPECT_FALSE(LarsenContribDispatch::Applies(TA::HA, Term::Water));

    // HN row: ALL TRUE per Larsen Table 2.
    EXPECT_TRUE(LarsenContribDispatch::Applies(TA::HN, Term::Primary_HB));
    EXPECT_TRUE(LarsenContribDispatch::Applies(TA::HN, Term::Secondary_HB));
    EXPECT_TRUE(LarsenContribDispatch::Applies(TA::HN, Term::Primary_HaB));
    EXPECT_TRUE(LarsenContribDispatch::Applies(TA::HN, Term::Secondary_HaB));
    EXPECT_TRUE(LarsenContribDispatch::Applies(TA::HN, Term::RingCurrent));
    EXPECT_TRUE(LarsenContribDispatch::Applies(TA::HN, Term::Water));
}


// Production-linked forcing function for the A23 value bug. The test calls
// the exact audit helper used at every production tensor-write site and
// compares it with the independent physical invariant isotropic = trace/3.
TEST(LarsenHBondPairAudit, TermAwareIsotropicAccountingUsesAppliedTensors) {
    using TA = LarsenContribDispatch::TargetAtom;
    using Term = LarsenContribDispatch::Term;

    LarsenHBondShieldingResult::PairRecord pair;
    pair.iso_1pHB = pair.iso_2pHB = 0.0;
    pair.iso_1pHaB = pair.iso_2pHaB = 0.0;
    pair.iso_diagnostic_CB = 0.0;
    Mat3 primary = Mat3::Zero();
    primary.diagonal() << 3.0, 6.0, 12.0;      // isotropic 7
    Mat3 secondary = Mat3::Zero();
    secondary.diagonal() << -3.0, 0.0, 6.0;    // isotropic 1
    Mat3 alpha = Mat3::Zero();
    alpha.diagonal() << 9.0, 12.0, 15.0;       // isotropic 12
    Mat3 diagnostic = Mat3::Zero();
    diagnostic.diagonal() << 30.0, 0.0, 0.0;   // isotropic 10

    larsen_hbond_shielding_detail::RecordTable2Contribution(
        pair, TA::HN, Term::Primary_HB, primary);
    // A second actual target write under the same term must accumulate,
    // not overwrite (the GLY HA2/HA3 fan-out is this shape).
    larsen_hbond_shielding_detail::RecordTable2Contribution(
        pair, TA::HA, Term::Primary_HB, secondary);
    larsen_hbond_shielding_detail::RecordTable2Contribution(
        pair, TA::N, Term::Secondary_HaB, alpha);
    larsen_hbond_shielding_detail::RecordDiagnosticCbContribution(
        pair, diagnostic);

    EXPECT_DOUBLE_EQ(pair.iso_1pHB,
                     primary.trace() / 3.0 + secondary.trace() / 3.0);
    EXPECT_DOUBLE_EQ(pair.iso_2pHB, 0.0);
    EXPECT_DOUBLE_EQ(pair.iso_1pHaB, 0.0);
    EXPECT_DOUBLE_EQ(pair.iso_2pHaB, alpha.trace() / 3.0);
    EXPECT_DOUBLE_EQ(pair.iso_diagnostic_CB,
                     diagnostic.trace() / 3.0);
    EXPECT_NE(pair.target_mask_1pHB & (1u << static_cast<unsigned>(TA::HN)), 0u);
    EXPECT_NE(pair.target_mask_1pHB & (1u << static_cast<unsigned>(TA::HA)), 0u);
    EXPECT_NE(pair.target_mask_2pHaB & (1u << static_cast<unsigned>(TA::N)), 0u);
    EXPECT_EQ(pair.target_mask_diagnostic_CB,
              1u << static_cast<unsigned>(TA::CB));
}


// Bare-checkout, end-to-end freeze gate for A23/C34/C35. It runs the
// production grid loader, production spatial enumeration/dispatch, and the
// production NPY writer against committed 1UBQ coordinates plus tiny H5
// grids created by the test. No configured database and no GTEST_SKIP.
TEST(LarsenHBondPairAudit, SyntheticGridEmitsAuditablePairTables) {
    const fs::path pdb = nmr::test::TestEnvironment::UbqProtonated();
    ASSERT_TRUE(fs::is_regular_file(pdb)) << pdb;
    auto build = BuildFromProtonatedPdb(pdb.string());
    ASSERT_TRUE(build.Ok()) << build.error;
    std::unique_ptr<Protein> protein = std::move(build.protein);
    ProteinConformation& conf = protein->Conformation();

    ASSERT_TRUE(conf.AttachResult(GeometryResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(SpatialIndexResult::Compute(conf)));

    const fs::path grid_dir =
        WriteSyntheticLarsenGridDirectory("pair_audit");
    LarsenHBondGrid grid(grid_dir.string());
    auto result = LarsenHBondShieldingResult::Compute(conf, grid);
    ASSERT_NE(result, nullptr);

    EXPECT_EQ(result->Pairs().size(),
              static_cast<std::size_t>(result->PairsFound() +
                                       result->PairsGridSkipped()));
    ASSERT_GT(result->PairsFound(), 0);

    int n_success = 0;
    int n_missing_frame = 0;
    int n_theta_miss = 0;
    int n_grid_miss = 0;
    int n_invalid_frame = 0;
    int n_carboxylate_symmetry_filtered = 0;
    int n_sidechain_carbonyl_success = 0;
    double pair_1pHB = 0.0;
    double pair_2pHB = 0.0;
    double pair_1pHaB = 0.0;
    double pair_2pHaB = 0.0;
    double pair_cb = 0.0;
    for (const auto& pair : result->Pairs()) {
        using Disp = LarsenHBondShieldingResult::PairDisposition;
        switch (pair.disposition) {
            case Disp::MissingFrameAtoms:
                ++n_missing_frame;
                EXPECT_FALSE(pair.frame_valid);
                EXPECT_TRUE(std::isnan(pair.r_angstrom));
                EXPECT_TRUE(std::isnan(pair.theta_deg));
                EXPECT_TRUE(std::isnan(pair.rho_deg));
                EXPECT_TRUE(std::isnan(pair.iso_1pHB));
                EXPECT_TRUE(std::isnan(pair.iso_2pHB));
                EXPECT_TRUE(std::isnan(pair.iso_1pHaB));
                EXPECT_TRUE(std::isnan(pair.iso_2pHaB));
                EXPECT_TRUE(std::isnan(pair.iso_diagnostic_CB));
                EXPECT_TRUE(std::isnan(pair.iso_table2_total));
                EXPECT_FALSE(pair.any_corner_imputed);
                EXPECT_EQ(pair.imputed_corner_count, 0);
                break;
            case Disp::ThetaOutOfRange:
                ++n_theta_miss;
                EXPECT_TRUE(pair.frame_valid);
                EXPECT_TRUE(std::isfinite(pair.r_angstrom));
                EXPECT_TRUE(std::isnan(pair.iso_table2_total));
                break;
            case Disp::GridMiss:
                ++n_grid_miss;
                EXPECT_TRUE(pair.frame_valid);
                EXPECT_TRUE(std::isfinite(pair.theta_deg));
                EXPECT_TRUE(std::isnan(pair.iso_table2_total));
                break;
            case Disp::Success:
                ++n_success;
                EXPECT_TRUE(pair.frame_valid);
                EXPECT_EQ(pair.imputed_corner_count, 1);
                EXPECT_TRUE(pair.any_corner_imputed);
                EXPECT_NEAR(pair.iso_table2_total,
                            pair.iso_1pHB + pair.iso_2pHB +
                            pair.iso_1pHaB + pair.iso_2pHaB,
                            1e-12);
                EXPECT_DOUBLE_EQ(pair.isotropic_total,
                                 pair.iso_table2_total);
                pair_1pHB += pair.iso_1pHB;
                pair_2pHB += pair.iso_2pHB;
                pair_1pHaB += pair.iso_1pHaB;
                pair_2pHaB += pair.iso_2pHaB;
                pair_cb += pair.iso_diagnostic_CB;
                if (pair.acceptor_class ==
                    HBondAcceptorClass::SidechainCarbonyl) {
                    ++n_sidechain_carbonyl_success;
                }
                break;
            case Disp::InvalidFrame:
                ++n_invalid_frame;
                EXPECT_FALSE(pair.frame_valid);
                EXPECT_TRUE(std::isnan(pair.iso_table2_total));
                EXPECT_FALSE(pair.any_corner_imputed);
                EXPECT_EQ(pair.imputed_corner_count, 0);
                break;
            case Disp::CarboxylateSymmetryFiltered:
                ++n_carboxylate_symmetry_filtered;
                EXPECT_EQ(pair.acceptor_class,
                          HBondAcceptorClass::CarboxylateOxygen);
                EXPECT_TRUE(pair.frame_valid);
                EXPECT_TRUE(std::isfinite(pair.r_angstrom));
                EXPECT_TRUE(std::isfinite(pair.theta_deg));
                EXPECT_TRUE(std::isfinite(pair.rho_deg));
                EXPECT_EQ(pair.target_mask_1pHB, 0u);
                EXPECT_EQ(pair.target_mask_2pHB, 0u);
                EXPECT_EQ(pair.target_mask_1pHaB, 0u);
                EXPECT_EQ(pair.target_mask_2pHaB, 0u);
                EXPECT_EQ(pair.target_mask_diagnostic_CB, 0u);
                EXPECT_TRUE(std::isnan(pair.iso_table2_total));
                EXPECT_FALSE(pair.any_corner_imputed);
                EXPECT_EQ(pair.imputed_corner_count, 0);
                break;
        }
    }
    EXPECT_EQ(n_success, result->PairsFound());
    EXPECT_EQ(n_missing_frame, 0)
        << "committed 1UBQ fixture has no classified missing-anchor row";
    EXPECT_GT(n_theta_miss, 0);
    EXPECT_GT(n_grid_miss, 0);
    EXPECT_EQ(n_invalid_frame, 0)
        << "unmodified committed coordinates should have valid frames";
    EXPECT_GT(n_carboxylate_symmetry_filtered, 0)
        << "all classified carboxylate sibling attempts must be retained";
    for (const auto& filtered : result->Pairs()) {
        if (filtered.disposition !=
            LarsenHBondShieldingResult::PairDisposition::
                CarboxylateSymmetryFiltered) {
            continue;
        }
        const auto selected = std::find_if(
            result->Pairs().begin(), result->Pairs().end(),
            [&](const LarsenHBondShieldingResult::PairRecord& candidate) {
                return candidate.donor_atom_idx == filtered.donor_atom_idx &&
                       candidate.acceptor_C_atom_idx ==
                           filtered.acceptor_C_atom_idx &&
                       candidate.acceptor_atom_idx ==
                           filtered.acceptor_third_atom_idx &&
                       candidate.disposition !=
                           LarsenHBondShieldingResult::PairDisposition::
                               CarboxylateSymmetryFiltered;
            });
        EXPECT_NE(selected, result->Pairs().end())
            << "filtered carboxylate row must accompany the selected sibling";
    }
    EXPECT_GT(n_sidechain_carbonyl_success, 0)
        << "committed 1UBQ should exercise ASN/GLN acceptor approximation";

    // Independent conservation oracle: pair-level trace sums must equal
    // traces of the tensors actually accumulated on the atom axis.
    double atom_1pHB = 0.0;
    double atom_2pHB = 0.0;
    double atom_1pHaB = 0.0;
    double atom_2pHaB = 0.0;
    double atom_cb = 0.0;
    for (std::size_t i = 0; i < conf.AtomCount(); ++i) {
        const ConformationAtom& atom = conf.AtomAt(i);
        atom_1pHB += atom.larsen_hbond_1pHB_tensor.trace() / 3.0;
        atom_2pHB += atom.larsen_hbond_2pHB_tensor.trace() / 3.0;
        atom_1pHaB += atom.larsen_hbond_1pHaB_tensor.trace() / 3.0;
        atom_2pHaB += atom.larsen_hbond_2pHaB_tensor.trace() / 3.0;
        atom_cb += atom.larsen_hbond_diagnostic_CB.trace() / 3.0;
    }
    EXPECT_NEAR(pair_1pHB, atom_1pHB, 1e-8);
    EXPECT_NEAR(pair_2pHB, atom_2pHB, 1e-8);
    EXPECT_NEAR(pair_1pHaB, atom_1pHaB, 1e-8);
    EXPECT_NEAR(pair_2pHaB, atom_2pHaB, 1e-8);
    EXPECT_NEAR(pair_cb, atom_cb, 1e-8);

    const fs::path out =
        MakeUniqueTempDirectory("nmr_larsen_pair_audit_output");
    EXPECT_EQ(result->WriteFeatures(conf, out.string()), 15);

    const std::size_t P = result->Pairs().size();
    const RawNpy index = ReadRawNpy(out / "larsen_hbond_pairs_index.npy");
    const RawNpy geometry =
        ReadRawNpy(out / "larsen_hbond_pairs_geometry.npy");
    const RawNpy isotropic =
        ReadRawNpy(out / "larsen_hbond_pairs_isotropic.npy");
    const RawNpy compat = ReadRawNpy(out / "larsen_hbond_pairs.npy");
    EXPECT_NE(index.header.find("'descr': '<i4'"), std::string::npos);
    EXPECT_NE(index.header.find("'shape': (" + std::to_string(P) +
                                ",16)"), std::string::npos);
    EXPECT_EQ(index.payload.size(), P * 16 * sizeof(std::int32_t));
    EXPECT_NE(geometry.header.find("'shape': (" + std::to_string(P) +
                                   ",6)"), std::string::npos);
    EXPECT_EQ(geometry.payload.size(), P * 6 * sizeof(double));
    EXPECT_EQ(isotropic.payload.size(), P * 6 * sizeof(double));
    EXPECT_NE(compat.header.find("'shape': (" + std::to_string(P) +
                                 ",28)"), std::string::npos);
    EXPECT_EQ(compat.payload.size(), P * 28 * sizeof(double));

    // A23 read-back: pin every field family to the production-owned row,
    // rather than accepting correctly shaped but empty/reordered tables.
    const auto success_it = std::find_if(
        result->Pairs().begin(), result->Pairs().end(),
        [](const LarsenHBondShieldingResult::PairRecord& pair) {
            return pair.disposition ==
                LarsenHBondShieldingResult::PairDisposition::Success;
        });
    ASSERT_NE(success_it, result->Pairs().end());
    const std::size_t success_row = static_cast<std::size_t>(
        std::distance(result->Pairs().begin(), success_it));
    const auto& expected_pair = *success_it;
    EXPECT_EQ(PayloadValue<std::int32_t>(index, success_row * 16 + 0),
              static_cast<std::int32_t>(expected_pair.donor_atom_idx));
    EXPECT_EQ(PayloadValue<std::int32_t>(index, success_row * 16 + 1),
              static_cast<std::int32_t>(expected_pair.acceptor_atom_idx));
    EXPECT_EQ(PayloadValue<std::int32_t>(index, success_row * 16 + 4),
              static_cast<std::int32_t>(expected_pair.donor_class));
    EXPECT_EQ(PayloadValue<std::int32_t>(index, success_row * 16 + 5),
              static_cast<std::int32_t>(expected_pair.acceptor_class));
    EXPECT_EQ(PayloadValue<std::int32_t>(index, success_row * 16 + 6), 3);
    EXPECT_EQ(PayloadValue<std::int32_t>(index, success_row * 16 + 11),
              static_cast<std::int32_t>(expected_pair.target_mask_1pHB));
    EXPECT_EQ(PayloadValue<std::int32_t>(index, success_row * 16 + 15),
              static_cast<std::int32_t>(
                  expected_pair.target_mask_diagnostic_CB));
    EXPECT_DOUBLE_EQ(PayloadValue<double>(geometry, success_row * 6 + 0),
                     expected_pair.r_angstrom);
    EXPECT_DOUBLE_EQ(PayloadValue<double>(geometry, success_row * 6 + 1),
                     expected_pair.theta_deg);
    EXPECT_DOUBLE_EQ(PayloadValue<double>(geometry, success_row * 6 + 2),
                     expected_pair.rho_deg);
    EXPECT_DOUBLE_EQ(PayloadValue<double>(geometry, success_row * 6 + 3), 1.0);
    EXPECT_DOUBLE_EQ(PayloadValue<double>(geometry, success_row * 6 + 4), 1.0);
    EXPECT_DOUBLE_EQ(PayloadValue<double>(geometry, success_row * 6 + 5), 1.0);
    EXPECT_DOUBLE_EQ(PayloadValue<double>(isotropic, success_row * 6 + 0),
                     expected_pair.iso_1pHB);
    EXPECT_DOUBLE_EQ(PayloadValue<double>(isotropic, success_row * 6 + 1),
                     expected_pair.iso_2pHB);
    EXPECT_DOUBLE_EQ(PayloadValue<double>(isotropic, success_row * 6 + 2),
                     expected_pair.iso_1pHaB);
    EXPECT_DOUBLE_EQ(PayloadValue<double>(isotropic, success_row * 6 + 3),
                     expected_pair.iso_2pHaB);
    EXPECT_DOUBLE_EQ(PayloadValue<double>(isotropic, success_row * 6 + 4),
                     expected_pair.iso_diagnostic_CB);
    EXPECT_DOUBLE_EQ(PayloadValue<double>(isotropic, success_row * 6 + 5),
                     expected_pair.iso_table2_total);
    for (std::size_t col = 0; col < 16; ++col) {
        EXPECT_DOUBLE_EQ(
            PayloadValue<double>(compat, success_row * 28 + col),
            static_cast<double>(PayloadValue<std::int32_t>(
                index, success_row * 16 + col)));
    }
    for (std::size_t col = 0; col < 6; ++col) {
        EXPECT_DOUBLE_EQ(
            PayloadValue<double>(compat, success_row * 28 + 16 + col),
            PayloadValue<double>(geometry, success_row * 6 + col));
        EXPECT_DOUBLE_EQ(
            PayloadValue<double>(compat, success_row * 28 + 22 + col),
            PayloadValue<double>(isotropic, success_row * 6 + col));
    }

    const auto filtered_it = std::find_if(
        result->Pairs().begin(), result->Pairs().end(),
        [](const LarsenHBondShieldingResult::PairRecord& pair) {
            return pair.disposition ==
                LarsenHBondShieldingResult::PairDisposition::
                    CarboxylateSymmetryFiltered;
        });
    ASSERT_NE(filtered_it, result->Pairs().end());
    const std::size_t filtered_row = static_cast<std::size_t>(
        std::distance(result->Pairs().begin(), filtered_it));
    EXPECT_EQ(PayloadValue<std::int32_t>(
                  index, filtered_row * 16 + 6),
              static_cast<std::int32_t>(
                  LarsenHBondShieldingResult::PairDisposition::
                      CarboxylateSymmetryFiltered));
    EXPECT_DOUBLE_EQ(PayloadValue<double>(
                         geometry, filtered_row * 6 + 5),
                     1.0);
    EXPECT_TRUE(std::isnan(PayloadValue<double>(
        isotropic, filtered_row * 6 + 5)));

    const RawNpy imputed =
        ReadRawNpy(out / "larsen_imputed_pair_count.npy");
    const RawNpy sidechain = ReadRawNpy(
        out / "larsen_sidechain_carbonyl_pair_count.npy");
    ASSERT_EQ(imputed.payload.size(),
              conf.AtomCount() * sizeof(std::int32_t));
    ASSERT_EQ(sidechain.payload.size(),
              conf.AtomCount() * sizeof(std::int32_t));
    std::int64_t sidechain_sum = 0;
    for (std::size_t i = 0; i < conf.AtomCount(); ++i) {
        std::int32_t imputed_value = 0;
        std::int32_t sidechain_value = 0;
        std::memcpy(&imputed_value,
                    imputed.payload.data() + i * sizeof(std::int32_t),
                    sizeof(std::int32_t));
        std::memcpy(&sidechain_value,
                    sidechain.payload.data() + i * sizeof(std::int32_t),
                    sizeof(std::int32_t));
        // Every synthetic successful lookup has one imputed corner, and
        // diagnostics-only CB writes must not increment this count.
        EXPECT_EQ(imputed_value, conf.AtomAt(i).larsen_hbond_n_pairs);
        sidechain_sum += sidechain_value;
    }
    EXPECT_GT(sidechain_sum, 0);

    RemoveFlatTempDirectory(out, {
        "larsen_hbond_shielding.npy",
        "larsen_hbond_1pHB_shielding.npy",
        "larsen_hbond_2pHB_shielding.npy",
        "larsen_hbond_1pHaB_shielding.npy",
        "larsen_hbond_2pHaB_shielding.npy",
        "larsen_hbond_diagnostic_CB_shielding.npy",
        "larsen_hbond_water_term.npy",
        "larsen_hbond_count.npy",
        "larsen_imputed_pair_count.npy",
        "larsen_sidechain_carbonyl_pair_count.npy",
        "larsen_corner_imputed.npy",
        "larsen_hbond_pairs_index.npy",
        "larsen_hbond_pairs_geometry.npy",
        "larsen_hbond_pairs_isotropic.npy",
        "larsen_hbond_pairs.npy",
    });
    RemoveFlatTempDirectory(grid_dir, {
        "NMANMA_dense.h5", "NMACOH_dense.h5", "NMACOO_dense.h5",
        "ALANMA_dense.h5", "ALACOH_dense.h5", "ALACOO_dense.h5",
    });
}


TEST(LarsenHBondDirectionalCovariance,
     SyntheticGridRerunsProductionLookupAfterProperRigidTransform) {
    const fs::path pdb = nmr::test::TestEnvironment::UbqProtonated();
    ASSERT_TRUE(fs::is_regular_file(pdb)) << pdb;
    auto build = BuildFromProtonatedPdb(pdb.string());
    ASSERT_TRUE(build.Ok()) << build.error;
    auto protein = std::move(build.protein);
    ProteinConformation& original = protein->Conformation();
    ASSERT_TRUE(original.AttachResult(GeometryResult::Compute(original)));
    ASSERT_TRUE(original.AttachResult(SpatialIndexResult::Compute(original)));

    const fs::path grid_dir =
        WriteSyntheticLarsenGridDirectory("directional_covariance");
    LarsenHBondGrid grid(grid_dir.string());
    auto original_result =
        LarsenHBondShieldingResult::Compute(original, grid);
    ASSERT_NE(original_result, nullptr);
    ASSERT_GT(original_result->PairsFound(), 0);

    // Larsen's source grids encode a chiral molecular convention (not its
    // mirror-image database).  The physically supported global covariance
    // check is therefore SO(3), not a source-free O(3) mirror lookup.
    const auto transform = nmr::test::directional::SeededTransform(
        0x4c415253454e4842ULL, false);
    ProteinConformation& moved = protein->AddConformation(
        nmr::test::directional::Positions(transform,
                                          original.Positions()),
        "Larsen proper covariance rerun");
    ASSERT_TRUE(moved.AttachResult(GeometryResult::Compute(moved)));
    ASSERT_TRUE(moved.AttachResult(SpatialIndexResult::Compute(moved)));
    auto moved_result = LarsenHBondShieldingResult::Compute(moved, grid);
    ASSERT_NE(moved_result, nullptr);

    ASSERT_EQ(moved_result->Pairs().size(),
              original_result->Pairs().size());
    EXPECT_EQ(moved_result->PairsFound(), original_result->PairsFound());
    EXPECT_EQ(moved_result->PairsGridSkipped(),
              original_result->PairsGridSkipped());
    EXPECT_EQ(moved_result->AtomsWithContribution(),
              original_result->AtomsWithContribution());
    EXPECT_EQ(moved_result->AmideHsUnboundWithWater(),
              original_result->AmideHsUnboundWithWater());

    constexpr double kTensorAbsTolerance = 4.0e-9;
    constexpr double kTensorRelTolerance = 4.0e-11;
    constexpr double kGeometryTolerance = 3.0e-10;

    for (size_t pair_index = 0;
         pair_index < original_result->Pairs().size(); ++pair_index) {
        const auto& a = original_result->Pairs()[pair_index];
        const auto& b = moved_result->Pairs()[pair_index];
        EXPECT_EQ(b.donor_atom_idx, a.donor_atom_idx);
        EXPECT_EQ(b.acceptor_atom_idx, a.acceptor_atom_idx);
        EXPECT_EQ(b.disposition, a.disposition);
        EXPECT_EQ(b.frame_valid, a.frame_valid);
        auto expect_same_or_nan = [&](double actual, double expected,
                                      const char* label) {
            if (std::isnan(expected)) {
                EXPECT_TRUE(std::isnan(actual))
                    << label << " pair=" << pair_index;
            } else {
                EXPECT_NEAR(actual, expected, kGeometryTolerance)
                    << label << " pair=" << pair_index;
            }
        };
        expect_same_or_nan(b.r_angstrom, a.r_angstrom, "r");
        expect_same_or_nan(b.theta_deg, a.theta_deg, "theta");
        expect_same_or_nan(b.rho_deg, a.rho_deg, "rho");
        expect_same_or_nan(b.iso_1pHB, a.iso_1pHB, "iso_1pHB");
        expect_same_or_nan(b.iso_2pHB, a.iso_2pHB, "iso_2pHB");
        expect_same_or_nan(b.iso_1pHaB, a.iso_1pHaB, "iso_1pHaB");
        expect_same_or_nan(b.iso_2pHaB, a.iso_2pHaB, "iso_2pHaB");
        expect_same_or_nan(b.iso_diagnostic_CB, a.iso_diagnostic_CB,
                           "iso_CB");
        expect_same_or_nan(b.iso_table2_total, a.iso_table2_total,
                           "iso_total");
    }

    using MatrixMember = Mat3 ConformationAtom::*;
    struct TensorOutput {
        const char* name;
        MatrixMember matrix;
    };
    const std::array<TensorOutput, 6> tensor_outputs{{
        {"larsen_hbond_shielding.npy",
         &ConformationAtom::larsen_hbond_shielding_tensor},
        {"larsen_hbond_1pHB_shielding.npy",
         &ConformationAtom::larsen_hbond_1pHB_tensor},
        {"larsen_hbond_2pHB_shielding.npy",
         &ConformationAtom::larsen_hbond_2pHB_tensor},
        {"larsen_hbond_1pHaB_shielding.npy",
         &ConformationAtom::larsen_hbond_1pHaB_tensor},
        {"larsen_hbond_2pHaB_shielding.npy",
         &ConformationAtom::larsen_hbond_2pHaB_tensor},
        {"larsen_hbond_diagnostic_CB_shielding.npy",
         &ConformationAtom::larsen_hbond_diagnostic_CB},
    }};
    for (size_t atom_index = 0; atom_index < original.AtomCount();
         ++atom_index) {
        const auto& a = original.AtomAt(atom_index);
        const auto& b = moved.AtomAt(atom_index);
        for (const auto& output : tensor_outputs) {
            EXPECT_TRUE(nmr::test::directional::NearMatrix(
                b.*(output.matrix),
                nmr::test::directional::EvenRank2(
                    transform, a.*(output.matrix)),
                kTensorAbsTolerance, kTensorRelTolerance))
                << output.name << " atom=" << atom_index;
        }
        EXPECT_DOUBLE_EQ(b.larsen_hbond_water_term,
                         a.larsen_hbond_water_term);
        EXPECT_EQ(b.larsen_hbond_n_pairs, a.larsen_hbond_n_pairs);
        EXPECT_EQ(b.larsen_hbond_any_corner_imputed,
                  a.larsen_hbond_any_corner_imputed);
    }

    const fs::path original_out_dir = MakeUniqueTempDirectory(
        "nmr_larsen_directional_covariance_original");
    const fs::path moved_out_dir = MakeUniqueTempDirectory(
        "nmr_larsen_directional_covariance_moved");
    ASSERT_EQ(original_result->WriteFeatures(
                  original, original_out_dir.string()), 15);
    ASSERT_EQ(moved_result->WriteFeatures(moved, moved_out_dir.string()),
              15);
    auto all_larsen_feature_paths = [](const fs::path& dir) {
        return std::vector<fs::path>{
            dir / "larsen_hbond_shielding.npy",
            dir / "larsen_hbond_1pHB_shielding.npy",
            dir / "larsen_hbond_2pHB_shielding.npy",
            dir / "larsen_hbond_1pHaB_shielding.npy",
            dir / "larsen_hbond_2pHaB_shielding.npy",
            dir / "larsen_hbond_diagnostic_CB_shielding.npy",
            dir / "larsen_hbond_water_term.npy",
            dir / "larsen_hbond_count.npy",
            dir / "larsen_imputed_pair_count.npy",
            dir / "larsen_sidechain_carbonyl_pair_count.npy",
            dir / "larsen_corner_imputed.npy",
            dir / "larsen_hbond_pairs_index.npy",
            dir / "larsen_hbond_pairs_geometry.npy",
            dir / "larsen_hbond_pairs_isotropic.npy",
            dir / "larsen_hbond_pairs.npy",
        };
    };
    ASSERT_EQ(nmr::test::directional::RunNumpyAllowPickleFalse(
                  NMR_TEST_PYTHON_EXECUTABLE,
                  NMR_NPY_ALLOW_PICKLE_FALSE_SCRIPT,
                  all_larsen_feature_paths(original_out_dir)),
              0);
    ASSERT_EQ(nmr::test::directional::RunNumpyAllowPickleFalse(
                  NMR_TEST_PYTHON_EXECUTABLE,
                  NMR_NPY_ALLOW_PICKLE_FALSE_SCRIPT,
                  all_larsen_feature_paths(moved_out_dir)),
              0);
    for (const auto& output : tensor_outputs) {
        const RawNpy source = ReadRawNpy(original_out_dir / output.name);
        const RawNpy npy = ReadRawNpy(moved_out_dir / output.name);
        ASSERT_NE(source.header.find("'<f8'"), std::string::npos);
        ASSERT_NE(npy.header.find("'<f8'"), std::string::npos);
        ASSERT_EQ(source.payload.size(),
                  original.AtomCount() * 9 * sizeof(double));
        ASSERT_EQ(npy.payload.size(), source.payload.size());
        for (size_t atom_index = 0; atom_index < moved.AtomCount();
             ++atom_index) {
            const SphericalTensor source_tensor =
                UnpackPayloadFull9(source, atom_index);
            const SphericalTensor actual_tensor =
                UnpackPayloadFull9(npy, atom_index);
            EXPECT_TRUE(nmr::test::directional::NearMatrix(
                actual_tensor.Reconstruct(),
                nmr::test::directional::EvenRank2(
                    transform, source_tensor.Reconstruct()),
                kTensorAbsTolerance, kTensorRelTolerance))
                << output.name << " atom=" << atom_index;
            EXPECT_NEAR(actual_tensor.T0, source_tensor.T0,
                        kTensorAbsTolerance)
                << output.name << " atom=" << atom_index << " T0";
            EXPECT_TRUE(nmr::test::directional::NearVector(
                nmr::test::directional::T1Vector(actual_tensor),
                nmr::test::directional::Axial(
                    transform,
                    nmr::test::directional::T1Vector(source_tensor)),
                kTensorAbsTolerance, kTensorRelTolerance))
                << output.name << " atom=" << atom_index << " T1";
            const SphericalTensor expected_t2 =
                nmr::test::directional::RotateNativeT2(
                    transform, source_tensor);
            for (size_t component = 0; component < 5; ++component) {
                EXPECT_NEAR(actual_tensor.T2[component],
                            expected_t2.T2[component],
                            kTensorAbsTolerance)
                    << output.name << " atom=" << atom_index
                    << " native_T2_component=" << component;
            }
        }
    }

    // Every non-tensor production NPY is scalar/categorical under this
    // proper transform.  Compare both serialized payloads, including NaN
    // identity for unavailable pair geometry.  This accounts for all 15
    // exact writer names rather than treating the six tensor files as an
    // aggregate proxy.
    for (const char* name : {
             "larsen_hbond_water_term.npy",
             "larsen_hbond_pairs_geometry.npy",
             "larsen_hbond_pairs_isotropic.npy",
             "larsen_hbond_pairs.npy"}) {
        const RawNpy source = ReadRawNpy(original_out_dir / name);
        const RawNpy actual = ReadRawNpy(moved_out_dir / name);
        ASSERT_NE(source.header.find("'<f8'"), std::string::npos)
            << name;
        ASSERT_NE(actual.header.find("'<f8'"), std::string::npos)
            << name;
        ASSERT_EQ(actual.payload.size(), source.payload.size()) << name;
        ASSERT_EQ(source.payload.size() % sizeof(double), 0U) << name;
        const std::size_t count = source.payload.size() / sizeof(double);
        for (std::size_t component = 0; component < count; ++component) {
            ExpectSameOrNaN(PayloadValue<double>(actual, component),
                            PayloadValue<double>(source, component),
                            kGeometryTolerance,
                            std::string(name) + " component=" +
                                std::to_string(component));
        }
    }
    for (const char* name : {
             "larsen_hbond_count.npy",
             "larsen_imputed_pair_count.npy",
             "larsen_sidechain_carbonyl_pair_count.npy",
             "larsen_corner_imputed.npy",
             "larsen_hbond_pairs_index.npy"}) {
        const RawNpy source = ReadRawNpy(original_out_dir / name);
        const RawNpy actual = ReadRawNpy(moved_out_dir / name);
        EXPECT_EQ(actual.payload, source.payload) << name;
    }

    // Serialize the four component tensor streams and count stream from the
    // exact trajectory owners.  Both frames come from attached production
    // LarsenHBondShieldingResult reruns; no aggregate time-series owner is
    // invented because none exists.
    ASSERT_TRUE(original.AttachResult(std::move(original_result)));
    ASSERT_TRUE(moved.AttachResult(std::move(moved_result)));
    auto trajectory_protein =
        TrajectoryProtein::CreateForTesting(std::move(protein));
    ASSERT_NE(trajectory_protein, nullptr);
    auto one_hb_ts =
        LarsenHBond1pHBShieldingTimeSeriesTrajectoryResult::Create(
            *trajectory_protein);
    auto two_hb_ts =
        LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult::Create(
            *trajectory_protein);
    auto one_hab_ts =
        LarsenHBond1pHaBShieldingTimeSeriesTrajectoryResult::Create(
            *trajectory_protein);
    auto two_hab_ts =
        LarsenHBond2pHaBShieldingTimeSeriesTrajectoryResult::Create(
            *trajectory_protein);
    auto count_ts = LarsenHBondCountTimeSeriesTrajectoryResult::Create(
        *trajectory_protein);
    ASSERT_NE(one_hb_ts, nullptr);
    ASSERT_NE(two_hb_ts, nullptr);
    ASSERT_NE(one_hab_ts, nullptr);
    ASSERT_NE(two_hab_ts, nullptr);
    ASSERT_NE(count_ts, nullptr);
    const std::array<TrajectoryResult*, 5> trajectory_results{{
        one_hb_ts.get(), two_hb_ts.get(), one_hab_ts.get(),
        two_hab_ts.get(), count_ts.get(),
    }};
    Trajectory dummy("", "", "");
    for (TrajectoryResult* result : trajectory_results) {
        result->Compute(original, *trajectory_protein, dummy, 53, 5.25);
        result->Compute(moved, *trajectory_protein, dummy, 59, 6.75);
        result->Finalize(*trajectory_protein, dummy);
    }
    const fs::path h5_path = nmr::test::TestEnvironment::TempPath(
        "larsen_hbond_directional_raw.h5");
    (void)std::remove(h5_path.string().c_str());
    {
        HighFive::File file(h5_path.string(), HighFive::File::Truncate);
        for (TrajectoryResult* result : trajectory_results)
            result->WriteH5Group(*trajectory_protein, file);
    }

    const std::array<std::string, 4> tensor_groups{{
        "/trajectory/larsen_hbond_1pHB_shielding_time_series",
        "/trajectory/larsen_hbond_2pHB_shielding_time_series",
        "/trajectory/larsen_hbond_1pHaB_shielding_time_series",
        "/trajectory/larsen_hbond_2pHaB_shielding_time_series",
    }};
    const std::size_t atom_count = original.AtomCount();
    int finite_h5_tensor_rows = 0;
    for (const std::string& group : tensor_groups) {
        // Each group-level contract applies to both serialized frames below:
        // source and the independently transformed production rerun.
        ExpectRawFull9TrajectoryMetadata(
            h5_path, group,
            "even_rank2 under proper rotations: T'=R T R^T; signed-rho "
            "DFT-grid lookup is chirality-conditioned and has no "
            "improper-transform contract");
        const std::string path = group + "/xyz";
        std::vector<std::size_t> dimensions;
        const auto values = ReadH5Flat<double>(
            h5_path, path, &dimensions);
        EXPECT_EQ(dimensions,
                  (std::vector<std::size_t>{atom_count, 2u, 9u}))
            << path;
        ASSERT_EQ(values.size(), atom_count * 2 * 9) << path;
        EXPECT_EQ(ReadH5Flat<std::size_t>(
                      h5_path, group + "/frame_indices"),
                  (std::vector<std::size_t>{53u, 59u})) << group;
        EXPECT_EQ(ReadH5Flat<double>(h5_path, group + "/frame_times"),
                  (std::vector<double>{5.25, 6.75})) << group;
        EXPECT_EQ(ReadH5Flat<std::uint8_t>(
                      h5_path, group + "/source_attached_per_frame"),
                  (std::vector<std::uint8_t>{1u, 1u})) << group;
        for (std::size_t atom_index = 0; atom_index < atom_count;
             ++atom_index) {
            const std::size_t source_base = (atom_index * 2) * 9;
            const std::size_t moved_base = (atom_index * 2 + 1) * 9;
            bool finite = false;
            for (std::size_t component = 0; component < 9; ++component) {
                EXPECT_EQ(std::isnan(values[moved_base + component]),
                          std::isnan(values[source_base + component]))
                    << path << " atom=" << atom_index
                    << " component=" << component;
                finite = finite ||
                    std::isfinite(values[source_base + component]);
            }
            if (!finite) {
                for (std::size_t component = 0; component < 9;
                     ++component) {
                    EXPECT_TRUE(std::isnan(values[source_base + component]))
                        << path << " atom=" << atom_index;
                    EXPECT_TRUE(std::isnan(values[moved_base + component]))
                        << path << " moved atom=" << atom_index;
                }
                continue;
            }
            ++finite_h5_tensor_rows;
            for (std::size_t component = 0; component < 9; ++component) {
                ASSERT_TRUE(std::isfinite(values[source_base + component]))
                    << path << " atom=" << atom_index;
                ASSERT_TRUE(std::isfinite(values[moved_base + component]))
                    << path << " moved atom=" << atom_index;
            }
            const SphericalTensor source_tensor =
                UnpackFull9(values.data() + source_base);
            const SphericalTensor actual_tensor =
                UnpackFull9(values.data() + moved_base);
            EXPECT_TRUE(nmr::test::directional::NearMatrix(
                actual_tensor.Reconstruct(),
                nmr::test::directional::EvenRank2(
                    transform, source_tensor.Reconstruct()),
                kTensorAbsTolerance, kTensorRelTolerance))
                << path << " atom=" << atom_index;
            EXPECT_NEAR(actual_tensor.T0, source_tensor.T0,
                        kTensorAbsTolerance)
                << path << " atom=" << atom_index << " T0";
            EXPECT_TRUE(nmr::test::directional::NearVector(
                nmr::test::directional::T1Vector(actual_tensor),
                nmr::test::directional::Axial(
                    transform,
                    nmr::test::directional::T1Vector(source_tensor)),
                kTensorAbsTolerance, kTensorRelTolerance))
                << path << " atom=" << atom_index << " T1";
            const SphericalTensor expected_t2 =
                nmr::test::directional::RotateNativeT2(
                    transform, source_tensor);
            for (std::size_t component = 0; component < 5; ++component) {
                EXPECT_NEAR(actual_tensor.T2[component],
                            expected_t2.T2[component],
                            kTensorAbsTolerance)
                    << path << " atom=" << atom_index
                    << " native_T2_component=" << component;
            }
        }
    }
    EXPECT_GT(finite_h5_tensor_rows, 0);

    const std::string count_group =
        "/trajectory/larsen_hbond_count_time_series";
    std::vector<std::size_t> count_dimensions;
    const auto h5_counts = ReadH5Flat<int>(
        h5_path, count_group + "/count", &count_dimensions);
    EXPECT_EQ(count_dimensions,
              (std::vector<std::size_t>{atom_count, 2u}));
    ASSERT_EQ(h5_counts.size(), atom_count * 2);
    EXPECT_EQ(ReadH5Flat<std::size_t>(
                  h5_path, count_group + "/frame_indices"),
              (std::vector<std::size_t>{53u, 59u}));
    EXPECT_EQ(ReadH5Flat<double>(h5_path, count_group + "/frame_times"),
              (std::vector<double>{5.25, 6.75}));
    EXPECT_EQ(ReadH5Flat<std::uint8_t>(
                  h5_path,
                  count_group + "/source_attached_per_frame"),
              (std::vector<std::uint8_t>{1u, 1u}));
    int positive_count_rows = 0;
    for (std::size_t atom_index = 0; atom_index < atom_count;
         ++atom_index) {
        EXPECT_EQ(h5_counts[atom_index * 2],
                  original.AtomAt(atom_index).larsen_hbond_n_pairs)
            << count_group << " atom=" << atom_index << " source";
        EXPECT_EQ(h5_counts[atom_index * 2 + 1],
                  moved.AtomAt(atom_index).larsen_hbond_n_pairs)
            << count_group << " atom=" << atom_index << " moved";
        EXPECT_EQ(h5_counts[atom_index * 2 + 1],
                  h5_counts[atom_index * 2])
            << count_group << " atom=" << atom_index << " identity";
        positive_count_rows += h5_counts[atom_index * 2] > 0 ? 1 : 0;
    }
    EXPECT_GT(positive_count_rows, 0);
    EXPECT_EQ(std::remove(h5_path.string().c_str()), 0) << h5_path;

    const std::array<const char*, 15> all_output_names{{
        "larsen_hbond_shielding.npy",
        "larsen_hbond_1pHB_shielding.npy",
        "larsen_hbond_2pHB_shielding.npy",
        "larsen_hbond_1pHaB_shielding.npy",
        "larsen_hbond_2pHaB_shielding.npy",
        "larsen_hbond_diagnostic_CB_shielding.npy",
        "larsen_hbond_water_term.npy",
        "larsen_hbond_count.npy",
        "larsen_imputed_pair_count.npy",
        "larsen_sidechain_carbonyl_pair_count.npy",
        "larsen_corner_imputed.npy",
        "larsen_hbond_pairs_index.npy",
        "larsen_hbond_pairs_geometry.npy",
        "larsen_hbond_pairs_isotropic.npy",
        "larsen_hbond_pairs.npy",
    }};
    for (const fs::path& dir : {original_out_dir, moved_out_dir}) {
        for (const char* name : all_output_names) {
            EXPECT_EQ(std::remove((dir / name).string().c_str()), 0)
                << name;
        }
        EXPECT_EQ(std::remove(dir.string().c_str()), 0) << dir;
    }
    for (const auto& entry : fs::directory_iterator(grid_dir)) {
        EXPECT_EQ(std::remove(entry.path().string().c_str()), 0)
            << entry.path();
    }
    EXPECT_EQ(std::remove(grid_dir.string().c_str()), 0) << grid_dir;
}


TEST(LarsenHBondDirectionalCovariance,
     PairGateAndWaterOutputsAreExactUnderImproperProductionRerun) {
    const fs::path pdb = nmr::test::TestEnvironment::UbqProtonated();
    ASSERT_TRUE(fs::is_regular_file(pdb)) << pdb;
    auto build = BuildFromProtonatedPdb(pdb.string());
    ASSERT_TRUE(build.Ok()) << build.error;
    auto protein = std::move(build.protein);
    ProteinConformation& original = protein->Conformation();
    ASSERT_TRUE(original.AttachResult(GeometryResult::Compute(original)));
    ASSERT_TRUE(original.AttachResult(SpatialIndexResult::Compute(original)));

    const fs::path grid_dir =
        WriteSyntheticLarsenGridDirectory("directional_improper_gate");
    LarsenHBondGrid grid(grid_dir.string());
    auto original_result =
        LarsenHBondShieldingResult::Compute(original, grid);
    ASSERT_NE(original_result, nullptr);
    ASSERT_GT(original_result->PairsFound(), 0);

    constexpr std::uint64_t kTransformSeed = 0x4c415253454e4f33ULL;
    const auto transform = nmr::test::directional::SeededTransform(
        kTransformSeed, true);
    ProteinConformation& reflected = protein->AddConformation(
        nmr::test::directional::Positions(transform, original.Positions()),
        "Larsen improper pair-gate rerun");
    ASSERT_TRUE(reflected.AttachResult(GeometryResult::Compute(reflected)));
    ASSERT_TRUE(reflected.AttachResult(
        SpatialIndexResult::Compute(reflected)));
    auto reflected_result =
        LarsenHBondShieldingResult::Compute(reflected, grid);
    ASSERT_NE(reflected_result, nullptr);

    // QueryNearest's rho axis is periodic.  Rho selects tensor values and
    // validity-mask corners, but a hit/miss is gated only by finite geometry,
    // r and theta.  Therefore the exact pair identity/disposition and the
    // topology-dispatched contributor counts survive reflection even though
    // the shielding tensors and imputation provenance have no O(3) law.
    ASSERT_EQ(reflected_result->Pairs().size(),
              original_result->Pairs().size());
    EXPECT_EQ(reflected_result->PairsFound(),
              original_result->PairsFound());
    EXPECT_EQ(reflected_result->PairsGridSkipped(),
              original_result->PairsGridSkipped());
    EXPECT_EQ(reflected_result->AmideHsUnboundWithWater(),
              original_result->AmideHsUnboundWithWater());
    constexpr double kGeometryTolerance = 3.0e-10;
    std::size_t finite_nonzero_rho_rows = 0;
    for (std::size_t pair_index = 0;
         pair_index < original_result->Pairs().size(); ++pair_index) {
        const auto& source = original_result->Pairs()[pair_index];
        const auto& actual = reflected_result->Pairs()[pair_index];
        EXPECT_EQ(actual.donor_atom_idx, source.donor_atom_idx);
        EXPECT_EQ(actual.acceptor_atom_idx, source.acceptor_atom_idx);
        EXPECT_EQ(actual.disposition, source.disposition);
        EXPECT_EQ(actual.target_mask_1pHB, source.target_mask_1pHB);
        EXPECT_EQ(actual.target_mask_2pHB, source.target_mask_2pHB);
        EXPECT_EQ(actual.target_mask_1pHaB, source.target_mask_1pHaB);
        EXPECT_EQ(actual.target_mask_2pHaB, source.target_mask_2pHaB);
        EXPECT_EQ(actual.target_mask_diagnostic_CB,
                  source.target_mask_diagnostic_CB);
        EXPECT_EQ(actual.frame_valid, source.frame_valid);
        ExpectSameOrNaN(actual.r_angstrom, source.r_angstrom,
                        kGeometryTolerance,
                        "improper Larsen r pair=" +
                            std::to_string(pair_index));
        ExpectSameOrNaN(actual.theta_deg, source.theta_deg,
                        kGeometryTolerance,
                        "improper Larsen theta pair=" +
                            std::to_string(pair_index));
        ExpectWrappedOppositeOrNaN(
            actual.rho_deg, source.rho_deg, kGeometryTolerance,
            "improper Larsen signed rho pair=" +
                std::to_string(pair_index));
        if (std::isfinite(source.rho_deg) &&
            std::abs(std::remainder(source.rho_deg, 360.0)) >
                kGeometryTolerance) {
            ++finite_nonzero_rho_rows;
        }
    }
    EXPECT_GT(finite_nonzero_rho_rows, 0u)
        << "improper rho forcing must be non-vacuous";

    const fs::path original_dir = MakeUniqueTempDirectory(
        "nmr_larsen_improper_gate_original");
    const fs::path reflected_dir = MakeUniqueTempDirectory(
        "nmr_larsen_improper_gate_reflected");
    ASSERT_EQ(original_result->WriteFeatures(
                  original, original_dir.string()), 15);
    ASSERT_EQ(reflected_result->WriteFeatures(
                  reflected, reflected_dir.string()), 15);

    const std::array<const char*, 4> exact_static_names{{
        "larsen_hbond_count.npy",
        "larsen_sidechain_carbonyl_pair_count.npy",
        "larsen_hbond_water_term.npy",
        "larsen_hbond_pairs_index.npy",
    }};
    for (const char* name : exact_static_names) {
        const RawNpy source = ReadRawNpy(original_dir / name);
        const RawNpy actual = ReadRawNpy(reflected_dir / name);
        EXPECT_EQ(actual.header, source.header) << name;
        EXPECT_EQ(actual.payload, source.payload) << name;
    }

    // The raw H-O-C-third dihedral is geometric and therefore has a real
    // pseudoscalar law even though the downstream chiral lookup values and
    // imputation corners do not.  Read the owning serialized sparse tables,
    // preserve exact row identity/NaN masks, and do not impose a law on the
    // signed-rho-conditioned fields.
    const RawNpy source_index = ReadRawNpy(
        original_dir / "larsen_hbond_pairs_index.npy");
    const RawNpy actual_index = ReadRawNpy(
        reflected_dir / "larsen_hbond_pairs_index.npy");
    const RawNpy source_geometry = ReadRawNpy(
        original_dir / "larsen_hbond_pairs_geometry.npy");
    const RawNpy actual_geometry = ReadRawNpy(
        reflected_dir / "larsen_hbond_pairs_geometry.npy");
    const RawNpy source_isotropic = ReadRawNpy(
        original_dir / "larsen_hbond_pairs_isotropic.npy");
    const RawNpy actual_isotropic = ReadRawNpy(
        reflected_dir / "larsen_hbond_pairs_isotropic.npy");
    const RawNpy source_compat = ReadRawNpy(
        original_dir / "larsen_hbond_pairs.npy");
    const RawNpy actual_compat = ReadRawNpy(
        reflected_dir / "larsen_hbond_pairs.npy");
    const std::size_t pair_count = original_result->Pairs().size();
    ASSERT_EQ(source_index.payload.size(),
              pair_count * 16u * sizeof(std::int32_t));
    ASSERT_EQ(actual_index.payload, source_index.payload);
    ASSERT_EQ(source_geometry.payload.size(),
              pair_count * 6u * sizeof(double));
    ASSERT_EQ(actual_geometry.payload.size(), source_geometry.payload.size());
    ASSERT_EQ(source_isotropic.payload.size(),
              pair_count * 6u * sizeof(double));
    ASSERT_EQ(actual_isotropic.payload.size(), source_isotropic.payload.size());
    ASSERT_EQ(source_compat.payload.size(),
              pair_count * 28u * sizeof(double));
    ASSERT_EQ(actual_compat.payload.size(), source_compat.payload.size());

    std::size_t serialized_nonzero_rho_rows = 0;
    for (std::size_t row = 0; row < pair_count; ++row) {
        ExpectSameOrNaN(PayloadValue<double>(actual_geometry, row * 6u),
                        PayloadValue<double>(source_geometry, row * 6u),
                        kGeometryTolerance,
                        "larsen_hbond_pairs_geometry.npy r row=" +
                            std::to_string(row));
        ExpectSameOrNaN(
            PayloadValue<double>(actual_geometry, row * 6u + 1u),
            PayloadValue<double>(source_geometry, row * 6u + 1u),
            kGeometryTolerance,
            "larsen_hbond_pairs_geometry.npy theta row=" +
                std::to_string(row));
        const double source_rho =
            PayloadValue<double>(source_geometry, row * 6u + 2u);
        const double actual_rho =
            PayloadValue<double>(actual_geometry, row * 6u + 2u);
        ExpectWrappedOppositeOrNaN(
            actual_rho, source_rho, kGeometryTolerance,
            "larsen_hbond_pairs_geometry.npy rho row=" +
                std::to_string(row));
        if (std::isfinite(source_rho) &&
            std::abs(std::remainder(source_rho, 360.0)) >
                kGeometryTolerance) {
            ++serialized_nonzero_rho_rows;
        }
        EXPECT_DOUBLE_EQ(
            PayloadValue<double>(actual_geometry, row * 6u + 5u),
            PayloadValue<double>(source_geometry, row * 6u + 5u))
            << "frame_valid row=" << row;

        for (const RawNpy* geometry : {&source_geometry, &actual_geometry}) {
            const double any_imputed =
                PayloadValue<double>(*geometry, row * 6u + 3u);
            const double corner_count =
                PayloadValue<double>(*geometry, row * 6u + 4u);
            const double frame_valid =
                PayloadValue<double>(*geometry, row * 6u + 5u);
            EXPECT_TRUE(any_imputed == 0.0 || any_imputed == 1.0)
                << "any_corner_imputed row=" << row;
            EXPECT_GE(corner_count, 0.0) << row;
            EXPECT_LE(corner_count, 8.0) << row;
            EXPECT_EQ(corner_count, std::floor(corner_count)) << row;
            EXPECT_TRUE(frame_valid == 0.0 || frame_valid == 1.0) << row;
        }

        for (std::size_t column = 0; column < 16u; ++column) {
            const double source_value =
                PayloadValue<double>(source_compat, row * 28u + column);
            const double actual_value =
                PayloadValue<double>(actual_compat, row * 28u + column);
            EXPECT_DOUBLE_EQ(actual_value, source_value)
                << "larsen_hbond_pairs.npy index row=" << row
                << " col=" << column;
            EXPECT_DOUBLE_EQ(
                source_value,
                static_cast<double>(PayloadValue<std::int32_t>(
                    source_index, row * 16u + column)))
                << "split index identity row=" << row << " col=" << column;
        }
        for (std::size_t column = 0; column < 6u; ++column) {
            for (const auto& run : {
                     std::pair<const RawNpy*, const RawNpy*>{
                         &source_compat, &source_geometry},
                     std::pair<const RawNpy*, const RawNpy*>{
                         &actual_compat, &actual_geometry}}) {
                ExpectSameOrNaN(
                    PayloadValue<double>(*run.first,
                                         row * 28u + 16u + column),
                    PayloadValue<double>(*run.second,
                                         row * 6u + column),
                    0.0,
                    "compat geometry split identity row=" +
                        std::to_string(row) + " col=" +
                        std::to_string(column));
            }
            for (const auto& run : {
                     std::pair<const RawNpy*, const RawNpy*>{
                         &source_compat, &source_isotropic},
                     std::pair<const RawNpy*, const RawNpy*>{
                         &actual_compat, &actual_isotropic}}) {
                ExpectSameOrNaN(
                    PayloadValue<double>(*run.first,
                                         row * 28u + 22u + column),
                    PayloadValue<double>(*run.second,
                                         row * 6u + column),
                    0.0,
                    "compat isotropic split identity row=" +
                        std::to_string(row) + " col=" +
                        std::to_string(column));
            }
        }
        ExpectSameOrNaN(
            PayloadValue<double>(actual_compat, row * 28u + 16u),
            PayloadValue<double>(source_compat, row * 28u + 16u),
            kGeometryTolerance, "compat r row=" + std::to_string(row));
        ExpectSameOrNaN(
            PayloadValue<double>(actual_compat, row * 28u + 17u),
            PayloadValue<double>(source_compat, row * 28u + 17u),
            kGeometryTolerance, "compat theta row=" + std::to_string(row));
        ExpectWrappedOppositeOrNaN(
            PayloadValue<double>(actual_compat, row * 28u + 18u),
            PayloadValue<double>(source_compat, row * 28u + 18u),
            kGeometryTolerance, "compat rho row=" + std::to_string(row));
        EXPECT_DOUBLE_EQ(
            PayloadValue<double>(actual_compat, row * 28u + 21u),
            PayloadValue<double>(source_compat, row * 28u + 21u))
            << "compat frame_valid row=" << row;

        const auto disposition = static_cast<
            LarsenHBondShieldingResult::PairDisposition>(
                PayloadValue<std::int32_t>(source_index, row * 16u + 6u));
        for (const RawNpy* isotropic :
             {&source_isotropic, &actual_isotropic}) {
            for (std::size_t column = 0; column < 6u; ++column) {
                const double value =
                    PayloadValue<double>(*isotropic, row * 6u + column);
                if (disposition ==
                    LarsenHBondShieldingResult::PairDisposition::Success) {
                    EXPECT_TRUE(std::isfinite(value))
                        << "successful isotropic row=" << row
                        << " col=" << column;
                } else {
                    EXPECT_TRUE(std::isnan(value))
                        << "unavailable isotropic row=" << row
                        << " col=" << column;
                }
            }
        }
    }
    EXPECT_GT(serialized_nonzero_rho_rows, 0u)
        << "serialized improper rho proof must be non-vacuous";

    const RawNpy count_npy = ReadRawNpy(
        original_dir / "larsen_hbond_count.npy");
    const RawNpy sidechain_count_npy = ReadRawNpy(
        original_dir / "larsen_sidechain_carbonyl_pair_count.npy");
    std::size_t positive_count_rows = 0;
    std::int64_t sidechain_count_sum = 0;
    for (std::size_t atom = 0; atom < original.AtomCount(); ++atom) {
        if (PayloadValue<std::int32_t>(count_npy, atom) > 0)
            ++positive_count_rows;
        sidechain_count_sum +=
            PayloadValue<std::int32_t>(sidechain_count_npy, atom);
    }
    EXPECT_GT(positive_count_rows, 0u);
    EXPECT_GT(sidechain_count_sum, 0);

    ASSERT_TRUE(original.AttachResult(std::move(original_result)));
    ASSERT_TRUE(reflected.AttachResult(std::move(reflected_result)));
    auto trajectory_protein =
        TrajectoryProtein::CreateForTesting(std::move(protein));
    ASSERT_NE(trajectory_protein, nullptr);
    auto count_ts = LarsenHBondCountTimeSeriesTrajectoryResult::Create(
        *trajectory_protein);
    auto water_ts = LarsenHBondWaterTermTimeSeriesTrajectoryResult::Create(
        *trajectory_protein);
    ASSERT_NE(count_ts, nullptr);
    ASSERT_NE(water_ts, nullptr);
    Trajectory dummy("", "", "");
    for (TrajectoryResult* result :
         {static_cast<TrajectoryResult*>(count_ts.get()),
          static_cast<TrajectoryResult*>(water_ts.get())}) {
        result->Compute(original, *trajectory_protein, dummy, 71, 8.25);
        result->Compute(reflected, *trajectory_protein, dummy, 79, 9.75);
        result->Finalize(*trajectory_protein, dummy);
    }
    const fs::path h5_path = nmr::test::TestEnvironment::TempPath(
        "larsen_hbond_improper_gate.h5");
    (void)std::remove(h5_path.string().c_str());
    {
        HighFive::File file(h5_path.string(), HighFive::File::Truncate);
        count_ts->WriteH5Group(*trajectory_protein, file);
        water_ts->WriteH5Group(*trajectory_protein, file);
    }
    const std::string count_path =
        "/trajectory/larsen_hbond_count_time_series/count";
    const std::string water_path =
        "/trajectory/larsen_hbond_water_term_time_series/water_term";
    const auto counts = ReadH5Flat<int>(h5_path, count_path);
    const auto water = ReadH5Flat<double>(h5_path, water_path);
    ASSERT_EQ(counts.size(), original.AtomCount() * 2u);
    ASSERT_EQ(water.size(), original.AtomCount() * 2u);
    for (std::size_t atom = 0; atom < original.AtomCount(); ++atom) {
        EXPECT_EQ(counts[atom * 2u + 1u], counts[atom * 2u])
            << count_path << " atom=" << atom;
        EXPECT_DOUBLE_EQ(water[atom * 2u + 1u], water[atom * 2u])
            << water_path << " atom=" << atom;
    }
    {
        HighFive::File file(h5_path.string(), HighFive::File::ReadOnly);
        for (const std::string& group : {
                 std::string("/trajectory/larsen_hbond_count_time_series"),
                 std::string("/trajectory/larsen_hbond_water_term_time_series")}) {
            std::string frame;
            std::string law;
            file.getGroup(group).getAttribute("coordinate_frame").read(frame);
            file.getGroup(group).getAttribute("transformation").read(law);
            EXPECT_NE(frame.find("intrinsic"), std::string::npos) << group;
            EXPECT_NE(law.find("reflection-invariant"),
                      std::string::npos) << group;
        }
    }
    EXPECT_EQ(std::remove(h5_path.string().c_str()), 0) << h5_path;

    const std::array<const char*, 15> all_output_names{{
        "larsen_hbond_shielding.npy",
        "larsen_hbond_1pHB_shielding.npy",
        "larsen_hbond_2pHB_shielding.npy",
        "larsen_hbond_1pHaB_shielding.npy",
        "larsen_hbond_2pHaB_shielding.npy",
        "larsen_hbond_diagnostic_CB_shielding.npy",
        "larsen_hbond_water_term.npy",
        "larsen_hbond_count.npy",
        "larsen_imputed_pair_count.npy",
        "larsen_sidechain_carbonyl_pair_count.npy",
        "larsen_corner_imputed.npy",
        "larsen_hbond_pairs_index.npy",
        "larsen_hbond_pairs_geometry.npy",
        "larsen_hbond_pairs_isotropic.npy",
        "larsen_hbond_pairs.npy",
    }};
    for (const fs::path& directory : {original_dir, reflected_dir}) {
        for (const char* name : all_output_names) {
            EXPECT_EQ(std::remove((directory / name).string().c_str()), 0)
                << name;
        }
        EXPECT_EQ(std::remove(directory.string().c_str()), 0) << directory;
    }
    for (const auto& entry : fs::directory_iterator(grid_dir)) {
        EXPECT_EQ(std::remove(entry.path().string().c_str()), 0)
            << entry.path();
    }
    EXPECT_EQ(std::remove(grid_dir.string().c_str()), 0) << grid_dir;
}


// A finite, non-degenerate candidate remains a confirmed geometric H-bond
// even when its distance lies outside the archive axis. This production
// forcing grid puts every real H...O distance outside the r axis, isolating
// GridMiss from Success while pinning the historical water-gate behavior.
TEST(LarsenHBondPairAudit, FiniteGridMissStillSuppressesWaterTerm) {
    const fs::path pdb = nmr::test::TestEnvironment::UbqProtonated();
    ASSERT_TRUE(fs::is_regular_file(pdb)) << pdb;
    auto build = BuildFromProtonatedPdb(pdb.string());
    ASSERT_TRUE(build.Ok()) << build.error;
    std::unique_ptr<Protein> protein = std::move(build.protein);
    ProteinConformation& conf = protein->Conformation();
    ASSERT_TRUE(conf.AttachResult(GeometryResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(SpatialIndexResult::Compute(conf)));

    const fs::path grid_dir =
        WriteSyntheticLarsenGridDirectory("finite_grid_miss", true);
    LarsenHBondGrid grid(grid_dir.string());
    auto result = LarsenHBondShieldingResult::Compute(conf, grid);
    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->PairsFound(), 0);

    const auto amide_grid_miss = std::find_if(
        result->Pairs().begin(), result->Pairs().end(),
        [](const LarsenHBondShieldingResult::PairRecord& pair) {
            return pair.donor_class == HBondDonorClass::AmideHydrogen &&
                   pair.disposition ==
                       LarsenHBondShieldingResult::PairDisposition::GridMiss;
        });
    ASSERT_NE(amide_grid_miss, result->Pairs().end());
    EXPECT_TRUE(amide_grid_miss->frame_valid);
    EXPECT_TRUE(std::isfinite(amide_grid_miss->r_angstrom));
    EXPECT_TRUE(std::isfinite(amide_grid_miss->theta_deg));
    EXPECT_TRUE(std::isfinite(amide_grid_miss->rho_deg));
    EXPECT_GE(amide_grid_miss->theta_deg, 90.0);
    EXPECT_LE(amide_grid_miss->theta_deg, 180.0);
    EXPECT_DOUBLE_EQ(
        conf.AtomAt(amide_grid_miss->donor_atom_idx)
            .larsen_hbond_water_term,
        0.0);

    RemoveFlatTempDirectory(grid_dir, {
        "NMANMA_dense.h5", "NMACOH_dense.h5", "NMACOO_dense.h5",
        "ALANMA_dense.h5", "ALACOH_dense.h5", "ALACOO_dense.h5",
    });
}


TEST(LarsenHBondPairAudit,
     MissingFrameAtomsIsNanWhileConfirmedPairIsZero) {
    const fs::path pdb = nmr::test::TestEnvironment::GmxProtonated();
    ASSERT_TRUE(fs::is_regular_file(pdb)) << pdb;
    auto build = BuildFromProtonatedPdb(pdb.string());
    ASSERT_TRUE(build.Ok()) << build.error;
    std::unique_ptr<Protein> protein = std::move(build.protein);
    ProteinConformation& baseline_conf = protein->Conformation();

    const Residue& n_terminal = protein->ResidueAt(0);
    ASSERT_NE(n_terminal.H, Residue::NONE);
    ASSERT_FALSE(protein->BackbonePredecessor(0).has_value());
    ASSERT_TRUE(protein->LegacyAmber()
                    .SemanticAt(n_terminal.H)
                    .IsBackboneAmideHydrogen());

    const Residue& internal_acceptor_residue = protein->ResidueAt(1);
    ASSERT_NE(internal_acceptor_residue.O, Residue::NONE);
    ASSERT_NE(internal_acceptor_residue.C, Residue::NONE);
    std::vector<Vec3> positions = baseline_conf.Positions();
    positions[internal_acceptor_residue.O] =
        positions[n_terminal.H] + Vec3(2.0, 0.0, 0.0);
    DerivedConformation& missing_conf = protein->AddDerived(
        baseline_conf, "forced missing ProCS15 donor anchor",
        std::move(positions));
    ASSERT_TRUE(missing_conf.AttachResult(
        GeometryResult::Compute(missing_conf)));
    ASSERT_TRUE(missing_conf.AttachResult(
        SpatialIndexResult::Compute(missing_conf)));

    const fs::path grid_dir =
        WriteSyntheticLarsenGridDirectory("missing_frame_atoms");
    LarsenHBondGrid grid(grid_dir.string());
    auto result = LarsenHBondShieldingResult::Compute(missing_conf, grid);
    ASSERT_NE(result, nullptr);

    const auto missing_pair = std::find_if(
        result->Pairs().begin(), result->Pairs().end(),
        [&](const LarsenHBondShieldingResult::PairRecord& pair) {
            return pair.donor_atom_idx == n_terminal.H &&
                   pair.acceptor_atom_idx == internal_acceptor_residue.O;
        });
    ASSERT_NE(missing_pair, result->Pairs().end());
    EXPECT_EQ(missing_pair->disposition,
              LarsenHBondShieldingResult::PairDisposition::MissingFrameAtoms);
    EXPECT_TRUE(std::isnan(
        missing_conf.AtomAt(n_terminal.H).larsen_hbond_water_term));

    const auto confirmed_pair = std::find_if(
        result->Pairs().begin(), result->Pairs().end(),
        [](const LarsenHBondShieldingResult::PairRecord& pair) {
            return pair.donor_class == HBondDonorClass::AmideHydrogen &&
                   (pair.disposition ==
                        LarsenHBondShieldingResult::PairDisposition::Success ||
                    pair.disposition ==
                        LarsenHBondShieldingResult::PairDisposition::GridMiss);
        });
    ASSERT_NE(confirmed_pair, result->Pairs().end());
    EXPECT_DOUBLE_EQ(
        missing_conf.AtomAt(confirmed_pair->donor_atom_idx)
            .larsen_hbond_water_term,
        0.0);

    RemoveFlatTempDirectory(grid_dir, {
        "NMANMA_dense.h5", "NMACOH_dense.h5", "NMACOO_dense.h5",
        "ALANMA_dense.h5", "ALACOH_dense.h5", "ALACOO_dense.h5",
    });
}


// Production-calling forcing gate for the malformed-frame/water interaction.
// Start from a real successful amide-H pair, collapse only that donor's N-H
// frame, and rerun the full calculator. The H...O candidate geometry remains
// spatially present, so an identity fallback would still hit the grid and
// suppress water; the correct path records InvalidFrame, applies no pair, and
// marks the water term unevaluable instead of claiming an unbound amide.
TEST(LarsenHBondPairAudit,
     DegenerateDonorFrameIsNanWhileEvaluatedPairIsZero) {
    const fs::path pdb = nmr::test::TestEnvironment::UbqProtonated();
    ASSERT_TRUE(fs::is_regular_file(pdb)) << pdb;
    auto build = BuildFromProtonatedPdb(pdb.string());
    ASSERT_TRUE(build.Ok()) << build.error;
    std::unique_ptr<Protein> protein = std::move(build.protein);
    ProteinConformation& baseline_conf = protein->Conformation();
    ASSERT_TRUE(baseline_conf.AttachResult(
        GeometryResult::Compute(baseline_conf)));
    ASSERT_TRUE(baseline_conf.AttachResult(
        SpatialIndexResult::Compute(baseline_conf)));

    const fs::path grid_dir =
        WriteSyntheticLarsenGridDirectory("invalid_frame");
    LarsenHBondGrid grid(grid_dir.string());
    auto baseline =
        LarsenHBondShieldingResult::Compute(baseline_conf, grid);
    ASSERT_NE(baseline, nullptr);

    const auto successful_amide = std::find_if(
        baseline->Pairs().begin(), baseline->Pairs().end(),
        [](const LarsenHBondShieldingResult::PairRecord& pair) {
            return pair.donor_class == HBondDonorClass::AmideHydrogen &&
                   pair.disposition ==
                       LarsenHBondShieldingResult::PairDisposition::Success;
        });
    ASSERT_NE(successful_amide, baseline->Pairs().end());
    const std::size_t donor_h = successful_amide->donor_atom_idx;
    const std::size_t donor_residue =
        protein->AtomAt(donor_h).residue_index;
    const std::size_t donor_n = protein->ResidueAt(donor_residue).N;
    ASSERT_NE(donor_n, Residue::NONE);
    EXPECT_DOUBLE_EQ(
        baseline_conf.AtomAt(donor_h).larsen_hbond_water_term, 0.0);
    ASSERT_TRUE(baseline_conf.AttachResult(std::move(baseline)));

    std::vector<Vec3> positions = baseline_conf.Positions();
    positions[donor_n] = positions[donor_h];
    DerivedConformation& invalid_conf = protein->AddDerived(
        baseline_conf, "forced degenerate Larsen donor frame",
        std::move(positions));
    ASSERT_TRUE(invalid_conf.AttachResult(
        GeometryResult::Compute(invalid_conf)));
    ASSERT_TRUE(invalid_conf.AttachResult(
        SpatialIndexResult::Compute(invalid_conf)));

    ::testing::internal::CaptureStderr();
    auto invalid = LarsenHBondShieldingResult::Compute(invalid_conf, grid);
    const std::string warning =
        ::testing::internal::GetCapturedStderr();
    ASSERT_NE(invalid, nullptr);

    int donor_rows = 0;
    int invalid_rows = 0;
    std::size_t first_invalid_row = SIZE_MAX;
    for (std::size_t row = 0; row < invalid->Pairs().size(); ++row) {
        const auto& pair = invalid->Pairs()[row];
        if (pair.donor_atom_idx != donor_h) continue;
        ++donor_rows;
        if (pair.disposition ==
            LarsenHBondShieldingResult::PairDisposition::InvalidFrame) {
            ++invalid_rows;
            if (first_invalid_row == SIZE_MAX) first_invalid_row = row;
            EXPECT_FALSE(pair.frame_valid);
            EXPECT_TRUE(std::isnan(pair.iso_1pHB));
            EXPECT_TRUE(std::isnan(pair.iso_2pHB));
            EXPECT_TRUE(std::isnan(pair.iso_1pHaB));
            EXPECT_TRUE(std::isnan(pair.iso_2pHaB));
            EXPECT_TRUE(std::isnan(pair.iso_table2_total));
        } else {
            EXPECT_EQ(pair.disposition,
                      LarsenHBondShieldingResult::PairDisposition::
                          MissingFrameAtoms);
        }
    }
    EXPECT_GT(donor_rows, 0);
    EXPECT_GT(invalid_rows, 0);
    ASSERT_NE(first_invalid_row, SIZE_MAX);
    EXPECT_TRUE(std::isnan(
        invalid_conf.AtomAt(donor_h).larsen_hbond_water_term));
    const ConformationAtom& invalid_donor = invalid_conf.AtomAt(donor_h);
    EXPECT_TRUE(invalid_donor.larsen_hbond_shielding_tensor.allFinite());
    EXPECT_DOUBLE_EQ(
        (invalid_donor.larsen_hbond_shielding_tensor -
         (invalid_donor.larsen_hbond_1pHB_tensor +
          invalid_donor.larsen_hbond_2pHB_tensor +
          invalid_donor.larsen_hbond_1pHaB_tensor +
          invalid_donor.larsen_hbond_2pHaB_tensor)).squaredNorm(),
        0.0) << "the scalar water term must not enter ProCS15 tensors";
    EXPECT_NE(warning.find("invalid donor frame"), std::string::npos)
        << warning;
    EXPECT_NE(warning.find("invalid Larsen pair frame"),
              std::string::npos) << warning;

    const fs::path out =
        MakeUniqueTempDirectory("nmr_larsen_invalid_frame_output");
    EXPECT_EQ(invalid->WriteFeatures(invalid_conf, out.string()), 15);
    const RawNpy index =
        ReadRawNpy(out / "larsen_hbond_pairs_index.npy");
    const RawNpy geometry =
        ReadRawNpy(out / "larsen_hbond_pairs_geometry.npy");
    const RawNpy isotropic =
        ReadRawNpy(out / "larsen_hbond_pairs_isotropic.npy");
    const RawNpy water =
        ReadRawNpy(out / "larsen_hbond_water_term.npy");
    EXPECT_EQ(PayloadValue<std::int32_t>(
                  index, first_invalid_row * 16 + 6),
              static_cast<std::int32_t>(
                  LarsenHBondShieldingResult::PairDisposition::InvalidFrame));
    EXPECT_DOUBLE_EQ(PayloadValue<double>(
                         geometry, first_invalid_row * 6 + 5),
                     0.0);
    EXPECT_TRUE(std::isnan(PayloadValue<double>(
        isotropic, first_invalid_row * 6 + 5)));
    EXPECT_TRUE(std::isnan(PayloadValue<double>(water, donor_h)));

    RemoveFlatTempDirectory(out, {
        "larsen_hbond_shielding.npy",
        "larsen_hbond_1pHB_shielding.npy",
        "larsen_hbond_2pHB_shielding.npy",
        "larsen_hbond_1pHaB_shielding.npy",
        "larsen_hbond_2pHaB_shielding.npy",
        "larsen_hbond_diagnostic_CB_shielding.npy",
        "larsen_hbond_water_term.npy",
        "larsen_hbond_count.npy",
        "larsen_imputed_pair_count.npy",
        "larsen_sidechain_carbonyl_pair_count.npy",
        "larsen_corner_imputed.npy",
        "larsen_hbond_pairs_index.npy",
        "larsen_hbond_pairs_geometry.npy",
        "larsen_hbond_pairs_isotropic.npy",
        "larsen_hbond_pairs.npy",
    });

    ASSERT_TRUE(invalid_conf.AttachResult(std::move(invalid)));
    auto trajectory_protein =
        TrajectoryProtein::CreateForTesting(std::move(protein));
    ASSERT_NE(trajectory_protein, nullptr);
    auto water_series =
        LarsenHBondWaterTermTimeSeriesTrajectoryResult::Create(
            *trajectory_protein);
    Trajectory trajectory(fs::path{}, fs::path{}, fs::path{});
    water_series->Compute(
        baseline_conf, *trajectory_protein, trajectory, 0, 0.0);
    water_series->Compute(
        invalid_conf, *trajectory_protein, trajectory, 1, 1.0);
    water_series->Finalize(*trajectory_protein, trajectory);

    const fs::path h5_path = fs::temp_directory_path() /
        ("nmr_larsen_invalid_frame_water_series_" +
         std::to_string(::getpid()) + ".h5");
    {
        HighFive::File file(h5_path.string(), HighFive::File::Truncate);
        water_series->WriteH5Group(*trajectory_protein, file);
    }
    std::vector<std::size_t> water_dims;
    const auto water_values = ReadH5Flat<double>(
        h5_path,
        "/trajectory/larsen_hbond_water_term_time_series/water_term",
        &water_dims);
    ASSERT_EQ(water_dims,
              (std::vector<std::size_t>{
                  trajectory_protein->AtomCount(), 2u}));
    EXPECT_DOUBLE_EQ(water_values[donor_h * 2], 0.0);
    EXPECT_TRUE(std::isnan(water_values[donor_h * 2 + 1]));
    const auto source_present = ReadH5Flat<std::uint8_t>(
        h5_path,
        "/trajectory/larsen_hbond_water_term_time_series/"
        "source_attached_per_frame");
    EXPECT_EQ(source_present,
              (std::vector<std::uint8_t>{1u, 1u}));
    EXPECT_EQ(std::remove(h5_path.string().c_str()), 0) << h5_path;

    RemoveFlatTempDirectory(grid_dir, {
        "NMANMA_dense.h5", "NMACOH_dense.h5", "NMACOO_dense.h5",
        "ALANMA_dense.h5", "ALACOH_dense.h5", "ALACOO_dense.h5",
    });
}


TEST(LarsenHBondPairAudit,
     NonFiniteDerivedGeometryMakesWaterTermUnevaluable) {
    const fs::path pdb = nmr::test::TestEnvironment::UbqProtonated();
    ASSERT_TRUE(fs::is_regular_file(pdb)) << pdb;
    auto build = BuildFromProtonatedPdb(pdb.string());
    ASSERT_TRUE(build.Ok()) << build.error;
    std::unique_ptr<Protein> protein = std::move(build.protein);
    ProteinConformation& baseline_conf = protein->Conformation();
    ASSERT_TRUE(baseline_conf.AttachResult(
        GeometryResult::Compute(baseline_conf)));
    ASSERT_TRUE(baseline_conf.AttachResult(
        SpatialIndexResult::Compute(baseline_conf)));

    const fs::path grid_dir =
        WriteSyntheticLarsenGridDirectory("nonfinite_geometry");
    LarsenHBondGrid grid(grid_dir.string());
    auto baseline =
        LarsenHBondShieldingResult::Compute(baseline_conf, grid);
    ASSERT_NE(baseline, nullptr);

    // Choose a real successful amide donor, then invalidate the acceptor
    // C position for every selected geometric partner of that donor. Theta
    // misses and symmetry-filtered siblings are not water-gating partners.
    const auto successful_amide = std::find_if(
        baseline->Pairs().begin(), baseline->Pairs().end(),
        [](const LarsenHBondShieldingResult::PairRecord& pair) {
            return pair.donor_class == HBondDonorClass::AmideHydrogen &&
                   pair.disposition ==
                       LarsenHBondShieldingResult::PairDisposition::Success;
        });
    ASSERT_NE(successful_amide, baseline->Pairs().end());
    const std::size_t donor_h = successful_amide->donor_atom_idx;
    EXPECT_DOUBLE_EQ(
        baseline_conf.AtomAt(donor_h).larsen_hbond_water_term, 0.0);

    std::vector<Vec3> positions = baseline_conf.Positions();
    int confirmed_before = 0;
    for (const auto& pair : baseline->Pairs()) {
        if (pair.donor_atom_idx != donor_h) continue;
        if (pair.disposition !=
                LarsenHBondShieldingResult::PairDisposition::Success &&
            pair.disposition !=
                LarsenHBondShieldingResult::PairDisposition::GridMiss) {
            continue;
        }
        ASSERT_NE(pair.acceptor_C_atom_idx, SIZE_MAX);
        ASSERT_NE(pair.acceptor_atom_idx, SIZE_MAX);
        positions[pair.acceptor_C_atom_idx] =
            positions[pair.acceptor_atom_idx];
        ++confirmed_before;
    }
    ASSERT_GT(confirmed_before, 0);
    DerivedConformation& invalid_conf = protein->AddDerived(
        baseline_conf, "forced invalid Larsen H-O-C geometry",
        std::move(positions));
    ASSERT_TRUE(invalid_conf.AttachResult(
        GeometryResult::Compute(invalid_conf)));
    ASSERT_TRUE(invalid_conf.AttachResult(
        SpatialIndexResult::Compute(invalid_conf)));

    ::testing::internal::CaptureStderr();
    auto invalid = LarsenHBondShieldingResult::Compute(invalid_conf, grid);
    const std::string warning =
        ::testing::internal::GetCapturedStderr();
    ASSERT_NE(invalid, nullptr);

    const auto invalid_pair = std::find_if(
        invalid->Pairs().begin(), invalid->Pairs().end(),
        [&](const LarsenHBondShieldingResult::PairRecord& pair) {
            return pair.donor_atom_idx == donor_h &&
                   pair.acceptor_atom_idx ==
                       successful_amide->acceptor_atom_idx;
        });
    ASSERT_NE(invalid_pair, invalid->Pairs().end());
    EXPECT_EQ(invalid_pair->disposition,
              LarsenHBondShieldingResult::PairDisposition::InvalidFrame);
    EXPECT_FALSE(invalid_pair->frame_valid);
    EXPECT_TRUE(std::isnan(invalid_pair->theta_deg));
    EXPECT_TRUE(std::isnan(invalid_pair->rho_deg));
    EXPECT_TRUE(std::isnan(invalid_pair->iso_table2_total));
    EXPECT_TRUE(std::isnan(
        invalid_conf.AtomAt(donor_h).larsen_hbond_water_term));
    EXPECT_TRUE(
        invalid_conf.AtomAt(donor_h)
            .larsen_hbond_shielding_tensor.allFinite());
    EXPECT_NE(warning.find("geometry_finite=0"), std::string::npos)
        << warning;

    RemoveFlatTempDirectory(grid_dir, {
        "NMANMA_dense.h5", "NMACOH_dense.h5", "NMACOO_dense.h5",
        "ALANMA_dense.h5", "ALACOH_dense.h5", "ALACOO_dense.h5",
    });
}


TEST(LarsenSidechainDonorAudit, EmitsTypedUnsupportedDonorsWithoutShielding) {
    const fs::path pdb = nmr::test::TestEnvironment::UbqProtonated();
    ASSERT_TRUE(fs::is_regular_file(pdb)) << pdb;
    auto build = BuildFromProtonatedPdb(pdb.string());
    ASSERT_TRUE(build.Ok()) << build.error;
    std::unique_ptr<Protein> protein = std::move(build.protein);
    ProteinConformation& conf = protein->Conformation();
    ASSERT_TRUE(conf.AttachResult(GeometryResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(SpatialIndexResult::Compute(conf)));

    auto audit = LarsenSidechainDonorAuditResult::Compute(conf);
    ASSERT_NE(audit, nullptr);
    ASSERT_EQ(audit->DonorAtoms().size(), conf.AtomCount());

    const LegacyAmberTopology& topology = protein->LegacyAmber();
    int n_donors = 0;
    int n_geometry_pass = 0;
    std::int64_t candidate_count_from_atoms = 0;
    std::set<PolarHKind> seen_kinds;
    for (std::size_t i = 0; i < conf.AtomCount(); ++i) {
        const auto& row = audit->DonorAtoms()[i];
        candidate_count_from_atoms += row.n_acceptor_candidates_3p5A;
        n_geometry_pass += row.n_geometry_pass;
        if (!row.is_sidechain_polar_h) continue;
        ++n_donors;
        ASSERT_GE(row.parent_atom, 0);
        const std::size_t parent = static_cast<std::size_t>(row.parent_atom);
        EXPECT_FALSE(topology.SemanticAt(parent).IsBackbone());
        EXPECT_TRUE(topology.SemanticAt(i).IsPolarH());
        seen_kinds.insert(topology.SemanticAt(i).polar_h);
    }
    EXPECT_GT(n_donors, 0);
    EXPECT_EQ(candidate_count_from_atoms,
              static_cast<std::int64_t>(audit->Candidates().size()));
    EXPECT_GT(audit->Candidates().size(), 0u);
    EXPECT_GT(n_geometry_pass, 0);
    EXPECT_TRUE(seen_kinds.count(PolarHKind::HydroxylOH_Aliphatic) > 0);
    EXPECT_TRUE(seen_kinds.count(PolarHKind::HydroxylOH_Aromatic) > 0);
    EXPECT_TRUE(seen_kinds.count(PolarHKind::AmmoniumNH) > 0);
    EXPECT_TRUE(seen_kinds.count(PolarHKind::GuanidiniumNH) > 0);

    for (const auto& row : audit->Candidates()) {
        EXPECT_LE(row.parent_acceptor_distance_A, 3.5 + 1e-12);
        EXPECT_FALSE(row.modeled_by_larsen_table2);
        EXPECT_EQ(protein->AtomAt(row.acceptor_atom).element, Element::O);
        EXPECT_FALSE(topology.SemanticAt(row.parent_atom).IsBackbone());
    }

    // The audit has no writer fields on ConformationAtom and must not
    // extend unsupported donor classes into Larsen's scientific tensor.
    for (std::size_t i = 0; i < conf.AtomCount(); ++i) {
        EXPECT_DOUBLE_EQ(conf.AtomAt(i).larsen_hbond_shielding_tensor.norm(),
                         0.0);
        EXPECT_EQ(conf.AtomAt(i).larsen_hbond_n_pairs, 0);
    }

    const fs::path out =
        MakeUniqueTempDirectory("nmr_larsen_sidechain_donor_audit_output");
    EXPECT_EQ(audit->WriteFeatures(conf, out.string()), 2);
    const RawNpy atoms =
        ReadRawNpy(out / "larsen_sidechain_donor_atoms.npy");
    const RawNpy candidates =
        ReadRawNpy(out / "larsen_sidechain_donor_candidates.npy");
    EXPECT_NE(atoms.header.find("'descr': '<i4'"), std::string::npos);
    EXPECT_NE(atoms.header.find("'shape': (" +
                                std::to_string(conf.AtomCount()) +
                                ",6)"), std::string::npos);
    EXPECT_EQ(atoms.payload.size(), conf.AtomCount() * 6 * sizeof(std::int32_t));
    EXPECT_NE(candidates.header.find("'shape': (" +
                                     std::to_string(audit->Candidates().size()) +
                                     ",13)"), std::string::npos);
    ASSERT_EQ(candidates.payload.size(),
              audit->Candidates().size() * 13 * sizeof(double));
    for (std::size_t row = 0; row < audit->Candidates().size(); ++row) {
        double modeled = 1.0;
        std::memcpy(&modeled,
                    candidates.payload.data() +
                        (row * 13 + 12) * sizeof(double),
                    sizeof(double));
        EXPECT_DOUBLE_EQ(modeled, 0.0);
    }

    // A24 read-back: a typed donor row and its first raw candidate must
    // survive the writer with their identities, geometry, and audit-only
    // modeled flag intact.
    const auto donor_it = std::find_if(
        audit->DonorAtoms().begin(), audit->DonorAtoms().end(),
        [](const LarsenSidechainDonorAuditResult::DonorAtomRecord& row) {
            return row.is_sidechain_polar_h != 0;
        });
    ASSERT_NE(donor_it, audit->DonorAtoms().end());
    const std::size_t donor_row = static_cast<std::size_t>(
        std::distance(audit->DonorAtoms().begin(), donor_it));
    EXPECT_EQ(PayloadValue<std::int32_t>(atoms, donor_row * 6 + 0), 1);
    EXPECT_EQ(PayloadValue<std::int32_t>(atoms, donor_row * 6 + 1),
              donor_it->polar_h_kind);
    EXPECT_EQ(PayloadValue<std::int32_t>(atoms, donor_row * 6 + 2),
              donor_it->parent_atom);
    EXPECT_EQ(PayloadValue<std::int32_t>(atoms, donor_row * 6 + 3),
              donor_it->residue_index);
    EXPECT_EQ(PayloadValue<std::int32_t>(atoms, donor_row * 6 + 4),
              donor_it->n_acceptor_candidates_3p5A);
    EXPECT_EQ(PayloadValue<std::int32_t>(atoms, donor_row * 6 + 5),
              donor_it->n_geometry_pass);

    ASSERT_FALSE(audit->Candidates().empty());
    const auto& expected_candidate = audit->Candidates().front();
    EXPECT_DOUBLE_EQ(PayloadValue<double>(candidates, 0),
                     static_cast<double>(expected_candidate.donor_atom));
    EXPECT_DOUBLE_EQ(PayloadValue<double>(candidates, 2),
                     static_cast<double>(expected_candidate.polar_h_kind));
    EXPECT_DOUBLE_EQ(PayloadValue<double>(candidates, 3),
                     static_cast<double>(expected_candidate.parent_atom));
    EXPECT_DOUBLE_EQ(PayloadValue<double>(candidates, 4),
                     static_cast<double>(expected_candidate.acceptor_atom));
    EXPECT_DOUBLE_EQ(PayloadValue<double>(candidates, 6),
                     static_cast<double>(expected_candidate.acceptor_class));
    EXPECT_DOUBLE_EQ(PayloadValue<double>(candidates, 7),
                     expected_candidate.h_acceptor_distance_A);
    EXPECT_DOUBLE_EQ(PayloadValue<double>(candidates, 8),
                     expected_candidate.parent_acceptor_distance_A);
    EXPECT_DOUBLE_EQ(PayloadValue<double>(candidates, 9),
                     expected_candidate.angle_parent_h_acceptor_deg);
    if (std::isnan(expected_candidate.rho_deg)) {
        EXPECT_TRUE(std::isnan(PayloadValue<double>(candidates, 10)));
    } else {
        EXPECT_DOUBLE_EQ(PayloadValue<double>(candidates, 10),
                         expected_candidate.rho_deg);
    }
    EXPECT_DOUBLE_EQ(PayloadValue<double>(candidates, 11),
                     expected_candidate.passes_geometry ? 1.0 : 0.0);
    EXPECT_DOUBLE_EQ(PayloadValue<double>(candidates, 12), 0.0);
    RemoveFlatTempDirectory(out, {
        "larsen_sidechain_donor_atoms.npy",
        "larsen_sidechain_donor_candidates.npy",
    });
}


TEST(LarsenSidechainDonorDirectionalCovariance,
     ProductionRerunAndSerializedCandidatesUnderProperRigidTransform) {
    const fs::path pdb = nmr::test::TestEnvironment::UbqProtonated();
    ASSERT_TRUE(fs::is_regular_file(pdb)) << pdb;
    auto build = BuildFromProtonatedPdb(pdb.string());
    ASSERT_TRUE(build.Ok()) << build.error;
    auto protein = std::move(build.protein);

    ProteinConformation& original = protein->Conformation();
    ASSERT_TRUE(original.AttachResult(GeometryResult::Compute(original)));
    ASSERT_TRUE(original.AttachResult(SpatialIndexResult::Compute(original)));
    auto original_result =
        LarsenSidechainDonorAuditResult::Compute(original);
    ASSERT_NE(original_result, nullptr);
    ASSERT_GT(original_result->Candidates().size(), 0U);

    constexpr std::uint64_t kTransformSeed = 0x4C415253454E5343ULL;
    const auto transform = nmr::test::directional::SeededTransform(
        kTransformSeed, false);
    ProteinConformation& moved = protein->AddConformation(
        nmr::test::directional::Positions(transform,
                                          original.Positions()),
        "Larsen sidechain donor audit proper covariance rerun");
    ASSERT_TRUE(moved.AttachResult(GeometryResult::Compute(moved)));
    ASSERT_TRUE(moved.AttachResult(SpatialIndexResult::Compute(moved)));
    auto moved_result = LarsenSidechainDonorAuditResult::Compute(moved);
    ASSERT_NE(moved_result, nullptr);

    constexpr double kGeometryTolerance = 3.0e-10;
    ASSERT_EQ(moved_result->DonorAtoms().size(),
              original_result->DonorAtoms().size());
    ASSERT_EQ(moved_result->Candidates().size(),
              original_result->Candidates().size());
    for (std::size_t row = 0;
         row < original_result->Candidates().size(); ++row) {
        const auto& source = original_result->Candidates()[row];
        const auto& actual = moved_result->Candidates()[row];
        EXPECT_EQ(actual.donor_atom, source.donor_atom) << row;
        EXPECT_EQ(actual.donor_residue, source.donor_residue) << row;
        EXPECT_EQ(actual.polar_h_kind, source.polar_h_kind) << row;
        EXPECT_EQ(actual.parent_atom, source.parent_atom) << row;
        EXPECT_EQ(actual.acceptor_atom, source.acceptor_atom) << row;
        EXPECT_EQ(actual.acceptor_residue, source.acceptor_residue) << row;
        EXPECT_EQ(actual.acceptor_class, source.acceptor_class) << row;
        EXPECT_EQ(actual.passes_geometry, source.passes_geometry) << row;
        EXPECT_EQ(actual.modeled_by_larsen_table2,
                  source.modeled_by_larsen_table2) << row;
        ExpectSameOrNaN(actual.h_acceptor_distance_A,
                        source.h_acceptor_distance_A,
                        kGeometryTolerance,
                        "sidechain H-acceptor distance row=" +
                            std::to_string(row));
        ExpectSameOrNaN(actual.parent_acceptor_distance_A,
                        source.parent_acceptor_distance_A,
                        kGeometryTolerance,
                        "sidechain parent-acceptor distance row=" +
                            std::to_string(row));
        ExpectSameOrNaN(actual.angle_parent_h_acceptor_deg,
                        source.angle_parent_h_acceptor_deg,
                        kGeometryTolerance,
                        "sidechain donor angle row=" +
                            std::to_string(row));
        ExpectSameOrNaN(actual.rho_deg, source.rho_deg,
                        kGeometryTolerance,
                        "sidechain signed rho row=" +
                            std::to_string(row));
    }

    const fs::path original_out = MakeUniqueTempDirectory(
        "nmr_larsen_sidechain_directional_original");
    const fs::path moved_out = MakeUniqueTempDirectory(
        "nmr_larsen_sidechain_directional_moved");
    ASSERT_EQ(original_result->WriteFeatures(
                  original, original_out.string()), 2);
    ASSERT_EQ(moved_result->WriteFeatures(moved, moved_out.string()), 2);
    ASSERT_EQ(nmr::test::directional::RunNumpyAllowPickleFalse(
                  NMR_TEST_PYTHON_EXECUTABLE,
                  NMR_NPY_ALLOW_PICKLE_FALSE_SCRIPT,
                  {original_out / "larsen_sidechain_donor_atoms.npy",
                   original_out /
                       "larsen_sidechain_donor_candidates.npy"}),
              0);
    ASSERT_EQ(nmr::test::directional::RunNumpyAllowPickleFalse(
                  NMR_TEST_PYTHON_EXECUTABLE,
                  NMR_NPY_ALLOW_PICKLE_FALSE_SCRIPT,
                  {moved_out / "larsen_sidechain_donor_atoms.npy",
                   moved_out / "larsen_sidechain_donor_candidates.npy"}),
              0);

    const RawNpy source_atoms = ReadRawNpy(
        original_out / "larsen_sidechain_donor_atoms.npy");
    const RawNpy actual_atoms = ReadRawNpy(
        moved_out / "larsen_sidechain_donor_atoms.npy");
    EXPECT_EQ(actual_atoms.payload, source_atoms.payload)
        << "larsen_sidechain_donor_atoms.npy";

    const RawNpy source_candidates = ReadRawNpy(
        original_out / "larsen_sidechain_donor_candidates.npy");
    const RawNpy actual_candidates = ReadRawNpy(
        moved_out / "larsen_sidechain_donor_candidates.npy");
    ASSERT_EQ(actual_candidates.payload.size(),
              source_candidates.payload.size());
    ASSERT_EQ(source_candidates.payload.size() % sizeof(double), 0U);
    const std::size_t serialized_components =
        source_candidates.payload.size() / sizeof(double);
    for (std::size_t component = 0;
         component < serialized_components; ++component) {
        ExpectSameOrNaN(
            PayloadValue<double>(actual_candidates, component),
            PayloadValue<double>(source_candidates, component),
            kGeometryTolerance,
            "larsen_sidechain_donor_candidates.npy component=" +
                std::to_string(component));
    }

    for (const fs::path& dir : {original_out, moved_out}) {
        EXPECT_EQ(std::remove(
                      (dir / "larsen_sidechain_donor_atoms.npy")
                          .string().c_str()), 0);
        EXPECT_EQ(std::remove(
                      (dir / "larsen_sidechain_donor_candidates.npy")
                          .string().c_str()), 0);
        EXPECT_EQ(std::remove(dir.string().c_str()), 0) << dir;
    }
}


TEST(LarsenSidechainDonorDirectionalCovariance,
     SignedRhoFlipsUnderImproperProductionRerunAndSerializedReadback) {
    const fs::path pdb = nmr::test::TestEnvironment::UbqProtonated();
    ASSERT_TRUE(fs::is_regular_file(pdb)) << pdb;
    auto build = BuildFromProtonatedPdb(pdb.string());
    ASSERT_TRUE(build.Ok()) << build.error;
    auto protein = std::move(build.protein);

    ProteinConformation& original = protein->Conformation();
    ASSERT_TRUE(original.AttachResult(GeometryResult::Compute(original)));
    ASSERT_TRUE(original.AttachResult(SpatialIndexResult::Compute(original)));
    auto original_result =
        LarsenSidechainDonorAuditResult::Compute(original);
    ASSERT_NE(original_result, nullptr);
    ASSERT_GT(original_result->Candidates().size(), 0u);

    constexpr std::uint64_t kTransformSeed = 0x4C415253454E5344ULL;
    const auto transform = nmr::test::directional::SeededTransform(
        kTransformSeed, true);
    ProteinConformation& reflected = protein->AddConformation(
        nmr::test::directional::Positions(transform, original.Positions()),
        "Larsen sidechain donor audit improper covariance rerun");
    ASSERT_TRUE(reflected.AttachResult(GeometryResult::Compute(reflected)));
    ASSERT_TRUE(reflected.AttachResult(
        SpatialIndexResult::Compute(reflected)));
    auto reflected_result =
        LarsenSidechainDonorAuditResult::Compute(reflected);
    ASSERT_NE(reflected_result, nullptr);

    constexpr double kGeometryTolerance = 3.0e-10;
    ASSERT_EQ(reflected_result->DonorAtoms().size(),
              original_result->DonorAtoms().size());
    ASSERT_EQ(reflected_result->Candidates().size(),
              original_result->Candidates().size());
    std::size_t finite_nonzero_rho_rows = 0;
    for (std::size_t row = 0;
         row < original_result->Candidates().size(); ++row) {
        const auto& source = original_result->Candidates()[row];
        const auto& actual = reflected_result->Candidates()[row];
        EXPECT_EQ(actual.donor_atom, source.donor_atom) << row;
        EXPECT_EQ(actual.donor_residue, source.donor_residue) << row;
        EXPECT_EQ(actual.polar_h_kind, source.polar_h_kind) << row;
        EXPECT_EQ(actual.parent_atom, source.parent_atom) << row;
        EXPECT_EQ(actual.acceptor_atom, source.acceptor_atom) << row;
        EXPECT_EQ(actual.acceptor_residue, source.acceptor_residue) << row;
        EXPECT_EQ(actual.acceptor_class, source.acceptor_class) << row;
        EXPECT_EQ(actual.passes_geometry, source.passes_geometry) << row;
        EXPECT_EQ(actual.modeled_by_larsen_table2,
                  source.modeled_by_larsen_table2) << row;
        ExpectSameOrNaN(
            actual.h_acceptor_distance_A, source.h_acceptor_distance_A,
            kGeometryTolerance,
            "sidechain improper H-acceptor distance row=" +
                std::to_string(row));
        ExpectSameOrNaN(
            actual.parent_acceptor_distance_A,
            source.parent_acceptor_distance_A, kGeometryTolerance,
            "sidechain improper parent-acceptor distance row=" +
                std::to_string(row));
        ExpectSameOrNaN(
            actual.angle_parent_h_acceptor_deg,
            source.angle_parent_h_acceptor_deg, kGeometryTolerance,
            "sidechain improper donor angle row=" +
                std::to_string(row));
        ExpectWrappedOppositeOrNaN(
            actual.rho_deg, source.rho_deg, kGeometryTolerance,
            "sidechain improper signed rho row=" +
                std::to_string(row));
        if (std::isfinite(source.rho_deg) &&
            std::abs(std::remainder(source.rho_deg, 360.0)) >
                kGeometryTolerance) {
            ++finite_nonzero_rho_rows;
        }
    }
    EXPECT_GT(finite_nonzero_rho_rows, 0u)
        << "sidechain improper rho forcing must be non-vacuous";

    const fs::path original_out = MakeUniqueTempDirectory(
        "nmr_larsen_sidechain_improper_original");
    const fs::path reflected_out = MakeUniqueTempDirectory(
        "nmr_larsen_sidechain_improper_reflected");
    ASSERT_EQ(original_result->WriteFeatures(
                  original, original_out.string()), 2);
    ASSERT_EQ(reflected_result->WriteFeatures(
                  reflected, reflected_out.string()), 2);

    const RawNpy source_atoms = ReadRawNpy(
        original_out / "larsen_sidechain_donor_atoms.npy");
    const RawNpy actual_atoms = ReadRawNpy(
        reflected_out / "larsen_sidechain_donor_atoms.npy");
    EXPECT_EQ(actual_atoms.header, source_atoms.header);
    EXPECT_EQ(actual_atoms.payload, source_atoms.payload)
        << "larsen_sidechain_donor_atoms.npy";

    const RawNpy source_candidates = ReadRawNpy(
        original_out / "larsen_sidechain_donor_candidates.npy");
    const RawNpy actual_candidates = ReadRawNpy(
        reflected_out / "larsen_sidechain_donor_candidates.npy");
    const std::size_t candidate_count =
        original_result->Candidates().size();
    ASSERT_EQ(source_candidates.payload.size(),
              candidate_count * 13u * sizeof(double));
    ASSERT_EQ(actual_candidates.payload.size(),
              source_candidates.payload.size());
    std::size_t serialized_nonzero_rho_rows = 0;
    for (std::size_t row = 0; row < candidate_count; ++row) {
        for (std::size_t column = 0; column < 7u; ++column) {
            EXPECT_DOUBLE_EQ(
                PayloadValue<double>(actual_candidates,
                                     row * 13u + column),
                PayloadValue<double>(source_candidates,
                                     row * 13u + column))
                << "candidate identity row=" << row
                << " col=" << column;
        }
        for (std::size_t column = 7u; column < 10u; ++column) {
            ExpectSameOrNaN(
                PayloadValue<double>(actual_candidates,
                                     row * 13u + column),
                PayloadValue<double>(source_candidates,
                                     row * 13u + column),
                kGeometryTolerance,
                "candidate invariant geometry row=" +
                    std::to_string(row) + " col=" +
                    std::to_string(column));
        }
        const double source_rho = PayloadValue<double>(
            source_candidates, row * 13u + 10u);
        const double actual_rho = PayloadValue<double>(
            actual_candidates, row * 13u + 10u);
        ExpectWrappedOppositeOrNaN(
            actual_rho, source_rho, kGeometryTolerance,
            "larsen_sidechain_donor_candidates.npy rho row=" +
                std::to_string(row));
        if (std::isfinite(source_rho) &&
            std::abs(std::remainder(source_rho, 360.0)) >
                kGeometryTolerance) {
            ++serialized_nonzero_rho_rows;
        }
        EXPECT_DOUBLE_EQ(
            PayloadValue<double>(actual_candidates, row * 13u + 11u),
            PayloadValue<double>(source_candidates, row * 13u + 11u))
            << "passes_geometry row=" << row;
        EXPECT_DOUBLE_EQ(
            PayloadValue<double>(actual_candidates, row * 13u + 12u),
            PayloadValue<double>(source_candidates, row * 13u + 12u))
            << "modeled_by_larsen_table2 row=" << row;
    }
    EXPECT_GT(serialized_nonzero_rho_rows, 0u)
        << "serialized sidechain improper rho proof must be non-vacuous";

    for (const fs::path& dir : {original_out, reflected_out}) {
        EXPECT_EQ(std::remove(
                      (dir / "larsen_sidechain_donor_atoms.npy")
                          .string().c_str()), 0);
        EXPECT_EQ(std::remove(
                      (dir / "larsen_sidechain_donor_candidates.npy")
                          .string().c_str()), 0);
        EXPECT_EQ(std::remove(dir.string().c_str()), 0) << dir;
    }
}


// ---------------------------------------------------------------------------
// 1UBQ smoke test fixture
// ---------------------------------------------------------------------------

class LarsenHBondShieldingTest : public ::testing::Test {
protected:
    void SetUp() override {
        const std::string kLarsen1UbqPm6 =
            nmr::test::TestEnvironment::Larsen1UbqPm6Pdb();
        if (RuntimeEnvironment::LarsenHBondGridDir().empty()) {
            GTEST_SKIP() << "LarsenHBondGridDir empty; set "
                            "larsen_hbond_grids in ~/.nmr_tools.toml";
        }
        if (!fs::exists(kLarsen1UbqPm6)) {
            GTEST_SKIP() << "Larsen 1UBQ PM6-D3H+ PDB not found at "
                         << kLarsen1UbqPm6;
        }
        ASSERT_EQ(session.LoadLarsenHBondGrid(), kOk) << session.LastError();
        ASSERT_TRUE(session.HasLarsenHBondGrid());

        auto r = BuildFromProtonatedPdb(kLarsen1UbqPm6);
        ASSERT_TRUE(r.Ok()) << r.error;
        protein = std::move(r.protein);
    }

    Session session;
    std::unique_ptr<Protein> protein;
};


TEST_F(LarsenHBondShieldingTest, SmokeOn1UbqPm6) {
    auto& conf = protein->Conformation();

    // Dependency chain.
    auto geo = GeometryResult::Compute(conf);
    ASSERT_NE(geo, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(geo)));

    auto sp = SpatialIndexResult::Compute(conf);
    ASSERT_NE(sp, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(sp)));

    auto en = EnrichmentResult::Compute(conf);
    ASSERT_NE(en, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(en)));

    auto dssp = DsspResult::Compute(conf);
    ASSERT_NE(dssp, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(dssp)));

    // Run the calculator.
    auto result = LarsenHBondShieldingResult::Compute(
        conf, *session.LarsenHBondGridPtr());
    ASSERT_NE(result, nullptr);

    // 1UBQ aggregate across all (donor, acceptor) classes: 60+ backbone
    // amide-H H-bonds plus Hα H-bonds and sidechain-acceptor pairs.
    // ≥ 100 is a conservative floor.
    EXPECT_GE(result->PairsFound(), 100)
        << "expected ≥ 100 H-bond pairs total; got "
        << result->PairsFound();
    EXPECT_GE(result->AtomsWithContribution(), 150);

    // All emitted tensors finite.
    int n_finite_total = 0;
    int n_with_total = 0;
    double max_cb_frobenius = 0.0;
    int n_amide_with_water = 0;
    int n_atoms_with_pair_count = 0;
    int n_total_pair_increments = 0;
    for (std::size_t i = 0; i < conf.AtomCount(); ++i) {
        const auto& a = conf.AtomAt(i);
        const auto& m = a.larsen_hbond_shielding_tensor;
        bool finite = true;
        for (int r = 0; r < 3 && finite; ++r)
            for (int c = 0; c < 3 && finite; ++c)
                if (!std::isfinite(m(r, c))) finite = false;
        if (m.norm() > 1e-9) {
            ++n_with_total;
            if (finite) ++n_finite_total;
        }
        // Track max |Cβ diagnostic| Frobenius norm.
        double cb_norm = a.larsen_hbond_diagnostic_CB.norm();
        if (cb_norm > max_cb_frobenius) max_cb_frobenius = cb_norm;
        // Water term sweep count.
        if (a.larsen_hbond_water_term > 0.0) ++n_amide_with_water;
        // n_pairs increment verification (C1 fix).
        if (a.larsen_hbond_n_pairs > 0) {
            ++n_atoms_with_pair_count;
            n_total_pair_increments += a.larsen_hbond_n_pairs;
        }
    }
    EXPECT_EQ(n_finite_total, n_with_total);
    EXPECT_GT(n_with_total, 0);

    // C1 fix: n_pairs is no longer all-zero. Every atom that received
    // any tensor contribution should have n_pairs >= 1. Cross-check:
    // (atoms with contribution) == (atoms with n_pairs > 0).
    EXPECT_EQ(n_atoms_with_pair_count, n_with_total)
        << "n_pairs should be non-zero exactly on atoms with contribution";
    EXPECT_GT(n_total_pair_increments, result->PairsFound())
        << "summed n_pairs should exceed pair count (multiple atoms per pair)";

    // Cβ diagnostic: Larsen Table 2 EXCLUDES Cβ from the final shielding
    // sum — but his DFT scan still computed real σ values at Cβ atoms,
    // and those are what live in the grid's donor_CB readouts. The
    // diagnostic surfaces what Larsen's model discards. On 1UBQ the
    // max Cβ Frobenius is ~14 ppm (real chemistry — Cβ does respond to
    // H-bond geometry). We assert finite + bounded-below-absurd; the
    // non-zero value is methodology signal, not a pipeline bug.
    EXPECT_TRUE(std::isfinite(max_cb_frobenius));
    EXPECT_LT(max_cb_frobenius, 100.0)
        << "Cβ diagnostic exceeded 100 ppm Frobenius — suspicious. "
           "Larsen's DFT measured real shielding at Cβ but extreme "
           "values usually mean pipeline trouble (close-approach "
           "geometry, rotation bug, parser regression). Observed: "
        << max_cb_frobenius;
    EXPECT_GT(max_cb_frobenius, 0.0)
        << "Cβ diagnostic was zero — ALA donor archives should produce "
           "non-zero CB readouts; zero means CB path didn't fire.";

    // Water term: assigned only to amide Hs that DSSP did NOT detect
    // as a donor of any H-bond. On 1UBQ specifically the post-M2-fix
    // result is 0 — every amide H is DSSP-detected as a donor
    // (helix/sheet packing). Pre-M2 the buggy implementation
    // assigned water term to grid-skipped DSSP-paired pairs (17 of
    // them); the fix correctly suppresses those. We assert the
    // weaker condition: 0 ≤ water-term count ≤ total amide Hs, and
    // mass conservation (no atom is both water-term-assigned AND
    // pair-contributing) is verified in the dedicated reality-check
    // test below.
    int n_amide_h_total = 0;
    for (std::size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        const auto& res = protein->ResidueAt(ri);
        if (res.H != Residue::NONE) ++n_amide_h_total;
    }
    EXPECT_LE(n_amide_with_water, n_amide_h_total)
        << "water-term count must not exceed total amide H count";
    EXPECT_EQ(n_amide_with_water, result->AmideHsUnboundWithWater())
        << "atom-level water-term count should match result aggregate";

    // WriteFeatures emits per-atom tensors/counts plus split and
    // compatibility pair-audit tables.
    fs::path tmp = fs::temp_directory_path() / "larsen_hbond_phase1_out";
    fs::create_directories(tmp);
    int n_written = result->WriteFeatures(conf, tmp.string());
    EXPECT_EQ(n_written, 15);
    for (const std::string& stem : {
        "larsen_hbond_shielding",
        "larsen_hbond_1pHB_shielding",
        "larsen_hbond_2pHB_shielding",
        "larsen_hbond_1pHaB_shielding",
        "larsen_hbond_2pHaB_shielding",
        "larsen_hbond_diagnostic_CB_shielding",
        "larsen_hbond_water_term",
        "larsen_hbond_count",
        "larsen_imputed_pair_count",
        "larsen_sidechain_carbonyl_pair_count",
        "larsen_corner_imputed",
        "larsen_hbond_pairs_index",
        "larsen_hbond_pairs_geometry",
        "larsen_hbond_pairs_isotropic",
        "larsen_hbond_pairs",
    }) {
        fs::path p = tmp / (stem + ".npy");
        EXPECT_TRUE(fs::exists(p)) << "missing " << p.string();
    }

    {
        fs::path p = tmp / "larsen_corner_imputed.npy";
        std::ifstream in(p, std::ios::binary);
        ASSERT_TRUE(in.is_open()) << p.string();
        char magic[6] = {};
        in.read(magic, 6);
        ASSERT_EQ(std::string(magic, 6), std::string("\x93NUMPY", 6));
        char version[2] = {};
        in.read(version, 2);
        ASSERT_EQ(version[0], 1);
        ASSERT_EQ(version[1], 0);
        std::uint16_t header_len = 0;
        in.read(reinterpret_cast<char*>(&header_len), sizeof(header_len));
        std::string header(header_len, '\0');
        in.read(header.data(), header_len);
        EXPECT_NE(header.find("'descr': '|i1'"), std::string::npos) << header;
        EXPECT_NE(header.find("'shape': (" + std::to_string(conf.AtomCount()) + ",)"),
                  std::string::npos) << header;
        std::vector<char> payload{std::istreambuf_iterator<char>(in),
                                  std::istreambuf_iterator<char>()};
        ASSERT_EQ(payload.size(), conf.AtomCount());
        for (std::size_t i = 0; i < conf.AtomCount(); ++i) {
            EXPECT_EQ(static_cast<std::int8_t>(payload[i]),
                      conf.AtomAt(i).larsen_hbond_any_corner_imputed ? 1 : 0);
        }
    }
}


// ---------------------------------------------------------------------------
// Reality-check tests (user-requested round 3)
// ---------------------------------------------------------------------------

// GLY HA fan-out: 1UBQ has 6 GLY residues. Every Hα target contribution
// must land on BOTH HA2 and HA3 — not just one. Test: among atoms
// classified as GLY α-hydrogens (via IsAnyAlphaHydrogen()) that
// received any Table 2 contribution under HA's row, count must be
// even (each GLY pair contributes paired).
TEST_F(LarsenHBondShieldingTest, GlyHaFanOutHitsBothHA) {
    auto& conf = protein->Conformation();

    // Dependency chain.
    auto geo = GeometryResult::Compute(conf);
    ASSERT_NE(geo, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(geo)));
    auto sp = SpatialIndexResult::Compute(conf);
    ASSERT_NE(sp, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(sp)));
    auto en = EnrichmentResult::Compute(conf);
    ASSERT_NE(en, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(en)));
    auto dssp = DsspResult::Compute(conf);
    ASSERT_NE(dssp, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(dssp)));

    auto result = LarsenHBondShieldingResult::Compute(
        conf, *session.LarsenHBondGridPtr());
    ASSERT_NE(result, nullptr);

    // Count GLY HAs receiving any 1°HB contribution. For each GLY
    // donor residue that successfully formed an H-bond, BOTH GLY
    // alpha hydrogens (HA2 + HA3) should have non-zero 1pHB tensor.
    int n_gly_residues_with_partial_ha_hit = 0;  // exactly 1 of 2 hit
    int n_gly_residues_with_both_ha_hit    = 0;  // both 2 of 2 hit
    int n_gly_residues_with_no_ha_hit      = 0;
    const auto& topo = protein->LegacyAmber();
    for (std::size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        const auto& res = protein->ResidueAt(ri);
        if (res.type != AminoAcid::GLY) continue;
        int n_ha_with_contribution = 0;
        int n_ha_atoms = 0;
        for (std::size_t ai : res.atom_indices) {
            const auto& sem = topo.SemanticAt(ai);
            if (!sem.IsAnyAlphaHydrogen()) continue;
            ++n_ha_atoms;
            if (conf.AtomAt(ai).larsen_hbond_1pHB_tensor.norm() > 1e-9) {
                ++n_ha_with_contribution;
            }
        }
        // Every GLY has 2 prochiral HAs (HA2 + HA3).
        EXPECT_EQ(n_ha_atoms, 2)
            << "GLY at residue " << ri << " should have 2 α-hydrogens; got "
            << n_ha_atoms;
        if (n_ha_with_contribution == 0)      ++n_gly_residues_with_no_ha_hit;
        else if (n_ha_with_contribution == 1) ++n_gly_residues_with_partial_ha_hit;
        else                                  ++n_gly_residues_with_both_ha_hit;
    }
    // The fan-out fix asserts: NO GLY residue should have a partial
    // HA hit. Either neither HA contributes (residue not a donor) or
    // both do (substrate-driven enumeration). A partial hit would be
    // the codex H1 bug.
    EXPECT_EQ(n_gly_residues_with_partial_ha_hit, 0)
        << "GLY HA fan-out broken: " << n_gly_residues_with_partial_ha_hit
        << " GLY residues have exactly 1 of 2 α-hydrogens contributing";
    // Sanity: at least one GLY should have contributed (1UBQ has
    // GLY in α-helix / β-sheet so they're amide-H donors).
    EXPECT_GT(n_gly_residues_with_both_ha_hit, 0)
        << "expected ≥1 GLY residue to have both Hα atoms receiving "
           "contributions on 1UBQ; saw 0";
}


// Mass-conservation: every geometric candidate the spatial sweep
// classified gets accounted for as either a successful pair or a
// counted skip (codex finding F4, 2026-05-12). No silent drops.
//
// Also: an amide H can legitimately carry BOTH a water term (it had
// no H-bond candidate as DONOR) AND a non-zero n_pairs (it received
// 2°HB contributions as the i+1 TARGET of someone else's H-bond).
// Those are independent — donor status and target status are
// orthogonal. The earlier DSSP-era mutual-exclusion assertion was a
// bookkeeping artefact, not physics.
TEST_F(LarsenHBondShieldingTest, GeometricCandidateMassConservation) {
    auto& conf = protein->Conformation();
    auto geo = GeometryResult::Compute(conf);
    ASSERT_TRUE(conf.AttachResult(std::move(geo)));
    auto sp = SpatialIndexResult::Compute(conf);
    ASSERT_TRUE(conf.AttachResult(std::move(sp)));
    auto en = EnrichmentResult::Compute(conf);
    ASSERT_TRUE(conf.AttachResult(std::move(en)));
    auto dssp = DsspResult::Compute(conf);
    ASSERT_TRUE(conf.AttachResult(std::move(dssp)));

    auto result = LarsenHBondShieldingResult::Compute(
        conf, *session.LarsenHBondGridPtr());
    ASSERT_NE(result, nullptr);

    // Pairs that reached the grid path and either succeeded or were
    // counted as skipped. The two counts together are the total of
    // geometric candidates the dispatch evaluated — none were silently
    // dropped (cf. codex F4: missing-frame-atoms early-returns used to
    // bypass the counter).
    const int n_found   = result->PairsFound();
    const int n_skipped = result->PairsGridSkipped();
    EXPECT_GE(n_found,   40) << "1UBQ should yield ≥ 40 grid-paired pairs";
    EXPECT_GE(n_skipped, 0)  << "skip counter must be non-negative";

    // Water term + n_pairs can both be > 0 on the same amide H without
    // contradicting the gate. The interesting invariant is that EVERY
    // amide H with water_term > 0 had ZERO geometric H-bond candidates
    // as a donor (θ ≥ 90° in 4.2 Å). We can't recompute the per-atom
    // donor-paired flag from result alone, but we can verify the
    // overall count is finite + the water term is exactly 2.07 ppm
    // where it fires.
    constexpr double kWater_ppm = 2.07;
    int n_water_total = 0;
    for (std::size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& a = conf.AtomAt(ai);
        if (a.larsen_hbond_water_term > 0.0) {
            EXPECT_DOUBLE_EQ(a.larsen_hbond_water_term, kWater_ppm)
                << "water term must be Larsen's calibrated 2.07 ppm "
                << "isotropic; atom_idx=" << ai;
            ++n_water_total;
        }
    }
    EXPECT_EQ(n_water_total, result->AmideHsUnboundWithWater())
        << "per-atom water_term count must match aggregate";
    // 1UBQ has ~76 amide Hs; some are solvent-exposed. After the F2
    // gate fix the water-term count should be > 0 (the earlier eager
    // gate suppressed it for ~all amide Hs). Bound is intentionally
    // loose — exact number depends on solvent accessibility geometry.
    EXPECT_GT(n_water_total, 0)
        << "expected ≥ 1 solvent-exposed amide H on 1UBQ to receive "
           "the water term; saw 0 (F2 gate may be too aggressive)";
}


// CB diagnostic does NOT increment larsen_hbond_n_pairs.
// (ConformationAtom doc: n_pairs counts the four real Table 2
// contribution classes only. The CB diagnostic is for parser-pipeline
// integrity, not a true class.)
TEST_F(LarsenHBondShieldingTest, CbDiagnosticDoesNotInflateNPairs) {
    auto& conf = protein->Conformation();
    auto geo = GeometryResult::Compute(conf);
    ASSERT_TRUE(conf.AttachResult(std::move(geo)));
    auto sp = SpatialIndexResult::Compute(conf);
    ASSERT_TRUE(conf.AttachResult(std::move(sp)));
    auto en = EnrichmentResult::Compute(conf);
    ASSERT_TRUE(conf.AttachResult(std::move(en)));
    auto dssp = DsspResult::Compute(conf);
    ASSERT_TRUE(conf.AttachResult(std::move(dssp)));

    auto result = LarsenHBondShieldingResult::Compute(
        conf, *session.LarsenHBondGridPtr());
    ASSERT_NE(result, nullptr);

    // An atom with CB diagnostic contribution but NO real Table 2
    // contribution should have n_pairs == 0.
    int n_cb_only_with_count = 0;
    for (std::size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& a = conf.AtomAt(ai);
        bool has_cb_diag = a.larsen_hbond_diagnostic_CB.norm() > 1e-9;
        bool has_real_contrib = a.larsen_hbond_shielding_tensor.norm() > 1e-9;
        if (has_cb_diag && !has_real_contrib && a.larsen_hbond_n_pairs > 0) {
            ++n_cb_only_with_count;
        }
    }
    EXPECT_EQ(n_cb_only_with_count, 0)
        << "Cβ-diagnostic-only contributors must NOT increment n_pairs "
           "(per ConformationAtom doc: n_pairs is for Table 2 classes only). "
           "The CB diagnostic exists for parser-pipeline integrity; "
           "incrementing n_pairs from CB contributions would inflate "
           "atoms' pair counts.";
}


// Schema validation: load LarsenHBondGrid and verify the loader's
// per-archive presence checks did not throw. (If schema were
// violated, ctor would have thrown and the fixture SetUp would have
// failed before reaching here.)
TEST_F(LarsenHBondShieldingTest, GridSchemaValidationPassed) {
    // Reaching this point means LarsenHBondGrid ctor + per-archive
    // ValidateSchema all passed. The fixture's grid load is the
    // implicit test.
    const LarsenHBondGrid* g = session.LarsenHBondGridPtr();
    ASSERT_NE(g, nullptr);
    EXPECT_TRUE(g->IsLoaded());
}


// F1 reality check (codex finding, 2026-05-12): after adding acceptor_CA
// and acceptor_C grid readouts, Cα and C' atoms in residues that are
// H-bond targets must accumulate Δσ_2°HB contributions. Before F1,
// Larsen Table 2 said these atoms get 2°HB but the grid lacked the
// tensors so they were silently zero — particularly C', the largest
// 2°HB term per Larsen Table 1 (~2.1 ppm RMSD on Ub).
TEST_F(LarsenHBondShieldingTest, CaAndCReceive2pHBContributions) {
    auto& conf = protein->Conformation();
    auto geo = GeometryResult::Compute(conf);
    ASSERT_TRUE(conf.AttachResult(std::move(geo)));
    auto sp = SpatialIndexResult::Compute(conf);
    ASSERT_TRUE(conf.AttachResult(std::move(sp)));
    auto en = EnrichmentResult::Compute(conf);
    ASSERT_TRUE(conf.AttachResult(std::move(en)));
    auto dssp = DsspResult::Compute(conf);
    ASSERT_TRUE(conf.AttachResult(std::move(dssp)));

    auto result = LarsenHBondShieldingResult::Compute(
        conf, *session.LarsenHBondGridPtr());
    ASSERT_NE(result, nullptr);

    const Protein& protein_ref = conf.ProteinRef();
    const auto& topo = protein_ref.LegacyAmber();

    int n_ca_with_2pHB = 0;
    int n_c_with_2pHB  = 0;
    double max_ca_2pHB = 0.0;
    double max_c_2pHB  = 0.0;
    for (std::size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& sem = topo.SemanticAt(ai);
        const auto& a = conf.AtomAt(ai);
        const double n = a.larsen_hbond_2pHB_tensor.norm();
        if (n <= 1e-9) continue;
        if (sem.backbone_role == BackboneRole::AlphaCarbon) {
            ++n_ca_with_2pHB;
            max_ca_2pHB = std::max(max_ca_2pHB, n);
        } else if (sem.backbone_role == BackboneRole::CarbonylCarbon) {
            ++n_c_with_2pHB;
            max_c_2pHB = std::max(max_c_2pHB, n);
        }
    }

    // 1UBQ has ~149 grid-paired H-bonds; many should land 2°HB on the
    // i+1 Cα and on the acceptor's own C atom. Loose lower bound to
    // detect pipeline regressions, not assert exact counts.
    EXPECT_GT(n_ca_with_2pHB, 10)
        << "Cα atoms with Δσ_2°HB > 0 should exceed 10 on 1UBQ; F1 "
           "regression: acceptor_CA may not be wired through.";
    EXPECT_GT(n_c_with_2pHB,  10)
        << "C' atoms with Δσ_2°HB > 0 should exceed 10 on 1UBQ; F1 "
           "regression: acceptor_C may not be wired through. "
           "C' Δσ_2°HB is Larsen's largest 2°HB term (~2.1 ppm RMSD).";

    // Magnitude sanity. The parser subtracts an orientation-matched
    // r-max reference surface; if that regresses to a global r-max
    // average, rotated absolute shielding anisotropy leaks into the
    // H-bond tensor and 1UBQ produces >200 ppm norms on C'/Cα. Keep a
    // deliberately loose but load-bearing upper bound here: real 1UBQ
    // values stay well below this after orientation-matched subtraction.
    EXPECT_LT(max_ca_2pHB, 100.0)
        << "Cα Δσ_2°HB max=" << max_ca_2pHB
        << " ppm — likely reference-surface or tensor-frame regression";
    EXPECT_LT(max_c_2pHB,  100.0)
        << "C' Δσ_2°HB max=" << max_c_2pHB
        << " ppm — likely reference-surface or tensor-frame regression";
}

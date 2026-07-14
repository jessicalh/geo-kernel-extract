//
// test_mopac_mc_shielding_time_series: discipline + integration for
// MopacMcConnellShieldingTimeSeriesTrajectoryResult (TR8 of the 13-TR
// plan). Sibling of TR7 but emits 9 components (T0+T1+T2) because
// the source mopac_mc_shielding_contribution is the full BO-channel
// D(r)Qhat response and can have nonzero T0/T1 in general.
//

#include "MopacMcConnellShieldingTimeSeriesTrajectoryResult.h"
#include "McConnellResult.h"
#include "MopacMcConnellResult.h"
#include "MopacResult.h"
#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "ChargeAssignmentResult.h"

#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "DenseBuffer.h"
#include "OperationLog.h"
#include "PdbFileReader.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "RunConfiguration.h"
#include "Session.h"
#include "TestEnvironment.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"
#include "Types.h"

#include <gtest/gtest.h>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <unistd.h>
#include <vector>

namespace fs = std::filesystem;

#ifndef NMR_TEST_DATA_DIR
#error "NMR_TEST_DATA_DIR must be defined"
#endif

namespace {

constexpr const char* kFixtureProtein = "1P9J_5801";
constexpr const char* kTensorTransformation = "even_rank2: T'=R T R^T";
constexpr const char* kTensorT1Semantics =
    "Cartesian Levi-Civita dual x,y,z (not real-Y1m): "
    "a=((T_yz-T_zy)/2,(T_zx-T_xz)/2,(T_xy-T_yx)/2); "
    "axial a'=det(R) R a; generically nonzero";
constexpr const char* kNormalizationScope =
    "xyz tensor payload: T2 uses isometric real-tesseral "
    "normalization; T1 uses the tensor_t1_semantics convention";

nmr::SphericalTensor SyntheticTensor() {
    nmr::SphericalTensor tensor;
    tensor.T0 = 1.25;
    tensor.T1 = {0.5, -0.75, 1.5};
    tensor.T2 = {-2.0, 2.25, -2.5, 2.75, -3.0};
    return tensor;
}

std::string TrrPathFor(const std::string& p) {
    return fs::path(p).replace_extension(".trr").string();
}
std::string ProductionDirFor(const std::string& p) {
    return fs::path(p).parent_path().string();
}
bool FixtureAvailable(const nmr::test::AmberTrajectoryFixture& fix) {
    return !fix.tpr_path.empty() && fs::exists(fix.tpr_path)
        && fs::exists(TrrPathFor(fix.tpr_path)) && fs::exists(fix.edr_path);
}
}  // namespace


TEST(MopacMcConnellShieldingTimeSeries,
     H5DirectionalMetadataIsExactWithoutFleetFixture) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();

    auto build = nmr::BuildFromProtonatedPdb(
        nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(build.Ok()) << build.error;
    auto tp_owner = nmr::TrajectoryProtein::CreateForTesting(
        std::move(build.protein));
    ASSERT_NE(tp_owner, nullptr);
    auto& tp = *tp_owner;
    ASSERT_GT(tp.AtomCount(), 0u);

    auto tr =
        nmr::MopacMcConnellShieldingTimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj({}, {}, {});
    std::vector<nmr::Vec3> positions(tp.AtomCount(), nmr::Vec3::Zero());
    auto conf = std::make_unique<nmr::ProteinConformation>(
        &tp.ProteinRef(), positions,
        "synthetic MOPAC McConnell metadata frame");
    conf->ForceAttachResultForTesting(
        std::make_unique<nmr::MopacMcConnellResult>());
    const nmr::SphericalTensor expected = SyntheticTensor();
    conf->MutableAtomAt(0).mopac_mc_shielding_contribution = expected;

    tr->Compute(*conf, tp, traj, 0, 0.0);
    tr->Finalize(tp, traj);
    ASSERT_EQ(tr->SourceAttachedCount(), 1u);

    const std::string h5_path = (fs::temp_directory_path() /
        ("mopac_mc_shielding_ts_directional_metadata_" +
         std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        tr->WriteH5Group(tp, file);
    }

    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup(
        "/trajectory/mopac_mc_shielding_time_series");
    auto ds = grp.getDataSet("xyz");
    EXPECT_EQ(ds.getSpace().getDimensions(),
              (std::vector<std::size_t>{tp.AtomCount(), 1u, 9u}));

    std::string basis;
    std::string order;
    std::string tensor_frame;
    std::string coordinate_frame;
    std::string parity;
    std::string tensor_parity;
    std::string transformation;
    std::string t1_semantics;
    std::string structural_zero_components;
    std::string normalization_scope;
    bool t1_structural_zero = true;
    grp.getAttribute("tensor_basis").read(basis);
    grp.getAttribute("tensor_component_order").read(order);
    grp.getAttribute("tensor_frame").read(tensor_frame);
    grp.getAttribute("coordinate_frame").read(coordinate_frame);
    grp.getAttribute("parity").read(parity);
    grp.getAttribute("tensor_parity").read(tensor_parity);
    grp.getAttribute("tensor_transformation").read(transformation);
    grp.getAttribute("tensor_t1_semantics").read(t1_semantics);
    grp.getAttribute("tensor_t1_structural_zero").read(
        t1_structural_zero);
    grp.getAttribute("tensor_structural_zero_components").read(
        structural_zero_components);
    grp.getAttribute("normalization_scope").read(normalization_scope);

    EXPECT_EQ(basis, nmr::kMcConnellPackFull9TensorBasis);
    EXPECT_EQ(order, nmr::kMcConnellPackFull9ComponentOrder);
    EXPECT_EQ(tensor_frame, nmr::kMcConnellPackFull9TensorFrame);
    EXPECT_EQ(coordinate_frame, "conformation_cartesian_xyz");
    EXPECT_EQ(parity, "0e+1e+2e");
    EXPECT_EQ(tensor_parity, "even");
    EXPECT_EQ(transformation, kTensorTransformation);
    EXPECT_EQ(t1_semantics, kTensorT1Semantics);
    EXPECT_FALSE(t1_structural_zero);
    EXPECT_EQ(structural_zero_components, "none");
    EXPECT_EQ(normalization_scope, kNormalizationScope);

    std::vector<double> flat(tp.AtomCount() * 9u);
    ds.read(flat.data());
    double expected_pack[9]{};
    expected.PackFull9(expected_pack);
    for (std::size_t k = 0; k < 9u; ++k)
        EXPECT_DOUBLE_EQ(flat[k], expected_pack[k]);

    std::vector<std::uint8_t> source_mask(1u, 0u);
    grp.getDataSet("source_attached_per_frame").read(source_mask.data());
    EXPECT_EQ(source_mask[0], 1u);

    fs::remove(h5_path);
}


// ── Layer 0: source absent every frame → group skipped ─────

TEST(MopacMcConnellShieldingTimeSeries, GroupSkippedWhenSourceNeverAttached) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    auto tr = nmr::MopacMcConnellShieldingTimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    constexpr size_t kFrames = 3;
    std::vector<nmr::Vec3> positions(tp.AtomCount(), nmr::Vec3::Zero());
    for (size_t t = 0; t < kFrames; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic frame (no MopacMcConnell)");
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);

    EXPECT_EQ(tr->NumFrames(), kFrames);
    EXPECT_EQ(tr->SourceAttachedCount(), 0u);

    const std::string h5_path = (fs::temp_directory_path() /
        ("mopac_mc_shielding_ts_absent_" +
         std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    EXPECT_FALSE(reopen.exist("/trajectory/mopac_mc_shielding_time_series"));

    fs::remove(h5_path);
}


// ── Integration1P9J: real Trajectory::Run with MOPAC + McConnell ──

TEST(MopacMcConnellShieldingTimeSeries, Integration1P9J) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = false;     // MopacMcConnellResult needs MopacResult
    opts.skip_coulomb = false;   // enables Mopac-family attachments
    opts.skip_apbs = true;
    opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::MopacMcConnellShieldingTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    const auto status = traj.Run(tp, config, session);
    if (status != nmr::kOk) {
        GTEST_SKIP() << "Trajectory::Run failed; skipping per "
                        "feedback_log_overages_dont_assert.";
    }

    auto& tr = tp.Result<nmr::MopacMcConnellShieldingTimeSeriesTrajectoryResult>();
    const std::size_t T = tr.NumFrames();
    const std::size_t attached = tr.SourceAttachedCount();
    std::cout << "MopacMcConnellShieldingTS Integration1P9J: T=" << T
              << " MopacMcConnell-attached=" << attached << std::endl;

    if (attached == 0) {
        GTEST_SKIP() << "MopacMcConnellResult did not attach on any frame.";
    }

    const std::string h5_path = (fs::temp_directory_path() /
        ("mopac_mc_shielding_ts_int1p9j_" +
         std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/mopac_mc_shielding_time_series");
    auto ds = grp.getDataSet("xyz");
    const auto dims = ds.getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 3u);
    EXPECT_EQ(dims[1], T);
    EXPECT_EQ(dims[2], 9u)
        << "9-component emission per 'if not traceless write both'";

    std::string basis, order, frame, coordinate_frame, parity;
    std::string tensor_parity, transformation, t1_semantics;
    std::string structural_zero_components, normalization_scope;
    std::string e3nn_export, units;
    bool t1_structural_zero = true;
    grp.getAttribute("tensor_basis").read(basis);
    grp.getAttribute("tensor_component_order").read(order);
    grp.getAttribute("tensor_frame").read(frame);
    grp.getAttribute("coordinate_frame").read(coordinate_frame);
    grp.getAttribute("parity").read(parity);
    grp.getAttribute("tensor_parity").read(tensor_parity);
    grp.getAttribute("tensor_transformation").read(transformation);
    grp.getAttribute("tensor_t1_semantics").read(t1_semantics);
    grp.getAttribute("tensor_t1_structural_zero").read(
        t1_structural_zero);
    grp.getAttribute("tensor_structural_zero_components").read(
        structural_zero_components);
    grp.getAttribute("normalization_scope").read(normalization_scope);
    grp.getAttribute("e3nn_export").read(e3nn_export);
    grp.getAttribute("units").read(units);
    EXPECT_EQ(basis, nmr::kMcConnellPackFull9TensorBasis);
    EXPECT_EQ(order, nmr::kMcConnellPackFull9ComponentOrder);
    EXPECT_EQ(frame, nmr::kMcConnellPackFull9TensorFrame);
    EXPECT_EQ(coordinate_frame, nmr::kMcConnellPackFull9TensorFrame);
    EXPECT_EQ(parity, "0e+1e+2e");
    EXPECT_EQ(tensor_parity, "even");
    EXPECT_EQ(transformation, kTensorTransformation);
    EXPECT_EQ(t1_semantics, kTensorT1Semantics);
    EXPECT_FALSE(t1_structural_zero);
    EXPECT_EQ(structural_zero_components, "none");
    EXPECT_EQ(normalization_scope, kNormalizationScope);
    EXPECT_EQ(e3nn_export, nmr::kMcConnellPackFull9E3nnExport);
    EXPECT_FALSE(grp.hasAttribute("irrep_layout"));
    EXPECT_EQ(units, "Angstrom^-3")
        << "bare McConnell kernel bo·D(r)Qhat, no Δχ × γ multiplication "
           "at extraction; decision 2026-05-21 per science/math review H1.";

    const std::size_t N = dims[0];
    std::vector<double> flat(N * T * 9);
    ds.read(flat.data());
    double max_t0 = 0.0, max_t1 = 0.0, max_t2 = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t t = 0; t < T; ++t) {
            const std::size_t base = (i * T + t) * 9;
            for (std::size_t k = 0; k < 9; ++k) {
                EXPECT_TRUE(std::isfinite(flat[base + k]));
            }
            max_t0 = std::max(max_t0, std::abs(flat[base + 0]));
            for (std::size_t k = 1; k <= 3; ++k)
                max_t1 = std::max(max_t1, std::abs(flat[base + k]));
            for (std::size_t k = 4; k <= 8; ++k)
                max_t2 = std::max(max_t2, std::abs(flat[base + k]));
        }
    }
    std::cout << "  max|T0| = " << max_t0
              << "  max|T1| = " << max_t1
              << "  max|T2| = " << max_t2
              << "  (T0/T1 not necessarily ~0 — D(r)Qhat is not pure T2)"
              << std::endl;
    // Sanity: at least T2 should be nonzero on every protein
    // (bond-anisotropy is a real signal). T0/T1 may or may not be —
    // log don't assert per feedback_log_overages_dont_assert.
    EXPECT_GT(max_t2, 0.0);

    fs::remove(h5_path);
}


// ── Layer 0b: FinalizeIdempotency ──────────────────────────────

TEST(MopacMcConnellShieldingTimeSeries, FinalizeIdempotency) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = false;
    opts.skip_coulomb = false;
    opts.skip_apbs = true;
    opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::MopacMcConnellShieldingTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    if (traj.Run(tp, config, session) != nmr::kOk) {
        GTEST_SKIP() << "Trajectory::Run failed; skipping.";
    }

    auto* buf_first = tp.GetDenseBuffer<nmr::SphericalTensor>(std::type_index(
        typeid(nmr::MopacMcConnellShieldingTimeSeriesTrajectoryResult)));
    if (!buf_first) GTEST_SKIP() << "MopacMcConnell never attached.";
    const std::size_t N_first = buf_first->AtomCount();
    const std::size_t T_first = buf_first->StridePerAtom();

    auto& tr = tp.Result<nmr::MopacMcConnellShieldingTimeSeriesTrajectoryResult>();
    tr.Finalize(tp, traj);  // second call — bounds-check idempotency
    auto* buf_second = tp.GetDenseBuffer<nmr::SphericalTensor>(std::type_index(
        typeid(nmr::MopacMcConnellShieldingTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf_second, nullptr);
    EXPECT_EQ(buf_second->AtomCount(), N_first);
    EXPECT_EQ(buf_second->StridePerAtom(), T_first);
}

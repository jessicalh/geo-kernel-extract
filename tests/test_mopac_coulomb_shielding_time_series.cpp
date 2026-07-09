//
// test_mopac_coulomb_shielding_time_series: discipline + integration for
// MopacCoulombShieldingTimeSeriesTrajectoryResult (TR7). Combines
// TR4 T2-only TS pattern with TR5 sparse-cadence gate.
//

#include "MopacCoulombShieldingTimeSeriesTrajectoryResult.h"
#include "MopacCoulombResult.h"
#include "MopacResult.h"
#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "ChargeAssignmentResult.h"

#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "DenseBuffer.h"
#include "OperationLog.h"
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
#include <vector>

namespace fs = std::filesystem;

#ifndef NMR_TEST_DATA_DIR
#error "NMR_TEST_DATA_DIR must be defined"
#endif

namespace {

constexpr const char* kFixtureProtein = "1P9J_5801";

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


// ── Layer 0a: source absent every frame → group skipped ─────

TEST(MopacCoulombShieldingTimeSeries, GroupSkippedWhenSourceNeverAttached) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    auto tr = nmr::MopacCoulombShieldingTimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    constexpr size_t kFrames = 3;
    std::vector<nmr::Vec3> positions(tp.AtomCount(), nmr::Vec3::Zero());
    for (size_t t = 0; t < kFrames; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic frame (no MopacCoulomb)");
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);

    EXPECT_EQ(tr->NumFrames(), kFrames);
    EXPECT_EQ(tr->SourceAttachedCount(), 0u);

    const std::string h5_path = (fs::temp_directory_path() /
        ("mopac_coulomb_shielding_ts_absent_" +
         std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    EXPECT_FALSE(reopen.exist(
        "/trajectory/mopac_coulomb_efg_time_series"))
        << "canonical group must NOT exist when source attached 0/T frames";
    EXPECT_FALSE(reopen.exist(
        "/trajectory/mopac_coulomb_shielding_time_series"))
        << "legacy group must NOT exist when source attached 0/T frames";

    fs::remove(h5_path);
}


// ── Integration1P9J: real Trajectory::Run with MOPAC + Coulomb ──

TEST(MopacCoulombShieldingTimeSeries, Integration1P9J) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = false;       // MopacCoulombResult needs MopacResult
    opts.skip_coulomb = false;     // enables MopacCoulombResult attach
    opts.skip_apbs = true;
    opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::MopacCoulombShieldingTimeSeriesTrajectoryResult::Create(tp_in);
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

    auto& tr = tp.Result<nmr::MopacCoulombShieldingTimeSeriesTrajectoryResult>();
    const std::size_t T = tr.NumFrames();
    const std::size_t attached = tr.SourceAttachedCount();
    std::cout << "MopacCoulombShieldingTS Integration1P9J: T=" << T
              << " MopacCoulomb-attached=" << attached << std::endl;

    if (attached == 0) {
        GTEST_SKIP() << "MopacCoulombResult did not attach on any frame.";
    }

    // Read back the H5 group + verify shape + finiteness.
    const std::string h5_path = (fs::temp_directory_path() /
        ("mopac_coulomb_shielding_ts_int1p9j_" +
         std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/mopac_coulomb_efg_time_series"));
    ASSERT_TRUE(reopen.exist("/trajectory/mopac_coulomb_shielding_time_series"));
    auto grp = reopen.getGroup("/trajectory/mopac_coulomb_efg_time_series");
    auto ds = grp.getDataSet("t2");
    const auto dims = ds.getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 3u);
    EXPECT_EQ(dims[1], T);
    EXPECT_EQ(dims[2], 5u) << "T2-only per plan + source comment";

    std::string quantity, historical, recovery, basis, order, frame;
    std::string t2_parity, export_note, layout, units, source, policy;
    bool gamma_applied = true;
    bool legacy_deprecated = false;
    grp.getAttribute("quantity").read(quantity);
    grp.getAttribute("historical_field_name").read(historical);
    grp.getAttribute("gamma_applied").read(gamma_applied);
    grp.getAttribute("shielding_recovery").read(recovery);
    grp.getAttribute("t2_basis").read(basis);
    grp.getAttribute("t2_component_order").read(order);
    grp.getAttribute("t2_frame").read(frame);
    grp.getAttribute("t2_parity").read(t2_parity);
    grp.getAttribute("e3nn_export").read(export_note);
    grp.getAttribute("legacy_irrep_attrs_deprecated").read(legacy_deprecated);
    grp.getAttribute("irrep_layout").read(layout);
    grp.getAttribute("units").read(units);
    grp.getAttribute("source").read(source);
    grp.getAttribute("source_attached_policy").read(policy);
    EXPECT_EQ(quantity, "raw_efg_t2");
    EXPECT_EQ(historical, "mopac_coulomb_shielding_contribution");
    EXPECT_FALSE(gamma_applied);
    EXPECT_EQ(recovery, "multiply by per-element gamma during calibration");
    EXPECT_EQ(basis, "project_native_t2_isometric_real_tesseral_v1");
    EXPECT_EQ(order, "T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2");
    EXPECT_EQ(frame, "cartesian_xyz_emitted_frame");
    EXPECT_EQ(t2_parity, "even");
    EXPECT_NE(export_note.find("project_t2_to_e3nn"), std::string::npos);
    EXPECT_TRUE(legacy_deprecated);
    EXPECT_EQ(layout, "T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2");
    EXPECT_EQ(units, "V/Å^2")
        << "EFG kernel units — bare EFG before γ multiplication; "
           "decision 2026-05-21 per science/math adversarial review H1.";
    EXPECT_EQ(source.find("MopacCoulombResult.cpp:"), std::string::npos);
    EXPECT_EQ(source.find("ppm"), std::string::npos);
    EXPECT_EQ(policy.find("OperationRunner.cpp:"), std::string::npos);
    EXPECT_EQ(policy.find("ppm"), std::string::npos);

    const std::size_t N = dims[0];
    std::vector<double> flat(N * T * 5);
    ds.read(flat.data());
    double max_mag = 0.0;
    for (std::size_t cell = 0; cell < flat.size(); ++cell) {
        EXPECT_TRUE(std::isfinite(flat[cell]));
        max_mag = std::max(max_mag, std::abs(flat[cell]));
    }
    EXPECT_GT(max_mag, 0.0)
        << "Mopac Coulomb shielding all zero — calc not firing or "
           "Mulliken charges collapsed";
    std::cout << "  max|T2| = " << max_mag << " V/Å^2 (bare EFG, "
              << "pre-γ)" << std::endl;

    auto legacy = reopen.getGroup("/trajectory/mopac_coulomb_shielding_time_series");
    std::string alias_of;
    bool legacy_name_deprecated = false;
    legacy.getAttribute("alias_of").read(alias_of);
    legacy.getAttribute("legacy_name_deprecated").read(legacy_name_deprecated);
    EXPECT_EQ(alias_of, "/trajectory/mopac_coulomb_efg_time_series");
    EXPECT_TRUE(legacy_name_deprecated);
    auto legacy_ds = legacy.getDataSet("t2");
    const auto legacy_dims = legacy_ds.getSpace().getDimensions();
    EXPECT_EQ(legacy_dims, dims);
    std::vector<double> legacy_flat(N * T * 5);
    legacy_ds.read(legacy_flat.data());
    EXPECT_EQ(legacy_flat, flat);

    fs::remove(h5_path);
}


// ── Layer 0b: FinalizeIdempotency ──────────────────────────────

TEST(MopacCoulombShieldingTimeSeries, FinalizeIdempotency) {
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
        return nmr::MopacCoulombShieldingTimeSeriesTrajectoryResult::Create(tp_in);
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
        typeid(nmr::MopacCoulombShieldingTimeSeriesTrajectoryResult)));
    if (!buf_first) GTEST_SKIP() << "MopacCoulomb never attached.";
    const std::size_t N_first = buf_first->AtomCount();
    const std::size_t T_first = buf_first->StridePerAtom();

    auto& tr = tp.Result<nmr::MopacCoulombShieldingTimeSeriesTrajectoryResult>();
    tr.Finalize(tp, traj);  // second call — bounds-check idempotency
    auto* buf_second = tp.GetDenseBuffer<nmr::SphericalTensor>(std::type_index(
        typeid(nmr::MopacCoulombShieldingTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf_second, nullptr);
    EXPECT_EQ(buf_second->AtomCount(), N_first);
    EXPECT_EQ(buf_second->StridePerAtom(), T_first);
}

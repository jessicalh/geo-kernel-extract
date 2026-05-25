//
// test_dssp8_time_series: discipline + integration for
// Dssp8TimeSeriesTrajectoryResult.
//

#include "CalculatorConfig.h"
#include "DsspResult.h"
#include "Dssp8TimeSeriesTrajectoryResult.h"
#include "EnrichmentResult.h"
#include "GeometryResult.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "RunConfiguration.h"
#include "Session.h"
#include "SpatialIndexResult.h"
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
nmr::RunConfiguration BuildConfig(unsigned stride) {
    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true;
    opts.skip_dssp = false;  // DSSP required for this TR's source
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::DsspResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::Dssp8TimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(stride);
    return config;
}

nmr::RunConfiguration BuildConfigNoDssp(unsigned stride) {
    auto cfg = BuildConfig(stride);
    cfg.MutablePerFrameRunOptions().skip_dssp = true;
    return cfg;
}

}  // namespace


TEST(Dssp8TimeSeries, Frame0Semantics) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildConfig(99999);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::Dssp8TimeSeriesTrajectoryResult>();
    EXPECT_EQ(tr.NumFrames(), 1u);
}


TEST(Dssp8TimeSeries, FinalizeIdempotency) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildConfig(99999);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    auto& tr = tp.Result<nmr::Dssp8TimeSeriesTrajectoryResult>();
    const std::size_t T = tr.NumFrames();
    tr.Finalize(tp, traj);
    EXPECT_EQ(tr.NumFrames(), T);
}


TEST(Dssp8TimeSeries, H5RoundTrip) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildConfig(99999);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::Dssp8TimeSeriesTrajectoryResult>();
    const std::size_t R = tp.ProteinRef().ResidueCount();

    const std::string h5_path = (fs::temp_directory_path() /
        ("dssp8_ts_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/dssp8_time_series"));
    auto grp = reopen.getGroup("/trajectory/dssp8_time_series");

    // ss8_code shape
    {
        const auto dims = grp.getDataSet("ss8_code").getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 2u);
        EXPECT_EQ(dims[0], R);
        EXPECT_EQ(dims[1], tr.NumFrames());
    }
    // H-bond datasets shape
    for (const std::string& name : {"hbond_acceptor_partner",
                                     "hbond_donor_partner"}) {
        const auto dims = grp.getDataSet(name).getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 3u);
        EXPECT_EQ(dims[0], R);
        EXPECT_EQ(dims[1], tr.NumFrames());
        EXPECT_EQ(dims[2], 2u);
    }
    for (const std::string& name : {"hbond_acceptor_energy",
                                     "hbond_donor_energy"}) {
        const auto dims = grp.getDataSet(name).getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 3u);
        EXPECT_EQ(dims[0], R);
        EXPECT_EQ(dims[2], 2u);
    }
    EXPECT_TRUE(grp.exist("residue_index_per_atom"));
    EXPECT_TRUE(grp.exist("source_attached_per_frame"));

    // Convention attrs
    std::string legend, units, policy;
    grp.getAttribute("ss8_legend").read(legend);
    grp.getAttribute("hbond_energy_units").read(units);
    grp.getAttribute("source_attached_policy").read(policy);
    EXPECT_NE(legend.find("alpha helix"), std::string::npos);
    EXPECT_EQ(units, "kcal/mol");
    EXPECT_NE(policy.find("conditional"), std::string::npos);

    fs::remove(h5_path);
}


TEST(Dssp8TimeSeries, Integration1P9J) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildConfig(300);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::Dssp8TimeSeriesTrajectoryResult>();
    const std::size_t R = tp.ProteinRef().ResidueCount();
    const std::size_t T = tr.NumFrames();
    EXPECT_GE(T, 2u);

    const std::string h5_path = (fs::temp_directory_path() /
        ("dssp8_ts_int_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/dssp8_time_series");

    // Range sanity: ss8_code in {0..7, 255}.
    std::vector<std::vector<std::uint8_t>> ss;
    grp.getDataSet("ss8_code").read(ss);
    ASSERT_EQ(ss.size(), R);
    std::size_t valid_codes = 0;
    for (const auto& row : ss) {
        for (auto c : row) {
            if (c <= 7) ++valid_codes;
            else EXPECT_EQ(c, 255u);
        }
    }
    EXPECT_GT(valid_codes, 0u);
    std::cout << "Dssp8TimeSeries: " << R << " residues, " << T
              << " frames; valid_codes=" << valid_codes << "/" << (R*T) << "\n";

    // H-bond energies <= 0 for stable bonds (kcal/mol negative); NaN
    // when no partner. Loose check on observed bonds.
    std::vector<std::vector<std::vector<double>>> acc_e;
    grp.getDataSet("hbond_acceptor_energy").read(acc_e);
    std::size_t observed_bonds = 0;
    for (const auto& row : acc_e) {
        for (const auto& frame : row) {
            for (double e : frame) {
                if (std::isfinite(e)) {
                    ++observed_bonds;
                    // libdssp writes the two best (lowest-energy)
                    // partner candidates to the slot regardless of the
                    // Kabsch-Sander 1983 strict -0.5 kcal/mol threshold
                    // (which gates SS-classification COUNTING, not
                    // slot-write). Empirically observed range on 1P9J
                    // includes E ∈ [-0.5, 0]: real but weak. Test:
                    // observed slots must be negative (attractive).
                    // Sign-flipped energy or positive would be a real
                    // libdssp serialisation bug.
                    EXPECT_LE(e, 0.0)
                        << "DSSP H-bond energy not attractive: " << e;
                }
            }
        }
    }
    EXPECT_GT(observed_bonds, 0u);
    std::cout << "  observed acceptor bonds: " << observed_bonds << "\n";

    fs::remove(h5_path);
}


TEST(Dssp8TimeSeries, SyntheticAllAbsentSkipsGroup) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    const std::size_t N = tp.AtomCount();
    auto tr = nmr::Dssp8TimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    std::vector<nmr::Vec3> positions(N, nmr::Vec3::Zero());
    for (std::size_t t = 0; t < 3; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic frame");
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);

    const std::string h5_path = (fs::temp_directory_path() /
        ("dssp8_ts_absent_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    EXPECT_FALSE(reopen.exist("/trajectory/dssp8_time_series"))
        << "All-absent run should skip group emission.";
    fs::remove(h5_path);
}

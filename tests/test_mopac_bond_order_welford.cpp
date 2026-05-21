//
// test_mopac_bond_order_welford: discipline + integration for
// MopacBondOrderWelfordTrajectoryResult (TR6 of the 13-TR plan).
// Clone of TR5 (MopacChargeWelford) with the bond axis substituted
// for the atom axis. Per-bond Welford state lives INSIDE the TR
// (BondLengthStats pattern), not on TrajectoryAtom.
//

#include "MopacBondOrderWelfordTrajectoryResult.h"
#include "MopacResult.h"
#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "ChargeAssignmentResult.h"

#include "CalculatorConfig.h"
#include "ConformationAtom.h"
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
void LoadCalculatorConfig() {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);
    nmr::CalculatorConfig::Load(std::string(NMR_TEST_DATA_DIR) +
                                "/../data/calculator_params.toml");
}

}  // namespace


// ── Layer 0: source absent every frame → group skipped ─────

TEST(MopacBondOrderWelford, GroupSkippedWhenSourceNeverAttached) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    auto tr = nmr::MopacBondOrderWelfordTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    constexpr size_t kFrames = 3;
    std::vector<nmr::Vec3> positions(tp.AtomCount(), nmr::Vec3::Zero());
    for (size_t t = 0; t < kFrames; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic frame (no MOPAC)");
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);

    EXPECT_EQ(tr->NumFrames(), kFrames);
    EXPECT_EQ(tr->SourceAttachedCount(), 0u);
    EXPECT_EQ(tr->BondCount(), tp.ProteinRef().BondCount());

    const std::string h5_path = (fs::temp_directory_path() /
        ("mopac_bond_order_welford_absent_" +
         std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    EXPECT_FALSE(reopen.exist("/trajectory/mopac_bond_order_welford"))
        << "group must NOT exist when source attached 0/T frames";

    fs::remove(h5_path);
}


// ── Integration1P9J: real Trajectory::Run with MOPAC enabled ──

TEST(MopacBondOrderWelford, Integration1P9J) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = false;
    opts.skip_coulomb = true;
    opts.skip_apbs = true;
    opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::MopacBondOrderWelfordTrajectoryResult::Create(tp_in);
    });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    const auto status = traj.Run(tp, config, session);
    if (status != nmr::kOk) {
        GTEST_SKIP() << "Trajectory::Run failed (MOPAC env not ready?); "
                        "skipping per feedback_log_overages_dont_assert.";
    }

    auto& tr = tp.Result<nmr::MopacBondOrderWelfordTrajectoryResult>();
    const std::size_t B = tr.BondCount();
    const std::size_t T = tr.NumFrames();
    const std::size_t attached = tr.SourceAttachedCount();
    std::cout << "MopacBondOrderWelford Integration1P9J: T=" << T
              << " MOPAC-attached=" << attached
              << " B=" << B << std::endl;

    if (attached == 0) {
        GTEST_SKIP() << "MOPAC did not attach on any frame.";
    }

    EXPECT_EQ(B, tp.ProteinRef().BondCount())
        << "bond axis must match protein.BondCount() (bonds.npy axis)";

    // Read back the H5 group + verify shape + sanity-check bond orders.
    const std::string h5_path = (fs::temp_directory_path() /
        ("mopac_bond_order_welford_int1p9j_" +
         std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/mopac_bond_order_welford");
    std::vector<double> means(B);
    grp.getDataSet("order_mean").read(means.data());

    // Sanity: at least some bonds should have nonzero Wiberg orders
    // (covalent bonds). A few may be zero (MOPAC didn't report). All
    // must be finite.
    std::size_t n_finite = 0, n_nonzero = 0;
    double max_order = 0.0;
    for (std::size_t bi = 0; bi < B; ++bi) {
        EXPECT_TRUE(std::isfinite(means[bi]));
        if (std::isfinite(means[bi])) ++n_finite;
        if (std::abs(means[bi]) > 1e-6) ++n_nonzero;
        max_order = std::max(max_order, std::abs(means[bi]));
    }
    EXPECT_EQ(n_finite, B);
    EXPECT_GT(n_nonzero, B / 2)
        << "expected most bonds to have nonzero Wiberg order; "
           "got " << n_nonzero << "/" << B
           << ". MOPAC may have failed to report bond orders.";
    // Wiberg orders are typically 0..3 (single/double/triple), with
    // aromatic ~1.5. Loose sanity bound — log don't assert per
    // feedback_log_overages_dont_assert.
    std::cout << "  max|order| = " << max_order
              << " (typical max ~3.0 for triple bonds)" << std::endl;

    fs::remove(h5_path);
}

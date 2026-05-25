//
// test_dihedral_bin_transition: discipline + integration for
// DihedralBinTransitionTrajectoryResult. AV companion to DihedralTS;
// stats-only (no per-frame bin labels — those are derivable from
// DihedralTS phi/psi/chi via the same binning function copied here).
//

#include "AminoAcidType.h"
#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "DihedralBinTransitionTrajectoryResult.h"
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
    opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::DihedralBinTransitionTrajectoryResult::Create(tp_in);
    });
    config.SetStride(stride);
    return config;
}

}  // namespace


TEST(DihedralBinTransition, Frame0Semantics) {
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

    const auto& tr = tp.Result<nmr::DihedralBinTransitionTrajectoryResult>();
    EXPECT_EQ(tr.NumFrames(), 1u);
}


TEST(DihedralBinTransition, FinalizeIdempotency) {
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

    auto& tr = tp.Result<nmr::DihedralBinTransitionTrajectoryResult>();
    const std::size_t T = tr.NumFrames();
    tr.Finalize(tp, traj);
    EXPECT_EQ(tr.NumFrames(), T);
}


TEST(DihedralBinTransition, H5RoundTrip) {
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

    const auto& tr = tp.Result<nmr::DihedralBinTransitionTrajectoryResult>();
    const std::size_t R = tp.ProteinRef().ResidueCount();

    const std::string h5_path = (fs::temp_directory_path() /
        ("dihedral_bin_trans_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/dihedral_bin_transition"));
    auto grp = reopen.getGroup("/trajectory/dihedral_bin_transition");

    // Per-residue stats datasets.
    for (const std::string& name : {"backbone_transition_count",
                                     "backbone_dominant_region",
                                     "n_frames_observed"}) {
        ASSERT_TRUE(grp.exist(name)) << name;
        const auto dims = grp.getDataSet(name).getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 1u);
        EXPECT_EQ(dims[0], R);
    }
    {
        ASSERT_TRUE(grp.exist("backbone_bin_occupancy"));
        const auto dims =
            grp.getDataSet("backbone_bin_occupancy").getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 2u);
        EXPECT_EQ(dims[0], R);
        EXPECT_EQ(dims[1], 6u);
    }
    for (const std::string& name : {"chi_transition_count",
                                     "chi_dominant_rotamer",
                                     "chi_n_frames_observed"}) {
        ASSERT_TRUE(grp.exist(name)) << name;
        const auto dims = grp.getDataSet(name).getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 2u);
        EXPECT_EQ(dims[0], R);
        EXPECT_EQ(dims[1], 4u);
    }
    {
        ASSERT_TRUE(grp.exist("chi_rotamer_occupancy"));
        const auto dims =
            grp.getDataSet("chi_rotamer_occupancy").getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 3u);
        EXPECT_EQ(dims[0], R);
        EXPECT_EQ(dims[1], 4u);
        EXPECT_EQ(dims[2], 3u);
    }
    EXPECT_TRUE(grp.exist("frame_indices"));
    EXPECT_TRUE(grp.exist("frame_times"));
    EXPECT_TRUE(grp.exist("source_attached_per_frame"));

    // Convention attrs
    std::string legend, gate, source_policy;
    grp.getAttribute("backbone_bin_legend").read(legend);
    grp.getAttribute("transition_gate").read(gate);
    grp.getAttribute("source_attached_policy").read(source_policy);
    EXPECT_NE(legend.find("alphaR"), std::string::npos);
    EXPECT_NE(gate.find("Both prev and curr"), std::string::npos);
    EXPECT_NE(source_policy.find("always_attached"), std::string::npos);

    fs::remove(h5_path);
}


TEST(DihedralBinTransition, Integration1P9J) {
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

    const auto& tr = tp.Result<nmr::DihedralBinTransitionTrajectoryResult>();
    const std::size_t R = tp.ProteinRef().ResidueCount();
    const std::size_t T = tr.NumFrames();
    EXPECT_GE(T, 2u);

    const std::string h5_path = (fs::temp_directory_path() /
        ("dihedral_bin_trans_int_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/dihedral_bin_transition");

    std::vector<std::uint32_t> bb_trans_count, n_obs;
    grp.getDataSet("backbone_transition_count").read(bb_trans_count);
    grp.getDataSet("n_frames_observed").read(n_obs);
    ASSERT_EQ(bb_trans_count.size(), R);
    ASSERT_EQ(n_obs.size(), R);

    std::vector<std::vector<std::uint32_t>> bb_occ;
    grp.getDataSet("backbone_bin_occupancy").read(bb_occ);
    ASSERT_EQ(bb_occ.size(), R);

    // Invariants:
    // (1) Sum of bin occupancy across non-unassigned bins equals
    //     n_frames_observed for each residue.
    // (2) Sum of ALL bin occupancy (including unassigned bin 0) equals
    //     T — every frame accounted for (codex-review-2026-05-19 fix).
    // (3) transition_count <= n_frames_observed (transitions can occur
    //     at most once per consecutive observed-frame pair).
    // (4) For non-N-term, non-C-term residues on 1P9J (linear single-
    //     chain monotonic), n_frames_observed should equal T.
    std::size_t residues_fully_observed = 0;
    std::size_t total_transitions = 0;
    for (std::size_t ri = 0; ri < R; ++ri) {
        std::uint32_t sum_non_unassigned = 0;
        for (std::size_t b = 1; b < 6; ++b) sum_non_unassigned += bb_occ[ri][b];
        EXPECT_EQ(sum_non_unassigned, n_obs[ri])
            << "non-unassigned occupancy sum != n_frames_observed at ri=" << ri;
        std::uint32_t sum_all = bb_occ[ri][0] + sum_non_unassigned;
        EXPECT_EQ(sum_all, T)
            << "total occupancy sum (including unassigned bin 0) != T at ri=" << ri;
        EXPECT_LE(bb_trans_count[ri], n_obs[ri])
            << "transition_count > n_frames_observed at ri=" << ri;
        if (n_obs[ri] == T) ++residues_fully_observed;
        total_transitions += bb_trans_count[ri];
    }
    EXPECT_GE(residues_fully_observed, R - 5u)
        << "Expected most residues to have phi+psi observed every frame";

    std::cout << "DihedralBinTransition: " << R << " residues, " << T
              << " frames; residues_fully_observed=" << residues_fully_observed
              << "; total_backbone_transitions=" << total_transitions << "\n";

    // Chi rotamer invariants: occupancy sum equals
    // chi_n_frames_observed per (residue, chi).
    std::vector<std::vector<std::vector<std::uint32_t>>> chi_occ;
    grp.getDataSet("chi_rotamer_occupancy").read(chi_occ);
    ASSERT_EQ(chi_occ.size(), R);
    std::vector<std::vector<std::uint32_t>> chi_n_obs;
    grp.getDataSet("chi_n_frames_observed").read(chi_n_obs);
    ASSERT_EQ(chi_n_obs.size(), R);

    std::size_t chi_observations = 0;
    for (std::size_t ri = 0; ri < R; ++ri) {
        for (std::size_t k = 0; k < 4; ++k) {
            const auto sum = chi_occ[ri][k][0] + chi_occ[ri][k][1]
                + chi_occ[ri][k][2];
            EXPECT_EQ(sum, chi_n_obs[ri][k])
                << "chi occupancy sum mismatch at ri=" << ri << " k=" << k;
            if (chi_n_obs[ri][k] > 0) ++chi_observations;
        }
    }
    EXPECT_GT(chi_observations, 0u);
    std::cout << "  chi observations (residue × chi pairs with >0 obs): "
              << chi_observations << "\n";

    // dominant_region == argmax of occupancy on residues with >0 obs.
    std::vector<std::uint8_t> dom_region;
    grp.getDataSet("backbone_dominant_region").read(dom_region);
    ASSERT_EQ(dom_region.size(), R);
    for (std::size_t ri = 0; ri < R; ++ri) {
        if (n_obs[ri] == 0) {
            EXPECT_EQ(dom_region[ri], 0u);  // unassigned
            continue;
        }
        // Verify argmax: dom_region[ri] is one of the bins with max count.
        std::uint32_t max_count = 0;
        for (std::size_t b = 0; b < 6; ++b)
            max_count = std::max(max_count, bb_occ[ri][b]);
        EXPECT_EQ(bb_occ[ri][dom_region[ri]], max_count)
            << "dominant_region not argmax at ri=" << ri;
    }

    fs::remove(h5_path);
}

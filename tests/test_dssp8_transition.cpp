//
// test_dssp8_transition: discipline + integration for
// Dssp8TransitionTrajectoryResult. AV companion to Dssp8TimeSeries.
//

#include "CalculatorConfig.h"
#include "DirectionalTestHelpers.h"
#include "DsspResult.h"
#include "Dssp8TransitionTrajectoryResult.h"
#include "GeometryResult.h"
#include "OperationLog.h"
#include "PdbFileReader.h"
#include "Protein.h"
#include "ProteinConformation.h"
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

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <iostream>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
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
    opts.skip_dssp = false;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::DsspResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::Dssp8TransitionTrajectoryResult::Create(tp_in);
    });
    config.SetStride(stride);
    return config;
}

template <typename T>
std::vector<T> ReadDsspTransitionFlat(const std::string& path,
                                      const std::string& dataset) {
    HighFive::File file(path, HighFive::File::ReadOnly);
    auto data_set = file.getDataSet(dataset);
    const auto dims = data_set.getSpace().getDimensions();
    const std::size_t count = std::accumulate(
        dims.begin(), dims.end(), std::size_t{1},
        std::multiplies<std::size_t>());
    std::vector<T> values(count);
    if (!values.empty()) data_set.read(values.data());
    return values;
}

std::string FreshDsspTransitionPath(const std::string& stem) {
    const std::string path = nmr::test::TestEnvironment::TempPath(stem);
    (void)std::remove(path.c_str());
    return path;
}

void ExpectDssp8TransitionDirectionalMetadata(const std::string& path) {
    HighFive::File file(path, HighFive::File::ReadOnly);
    const std::string group_path = "/trajectory/dssp8_transition";
    auto group = file.getGroup(group_path);
    std::string serialization;
    group.getAttribute("dssp_coordinate_serialization").read(serialization);
    EXPECT_EQ(
        serialization,
        "DsspResult writes a temporary PDB with coordinates rounded to "
        "0.001 Angstrom (%8.3f) before cif++/libdssp. The SS8 statistics "
        "have an O(3)-invariant continuum law, but a transformed production "
        "rerun can cross a libdssp category boundary because of this "
        "serialization.")
        << group_path << " @dssp_coordinate_serialization";

    auto expect_strings = [&](const char* dataset_name,
                              const std::vector<std::pair<
                                  std::string, std::string>>& expected) {
        const std::string dataset_path =
            group_path + "/" + dataset_name;
        auto dataset = file.getDataSet(dataset_path);
        for (const auto& [name, value] : expected) {
            std::string actual;
            dataset.getAttribute(name).read(actual);
            EXPECT_EQ(actual, value) << dataset_path << " @" << name;
        }
    };

    expect_strings("source_attached_per_frame", {
        {"units", "dimensionless"},
        {"coordinate_frame", "intrinsic_source_availability"},
        {"parity", "even"},
        {"transformation",
         "exact rotation/translation/reflection-invariant source-attachment "
         "mask"},
        {"validity",
         "1 means DsspResult attached for this frame; 0 breaks every "
         "residue's transition chain and contributes no occupancy"},
    });

    const std::string ss8_transformation =
        "physical O(3)-invariant H/G/I/E/B/T/S/C statistic after the "
        "explicit libdssp PPII-to-coil collapse; translation invariant. "
        "A transformed production rerun may change the statistic when "
        "0.001-A temporary-PDB rounding crosses a libdssp classification "
        "boundary; that is serialization sensitivity, not parity";
    auto ss8_contract = [&](const std::string& units,
                            const std::string& validity) {
        return std::vector<std::pair<std::string, std::string>>{
            {"units", units},
            {"coordinate_frame", "intrinsic_dssp8_category_statistics"},
            {"parity", "even"},
            {"transformation", ss8_transformation},
            {"validity", validity},
        };
    };
    expect_strings("ss8_transition_count", ss8_contract(
        "transitions",
        "physical uint32 count of changes between consecutive observed "
        "SS8 labels; zero is a valid no-transition value. Source-absent "
        "or unmapped frames contribute nothing and break the chain"));
    expect_strings("ss8_dominant", ss8_contract(
        "category",
        "255 iff n_frames_observed is zero; otherwise 0..7 is the "
        "highest-occupancy SS8 code, with ties resolved to the lowest "
        "numeric code"));
    expect_strings("ss8_occupancy", ss8_contract(
        "frame_count",
        "every uint32 entry is a physical per-state frame count; sum "
        "over the state axis equals n_frames_observed for that residue"));
    auto matrix_contract = ss8_contract(
        "transitions",
        "every uint32 entry is a physical consecutive-observed-pair "
        "transition count; diagonal entries are structural zeros because "
        "self-edges are excluded, not unavailable values");
    matrix_contract.emplace_back(
        "structural_zero_entries",
        "all [ri,state,state] diagonal entries are identically zero by "
        "the transition-owner self-edge exclusion policy");
    expect_strings("ss8_transition_matrix", matrix_contract);

    expect_strings("n_frames_observed", {
        {"units", "frame_count"},
        {"coordinate_frame", "intrinsic_observation_count"},
        {"parity", "even"},
        {"transformation",
         "exact rotation/translation/reflection-invariant count of mapped "
         "DSSP residue rows for a fixed typed topology; independent of the "
         "SS8 category value"},
        {"validity",
         "every uint32 entry is a physical count; zero means no mapped "
         "DSSP observation for that residue in any attached frame"},
    });
}

}  // namespace


TEST(Dssp8Transition,
     DirectionalEightStateCollapseRerunProductionDsspO3SerializedH5) {
    using nmr::test::directional::Positions;
    using nmr::test::directional::SeededTransform;

    nmr::test::TestEnvironment::LoadCalculatorConfig();
    const std::string pdb = nmr::test::TestEnvironment::UbqProtonated();
    if (pdb.empty() || !fs::exists(pdb))
        GTEST_SKIP() << "1UBQ protonated fixture unavailable";

    auto source_loaded = nmr::BuildFromProtonatedPdb(pdb);
    ASSERT_TRUE(source_loaded.Ok()) << source_loaded.error;
    auto source_tp = nmr::TrajectoryProtein::CreateForTesting(
        std::move(source_loaded.protein));
    ASSERT_NE(source_tp, nullptr);
    const auto source_positions =
        source_tp->CanonicalConformation().Positions();

    // A second, deterministic molecular frame displaces every backbone
    // carbonyl oxygen far from its source position.  This leaves the typed
    // protein/topology and peptide-chain coordinates intact while removing
    // the distance-derived carbonyl H-bond pattern that owns DSSP's H/E
    // assignments.  Both frames below still run through production
    // DsspResult::Compute (including its PDB/libdssp boundary); no synthetic
    // DsspResidue is injected into the transition owner.
    auto disrupted_positions = source_positions;
    constexpr double kCarbonylDxA = 80.0;
    constexpr double kCarbonylDyA = -60.0;
    constexpr double kCarbonylDzA = 40.0;
    std::size_t displaced_carbonyls = 0;
    for (std::size_t ai = 0; ai < source_tp->ProteinRef().AtomCount(); ++ai) {
        if (source_tp->ProteinRef().AtomAt(ai).pdb_atom_name != "O") continue;
        disrupted_positions[ai] +=
            nmr::Vec3(kCarbonylDxA, kCarbonylDyA, kCarbonylDzA);
        ++displaced_carbonyls;
    }
    ASSERT_GT(displaced_carbonyls, 50u);
    const std::vector<std::vector<nmr::Vec3>> source_frames = {
        source_positions, disrupted_positions};

    auto source_result =
        nmr::Dssp8TransitionTrajectoryResult::Create(*source_tp);
    ASSERT_NE(source_result, nullptr);
    nmr::Trajectory dummy("", "", "");
    for (std::size_t frame = 0; frame < source_frames.size(); ++frame) {
        auto source_conf = source_tp->TickConformation(source_frames[frame]);
        auto source_dssp = nmr::DsspResult::Compute(*source_conf);
        ASSERT_NE(source_dssp, nullptr) << "source frame=" << frame;
        ASSERT_TRUE(source_conf->AttachResult(std::move(source_dssp)));
        source_result->Compute(
            *source_conf, *source_tp, dummy, 83 + frame, 5.25 + frame);
    }
    source_result->Finalize(*source_tp, dummy);
    const std::string source_path =
        FreshDsspTransitionPath("dssp8_transition_directional_source.h5");
    {
        HighFive::File file(source_path, HighFive::File::Truncate);
        source_result->WriteH5Group(*source_tp, file);
    }
    ExpectDssp8TransitionDirectionalMetadata(source_path);

    const std::string group = "/trajectory/dssp8_transition/";
    const auto source_count = ReadDsspTransitionFlat<std::uint32_t>(
        source_path, group + "ss8_transition_count");
    const auto source_dominant = ReadDsspTransitionFlat<std::uint8_t>(
        source_path, group + "ss8_dominant");
    const auto source_observed = ReadDsspTransitionFlat<std::uint32_t>(
        source_path, group + "n_frames_observed");
    const auto source_occupancy = ReadDsspTransitionFlat<std::uint32_t>(
        source_path, group + "ss8_occupancy");
    const auto source_matrix = ReadDsspTransitionFlat<std::uint32_t>(
        source_path, group + "ss8_transition_matrix");
    ASSERT_EQ(source_count.size(), source_tp->ProteinRef().ResidueCount());
    ASSERT_GT(std::count(source_dominant.begin(), source_dominant.end(), 7u),
              0);
    const std::uint64_t source_transition_total = std::accumulate(
        source_count.begin(), source_count.end(), std::uint64_t{0});
    const std::uint64_t source_matrix_total = std::accumulate(
        source_matrix.begin(), source_matrix.end(), std::uint64_t{0});
    ASSERT_GT(source_transition_total, 0u)
        << "two production libdssp frames did not produce an SS8 change";
    ASSERT_EQ(source_matrix_total, source_transition_total);

    constexpr std::uint64_t kTransformSeed = 0x4453535038545201ULL;
    for (const bool improper : {false, true}) {
        auto moved_loaded = nmr::BuildFromProtonatedPdb(pdb);
        ASSERT_TRUE(moved_loaded.Ok()) << moved_loaded.error;
        auto moved_tp = nmr::TrajectoryProtein::CreateForTesting(
            std::move(moved_loaded.protein));
        ASSERT_NE(moved_tp, nullptr);
        const auto transform = SeededTransform(kTransformSeed, improper);
        auto moved_result =
            nmr::Dssp8TransitionTrajectoryResult::Create(*moved_tp);
        ASSERT_NE(moved_result, nullptr);
        for (std::size_t frame = 0; frame < source_frames.size(); ++frame) {
            auto moved_conf = moved_tp->TickConformation(
                Positions(transform, source_frames[frame]));
            auto moved_dssp = nmr::DsspResult::Compute(*moved_conf);
            ASSERT_NE(moved_dssp, nullptr)
                << "improper=" << improper << " frame=" << frame;
            ASSERT_TRUE(moved_conf->AttachResult(std::move(moved_dssp)));
            moved_result->Compute(
                *moved_conf, *moved_tp, dummy,
                83 + frame, 5.25 + frame);
        }
        moved_result->Finalize(*moved_tp, dummy);
        const std::string moved_path = FreshDsspTransitionPath(
            improper ? "dssp8_transition_directional_improper.h5" :
                       "dssp8_transition_directional_proper.h5");
        {
            HighFive::File file(moved_path, HighFive::File::Truncate);
            moved_result->WriteH5Group(*moved_tp, file);
        }
        ExpectDssp8TransitionDirectionalMetadata(moved_path);

        // The transition owner consumes the eight-state code, not the
        // reflection-sensitive PPII flag. Its documented P -> C=7 collapse
        // therefore makes every serialized transition/occupancy category
        // O(3)-invariant, including under the real improper libdssp rerun.
        EXPECT_EQ(ReadDsspTransitionFlat<std::uint32_t>(
                      moved_path, group + "ss8_transition_count"),
                  source_count);
        EXPECT_EQ(ReadDsspTransitionFlat<std::uint8_t>(
                      moved_path, group + "ss8_dominant"),
                  source_dominant);
        EXPECT_EQ(ReadDsspTransitionFlat<std::uint32_t>(
                      moved_path, group + "n_frames_observed"),
                  source_observed);
        EXPECT_EQ(ReadDsspTransitionFlat<std::uint32_t>(
                      moved_path, group + "ss8_occupancy"),
                  source_occupancy);
        EXPECT_EQ(ReadDsspTransitionFlat<std::uint32_t>(
                      moved_path, group + "ss8_transition_matrix"),
                  source_matrix);
        EXPECT_EQ(ReadDsspTransitionFlat<std::size_t>(
                      moved_path, group + "frame_indices"),
                  ReadDsspTransitionFlat<std::size_t>(
                      source_path, group + "frame_indices"));
        EXPECT_EQ(ReadDsspTransitionFlat<double>(
                      moved_path, group + "frame_times"),
                  ReadDsspTransitionFlat<double>(
                      source_path, group + "frame_times"));
        EXPECT_EQ(ReadDsspTransitionFlat<std::uint8_t>(
                      moved_path, group + "source_attached_per_frame"),
                  ReadDsspTransitionFlat<std::uint8_t>(
                      source_path, group + "source_attached_per_frame"));
        EXPECT_EQ(std::remove(moved_path.c_str()), 0);
    }
    EXPECT_EQ(std::remove(source_path.c_str()), 0);
}


TEST(Dssp8Transition, H5RoundTripSingleFrame) {
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

    const auto& tr = tp.Result<nmr::Dssp8TransitionTrajectoryResult>();
    const std::size_t R = tp.ProteinRef().ResidueCount();
    EXPECT_EQ(tr.NumFrames(), 1u);

    const std::string h5_path = (fs::temp_directory_path() /
        ("dssp8_trans_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/dssp8_transition"));
    auto grp = reopen.getGroup("/trajectory/dssp8_transition");

    for (const std::string& name : {"ss8_transition_count",
                                     "ss8_dominant",
                                     "n_frames_observed"}) {
        ASSERT_TRUE(grp.exist(name)) << name;
        const auto dims = grp.getDataSet(name).getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 1u);
        EXPECT_EQ(dims[0], R);
    }
    {
        const auto dims =
            grp.getDataSet("ss8_occupancy").getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 2u);
        EXPECT_EQ(dims[0], R);
        EXPECT_EQ(dims[1], 8u);
    }
    {
        const auto dims =
            grp.getDataSet("ss8_transition_matrix").getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 3u);
        EXPECT_EQ(dims[0], R);
        EXPECT_EQ(dims[1], 8u);
        EXPECT_EQ(dims[2], 8u);
    }

    // Single-frame run: zero transitions, occupancy sums == n_obs.
    std::vector<std::uint32_t> trans_count;
    grp.getDataSet("ss8_transition_count").read(trans_count);
    for (auto c : trans_count) EXPECT_EQ(c, 0u);

    fs::remove(h5_path);
}


TEST(Dssp8Transition, Integration1P9J) {
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

    const auto& tr = tp.Result<nmr::Dssp8TransitionTrajectoryResult>();
    const std::size_t R = tp.ProteinRef().ResidueCount();
    const std::size_t T = tr.NumFrames();
    EXPECT_GE(T, 2u);

    const std::string h5_path = (fs::temp_directory_path() /
        ("dssp8_trans_int_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/dssp8_transition");

    // Invariants:
    // (1) Sum of ss8_occupancy across all 8 states == n_frames_observed.
    // (2) Transition matrix diagonal is zero (self-transitions excluded).
    // (3) Sum of off-diagonal transition matrix == ss8_transition_count.
    std::vector<std::uint32_t> trans_count, n_obs;
    grp.getDataSet("ss8_transition_count").read(trans_count);
    grp.getDataSet("n_frames_observed").read(n_obs);
    std::vector<std::vector<std::uint32_t>> occ;
    grp.getDataSet("ss8_occupancy").read(occ);
    std::vector<std::vector<std::vector<std::uint32_t>>> mat;
    grp.getDataSet("ss8_transition_matrix").read(mat);

    std::size_t total_trans = 0, total_obs = 0;
    for (std::size_t ri = 0; ri < R; ++ri) {
        std::uint32_t occ_sum = 0;
        for (std::size_t s = 0; s < 8; ++s) occ_sum += occ[ri][s];
        EXPECT_EQ(occ_sum, n_obs[ri])
            << "occupancy sum != n_frames_observed at ri=" << ri;

        std::uint32_t mat_off_diag_sum = 0;
        for (std::size_t p = 0; p < 8; ++p) {
            for (std::size_t c = 0; c < 8; ++c) {
                if (p == c) {
                    EXPECT_EQ(mat[ri][p][c], 0u)
                        << "diagonal nonzero at ri=" << ri
                        << " state=" << p;
                } else {
                    mat_off_diag_sum += mat[ri][p][c];
                }
            }
        }
        EXPECT_EQ(mat_off_diag_sum, trans_count[ri])
            << "matrix off-diagonal sum != transition_count at ri=" << ri;

        total_trans += trans_count[ri];
        total_obs   += n_obs[ri];
    }
    std::cout << "Dssp8Transition: " << R << " residues, " << T
              << " frames; total_obs=" << total_obs
              << " total_transitions=" << total_trans << "\n";

    fs::remove(h5_path);
}

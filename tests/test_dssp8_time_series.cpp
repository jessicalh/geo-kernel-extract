//
// test_dssp8_time_series: discipline + integration for
// Dssp8TimeSeriesTrajectoryResult.
//

#include "CalculatorConfig.h"
#include "DirectionalTestHelpers.h"
#include "DsspResult.h"
#include "Dssp8TimeSeriesTrajectoryResult.h"
#include "EnrichmentResult.h"
#include "GeometryResult.h"
#include "OperationLog.h"
#include "PdbFileReader.h"
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
#include <unistd.h>
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

template <typename T>
std::vector<T> ReadDsspH5Flat(
        const std::string& path,
        const std::string& dataset,
        std::vector<std::size_t>* dimensions = nullptr) {
    HighFive::File file(path, HighFive::File::ReadOnly);
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

std::string FreshDsspH5Path(const std::string& stem) {
    const std::string path = nmr::test::TestEnvironment::TempPath(stem);
    (void)std::remove(path.c_str());
    return path;
}

void ExpectDsspSameIncludingNaN(const std::vector<double>& actual,
                                const std::vector<double>& expected,
                                double abs_tolerance,
                                const std::string& dataset) {
    ASSERT_EQ(actual.size(), expected.size()) << dataset;
    for (std::size_t i = 0; i < expected.size(); ++i) {
        SCOPED_TRACE(dataset + " index=" + std::to_string(i));
        if (std::isnan(expected[i])) {
            EXPECT_TRUE(std::isnan(actual[i]));
        } else {
            ASSERT_TRUE(std::isfinite(actual[i]));
            EXPECT_NEAR(actual[i], expected[i], abs_tolerance);
        }
    }
}

void ExpectDssp8DirectionalMetadata(const std::string& path) {
    HighFive::File file(path, HighFive::File::ReadOnly);
    const std::string group_path = "/trajectory/dssp8_time_series";
    auto group = file.getGroup(group_path);
    std::string serialization;
    group.getAttribute("dssp_coordinate_serialization").read(serialization);
    EXPECT_EQ(
        serialization,
        "DsspResult writes a temporary PDB with coordinates rounded to "
        "0.001 Angstrom (%8.3f) before cif++/libdssp. Physical O(3) laws "
        "below are continuum laws; transformed production reruns can cross "
        "a libdssp decision boundary or show bounded continuous-value drift "
        "because of this serialization.")
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
         "1 means DsspResult attached for this frame; 0 means all "
         "per-residue payloads in that frame carry unavailable sentinels"},
    });
    expect_strings("ss8_code", {
        {"units", "category"},
        {"coordinate_frame", "intrinsic_dssp8_category"},
        {"parity", "even"},
        {"transformation",
         "physical O(3)-invariant H/G/I/E/B/T/S/C category after the "
         "explicit libdssp PPII-to-coil collapse; translation invariant. "
         "A transformed production rerun may change a boundary category "
         "when 0.001-A temporary-PDB rounding crosses a libdssp decision "
         "threshold; that is serialization sensitivity, not parity"},
        {"validity",
         "255 means no DSSP observation for this residue/frame; "
         "source_attached_per_frame gates whole-frame availability and "
         "ss8_code != 255 is the per-residue observation mask"},
    });
    expect_strings("ppii_flag", {
        {"units", "category"},
        {"coordinate_frame", "intrinsic_signed_dssp_chirality_category"},
        {"parity", "mixed"},
        {"transformation",
         "translation/proper-rotation invariant PPII category; reflection "
         "negates the signed phi/psi geometry used by libdssp, so the "
         "chirality-conditioned predicate has no fixed improper-transform "
         "map. The 0.001-A temporary-PDB rounding can affect values at a "
         "classification boundary"},
        {"validity",
         "1=PPII, 0=observed non-PPII, 255=no DSSP observation; "
         "source_attached_per_frame gates whole-frame availability"},
    });

    const std::vector<std::pair<std::string, std::string>> partner_contract = {
        {"units", "residue_index"},
        {"coordinate_frame", "intrinsic_protein_residue_identity"},
        {"parity", "even"},
        {"transformation",
         "O(3)-invariant residue identity of the distance-derived DSSP "
         "hydrogen-bond partner; translation invariant"},
        {"validity",
         "nonnegative values are protein residue indices; -1 means no "
         "mapped partner or no DSSP observation. Use ss8_code != 255 for "
         "the per-residue observation mask and source_attached_per_frame "
         "for whole-frame availability"},
    };
    expect_strings("hbond_acceptor_partner", partner_contract);
    expect_strings("hbond_donor_partner", partner_contract);

    const std::vector<std::pair<std::string, std::string>> energy_contract = {
        {"units", "kcal/mol"},
        {"coordinate_frame", "intrinsic_dssp_hbond_energy"},
        {"parity", "even"},
        {"transformation",
         "continuum rotation/translation/reflection-invariant DSSP "
         "electrostatic H-bond energy scalar; transformed production "
         "reruns use an absolute 5e-3 kcal/mol envelope because the "
         "temporary PDB rounds coordinates to 0.001 Angstrom"},
        {"covariance_tolerance",
         "absolute 5e-3 kcal/mol at the serialized DsspResult/libdssp "
         "boundary; recorded forcing-transform maximum 4e-3 kcal/mol"},
        {"validity",
         "finite only for a mapped partner on an observed residue; NaN "
         "means no mapped partner or no DSSP observation. Partner >= 0 "
         "and ss8_code != 255 are the corresponding validity gates"},
    };
    expect_strings("hbond_acceptor_energy", energy_contract);
    expect_strings("hbond_donor_energy", energy_contract);

    expect_strings("residue_index_per_atom", {
        {"units", "residue_index"},
        {"coordinate_frame", "intrinsic_topology_lookup"},
        {"parity", "even"},
        {"transformation",
         "exact rotation/translation/reflection-invariant atom-to-residue "
         "topology lookup"},
        {"validity",
         "ordinary typed proteins populate every atom row with a "
         "nonnegative protein residue index"},
    });
}

}  // namespace


TEST(Dssp8TimeSeries,
     DirectionalCategoricalRerunProductionDsspO3SerializedH5) {
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
    auto source_conf = source_tp->TickConformation(source_positions);
    auto source_dssp = nmr::DsspResult::Compute(*source_conf);
    ASSERT_NE(source_dssp, nullptr);
    ASSERT_TRUE(source_conf->AttachResult(std::move(source_dssp)));
    auto source_result =
        nmr::Dssp8TimeSeriesTrajectoryResult::Create(*source_tp);
    ASSERT_NE(source_result, nullptr);
    nmr::Trajectory dummy("", "", "");
    source_result->Compute(*source_conf, *source_tp, dummy, 73, 4.5);
    source_result->Finalize(*source_tp, dummy);
    const std::string source_path =
        FreshDsspH5Path("dssp8_directional_source.h5");
    {
        HighFive::File file(source_path, HighFive::File::Truncate);
        source_result->WriteH5Group(*source_tp, file);
    }
    ExpectDssp8DirectionalMetadata(source_path);

    const std::string group = "/trajectory/dssp8_time_series/";
    const std::size_t residue_count = source_tp->ProteinRef().ResidueCount();
    std::vector<std::size_t> ss_dims;
    const auto source_ss = ReadDsspH5Flat<std::uint8_t>(
        source_path, group + "ss8_code", &ss_dims);
    const auto source_ppii = ReadDsspH5Flat<std::uint8_t>(
        source_path, group + "ppii_flag");
    const auto source_acceptor_partner = ReadDsspH5Flat<std::int32_t>(
        source_path, group + "hbond_acceptor_partner");
    const auto source_donor_partner = ReadDsspH5Flat<std::int32_t>(
        source_path, group + "hbond_donor_partner");
    const auto source_acceptor_energy = ReadDsspH5Flat<double>(
        source_path, group + "hbond_acceptor_energy");
    const auto source_donor_energy = ReadDsspH5Flat<double>(
        source_path, group + "hbond_donor_energy");
    const auto source_frame_indices = ReadDsspH5Flat<std::size_t>(
        source_path, group + "frame_indices");
    const auto source_frame_times = ReadDsspH5Flat<double>(
        source_path, group + "frame_times");
    ASSERT_EQ(ss_dims,
              (std::vector<std::size_t>{residue_count, 1u}));
    ASSERT_EQ(source_ppii.size(), residue_count);
    ASSERT_GT(std::count(source_ppii.begin(), source_ppii.end(), 1u), 0)
        << "1UBQ must exercise the reflection-sensitive PPII category";
    ASSERT_EQ(source_frame_indices, (std::vector<std::size_t>{73u}));
    ASSERT_EQ(source_frame_times, (std::vector<double>{4.5}));

    constexpr std::uint64_t kTransformSeed = 0x4453535038545301ULL;
    // libdssp's electrostatic H-bond energy depends only on interatomic
    // distances, but the production boundary writes a temporary PDB with
    // 0.001-A coordinate quantization before libdssp reads it. The recorded
    // maximum energy drift for this forcing transform is 0.004 kcal/mol;
    // 0.005 is the explicit serialization-bound tolerance.
    constexpr double kEnergyAbsTolerance = 5.0e-3;
    std::size_t improper_ppii_changed = 0;
    std::size_t improper_ppii_unchanged_observed = 0;
    std::size_t improper_ss_changed = 0;
    std::size_t improper_ss_unchanged_observed = 0;

    for (const bool improper : {false, true}) {
        auto moved_loaded = nmr::BuildFromProtonatedPdb(pdb);
        ASSERT_TRUE(moved_loaded.Ok()) << moved_loaded.error;
        auto moved_tp = nmr::TrajectoryProtein::CreateForTesting(
            std::move(moved_loaded.protein));
        ASSERT_NE(moved_tp, nullptr);
        const auto transform = SeededTransform(kTransformSeed, improper);
        auto moved_conf = moved_tp->TickConformation(
            Positions(transform, source_positions));
        auto moved_dssp = nmr::DsspResult::Compute(*moved_conf);
        ASSERT_NE(moved_dssp, nullptr);
        ASSERT_TRUE(moved_conf->AttachResult(std::move(moved_dssp)));
        auto moved_result =
            nmr::Dssp8TimeSeriesTrajectoryResult::Create(*moved_tp);
        ASSERT_NE(moved_result, nullptr);
        moved_result->Compute(*moved_conf, *moved_tp, dummy, 73, 4.5);
        moved_result->Finalize(*moved_tp, dummy);
        const std::string moved_path = FreshDsspH5Path(
            improper ? "dssp8_directional_improper.h5" :
                       "dssp8_directional_proper.h5");
        {
            HighFive::File file(moved_path, HighFive::File::Truncate);
            moved_result->WriteH5Group(*moved_tp, file);
        }
        ExpectDssp8DirectionalMetadata(moved_path);

        const auto moved_ss = ReadDsspH5Flat<std::uint8_t>(
            moved_path, group + "ss8_code");
        const auto moved_ppii = ReadDsspH5Flat<std::uint8_t>(
            moved_path, group + "ppii_flag");
        ASSERT_EQ(moved_ss.size(), source_ss.size());
        ASSERT_EQ(moved_ppii.size(), source_ppii.size());
        // The ordinary DSSP8 code collapses libdssp's P extension into the
        // coil code (7), so its serialized eight-state category remains O(3)
        // invariant. The separate PPII flag retains the signed distinction.
        EXPECT_EQ(moved_ss, source_ss);
        if (!improper) {
            EXPECT_EQ(moved_ppii, source_ppii);
        } else {
            for (std::size_t ri = 0; ri < residue_count; ++ri) {
                const bool source_observed =
                    source_ss[ri] !=
                    nmr::Dssp8TimeSeriesTrajectoryResult::kSSUnassigned;
                const bool moved_observed =
                    moved_ss[ri] !=
                    nmr::Dssp8TimeSeriesTrajectoryResult::kSSUnassigned;
                if (source_observed && moved_observed) {
                    if (moved_ss[ri] == source_ss[ri])
                        ++improper_ss_unchanged_observed;
                    else
                        ++improper_ss_changed;
                    if (moved_ppii[ri] == source_ppii[ri])
                        ++improper_ppii_unchanged_observed;
                    else
                        ++improper_ppii_changed;
                }
            }
        }

        // DSSP H-bond partners and energies are distance-derived and remain
        // invariant under both proper and improper rigid transforms.
        EXPECT_EQ(ReadDsspH5Flat<std::int32_t>(
                      moved_path, group + "hbond_acceptor_partner"),
                  source_acceptor_partner);
        EXPECT_EQ(ReadDsspH5Flat<std::int32_t>(
                      moved_path, group + "hbond_donor_partner"),
                  source_donor_partner);
        ExpectDsspSameIncludingNaN(
            ReadDsspH5Flat<double>(
                moved_path, group + "hbond_acceptor_energy"),
            source_acceptor_energy, kEnergyAbsTolerance,
            "hbond_acceptor_energy");
        ExpectDsspSameIncludingNaN(
            ReadDsspH5Flat<double>(
                moved_path, group + "hbond_donor_energy"),
            source_donor_energy, kEnergyAbsTolerance,
            "hbond_donor_energy");
        EXPECT_EQ(ReadDsspH5Flat<std::int32_t>(
                      moved_path, group + "residue_index_per_atom"),
                  ReadDsspH5Flat<std::int32_t>(
                      source_path, group + "residue_index_per_atom"));
        EXPECT_EQ(ReadDsspH5Flat<std::uint8_t>(
                      moved_path, group + "source_attached_per_frame"),
                  ReadDsspH5Flat<std::uint8_t>(
                      source_path, group + "source_attached_per_frame"));
        EXPECT_EQ(ReadDsspH5Flat<std::size_t>(
                      moved_path, group + "frame_indices"),
                  source_frame_indices);
        EXPECT_EQ(ReadDsspH5Flat<double>(
                      moved_path, group + "frame_times"),
                  source_frame_times);
        EXPECT_EQ(std::remove(moved_path.c_str()), 0);
    }

    // Proper rotations preserve libdssp's signed torsions and categories.
    // Reflection negates the signed torsions used by the PPII predicate:
    // some categories change while others remain fixed. This directly
    // refutes both an invariant claim and a global-complement claim for the
    // serialized PPII flag. ss8_code is invariant because P and C deliberately
    // share code 7 at this eight-state boundary.
    EXPECT_GT(improper_ppii_changed, 0u);
    EXPECT_GT(improper_ppii_unchanged_observed, 0u);
    EXPECT_EQ(improper_ss_changed, 0u);
    EXPECT_GT(improper_ss_unchanged_observed, 0u);
    EXPECT_EQ(std::remove(source_path.c_str()), 0);
}


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


TEST(Dssp8TimeSeries, SyntheticPpiiFrameThenSourceAbsent) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    const std::size_t N = tp.AtomCount();
    const std::size_t R = tp.ProteinRef().ResidueCount();
    ASSERT_GT(R, 0u);
    auto tr = nmr::Dssp8TimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    std::vector<nmr::Vec3> positions(N, nmr::Vec3::Zero());

    const std::size_t p_residue = 0;
    {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic PPII frame");
        std::vector<nmr::DsspResidue> residues(R);
        residues[p_residue].observed = true;
        residues[p_residue].secondary_structure = 'P';
        ASSERT_TRUE(conf->AttachResult(
            nmr::DsspResult::CreateForTesting(std::move(residues))));
        tr->Compute(*conf, tp, traj, 0, 0.0);
    }
    {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic source-absent frame");
        tr->Compute(*conf, tp, traj, 1, 1.0);
    }
    tr->Finalize(tp, traj);

    const std::string h5_path = (fs::temp_directory_path() /
        ("dssp8_ts_ppii_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/dssp8_time_series"));
    auto grp = reopen.getGroup("/trajectory/dssp8_time_series");

    std::vector<std::vector<std::uint8_t>> ss;
    std::vector<std::vector<std::uint8_t>> ppii;
    grp.getDataSet("ss8_code").read(ss);
    grp.getDataSet("ppii_flag").read(ppii);
    ASSERT_EQ(ss.size(), R);
    ASSERT_EQ(ppii.size(), R);
    ASSERT_EQ(ss[p_residue].size(), 2u);
    ASSERT_EQ(ppii[p_residue].size(), 2u);

    EXPECT_EQ(ss[p_residue][0], 7u);
    EXPECT_EQ(ppii[p_residue][0], 1u);
    EXPECT_EQ(ss[p_residue][1],
              nmr::Dssp8TimeSeriesTrajectoryResult::kSSUnassigned);
    EXPECT_EQ(ppii[p_residue][1],
              nmr::Dssp8TimeSeriesTrajectoryResult::kSSUnassigned);

    fs::remove(h5_path);
}

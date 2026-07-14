//
// test_sasa_time_series: discipline + integration for
// SasaTimeSeriesTrajectoryResult. Scalar-double DenseBuffer pattern.
//

#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "SasaResult.h"

#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "DenseBuffer.h"
#include "SasaTimeSeriesTrajectoryResult.h"
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
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <unistd.h>
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

nmr::ProteinConformation& BuildOneAtomSasaFixture(nmr::Protein& protein) {
    nmr::Residue residue;
    residue.type = nmr::AminoAcid::ALA;
    residue.sequence_number = 1;
    const size_t ri = protein.AddResidue(std::move(residue));
    auto atom = nmr::Atom::Create(nmr::Element::C);
    atom->residue_index = ri;
    const size_t ai = protein.AddAtom(std::move(atom));
    protein.MutableResidueAt(ri).atom_indices.push_back(ai);
    return protein.AddConformation({nmr::Vec3(0.0, 0.0, 0.0)},
                                   "one atom sasa");
}

bool AttachGeometryAndSpatial(nmr::ProteinConformation& conf) {
    auto geometry = nmr::GeometryResult::Compute(conf);
    if (!geometry || !conf.AttachResult(std::move(geometry))) return false;
    auto spatial = nmr::SpatialIndexResult::Compute(conf);
    if (!spatial || !conf.AttachResult(std::move(spatial))) return false;
    return true;
}

void RunBadSasaPointConfigInChild(const char* value) {
    const fs::path toml = fs::temp_directory_path() /
        ("sasa_bad_points_" + std::to_string(::getpid()) + ".toml");
    {
        std::ofstream out(toml);
        out << "sasa_n_points = " << value << "\n";
    }
    nmr::CalculatorConfig::Load(toml.string());
    nmr::Protein protein;
    auto& conf = BuildOneAtomSasaFixture(protein);
    if (!AttachGeometryAndSpatial(conf)) {
        fs::remove(toml);
        _exit(2);
    }
    auto sasa = nmr::SasaResult::Compute(conf);
    fs::remove(toml);
    _exit(sasa == nullptr ? 0 : 1);
}
}  // namespace


TEST(SasaTimeSeries, sasa_n_points_validation) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    {
        nmr::Protein protein;
        auto& conf = BuildOneAtomSasaFixture(protein);
        ASSERT_TRUE(AttachGeometryAndSpatial(conf));
        auto sasa = nmr::SasaResult::Compute(conf);
        ASSERT_NE(sasa, nullptr);
        EXPECT_TRUE(std::isfinite(conf.AtomAt(0).atom_sasa));
        EXPECT_GE(conf.AtomAt(0).atom_sasa, 0.0);
    }

    EXPECT_EXIT({ RunBadSasaPointConfigInChild("0"); },
                ::testing::ExitedWithCode(0), "");
    EXPECT_EXIT({ RunBadSasaPointConfigInChild("-1"); },
                ::testing::ExitedWithCode(0), "");
    EXPECT_EXIT({ RunBadSasaPointConfigInChild("92.5"); },
                ::testing::ExitedWithCode(0), "");
    EXPECT_EXIT({ RunBadSasaPointConfigInChild("nan"); },
                ::testing::ExitedWithCode(0), "");
}


TEST(SasaTimeSeries, SyntheticFourFrames) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    const size_t N = tp.AtomCount();
    auto tr = nmr::SasaTimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    constexpr size_t kFrames = 4;
    std::vector<nmr::Vec3> positions(N, nmr::Vec3::Zero());
    for (size_t t = 0; t < kFrames; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic frame");
        for (size_t i = 0; i < N; ++i) {
            conf->MutableAtomAt(i).atom_sasa = static_cast<double>(i) + t * 100.0;
        }
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);

    auto* buf = tp.GetDenseBuffer<double>(std::type_index(typeid(
        nmr::SasaTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);
    EXPECT_EQ(buf->StridePerAtom(), kFrames);
    for (size_t i : {size_t(0), N / 2, N - 1})
        for (size_t t = 0; t < kFrames; ++t)
            EXPECT_DOUBLE_EQ(buf->At(i, t), static_cast<double>(i) + t * 100.0);

    const std::string h5_path = (fs::temp_directory_path() /
        ("sasa_ts_unit_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/sasa_time_series"));
    fs::remove(h5_path);
}


TEST(SasaTimeSeries, FiniteGridDirectionalMetadataWithoutFleetFixture) {
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

    auto tr = nmr::SasaTimeSeriesTrajectoryResult::Create(tp);
    ASSERT_NE(tr, nullptr);
    nmr::Trajectory dummy("", "", "");
    std::vector<nmr::Vec3> positions(tp.AtomCount(), nmr::Vec3::Zero());
    nmr::ProteinConformation conf(
        &tp.ProteinRef(), positions, "synthetic SASA metadata frame");
    for (std::size_t atom = 0; atom < tp.AtomCount(); ++atom)
        conf.MutableAtomAt(atom).atom_sasa = 0.25 * atom;
    tr->Compute(conf, tp, dummy, 9u, 1.5);
    tr->Finalize(tp, dummy);

    const std::string h5_path = (fs::temp_directory_path() /
        ("sasa_ts_directional_metadata_" + std::to_string(::getpid()) +
         ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        tr->WriteH5Group(tp, file);
    }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/sasa_time_series");
    std::string layout, parity, frame, transformation, scope, units;
    grp.getAttribute("irrep_layout").read(layout);
    grp.getAttribute("parity").read(parity);
    grp.getAttribute("coordinate_frame").read(frame);
    grp.getAttribute("transformation").read(transformation);
    grp.getAttribute("directional_metadata_scope").read(scope);
    grp.getDataSet("sasa").getAttribute("units").read(units);
    EXPECT_EQ(layout, "raw_scalar_no_exact_o3_irrep");
    EXPECT_EQ(parity, "mixed");
    EXPECT_EQ(frame, "lab_fixed_fibonacci_sampling_grid");
    EXPECT_EQ(transformation,
              "continuum rotation-invariant scalar; live finite lab-fixed Fibonacci "
              "estimator has no exact O(3) law and is only approximately invariant "
              "within the recorded covariance-test envelope");
    EXPECT_EQ(scope, "sasa dataset only");
    EXPECT_EQ(units, "Angstrom^2");
    fs::remove(h5_path);
}


TEST(SasaTimeSeries, Frame0Semantics) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true; opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::SasaResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::SasaTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    EXPECT_EQ(traj.FrameCount(), 1u);
}


TEST(SasaTimeSeries, FinalizeIdempotency) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true; opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::SasaResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::SasaTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    auto* buf_first = tp.GetDenseBuffer<double>(std::type_index(typeid(
        nmr::SasaTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf_first, nullptr);
    const std::size_t T_first = buf_first->StridePerAtom();
    auto& tr = tp.Result<nmr::SasaTimeSeriesTrajectoryResult>();
    tr.Finalize(tp, traj);
    auto* buf_second = tp.GetDenseBuffer<double>(std::type_index(typeid(
        nmr::SasaTimeSeriesTrajectoryResult)));
    EXPECT_EQ(buf_second->StridePerAtom(), T_first);
}


TEST(SasaTimeSeries, H5RoundTrip) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true; opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::SasaResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::SasaTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::SasaTimeSeriesTrajectoryResult>();
    const std::string h5_path = (fs::temp_directory_path() /
        ("sasa_ts_h5_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/sasa_time_series");
    auto ds = grp.getDataSet("sasa");
    const auto dims = ds.getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 2u);
    EXPECT_EQ(dims[0], tp.AtomCount());
    EXPECT_EQ(dims[1], 1u);

    std::string parity, units, layout;
    grp.getAttribute("parity").read(parity);
    grp.getAttribute("units").read(units);
    grp.getAttribute("irrep_layout").read(layout);
    EXPECT_EQ(parity, "mixed");
    EXPECT_EQ(units, "Angstrom^2");
    EXPECT_EQ(layout, "raw_scalar_no_exact_o3_irrep");
    fs::remove(h5_path);
}


TEST(SasaTimeSeries, Integration1P9J) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true; opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::SasaResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::SasaTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(300);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    auto* buf = tp.GetDenseBuffer<double>(std::type_index(typeid(
        nmr::SasaTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);

    double max_sasa = 0.0;
    size_t populated = 0;
    for (size_t i = 0; i < buf->AtomCount(); ++i) {
        for (size_t t = 0; t < buf->StridePerAtom(); ++t) {
            const double v = buf->At(i, t);
            EXPECT_TRUE(std::isfinite(v));
            EXPECT_GE(v, 0.0);  // SASA cannot be negative
            max_sasa = std::max(max_sasa, v);
        }
        if (buf->At(i, 0) > 1e-12) ++populated;
    }
    EXPECT_GT(max_sasa, 1e-3) << "SASA all zero — calc not firing";
    // SASA per atom is bounded by 4πr² of a ~3 Å sphere (r_vdW + r_probe) → ~113 Å²/atom max.
    EXPECT_LT(max_sasa, 200.0) << "SASA exceeds physical bound";
    std::cout << "SasaTimeSeries max=" << max_sasa << " A^2, populated=" << populated
              << "/" << buf->AtomCount() << "\n";
}

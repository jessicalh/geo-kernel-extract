//
// test_dihedral_time_series: discipline + integration for
// DihedralTimeSeriesTrajectoryResult. Per-residue per-frame
// phi/psi/omega/chi timelines + Ramachandran-region classification +
// static per-residue masks. Self-contained TR (positions + topology
// only — no source ConformationResult dependency).
//

#include "AminoAcidType.h"
#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "DihedralTimeSeriesTrajectoryResult.h"
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
#include <limits>
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

nmr::RunConfiguration BuildConfig(unsigned stride) {
    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    // DihedralTS reads positions only — strip everything heavy.
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true;
    opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::DihedralTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(stride);
    return config;
}

}  // namespace


// ── Frame 0 smoke ────────────────────────────────────────────────────

TEST(DihedralTimeSeries, Frame0Semantics) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildConfig(99999);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    EXPECT_EQ(traj.FrameCount(), 1u);

    const auto& tr = tp.Result<nmr::DihedralTimeSeriesTrajectoryResult>();
    EXPECT_EQ(tr.NumFrames(), 1u);
}


// ── Finalize idempotency (data-flow short-circuit) ───────────────────

TEST(DihedralTimeSeries, FinalizeIdempotency) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildConfig(99999);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    auto& tr = tp.Result<nmr::DihedralTimeSeriesTrajectoryResult>();
    const std::size_t T = tr.NumFrames();
    tr.Finalize(tp, traj);
    EXPECT_EQ(tr.NumFrames(), T);
}


// ── H5 round-trip: schema + attrs + dataset shapes ───────────────────

TEST(DihedralTimeSeries, H5RoundTrip) {
    LoadCalculatorConfig();
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

    const auto& tr = tp.Result<nmr::DihedralTimeSeriesTrajectoryResult>();
    const std::size_t R = tp.ProteinRef().ResidueCount();
    const std::size_t T = tr.NumFrames();

    const std::string h5_path = (fs::temp_directory_path() /
        ("dihedral_ts_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/dihedral_time_series"));
    auto grp = reopen.getGroup("/trajectory/dihedral_time_series");

    // Convention attrs.
    std::string units, periodicity, convention;
    grp.getAttribute("angle_units").read(units);
    grp.getAttribute("periodicity").read(periodicity);
    grp.getAttribute("angle_convention").read(convention);
    EXPECT_EQ(units, "radians");
    EXPECT_EQ(periodicity, "2pi");
    EXPECT_NE(convention.find("IUPAC"), std::string::npos);

    // Per-frame datasets.
    for (const std::string& name : {"phi", "psi", "omega", "omega_deviation"}) {
        ASSERT_TRUE(grp.exist(name)) << name;
        const auto dims = grp.getDataSet(name).getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 2u);
        EXPECT_EQ(dims[0], R);
        EXPECT_EQ(dims[1], T);
    }
    {
        ASSERT_TRUE(grp.exist("chi"));
        const auto dims = grp.getDataSet("chi").getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 3u);
        EXPECT_EQ(dims[0], R);
        EXPECT_EQ(dims[1], T);
        EXPECT_EQ(dims[2], 4u);
    }
    {
        ASSERT_TRUE(grp.exist("rama_region"));
        const auto dims = grp.getDataSet("rama_region").getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 2u);
        EXPECT_EQ(dims[0], R);
        EXPECT_EQ(dims[1], T);
    }

    // Static masks.
    EXPECT_TRUE(grp.exist("chi_exists"));
    EXPECT_TRUE(grp.exist("omega_is_xpro"));
    EXPECT_TRUE(grp.exist("is_glycine"));
    EXPECT_TRUE(grp.exist("is_proline"));
    EXPECT_TRUE(grp.exist("is_pre_proline"));
    EXPECT_TRUE(grp.exist("residue_terminal_state"));
    EXPECT_TRUE(grp.exist("chain_id_per_residue"));
    EXPECT_TRUE(grp.exist("residue_index_per_atom"));
    EXPECT_TRUE(grp.exist("frame_indices"));
    EXPECT_TRUE(grp.exist("frame_times"));
    EXPECT_TRUE(grp.exist("source_attached_per_frame"));

    fs::remove(h5_path);
}


// ── Integration on 1P9J: real data, real boundaries ──────────────────

TEST(DihedralTimeSeries, Integration1P9J) {
    LoadCalculatorConfig();
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

    const auto& tr = tp.Result<nmr::DihedralTimeSeriesTrajectoryResult>();
    const std::size_t R = tp.ProteinRef().ResidueCount();
    const std::size_t T = tr.NumFrames();
    EXPECT_GE(T, 2u);

    const std::string h5_path = (fs::temp_directory_path() /
        ("dihedral_ts_int_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/dihedral_time_series");

    // Boundary-NaN discipline: N-terminus residue must have NaN phi.
    // C-terminus must have NaN psi and omega.
    std::vector<std::vector<double>> phi_data, psi_data, omega_data;
    grp.getDataSet("phi").read(phi_data);
    grp.getDataSet("psi").read(psi_data);
    grp.getDataSet("omega").read(omega_data);
    ASSERT_EQ(phi_data.size(), R);
    ASSERT_EQ(phi_data[0].size(), T);

    std::vector<std::uint8_t> terminal_state;
    grp.getDataSet("residue_terminal_state").read(terminal_state);
    ASSERT_EQ(terminal_state.size(), R);

    std::size_t n_term_seen = 0, c_term_seen = 0;
    for (std::size_t ri = 0; ri < R; ++ri) {
        const std::uint8_t ts = terminal_state[ri];
        // 1 = NTerminus, 2 = CTerminus, 3 = NAndCTerminus
        if (ts == 1u || ts == 3u) {
            EXPECT_TRUE(std::isnan(phi_data[ri][0]))
                << "N-terminus residue " << ri << " has finite phi";
            ++n_term_seen;
        }
        if (ts == 2u || ts == 3u) {
            EXPECT_TRUE(std::isnan(psi_data[ri][0]))
                << "C-terminus residue " << ri << " has finite psi";
            EXPECT_TRUE(std::isnan(omega_data[ri][0]))
                << "C-terminus residue " << ri << " has finite omega";
            ++c_term_seen;
        }
    }
    std::cout << "DihedralTimeSeries: " << R << " residues, " << T << " frames; "
              << "termini observed: N=" << n_term_seen
              << " C=" << c_term_seen << "\n";

    // chi_exists static: GLY rows should be all-zero, PHE/TYR rows
    // should be 1 in slots 0+1.
    std::vector<std::vector<std::uint8_t>> chi_exists_data;
    grp.getDataSet("chi_exists").read(chi_exists_data);
    ASSERT_EQ(chi_exists_data.size(), R);
    std::size_t gly_count = 0, two_chi_count = 0;
    for (std::size_t ri = 0; ri < R; ++ri) {
        const auto& row = chi_exists_data[ri];
        ASSERT_EQ(row.size(), 4u);
        const nmr::Residue& res = tp.ProteinRef().ResidueAt(ri);
        if (res.type == nmr::AminoAcid::GLY) {
            EXPECT_EQ(row[0] + row[1] + row[2] + row[3], 0)
                << "GLY residue " << ri << " has chi_exists != all-zero";
            ++gly_count;
        }
        if (res.type == nmr::AminoAcid::PHE ||
            res.type == nmr::AminoAcid::TYR ||
            res.type == nmr::AminoAcid::HIS ||
            res.type == nmr::AminoAcid::TRP) {
            EXPECT_EQ(row[0], 1u);
            EXPECT_EQ(row[1], 1u);
            ++two_chi_count;
        }
    }
    std::cout << "  GLY residues (all-zero chi_exists): " << gly_count
              << "; aromatic residues (>=2 chi): " << two_chi_count << "\n";

    // Range sanity on emitted dihedrals.
    std::size_t phi_finite = 0;
    for (std::size_t ri = 0; ri < R; ++ri) {
        for (std::size_t t = 0; t < T; ++t) {
            const double v = phi_data[ri][t];
            if (std::isfinite(v)) {
                ++phi_finite;
                EXPECT_GE(v, -M_PI) << "phi out of (-pi, pi] at (" << ri << "," << t << ")";
                EXPECT_LE(v,  M_PI) << "phi out of (-pi, pi] at (" << ri << "," << t << ")";
            }
        }
    }
    EXPECT_GE(phi_finite, R - 5u)
        << "Expected most rows to have at least one finite phi sample";

    fs::remove(h5_path);
}

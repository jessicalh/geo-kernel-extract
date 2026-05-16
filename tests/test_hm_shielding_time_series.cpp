//
// test_hm_shielding_time_series: discipline + integration tests for
// HmShieldingTimeSeriesTrajectoryResult.
//
// Mirrors test_tripeptide_bb_shielding_time_series.cpp but the source
// (HaighMallionResult) is unconditionally attached in
// PerFrameExtractionSet — no source-attached gate, no DSN gating.
//

// Nanoflann-using headers must come before GROMACS headers to keep
// gromacs vectypes.h #define DIM 3 from colliding with nanoflann's
// DIM template parameter.
#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "BiotSavartResult.h"
#include "HaighMallionResult.h"

#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "DenseBuffer.h"
#include "HmShieldingTimeSeriesTrajectoryResult.h"
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
#include <vector>

namespace fs = std::filesystem;

#ifndef NMR_TEST_DATA_DIR
#error "NMR_TEST_DATA_DIR must be defined"
#endif

namespace {

constexpr const char* kFixtureProtein = "1P9J_5801";

constexpr int kFingerprintResidueNumber = 28;
const std::string kFingerprintChainId = "A";

// HM shielding magnitude band: literature ring-current shielding at
// backbone amide H near aromatic rings runs ±5-10 ppm. 50 ppm is a
// generous bound that catches a runaway calculation without tripping
// on legitimate ring-proximity values.
constexpr double kT0SanityBoundPpm = 50.0;


std::string TrrPathFor(const std::string& tpr_path) {
    return fs::path(tpr_path).replace_extension(".trr").string();
}

std::string ProductionDirFor(const std::string& tpr_path) {
    return fs::path(tpr_path).parent_path().string();
}

bool FixtureAvailable(const nmr::test::AmberTrajectoryFixture& fix) {
    return !fix.tpr_path.empty() &&
           fs::exists(fix.tpr_path) &&
           fs::exists(TrrPathFor(fix.tpr_path)) &&
           fs::exists(fix.edr_path);
}

void LoadCalculatorConfig() {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);
    nmr::CalculatorConfig::Load(
        std::string(NMR_TEST_DATA_DIR) + "/../data/calculator_params.toml");
}


nmr::SphericalTensor SyntheticTensor(size_t atom_i, size_t frame_t) {
    nmr::SphericalTensor s;
    s.T0 = static_cast<double>(atom_i) + static_cast<double>(frame_t) * 100.0;
    for (size_t k = 0; k < 3; ++k) {
        s.T1[k] = static_cast<double>(atom_i)
                  + static_cast<double>(frame_t) * 100.0
                  + static_cast<double>(k + 1) * 1e-2;
    }
    for (size_t k = 0; k < 5; ++k) {
        s.T2[k] = static_cast<double>(atom_i)
                  + static_cast<double>(frame_t) * 100.0
                  + static_cast<double>(k + 1) * 1e-3;
    }
    return s;
}


bool SphericalEqual(const nmr::SphericalTensor& a,
                    const nmr::SphericalTensor& b,
                    double tol) {
    if (std::abs(a.T0 - b.T0) > tol) return false;
    for (size_t k = 0; k < 3; ++k)
        if (std::abs(a.T1[k] - b.T1[k]) > tol) return false;
    for (size_t k = 0; k < 5; ++k)
        if (std::abs(a.T2[k] - b.T2[k]) > tol) return false;
    return true;
}

}  // namespace


// ============================================================================
// UNIT: 4-frame synthetic round-trip.
// ============================================================================

TEST(HmShieldingTimeSeries, SyntheticFourFrames) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();

    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix))
        GTEST_SKIP() << "fleet_amber " << kFixtureProtein
                     << " fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    const size_t Ntp = tp.AtomCount();
    ASSERT_GT(Ntp, 0u);

    auto tr = nmr::HmShieldingTimeSeriesTrajectoryResult::Create(tp);

    nmr::Trajectory traj(TrrPathFor(fix.tpr_path),
                         fix.tpr_path, fix.edr_path);

    constexpr size_t kFrames = 4;
    const auto& protein_ref = tp.ProteinRef();
    std::vector<nmr::Vec3> positions(Ntp, nmr::Vec3::Zero());

    for (size_t t = 0; t < kFrames; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &protein_ref, positions, "synthetic frame");
        for (size_t i = 0; i < Ntp; ++i) {
            conf->MutableAtomAt(i).hm_shielding_contribution =
                SyntheticTensor(i, t);
        }
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    EXPECT_EQ(tr->NumFrames(), kFrames);

    tr->Finalize(tp, traj);

    auto* buf =
        tp.GetDenseBuffer<nmr::SphericalTensor>(std::type_index(typeid(
            nmr::HmShieldingTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);
    EXPECT_EQ(buf->AtomCount(), Ntp);
    EXPECT_EQ(buf->StridePerAtom(), kFrames);

    for (size_t i : {size_t(0), Ntp / 2, Ntp - 1}) {
        for (size_t t = 0; t < kFrames; ++t) {
            const auto expected = SyntheticTensor(i, t);
            const auto& got = buf->At(i, t);
            EXPECT_TRUE(SphericalEqual(got, expected, 1e-12))
                << "buffer mismatch at atom " << i << " frame " << t;
        }
    }

    const std::string h5_path = (fs::temp_directory_path() /
        ("hm_shielding_ts_unit_" +
         std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        tr->WriteH5Group(tp, file);
    }
    ASSERT_TRUE(fs::exists(h5_path));

    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist(
        "/trajectory/hm_shielding_time_series"));
    auto grp = reopen.getGroup(
        "/trajectory/hm_shielding_time_series");
    auto ds = grp.getDataSet("xyz");
    const auto dims = ds.getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 3u);
    EXPECT_EQ(dims[0], Ntp);
    EXPECT_EQ(dims[1], kFrames);
    EXPECT_EQ(dims[2], 9u);

    std::vector<double> flat(Ntp * kFrames * 9);
    ds.read(flat.data());
    const size_t i = Ntp / 2;
    const size_t t = 2;
    const size_t base = (i * kFrames + t) * 9;
    const auto expected = SyntheticTensor(i, t);
    EXPECT_DOUBLE_EQ(flat[base + 0], expected.T0);
    EXPECT_DOUBLE_EQ(flat[base + 1], expected.T1[0]);
    EXPECT_DOUBLE_EQ(flat[base + 2], expected.T1[1]);
    EXPECT_DOUBLE_EQ(flat[base + 3], expected.T1[2]);
    EXPECT_DOUBLE_EQ(flat[base + 4], expected.T2[0]);
    EXPECT_DOUBLE_EQ(flat[base + 5], expected.T2[1]);
    EXPECT_DOUBLE_EQ(flat[base + 6], expected.T2[2]);
    EXPECT_DOUBLE_EQ(flat[base + 7], expected.T2[3]);
    EXPECT_DOUBLE_EQ(flat[base + 8], expected.T2[4]);

    fs::remove(h5_path);
}


// ============================================================================
// DISCIPLINE: Frame-0 semantics (stride ≥ fixture length → only frame 0).
// ============================================================================

TEST(HmShieldingTimeSeries, Frame0Semantics) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix))
        GTEST_SKIP() << "fleet_amber " << kFixtureProtein
                     << " fixture not on disk";

    nmr::RunConfiguration config;
    config.SetName("HmShieldingTimeSeriesFrame0Semantics");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac   = true;
    opts.skip_coulomb = true;
    opts.skip_apbs    = true;
    opts.skip_dssp    = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::HaighMallionResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
            -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::HmShieldingTimeSeriesTrajectoryResult::Create(tp_in);
        });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path),
                         fix.tpr_path, fix.edr_path);
    nmr::Session session;
    const nmr::Status s = traj.Run(tp, config, session);
    ASSERT_EQ(s, nmr::kOk);
    ASSERT_EQ(traj.FrameCount(), 1u);

    ASSERT_TRUE(tp.HasResult<nmr::HmShieldingTimeSeriesTrajectoryResult>());

    auto* buf = tp.GetDenseBuffer<nmr::SphericalTensor>(std::type_index(
        typeid(nmr::HmShieldingTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);
    EXPECT_EQ(buf->AtomCount(), tp.AtomCount());
    EXPECT_EQ(buf->StridePerAtom(), 1u);
}


// ============================================================================
// DISCIPLINE: Finalize idempotency.
// ============================================================================

TEST(HmShieldingTimeSeries, FinalizeIdempotency) {
    LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    config.SetName("HmShieldingTimeSeriesFinalizeIdempotency");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac   = true;
    opts.skip_coulomb = true;
    opts.skip_apbs    = true;
    opts.skip_dssp    = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::HaighMallionResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
            -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::HmShieldingTimeSeriesTrajectoryResult::Create(tp_in);
        });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path),
                         fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    auto* buf_first = tp.GetDenseBuffer<nmr::SphericalTensor>(
        std::type_index(typeid(nmr::HmShieldingTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf_first, nullptr);
    const std::size_t N_first = buf_first->AtomCount();
    const std::size_t T_first = buf_first->StridePerAtom();

    auto& tr = tp.Result<nmr::HmShieldingTimeSeriesTrajectoryResult>();
    tr.Finalize(tp, traj);

    auto* buf_second = tp.GetDenseBuffer<nmr::SphericalTensor>(
        std::type_index(typeid(nmr::HmShieldingTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf_second, nullptr);
    EXPECT_EQ(buf_second->AtomCount(), N_first);
    EXPECT_EQ(buf_second->StridePerAtom(), T_first);
    EXPECT_EQ(tr.NumFrames(), T_first);
}


// ============================================================================
// DISCIPLINE: H5 round-trip via a temp file.
// ============================================================================

TEST(HmShieldingTimeSeries, H5RoundTrip) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    config.SetName("HmShieldingTimeSeriesH5RoundTrip");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac   = true;
    opts.skip_coulomb = true;
    opts.skip_apbs    = true;
    opts.skip_dssp    = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::HaighMallionResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
            -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::HmShieldingTimeSeriesTrajectoryResult::Create(tp_in);
        });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path),
                         fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::HmShieldingTimeSeriesTrajectoryResult>();

    const std::string h5_path = (fs::temp_directory_path() /
        ("hm_shielding_ts_h5_roundtrip_" +
         std::to_string(::getpid()) + ".h5")).string();
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        tr.WriteH5Group(tp, file);
    }
    ASSERT_TRUE(fs::exists(h5_path));

    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/hm_shielding_time_series"));

    auto grp = reopen.getGroup("/trajectory/hm_shielding_time_series");
    auto ds = grp.getDataSet("xyz");
    const auto dims = ds.getSpace().getDimensions();
    ASSERT_EQ(dims.size(), 3u);
    EXPECT_EQ(dims[0], tp.AtomCount());
    EXPECT_EQ(dims[1], 1u);
    EXPECT_EQ(dims[2], 9u);

    std::string parity, normalization, units, layout;
    grp.getAttribute("parity").read(parity);
    grp.getAttribute("normalization").read(normalization);
    grp.getAttribute("units").read(units);
    grp.getAttribute("irrep_layout").read(layout);
    EXPECT_EQ(parity, "0e+1o+2e");
    EXPECT_EQ(normalization, "isometric_real_sph");
    EXPECT_EQ(units, "ppm");
    EXPECT_EQ(layout,
        "T0,T1_m-1,T1_m0,T1_m+1,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2");

    fs::remove(h5_path);
}


// ============================================================================
// INTEGRATION: Runs Trajectory::Run with HaighMallionResult attached;
// verifies the per-atom HM shielding TR captures finite values across
// all atoms / frames and that at least one atom has non-trivial
// shielding (catches all-zero emission).
// ============================================================================

TEST(HmShieldingTimeSeries, Integration1P9J) {
    LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix))
        GTEST_SKIP() << "fleet_amber " << kFixtureProtein
                     << " fixture not on disk";

    nmr::RunConfiguration config;
    config.SetName("HmShieldingTimeSeriesIntegration");
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac   = true;
    opts.skip_coulomb = true;
    opts.skip_apbs    = true;
    opts.skip_dssp    = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.RequireConformationResult(typeid(nmr::HaighMallionResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
            -> std::unique_ptr<nmr::TrajectoryResult> {
            return nmr::HmShieldingTimeSeriesTrajectoryResult::Create(tp_in);
        });
    config.SetStride(300);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path),
                         fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    EXPECT_GE(traj.FrameCount(), 1u);
    EXPECT_LE(traj.FrameCount(), 10u);

    // Anchor: any CA atom in a representative residue — HM may be zero
    // there if no ring is nearby, but finiteness is required.
    const auto& protein = tp.ProteinRef();
    size_t fingerprint_atom = nmr::Residue::NONE;
    for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
        const auto& res = protein.ResidueAt(ri);
        if (res.chain_id != kFingerprintChainId) continue;
        if (res.sequence_number != kFingerprintResidueNumber) continue;
        fingerprint_atom = res.CA;
        break;
    }
    ASSERT_NE(fingerprint_atom, nmr::Residue::NONE);

    auto* buf = tp.GetDenseBuffer<nmr::SphericalTensor>(std::type_index(
        typeid(nmr::HmShieldingTimeSeriesTrajectoryResult)));
    ASSERT_NE(buf, nullptr);
    EXPECT_EQ(buf->AtomCount(), tp.AtomCount());
    EXPECT_EQ(buf->StridePerAtom(), traj.FrameCount());

    const nmr::SphericalTensor& cell = buf->At(fingerprint_atom, 0);
    EXPECT_TRUE(std::isfinite(cell.T0)
             && std::isfinite(cell.T1[0]) && std::isfinite(cell.T1[1])
             && std::isfinite(cell.T1[2])
             && std::isfinite(cell.T2[0]) && std::isfinite(cell.T2[1])
             && std::isfinite(cell.T2[2]) && std::isfinite(cell.T2[3])
             && std::isfinite(cell.T2[4]))
        << "non-finite component at fingerprint atom " << fingerprint_atom;
    EXPECT_LT(std::abs(cell.T0), kT0SanityBoundPpm)
        << "T0 outside HM sanity band: " << cell.T0 << " ppm at atom "
        << fingerprint_atom;

    // Full-buffer 9-component isfinite sweep + population floor. HM is
    // sparse — only atoms within ~15 Å of an aromatic ring contribute.
    // On 1P9J_5801 (56 residues, multiple aromatic residues), expect
    // at least one atom with measurable shielding; a floor of 1 is a
    // catastrophic-all-zero check, not a coverage assertion.
    size_t populated = 0;
    for (size_t i = 0; i < buf->AtomCount(); ++i) {
        for (size_t t = 0; t < buf->StridePerAtom(); ++t) {
            const nmr::SphericalTensor& st = buf->At(i, t);
            EXPECT_TRUE(std::isfinite(st.T0)
                     && std::isfinite(st.T1[0]) && std::isfinite(st.T1[1])
                     && std::isfinite(st.T1[2])
                     && std::isfinite(st.T2[0]) && std::isfinite(st.T2[1])
                     && std::isfinite(st.T2[2]) && std::isfinite(st.T2[3])
                     && std::isfinite(st.T2[4]))
                << "non-finite component at atom " << i << " frame " << t;
        }
        if (std::abs(buf->At(i, 0).T0) > 1e-12) ++populated;
    }
    EXPECT_GT(populated, 0u)
        << "HM shielding populated zero atoms — no aromatic ring contribution "
        << "detected on a fixture with aromatic residues, suggests calculator "
        << "regression";

    std::cout << "HmShieldingTimeSeries CA-anchor diagnostic: "
              << "chain=" << kFingerprintChainId
              << " residue=" << kFingerprintResidueNumber
              << " atom_idx=" << fingerprint_atom
              << " frame 0 T0=" << cell.T0 << " ppm, populated="
              << populated << "/" << buf->AtomCount() << "\n";
}

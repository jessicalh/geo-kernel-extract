//
// test_ring_neighbourhood_trajectory_stats: discipline + integration
// for RingNeighbourhoodTrajectoryStats (TR10). Static (atom, ring)
// snapshot at frame 0 + per-frame geometric residual on those pairs.
// Self-contained TR (positions + GeometryResult + SpatialIndexResult,
// no source ConformationResult dependency).
//

#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "EnrichmentResult.h"
#include "GeometryResult.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "Ring.h"
#include "RingNeighbourhoodTrajectoryStats.h"
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
#include <filesystem>
#include <limits>
#include <memory>
#include <set>
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
    // TR10 reads positions + ring_geometries + spatial index only.
    // Strip everything heavy.
    opts.skip_mopac = true;
    opts.skip_coulomb = true;
    opts.skip_apbs = true;
    opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::RingNeighbourhoodTrajectoryStats::Create(tp_in);
    });
    config.SetStride(stride);
    return config;
}

}  // namespace


// ── Frame 0 smoke ───────────────────────────────────────────────────

TEST(RingNeighbourhoodTrajectoryStats, Frame0Semantics) {
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
    EXPECT_EQ(traj.FrameCount(), 1u);

    const auto& tr = tp.Result<nmr::RingNeighbourhoodTrajectoryStats>();
    EXPECT_EQ(tr.NumFrames(), 1u);
    EXPECT_EQ(tr.NumAtoms(), tp.AtomCount());
    // 1P9J has aromatic residues; per_atom_ring_list should be non-empty
    // for at least some atoms (atoms near PHE/TYR/TRP/HIS rings).
    EXPECT_GT(tr.RPerAtomMax(), 0u)
        << "1P9J should have at least one (atom, ring) pair within cutoff";
}


// ── Finalize idempotency ────────────────────────────────────────────

TEST(RingNeighbourhoodTrajectoryStats, FinalizeIdempotency) {
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

    auto& tr = tp.Result<nmr::RingNeighbourhoodTrajectoryStats>();
    const std::size_t T = tr.NumFrames();
    const std::size_t R = tr.RPerAtomMax();
    tr.Finalize(tp, traj);  // second call should be safe no-op-ish
    EXPECT_EQ(tr.NumFrames(), T);
    EXPECT_EQ(tr.RPerAtomMax(), R);
}


// ── H5 round-trip: schema + attrs + dataset shapes ──────────────────

TEST(RingNeighbourhoodTrajectoryStats, H5RoundTrip) {
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

    const auto& tr = tp.Result<nmr::RingNeighbourhoodTrajectoryStats>();
    const std::size_t N = tp.AtomCount();
    const std::size_t T = tr.NumFrames();
    const std::size_t R = tr.RPerAtomMax();
    if (R == 0) GTEST_SKIP() << "no aromatic-ring pairs in protein";

    const std::string h5_path = (fs::temp_directory_path() /
        ("rnts_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/ring_neighbourhood_trajectory_stats"));
    auto grp = reopen.getGroup("/trajectory/ring_neighbourhood_trajectory_stats");

    // Convention attrs.
    std::string channel_layout, units, policy;
    grp.getAttribute("channel_layout").read(channel_layout);
    grp.getAttribute("units").read(units);
    grp.getAttribute("source_attached_policy").read(policy);
    EXPECT_EQ(channel_layout, "distance,rho,z,in_plane_angle");
    EXPECT_NE(units.find("Angstrom"), std::string::npos);
    EXPECT_NE(units.find("radians"), std::string::npos);
    EXPECT_NE(policy.find("always_attached"), std::string::npos);

    // Cutoff attr matches CalculatorConfig
    double cutoff_attr = 0.0;
    grp.getAttribute("ring_current_spatial_cutoff_A").read(cutoff_attr);
    EXPECT_DOUBLE_EQ(cutoff_attr,
        nmr::CalculatorConfig::Get("ring_current_spatial_cutoff"));

    // Dataset shapes
    {
        ASSERT_TRUE(grp.exist("geometry"));
        const auto dims = grp.getDataSet("geometry").getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 4u);
        EXPECT_EQ(dims[0], N);
        EXPECT_EQ(dims[1], T);
        EXPECT_EQ(dims[2], R);
        EXPECT_EQ(dims[3], nmr::RingNeighbourhoodTrajectoryStats::kChannelCount);
    }
    {
        ASSERT_TRUE(grp.exist("ring_membership_per_atom"));
        const auto dims = grp.getDataSet("ring_membership_per_atom")
            .getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 2u);
        EXPECT_EQ(dims[0], N);
        EXPECT_EQ(dims[1], R);
    }
    EXPECT_TRUE(grp.exist("frame_indices"));
    EXPECT_TRUE(grp.exist("frame_times"));
    EXPECT_TRUE(grp.exist("source_attached_per_frame"));

    // ring_membership_per_atom invariants:
    //   - all values in [-1, n_aromatic-1]
    //   - per-row: live indices come first (sorted asc), then -1 padding
    //   - sentinel slot count consistent with R_max
    std::vector<std::vector<std::int32_t>> mem2d;
    grp.getDataSet("ring_membership_per_atom").read(mem2d);
    ASSERT_EQ(mem2d.size(), N);
    std::vector<std::int32_t> mem(N * R);
    for (std::size_t i = 0; i < N; ++i) {
        ASSERT_EQ(mem2d[i].size(), R);
        for (std::size_t r = 0; r < R; ++r) mem[i * R + r] = mem2d[i][r];
    }
    const std::size_t n_aromatic = tp.ProteinRef().RingCount();
    for (std::size_t i = 0; i < N; ++i) {
        bool sentinel_seen = false;
        std::int32_t prev_idx = -1;
        for (std::size_t r = 0; r < R; ++r) {
            const std::int32_t v = mem[i * R + r];
            if (v == -1) {
                sentinel_seen = true;
                continue;
            }
            // No live index after a sentinel (live slots are first)
            EXPECT_FALSE(sentinel_seen)
                << "atom " << i << " slot " << r << " is live (" << v
                << ") after sentinel";
            EXPECT_GE(v, 0);
            EXPECT_LT(static_cast<std::size_t>(v), n_aromatic)
                << "atom " << i << " slot " << r
                << " references non-aromatic ring " << v;
            // Sorted ascending: ring indices monotonically increasing
            if (prev_idx >= 0) {
                EXPECT_GT(v, prev_idx)
                    << "atom " << i << " slot " << r
                    << " membership not strictly ascending";
            }
            prev_idx = v;
        }
    }
}


// ── Geometric content sanity: finite + range checks ─────────────────

TEST(RingNeighbourhoodTrajectoryStats, GeometryRanges1P9J) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    // Stride 300 over 1P9J's 1501-frame fixture = ~5 frames.
    auto config = BuildConfig(300);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::RingNeighbourhoodTrajectoryStats>();
    const std::size_t N = tp.AtomCount();
    const std::size_t T = tr.NumFrames();
    const std::size_t R = tr.RPerAtomMax();
    if (R == 0) GTEST_SKIP() << "no aromatic-ring pairs in protein";

    const std::string h5_path = (fs::temp_directory_path() /
        ("rnts_geom_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/ring_neighbourhood_trajectory_stats");
    std::vector<std::vector<std::vector<std::vector<double>>>> geom4d;
    grp.getDataSet("geometry").read(geom4d);
    ASSERT_EQ(geom4d.size(), N);
    std::vector<double> geom(N * T * R * 4u);
    for (std::size_t i = 0; i < N; ++i) {
        ASSERT_EQ(geom4d[i].size(), T);
        for (std::size_t t = 0; t < T; ++t) {
            ASSERT_EQ(geom4d[i][t].size(), R);
            for (std::size_t r = 0; r < R; ++r) {
                ASSERT_EQ(geom4d[i][t][r].size(), 4u);
                for (std::size_t ch = 0; ch < 4; ++ch)
                    geom[((i * T + t) * R + r) * 4 + ch] = geom4d[i][t][r][ch];
            }
        }
    }
    std::vector<std::vector<std::int32_t>> mem2d;
    grp.getDataSet("ring_membership_per_atom").read(mem2d);
    std::vector<std::int32_t> mem(N * R);
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t r = 0; r < R; ++r) mem[i * R + r] = mem2d[i][r];

    // Per-channel sanity, restricted to live (atom, ring) slots.
    std::size_t n_live = 0;
    std::size_t n_overage = 0;  // distance > cutoff after drift -- counted, not asserted
    const double cutoff = nmr::CalculatorConfig::Get("ring_current_spatial_cutoff");
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t t = 0; t < T; ++t) {
            for (std::size_t r = 0; r < R; ++r) {
                if (mem[i * R + r] == -1) {
                    // Padded slot: all 4 channels must be NaN.
                    for (std::size_t ch = 0; ch < 4; ++ch) {
                        const double v = geom[((i * T + t) * R + r) * 4 + ch];
                        EXPECT_TRUE(std::isnan(v))
                            << "atom " << i << " frame " << t << " slot " << r
                            << " channel " << ch
                            << " in padded slot must be NaN, got " << v;
                    }
                    continue;
                }
                ++n_live;
                const std::size_t base = ((i * T + t) * R + r) * 4;
                const double dist = geom[base + 0];
                const double rho  = geom[base + 1];
                const double z    = geom[base + 2];
                const double phi  = geom[base + 3];

                ASSERT_TRUE(std::isfinite(dist))
                    << "non-finite distance at (" << i << "," << t << "," << r << ")";
                ASSERT_TRUE(std::isfinite(rho))
                    << "non-finite rho";
                ASSERT_TRUE(std::isfinite(z))
                    << "non-finite z";
                // phi may be NaN at rho ≈ 0 by design (singular case).

                EXPECT_GE(dist, 0.0);
                EXPECT_GE(rho, 0.0);
                // distance^2 ≈ rho^2 + z^2 (within roundoff)
                const double check = std::sqrt(rho*rho + z*z);
                EXPECT_NEAR(dist, check, 1e-9)
                    << "rho/z decomposition inconsistent at (" << i << ","
                    << t << "," << r << ")";

                if (std::isfinite(phi)) {
                    EXPECT_GE(phi, 0.0);
                    EXPECT_LT(phi, 2.0 * M_PI);
                }

                if (dist > cutoff) ++n_overage;
            }
        }
    }
    EXPECT_GT(n_live, 0u) << "no live (atom, ring) slots populated";
    nmr::OperationLog::Info(
        "RingNeighbourhoodTrajectoryStatsTest::GeometryRanges1P9J",
        "n_live=" + std::to_string(n_live) +
        " n_overage_past_cutoff=" + std::to_string(n_overage) +
        " (drift past 15A acceptable; consumer filters on distance)");
}


// ── Literature-anchored probe (PHE ring membership population) ──────
//
// 1P9J is an SH3 domain with multiple PHE rings (canonical ring-current
// system). For the FIRST PHE ring found, confirm the static membership
// captures the expected pattern: the ring's own vertex atoms always sit
// within cutoff, plus many adjacent-residue heavy atoms (15A includes
// neighbouring residues' backbone + nearby sidechains).
//
// Robust to PDB-vs-BMRB residue numbering differences -- doesn't pin to
// PHE33 specifically.

TEST(RingNeighbourhoodTrajectoryStats, PheRingMembershipProbe) {
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

    const auto& tr = tp.Result<nmr::RingNeighbourhoodTrajectoryStats>();
    if (tr.RPerAtomMax() == 0) GTEST_SKIP() << "no aromatic rings in protein";

    // Find the first PHE ring in the protein.
    std::size_t phe_ring_idx = SIZE_MAX;
    int phe_seq = -1;
    for (std::size_t ri = 0; ri < tp.ProteinRef().RingCount(); ++ri) {
        const auto& ring = tp.ProteinRef().RingAt(ri);
        const auto& res = tp.ProteinRef().ResidueAt(ring.parent_residue_index);
        if (res.type == nmr::AminoAcid::PHE) {
            phe_ring_idx = ri;
            phe_seq = res.sequence_number;
            break;
        }
    }
    if (phe_ring_idx == SIZE_MAX) {
        GTEST_SKIP() << "no PHE ring in this protein";
    }

    // Vertex-atom invariant: each of the ring's own vertices MUST appear
    // with this ring index in its membership (distance to center is the
    // ring radius < cutoff, by construction).
    const auto& phe_ring = tp.ProteinRef().RingAt(phe_ring_idx);
    for (std::size_t v : phe_ring.atom_indices) {
        const auto& rings = tr.RingListForAtom(v);
        EXPECT_NE(std::find(rings.begin(), rings.end(), phe_ring_idx),
                   rings.end())
            << "PHE ring " << phe_ring_idx << " (resnum " << phe_seq
            << ") missing from its own vertex atom " << v
            << "'s membership snapshot";
    }

    // Population-level probe: the PHE ring's index should appear for
    // many atoms within 15A. SH3-domain proteins are compact (~80
    // residues); a central PHE ring's 15A neighbourhood typically
    // includes 50-200 atoms.
    std::size_t n_atoms_with_phe = 0;
    for (std::size_t ai = 0; ai < tp.AtomCount(); ++ai) {
        const auto& rings = tr.RingListForAtom(ai);
        if (std::find(rings.begin(), rings.end(), phe_ring_idx)
            != rings.end()) {
            ++n_atoms_with_phe;
        }
    }
    EXPECT_GT(n_atoms_with_phe, 10u)
        << "PHE ring (resnum " << phe_seq << ") should appear in many "
        << "atoms' membership (15A cutoff is wide for SH3 packing)";
    nmr::OperationLog::Info(
        "RingNeighbourhoodTrajectoryStatsTest::PheRingMembershipProbe",
        "phe_ring_idx=" + std::to_string(phe_ring_idx) +
        " phe_seq=" + std::to_string(phe_seq) +
        " n_atoms_with_phe_in_static_snapshot=" +
        std::to_string(n_atoms_with_phe));
}

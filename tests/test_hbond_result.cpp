#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include <array>
#include <cstdint>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <unistd.h>

#include <highfive/H5DataSet.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include "HBondResult.h"
#include "HBondCountWelfordTrajectoryResult.h"
#include "McConnellResult.h"
#include "CoulombResult.h"
#include "RingSusceptibilityResult.h"
#include "DsspResult.h"
#include "ChargeAssignmentResult.h"
#include "ChargeSource.h"
#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "OrcaRunLoader.h"
#include "OrcaShieldingResult.h"
#include "MutationDeltaResult.h"
#include "OperationLog.h"
#include "PhysicalConstants.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"

#include <filesystem>
namespace fs = std::filesystem;
using namespace nmr;

namespace {

struct SyntheticHBondFixture {
    std::unique_ptr<Protein> protein;
    size_t donor_n = SIZE_MAX;
    size_t donor_h = SIZE_MAX;
    size_t acceptor_o = SIZE_MAX;
    size_t local_target = SIZE_MAX;
    size_t remote_target = SIZE_MAX;
};

SyntheticHBondFixture BuildNonCollinearExplicitHydrogenFixture() {
    SyntheticHBondFixture f;
    f.protein = std::make_unique<Protein>();

    // Residue-vector order deliberately crosses a chain boundary at 1 -> 2.
    // The donor is the final residue of chain A and the acceptor is the first
    // residue of chain B, so their storage indices differ by one even though
    // there is no peptide path between them.  Residue 0 is a real peptide
    // predecessor of the donor; residue 3 is a disconnected control.
    const std::array<const char*, 4> chains{{"A", "A", "B", "C"}};
    for (int ri = 0; ri < 4; ++ri) {
        Residue r;
        r.type = AminoAcid::ALA;
        r.sequence_number = ri + 1;
        r.chain_id = chains[static_cast<size_t>(ri)];
        f.protein->AddResidue(std::move(r));
    }

    std::vector<Vec3> positions;
    auto add_atom = [&](size_t residue, const char* name, Element element,
                        const Vec3& pos) {
        auto atom = Atom::Create(element);
        atom->pdb_atom_name = name;
        atom->residue_index = residue;
        const size_t ai = f.protein->AddAtom(std::move(atom));
        f.protein->MutableResidueAt(residue).atom_indices.push_back(ai);
        positions.push_back(pos);
        return ai;
    };

    // C(0)-N(1) is a real PeptideCN edge.  The H-bond itself spans the
    // A/B chain boundary (residue 1 -> residue 2), and is deliberately
    // non-collinear in N-H...O so the explicit-H source convention remains
    // independently pinned by the same fixture.
    f.local_target = add_atom(
        0, "C", Element::C, Vec3(-1.33, 0.0, 0.0));
    f.donor_n = add_atom(1, "N", Element::N, Vec3(0.0, 0.0, 0.0));
    f.donor_h = add_atom(1, "H", Element::H, Vec3(1.0, 0.0, 0.0));
    f.acceptor_o = add_atom(2, "O", Element::O, Vec3(1.0, 2.0, 0.0));
    f.remote_target = add_atom(
        3, "C", Element::C, Vec3(2.0, 2.0, 0.0));

    f.protein->FinalizeConstruction(positions);
    f.protein->AddConformation(std::move(positions), "explicit-H forcing fixture");
    return f;
}

template <typename T>
std::vector<T> ReadNpyPayload(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    EXPECT_TRUE(in.is_open()) << path;
    if (!in.is_open()) return {};

    char magic[6] = {};
    in.read(magic, sizeof(magic));
    EXPECT_EQ(std::string(magic, sizeof(magic)),
              std::string("\x93NUMPY", sizeof(magic)));
    char version[2] = {};
    in.read(version, sizeof(version));
    EXPECT_EQ(version[0], 1);
    EXPECT_EQ(version[1], 0);
    std::uint16_t header_length = 0;
    in.read(reinterpret_cast<char*>(&header_length), sizeof(header_length));
    in.seekg(header_length, std::ios::cur);

    const std::streampos payload_begin = in.tellg();
    in.seekg(0, std::ios::end);
    const std::streamoff payload_bytes = in.tellg() - payload_begin;
    EXPECT_EQ(payload_bytes % static_cast<std::streamoff>(sizeof(T)), 0);
    in.seekg(payload_begin);

    std::vector<T> values(
        static_cast<size_t>(payload_bytes / sizeof(T)));
    if (!values.empty()) {
        in.read(reinterpret_cast<char*>(values.data()), payload_bytes);
    }
    return values;
}

}  // namespace


TEST(HBondGeometryKernel, UsesExplicitHydrogenAndAppliesTargetSequenceFilter) {
    auto f = BuildNonCollinearExplicitHydrogenFixture();
    auto& conf = f.protein->Conformation();

    ASSERT_TRUE(conf.AttachResult(GeometryResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(SpatialIndexResult::Compute(conf)));

    std::vector<DsspResidue> residues(f.protein->ResidueCount());
    residues[1].observed = true;
    residues[1].acceptors[0].residue_index = 2;
    ASSERT_TRUE(conf.AttachResult(DsspResult::CreateForTesting(std::move(residues))));

    // The graph, not vector adjacency, owns all three answers.
    EXPECT_EQ(hbond_result_detail::BackboneResidueSeparation(
                  *f.protein, 1, 2),
              -1) << "adjacent residue slots on different chains are disconnected";
    EXPECT_EQ(hbond_result_detail::BackboneResidueSeparation(
                  *f.protein, 0, 1),
              1) << "the explicit C(0)-N(1) peptide bond is one step";
    EXPECT_EQ(hbond_result_detail::BackboneResidueSeparation(
                  *f.protein, 3, 1),
              -1);

    auto hbond = HBondResult::Compute(conf);
    ASSERT_NE(hbond, nullptr);
    EXPECT_EQ(hbond->HBondCount(), 1u)
        << "cross-chain H-bond must not be discarded because its residue "
           "storage indices are adjacent";
    ASSERT_TRUE(conf.AttachResult(std::move(hbond)));

    // Blessed values for donor H=(1,0,0), acceptor O=(1,2,0), and remote
    // target=(2,2,0).  These numbers are pinned independently, rather than
    // recomputing the production formula in the test.
    const auto& remote = conf.AtomAt(f.remote_target);
    EXPECT_EQ(remote.hbond_count_within_3_5A, 1)
        << "the disconnected-chain H-bond is inside the 3.5 A count shell";
    EXPECT_NEAR(remote.hbond_nearest_dist, 2.236067977499790, 1e-12);
    EXPECT_NEAR(remote.hbond_inv_d3, 0.089442719099992, 1e-12);
    EXPECT_NEAR(remote.hbond_mcconnell_scalar, 0.125219806739988, 1e-12);
    EXPECT_NEAR(remote.hbond_nearest_dir.x(), 0.447213595499958, 1e-12);
    EXPECT_NEAR(remote.hbond_nearest_dir.y(), 0.894427190999916, 1e-12);
    EXPECT_NEAR(remote.hbond_nearest_dir.z(), 0.0, 1e-12);

    // Residue 0 is the donor's true peptide predecessor and therefore lies
    // within the configured two-residue endpoint exclusion.
    // It is geometrically valid, so a zero here specifically freezes the
    // production SequentialExclusionFilter wiring (C18).
    const auto& local = conf.AtomAt(f.local_target);
    EXPECT_EQ(local.hbond_count_within_3_5A, 0);
    EXPECT_DOUBLE_EQ(local.hbond_mcconnell_scalar, 0.0);
    EXPECT_DOUBLE_EQ(local.hbond_inv_d3, 0.0);

    // Directly pin the production helper's explicit-H convention too.  It
    // has external linkage in HBondResult's per-file named namespace.
    const auto kernel = hbond_result_detail::ComputeKernel(
        Vec3(2.0, 2.0, 0.0), Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0));
    EXPECT_NEAR(kernel.distance, 2.236067977499790, 1e-12);
    EXPECT_NEAR(kernel.f, 0.125219806739988, 1e-12);

    // Pin the changed producer arrays, not merely their in-memory source.
    // The expected count/flags come from the independently constructed
    // chain graph and DSSP pair above; WriteFeatures is pure read-back.
    const fs::path output_dir = fs::temp_directory_path() /
        ("hbond_multichain_forcing_" + std::to_string(::getpid()));
    ASSERT_TRUE(fs::create_directories(output_dir));
    ASSERT_EQ(conf.Result<HBondResult>().WriteFeatures(
                  conf, output_dir.string()),
              6);
    const auto scalars = ReadNpyPayload<double>(
        output_dir / "hbond_scalars.npy");
    const auto flags = ReadNpyPayload<std::int8_t>(
        output_dir / "hbond_flags.npy");
    ASSERT_EQ(scalars.size(), conf.AtomCount() * 4);
    ASSERT_EQ(flags.size(), conf.AtomCount() * 3);
    EXPECT_DOUBLE_EQ(scalars[f.remote_target * 4 + 2], 1.0);
    EXPECT_DOUBLE_EQ(scalars[f.local_target * 4 + 2], 0.0);
    EXPECT_EQ(flags[f.remote_target * 3 + 0], 1);
    EXPECT_EQ(flags[f.donor_n * 3 + 1], 1);
    EXPECT_EQ(flags[f.acceptor_o * 3 + 2], 1);
    EXPECT_EQ(flags[f.local_target * 3 + 0], 0);
    EXPECT_EQ(flags[f.local_target * 3 + 1], 0);
    EXPECT_EQ(flags[f.local_target * 3 + 2], 0);
    // Avoid std::filesystem::remove_all: libtorch exports a broken copy of
    // that symbol in this test executable.  This is the established cleanup
    // pattern used by the other producer-emission tests.
    std::error_code ec;
    fs::remove(output_dir / "hbond_scalars.npy", ec);
    fs::remove(output_dir / "hbond_flags.npy", ec);
    fs::remove(output_dir / "hbond_nearest_dir.npy", ec);
    fs::remove(output_dir, ec);
}


TEST(HBondGeometryKernel,
     CrossChainProductionCountFeedsProductionWelfordConsumer) {
    auto f = BuildNonCollinearExplicitHydrogenFixture();
    ProteinConformation* conf = &f.protein->Conformation();

    // A second, independently computed frame moves only the disconnected
    // target outside the count shell.  The cross-chain DSSP pair remains
    // accepted, giving the Welford consumer the exact sequence [1, 0].
    std::vector<Vec3> far_positions = conf->Positions();
    far_positions[f.remote_target] = Vec3(20.0, 20.0, 0.0);
    auto far_conf = std::make_unique<ProteinConformation>(
        f.protein.get(), std::move(far_positions),
        "cross-chain target outside count shell");

    ASSERT_TRUE(conf->AttachResult(GeometryResult::Compute(*conf)));
    ASSERT_TRUE(conf->AttachResult(SpatialIndexResult::Compute(*conf)));
    std::vector<DsspResidue> residues(f.protein->ResidueCount());
    residues[1].observed = true;
    residues[1].acceptors[0].residue_index = 2;
    ASSERT_TRUE(conf->AttachResult(
        DsspResult::CreateForTesting(std::move(residues))));
    ASSERT_TRUE(conf->AttachResult(HBondResult::Compute(*conf)));
    ASSERT_EQ(conf->AtomAt(f.remote_target).hbond_count_within_3_5A, 1);

    ASSERT_TRUE(far_conf->AttachResult(GeometryResult::Compute(*far_conf)));
    ASSERT_TRUE(far_conf->AttachResult(
        SpatialIndexResult::Compute(*far_conf)));
    std::vector<DsspResidue> far_residues(f.protein->ResidueCount());
    far_residues[1].observed = true;
    far_residues[1].acceptors[0].residue_index = 2;
    ASSERT_TRUE(far_conf->AttachResult(
        DsspResult::CreateForTesting(std::move(far_residues))));
    ASSERT_TRUE(far_conf->AttachResult(HBondResult::Compute(*far_conf)));
    ASSERT_EQ(far_conf->AtomAt(f.remote_target).hbond_count_within_3_5A, 0);

    // The synthetic Protein is already finalized and owns the production
    // HBondResult. Seat that same object in a trajectory buffer without an
    // external TPR/TRR dependency, then dispatch the real Welford Compute.
    auto tp = TrajectoryProtein::CreateForTesting(std::move(f.protein));
    ASSERT_NE(tp, nullptr);
    ASSERT_TRUE(tp->AttachResult(
        HBondCountWelfordTrajectoryResult::Create(*tp)));
    Trajectory trajectory(fs::path{}, fs::path{}, fs::path{});
    tp->DispatchCompute(*conf, trajectory, 17, 4.0);
    tp->DispatchCompute(*far_conf, trajectory, 18, 5.0);
    tp->FinalizeAllResults(trajectory);

    const auto& remote =
        tp->AtomAt(f.remote_target).hbond_count_welford;
    EXPECT_EQ(remote.n_frames, 2u);
    EXPECT_DOUBLE_EQ(remote.count.mean, 0.5);
    EXPECT_DOUBLE_EQ(remote.count.m2, 0.5);
    EXPECT_DOUBLE_EQ(remote.count.std, std::sqrt(0.5));
    EXPECT_DOUBLE_EQ(remote.count.min, 0.0);
    EXPECT_DOUBLE_EQ(remote.count.max, 1.0);
    EXPECT_EQ(remote.count.min_frame, 18u);
    EXPECT_EQ(remote.count.max_frame, 17u);
    EXPECT_DOUBLE_EQ(remote.occupancy_fraction.mean, 0.5);
    EXPECT_DOUBLE_EQ(remote.count_delta.mean, -1.0);
    EXPECT_DOUBLE_EQ(remote.count_abs_delta.mean, 1.0);
    EXPECT_DOUBLE_EQ(remote.count_delta_squared.mean, 1.0);
    EXPECT_DOUBLE_EQ(remote.count_dxdt.mean, -1.0);
    EXPECT_DOUBLE_EQ(remote.count_rms_delta, 1.0);
    EXPECT_EQ(remote.delta_n, 1u);
    EXPECT_EQ(remote.dxdt_n, 1u);

    const auto& local =
        tp->AtomAt(f.local_target).hbond_count_welford;
    EXPECT_EQ(local.n_frames, 2u);
    EXPECT_DOUBLE_EQ(local.count.mean, 0.0);
    EXPECT_DOUBLE_EQ(local.count.std, 0.0);
    EXPECT_DOUBLE_EQ(local.occupancy_fraction.mean, 0.0);
    EXPECT_DOUBLE_EQ(local.count_delta.mean, 0.0);

    // Pin the emitted consumer contract, including canonical datasets and
    // the legacy aliases that downstream SDKs still expose.
    const fs::path h5_path = fs::temp_directory_path() /
        ("hbond_multichain_welford_" + std::to_string(::getpid()) + ".h5");
    {
        HighFive::File file(h5_path.string(), HighFive::File::Truncate);
        tp->WriteH5(file);
    }
    {
        HighFive::File file(h5_path.string(), HighFive::File::ReadOnly);
        ASSERT_TRUE(file.exist("/trajectory/hbond_count_welford"));
        auto group = file.getGroup("/trajectory/hbond_count_welford");

        bool finalized = false;
        size_t n_frames = 0;
        group.getAttribute("finalized").read(finalized);
        group.getAttribute("n_frames").read(n_frames);
        EXPECT_TRUE(finalized);
        EXPECT_EQ(n_frames, 2u);

        auto read_double = [&](const char* name) {
            std::vector<double> values;
            group.getDataSet(name).read(values);
            EXPECT_EQ(values.size(), tp->AtomCount()) << name;
            return values;
        };
        auto read_size = [&](const char* name) {
            std::vector<size_t> values;
            group.getDataSet(name).read(values);
            EXPECT_EQ(values.size(), tp->AtomCount()) << name;
            return values;
        };

        const auto count_mean = read_double("count_mean");
        const auto count_m2 = read_double("count_m2");
        const auto count_std = read_double("count_std");
        const auto count_min = read_double("count_min");
        const auto count_max = read_double("count_max");
        const auto count_min_frame = read_size("count_min_frame");
        const auto count_max_frame = read_size("count_max_frame");
        const auto occupancy_mean =
            read_double("occupancy_fraction_mean");
        const auto delta_mean = read_double("count_delta_mean");
        const auto abs_delta_mean =
            read_double("count_abs_delta_mean");
        const auto delta_squared_mean =
            read_double("count_delta_squared_mean");
        const auto dxdt_mean = read_double("count_dxdt_mean");
        const auto rms_delta = read_double("rms_delta");
        const auto n_frames_per_atom = read_size("n_frames_per_atom");
        const auto delta_n_per_atom = read_size("delta_n_per_atom");
        const auto dxdt_n_per_atom = read_size("dxdt_n_per_atom");
        const auto legacy_mean = read_double("mean");
        const auto legacy_std = read_double("std");
        const auto legacy_min = read_double("min");
        const auto legacy_max = read_double("max");
        const auto legacy_delta_mean = read_double("delta_mean");

        const size_t ri = f.remote_target;
        EXPECT_DOUBLE_EQ(count_mean[ri], 0.5);
        EXPECT_DOUBLE_EQ(count_m2[ri], 0.5);
        EXPECT_DOUBLE_EQ(count_std[ri], std::sqrt(0.5));
        EXPECT_DOUBLE_EQ(count_min[ri], 0.0);
        EXPECT_DOUBLE_EQ(count_max[ri], 1.0);
        EXPECT_EQ(count_min_frame[ri], 18u);
        EXPECT_EQ(count_max_frame[ri], 17u);
        EXPECT_DOUBLE_EQ(occupancy_mean[ri], 0.5);
        EXPECT_DOUBLE_EQ(delta_mean[ri], -1.0);
        EXPECT_DOUBLE_EQ(abs_delta_mean[ri], 1.0);
        EXPECT_DOUBLE_EQ(delta_squared_mean[ri], 1.0);
        EXPECT_DOUBLE_EQ(dxdt_mean[ri], -1.0);
        EXPECT_DOUBLE_EQ(rms_delta[ri], 1.0);
        EXPECT_EQ(n_frames_per_atom[ri], 2u);
        EXPECT_EQ(delta_n_per_atom[ri], 1u);
        EXPECT_EQ(dxdt_n_per_atom[ri], 1u);
        EXPECT_DOUBLE_EQ(legacy_mean[ri], count_mean[ri]);
        EXPECT_DOUBLE_EQ(legacy_std[ri], count_std[ri]);
        EXPECT_DOUBLE_EQ(legacy_min[ri], count_min[ri]);
        EXPECT_DOUBLE_EQ(legacy_max[ri], count_max[ri]);
        EXPECT_DOUBLE_EQ(legacy_delta_mean[ri], delta_mean[ri]);

        const size_t li = f.local_target;
        EXPECT_DOUBLE_EQ(count_mean[li], 0.0);
        EXPECT_DOUBLE_EQ(occupancy_mean[li], 0.0);
        EXPECT_DOUBLE_EQ(delta_mean[li], 0.0);
    }
    std::error_code ec;
    fs::remove(h5_path, ec);
}



// ============================================================================
// Single protein test on ORCA test protein (A0A7C5FAR6)
// ============================================================================

TEST(HBondOrcaTest, RunOnProtonatedProtein) {
    OrcaRunFiles files;
    files.pdb_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.pdb";
    files.xyz_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.xyz";
    files.prmtop_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.prmtop";

    if (!fs::exists(files.xyz_path) || !fs::exists(files.prmtop_path))
        GTEST_SKIP() << "ORCA test data not found";

    auto load = BuildFromOrca(files);
    ASSERT_TRUE(load.Ok()) << load.error;

    auto& conf = load.protein->Conformation();
    conf.AttachResult(GeometryResult::Compute(conf));
    conf.AttachResult(SpatialIndexResult::Compute(conf));

    auto dssp = DsspResult::Compute(conf);
    ASSERT_NE(dssp, nullptr) << "DSSP must succeed";
    conf.AttachResult(std::move(dssp));

    auto hbond = HBondResult::Compute(conf);
    ASSERT_NE(hbond, nullptr);
    conf.AttachResult(std::move(hbond));

    int checked = 0;
    int has_hbond = 0;
    double max_scalar = 0;
    int donors = 0, acceptors = 0;

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& ca = conf.AtomAt(ai);

        if (ca.hbond_nearest_dist > 0 && ca.hbond_nearest_dist < NO_DATA_SENTINEL) {
            has_hbond++;
            EXPECT_GT(ca.hbond_inv_d3, 0.0);
            EXPECT_NEAR(ca.hbond_nearest_dir.norm(), 1.0, 1e-8);
            checked++;
        }

        max_scalar = std::max(max_scalar, std::abs(ca.hbond_mcconnell_scalar));
        if (ca.hbond_is_donor) donors++;
        if (ca.hbond_is_acceptor) acceptors++;
    }

    std::cout << "  ORCA protein HBond summary:\n"
              << "    atoms=" << conf.AtomCount()
              << " residues=" << load.protein->ResidueCount() << "\n"
              << "    atoms with H-bond neighbours: " << has_hbond << "\n"
              << "    donors: " << donors << " acceptors: " << acceptors << "\n"
              << "    nearest geometry verified on " << checked << " atoms\n"
              << "    max |scalar| = " << max_scalar << " A^-3\n";

    EXPECT_GT(has_hbond, 0) << "Some atoms should have H-bond neighbours";
    EXPECT_GT(donors, 0) << "Should have donor atoms";
    EXPECT_GT(acceptors, 0) << "Should have acceptor atoms";
    EXPECT_GT(max_scalar, 0.001) << "Scalar should be non-zero";
}


// ============================================================================
// Batch test on 465 clean pairs: full analytical process (lesson 21)
// ============================================================================

static std::string FindNmrOutput(const std::string& dir,
                                  const std::string& prefix) {
    std::string exact = dir + prefix + "_nmr.out";
    if (fs::exists(exact)) return exact;
    for (const auto& entry : fs::directory_iterator(dir)) {
        std::string name = entry.path().filename().string();
        if (name.find(prefix) == 0 && name.find("_nmr.out") != std::string::npos)
            return entry.path().string();
    }
    return "";
}


TEST(BatchHBond, AllCleanPairs) {
    if (!fs::exists(nmr::test::TestEnvironment::Consolidated()))
        GTEST_SKIP() << "Consolidated directory not found";

    uint32_t saved_mask = OperationLog::GetChannelMask();
    OperationLog::SetChannelMask(0);

    int processed = 0, skipped = 0, failed = 0;
    double grand_max_scalar = 0;
    double grand_sum_scalar = 0;
    int grand_scalar_count = 0;

    // Per-secondary-structure: H-bond count in helix vs sheet vs coil
    int hb_in_helix = 0, hb_in_sheet = 0, hb_in_coil = 0;

    // DFT proximity: H-bond signal near vs far from mutation sites
    double sum_hb_near = 0, sum_hb_far = 0;
    int total_near = 0, total_far = 0;

    for (const auto& entry : fs::directory_iterator(nmr::test::TestEnvironment::Consolidated())) {
        if (!entry.is_directory()) continue;
        std::string protein_id = entry.path().filename().string();
        std::string dir = entry.path().string() + "/";

        // Require complete pair: WT + ALA with prmtop + xyz + nmr
        std::string wt_prmtop = dir + protein_id + "_WT.prmtop";
        std::string wt_xyz = dir + protein_id + "_WT.xyz";
        std::string ala_prmtop = dir + protein_id + "_ALA.prmtop";
        std::string ala_xyz = dir + protein_id + "_ALA.xyz";
        if (!fs::exists(wt_prmtop) || !fs::exists(wt_xyz) ||
            !fs::exists(ala_prmtop) || !fs::exists(ala_xyz)) {
            skipped++;
            continue;
        }

        std::string wt_nmr = FindNmrOutput(dir, protein_id + "_WT");
        std::string ala_nmr = FindNmrOutput(dir, protein_id + "_ALA");
        if (wt_nmr.empty() || ala_nmr.empty()) { skipped++; continue; }

        OrcaRunFiles files;
        files.pdb_path = dir + protein_id + "_WT.pdb";
        files.xyz_path = wt_xyz;
        files.prmtop_path = wt_prmtop;

        auto load = BuildFromOrca(files);
        if (!load.Ok()) { failed++; continue; }

        auto& conf = load.protein->Conformation();
        conf.AttachResult(GeometryResult::Compute(conf));

        PrmtopChargeSource cs(files.prmtop_path);
        conf.AttachResult(ChargeAssignmentResult::Compute(conf, cs));
        conf.AttachResult(SpatialIndexResult::Compute(conf));

        auto dssp = DsspResult::Compute(conf);
        if (!dssp) { failed++; continue; }
        conf.AttachResult(std::move(dssp));

        // All four calculators for T2 independence
        conf.AttachResult(McConnellResult::Compute(conf));
        conf.AttachResult(CoulombResult::Compute(conf));
        conf.AttachResult(RingSusceptibilityResult::Compute(conf));

        auto hbond = HBondResult::Compute(conf);
        if (!hbond) { failed++; continue; }
        conf.AttachResult(std::move(hbond));

        // ORCA shielding
        auto orca = OrcaShieldingResult::Compute(conf, wt_nmr);
        if (orca) conf.AttachResult(std::move(orca));

        const auto& dssp_data = conf.Result<DsspResult>();

        for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
            const auto& ca = conf.AtomAt(ai);

            double scalar = std::abs(ca.hbond_mcconnell_scalar);

            if (scalar > 1e-10) {
                grand_sum_scalar += scalar;
                grand_scalar_count++;
                grand_max_scalar = std::max(grand_max_scalar, scalar);
            }

            // Per-SS H-bond counting: for atoms that are donors
            if (ca.hbond_is_donor) {
                size_t ri = load.protein->AtomAt(ai).residue_index;
                char ss = dssp_data.SecondaryStructure(ri);
                if (ss == 'H' || ss == 'G' || ss == 'I') hb_in_helix++;
                else if (ss == 'E' || ss == 'B') hb_in_sheet++;
                else hb_in_coil++;
            }

        }

        // ------ DFT proximity: load ALA, compute MutationDelta ------
        // H-bonds are backbone — they persist across aromatic→ALA mutations.
        // So H-bond T0 should be spatially distributed (wherever there are
        // backbone H-bonds), NOT concentrated near mutation sites. This is
        // the complement to ring chi and Coulomb aromatic, which concentrate
        // near mutations.

        OrcaRunFiles ala_files;
        ala_files.pdb_path = dir + protein_id + "_ALA.pdb";
        ala_files.xyz_path = ala_xyz;
        ala_files.prmtop_path = ala_prmtop;

        auto ala_load = BuildFromOrca(ala_files);
        if (ala_load.Ok()) {
            auto& ala_conf = ala_load.protein->Conformation();
            ala_conf.AttachResult(GeometryResult::Compute(ala_conf));
            PrmtopChargeSource ala_cs(ala_files.prmtop_path);
            ala_conf.AttachResult(ChargeAssignmentResult::Compute(ala_conf, ala_cs));
            ala_conf.AttachResult(SpatialIndexResult::Compute(ala_conf));

            auto ala_orca = OrcaShieldingResult::Compute(ala_conf, ala_nmr);
            if (ala_orca) ala_conf.AttachResult(std::move(ala_orca));

            if (conf.HasResult<OrcaShieldingResult>() &&
                ala_conf.HasResult<OrcaShieldingResult>()) {

                auto delta = MutationDeltaResult::Compute(conf, ala_conf);
                if (delta) {
                    const auto& mut_sites = delta->MutationSites();
                    constexpr double TEST_NEAR_MUTATION_DIST = 8.0;

                    double sum_near = 0, sum_far = 0;
                    int nc = 0, fc = 0;

                    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
                        if (!delta->HasMatch(ai)) continue;
                        double hb_t0 = std::abs(
                            conf.AtomAt(ai).hbond_mcconnell_scalar);

                        bool near_mut = false;
                        Vec3 pos = conf.PositionAt(ai);
                        for (const auto& ms : mut_sites) {
                            const Residue& mr =
                                load.protein->ResidueAt(ms.residue_index);
                            if (mr.CA != Residue::NONE) {
                                double d = (pos - conf.PositionAt(mr.CA)).norm();
                                if (d < TEST_NEAR_MUTATION_DIST) {
                                    near_mut = true; break;
                                }
                            }
                        }

                        if (near_mut) { sum_near += hb_t0; nc++; }
                        else          { sum_far += hb_t0; fc++; }
                    }

                    if (nc > 0) { sum_hb_near += sum_near / nc; total_near++; }
                    if (fc > 0) { sum_hb_far += sum_far / fc; total_far++; }
                }
            }
        }

        processed++;
    }

    OperationLog::SetChannelMask(saved_mask);

    std::cout << "\n=== Batch HBondResult Summary ===\n"
              << "  Processed: " << processed
              << "  Skipped: " << skipped
              << "  Failed: " << failed << "\n\n";

    if (processed == 0) GTEST_SKIP() << "No proteins processed";

    double mean_scalar = grand_sum_scalar / std::max(grand_scalar_count, 1);

    std::cout << "  MAGNITUDES:\n"
              << "    Atoms with non-zero HBond scalar: " << grand_scalar_count << "\n"
              << "    Mean |scalar|: " << mean_scalar << " A^-3\n"
              << "    Max |scalar|:  " << grand_max_scalar << " A^-3\n\n";

    std::cout << "  H-BOND DONORS BY SECONDARY STRUCTURE:\n"
              << "    Helix: " << hb_in_helix << "\n"
              << "    Sheet: " << hb_in_sheet << "\n"
              << "    Coil:  " << hb_in_coil << "\n\n";

    // Assertions
    EXPECT_EQ(failed, 0);
    EXPECT_GT(processed, 100);

    // Physical magnitudes
    EXPECT_GT(grand_max_scalar, 0.001) << "H-bond scalar should be non-zero";

    // H-bonds should exist in all SS types
    EXPECT_GT(hb_in_helix, 0) << "Should have helix H-bond donors";
    EXPECT_GT(hb_in_sheet, 0) << "Should have sheet H-bond donors";

    // DFT proximity: H-bond signal near vs far from aromatic mutation sites.
    // Unlike ring chi and Coulomb aromatic (which concentrate near mutations),
    // H-bonds are backbone features that persist across aromatic→ALA mutations.
    // The near/far ratio should be close to 1 — H-bond signal is spatially
    // distributed, not concentrated where rings were removed.
    double grand_hb_near = (total_near > 0) ? sum_hb_near / total_near : 0;
    double grand_hb_far = (total_far > 0) ? sum_hb_far / total_far : 0;
    double hb_near_far_ratio = (grand_hb_far > 1e-10)
        ? grand_hb_near / grand_hb_far : 0;

    std::cout << "  DFT PROXIMITY (HBond |scalar| near vs far from mutation sites):\n"
              << "    Near mutation (<8A): " << grand_hb_near << " A^-3\n"
              << "    Far from mutation:   " << grand_hb_far << " A^-3\n"
              << "    Ratio near/far: " << hb_near_far_ratio
              << " (" << total_near << " proteins)\n\n";

    // H-bond near/far ratio should be much smaller than ring chi near/far
    // (which was 6.9x). H-bonds don't concentrate near aromatic mutations.
    // A ratio between 0.5 and 3.0 is reasonable — some variation from
    // packing effects but not the 7x signal that ring-specific calculators show.
    if (total_near > 10 && total_far > 10) {
        EXPECT_LT(hb_near_far_ratio, 5.0)
            << "H-bond near/far ratio should be much less than ring chi (6.9x) "
               "— H-bonds are not ring-specific";
    }
}

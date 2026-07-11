#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "HBondResult.h"
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

    // Seven residue slots give the remote control a minimum sequence
    // separation of three from the donor/acceptor endpoints.  Atom names are
    // intentionally empty: this is a geometry-kernel fixture, so the typed
    // semantic substrate's documented stub path applies.
    for (int ri = 0; ri < 7; ++ri) {
        Residue r;
        r.type = AminoAcid::Unknown;
        r.sequence_number = ri + 1;
        r.chain_id = "A";
        f.protein->AddResidue(std::move(r));
    }

    std::vector<Vec3> positions;
    auto add_atom = [&](size_t residue, Element element, const Vec3& pos) {
        auto atom = Atom::Create(element);
        atom->residue_index = residue;
        const size_t ai = f.protein->AddAtom(std::move(atom));
        f.protein->MutableResidueAt(residue).atom_indices.push_back(ai);
        positions.push_back(pos);
        return ai;
    };

    // Deliberately non-collinear N-H...O: the old N...O-midpoint proxy
    // cannot accidentally agree with the explicit-H production geometry.
    f.donor_n = add_atom(0, Element::N, Vec3(0.0, 0.0, 0.0));
    f.donor_h = add_atom(0, Element::H, Vec3(1.0, 0.0, 0.0));
    f.local_target = add_atom(1, Element::C, Vec3(3.0, 1.0, 0.0));
    f.acceptor_o = add_atom(3, Element::O, Vec3(1.0, 2.0, 0.0));
    f.remote_target = add_atom(6, Element::C, Vec3(2.0, 2.0, 0.0));

    auto& donor = f.protein->MutableResidueAt(0);
    donor.N = f.donor_n;
    donor.H = f.donor_h;
    f.protein->MutableResidueAt(3).O = f.acceptor_o;

    f.protein->FinalizeConstruction(positions);
    f.protein->AddConformation(std::move(positions), "explicit-H forcing fixture");
    return f;
}

}  // namespace


TEST(HBondGeometryKernel, UsesExplicitHydrogenAndAppliesTargetSequenceFilter) {
    auto f = BuildNonCollinearExplicitHydrogenFixture();
    auto& conf = f.protein->Conformation();

    ASSERT_TRUE(conf.AttachResult(GeometryResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(SpatialIndexResult::Compute(conf)));

    std::vector<DsspResidue> residues(f.protein->ResidueCount());
    residues[0].observed = true;
    residues[0].acceptors[0].residue_index = 3;
    ASSERT_TRUE(conf.AttachResult(DsspResult::CreateForTesting(std::move(residues))));

    auto hbond = HBondResult::Compute(conf);
    ASSERT_NE(hbond, nullptr);
    EXPECT_EQ(hbond->HBondCount(), 1u);
    ASSERT_TRUE(conf.AttachResult(std::move(hbond)));

    // Blessed values for donor H=(1,0,0), acceptor O=(1,2,0), and remote
    // target=(2,2,0).  These numbers are pinned independently, rather than
    // recomputing the production formula in the test.
    const auto& remote = conf.AtomAt(f.remote_target);
    EXPECT_NEAR(remote.hbond_nearest_dist, 2.236067977499790, 1e-12);
    EXPECT_NEAR(remote.hbond_inv_d3, 0.089442719099992, 1e-12);
    EXPECT_NEAR(remote.hbond_mcconnell_scalar, 0.125219806739988, 1e-12);
    EXPECT_NEAR(remote.hbond_nearest_dir.x(), 0.447213595499958, 1e-12);
    EXPECT_NEAR(remote.hbond_nearest_dir.y(), 0.894427190999916, 1e-12);
    EXPECT_NEAR(remote.hbond_nearest_dir.z(), 0.0, 1e-12);

    // Residue 1 is within the configured two-residue endpoint exclusion.
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

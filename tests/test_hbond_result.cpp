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

    // T0 = f identity: verify on nearest H-bond tensor
    int checked = 0;
    double max_diff = 0.0;
    int has_hbond = 0;
    double max_t0 = 0, max_t2 = 0;
    int donors = 0, acceptors = 0;

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& ca = conf.AtomAt(ai);

        if (ca.hbond_nearest_dist > 0 && ca.hbond_nearest_dist < NO_DATA_SENTINEL) {
            has_hbond++;

            // Verify T0 = f for nearest H-bond tensor
            double t0 = ca.hbond_nearest_spherical.T0;
            // Recompute f from the tensor trace
            double trace = ca.hbond_nearest_tensor.trace();
            double f_from_trace = trace / 3.0;
            double diff = std::abs(t0 - f_from_trace);
            max_diff = std::max(max_diff, diff);
            checked++;
        }

        max_t0 = std::max(max_t0, std::abs(ca.hbond_shielding_contribution.T0));
        max_t2 = std::max(max_t2, ca.hbond_shielding_contribution.T2Magnitude());
        if (ca.hbond_is_donor) donors++;
        if (ca.hbond_is_acceptor) acceptors++;
    }

    std::cout << "  ORCA protein HBond summary:\n"
              << "    atoms=" << conf.AtomCount()
              << " residues=" << load.protein->ResidueCount() << "\n"
              << "    atoms with H-bond neighbours: " << has_hbond << "\n"
              << "    donors: " << donors << " acceptors: " << acceptors << "\n"
              << "    T0=Trace/3 verified on " << checked
              << " nearest tensors, max diff = " << max_diff << "\n"
              << "    max |T0| = " << max_t0 << " A^-3\n"
              << "    max |T2| = " << max_t2 << " A^-3\n";

    EXPECT_GT(has_hbond, 0) << "Some atoms should have H-bond neighbours";
    EXPECT_GT(donors, 0) << "Should have donor atoms";
    EXPECT_GT(acceptors, 0) << "Should have acceptor atoms";
    EXPECT_LT(max_diff, 1e-10) << "T0 must equal Trace/3 at machine precision";
    EXPECT_GT(max_t0, 0.001) << "T0 should be non-zero";
    EXPECT_GT(max_t2, 0.001) << "T2 should be non-zero";
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

    // T2 independence accumulators
    auto t2_cos_sim = [](const std::array<double,5>& a,
                         const std::array<double,5>& b) -> double {
        double dot = 0, na = 0, nb = 0;
        for (int m = 0; m < 5; ++m) {
            dot += a[m] * b[m];
            na += a[m] * a[m];
            nb += b[m] * b[m];
        }
        double denom = std::sqrt(na * nb);
        return (denom < 1e-20) ? 0.0 : dot / denom;
    };

    constexpr double TEST_T2_MIN = 1e-4;

    struct T2Pair { double sum_abs_cos = 0; int count = 0; };
    T2Pair hb_vs_mc, hb_vs_co, hb_vs_rc;

    int processed = 0, skipped = 0, failed = 0;
    int total_hbonds = 0;
    double grand_max_t0 = 0, grand_max_t2 = 0;
    double grand_sum_t0 = 0, grand_sum_t2 = 0;
    int grand_t0_count = 0;

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

            double t0 = std::abs(ca.hbond_shielding_contribution.T0);
            double t2 = ca.hbond_shielding_contribution.T2Magnitude();

            if (t0 > 1e-10) {
                grand_sum_t0 += t0;
                grand_sum_t2 += t2;
                grand_t0_count++;
                grand_max_t0 = std::max(grand_max_t0, t0);
                grand_max_t2 = std::max(grand_max_t2, t2);
            }

            // Per-SS H-bond counting: for atoms that are donors
            if (ca.hbond_is_donor) {
                size_t ri = load.protein->AtomAt(ai).residue_index;
                char ss = dssp_data.SecondaryStructure(ri);
                if (ss == 'H' || ss == 'G' || ss == 'I') hb_in_helix++;
                else if (ss == 'E' || ss == 'B') hb_in_sheet++;
                else hb_in_coil++;
            }

            // T2 independence: HBond vs other three calculators
            double hb_t2_mag = ca.hbond_shielding_contribution.T2Magnitude();
            double mc_t2_mag = ca.mc_shielding_contribution.T2Magnitude();
            double co_t2_mag = ca.coulomb_shielding_contribution.T2Magnitude();
            double rc_t2_mag = ca.ringchi_shielding_contribution.T2Magnitude();

            if (hb_t2_mag > TEST_T2_MIN && mc_t2_mag > TEST_T2_MIN) {
                double c = t2_cos_sim(ca.hbond_shielding_contribution.T2,
                                      ca.mc_shielding_contribution.T2);
                hb_vs_mc.sum_abs_cos += std::abs(c);
                hb_vs_mc.count++;
            }
            if (hb_t2_mag > TEST_T2_MIN && co_t2_mag > TEST_T2_MIN) {
                double c = t2_cos_sim(ca.hbond_shielding_contribution.T2,
                                      ca.coulomb_shielding_contribution.T2);
                hb_vs_co.sum_abs_cos += std::abs(c);
                hb_vs_co.count++;
            }
            if (hb_t2_mag > TEST_T2_MIN && rc_t2_mag > TEST_T2_MIN) {
                double c = t2_cos_sim(ca.hbond_shielding_contribution.T2,
                                      ca.ringchi_shielding_contribution.T2);
                hb_vs_rc.sum_abs_cos += std::abs(c);
                hb_vs_rc.count++;
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
                            conf.AtomAt(ai).hbond_shielding_contribution.T0);

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

    double mean_t0 = grand_sum_t0 / std::max(grand_t0_count, 1);
    double mean_t2 = grand_sum_t2 / std::max(grand_t0_count, 1);

    std::cout << "  MAGNITUDES:\n"
              << "    Atoms with non-zero HBond T0: " << grand_t0_count << "\n"
              << "    Mean |T0|: " << mean_t0 << " A^-3\n"
              << "    Max |T0|:  " << grand_max_t0 << " A^-3\n"
              << "    Mean |T2|: " << mean_t2 << " A^-3\n"
              << "    Max |T2|:  " << grand_max_t2 << " A^-3\n"
              << "    T2/T0 ratio: " << ((mean_t0 > 1e-10) ? mean_t2/mean_t0 : 0) << "\n\n";

    std::cout << "  H-BOND DONORS BY SECONDARY STRUCTURE:\n"
              << "    Helix: " << hb_in_helix << "\n"
              << "    Sheet: " << hb_in_sheet << "\n"
              << "    Coil:  " << hb_in_coil << "\n\n";

    double mean_hb_mc = (hb_vs_mc.count > 0)
        ? hb_vs_mc.sum_abs_cos / hb_vs_mc.count : 0;
    double mean_hb_co = (hb_vs_co.count > 0)
        ? hb_vs_co.sum_abs_cos / hb_vs_co.count : 0;
    double mean_hb_rc = (hb_vs_rc.count > 0)
        ? hb_vs_rc.sum_abs_cos / hb_vs_rc.count : 0;

    std::cout << "  T2 INDEPENDENCE (mean |cos| in 5D):\n"
              << "    HBond vs McConnell:     " << mean_hb_mc
              << " (" << hb_vs_mc.count << " atoms)\n"
              << "    HBond vs Coulomb EFG:   " << mean_hb_co
              << " (" << hb_vs_co.count << " atoms)\n"
              << "    HBond vs Ring Chi:      " << mean_hb_rc
              << " (" << hb_vs_rc.count << " atoms)\n\n";

    // Assertions
    EXPECT_EQ(failed, 0);
    EXPECT_GT(processed, 100);

    // Physical magnitudes
    EXPECT_GT(grand_max_t0, 0.001) << "H-bond T0 should be non-zero";
    EXPECT_GT(grand_max_t2, 0.001) << "H-bond T2 should be non-zero";

    // T2 should exceed T0 (dipolar kernel)
    EXPECT_GT(mean_t2, mean_t0)
        << "T2 should exceed T0 for dipolar kernel";

    // H-bonds should exist in all SS types
    EXPECT_GT(hb_in_helix, 0) << "Should have helix H-bond donors";
    EXPECT_GT(hb_in_sheet, 0) << "Should have sheet H-bond donors";

    // T2 independence: HBond should not be parallel to other calculators
    if (hb_vs_mc.count > 1000) {
        EXPECT_LT(mean_hb_mc, 0.9)
            << "HBond and McConnell T2 should not be parallel";
    }

    // DFT proximity: H-bond signal near vs far from aromatic mutation sites.
    // Unlike ring chi and Coulomb aromatic (which concentrate near mutations),
    // H-bonds are backbone features that persist across aromatic→ALA mutations.
    // The near/far ratio should be close to 1 — H-bond signal is spatially
    // distributed, not concentrated where rings were removed.
    double grand_hb_near = (total_near > 0) ? sum_hb_near / total_near : 0;
    double grand_hb_far = (total_far > 0) ? sum_hb_far / total_far : 0;
    double hb_near_far_ratio = (grand_hb_far > 1e-10)
        ? grand_hb_near / grand_hb_far : 0;

    std::cout << "  DFT PROXIMITY (HBond |T0| near vs far from mutation sites):\n"
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

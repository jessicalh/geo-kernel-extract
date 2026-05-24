#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>

#include "McConnellResult.h"
#include "MutationDeltaResult.h"
#include "OrcaRunLoader.h"
#include "OrcaShieldingResult.h"
#include "ChargeAssignmentResult.h"
#include "ChargeSource.h"
#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "OperationLog.h"

namespace fs = std::filesystem;
using namespace nmr;



// ============================================================================
// Find NMR output file — handles both "ID_WT_nmr.out" and
// "ID_WT_20260327_120049_nmr.out" naming conventions.
// ============================================================================

static std::string FindNmrOutput(const std::string& dir,
                                  const std::string& prefix) {
    // Try exact name first
    std::string exact = dir + prefix + "_nmr.out";
    if (fs::exists(exact)) return exact;

    // Try dated pattern: prefix_YYYYMMDD_HHMMSS_nmr.out
    for (const auto& entry : fs::directory_iterator(dir)) {
        std::string name = entry.path().filename().string();
        if (name.find(prefix) == 0 && name.find("_nmr.out") != std::string::npos)
            return entry.path().string();
    }
    return "";
}


// ============================================================================
// Load a protein from consolidated directory (prmtop + xyz path)
// ============================================================================

struct BatchProtein {
    std::unique_ptr<Protein> protein;
    bool ok = false;
    std::string error;
};

static BatchProtein LoadFromConsolidated(
        const std::string& dir,
        const std::string& protein_id,
        const std::string& variant) {

    std::string prefix = protein_id + "_" + variant;

    OrcaRunFiles files;
    files.pdb_path = dir + prefix + ".pdb";
    files.xyz_path = dir + prefix + ".xyz";
    files.prmtop_path = dir + prefix + ".prmtop";
    files.tleap_script_path = dir + prefix + "_tleap.in";

    if (!fs::exists(files.xyz_path))
        return {nullptr, false, "xyz not found: " + files.xyz_path};
    if (!fs::exists(files.prmtop_path))
        return {nullptr, false, "prmtop not found: " + files.prmtop_path};

    auto load = BuildFromOrca(files);
    if (!load.Ok())
        return {nullptr, false, "BuildFromOrca failed: " + load.error};

    auto& conf = load.protein->Conformation();

    // Attach Layer 0 results needed for McConnell + MutationDelta
    conf.AttachResult(GeometryResult::Compute(conf));

    PrmtopChargeSource charge_source(files.prmtop_path);
    conf.AttachResult(ChargeAssignmentResult::Compute(conf, charge_source));

    conf.AttachResult(SpatialIndexResult::Compute(conf));

    // McConnell
    auto mc = McConnellResult::Compute(conf);
    if (!mc) return {nullptr, false, "McConnellResult failed"};
    conf.AttachResult(std::move(mc));

    // ORCA shielding
    std::string nmr_path = FindNmrOutput(dir, prefix);
    if (!nmr_path.empty()) {
        auto orca = OrcaShieldingResult::Compute(conf, nmr_path);
        if (orca) conf.AttachResult(std::move(orca));
    }

    return {std::move(load.protein), true, ""};
}


// ============================================================================
// Batch test: load all clean pairs, compute McConnell + MutationDelta,
// collect statistics.
// ============================================================================

TEST(BatchMcConnell, AllCleanPairs) {
    if (!fs::exists(nmr::test::TestEnvironment::Consolidated())) {
        GTEST_SKIP() << "Consolidated directory not found";
    }

    // Suppress per-atom logging for batch run
    uint32_t saved_mask = OperationLog::GetChannelMask();
    OperationLog::SetChannelMask(0);

    struct PairResult {
        std::string id;
        int wt_atoms = 0;
        int ala_atoms = 0;
        int matched = 0;
        int mutations = 0;
        double mean_abs_delta_t0 = 0;
        double max_abs_delta_t0 = 0;
        double mean_t2_magnitude = 0;
        double wt_max_mc_t0 = 0;
        double wt_max_mc_t2 = 0;
    };

    std::vector<PairResult> results;
    int skipped = 0;
    int failed = 0;

    for (const auto& entry : fs::directory_iterator(nmr::test::TestEnvironment::Consolidated())) {
        if (!entry.is_directory()) continue;
        std::string protein_id = entry.path().filename().string();
        std::string dir = entry.path().string() + "/";

        // Check clean path exists
        std::string wt_prmtop = dir + protein_id + "_WT.prmtop";
        std::string wt_xyz = dir + protein_id + "_WT.xyz";
        std::string ala_prmtop = dir + protein_id + "_ALA.prmtop";
        std::string ala_xyz = dir + protein_id + "_ALA.xyz";

        if (!fs::exists(wt_prmtop) || !fs::exists(wt_xyz) ||
            !fs::exists(ala_prmtop) || !fs::exists(ala_xyz)) {
            skipped++;
            continue;
        }

        // Check NMR output exists for both
        std::string wt_nmr = FindNmrOutput(dir, protein_id + "_WT");
        std::string ala_nmr = FindNmrOutput(dir, protein_id + "_ALA");
        if (wt_nmr.empty() || ala_nmr.empty()) {
            skipped++;
            continue;
        }

        // Load both
        auto wt = LoadFromConsolidated(dir, protein_id, "WT");
        auto ala = LoadFromConsolidated(dir, protein_id, "ALA");

        if (!wt.ok || !ala.ok) {
            failed++;
            if (failed <= 5) {
                std::cerr << "  FAIL " << protein_id << ": "
                          << (wt.ok ? "" : "WT: " + wt.error + " ")
                          << (ala.ok ? "" : "ALA: " + ala.error) << "\n";
            }
            continue;
        }

        auto& wt_conf = wt.protein->Conformation();
        const auto& ala_conf = ala.protein->Conformation();

        // Both must have ORCA shielding for MutationDelta
        if (!wt_conf.HasResult<OrcaShieldingResult>() ||
            !ala_conf.HasResult<OrcaShieldingResult>()) {
            failed++;
            continue;
        }

        // Compute MutationDelta
        auto delta = MutationDeltaResult::Compute(wt_conf, ala_conf);
        if (!delta) { failed++; continue; }

        // Collect statistics
        PairResult pr;
        pr.id = protein_id;
        pr.wt_atoms = static_cast<int>(wt_conf.AtomCount());
        pr.ala_atoms = static_cast<int>(ala_conf.AtomCount());
        pr.matched = static_cast<int>(delta->MatchedAtomCount());
        pr.mutations = static_cast<int>(delta->MutationSites().size());

        double sum_abs_t0 = 0, sum_t2 = 0;
        double max_abs = 0;
        int count = 0;
        for (size_t ai = 0; ai < wt_conf.AtomCount(); ++ai) {
            if (!delta->HasMatch(ai)) continue;
            double t0 = std::abs(delta->DeltaT0At(ai));
            sum_abs_t0 += t0;
            max_abs = std::max(max_abs, t0);
            sum_t2 += delta->DeltaShieldingSphericalAt(ai).T2Magnitude();
            count++;
        }
        if (count > 0) {
            pr.mean_abs_delta_t0 = sum_abs_t0 / count;
            pr.mean_t2_magnitude = sum_t2 / count;
        }
        pr.max_abs_delta_t0 = max_abs;

        // McConnell stats on WT
        double mc_max_t0 = 0, mc_max_t2 = 0;
        for (size_t ai = 0; ai < wt_conf.AtomCount(); ++ai) {
            const auto& sc = wt_conf.AtomAt(ai).mc_shielding_contribution;
            mc_max_t0 = std::max(mc_max_t0, std::abs(sc.T0));
            mc_max_t2 = std::max(mc_max_t2, sc.T2Magnitude());
        }
        pr.wt_max_mc_t0 = mc_max_t0;
        pr.wt_max_mc_t2 = mc_max_t2;

        results.push_back(pr);
    }

    OperationLog::SetChannelMask(saved_mask);

    // Print summary
    std::cout << "\n=== Batch McConnell + MutationDelta Summary ===\n";
    std::cout << "  Processed: " << results.size()
              << "  Skipped (missing files): " << skipped
              << "  Failed: " << failed << "\n\n";

    if (results.empty()) {
        GTEST_SKIP() << "No clean pairs processed";
    }

    // Aggregate statistics
    double total_mean_t0 = 0, total_mean_t2 = 0;
    double global_max_t0 = 0;
    int total_matched = 0;
    for (const auto& pr : results) {
        total_mean_t0 += pr.mean_abs_delta_t0;
        total_mean_t2 += pr.mean_t2_magnitude;
        global_max_t0 = std::max(global_max_t0, pr.max_abs_delta_t0);
        total_matched += pr.matched;
    }
    int n = static_cast<int>(results.size());

    std::cout << "  Grand mean |delta T0|: " << total_mean_t0 / n << " ppm\n";
    std::cout << "  Grand mean |T2|: " << total_mean_t2 / n << "\n";
    std::cout << "  Global max |delta T0|: " << global_max_t0 << " ppm\n";
    std::cout << "  Total matched atoms: " << total_matched << "\n";

    // McConnell stats
    double mc_mean_t0 = 0, mc_mean_t2 = 0;
    for (const auto& pr : results) {
        mc_mean_t0 += pr.wt_max_mc_t0;
        mc_mean_t2 += pr.wt_max_mc_t2;
    }
    std::cout << "  Mean max |McConnell T0|: " << mc_mean_t0 / n << " A^-3\n";
    std::cout << "  Mean max |McConnell T2|: " << mc_mean_t2 / n << " A^-3\n";

    // Print first 10 proteins
    std::cout << "\n  First 10 proteins:\n";
    for (int i = 0; i < std::min(10, n); ++i) {
        const auto& pr = results[i];
        std::cout << "    " << pr.id
                  << " WT=" << pr.wt_atoms << " ALA=" << pr.ala_atoms
                  << " matched=" << pr.matched
                  << " mutations=" << pr.mutations
                  << " mean|dT0|=" << pr.mean_abs_delta_t0
                  << " max|dT0|=" << pr.max_abs_delta_t0
                  << " MC_T0=" << pr.wt_max_mc_t0
                  << " MC_T2=" << pr.wt_max_mc_t2 << "\n";
    }

    // Assertions
    EXPECT_GT(results.size(), 100u)
        << "Should process at least 100 clean pairs";
    EXPECT_GT(total_matched, 10000)
        << "Should have many matched atoms across all pairs";
    EXPECT_GT(total_mean_t0 / n, 0.01)
        << "Grand mean delta T0 should be detectable";
}

#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include <filesystem>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <numeric>
#include <iomanip>

#include "CoulombResult.h"
#include "RingSusceptibilityResult.h"
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
// Find NMR output file
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


// ============================================================================
// Load and enrich a protein from consolidated directory
// ============================================================================

struct BatchProtein {
    std::unique_ptr<Protein> protein;
    bool ok = false;
    std::string error;
};

static BatchProtein LoadAndEnrich(const std::string& dir,
                                   const std::string& protein_id,
                                   const std::string& variant) {

    std::string prefix = protein_id + "_" + variant;

    OrcaRunFiles files;
    files.pdb_path = dir + prefix + ".pdb";
    files.xyz_path = dir + prefix + ".xyz";
    files.prmtop_path = dir + prefix + ".prmtop";

    if (!fs::exists(files.xyz_path))
        return {nullptr, false, "xyz not found"};
    if (!fs::exists(files.prmtop_path))
        return {nullptr, false, "prmtop not found"};

    auto load = BuildFromOrca(files);
    if (!load.Ok())
        return {nullptr, false, "BuildFromOrca: " + load.error};

    auto& conf = load.protein->Conformation();

    conf.AttachResult(GeometryResult::Compute(conf));

    PrmtopChargeSource charge_source(files.prmtop_path);
    conf.AttachResult(ChargeAssignmentResult::Compute(conf, charge_source));

    conf.AttachResult(SpatialIndexResult::Compute(conf));

    // McConnell (already batch-validated)
    auto mc = McConnellResult::Compute(conf);
    if (!mc) return {nullptr, false, "McConnellResult failed"};
    conf.AttachResult(std::move(mc));

    // Coulomb
    auto coulomb = CoulombResult::Compute(conf);
    if (!coulomb) return {nullptr, false, "CoulombResult failed"};
    conf.AttachResult(std::move(coulomb));

    // Ring susceptibility
    auto rchi = RingSusceptibilityResult::Compute(conf);
    if (!rchi) return {nullptr, false, "RingSusceptibilityResult failed"};
    conf.AttachResult(std::move(rchi));

    // ORCA shielding
    std::string nmr_path = FindNmrOutput(dir, prefix);
    if (!nmr_path.empty()) {
        auto orca = OrcaShieldingResult::Compute(conf, nmr_path);
        if (orca) conf.AttachResult(std::move(orca));
    }

    return {std::move(load.protein), true, ""};
}


// ============================================================================
// Batch test: all 465 clean pairs
// ============================================================================

TEST(BatchCoulombRingChi, AllCleanPairs) {
    if (!fs::exists(nmr::test::TestEnvironment::Consolidated())) {
        GTEST_SKIP() << "Consolidated directory not found";
    }

    uint32_t saved_mask = OperationLog::GetChannelMask();
    OperationLog::SetChannelMask(0);

    // Per-ring-type counters across the whole batch
    struct RingTypeStats {
        int count = 0;                 // number of rings of this type
        int atom_pairs = 0;            // total atom-ring pairs
        double sum_t0 = 0;            // sum of |T0| for mean
        double sum_t2 = 0;            // sum of |T2| for mean
        double max_t0 = 0;
        double max_t2 = 0;
    };
    std::array<RingTypeStats, 8> ring_type_stats = {};

    struct PairStats {
        std::string id;
        int wt_atoms = 0;
        int wt_rings = 0;
        int ala_rings = 0;

        // Coulomb stats (WT)
        double total_charge = 0;
        double mean_E_mag = 0;
        double max_E_mag = 0;
        double max_efg_trace = 0;
        double max_E_decomp_err = 0;
        double mean_backbone_proj = 0;

        // Ring susceptibility stats (WT)
        int ring_neighbour_pairs = 0;
        double max_t0_f_diff = 0;
        double max_rchi_t0 = 0;
        double max_rchi_t2 = 0;

        // DFT proximity: signal near mutation sites vs far
        double mean_rchi_t0_near_mutation = 0;
        double mean_rchi_t0_far = 0;
        double mean_coulomb_arom_E_near = 0;
        double mean_coulomb_arom_E_far = 0;
        int near_count = 0;
        int far_count = 0;

        // Comparison: WT vs ALA
        int wt_atoms_with_ring_neighbours = 0;
        int ala_atoms_with_ring_neighbours = 0;
    };

    // T2 independence analysis: cosine similarity between calculator T2
    // components in 5D space. If two calculators produce the same angular
    // pattern (cos ~ ±1), the model cannot distinguish them and one is
    // redundant. If cos ~ 0, they provide independent angular information.
    struct T2PairAccum {
        double sum_abs_cos = 0;  // sum of |cos(angle)| for mean
        int count = 0;
    };
    T2PairAccum mc_vs_coulomb, mc_vs_rchi, coulomb_vs_rchi;

    auto t2_cos_sim = [](const std::array<double,5>& a,
                         const std::array<double,5>& b) -> double {
        double dot = 0, na = 0, nb = 0;
        for (int m = 0; m < 5; ++m) {
            dot += a[m] * b[m];
            na += a[m] * a[m];
            nb += b[m] * b[m];
        }
        double denom = std::sqrt(na * nb);
        if (denom < 1e-20) return 0.0;  // one or both are zero
        return dot / denom;
    };

    // Test threshold: minimum T2 magnitude to include in independence check.
    // Below this, the T2 direction is dominated by numerical noise.
    constexpr double TEST_T2_MIN_FOR_INDEPENDENCE = 1e-4;

    std::vector<PairStats> results;
    int skipped = 0;
    int failed = 0;

    for (const auto& entry : fs::directory_iterator(nmr::test::TestEnvironment::Consolidated())) {
        if (!entry.is_directory()) continue;
        std::string protein_id = entry.path().filename().string();
        std::string dir = entry.path().string() + "/";

        // Check clean path exists (prmtop + xyz + nmr for both)
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
        if (wt_nmr.empty() || ala_nmr.empty()) {
            skipped++;
            continue;
        }

        // Load both
        auto wt = LoadAndEnrich(dir, protein_id, "WT");
        auto ala = LoadAndEnrich(dir, protein_id, "ALA");

        if (!wt.ok || !ala.ok) {
            failed++;
            if (failed <= 5)
                std::cerr << "  FAIL " << protein_id << ": "
                          << (wt.ok ? "" : "WT: " + wt.error + " ")
                          << (ala.ok ? "" : "ALA: " + ala.error) << "\n";
            continue;
        }

        auto& wt_conf = wt.protein->Conformation();
        auto& ala_conf = ala.protein->Conformation();

        PairStats ps;
        ps.id = protein_id;
        ps.wt_atoms = static_cast<int>(wt_conf.AtomCount());
        ps.wt_rings = static_cast<int>(wt.protein->RingCount());
        ps.ala_rings = static_cast<int>(ala.protein->RingCount());

        // ------ Coulomb validation on WT ------
        ps.total_charge = wt_conf.Result<ChargeAssignmentResult>().TotalCharge();

        double sum_E = 0, max_E = 0, max_trace = 0, max_decomp = 0;
        double sum_bb_proj = 0;
        for (size_t ai = 0; ai < wt_conf.AtomCount(); ++ai) {
            const auto& ca = wt_conf.AtomAt(ai);

            double E_mag = ca.coulomb_E_magnitude;
            sum_E += E_mag;
            max_E = std::max(max_E, E_mag);

            double trace = std::abs(ca.coulomb_EFG_total.trace());
            max_trace = std::max(max_trace, trace);

            // Decomposition check: |E_total - (E_bb + E_sc + E_arom)|
            Vec3 E_sum = ca.coulomb_E_backbone + ca.coulomb_E_sidechain
                       + ca.coulomb_E_aromatic;
            double decomp_err = (ca.coulomb_E_total - E_sum).norm();
            max_decomp = std::max(max_decomp, decomp_err);

            sum_bb_proj += ca.coulomb_E_backbone_frac;
        }
        ps.mean_E_mag = sum_E / wt_conf.AtomCount();
        ps.max_E_mag = max_E;
        ps.max_efg_trace = max_trace;
        ps.max_E_decomp_err = max_decomp;
        ps.mean_backbone_proj = sum_bb_proj / wt_conf.AtomCount();

        // ------ Ring susceptibility validation on WT ------
        double max_tf_diff = 0;
        double max_rchi_t0 = 0, max_rchi_t2 = 0;
        int ring_pairs = 0;
        int wt_with_rn = 0, ala_with_rn = 0;

        // Per-ring-type stats: count rings in this protein by type
        for (size_t ri = 0; ri < wt.protein->RingCount(); ++ri) {
            int type_idx = wt.protein->RingAt(ri).TypeIndexAsInt();
            if (type_idx >= 0 && type_idx < 8)
                ring_type_stats[type_idx].count++;
        }

        for (size_t ai = 0; ai < wt_conf.AtomCount(); ++ai) {
            const auto& ca = wt_conf.AtomAt(ai);
            if (!ca.ring_neighbours.empty()) wt_with_rn++;

            for (const auto& rn : ca.ring_neighbours) {
                if (std::abs(rn.chi_scalar) < 1e-15) continue;
                double diff = std::abs(rn.chi_spherical.T0 - rn.chi_scalar);
                max_tf_diff = std::max(max_tf_diff, diff);
                ring_pairs++;

                // Per-ring-type accumulation
                int ti = static_cast<int>(rn.ring_type);
                if (ti >= 0 && ti < 8) {
                    ring_type_stats[ti].atom_pairs++;
                    ring_type_stats[ti].sum_t0 += std::abs(rn.chi_spherical.T0);
                    ring_type_stats[ti].sum_t2 += rn.chi_spherical.T2Magnitude();
                    ring_type_stats[ti].max_t0 = std::max(
                        ring_type_stats[ti].max_t0, std::abs(rn.chi_spherical.T0));
                    ring_type_stats[ti].max_t2 = std::max(
                        ring_type_stats[ti].max_t2, rn.chi_spherical.T2Magnitude());
                }
            }

            max_rchi_t0 = std::max(max_rchi_t0,
                std::abs(ca.ringchi_shielding_contribution.T0));
            max_rchi_t2 = std::max(max_rchi_t2,
                ca.ringchi_shielding_contribution.T2Magnitude());
        }

        for (size_t ai = 0; ai < ala_conf.AtomCount(); ++ai) {
            if (!ala_conf.AtomAt(ai).ring_neighbours.empty()) ala_with_rn++;
        }

        ps.ring_neighbour_pairs = ring_pairs;
        ps.max_t0_f_diff = max_tf_diff;
        ps.max_rchi_t0 = max_rchi_t0;
        ps.max_rchi_t2 = max_rchi_t2;
        ps.wt_atoms_with_ring_neighbours = wt_with_rn;
        ps.ala_atoms_with_ring_neighbours = ala_with_rn;

        // ------ DFT comparison: ring chi signal vs DFT delta ------
        // Compute MutationDelta to get DFT deltas, then check whether
        // ring chi T0 is larger near mutation sites (where rings were removed)
        if (wt_conf.HasResult<OrcaShieldingResult>() &&
            ala_conf.HasResult<OrcaShieldingResult>()) {

            auto delta = MutationDeltaResult::Compute(wt_conf, ala_conf);
            if (delta) {
                // Identify atoms near mutation sites (within 5A of any
                // mutation site backbone atom)
                const auto& mut_sites = delta->MutationSites();

                double sum_rchi_near = 0, sum_rchi_far = 0;
                double sum_arom_E_near = 0, sum_arom_E_far = 0;
                int near_count = 0, far_count = 0;

                // Test threshold: atoms within 8A of mutation site CA
                // are "near." This is a test classification distance,
                // not a physics cutoff.
                constexpr double TEST_NEAR_MUTATION_DIST = 8.0;

                for (size_t ai = 0; ai < wt_conf.AtomCount(); ++ai) {
                    if (!delta->HasMatch(ai)) continue;

                    bool near_mutation = false;
                    Vec3 pos = wt_conf.PositionAt(ai);
                    for (const auto& ms : mut_sites) {
                        const Residue& mut_res =
                            wt.protein->ResidueAt(ms.residue_index);
                        if (mut_res.CA != Residue::NONE) {
                            double d = (pos - wt_conf.PositionAt(mut_res.CA)).norm();
                            if (d < TEST_NEAR_MUTATION_DIST) {
                                near_mutation = true; break;
                            }
                        }
                    }

                    double rchi_t0 = std::abs(
                        wt_conf.AtomAt(ai).ringchi_shielding_contribution.T0);
                    double arom_E = wt_conf.AtomAt(ai).aromatic_E_magnitude;

                    if (near_mutation) {
                        sum_rchi_near += rchi_t0;
                        sum_arom_E_near += arom_E;
                        near_count++;
                    } else {
                        sum_rchi_far += rchi_t0;
                        sum_arom_E_far += arom_E;
                        far_count++;
                    }
                }

                ps.near_count = near_count;
                ps.far_count = far_count;
                if (near_count > 0) {
                    ps.mean_rchi_t0_near_mutation = sum_rchi_near / near_count;
                    ps.mean_coulomb_arom_E_near = sum_arom_E_near / near_count;
                }
                if (far_count > 0) {
                    ps.mean_rchi_t0_far = sum_rchi_far / far_count;
                    ps.mean_coulomb_arom_E_far = sum_arom_E_far / far_count;
                }
            }
        }

        // ------ T2 independence: pairwise cosine similarity ------
        for (size_t ai = 0; ai < wt_conf.AtomCount(); ++ai) {
            const auto& ca = wt_conf.AtomAt(ai);

            double mc_t2_mag = ca.mc_shielding_contribution.T2Magnitude();
            double co_t2_mag = ca.coulomb_shielding_contribution.T2Magnitude();
            double rc_t2_mag = ca.ringchi_shielding_contribution.T2Magnitude();

            if (mc_t2_mag > TEST_T2_MIN_FOR_INDEPENDENCE &&
                co_t2_mag > TEST_T2_MIN_FOR_INDEPENDENCE) {
                double c = t2_cos_sim(ca.mc_shielding_contribution.T2,
                                      ca.coulomb_shielding_contribution.T2);
                mc_vs_coulomb.sum_abs_cos += std::abs(c);
                mc_vs_coulomb.count++;
            }

            if (mc_t2_mag > TEST_T2_MIN_FOR_INDEPENDENCE &&
                rc_t2_mag > TEST_T2_MIN_FOR_INDEPENDENCE) {
                double c = t2_cos_sim(ca.mc_shielding_contribution.T2,
                                      ca.ringchi_shielding_contribution.T2);
                mc_vs_rchi.sum_abs_cos += std::abs(c);
                mc_vs_rchi.count++;
            }

            if (co_t2_mag > TEST_T2_MIN_FOR_INDEPENDENCE &&
                rc_t2_mag > TEST_T2_MIN_FOR_INDEPENDENCE) {
                double c = t2_cos_sim(ca.coulomb_shielding_contribution.T2,
                                      ca.ringchi_shielding_contribution.T2);
                coulomb_vs_rchi.sum_abs_cos += std::abs(c);
                coulomb_vs_rchi.count++;
            }
        }

        results.push_back(ps);
    }

    OperationLog::SetChannelMask(saved_mask);

    // ======================================================================
    // Report and validate
    // ======================================================================

    int n = static_cast<int>(results.size());
    std::cout << "\n=== Batch Coulomb + RingSusceptibility Summary ===\n";
    std::cout << "  Processed: " << n
              << "  Skipped: " << skipped
              << "  Failed: " << failed << "\n\n";

    if (results.empty()) {
        GTEST_SKIP() << "No clean pairs processed";
    }

    // --- Coulomb aggregates ---
    double grand_mean_E = 0, grand_max_E = 0;
    double global_max_trace = 0, global_max_decomp = 0;
    double grand_mean_bb_proj = 0;
    int non_integer_charge = 0;

    for (const auto& ps : results) {
        grand_mean_E += ps.mean_E_mag;
        grand_max_E = std::max(grand_max_E, ps.max_E_mag);
        global_max_trace = std::max(global_max_trace, ps.max_efg_trace);
        global_max_decomp = std::max(global_max_decomp, ps.max_E_decomp_err);
        grand_mean_bb_proj += ps.mean_backbone_proj;

        double charge_frac = std::abs(ps.total_charge - std::round(ps.total_charge));
        if (charge_frac > 0.01) non_integer_charge++;
    }
    grand_mean_E /= n;
    grand_mean_bb_proj /= n;

    std::cout << "  COULOMB:\n"
              << "    Grand mean |E|: " << grand_mean_E << " V/A\n"
              << "    Global max |E|: " << grand_max_E << " V/A\n"
              << "    Max |EFG trace|: " << global_max_trace << "\n"
              << "    Max |E decomp err|: " << global_max_decomp << "\n"
              << "    Mean backbone proj: " << grand_mean_bb_proj << " V/A\n"
              << "    Non-integer charge: " << non_integer_charge << " / " << n << "\n\n";

    // --- Ring susceptibility aggregates ---
    double global_max_tf = 0;
    double grand_mean_rchi_t0 = 0, grand_mean_rchi_t2 = 0;
    int total_ring_pairs = 0;
    int wt_more_rn = 0;  // proteins where WT has more ring neighbours than ALA

    for (const auto& ps : results) {
        global_max_tf = std::max(global_max_tf, ps.max_t0_f_diff);
        grand_mean_rchi_t0 += ps.max_rchi_t0;
        grand_mean_rchi_t2 += ps.max_rchi_t2;
        total_ring_pairs += ps.ring_neighbour_pairs;
        if (ps.wt_atoms_with_ring_neighbours > ps.ala_atoms_with_ring_neighbours)
            wt_more_rn++;
    }
    grand_mean_rchi_t0 /= n;
    grand_mean_rchi_t2 /= n;

    std::cout << "  RING SUSCEPTIBILITY:\n"
              << "    Total ring-atom pairs: " << total_ring_pairs << "\n"
              << "    Max |T0 - f| across all: " << global_max_tf << "\n"
              << "    Mean max |T0|: " << grand_mean_rchi_t0 << " A^-3\n"
              << "    Mean max |T2|: " << grand_mean_rchi_t2 << " A^-3\n"
              << "    WT has more ring neighbours than ALA: "
              << wt_more_rn << " / " << n << "\n\n";

    // --- Per-ring-type breakdown ---
    std::cout << "  PER-RING-TYPE BREAKDOWN:\n";
    std::cout << "    Type       Rings   Pairs    Mean|T0|   Max|T0|   Mean|T2|   Max|T2|\n";
    for (int ti = 0; ti < 8; ++ti) {
        const auto& rts = ring_type_stats[ti];
        if (rts.count == 0) continue;
        double mean_t0 = (rts.atom_pairs > 0) ? rts.sum_t0 / rts.atom_pairs : 0;
        double mean_t2 = (rts.atom_pairs > 0) ? rts.sum_t2 / rts.atom_pairs : 0;
        std::cout << "    " << std::setw(8) << std::left
                  << RingTypeName(static_cast<RingTypeIndex>(ti))
                  << " " << std::setw(6) << std::right << rts.count
                  << " " << std::setw(8) << rts.atom_pairs
                  << "  " << std::setw(10) << std::fixed << std::setprecision(6) << mean_t0
                  << "  " << std::setw(8) << std::setprecision(4) << rts.max_t0
                  << "  " << std::setw(10) << std::setprecision(6) << mean_t2
                  << "  " << std::setw(8) << std::setprecision(4) << rts.max_t2
                  << "\n";
    }
    std::cout << "\n";

    // --- DFT proximity analysis ---
    // Calculator signals should be LARGER near mutation sites (where
    // aromatic sidechains were present in WT but removed in ALA).
    double sum_rchi_near = 0, sum_rchi_far = 0;
    double sum_arom_near = 0, sum_arom_far = 0;
    int total_near = 0, total_far = 0;
    for (const auto& ps : results) {
        if (ps.near_count > 0) {
            sum_rchi_near += ps.mean_rchi_t0_near_mutation;
            sum_arom_near += ps.mean_coulomb_arom_E_near;
            total_near++;
        }
        if (ps.far_count > 0) {
            sum_rchi_far += ps.mean_rchi_t0_far;
            sum_arom_far += ps.mean_coulomb_arom_E_far;
            total_far++;
        }
    }
    double grand_rchi_near = (total_near > 0) ? sum_rchi_near / total_near : 0;
    double grand_rchi_far = (total_far > 0) ? sum_rchi_far / total_far : 0;
    double grand_arom_near = (total_near > 0) ? sum_arom_near / total_near : 0;
    double grand_arom_far = (total_far > 0) ? sum_arom_far / total_far : 0;

    std::cout << "  DFT PROXIMITY (signal near vs far from mutation sites, <8A test threshold):\n"
              << "    Ring Chi |T0|:\n"
              << "      Near mutation: " << grand_rchi_near << " A^-3\n"
              << "      Far:           " << grand_rchi_far << " A^-3\n"
              << "      Ratio near/far: "
              << ((grand_rchi_far > 1e-10) ? grand_rchi_near / grand_rchi_far : 0) << "\n"
              << "    Coulomb aromatic |E|:\n"
              << "      Near mutation: " << grand_arom_near << " V/A\n"
              << "      Far:           " << grand_arom_far << " V/A\n"
              << "      Ratio near/far: "
              << ((grand_arom_far > 1e-10) ? grand_arom_near / grand_arom_far : 0) << "\n"
              << "    (" << total_near << " proteins with near/far data)\n\n";

    // Print first 10
    std::cout << "  First 10 proteins:\n";
    for (int i = 0; i < std::min(10, n); ++i) {
        const auto& ps = results[i];
        std::cout << "    " << ps.id
                  << " atoms=" << ps.wt_atoms
                  << " rings=" << ps.wt_rings << "/" << ps.ala_rings
                  << " mean|E|=" << ps.mean_E_mag
                  << " max|E|=" << ps.max_E_mag
                  << " q=" << ps.total_charge
                  << " rchi_T0=" << ps.max_rchi_t0
                  << " rchi_T2=" << ps.max_rchi_t2
                  << " rn_wt=" << ps.wt_atoms_with_ring_neighbours
                  << " rn_ala=" << ps.ala_atoms_with_ring_neighbours
                  << "\n";
    }

    // ======================================================================
    // Hard assertions
    // ======================================================================

    EXPECT_EQ(failed, 0) << failed << " proteins failed to process";

    EXPECT_GT(results.size(), 100u)
        << "Should process at least 100 clean pairs";

    // Coulomb: E-field magnitude in physically reasonable range
    // Case (1995): backbone amide H sees ~1-10 V/A. Mean across all
    // atoms (including distant ones) should be lower.
    EXPECT_GT(grand_mean_E, 0.1)
        << "Grand mean |E| too small — check ke multiplication";
    EXPECT_LT(grand_mean_E, 50.0)
        << "Grand mean |E| too large — check units";

    // EFG tracelessness: Gauss's law at machine precision
    EXPECT_LT(global_max_trace, 1e-8)
        << "EFG must be traceless everywhere";

    // Decomposition: exact vector sum
    EXPECT_LT(global_max_decomp, 1e-8)
        << "E decomposition must sum to total exactly";

    // Charges: prmtop gives integer charges
    EXPECT_EQ(non_integer_charge, 0)
        << "All prmtop charges should be integer total";

    // Ring susceptibility: T0 = f identity at machine precision
    EXPECT_LT(global_max_tf, 1e-10)
        << "T0 must equal f everywhere";

    // WT should have more ring neighbours than ALA (aromatic→ALA removes rings)
    EXPECT_GT(wt_more_rn, n / 2)
        << "Most WT proteins should have more ring neighbours than ALA mutant";

    // Ring susceptibility magnitudes physically reasonable
    // max |T0| ~ 0.01-5 A^-3 (depends on distance to ring)
    EXPECT_GT(grand_mean_rchi_t0, 0.001)
        << "Ring chi T0 too small";
    EXPECT_LT(grand_mean_rchi_t0, 50.0)
        << "Ring chi T0 too large";

    // Ring chi signal should be larger near mutation sites (where rings
    // were removed) than far away. If this fails, the 1/r³ spatial
    // decay is not working or the mutation site identification is wrong.
    if (total_near > 10 && total_far > 10) {
        EXPECT_GT(grand_rchi_near, grand_rchi_far)
            << "Ring chi signal should be stronger near mutation sites "
               "(near=" << grand_rchi_near << " far=" << grand_rchi_far << ")";

        // Coulomb aromatic E-field should also be larger near mutation
        // sites: the aromatic sidechain charges are the source, and
        // they're only present in WT (removed in ALA).
        EXPECT_GT(grand_arom_near, grand_arom_far)
            << "Aromatic E-field should be stronger near mutation sites "
               "(near=" << grand_arom_near << " far=" << grand_arom_far << ")";
    }

    // T2/T0 ratio: the dipolar kernel is pure T2 (traceless symmetric).
    // The T0 comes from asymmetric coupling with ring normal direction.
    // T2 should be larger than T0 for most atoms.
    EXPECT_GT(grand_mean_rchi_t2, grand_mean_rchi_t0)
        << "T2 should exceed T0 (dipolar kernel is predominantly T2)";

    // ======================================================================
    // T2 independence analysis
    // ======================================================================

    double mean_mc_co = (mc_vs_coulomb.count > 0)
        ? mc_vs_coulomb.sum_abs_cos / mc_vs_coulomb.count : 0;
    double mean_mc_rc = (mc_vs_rchi.count > 0)
        ? mc_vs_rchi.sum_abs_cos / mc_vs_rchi.count : 0;
    double mean_co_rc = (coulomb_vs_rchi.count > 0)
        ? coulomb_vs_rchi.sum_abs_cos / coulomb_vs_rchi.count : 0;

    std::cout << "  T2 INDEPENDENCE (mean |cos| in 5D, 0=orthogonal, 1=parallel):\n"
              << "    McConnell vs Coulomb EFG:     " << mean_mc_co
              << " (" << mc_vs_coulomb.count << " atoms)\n"
              << "    McConnell vs Ring Chi:        " << mean_mc_rc
              << " (" << mc_vs_rchi.count << " atoms)\n"
              << "    Coulomb EFG vs Ring Chi:      " << mean_co_rc
              << " (" << coulomb_vs_rchi.count << " atoms)\n\n";

    // If mean |cos| is near 1, the two calculators are producing the same
    // angular pattern and the model cannot distinguish them — one is
    // redundant. If near 0, they are orthogonal and maximally useful.
    // Random 5D vectors have expected |cos| = sqrt(2/(5*pi)) ≈ 0.36.
    // Values significantly below 0.9 indicate useful independence.
    if (mc_vs_coulomb.count > 1000) {
        EXPECT_LT(mean_mc_co, 0.9)
            << "McConnell and Coulomb T2 should not be parallel";
    }
    if (mc_vs_rchi.count > 1000) {
        EXPECT_LT(mean_mc_rc, 0.9)
            << "McConnell and RingChi T2 should not be parallel";
    }
    if (coulomb_vs_rchi.count > 1000) {
        EXPECT_LT(mean_co_rc, 0.9)
            << "Coulomb and RingChi T2 should not be parallel";
    }
}

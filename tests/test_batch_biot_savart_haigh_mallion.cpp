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
#include <array>

#include "BiotSavartResult.h"
#include "HaighMallionResult.h"
#include "McConnellResult.h"
#include "CoulombResult.h"
#include "RingSusceptibilityResult.h"
#include "HBondResult.h"
#include "DsspResult.h"
#include "MutationDeltaResult.h"
#include "OrcaRunLoader.h"
#include "OrcaShieldingResult.h"
#include "ChargeAssignmentResult.h"
#include "ChargeSource.h"
#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "OperationLog.h"
#include "PhysicalConstants.h"

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
// Load and run full pipeline in dependency order.
// Every calculator attached before BS and HM so the conformation has
// all prior results visible.
// ============================================================================

struct BatchProtein {
    std::unique_ptr<Protein> protein;
    bool ok = false;
    std::string error;
};

static BatchProtein LoadFullPipeline(const std::string& dir,
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

    // --- Foundation results ---
    conf.AttachResult(GeometryResult::Compute(conf));

    PrmtopChargeSource charge_source(files.prmtop_path);
    conf.AttachResult(ChargeAssignmentResult::Compute(conf, charge_source));

    conf.AttachResult(SpatialIndexResult::Compute(conf));

    auto dssp = DsspResult::Compute(conf);
    if (dssp) conf.AttachResult(std::move(dssp));

    // --- All classical calculators in order ---
    auto mc = McConnellResult::Compute(conf);
    if (!mc) return {nullptr, false, "McConnellResult failed"};
    conf.AttachResult(std::move(mc));

    auto coulomb = CoulombResult::Compute(conf);
    if (!coulomb) return {nullptr, false, "CoulombResult failed"};
    conf.AttachResult(std::move(coulomb));

    auto rchi = RingSusceptibilityResult::Compute(conf);
    if (!rchi) return {nullptr, false, "RingSusceptibilityResult failed"};
    conf.AttachResult(std::move(rchi));

    if (conf.HasResult<DsspResult>()) {
        auto hbond = HBondResult::Compute(conf);
        if (hbond) conf.AttachResult(std::move(hbond));
    }

    auto bs = BiotSavartResult::Compute(conf);
    if (!bs) return {nullptr, false, "BiotSavartResult failed"};
    conf.AttachResult(std::move(bs));

    auto hm = HaighMallionResult::Compute(conf);
    if (!hm) return {nullptr, false, "HaighMallionResult failed"};
    conf.AttachResult(std::move(hm));

    // --- ORCA shielding ---
    std::string nmr_path = FindNmrOutput(dir, prefix);
    if (!nmr_path.empty()) {
        auto orca = OrcaShieldingResult::Compute(conf, nmr_path);
        if (orca) conf.AttachResult(std::move(orca));
    }

    return {std::move(load.protein), true, ""};
}


// ============================================================================
// T2 cosine similarity in 5D
// ============================================================================

static double T2CosSim(const std::array<double,5>& a,
                       const std::array<double,5>& b) {
    double dot = 0, na = 0, nb = 0;
    for (int m = 0; m < 5; ++m) {
        dot += a[m] * b[m];
        na += a[m] * a[m];
        nb += b[m] * b[m];
    }
    double denom = std::sqrt(na * nb);
    if (denom < 1e-20) return 0.0;
    return dot / denom;
}

constexpr double T2_MIN_FOR_INDEPENDENCE = 1e-4;


// ============================================================================
// Distance bins for BS-HM convergence analysis
// ============================================================================

struct DistanceBin {
    double lo, hi;
    const char* label;
    int count = 0;
    double sum_bs_t0 = 0, sum_hm_t0 = 0;
    double sum_bs_t2 = 0, sum_hm_t2 = 0;
    double sum_abs_cos_t2 = 0;   // BS vs HM T2 cosine similarity
    int cos_count = 0;
    double sum_t0_ratio = 0;     // |HM_T0 / BS_T0| where both nonzero
    double sum_t0_ratio_sq = 0;
    int ratio_count = 0;
};


// ============================================================================
// Batch test: all 465 clean pairs, full pipeline, BS+HM analysis
// ============================================================================

TEST(BatchBiotSavartHaighMallion, AllCleanPairs) {
    if (!fs::exists(nmr::test::TestEnvironment::Consolidated())) {
        GTEST_SKIP() << "Consolidated directory not found";
    }

    uint32_t saved_mask = OperationLog::GetChannelMask();
    OperationLog::SetChannelMask(0);

    // --- Per-ring-type stats for BS and HM ---
    struct RingTypeStats {
        int count = 0;
        int atom_pairs = 0;
        double sum_bs_t0 = 0, sum_hm_t0 = 0;
        double sum_bs_t2 = 0, sum_hm_t2 = 0;
        double max_bs_t0 = 0, max_hm_t0 = 0;
        double max_bs_t2 = 0, max_hm_t2 = 0;
        // HM raw integral properties
        double max_hm_trace = 0;
        double max_hm_asym = 0;
    };
    std::array<RingTypeStats, 8> ring_stats = {};

    // --- Distance-binned BS-HM convergence ---
    std::array<DistanceBin, 5> dist_bins = {{
        {2.0, 4.0, "2-4A"},
        {4.0, 6.0, "4-6A"},
        {6.0, 8.0, "6-8A"},
        {8.0, 12.0, "8-12A"},
        {12.0, 15.0, "12-15A"}
    }};

    // --- T2 independence: all calculator pairs involving BS and HM ---
    struct T2PairAccum {
        double sum_abs_cos = 0;
        int count = 0;
    };
    T2PairAccum bs_vs_hm, bs_vs_mc, bs_vs_co, bs_vs_rchi, bs_vs_hb;
    T2PairAccum hm_vs_mc, hm_vs_co, hm_vs_rchi, hm_vs_hb;

    // --- Fused ring analysis (TRP): TRP5+TRP6 vs TRP9 ---
    struct FusedRingAccum {
        int trp_count = 0;           // TRP residues with all three rings
        double sum_t0_sum = 0;       // |T0(TRP5) + T0(TRP6)| per atom
        double sum_t0_perim = 0;     // |T0(TRP9)| per atom
        int atom_count = 0;
    };
    FusedRingAccum fused_bs, fused_hm;

    // --- DFT proximity ---
    double sum_bs_near = 0, sum_bs_far = 0;
    double sum_hm_near = 0, sum_hm_far = 0;
    int total_near = 0, total_far = 0;
    int proteins_with_prox = 0;

    int n_processed = 0, n_skipped = 0, n_failed = 0;

    for (const auto& entry : fs::directory_iterator(nmr::test::TestEnvironment::Consolidated())) {
        if (!entry.is_directory()) continue;
        std::string protein_id = entry.path().filename().string();
        std::string dir = entry.path().string() + "/";

        // Check all required files exist
        std::string wt_xyz = dir + protein_id + "_WT.xyz";
        std::string wt_prmtop = dir + protein_id + "_WT.prmtop";
        std::string ala_xyz = dir + protein_id + "_ALA.xyz";
        std::string ala_prmtop = dir + protein_id + "_ALA.prmtop";

        if (!fs::exists(wt_xyz) || !fs::exists(wt_prmtop) ||
            !fs::exists(ala_xyz) || !fs::exists(ala_prmtop)) {
            n_skipped++;
            continue;
        }

        std::string wt_nmr = FindNmrOutput(dir, protein_id + "_WT");
        std::string ala_nmr = FindNmrOutput(dir, protein_id + "_ALA");
        if (wt_nmr.empty() || ala_nmr.empty()) {
            n_skipped++;
            continue;
        }

        auto wt = LoadFullPipeline(dir, protein_id, "WT");
        auto ala = LoadFullPipeline(dir, protein_id, "ALA");

        if (!wt.ok || !ala.ok) {
            n_failed++;
            if (n_failed <= 5)
                std::cerr << "  FAIL " << protein_id << ": "
                          << (wt.ok ? "" : "WT: " + wt.error + " ")
                          << (ala.ok ? "" : "ALA: " + ala.error) << "\n";
            continue;
        }

        auto& wt_conf = wt.protein->Conformation();
        n_processed++;

        // ================================================================
        // Per-ring-type and distance-binned analysis on WT
        // ================================================================

        // Count rings by type
        for (size_t ri = 0; ri < wt.protein->RingCount(); ++ri) {
            int ti = wt.protein->RingAt(ri).TypeIndexAsInt();
            if (ti >= 0 && ti < 8) ring_stats[ti].count++;
        }

        // Build TRP fused ring map: residue -> {trp5_ri, trp6_ri, trp9_ri}
        struct TrpRings { size_t r5 = SIZE_MAX, r6 = SIZE_MAX, r9 = SIZE_MAX; };
        std::map<size_t, TrpRings> trp_map;
        for (size_t ri = 0; ri < wt.protein->RingCount(); ++ri) {
            const Ring& ring = wt.protein->RingAt(ri);
            size_t res = ring.parent_residue_index;
            if (ring.type_index == RingTypeIndex::TrpPyrrole)   trp_map[res].r5 = ri;
            if (ring.type_index == RingTypeIndex::TrpBenzene)   trp_map[res].r6 = ri;
            if (ring.type_index == RingTypeIndex::TrpPerimeter) trp_map[res].r9 = ri;
        }

        for (size_t ai = 0; ai < wt_conf.AtomCount(); ++ai) {
            const auto& ca = wt_conf.AtomAt(ai);

            for (const auto& rn : ca.ring_neighbours) {
                int ti = static_cast<int>(rn.ring_type);
                double dist = rn.distance_to_center;
                if (ti < 0 || ti >= 8) continue;

                double bs_t0 = std::abs(rn.G_spherical.T0);
                double bs_t2 = rn.G_spherical.T2Magnitude();
                double hm_t0_raw = std::abs(rn.hm_spherical.T0);  // should be ~0

                // Reconstruct HM full kernel G from stored data
                const RingGeometry& geom = wt_conf.ring_geometries[rn.ring_index];
                Vec3 V = rn.hm_tensor * geom.normal;
                double hm_full_t0 = std::abs(geom.normal.dot(V) / 3.0);

                SphericalTensor hm_G_st;
                {
                    Mat3 G_hm;
                    for (int a = 0; a < 3; ++a)
                        for (int b = 0; b < 3; ++b)
                            G_hm(a, b) = geom.normal(b) * V(a);
                    hm_G_st = SphericalTensor::Decompose(G_hm);
                }
                double hm_full_t2 = hm_G_st.T2Magnitude();

                // Per-ring-type accumulation
                ring_stats[ti].atom_pairs++;
                ring_stats[ti].sum_bs_t0 += bs_t0;
                ring_stats[ti].sum_hm_t0 += hm_full_t0;
                ring_stats[ti].sum_bs_t2 += bs_t2;
                ring_stats[ti].sum_hm_t2 += hm_full_t2;
                ring_stats[ti].max_bs_t0 = std::max(ring_stats[ti].max_bs_t0, bs_t0);
                ring_stats[ti].max_hm_t0 = std::max(ring_stats[ti].max_hm_t0, hm_full_t0);
                ring_stats[ti].max_bs_t2 = std::max(ring_stats[ti].max_bs_t2, bs_t2);
                ring_stats[ti].max_hm_t2 = std::max(ring_stats[ti].max_hm_t2, hm_full_t2);
                ring_stats[ti].max_hm_trace = std::max(
                    ring_stats[ti].max_hm_trace, std::abs(rn.hm_tensor.trace()));
                ring_stats[ti].max_hm_asym = std::max(
                    ring_stats[ti].max_hm_asym,
                    (rn.hm_tensor - rn.hm_tensor.transpose()).norm());

                // Distance-binned BS-HM convergence
                for (auto& bin : dist_bins) {
                    if (dist >= bin.lo && dist < bin.hi) {
                        bin.count++;
                        bin.sum_bs_t0 += bs_t0;
                        bin.sum_hm_t0 += hm_full_t0;
                        bin.sum_bs_t2 += bs_t2;
                        bin.sum_hm_t2 += hm_full_t2;

                        // T2 cosine similarity between BS and HM
                        if (bs_t2 > T2_MIN_FOR_INDEPENDENCE &&
                            hm_full_t2 > T2_MIN_FOR_INDEPENDENCE) {
                            double c = T2CosSim(rn.G_spherical.T2, hm_G_st.T2);
                            bin.sum_abs_cos_t2 += std::abs(c);
                            bin.cos_count++;
                        }

                        // T0 ratio where both are nonzero
                        if (bs_t0 > 1e-10 && hm_full_t0 > 1e-10) {
                            double ratio = hm_full_t0 / bs_t0;
                            bin.sum_t0_ratio += ratio;
                            bin.sum_t0_ratio_sq += ratio * ratio;
                            bin.ratio_count++;
                        }
                        break;
                    }
                }
            }

            // --- Fused ring analysis: TRP5+TRP6 vs TRP9 ---
            for (const auto& [res_idx, trp] : trp_map) {
                if (trp.r5 == SIZE_MAX || trp.r6 == SIZE_MAX || trp.r9 == SIZE_MAX)
                    continue;

                // Find this atom's RingNeighbourhood entries for each TRP ring
                const RingNeighbourhood* rn5 = nullptr;
                const RingNeighbourhood* rn6 = nullptr;
                const RingNeighbourhood* rn9 = nullptr;
                for (const auto& rn : ca.ring_neighbours) {
                    if (rn.ring_index == trp.r5) rn5 = &rn;
                    if (rn.ring_index == trp.r6) rn6 = &rn;
                    if (rn.ring_index == trp.r9) rn9 = &rn;
                }

                if (rn5 && rn6 && rn9) {
                    // BS: T0(TRP5) + T0(TRP6) vs T0(TRP9)
                    double bs_sum = rn5->G_spherical.T0 + rn6->G_spherical.T0;
                    fused_bs.sum_t0_sum += std::abs(bs_sum);
                    fused_bs.sum_t0_perim += std::abs(rn9->G_spherical.T0);
                    fused_bs.atom_count++;

                    // HM: same analysis using full kernel
                    auto hm_t0_for = [&](const RingNeighbourhood* rn) {
                        const RingGeometry& g = wt_conf.ring_geometries[rn->ring_index];
                        Vec3 v = rn->hm_tensor * g.normal;
                        return g.normal.dot(v) / 3.0;
                    };
                    double hm_sum = hm_t0_for(rn5) + hm_t0_for(rn6);
                    fused_hm.sum_t0_sum += std::abs(hm_sum);
                    fused_hm.sum_t0_perim += std::abs(hm_t0_for(rn9));
                    fused_hm.atom_count++;
                }
            }

            // --- T2 independence: BS and HM vs all other calculators ---
            double bs_t2_mag = ca.bs_shielding_contribution.T2Magnitude();
            double hm_t2_mag = ca.hm_shielding_contribution.T2Magnitude();
            double mc_t2_mag = ca.mc_shielding_contribution.T2Magnitude();
            double co_t2_mag = ca.coulomb_shielding_contribution.T2Magnitude();
            double rc_t2_mag = ca.ringchi_shielding_contribution.T2Magnitude();
            double hb_t2_mag = ca.hbond_shielding_contribution.T2Magnitude();

            auto accum = [&](T2PairAccum& acc, const std::array<double,5>& a,
                             const std::array<double,5>& b, double a_mag, double b_mag) {
                if (a_mag > T2_MIN_FOR_INDEPENDENCE &&
                    b_mag > T2_MIN_FOR_INDEPENDENCE) {
                    acc.sum_abs_cos += std::abs(T2CosSim(a, b));
                    acc.count++;
                }
            };

            // BS vs everything
            accum(bs_vs_hm, ca.bs_shielding_contribution.T2,
                  ca.hm_shielding_contribution.T2, bs_t2_mag, hm_t2_mag);
            accum(bs_vs_mc, ca.bs_shielding_contribution.T2,
                  ca.mc_shielding_contribution.T2, bs_t2_mag, mc_t2_mag);
            accum(bs_vs_co, ca.bs_shielding_contribution.T2,
                  ca.coulomb_shielding_contribution.T2, bs_t2_mag, co_t2_mag);
            accum(bs_vs_rchi, ca.bs_shielding_contribution.T2,
                  ca.ringchi_shielding_contribution.T2, bs_t2_mag, rc_t2_mag);
            accum(bs_vs_hb, ca.bs_shielding_contribution.T2,
                  ca.hbond_shielding_contribution.T2, bs_t2_mag, hb_t2_mag);

            // HM vs everything
            accum(hm_vs_mc, ca.hm_shielding_contribution.T2,
                  ca.mc_shielding_contribution.T2, hm_t2_mag, mc_t2_mag);
            accum(hm_vs_co, ca.hm_shielding_contribution.T2,
                  ca.coulomb_shielding_contribution.T2, hm_t2_mag, co_t2_mag);
            accum(hm_vs_rchi, ca.hm_shielding_contribution.T2,
                  ca.ringchi_shielding_contribution.T2, hm_t2_mag, rc_t2_mag);
            accum(hm_vs_hb, ca.hm_shielding_contribution.T2,
                  ca.hbond_shielding_contribution.T2, hm_t2_mag, hb_t2_mag);
        }

        // --- DFT proximity analysis ---
        auto& ala_conf = ala.protein->Conformation();
        if (wt_conf.HasResult<OrcaShieldingResult>() &&
            ala_conf.HasResult<OrcaShieldingResult>()) {

            auto delta = MutationDeltaResult::Compute(wt_conf, ala_conf);
            if (delta) {
                const auto& mut_sites = delta->MutationSites();
                constexpr double TEST_NEAR_DIST = 8.0;

                for (size_t ai = 0; ai < wt_conf.AtomCount(); ++ai) {
                    if (!delta->HasMatch(ai)) continue;

                    bool near = false;
                    Vec3 pos = wt_conf.PositionAt(ai);
                    for (const auto& ms : mut_sites) {
                        const Residue& mut_res =
                            wt.protein->ResidueAt(ms.residue_index);
                        if (mut_res.CA != Residue::NONE) {
                            double d = (pos - wt_conf.PositionAt(mut_res.CA)).norm();
                            if (d < TEST_NEAR_DIST) { near = true; break; }
                        }
                    }

                    double bs_t0 = std::abs(
                        wt_conf.AtomAt(ai).bs_shielding_contribution.T0);
                    double hm_t0 = std::abs(
                        wt_conf.AtomAt(ai).hm_shielding_contribution.T0);

                    if (near) {
                        sum_bs_near += bs_t0;
                        sum_hm_near += hm_t0;
                        total_near++;
                    } else {
                        sum_bs_far += bs_t0;
                        sum_hm_far += hm_t0;
                        total_far++;
                    }
                }
                proteins_with_prox++;
            }
        }
    }

    OperationLog::SetChannelMask(saved_mask);

    // ======================================================================
    // Report
    // ======================================================================

    int n = n_processed;
    std::cout << "\n=== Batch BiotSavart + HaighMallion Summary ===\n";
    std::cout << "  Processed: " << n
              << "  Skipped: " << n_skipped
              << "  Failed: " << n_failed << "\n\n";

    if (n == 0) {
        GTEST_SKIP() << "No clean pairs processed";
    }

    // --- Per-ring-type breakdown ---
    std::cout << "  PER-RING-TYPE BREAKDOWN:\n";
    std::cout << "    Type     Rings   Pairs    BS mean|T0|  HM mean|T0|  "
              << "BS mean|T2|  HM mean|T2|  HM max|Tr|   HM max asym\n";
    for (int ti = 0; ti < 8; ++ti) {
        const auto& s = ring_stats[ti];
        if (s.count == 0) continue;
        double bs_mt0 = (s.atom_pairs > 0) ? s.sum_bs_t0 / s.atom_pairs : 0;
        double hm_mt0 = (s.atom_pairs > 0) ? s.sum_hm_t0 / s.atom_pairs : 0;
        double bs_mt2 = (s.atom_pairs > 0) ? s.sum_bs_t2 / s.atom_pairs : 0;
        double hm_mt2 = (s.atom_pairs > 0) ? s.sum_hm_t2 / s.atom_pairs : 0;
        std::cout << "    " << std::setw(6) << std::left
                  << RingTypeName(static_cast<RingTypeIndex>(ti))
                  << " " << std::setw(6) << std::right << s.count
                  << " " << std::setw(8) << s.atom_pairs
                  << "  " << std::scientific << std::setprecision(3)
                  << std::setw(11) << bs_mt0
                  << "  " << std::setw(11) << hm_mt0
                  << "  " << std::setw(11) << bs_mt2
                  << "  " << std::setw(11) << hm_mt2
                  << "  " << std::setw(11) << s.max_hm_trace
                  << "  " << std::setw(11) << s.max_hm_asym
                  << "\n";
    }
    std::cout << "\n";

    // --- Distance-binned BS-HM convergence ---
    std::cout << "  BS-HM CONVERGENCE BY DISTANCE:\n";
    std::cout << "    Bin      Pairs    BS <|T0|>   HM <|T0|>  "
              << "T0 ratio(mean+/-sd)  T2 |cos|(BS,HM)\n";
    for (const auto& bin : dist_bins) {
        if (bin.count == 0) continue;
        double bs_mt0 = bin.sum_bs_t0 / bin.count;
        double hm_mt0 = bin.sum_hm_t0 / bin.count;
        double t2_cos = (bin.cos_count > 0)
            ? bin.sum_abs_cos_t2 / bin.cos_count : 0;
        double ratio_mean = (bin.ratio_count > 0)
            ? bin.sum_t0_ratio / bin.ratio_count : 0;
        double ratio_var = (bin.ratio_count > 1)
            ? (bin.sum_t0_ratio_sq / bin.ratio_count - ratio_mean * ratio_mean) : 0;
        double ratio_sd = (ratio_var > 0) ? std::sqrt(ratio_var) : 0;

        std::cout << "    " << std::setw(6) << std::left << bin.label
                  << " " << std::setw(8) << std::right << bin.count
                  << std::scientific << std::setprecision(3)
                  << "  " << std::setw(10) << bs_mt0
                  << "  " << std::setw(10) << hm_mt0
                  << std::fixed << std::setprecision(2)
                  << "  " << std::setw(8) << ratio_mean
                  << " +/- " << std::setw(6) << ratio_sd
                  << "  " << std::setprecision(4) << std::setw(8) << t2_cos
                  << " (" << bin.cos_count << ")"
                  << "\n";
    }
    std::cout << "\n";

    // --- Fused ring analysis ---
    if (fused_bs.atom_count > 0) {
        double bs_sum_mean = fused_bs.sum_t0_sum / fused_bs.atom_count;
        double bs_per_mean = fused_bs.sum_t0_perim / fused_bs.atom_count;
        double hm_sum_mean = fused_hm.sum_t0_sum / fused_hm.atom_count;
        double hm_per_mean = fused_hm.sum_t0_perim / fused_hm.atom_count;

        std::cout << "  FUSED RING (TRP): T0(TRP5)+T0(TRP6) vs T0(TRP9)\n"
                  << "    BS:  <|T0(5)+T0(6)|> = " << bs_sum_mean
                  << "  <|T0(9)|> = " << bs_per_mean
                  << "  ratio = " << (bs_per_mean > 1e-15 ? bs_sum_mean / bs_per_mean : 0) << "\n"
                  << "    HM:  <|T0(5)+T0(6)|> = " << hm_sum_mean
                  << "  <|T0(9)|> = " << hm_per_mean
                  << "  ratio = " << (hm_per_mean > 1e-15 ? hm_sum_mean / hm_per_mean : 0) << "\n"
                  << "    (" << fused_bs.atom_count << " atom-TRP triplets)\n\n";
    }

    // --- T2 independence ---
    auto report_t2 = [](const char* label, const T2PairAccum& acc) {
        double mean = (acc.count > 0) ? acc.sum_abs_cos / acc.count : 0;
        std::cout << "    " << std::setw(28) << std::left << label
                  << std::fixed << std::setprecision(4) << mean
                  << " (" << acc.count << " atoms)\n";
    };

    std::cout << "  T2 INDEPENDENCE (mean |cos| in 5D, random=0.36, parallel=1.0):\n";
    report_t2("BiotSavart vs HaighMallion:", bs_vs_hm);
    report_t2("BiotSavart vs McConnell:", bs_vs_mc);
    report_t2("BiotSavart vs Coulomb:", bs_vs_co);
    report_t2("BiotSavart vs RingSuscept:", bs_vs_rchi);
    report_t2("BiotSavart vs HBond:", bs_vs_hb);
    report_t2("HaighMallion vs McConnell:", hm_vs_mc);
    report_t2("HaighMallion vs Coulomb:", hm_vs_co);
    report_t2("HaighMallion vs RingSuscept:", hm_vs_rchi);
    report_t2("HaighMallion vs HBond:", hm_vs_hb);
    std::cout << "\n";

    // --- DFT proximity ---
    if (total_near > 0 && total_far > 0) {
        double bs_near = sum_bs_near / total_near;
        double bs_far = sum_bs_far / total_far;
        double hm_near = sum_hm_near / total_near;
        double hm_far = sum_hm_far / total_far;

        std::cout << "  DFT PROXIMITY (signal near vs far from mutation sites, <8A):\n"
                  << "    BiotSavart |T0|:\n"
                  << "      Near: " << bs_near << "  Far: " << bs_far
                  << "  Ratio: " << (bs_far > 1e-15 ? bs_near / bs_far : 0) << "\n"
                  << "    HaighMallion |T0|:\n"
                  << "      Near: " << hm_near << "  Far: " << hm_far
                  << "  Ratio: " << (hm_far > 1e-15 ? hm_near / hm_far : 0) << "\n"
                  << "    (" << proteins_with_prox << " proteins, "
                  << total_near << " near / " << total_far << " far atoms)\n\n";
    }

    // ======================================================================
    // Hard assertions
    // ======================================================================

    EXPECT_EQ(n_failed, 0) << n_failed << " proteins failed to process";
    EXPECT_GT(n, 100) << "Should process at least 100 clean pairs";

    // HM raw integral must be traceless everywhere
    for (int ti = 0; ti < 8; ++ti) {
        EXPECT_LT(ring_stats[ti].max_hm_trace, 1e-8)
            << "HM raw integral must be traceless for ring type " << ti;
        EXPECT_LT(ring_stats[ti].max_hm_asym, 1e-8)
            << "HM raw integral must be symmetric for ring type " << ti;
    }

    // BS and HM signal stronger near mutation sites
    if (total_near > 100 && total_far > 100) {
        double bs_near = sum_bs_near / total_near;
        double bs_far = sum_bs_far / total_far;
        double hm_near = sum_hm_near / total_near;
        double hm_far = sum_hm_far / total_far;

        EXPECT_GT(bs_near, bs_far)
            << "BS signal should be stronger near mutation sites";
        EXPECT_GT(hm_near, hm_far)
            << "HM signal should be stronger near mutation sites";
    }

    // T2 independence: BS and HM should not be parallel to existing calcs
    auto check_independent = [](const char* label, const T2PairAccum& acc) {
        if (acc.count > 1000) {
            double mean = acc.sum_abs_cos / acc.count;
            EXPECT_LT(mean, 0.9) << label << " T2 should not be parallel";
        }
    };
    check_independent("BS vs MC", bs_vs_mc);
    check_independent("BS vs Coulomb", bs_vs_co);
    check_independent("BS vs RingSuscept", bs_vs_rchi);
    check_independent("HM vs MC", hm_vs_mc);
    check_independent("HM vs Coulomb", hm_vs_co);
    check_independent("HM vs RingSuscept", hm_vs_rchi);

    // BS and HM model the same physics -- they SHOULD be correlated.
    // But the interesting question is how much. Report, don't assert.
}

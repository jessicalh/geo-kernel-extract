#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include <filesystem>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <iomanip>
#include <array>

#include "PiQuadrupoleResult.h"
#include "DispersionResult.h"
#include "McConnellResult.h"
#include "RingSusceptibilityResult.h"
#include "BiotSavartResult.h"
#include "CoulombResult.h"
#include "MutationDeltaResult.h"
#include "OrcaShieldingResult.h"
#include "ChargeAssignmentResult.h"
#include "ChargeSource.h"
#include "OrcaRunLoader.h"
#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "OperationLog.h"
#include "PhysicalConstants.h"

namespace fs = std::filesystem;
using namespace nmr;



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
// Load and enrich a protein with all 8 calculators + ORCA shielding
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

    // Prior calculators (needed for T2 independence comparison)
    conf.AttachResult(McConnellResult::Compute(conf));
    conf.AttachResult(CoulombResult::Compute(conf));
    conf.AttachResult(RingSusceptibilityResult::Compute(conf));
    conf.AttachResult(BiotSavartResult::Compute(conf));

    // New calculators
    auto pq = PiQuadrupoleResult::Compute(conf);
    if (!pq) return {nullptr, false, "PiQuadrupoleResult failed"};
    conf.AttachResult(std::move(pq));

    auto disp = DispersionResult::Compute(conf);
    if (!disp) return {nullptr, false, "DispersionResult failed"};
    conf.AttachResult(std::move(disp));

    // ORCA shielding (optional — needed for DFT proximity analysis)
    std::string nmr_path = FindNmrOutput(dir, prefix);
    if (!nmr_path.empty()) {
        auto orca = OrcaShieldingResult::Compute(conf, nmr_path);
        if (orca) conf.AttachResult(std::move(orca));
    }

    return {std::move(load.protein), true, ""};
}


// ============================================================================
// Batch test: PiQuadrupole + Dispersion on all 465 clean pairs
// ============================================================================

TEST(BatchPiQuadDisp, AllCleanPairs) {
    if (!fs::exists(nmr::test::TestEnvironment::Consolidated()))
        GTEST_SKIP() << "Consolidated directory not found";

    uint32_t saved_mask = OperationLog::GetChannelMask();
    OperationLog::SetChannelMask(0);

    // T2 independence: cosine similarity in 5D space
    struct T2PairAccum {
        double sum_abs_cos = 0;
        int count = 0;
    };
    // PQ vs: McConnell, Coulomb, RingChi, BS
    T2PairAccum pq_vs_mc, pq_vs_coulomb, pq_vs_rchi, pq_vs_bs;
    // Disp vs: McConnell, Coulomb, RingChi, BS, PQ
    T2PairAccum disp_vs_mc, disp_vs_coulomb, disp_vs_rchi, disp_vs_bs, disp_vs_pq;

    auto t2_cos_sim = [](const std::array<double,5>& a,
                         const std::array<double,5>& b) -> double {
        double dot = 0, na = 0, nb = 0;
        for (int m = 0; m < 5; ++m) {
            dot += a[m] * b[m];
            na += a[m] * a[m];
            nb += b[m] * b[m];
        }
        double denom = std::sqrt(na * nb);
        if (denom < 1e-20) return 0.0;
        return dot / denom;
    };

    constexpr double T2_MIN = 1e-4;

    // Per-ring-type stats for PQ
    struct RingTypeStats {
        int pairs = 0;
        double sum_scalar = 0;
        double sum_t2 = 0;
        double max_scalar = 0;
        double max_t2 = 0;
        double max_trace = 0;
    };
    std::array<RingTypeStats, 8> pq_type_stats = {};

    // DFT proximity accumulators
    constexpr double NEAR_MUTATION_DIST = 8.0;  // test threshold (A)
    double sum_pq_near = 0, sum_pq_far = 0;
    double sum_disp_near = 0, sum_disp_far = 0;
    int near_count = 0, far_count = 0;
    int pairs_with_dft = 0;

    int processed = 0, skipped = 0, failed = 0;
    int total_pq_pairs = 0, total_disp_pairs = 0, total_disp_contacts = 0;
    double global_max_pq_trace = 0;
    double global_max_pq_t2 = 0, global_max_disp_t2 = 0;

    for (const auto& entry : fs::directory_iterator(nmr::test::TestEnvironment::Consolidated())) {
        if (!entry.is_directory()) continue;
        std::string protein_id = entry.path().filename().string();
        std::string dir = entry.path().string() + "/";

        std::string wt_prmtop = dir + protein_id + "_WT.prmtop";
        std::string wt_xyz = dir + protein_id + "_WT.xyz";
        std::string ala_prmtop = dir + protein_id + "_ALA.prmtop";
        std::string ala_xyz = dir + protein_id + "_ALA.xyz";
        if (!fs::exists(wt_prmtop) || !fs::exists(wt_xyz)) {
            skipped++;
            continue;
        }

        auto wt = LoadAndEnrich(dir, protein_id, "WT");
        if (!wt.ok) {
            failed++;
            if (failed <= 5)
                std::cerr << "  FAIL " << protein_id << ": " << wt.error << "\n";
            continue;
        }

        auto& conf = wt.protein->Conformation();

        // --- PQ validation ---
        int pq_pairs = 0;
        for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
            const auto& ca = conf.AtomAt(ai);

            for (const auto& rn : ca.ring_neighbours) {
                if (rn.quad_tensor.isZero(1e-20)) continue;

                double trace = std::abs(rn.quad_tensor.trace());
                global_max_pq_trace = std::max(global_max_pq_trace, trace);

                int ti = static_cast<int>(rn.ring_type);
                if (ti >= 0 && ti < 8) {
                    pq_type_stats[ti].pairs++;
                    pq_type_stats[ti].sum_scalar += std::abs(rn.quad_scalar);
                    pq_type_stats[ti].sum_t2 += rn.quad_spherical.T2Magnitude();
                    pq_type_stats[ti].max_scalar = std::max(
                        pq_type_stats[ti].max_scalar, std::abs(rn.quad_scalar));
                    pq_type_stats[ti].max_t2 = std::max(
                        pq_type_stats[ti].max_t2, rn.quad_spherical.T2Magnitude());
                    pq_type_stats[ti].max_trace = std::max(
                        pq_type_stats[ti].max_trace, trace);
                }
                pq_pairs++;
            }

            global_max_pq_t2 = std::max(global_max_pq_t2,
                ca.piquad_shielding_contribution.T2Magnitude());

            // --- Disp stats ---
            for (const auto& rn : ca.ring_neighbours) {
                if (rn.disp_contacts > 0) {
                    total_disp_pairs++;
                    total_disp_contacts += rn.disp_contacts;
                }
            }
            global_max_disp_t2 = std::max(global_max_disp_t2,
                ca.disp_shielding_contribution.T2Magnitude());

            // --- T2 independence ---
            auto pq_t2 = ca.piquad_shielding_contribution.T2;
            auto disp_t2 = ca.disp_shielding_contribution.T2;
            auto mc_t2 = ca.mc_shielding_contribution.T2;
            auto coulomb_t2 = ca.coulomb_EFG_total_spherical.T2;
            auto rchi_t2 = ca.ringchi_shielding_contribution.T2;
            auto bs_t2 = ca.bs_shielding_contribution.T2;

            double pq_mag = ca.piquad_shielding_contribution.T2Magnitude();
            double disp_mag = ca.disp_shielding_contribution.T2Magnitude();
            double mc_mag = ca.mc_shielding_contribution.T2Magnitude();
            double coulomb_mag = ca.coulomb_EFG_total_spherical.T2Magnitude();
            double rchi_mag = ca.ringchi_shielding_contribution.T2Magnitude();
            double bs_mag = ca.bs_shielding_contribution.T2Magnitude();

            if (pq_mag > T2_MIN) {
                if (mc_mag > T2_MIN)      { pq_vs_mc.sum_abs_cos += std::abs(t2_cos_sim(pq_t2, mc_t2)); pq_vs_mc.count++; }
                if (coulomb_mag > T2_MIN)  { pq_vs_coulomb.sum_abs_cos += std::abs(t2_cos_sim(pq_t2, coulomb_t2)); pq_vs_coulomb.count++; }
                if (rchi_mag > T2_MIN)     { pq_vs_rchi.sum_abs_cos += std::abs(t2_cos_sim(pq_t2, rchi_t2)); pq_vs_rchi.count++; }
                if (bs_mag > T2_MIN)       { pq_vs_bs.sum_abs_cos += std::abs(t2_cos_sim(pq_t2, bs_t2)); pq_vs_bs.count++; }
            }
            if (disp_mag > T2_MIN) {
                if (mc_mag > T2_MIN)      { disp_vs_mc.sum_abs_cos += std::abs(t2_cos_sim(disp_t2, mc_t2)); disp_vs_mc.count++; }
                if (coulomb_mag > T2_MIN)  { disp_vs_coulomb.sum_abs_cos += std::abs(t2_cos_sim(disp_t2, coulomb_t2)); disp_vs_coulomb.count++; }
                if (rchi_mag > T2_MIN)     { disp_vs_rchi.sum_abs_cos += std::abs(t2_cos_sim(disp_t2, rchi_t2)); disp_vs_rchi.count++; }
                if (bs_mag > T2_MIN)       { disp_vs_bs.sum_abs_cos += std::abs(t2_cos_sim(disp_t2, bs_t2)); disp_vs_bs.count++; }
                if (pq_mag > T2_MIN)       { disp_vs_pq.sum_abs_cos += std::abs(t2_cos_sim(disp_t2, pq_t2)); disp_vs_pq.count++; }
            }
        }

        total_pq_pairs += pq_pairs;

        // --- DFT proximity analysis ---
        // Load ALA if available, compute MutationDelta, check whether
        // PQ/Disp signal is stronger near mutation sites.
        if (fs::exists(ala_prmtop) && fs::exists(ala_xyz) &&
            conf.HasResult<OrcaShieldingResult>()) {

            auto ala = LoadAndEnrich(dir, protein_id, "ALA");
            if (ala.ok && ala.protein->Conformation().HasResult<OrcaShieldingResult>()) {
                auto& ala_conf = ala.protein->Conformation();
                auto delta = MutationDeltaResult::Compute(conf, ala_conf);
                if (delta) {
                    const auto& mut_sites = delta->MutationSites();

                    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
                        if (!delta->HasMatch(ai)) continue;

                        bool near_mut = false;
                        Vec3 pos = conf.PositionAt(ai);
                        for (const auto& ms : mut_sites) {
                            const Residue& r = wt.protein->ResidueAt(ms.residue_index);
                            if (r.CA != Residue::NONE) {
                                if ((pos - conf.PositionAt(r.CA)).norm() < NEAR_MUTATION_DIST) {
                                    near_mut = true; break;
                                }
                            }
                        }

                        double pq_t2 = conf.AtomAt(ai).piquad_shielding_contribution.T2Magnitude();
                        double disp_t2 = conf.AtomAt(ai).disp_shielding_contribution.T2Magnitude();

                        if (near_mut) {
                            sum_pq_near += pq_t2;
                            sum_disp_near += disp_t2;
                            near_count++;
                        } else {
                            sum_pq_far += pq_t2;
                            sum_disp_far += disp_t2;
                            far_count++;
                        }
                    }
                    pairs_with_dft++;
                }
            }
        }

        processed++;
    }

    OperationLog::SetChannelMask(saved_mask);

    ASSERT_GT(processed, 400) << "Should process 400+ proteins";
    EXPECT_EQ(failed, 0) << "Zero failures expected";

    // PQ tracelessness across all proteins
    EXPECT_LT(global_max_pq_trace, 1e-10)
        << "PQ EFG must be traceless across all proteins";

    // --- Summary ---
    std::cout << "\n  === Batch PiQuadrupole + Dispersion ===\n"
              << "  Processed: " << processed
              << " Skipped: " << skipped
              << " Failed: " << failed << "\n\n";

    // PQ summary
    std::cout << "  PiQuadrupole:\n"
              << "    Total atom-ring pairs: " << total_pq_pairs << "\n"
              << "    Global max |Tr|: " << global_max_pq_trace << "\n"
              << "    Global max |T2|: " << global_max_pq_t2 << " A^-5\n\n";

    // PQ per-ring-type table
    const char* type_names[] = {
        "PHE", "TYR", "TRP6", "TRP5", "TRP9", "HIS", "HID", "HIE"
    };
    std::cout << "    Per-ring-type (PQ scalar and T2):\n"
              << "    " << std::setw(6) << "Type"
              << std::setw(10) << "Pairs"
              << std::setw(14) << "Mean|scalar|"
              << std::setw(14) << "Max|scalar|"
              << std::setw(14) << "Mean|T2|"
              << std::setw(14) << "Max|T2|"
              << std::setw(14) << "Max|Tr|" << "\n";
    for (int ti = 0; ti < 8; ++ti) {
        const auto& s = pq_type_stats[ti];
        if (s.pairs == 0) continue;
        std::cout << "    " << std::setw(6) << type_names[ti]
                  << std::setw(10) << s.pairs
                  << std::setw(14) << std::setprecision(6) << s.sum_scalar / s.pairs
                  << std::setw(14) << s.max_scalar
                  << std::setw(14) << s.sum_t2 / s.pairs
                  << std::setw(14) << s.max_t2
                  << std::setw(14) << std::scientific << s.max_trace
                  << std::fixed << "\n";
    }

    // Disp summary
    std::cout << "\n  Dispersion:\n"
              << "    Total atom-ring pairs with contacts: " << total_disp_pairs << "\n"
              << "    Total vertex contacts: " << total_disp_contacts << "\n"
              << "    Global max |T2|: " << global_max_disp_t2 << " A^-6\n\n";

    // T2 independence
    auto print_t2 = [](const char* name, const T2PairAccum& acc) {
        if (acc.count > 0)
            std::cout << "    " << std::setw(22) << name
                      << ": mean|cos| = " << std::setprecision(3)
                      << acc.sum_abs_cos / acc.count
                      << " (n=" << acc.count << ")\n";
    };

    std::cout << "  T2 independence (mean |cos| in 5D, random~0.36):\n";
    print_t2("PQ vs McConnell", pq_vs_mc);
    print_t2("PQ vs Coulomb", pq_vs_coulomb);
    print_t2("PQ vs RingChi", pq_vs_rchi);
    print_t2("PQ vs BiotSavart", pq_vs_bs);
    print_t2("Disp vs McConnell", disp_vs_mc);
    print_t2("Disp vs Coulomb", disp_vs_coulomb);
    print_t2("Disp vs RingChi", disp_vs_rchi);
    print_t2("Disp vs BiotSavart", disp_vs_bs);
    print_t2("Disp vs PiQuad", disp_vs_pq);

    // DFT proximity report
    if (near_count > 0 && far_count > 0) {
        double mean_pq_near = sum_pq_near / near_count;
        double mean_pq_far = sum_pq_far / far_count;
        double mean_disp_near = sum_disp_near / near_count;
        double mean_disp_far = sum_disp_far / far_count;

        std::cout << "\n  DFT proximity analysis (near < " << NEAR_MUTATION_DIST
                  << "A from mutation site CA):\n"
                  << "    Pairs with DFT: " << pairs_with_dft
                  << " Near atoms: " << near_count
                  << " Far atoms: " << far_count << "\n"
                  << "    PQ  mean|T2|: near=" << std::setprecision(4) << mean_pq_near
                  << " far=" << mean_pq_far
                  << " ratio=" << (mean_pq_far > 0 ? mean_pq_near / mean_pq_far : 0)
                  << "\n"
                  << "    Disp mean|T2|: near=" << mean_disp_near
                  << " far=" << mean_disp_far
                  << " ratio=" << (mean_disp_far > 0 ? mean_disp_near / mean_disp_far : 0)
                  << "\n";

        // PQ should be stronger near mutation sites (aromatic rings removed)
        // because the EFG from the ring quadrupole decays steeply (1/r^5).
        // Expect ratio > 1 (signal stronger where the rings were).
        EXPECT_GT(mean_pq_near / mean_pq_far, 1.0)
            << "PQ T2 should be stronger near mutation sites";

        // Dispersion should also be stronger near (1/r^6, even steeper).
        EXPECT_GT(mean_disp_near / mean_disp_far, 1.0)
            << "Disp T2 should be stronger near mutation sites";
    }
    std::cout << "\n";
}

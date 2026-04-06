#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include "PdbFileReader.h"
#include "Protein.h"
#include "ProtonationDetectionResult.h"
#include "GeometryResult.h"
#include "EnrichmentResult.h"
#include "ChargeAssignmentResult.h"
#include "SpatialIndexResult.h"
#include "DsspResult.h"
#include "ApbsFieldResult.h"
#include "MolecularGraphResult.h"
#include "MopacResult.h"
#include "RuntimeEnvironment.h"
#include "OperationLog.h"
#include <filesystem>
#include <cmath>
#include <map>

using namespace nmr;

// ============================================================================
// Test data paths
// ============================================================================




// ============================================================================
// Full pipeline integration test: loads 1UBQ, attaches ALL 9 results
// in dependency order, verifies each one.
// ============================================================================

class FullPipelineTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!std::filesystem::exists(nmr::test::TestEnvironment::UbqProtonated())) {
            GTEST_SKIP() << "1UBQ PDB not found at " << nmr::test::TestEnvironment::UbqProtonated();
        }
        auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
        if (!r.Ok()) GTEST_SKIP() << r.error;
        protein = std::move(r.protein);

        RuntimeEnvironment::Load();
    }
    std::unique_ptr<Protein> protein;
    int results_attached = 0;
    int results_failed = 0;
};


TEST_F(FullPipelineTest, EntirePipelineEndToEnd) {
    auto& conf = protein->Conformation();

    OperationLog::Info("FullPipelineTest",
        "Starting full pipeline on 1UBQ: " +
        std::to_string(conf.AtomCount()) + " atoms, " +
        std::to_string(protein->ResidueCount()) + " residues, " +
        std::to_string(protein->RingCount()) + " rings, " +
        std::to_string(protein->BondCount()) + " bonds");

    // ================================================================
    // Step 1: ProtonationDetectionResult (no dependencies)
    // 1UBQ crystal has no H atoms, so all titratable residues unresolved.
    // ================================================================
    {
        auto prot = ProtonationDetectionResult::Compute(conf);
        if (prot && conf.AttachResult(std::move(prot))) {
            results_attached++;
            OperationLog::Info("FullPipelineTest",
                "ProtonationDetectionResult: attached");
        } else {
            results_failed++;
            OperationLog::Error("FullPipelineTest",
                "ProtonationDetectionResult: FAILED to attach");
        }
    }

    // ================================================================
    // Step 2: GeometryResult (no dependencies)
    // ================================================================
    {
        auto geo = GeometryResult::Compute(conf);
        if (geo && conf.AttachResult(std::move(geo))) {
            results_attached++;
            OperationLog::Info("FullPipelineTest",
                "GeometryResult: attached");
        } else {
            results_failed++;
            OperationLog::Error("FullPipelineTest",
                "GeometryResult: FAILED to attach");
        }
    }

    // ================================================================
    // Step 3: EnrichmentResult (no dependencies)
    // ================================================================
    {
        auto enrich = EnrichmentResult::Compute(conf);
        if (enrich && conf.AttachResult(std::move(enrich))) {
            results_attached++;
            OperationLog::Info("FullPipelineTest",
                "EnrichmentResult: attached (unknowns=" +
                std::to_string(
                    conf.Result<EnrichmentResult>().UnknownCount()) + ")");
        } else {
            results_failed++;
            OperationLog::Error("FullPipelineTest",
                "EnrichmentResult: FAILED to attach");
        }
    }

    // ================================================================
    // Step 4: ChargeAssignmentResult (no dependencies, uses ff14SB file)
    // ================================================================
    bool charges_ok = false;
    {
        if (!std::filesystem::exists(nmr::test::TestEnvironment::Ff14sbParams())) {
            OperationLog::Warn("FullPipelineTest",
                "ff14sb_params.dat not found, using stub charges");
            auto charges = ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams().c_str());
            if (charges && conf.AttachResult(std::move(charges))) {
                results_attached++;
                charges_ok = true;
            } else {
                results_failed++;
            }
        } else {
            auto charges = ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams());
            if (charges && conf.AttachResult(std::move(charges))) {
                results_attached++;
                charges_ok = true;
                OperationLog::Info("FullPipelineTest",
                    "ChargeAssignmentResult: attached (source=" +
                    conf.Result<ChargeAssignmentResult>().Source() +
                    ", total_charge=" +
                    std::to_string(
                        conf.Result<ChargeAssignmentResult>().TotalCharge()) +
                    ")");
            } else {
                results_failed++;
                OperationLog::Error("FullPipelineTest",
                    "ChargeAssignmentResult: FAILED to attach");
            }
        }
    }

    // ================================================================
    // Step 5: SpatialIndexResult (depends on GeometryResult)
    // ================================================================
    {
        if (conf.HasResult<GeometryResult>()) {
            auto spatial = SpatialIndexResult::Compute(conf);
            if (spatial && conf.AttachResult(std::move(spatial))) {
                results_attached++;
                OperationLog::Info("FullPipelineTest",
                    "SpatialIndexResult: attached");
            } else {
                results_failed++;
                OperationLog::Error("FullPipelineTest",
                    "SpatialIndexResult: FAILED to attach");
            }
        } else {
            results_failed++;
            OperationLog::Error("FullPipelineTest",
                "SpatialIndexResult: SKIPPED (GeometryResult missing)");
        }
    }

    // ================================================================
    // Step 6: DsspResult (no dependencies)
    // ================================================================
    {
        auto dssp = DsspResult::Compute(conf);
        if (dssp && conf.AttachResult(std::move(dssp))) {
            results_attached++;
            OperationLog::Info("FullPipelineTest",
                "DsspResult: attached");
        } else {
            results_failed++;
            OperationLog::Error("FullPipelineTest",
                "DsspResult: FAILED to attach");
        }
    }

    // ================================================================
    // Step 7: ApbsFieldResult (depends on ChargeAssignmentResult)
    // ================================================================
    {
        if (conf.HasResult<ChargeAssignmentResult>()) {
            auto apbs = ApbsFieldResult::Compute(conf);
            if (apbs && conf.AttachResult(std::move(apbs))) {
                results_attached++;
                OperationLog::Info("FullPipelineTest",
                    "ApbsFieldResult: attached");
            } else {
                results_failed++;
                OperationLog::Error("FullPipelineTest",
                    "ApbsFieldResult: FAILED to attach");
            }
        } else {
            results_failed++;
            OperationLog::Error("FullPipelineTest",
                "ApbsFieldResult: SKIPPED (ChargeAssignmentResult missing)");
        }
    }

    // ================================================================
    // Step 8: MolecularGraphResult (depends on SpatialIndexResult)
    // ================================================================
    {
        if (conf.HasResult<SpatialIndexResult>()) {
            auto graph = MolecularGraphResult::Compute(conf);
            if (graph && conf.AttachResult(std::move(graph))) {
                results_attached++;
                OperationLog::Info("FullPipelineTest",
                    "MolecularGraphResult: attached");
            } else {
                results_failed++;
                OperationLog::Error("FullPipelineTest",
                    "MolecularGraphResult: FAILED to attach");
            }
        } else {
            results_failed++;
            OperationLog::Error("FullPipelineTest",
                "MolecularGraphResult: SKIPPED (SpatialIndexResult missing)");
        }
    }

    // ================================================================
    // Step 9: MopacResult on full protein (no result dependency)
    // ================================================================
    {
        auto mopac = MopacResult::Compute(conf, 0);
        if (mopac && conf.AttachResult(std::move(mopac))) {
            results_attached++;
            OperationLog::Info("FullPipelineTest",
                "MopacResult: attached (" +
                std::to_string(conf.AtomCount()) + " atoms)");
        } else {
            results_failed++;
            OperationLog::Error("FullPipelineTest",
                "MopacResult: FAILED to attach");
        }
    }

    // ================================================================
    // VERIFICATION PHASE
    // ================================================================

    fprintf(stderr,
        "\n"
        "================================================================\n"
        "  FULL PIPELINE VERIFICATION: 1UBQ\n"
        "================================================================\n");

    // --- Check 1: All 9 results are attached ---
    fprintf(stderr, "\n--- Check 1: Result attachment ---\n");
    fprintf(stderr, "  Results attached: %d / 9\n", results_attached);
    fprintf(stderr, "  Results failed:   %d\n", results_failed);

    EXPECT_TRUE(conf.HasResult<ProtonationDetectionResult>())
        << "ProtonationDetectionResult not attached";
    EXPECT_TRUE(conf.HasResult<GeometryResult>())
        << "GeometryResult not attached";
    EXPECT_TRUE(conf.HasResult<EnrichmentResult>())
        << "EnrichmentResult not attached";
    EXPECT_TRUE(conf.HasResult<ChargeAssignmentResult>())
        << "ChargeAssignmentResult not attached";
    EXPECT_TRUE(conf.HasResult<SpatialIndexResult>())
        << "SpatialIndexResult not attached";
    EXPECT_TRUE(conf.HasResult<DsspResult>())
        << "DsspResult not attached";
    EXPECT_TRUE(conf.HasResult<ApbsFieldResult>())
        << "ApbsFieldResult not attached";
    EXPECT_TRUE(conf.HasResult<MolecularGraphResult>())
        << "MolecularGraphResult not attached";
    EXPECT_TRUE(conf.HasResult<MopacResult>())
        << "MopacResult not attached";

    // --- Check 2: Every atom has role != Unknown (from EnrichmentResult) ---
    fprintf(stderr, "\n--- Check 2: Atom roles ---\n");
    if (conf.HasResult<EnrichmentResult>()) {
        int unknown_count = 0;
        for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
            if (conf.AtomAt(ai).role == AtomRole::Unknown) unknown_count++;
        }
        fprintf(stderr, "  Unknown roles: %d / %zu atoms\n",
            unknown_count, conf.AtomCount());
        // Allow very few unknowns (possible edge cases at termini)
        EXPECT_LT(unknown_count, 5)
            << unknown_count << " atoms have Unknown role";
    }

    // --- Check 3: Charges assigned ---
    fprintf(stderr, "\n--- Check 3: Charges ---\n");
    if (conf.HasResult<ChargeAssignmentResult>()) {
        auto& ca_result = conf.Result<ChargeAssignmentResult>();
        double sum_charge = ca_result.TotalCharge();
        fprintf(stderr, "  Assigned: %zu\n", ca_result.AssignedCount());
        fprintf(stderr, "  Unassigned: %zu\n", ca_result.UnassignedCount());
        fprintf(stderr, "  Total charge: %.2f\n", sum_charge);

        EXPECT_GT(ca_result.AssignedCount(),
                  conf.AtomCount() * 9 / 10)
            << "Too many unassigned atoms";
    }

    // --- Check 4: Every atom has spatial neighbours ---
    fprintf(stderr, "\n--- Check 4: Spatial neighbours ---\n");
    if (conf.HasResult<SpatialIndexResult>()) {
        int no_neighbours = 0;
        double sum_neighbours = 0.0;
        for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
            size_t nc = conf.AtomAt(ai).spatial_neighbours.size();
            if (nc == 0) no_neighbours++;
            sum_neighbours += static_cast<double>(nc);
        }
        double mean_neighbours = sum_neighbours /
            static_cast<double>(conf.AtomCount());
        fprintf(stderr, "  Atoms with no neighbours: %d\n", no_neighbours);
        fprintf(stderr, "  Mean neighbour count: %.1f\n", mean_neighbours);
        EXPECT_EQ(no_neighbours, 0)
            << no_neighbours << " atoms have no spatial neighbours";
    }

    // --- Check 5: At least one residue has SS assignment ---
    fprintf(stderr, "\n--- Check 5: Secondary structure ---\n");
    if (conf.HasResult<DsspResult>()) {
        const auto& dssp = conf.Result<DsspResult>();
        int n_helix = 0, n_sheet = 0, n_coil = 0, n_other = 0;
        for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
            char ss = dssp.SecondaryStructure(ri);
            if (ss == 'H' || ss == 'G' || ss == 'I') n_helix++;
            else if (ss == 'E' || ss == 'B') n_sheet++;
            else if (ss == 'C' || ss == ' ') n_coil++;
            else n_other++;
        }
        int total_res = static_cast<int>(protein->ResidueCount());
        double pct_helix = 100.0 * n_helix / total_res;
        double pct_sheet = 100.0 * n_sheet / total_res;
        double pct_coil  = 100.0 * (n_coil + n_other) / total_res;

        fprintf(stderr, "  Helix:  %d (%.1f%%)\n", n_helix, pct_helix);
        fprintf(stderr, "  Sheet:  %d (%.1f%%)\n", n_sheet, pct_sheet);
        fprintf(stderr, "  Coil:   %d (%.1f%%)\n", n_coil + n_other, pct_coil);

        EXPECT_GT(n_helix + n_sheet, 0)
            << "No residues have helix or sheet assignment";
    }

    // --- Check 6: APBS E-field is non-zero for most atoms ---
    fprintf(stderr, "\n--- Check 6: APBS E-field ---\n");
    if (conf.HasResult<ApbsFieldResult>()) {
        const auto& apbs = conf.Result<ApbsFieldResult>();
        int nonzero_E = 0;
        double sum_E_mag = 0.0;
        for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
            double E_mag = apbs.ElectricFieldAt(ai).norm();
            if (E_mag > 1e-10) nonzero_E++;
            sum_E_mag += E_mag;
        }
        double mean_E = sum_E_mag / static_cast<double>(conf.AtomCount());
        fprintf(stderr, "  Non-zero E-field: %d / %zu\n",
            nonzero_E, conf.AtomCount());
        fprintf(stderr, "  Mean |E|: %.6e\n", mean_E);

        EXPECT_GT(nonzero_E, static_cast<int>(conf.AtomCount()) / 2)
            << "Too few atoms with non-zero E-field";
    }

    // --- Check 7: Graph distances are assigned ---
    fprintf(stderr, "\n--- Check 7: Graph distances ---\n");
    if (conf.HasResult<MolecularGraphResult>()) {
        int has_ring_dist = 0;
        int has_N_dist = 0;
        int has_O_dist = 0;
        for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
            if (conf.AtomAt(ai).graph_dist_ring >= 0) has_ring_dist++;
            if (conf.AtomAt(ai).graph_dist_N >= 0) has_N_dist++;
            if (conf.AtomAt(ai).graph_dist_O >= 0) has_O_dist++;
        }
        fprintf(stderr, "  ring dist assigned: %d / %zu\n",
            has_ring_dist, conf.AtomCount());
        fprintf(stderr, "  N dist assigned:    %d / %zu\n",
            has_N_dist, conf.AtomCount());
        fprintf(stderr, "  O dist assigned:    %d / %zu\n",
            has_O_dist, conf.AtomCount());

        EXPECT_GT(has_ring_dist, 0)
            << "No atoms have graph_dist_ring assigned";
    }

    // --- SUMMARY ---
    fprintf(stderr,
        "\n"
        "================================================================\n"
        "  PIPELINE SUMMARY\n"
        "================================================================\n"
        "  Atom count:     %zu\n"
        "  Residue count:  %zu\n"
        "  Ring count:     %zu\n"
        "  Bond count:     %zu\n"
        "  Results:        %d attached, %d failed\n"
        "================================================================\n\n",
        conf.AtomCount(),
        protein->ResidueCount(),
        protein->RingCount(),
        protein->BondCount(),
        results_attached,
        results_failed);
}

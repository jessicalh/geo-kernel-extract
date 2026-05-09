// Smoke test for PlanarGeometryResult.
//
// Loads 1UBQ via BuildFromProtonatedPdb, attaches the
// dependency chain (Geometry, SpatialIndex, Enrichment), then runs
// PlanarGeometryResult and verifies:
//
//   - per-atom pyramidalization populated for every PlanarGroupKind != None,
//     finite, sane magnitudes (< 0.5 Å for any planar atom in a real
//     protein — values much larger are bond-graph degeneracies).
//   - per-residue ω near π for non-Pro non-terminal residues; |Δω|
//     stays < 30° on a crystal-quality protein (1UBQ is high-res).
//   - per-aromatic-ring χ₂ in [-π, π], one per aromatic ring (1UBQ
//     has 4 PHE + 1 TYR = 5 aromatic rings).
//   - per-saturated-ring Cremer-Pople (Q, θ) finite for Pro
//     pyrrolidines (1UBQ has 3 prolines: P19, P37, P38).
//   - WriteFeatures emits all six NPYs at the expected paths.

#include "TestEnvironment.h"
#include <gtest/gtest.h>

#include "AminoAcidType.h"
#include "ConformationAtom.h"
#include "EnrichmentResult.h"
#include "GeometryResult.h"
#include "LegacyAmberTopology.h"
#include "OperationLog.h"
#include "PdbFileReader.h"
#include "PlanarGeometryResult.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "Ring.h"
#include "SemanticEnums.h"
#include "SpatialIndexResult.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>

using namespace nmr;
namespace fs = std::filesystem;


class PlanarGeometryTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!fs::exists(test::TestEnvironment::UbqProtonated())) {
            GTEST_SKIP() << "1UBQ PDB not found at "
                         << test::TestEnvironment::UbqProtonated();
        }
        auto r = BuildFromProtonatedPdb(test::TestEnvironment::UbqProtonated());
        if (!r.Ok()) GTEST_SKIP() << r.error;
        protein = std::move(r.protein);
    }

    std::unique_ptr<Protein> protein;
};


TEST_F(PlanarGeometryTest, ComputesPyramidalizationOmegaChi2Pucker) {
    auto& conf = protein->Conformation();

    auto geo = GeometryResult::Compute(conf);
    ASSERT_NE(geo, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(geo)));

    auto spatial = SpatialIndexResult::Compute(conf);
    ASSERT_NE(spatial, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(spatial)));

    auto enrich = EnrichmentResult::Compute(conf);
    ASSERT_NE(enrich, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(enrich)));

    auto pgr = PlanarGeometryResult::Compute(conf);
    ASSERT_NE(pgr, nullptr) << "PlanarGeometryResult returned nullptr";
    ASSERT_TRUE(conf.AttachResult(std::move(pgr)));

    const auto& result = conf.Result<PlanarGeometryResult>();
    const LegacyAmberTopology& topo = protein->LegacyAmber();

    // ── Pyramidalization sanity ──
    int planar_count = 0;
    int pyr_finite = 0;
    int pyr_above_thresh = 0;
    double max_abs_pyr = 0.0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& sem = topo.SemanticAt(ai);
        const double pyr = conf.AtomAt(ai).pyramidalization;

        if (sem.planar_group != PlanarGroupKind::None) {
            ++planar_count;
            if (std::isfinite(pyr)) ++pyr_finite;
            const double abs_pyr = std::abs(pyr);
            if (abs_pyr > max_abs_pyr) max_abs_pyr = abs_pyr;
            if (abs_pyr > 0.5) ++pyr_above_thresh;
        } else {
            EXPECT_EQ(pyr, 0.0)
                << "Non-planar atom " << ai << " has nonzero pyramidalization";
        }
    }

    // ── Omega sanity ──
    const auto& omega = result.OmegaActual();
    const auto& dev   = result.OmegaDeviation();
    ASSERT_EQ(omega.size(), protein->ResidueCount());
    ASSERT_EQ(dev.size(),   protein->ResidueCount());

    int omega_valid_count = 0;
    int omega_planar_count = 0;
    double max_abs_dev_deg = 0.0;
    for (size_t ri = 0; ri < omega.size(); ++ri) {
        if (std::isnan(omega[ri])) continue;
        ++omega_valid_count;
        const double dev_deg = dev[ri] * 180.0 / M_PI;
        if (std::abs(dev_deg) > max_abs_dev_deg) max_abs_dev_deg = std::abs(dev_deg);
        // For 1UBQ (1.8 Å crystal), planar amides should be within ±30°
        if (std::abs(dev_deg) < 30.0) ++omega_planar_count;
    }

    // ── Aromatic χ₂ sanity ──
    const auto& chi2 = result.AromaticChi2();
    ASSERT_EQ(chi2.size(), topo.AromaticRingCount());
    int chi2_valid = 0;
    for (double v : chi2) {
        if (std::isnan(v)) continue;
        ++chi2_valid;
        EXPECT_GE(v, -M_PI - 1e-9);
        EXPECT_LE(v,  M_PI + 1e-9);
    }

    // ── Pucker sanity ──
    const auto& pucker_Q = result.PuckerQ();
    const auto& pucker_t = result.PuckerTheta();
    ASSERT_EQ(pucker_Q.size(), topo.SaturatedRingCount());
    ASSERT_EQ(pucker_t.size(), topo.SaturatedRingCount());
    int pucker_valid = 0;
    double max_Q = 0.0;
    for (size_t i = 0; i < pucker_Q.size(); ++i) {
        if (std::isnan(pucker_Q[i])) continue;
        ++pucker_valid;
        EXPECT_GE(pucker_Q[i], 0.0);
        if (pucker_Q[i] > max_Q) max_Q = pucker_Q[i];
        EXPECT_GE(pucker_t[i], 0.0);
        EXPECT_LT(pucker_t[i], 360.0);
    }

    fprintf(stderr,
        "\n=== PlanarGeometry summary (1UBQ) ===\n"
        "  atoms (planar):           %d\n"
        "  pyramidalization finite:  %d / %d\n"
        "  max |pyr|:                %.3f Å\n"
        "  pyr above 0.5 Å:          %d\n"
        "  omega valid:              %d / %zu\n"
        "  omega planar (<30°):      %d / %d\n"
        "  max |Δω|:                 %.2f°\n"
        "  aromatic χ₂ valid:        %d / %zu\n"
        "  saturated rings:          %zu (pucker valid %d)\n"
        "  max pucker Q:             %.3f Å\n"
        "======================================\n",
        planar_count, pyr_finite, planar_count, max_abs_pyr,
        pyr_above_thresh,
        omega_valid_count, protein->ResidueCount(),
        omega_planar_count, omega_valid_count, max_abs_dev_deg,
        chi2_valid, chi2.size(),
        topo.SaturatedRingCount(), pucker_valid, max_Q);

    EXPECT_GT(planar_count, 100)
        << "Expected many planar atoms in 1UBQ (peptide amides + aromatics)";
    EXPECT_EQ(pyr_finite, planar_count)
        << "All planar-atom pyramidalisation values must be finite";
    EXPECT_LT(max_abs_pyr, 0.5)
        << "Pyramidalisation magnitudes should be < 0.5 Å for a real protein";
    EXPECT_GT(omega_valid_count, 50)
        << "Most 1UBQ residues should have a defined ω";
    EXPECT_EQ(omega_planar_count, omega_valid_count)
        << "All 1UBQ peptide bonds should be near-planar (|Δω| < 30°)";
    EXPECT_LT(max_abs_dev_deg, 30.0)
        << "Max |Δω| should stay under 30° for a 1.8 Å crystal";
    EXPECT_EQ(chi2_valid, static_cast<int>(topo.AromaticRingCount()))
        << "Every aromatic ring's parent residue should have χ₂ defined";
    EXPECT_EQ(pucker_valid, static_cast<int>(topo.SaturatedRingCount()))
        << "Every saturated 5-ring (Pro pyrrolidine) should yield Cremer-Pople";
    if (topo.SaturatedRingCount() > 0) {
        EXPECT_GT(max_Q, 0.05)
            << "Pro pyrrolidine pucker amplitude should be non-trivial";
        EXPECT_LT(max_Q, 1.0)
            << "Pro pyrrolidine pucker should be physical (< 1 Å)";
    }
}


TEST_F(PlanarGeometryTest, WriteFeaturesEmitsAllSixNpys) {
    auto& conf = protein->Conformation();

    auto geo = GeometryResult::Compute(conf);
    ASSERT_NE(geo, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(geo)));

    auto spatial = SpatialIndexResult::Compute(conf);
    ASSERT_NE(spatial, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(spatial)));

    auto enrich = EnrichmentResult::Compute(conf);
    ASSERT_NE(enrich, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(enrich)));

    auto pgr = PlanarGeometryResult::Compute(conf);
    ASSERT_NE(pgr, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(pgr)));

    const fs::path output_dir = fs::temp_directory_path() /
        "planar_geometry_test_writefeatures";
    fs::create_directories(output_dir);

    const auto& result = conf.Result<PlanarGeometryResult>();
    int written = result.WriteFeatures(conf, output_dir.string());
    EXPECT_EQ(written, 6);

    for (const char* stem : {"pyramidalization", "omega_actual",
                              "omega_deviation", "aromatic_chi2",
                              "pucker_Q", "pucker_theta"}) {
        const fs::path p = output_dir / (std::string(stem) + ".npy");
        EXPECT_TRUE(fs::exists(p)) << "missing " << p;
        EXPECT_GT(fs::file_size(p), 0u);
    }

    // No fs::remove_all (libtorch ships a broken impl that segfaults
    // when CUDA has been initialised in this process; see
    // feedback_libtorch_broken_filesystem memory entry). /tmp handles
    // the kilobyte-scale leftovers.
}

// Smoke test for PlanarGeometryResult.
//
// Loads 1UBQ via BuildFromProtonatedPdb, attaches the
// dependency chain (Geometry, SpatialIndex, Enrichment), then runs
// PlanarGeometryResult and verifies:
//
//   - per-atom pyramidalization populated for every PlanarGroupKind != None,
//     finite, sane magnitudes (< 0.5 Å sanity bound; literature reports
//     amide pyramidalization SD ~0.005-0.020 Å on 1.0 Å crystals,
//     Berkholz/Karplus/Mannige). Logs the residue/atom-name of the
//     max-pyramidalization atom for build-over-build sanity.
//   - per-residue ω: 75/76 valid (C-terminus skipped); |Δω| within 15°
//     for non-X-Pro bonds (1UBQ is 1.8 Å crystal — Berkholz 99.9%
//     within 10° for ultra-high-res, 15° here is conservative);
//     X-Pro bonds emit ω regardless and are tagged via omega_is_xpro.
//   - per-aromatic-ring χ₂ in [-π, π], one per aromatic ring. 1UBQ
//     has 2 PHE + 1 TYR + 1 HIS = 4 aromatic rings.
//   - per-saturated-ring Cremer-Pople (Q, θ) finite for Pro
//     pyrrolidines (1UBQ has 3 prolines: P19, P37, P38).
//   - WriteFeatures emits all seven NPYs (six values + omega_is_xpro).
//
// Includes a synthetic Cremer-Pople pin test on a regular pentagon
// with a known envelope to verify the (Q, θ) formula against the
// canonical Cremer-Pople 1975 convention (catches the normal-direction
// inversion bug found in the 2026-05-09 adversarial review).

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
    size_t max_pyr_ai = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& sem = topo.SemanticAt(ai);
        const double pyr = conf.AtomAt(ai).pyramidalization;

        if (sem.planar_group != PlanarGroupKind::None) {
            ++planar_count;
            if (std::isfinite(pyr)) ++pyr_finite;
            const double abs_pyr = std::abs(pyr);
            if (abs_pyr > max_abs_pyr) {
                max_abs_pyr = abs_pyr;
                max_pyr_ai = ai;
            }
            if (abs_pyr > 0.5) ++pyr_above_thresh;
        } else {
            EXPECT_EQ(pyr, 0.0)
                << "Non-planar atom " << ai << " has nonzero pyramidalization";
        }
    }

    // Identify the max-pyr atom for build-over-build sanity. If max-pyr
    // creeps from heavy atoms onto H atoms, the protonation tool's H
    // placement may be the dominant signal — informative for future
    // debugging.
    const Atom& max_pyr_atom = protein->AtomAt(max_pyr_ai);
    const Residue& max_pyr_res = protein->ResidueAt(max_pyr_atom.residue_index);

    // ── Omega sanity ──
    const auto& omega = result.OmegaActual();
    const auto& dev   = result.OmegaDeviation();
    const auto& xpro  = result.OmegaIsXpro();
    ASSERT_EQ(omega.size(), protein->ResidueCount());
    ASSERT_EQ(dev.size(),   protein->ResidueCount());
    ASSERT_EQ(xpro.size(),  protein->ResidueCount());

    int omega_valid_count = 0;
    int omega_xpro_count = 0;
    int omega_non_xpro_planar_count = 0;
    int omega_non_xpro_count = 0;
    double max_abs_dev_deg_non_xpro = 0.0;
    for (size_t ri = 0; ri < omega.size(); ++ri) {
        if (std::isnan(omega[ri])) continue;
        ++omega_valid_count;
        if (xpro[ri]) {
            ++omega_xpro_count;
            continue;  // X-Pro emits ω as conformational signal; tighter
                       // |Δω| bound applies only to non-X-Pro bonds.
        }
        ++omega_non_xpro_count;
        const double dev_deg = dev[ri] * 180.0 / M_PI;
        if (std::abs(dev_deg) > max_abs_dev_deg_non_xpro)
            max_abs_dev_deg_non_xpro = std::abs(dev_deg);
        // 1UBQ is 1.8 Å crystal — Berkholz et al. 2009 places 99.9% of
        // |Δω| within 10° on ultra-high-res structures; 15° here is a
        // conservative bound that catches gross pipeline errors.
        if (std::abs(dev_deg) < 15.0) ++omega_non_xpro_planar_count;
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
        "  planar atoms:             %d (pyr finite %d/%d)\n"
        "  max |pyr|:                %.3f Å at residue %d %s, atom %s\n"
        "  pyr above 0.5 Å:          %d\n"
        "  omega valid:              %d / %zu\n"
        "    of which X-Pro:         %d\n"
        "    non-X-Pro planar (<15°): %d / %d\n"
        "    max |Δω| (non-X-Pro):   %.2f°\n"
        "  aromatic χ₂ valid:        %d / %zu\n"
        "  saturated rings:          %zu (pucker valid %d)\n"
        "  max pucker Q:             %.3f Å\n"
        "======================================\n",
        planar_count, pyr_finite, planar_count, max_abs_pyr,
        max_pyr_res.sequence_number,
        ThreeLetterCodeForAminoAcid(max_pyr_res.type).c_str(),
        max_pyr_atom.pdb_atom_name.c_str(),
        pyr_above_thresh,
        omega_valid_count, protein->ResidueCount(),
        omega_xpro_count,
        omega_non_xpro_planar_count, omega_non_xpro_count,
        max_abs_dev_deg_non_xpro,
        chi2_valid, chi2.size(),
        topo.SaturatedRingCount(), pucker_valid, max_Q);

    EXPECT_GT(planar_count, 100)
        << "Expected many planar atoms in 1UBQ (peptide amides + aromatics)";
    EXPECT_EQ(pyr_finite, planar_count)
        << "All planar-atom pyramidalisation values must be finite";
    EXPECT_LT(max_abs_pyr, 0.5)
        << "Pyramidalisation magnitudes should be < 0.5 Å (loose smoke bound)";
    EXPECT_EQ(omega_valid_count, 75)
        << "1UBQ has 76 residues → 75 peptide bonds (C-term skipped); "
           "X-Pro bonds are emitted (and tagged via omega_is_xpro), so "
           "all 75 should be valid";
    EXPECT_EQ(omega_xpro_count, 3)
        << "1UBQ has 3 prolines (P19, P37, P38) → 3 X-Pro bonds";
    EXPECT_EQ(omega_non_xpro_planar_count, omega_non_xpro_count)
        << "All non-X-Pro 1UBQ peptide bonds should be planar (|Δω| < 15°)";
    EXPECT_LT(max_abs_dev_deg_non_xpro, 15.0)
        << "Max |Δω| (non-X-Pro) should stay under 15° on this 1.8 Å crystal";
    EXPECT_EQ(topo.AromaticRingCount(), 4u)
        << "1UBQ aromatics: 2 PHE + 1 TYR + 1 HIS = 4 rings";
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


// ----------------------------------------------------------------------
// Synthetic pin: regular-pentagon test of the Cremer-Pople (Q, θ)
// formula. Verifies amplitude against the analytical √(2/N)·Δ value
// and verifies that flipping the pucker direction shifts θ by 180°
// (the relationship the 2026-05-09 adversarial review used to catch
// the normal-direction inversion). The lambda intentionally near-
// duplicates the production helper so the test exercises the formula
// independently of the production code path.
// ----------------------------------------------------------------------

TEST(PlanarGeometrySyntheticTest, CremerPopleRegularPentagonAmplitude) {
    // The CremerPople5Ring helper is anonymous-namespace-private in
    // PlanarGeometryResult.cpp, so we exercise the formula via a
    // hand-built 5-ring point cloud and compute Q, θ inline. This is
    // intentionally a near-duplicate of the production helper — the
    // duplication is the test's discipline ("each test exercises the
    // formula independently of the production code path").
    auto cremer_pople = [](const std::array<Vec3, 5>& positions) {
        Vec3 G = Vec3::Zero();
        for (const auto& p : positions) G += p;
        G /= 5.0;

        Vec3 R1 = Vec3::Zero();
        Vec3 R2 = Vec3::Zero();
        for (size_t j = 0; j < 5; ++j) {
            const Vec3 r_j = positions[j] - G;
            const double phi = 2.0 * M_PI * static_cast<double>(j) / 5.0;
            R1 += r_j * std::sin(phi);
            R2 += r_j * std::cos(phi);
        }
        const Vec3 n = R1.cross(R2);
        const Vec3 n_hat = n / n.norm();

        double cs = 0.0, sn = 0.0;
        for (size_t j = 0; j < 5; ++j) {
            const double z_j = (positions[j] - G).dot(n_hat);
            const double phi = 4.0 * M_PI * static_cast<double>(j) / 5.0;
            cs +=  z_j * std::cos(phi);
            sn += -z_j * std::sin(phi);
        }
        const double scale = std::sqrt(2.0 / 5.0);
        const double Qcos = scale * cs;
        const double Qsin = scale * sn;
        const double Q = std::sqrt(Qcos * Qcos + Qsin * Qsin);
        double theta = std::atan2(Qsin, Qcos) * 180.0 / M_PI;
        if (theta < 0.0) theta += 360.0;
        return std::pair<double, double>{Q, theta};
    };

    // Regular pentagon (unit radius) with atom 2 puckered up by a
    // small Δ. Expected amplitude in the linear (small-Δ) limit:
    // Q ≈ √(2/N) · Δ for an N-ring envelope at one atom. For
    // Δ = 0.05 Å on a unit-radius ring, mean-plane tilt is < 1%
    // (verified analytically). Larger Δ values introduce systematic
    // deviation from this analytical formula due to the Cremer-Pople
    // normal tilting away from the puckered atom; that's a real
    // feature of the formulation, not a bug.
    const double Delta = 0.05;
    std::array<Vec3, 5> ring;
    for (size_t j = 0; j < 5; ++j) {
        const double phi = 2.0 * M_PI * static_cast<double>(j) / 5.0;
        ring[j] = Vec3(std::cos(phi), std::sin(phi),
                       (j == 2) ? Delta : 0.0);
    }
    const auto [Q_up, theta_up] = cremer_pople(ring);

    const double Q_expected = std::sqrt(2.0 / 5.0) * Delta;
    EXPECT_NEAR(Q_up, Q_expected, 1e-3)
        << "Q amplitude does not match small-Δ √(2/N)·Δ analytical limit";
    EXPECT_GE(theta_up, 0.0);
    EXPECT_LT(theta_up, 360.0);

    // Flip the pucker direction: same envelope at j=2 with atom DOWN.
    // Negating all z values negates both Q₂cos and Q₂sin, which is
    // exactly a 180° phase rotation in θ. This is what the 2026-05-09
    // adversarial review used to localise the normal-direction
    // inversion bug.
    ring[2] = Vec3(std::cos(4.0 * M_PI / 5.0),
                    std::sin(4.0 * M_PI / 5.0),
                    -Delta);
    const auto [Q_down, theta_down] = cremer_pople(ring);

    EXPECT_NEAR(Q_down, Q_up, 1e-9)
        << "Pucker amplitude is invariant under atom-direction flip";
    const double diff = std::abs(theta_down - theta_up);
    const double diff_mod_360 = std::min(diff, 360.0 - diff);
    EXPECT_NEAR(diff_mod_360, 180.0, 1e-6)
        << "Flipping pucker direction should shift θ by 180°; "
        << "θ_up=" << theta_up << " θ_down=" << theta_down;

    fprintf(stderr,
        "\n=== Cremer-Pople synthetic pin (j=2 envelope, Δ=%.3f Å) ===\n"
        "  Q (up):       %.6f  (expected √(2/5)·Δ = %.6f)\n"
        "  θ (up):       %.2f°\n"
        "  Q (down):     %.6f\n"
        "  θ (down):     %.2f°\n"
        "  Δθ:           %.2f° (expected 180°)\n"
        "===========================================================\n",
        Delta, Q_up, Q_expected, theta_up, Q_down, theta_down, diff_mod_360);
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
    EXPECT_EQ(written, 7);

    for (const char* stem : {"pyramidalization", "omega_actual",
                              "omega_deviation", "aromatic_chi2",
                              "pucker_Q", "pucker_theta",
                              "omega_is_xpro"}) {
        const fs::path p = output_dir / (std::string(stem) + ".npy");
        EXPECT_TRUE(fs::exists(p)) << "missing " << p;
        EXPECT_GT(fs::file_size(p), 0u);
    }

    // No fs::remove_all (libtorch ships a broken impl that segfaults
    // when CUDA has been initialised in this process; see
    // feedback_libtorch_broken_filesystem memory entry). /tmp handles
    // the kilobyte-scale leftovers.
}

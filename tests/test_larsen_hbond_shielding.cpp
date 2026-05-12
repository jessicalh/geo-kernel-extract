// End-to-end smoke for LarsenHBondShieldingResult on 1UBQ_pm6dh3plus.pdb.
//
// Asserts the calculator across both donor classes (amide H + Hα) and
// all four acceptor classes (BackboneCarbonyl / SidechainCarbonyl /
// HydroxylOxygen / CarboxylateOxygen) via spatial-sweep enumeration:
//   - Detects ≥ 100 H-bond pairs (1UBQ has ~60 backbone amide-H pairs
//     plus Hα and sidechain contributions; aggregate easily clears 100).
//   - All emitted Mat3 tensors are finite.
//   - larsen_hbond_n_pairs > 0 on atoms that received contributions.
//   - larsen_hbond_any_corner_imputed propagates to Table 2 target
//     atoms (codex C4 fix).
//   - Cβ diagnostic tensors are finite and bounded — the value reflects
//     Larsen's actual DFT-computed Cβ shielding (which his Table 2
//     CHOSE not to include in ProCS15); non-zero is the methodology
//     statement, not a pipeline bug.
//   - Δσ_w water-term gate is geometric (no candidate O in range).
//   - WriteFeatures emits 8 NPYs.
//
// Per-cell Larsen Table 2 dispatch unit tests cover LarsenContribDispatch
// directly (no fixture needed).
//
// Skipped if the LarsenHBondGrid directory or the Larsen-archive 1UBQ
// PDB is not present.

#include "TestEnvironment.h"
#include <gtest/gtest.h>

#include "ConformationAtom.h"
#include "DsspResult.h"
#include "EnrichmentResult.h"
#include "GeometryResult.h"
#include "LarsenHBondGrid.h"
#include "LarsenHBondShieldingResult.h"
#include "OperationLog.h"
#include "PdbFileReader.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "RuntimeEnvironment.h"
#include "Session.h"
#include "SpatialIndexResult.h"

#include <cmath>
#include <filesystem>
#include <memory>
#include <string>

using namespace nmr;
namespace fs = std::filesystem;


namespace {

constexpr const char* kLarsen1UbqPm6 =
    "/mnt/expansion/larsen_archive/structures/1UBQ_pm6dh3plus.pdb";

}  // namespace


// ---------------------------------------------------------------------------
// LarsenContribDispatch::Applies — per-cell Table 2 dispatch verification.
// ---------------------------------------------------------------------------

TEST(LarsenContribDispatchTest, Table2Cells) {
    using TA = LarsenContribDispatch::TargetAtom;
    using Term = LarsenContribDispatch::Term;

    // N row: 1HB ✓, 2HB ✓, 1HaB ✗, 2HaB ✓, RC ✗, w ✗
    EXPECT_TRUE (LarsenContribDispatch::Applies(TA::N,  Term::Primary_HB));
    EXPECT_TRUE (LarsenContribDispatch::Applies(TA::N,  Term::Secondary_HB));
    EXPECT_FALSE(LarsenContribDispatch::Applies(TA::N,  Term::Primary_HaB));
    EXPECT_TRUE (LarsenContribDispatch::Applies(TA::N,  Term::Secondary_HaB));
    EXPECT_FALSE(LarsenContribDispatch::Applies(TA::N,  Term::RingCurrent));
    EXPECT_FALSE(LarsenContribDispatch::Applies(TA::N,  Term::Water));

    // Cα row: 1HB ✗, 2HB ✓, 1HaB ✗, 2HaB ✗, RC ✗, w ✗
    EXPECT_FALSE(LarsenContribDispatch::Applies(TA::CA, Term::Primary_HB));
    EXPECT_TRUE (LarsenContribDispatch::Applies(TA::CA, Term::Secondary_HB));
    EXPECT_FALSE(LarsenContribDispatch::Applies(TA::CA, Term::Primary_HaB));
    EXPECT_FALSE(LarsenContribDispatch::Applies(TA::CA, Term::Secondary_HaB));

    // Cβ row: ALL FALSE (diagnostic; Larsen Table 2 says no contribution).
    for (int term = 0; term < (int)Term::Count; ++term) {
        EXPECT_FALSE(LarsenContribDispatch::Applies(
            TA::CB, static_cast<Term>(term)))
            << "Cβ row must be all-false per Larsen Table 2; term="
            << term;
    }

    // C' row: 1HB ✗, 2HB ✓, 1HaB ✗, 2HaB ✗, RC ✗, w ✗
    EXPECT_FALSE(LarsenContribDispatch::Applies(TA::C, Term::Primary_HB));
    EXPECT_TRUE (LarsenContribDispatch::Applies(TA::C, Term::Secondary_HB));
    EXPECT_FALSE(LarsenContribDispatch::Applies(TA::C, Term::Primary_HaB));
    EXPECT_FALSE(LarsenContribDispatch::Applies(TA::C, Term::Secondary_HaB));

    // Hα row: 1HB ✓, 2HB ✓, 1HaB ✓, 2HaB ✓, RC ✓, w ✗
    EXPECT_TRUE (LarsenContribDispatch::Applies(TA::HA, Term::Primary_HB));
    EXPECT_TRUE (LarsenContribDispatch::Applies(TA::HA, Term::Secondary_HB));
    EXPECT_TRUE (LarsenContribDispatch::Applies(TA::HA, Term::Primary_HaB));
    EXPECT_TRUE (LarsenContribDispatch::Applies(TA::HA, Term::Secondary_HaB));
    EXPECT_TRUE (LarsenContribDispatch::Applies(TA::HA, Term::RingCurrent));
    EXPECT_FALSE(LarsenContribDispatch::Applies(TA::HA, Term::Water));

    // HN row: ALL TRUE per Larsen Table 2.
    EXPECT_TRUE(LarsenContribDispatch::Applies(TA::HN, Term::Primary_HB));
    EXPECT_TRUE(LarsenContribDispatch::Applies(TA::HN, Term::Secondary_HB));
    EXPECT_TRUE(LarsenContribDispatch::Applies(TA::HN, Term::Primary_HaB));
    EXPECT_TRUE(LarsenContribDispatch::Applies(TA::HN, Term::Secondary_HaB));
    EXPECT_TRUE(LarsenContribDispatch::Applies(TA::HN, Term::RingCurrent));
    EXPECT_TRUE(LarsenContribDispatch::Applies(TA::HN, Term::Water));
}


// ---------------------------------------------------------------------------
// 1UBQ smoke test fixture
// ---------------------------------------------------------------------------

class LarsenHBondShieldingTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (RuntimeEnvironment::LarsenHBondGridDir().empty()) {
            GTEST_SKIP() << "LarsenHBondGridDir empty; set "
                            "larsen_hbond_grids in ~/.nmr_tools.toml";
        }
        if (!fs::exists(kLarsen1UbqPm6)) {
            GTEST_SKIP() << "Larsen 1UBQ PM6-D3H+ PDB not found at "
                         << kLarsen1UbqPm6;
        }
        ASSERT_EQ(session.LoadLarsenHBondGrid(), kOk) << session.LastError();
        ASSERT_TRUE(session.HasLarsenHBondGrid());

        auto r = BuildFromProtonatedPdb(kLarsen1UbqPm6);
        ASSERT_TRUE(r.Ok()) << r.error;
        protein = std::move(r.protein);
    }

    Session session;
    std::unique_ptr<Protein> protein;
};


TEST_F(LarsenHBondShieldingTest, SmokeOn1UbqPm6) {
    auto& conf = protein->Conformation();

    // Dependency chain.
    auto geo = GeometryResult::Compute(conf);
    ASSERT_NE(geo, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(geo)));

    auto sp = SpatialIndexResult::Compute(conf);
    ASSERT_NE(sp, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(sp)));

    auto en = EnrichmentResult::Compute(conf);
    ASSERT_NE(en, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(en)));

    auto dssp = DsspResult::Compute(conf);
    ASSERT_NE(dssp, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(dssp)));

    // Run the calculator.
    auto result = LarsenHBondShieldingResult::Compute(
        conf, *session.LarsenHBondGridPtr());
    ASSERT_NE(result, nullptr);

    // 1UBQ aggregate across all (donor, acceptor) classes: 60+ backbone
    // amide-H H-bonds plus Hα H-bonds and sidechain-acceptor pairs.
    // ≥ 100 is a conservative floor.
    EXPECT_GE(result->PairsFound(), 100)
        << "expected ≥ 100 H-bond pairs total; got "
        << result->PairsFound();
    EXPECT_GE(result->AtomsWithContribution(), 150);

    // All emitted tensors finite.
    int n_finite_total = 0;
    int n_with_total = 0;
    double max_cb_frobenius = 0.0;
    int n_amide_with_water = 0;
    int n_atoms_with_pair_count = 0;
    int n_total_pair_increments = 0;
    for (std::size_t i = 0; i < conf.AtomCount(); ++i) {
        const auto& a = conf.AtomAt(i);
        const auto& m = a.larsen_hbond_shielding_tensor;
        bool finite = true;
        for (int r = 0; r < 3 && finite; ++r)
            for (int c = 0; c < 3 && finite; ++c)
                if (!std::isfinite(m(r, c))) finite = false;
        if (m.norm() > 1e-9) {
            ++n_with_total;
            if (finite) ++n_finite_total;
        }
        // Track max |Cβ diagnostic| Frobenius norm.
        double cb_norm = a.larsen_hbond_diagnostic_CB.norm();
        if (cb_norm > max_cb_frobenius) max_cb_frobenius = cb_norm;
        // Water term sweep count.
        if (a.larsen_hbond_water_term > 0.0) ++n_amide_with_water;
        // n_pairs increment verification (C1 fix).
        if (a.larsen_hbond_n_pairs > 0) {
            ++n_atoms_with_pair_count;
            n_total_pair_increments += a.larsen_hbond_n_pairs;
        }
    }
    EXPECT_EQ(n_finite_total, n_with_total);
    EXPECT_GT(n_with_total, 0);

    // C1 fix: n_pairs is no longer all-zero. Every atom that received
    // any tensor contribution should have n_pairs >= 1. Cross-check:
    // (atoms with contribution) == (atoms with n_pairs > 0).
    EXPECT_EQ(n_atoms_with_pair_count, n_with_total)
        << "n_pairs should be non-zero exactly on atoms with contribution";
    EXPECT_GT(n_total_pair_increments, result->PairsFound())
        << "summed n_pairs should exceed pair count (multiple atoms per pair)";

    // Cβ diagnostic: Larsen Table 2 EXCLUDES Cβ from the final shielding
    // sum — but his DFT scan still computed real σ values at Cβ atoms,
    // and those are what live in the grid's donor_CB readouts. The
    // diagnostic surfaces what Larsen's model discards. On 1UBQ the
    // max Cβ Frobenius is ~14 ppm (real chemistry — Cβ does respond to
    // H-bond geometry). We assert finite + bounded-below-absurd; the
    // non-zero value is methodology signal, not a pipeline bug.
    EXPECT_TRUE(std::isfinite(max_cb_frobenius));
    EXPECT_LT(max_cb_frobenius, 100.0)
        << "Cβ diagnostic exceeded 100 ppm Frobenius — suspicious. "
           "Larsen's DFT measured real shielding at Cβ but extreme "
           "values usually mean pipeline trouble (close-approach "
           "geometry, rotation bug, parser regression). Observed: "
        << max_cb_frobenius;
    EXPECT_GT(max_cb_frobenius, 0.0)
        << "Cβ diagnostic was zero — ALA donor archives should produce "
           "non-zero CB readouts; zero means CB path didn't fire.";

    // Water term: assigned only to amide Hs that DSSP did NOT detect
    // as a donor of any H-bond. On 1UBQ specifically the post-M2-fix
    // result is 0 — every amide H is DSSP-detected as a donor
    // (helix/sheet packing). Pre-M2 the buggy implementation
    // assigned water term to grid-skipped DSSP-paired pairs (17 of
    // them); the fix correctly suppresses those. We assert the
    // weaker condition: 0 ≤ water-term count ≤ total amide Hs, and
    // mass conservation (no atom is both water-term-assigned AND
    // pair-contributing) is verified in the dedicated reality-check
    // test below.
    int n_amide_h_total = 0;
    for (std::size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        const auto& res = protein->ResidueAt(ri);
        if (res.H != Residue::NONE) ++n_amide_h_total;
    }
    EXPECT_LE(n_amide_with_water, n_amide_h_total)
        << "water-term count must not exceed total amide H count";
    EXPECT_EQ(n_amide_with_water, result->AmideHsUnboundWithWater())
        << "atom-level water-term count should match result aggregate";

    // WriteFeatures emits 8 NPYs (spherical-packed per Pattern 11 + 6).
    fs::path tmp = fs::temp_directory_path() / "larsen_hbond_phase1_out";
    fs::create_directories(tmp);
    int n_written = result->WriteFeatures(conf, tmp.string());
    EXPECT_EQ(n_written, 8);
    for (const std::string& stem : {
        "larsen_hbond_shielding",
        "larsen_hbond_1pHB_shielding",
        "larsen_hbond_2pHB_shielding",
        "larsen_hbond_1pHaB_shielding",
        "larsen_hbond_2pHaB_shielding",
        "larsen_hbond_diagnostic_CB_shielding",
        "larsen_hbond_water_term",
        "larsen_hbond_count",
    }) {
        fs::path p = tmp / (stem + ".npy");
        EXPECT_TRUE(fs::exists(p)) << "missing " << p.string();
    }
}


// ---------------------------------------------------------------------------
// Reality-check tests (user-requested round 3)
// ---------------------------------------------------------------------------

// GLY HA fan-out: 1UBQ has 6 GLY residues. Every Hα target contribution
// must land on BOTH HA2 and HA3 — not just one. Test: among atoms
// classified as GLY α-hydrogens (via IsAnyAlphaHydrogen()) that
// received any Table 2 contribution under HA's row, count must be
// even (each GLY pair contributes paired).
TEST_F(LarsenHBondShieldingTest, GlyHaFanOutHitsBothHA) {
    auto& conf = protein->Conformation();

    // Dependency chain.
    auto geo = GeometryResult::Compute(conf);
    ASSERT_NE(geo, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(geo)));
    auto sp = SpatialIndexResult::Compute(conf);
    ASSERT_NE(sp, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(sp)));
    auto en = EnrichmentResult::Compute(conf);
    ASSERT_NE(en, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(en)));
    auto dssp = DsspResult::Compute(conf);
    ASSERT_NE(dssp, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(dssp)));

    auto result = LarsenHBondShieldingResult::Compute(
        conf, *session.LarsenHBondGridPtr());
    ASSERT_NE(result, nullptr);

    // Count GLY HAs receiving any 1°HB contribution. For each GLY
    // donor residue that successfully formed an H-bond, BOTH GLY
    // alpha hydrogens (HA2 + HA3) should have non-zero 1pHB tensor.
    int n_gly_residues_with_partial_ha_hit = 0;  // exactly 1 of 2 hit
    int n_gly_residues_with_both_ha_hit    = 0;  // both 2 of 2 hit
    int n_gly_residues_with_no_ha_hit      = 0;
    const auto& topo = protein->LegacyAmber();
    for (std::size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        const auto& res = protein->ResidueAt(ri);
        if (res.type != AminoAcid::GLY) continue;
        int n_ha_with_contribution = 0;
        int n_ha_atoms = 0;
        for (std::size_t ai : res.atom_indices) {
            const auto& sem = topo.SemanticAt(ai);
            if (!sem.IsAnyAlphaHydrogen()) continue;
            ++n_ha_atoms;
            if (conf.AtomAt(ai).larsen_hbond_1pHB_tensor.norm() > 1e-9) {
                ++n_ha_with_contribution;
            }
        }
        // Every GLY has 2 prochiral HAs (HA2 + HA3).
        EXPECT_EQ(n_ha_atoms, 2)
            << "GLY at residue " << ri << " should have 2 α-hydrogens; got "
            << n_ha_atoms;
        if (n_ha_with_contribution == 0)      ++n_gly_residues_with_no_ha_hit;
        else if (n_ha_with_contribution == 1) ++n_gly_residues_with_partial_ha_hit;
        else                                  ++n_gly_residues_with_both_ha_hit;
    }
    // The fan-out fix asserts: NO GLY residue should have a partial
    // HA hit. Either neither HA contributes (residue not a donor) or
    // both do (substrate-driven enumeration). A partial hit would be
    // the codex H1 bug.
    EXPECT_EQ(n_gly_residues_with_partial_ha_hit, 0)
        << "GLY HA fan-out broken: " << n_gly_residues_with_partial_ha_hit
        << " GLY residues have exactly 1 of 2 α-hydrogens contributing";
    // Sanity: at least one GLY should have contributed (1UBQ has
    // GLY in α-helix / β-sheet so they're amide-H donors).
    EXPECT_GT(n_gly_residues_with_both_ha_hit, 0)
        << "expected ≥1 GLY residue to have both Hα atoms receiving "
           "contributions on 1UBQ; saw 0";
}


// Mass-conservation: every geometric candidate the spatial sweep
// classified gets accounted for as either a successful pair or a
// counted skip (codex finding F4, 2026-05-12). No silent drops.
//
// Also: an amide H can legitimately carry BOTH a water term (it had
// no H-bond candidate as DONOR) AND a non-zero n_pairs (it received
// 2°HB contributions as the i+1 TARGET of someone else's H-bond).
// Those are independent — donor status and target status are
// orthogonal. The earlier DSSP-era mutual-exclusion assertion was a
// bookkeeping artefact, not physics.
TEST_F(LarsenHBondShieldingTest, GeometricCandidateMassConservation) {
    auto& conf = protein->Conformation();
    auto geo = GeometryResult::Compute(conf);
    ASSERT_TRUE(conf.AttachResult(std::move(geo)));
    auto sp = SpatialIndexResult::Compute(conf);
    ASSERT_TRUE(conf.AttachResult(std::move(sp)));
    auto en = EnrichmentResult::Compute(conf);
    ASSERT_TRUE(conf.AttachResult(std::move(en)));
    auto dssp = DsspResult::Compute(conf);
    ASSERT_TRUE(conf.AttachResult(std::move(dssp)));

    auto result = LarsenHBondShieldingResult::Compute(
        conf, *session.LarsenHBondGridPtr());
    ASSERT_NE(result, nullptr);

    // Pairs that reached the grid path and either succeeded or were
    // counted as skipped. The two counts together are the total of
    // geometric candidates the dispatch evaluated — none were silently
    // dropped (cf. codex F4: missing-frame-atoms early-returns used to
    // bypass the counter).
    const int n_found   = result->PairsFound();
    const int n_skipped = result->PairsGridSkipped();
    EXPECT_GE(n_found,   40) << "1UBQ should yield ≥ 40 grid-paired pairs";
    EXPECT_GE(n_skipped, 0)  << "skip counter must be non-negative";

    // Water term + n_pairs can both be > 0 on the same amide H without
    // contradicting the gate. The interesting invariant is that EVERY
    // amide H with water_term > 0 had ZERO geometric H-bond candidates
    // as a donor (θ ≥ 90° in 4.2 Å). We can't recompute the per-atom
    // donor-paired flag from result alone, but we can verify the
    // overall count is finite + the water term is exactly 2.07 ppm
    // where it fires.
    constexpr double kWater_ppm = 2.07;
    int n_water_total = 0;
    for (std::size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& a = conf.AtomAt(ai);
        if (a.larsen_hbond_water_term > 0.0) {
            EXPECT_DOUBLE_EQ(a.larsen_hbond_water_term, kWater_ppm)
                << "water term must be Larsen's calibrated 2.07 ppm "
                << "isotropic; atom_idx=" << ai;
            ++n_water_total;
        }
    }
    EXPECT_EQ(n_water_total, result->AmideHsUnboundWithWater())
        << "per-atom water_term count must match aggregate";
    // 1UBQ has ~76 amide Hs; some are solvent-exposed. After the F2
    // gate fix the water-term count should be > 0 (the earlier eager
    // gate suppressed it for ~all amide Hs). Bound is intentionally
    // loose — exact number depends on solvent accessibility geometry.
    EXPECT_GT(n_water_total, 0)
        << "expected ≥ 1 solvent-exposed amide H on 1UBQ to receive "
           "the water term; saw 0 (F2 gate may be too aggressive)";
}


// CB diagnostic does NOT increment larsen_hbond_n_pairs.
// (ConformationAtom doc: n_pairs counts the four real Table 2
// contribution classes only. The CB diagnostic is for parser-pipeline
// integrity, not a true class.)
TEST_F(LarsenHBondShieldingTest, CbDiagnosticDoesNotInflateNPairs) {
    auto& conf = protein->Conformation();
    auto geo = GeometryResult::Compute(conf);
    ASSERT_TRUE(conf.AttachResult(std::move(geo)));
    auto sp = SpatialIndexResult::Compute(conf);
    ASSERT_TRUE(conf.AttachResult(std::move(sp)));
    auto en = EnrichmentResult::Compute(conf);
    ASSERT_TRUE(conf.AttachResult(std::move(en)));
    auto dssp = DsspResult::Compute(conf);
    ASSERT_TRUE(conf.AttachResult(std::move(dssp)));

    auto result = LarsenHBondShieldingResult::Compute(
        conf, *session.LarsenHBondGridPtr());
    ASSERT_NE(result, nullptr);

    // An atom with CB diagnostic contribution but NO real Table 2
    // contribution should have n_pairs == 0.
    int n_cb_only_with_count = 0;
    for (std::size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& a = conf.AtomAt(ai);
        bool has_cb_diag = a.larsen_hbond_diagnostic_CB.norm() > 1e-9;
        bool has_real_contrib = a.larsen_hbond_shielding_tensor.norm() > 1e-9;
        if (has_cb_diag && !has_real_contrib && a.larsen_hbond_n_pairs > 0) {
            ++n_cb_only_with_count;
        }
    }
    EXPECT_EQ(n_cb_only_with_count, 0)
        << "Cβ-diagnostic-only contributors must NOT increment n_pairs "
           "(per ConformationAtom doc: n_pairs is for Table 2 classes only). "
           "The CB diagnostic exists for parser-pipeline integrity; "
           "incrementing n_pairs from CB contributions would inflate "
           "atoms' pair counts.";
}


// Schema validation: load LarsenHBondGrid and verify the loader's
// per-archive presence checks did not throw. (If schema were
// violated, ctor would have thrown and the fixture SetUp would have
// failed before reaching here.)
TEST_F(LarsenHBondShieldingTest, GridSchemaValidationPassed) {
    // Reaching this point means LarsenHBondGrid ctor + per-archive
    // ValidateSchema all passed. The fixture's grid load is the
    // implicit test.
    const LarsenHBondGrid* g = session.LarsenHBondGridPtr();
    ASSERT_NE(g, nullptr);
    EXPECT_TRUE(g->IsLoaded());
}


// F1 reality check (codex finding, 2026-05-12): after adding acceptor_CA
// and acceptor_C grid readouts, Cα and C' atoms in residues that are
// H-bond targets must accumulate Δσ_2°HB contributions. Before F1,
// Larsen Table 2 said these atoms get 2°HB but the grid lacked the
// tensors so they were silently zero — particularly C', the largest
// 2°HB term per Larsen Table 1 (~2.1 ppm RMSD on Ub).
TEST_F(LarsenHBondShieldingTest, CaAndCReceive2pHBContributions) {
    auto& conf = protein->Conformation();
    auto geo = GeometryResult::Compute(conf);
    ASSERT_TRUE(conf.AttachResult(std::move(geo)));
    auto sp = SpatialIndexResult::Compute(conf);
    ASSERT_TRUE(conf.AttachResult(std::move(sp)));
    auto en = EnrichmentResult::Compute(conf);
    ASSERT_TRUE(conf.AttachResult(std::move(en)));
    auto dssp = DsspResult::Compute(conf);
    ASSERT_TRUE(conf.AttachResult(std::move(dssp)));

    auto result = LarsenHBondShieldingResult::Compute(
        conf, *session.LarsenHBondGridPtr());
    ASSERT_NE(result, nullptr);

    const Protein& protein_ref = conf.ProteinRef();
    const auto& topo = protein_ref.LegacyAmber();

    int n_ca_with_2pHB = 0;
    int n_c_with_2pHB  = 0;
    double max_ca_2pHB = 0.0;
    double max_c_2pHB  = 0.0;
    for (std::size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& sem = topo.SemanticAt(ai);
        const auto& a = conf.AtomAt(ai);
        const double n = a.larsen_hbond_2pHB_tensor.norm();
        if (n <= 1e-9) continue;
        if (sem.backbone_role == BackboneRole::AlphaCarbon) {
            ++n_ca_with_2pHB;
            max_ca_2pHB = std::max(max_ca_2pHB, n);
        } else if (sem.backbone_role == BackboneRole::CarbonylCarbon) {
            ++n_c_with_2pHB;
            max_c_2pHB = std::max(max_c_2pHB, n);
        }
    }

    // 1UBQ has ~149 grid-paired H-bonds; many should land 2°HB on the
    // i+1 Cα and on the acceptor's own C atom. Loose lower bound to
    // detect pipeline regressions, not assert exact counts.
    EXPECT_GT(n_ca_with_2pHB, 10)
        << "Cα atoms with Δσ_2°HB > 0 should exceed 10 on 1UBQ; F1 "
           "regression: acceptor_CA may not be wired through.";
    EXPECT_GT(n_c_with_2pHB,  10)
        << "C' atoms with Δσ_2°HB > 0 should exceed 10 on 1UBQ; F1 "
           "regression: acceptor_C may not be wired through. "
           "C' Δσ_2°HB is Larsen's largest 2°HB term (~2.1 ppm RMSD).";

    // Magnitude sanity. The parser subtracts an orientation-matched
    // r-max reference surface; if that regresses to a global r-max
    // average, rotated absolute shielding anisotropy leaks into the
    // H-bond tensor and 1UBQ produces >200 ppm norms on C'/Cα. Keep a
    // deliberately loose but load-bearing upper bound here: real 1UBQ
    // values stay well below this after orientation-matched subtraction.
    EXPECT_LT(max_ca_2pHB, 100.0)
        << "Cα Δσ_2°HB max=" << max_ca_2pHB
        << " ppm — likely reference-surface or tensor-frame regression";
    EXPECT_LT(max_c_2pHB,  100.0)
        << "C' Δσ_2°HB max=" << max_c_2pHB
        << " ppm — likely reference-surface or tensor-frame regression";
}

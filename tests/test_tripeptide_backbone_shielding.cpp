// Smoke test for TripeptideBackboneShieldingResult.
//
// Loads 1UBQ_pm6dh3plus.pdb (Larsen's PM6-D3H+ optimised geometry —
// /mnt/expansion/larsen_archive/structures/, the published-RMSD
// validation target from Larsen 2015) via BuildFromProtonatedPdb,
// attaches the dependency chain (Geometry, SpatialIndex, Enrichment),
// then runs TripeptideBackboneShieldingResult against the local
// tensorcs15 replica and verifies:
//
//   - residues attempted == 76 (1UBQ length); matched is most of them
//     (terminal pair always skipped because phi/psi requires both
//     flanking residues; LYS/ARG/etc. with chi3/chi4 set should hit
//     the chi3+chi4 lookup path)
//   - per-atom tensors are finite where has_match is true
//   - frame_type discriminator: any SER residue should hit
//     orca_input_orientation rows (project SER PBE regen)
//   - WriteFeatures emits 3 NPYs
//
// Skipped if the tensorcs15 DSN is not configured (test machine
// doesn't have the local DB), or if the Larsen archive PDB is not at
// the expected path.

#include "TestEnvironment.h"
#include <gtest/gtest.h>

#include "AminoAcidType.h"
#include "ConformationAtom.h"
#include "EnrichmentResult.h"
#include "GeometryResult.h"
#include "OperationLog.h"
#include "PdbFileReader.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "RuntimeEnvironment.h"
#include "Session.h"
#include "SpatialIndexResult.h"
#include "TripeptideBackboneShieldingResult.h"
#include "TripeptideDftTable.h"

#include <cmath>
#include <filesystem>
#include <string>

using namespace nmr;
namespace fs = std::filesystem;


namespace {

constexpr const char* kLarsen1UbqPm6 =
    "/mnt/expansion/larsen_archive/structures/1UBQ_pm6dh3plus.pdb";

}  // namespace


class TripeptideBackboneShieldingTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (RuntimeEnvironment::TensorCs15Dsn().empty()) {
            GTEST_SKIP() << "tensorcs15 DSN not configured "
                            "(set [databases].tensorcs15 in "
                            "~/.nmr_tools.toml)";
        }
        if (!fs::exists(kLarsen1UbqPm6)) {
            GTEST_SKIP() << "Larsen 1UBQ PM6-D3H+ PDB not found at "
                         << kLarsen1UbqPm6
                         << " (download via larsen ERDA archive — see "
                            "/mnt/expansion/larsen_archive/README.md)";
        }
        ASSERT_EQ(session.LoadTripeptideDftTable(), kOk)
            << session.LastError();
        ASSERT_TRUE(session.HasTripeptideDftTable());

        auto r = BuildFromProtonatedPdb(kLarsen1UbqPm6);
        ASSERT_TRUE(r.Ok()) << r.error;
        protein = std::move(r.protein);
    }

    Session session;
    std::unique_ptr<Protein> protein;
};


TEST_F(TripeptideBackboneShieldingTest, RunsOn1UbqPm6) {
    auto& conf = protein->Conformation();

    auto geo = GeometryResult::Compute(conf);
    ASSERT_NE(geo, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(geo)));

    auto sp = SpatialIndexResult::Compute(conf);
    ASSERT_NE(sp, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(sp)));

    auto en = EnrichmentResult::Compute(conf);
    ASSERT_NE(en, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(en)));

    auto tbb = TripeptideBackboneShieldingResult::Compute(
        conf, *session.TripeptideDftTablePtr());
    ASSERT_NE(tbb, nullptr);

    // 1UBQ has 76 residues; ResiduesAttempted counts standard residues
    // (Unknown skipped). Termini fail phi/psi.
    EXPECT_EQ(tbb->ResiduesAttempted(),
              static_cast<int>(protein->ResidueCount()));
    EXPECT_GE(tbb->ResiduesMatched(), 60);  // most internals match
    EXPECT_GT(tbb->AtomsAssigned(), 0);

    // Backbone Kabsch RMSD should be small — the canonical tripeptide
    // and the protein backbone share the same N/CA/C bond geometry by
    // construction. ~0.5 Å is the empirical rough upper bound; if mean
    // exceeds 1 Å something is structurally wrong with the alignment.
    EXPECT_LT(tbb->MeanBackboneRmsd(), 1.0);

    // Spot check: every matched atom has finite tensor.
    int n_finite = 0, n_total = 0;
    for (size_t i = 0; i < conf.AtomCount(); ++i) {
        const auto& ca = conf.AtomAt(i);
        if (!ca.tripeptide_bb_has_match) continue;
        ++n_total;
        const auto& m = ca.tripeptide_bb_shielding_tensor;
        bool finite = true;
        for (int r = 0; r < 3 && finite; ++r)
            for (int c = 0; c < 3 && finite; ++c)
                if (!std::isfinite(m(r, c))) finite = false;
        if (finite) ++n_finite;
    }
    EXPECT_EQ(n_finite, n_total)
        << "non-finite tensor at " << (n_total - n_finite) << " atoms";

    // Spot check: SER residues that hit the DB should carry the
    // orca_input_orientation tag (project SER PBE regen). 1UBQ has
    // 4 SER residues. Not all may hit (e.g. terminal SER skipped).
    int n_ser_attempted = 0, n_ser_matched_pbe = 0;
    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        if (protein->ResidueAt(ri).type != AminoAcid::SER) continue;
        ++n_ser_attempted;
        const auto& m = tbb->ResidueMatches()[ri];
        if (m.calc_id == 0) continue;
        if (m.frame_type == "orca_input_orientation") {
            ++n_ser_matched_pbe;
        }
    }
    EXPECT_GE(n_ser_attempted, 1) << "1UBQ should have at least 1 SER";
    EXPECT_GT(n_ser_matched_pbe, 0)
        << "SER residues that match should hit ORCA PBE rows "
        << "(frame_type=orca_input_orientation), "
        << "got " << n_ser_matched_pbe << "/" << n_ser_attempted;

    // Test NPY emission.
    const std::string out_dir = "/tmp/tripeptide_bb_smoke_out";
    fs::create_directories(out_dir);
    int n_npy = tbb->WriteFeatures(conf, out_dir);
    EXPECT_EQ(n_npy, 4);
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_bb_shielding.npy"));
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_bb_residual_vec.npy"));
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_bb_match_distance.npy"));
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_bb_method_tag.npy"));

    // Residual sanity (post-Kabsch). The 3-point N/CA/C Kabsch has a
    // mean ~0.5 Å fit residual (chi-grid coarseness on flanking
    // residues' contribution to backbone C/O positions); individual
    // BB atoms can be 1-2 Å off. Above 3 Å on a backbone atom
    // indicates a Kabsch failure or mismapping — fail loud.
    //
    // The residual itself IS the ML feature (residual_vec on
    // ConformationAtom) — we don't gate ML emission on it. This test
    // just sanity-checks the distribution.
    int n_bb_extreme = 0;
    double max_bb_residual = 0.0;
    int n_bb_total = 0;
    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        const auto& res = protein->ResidueAt(ri);
        for (size_t slot : {res.N, res.CA, res.C, res.O}) {
            if (slot == Residue::NONE) continue;
            const auto& ca = conf.AtomAt(slot);
            if (!ca.tripeptide_bb_has_match) continue;
            ++n_bb_total;
            const double d = ca.tripeptide_bb_match_distance;
            if (d > max_bb_residual) max_bb_residual = d;
            if (d > 3.0) ++n_bb_extreme;
        }
    }
    EXPECT_GT(n_bb_total, 100)
        << "expected backbone atoms to mostly get tensors";
    EXPECT_LT(n_bb_extreme, 5)
        << "more than 5 backbone atoms have post-Kabsch residual > 3 Å "
        << "— Kabsch failure or substrate disagreement (max="
        << max_bb_residual << " Å)";
}

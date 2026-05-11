// Smoke test for TripeptideNeighborShieldingResult.
//
// Per Larsen 2015 Eq 3: each residue i receives Δσ_BB^{i-1}(i) +
// Δσ_BB^{i+1}(i) — read at the flanking ALA cap atoms of the
// (i±1)-centered tripeptides, with the AAA reference at standard
// angles (φ_std=-120°, ψ_std=140°) subtracted.
//
// On 1UBQ_pm6dh3plus.pdb we expect:
//   - residues_any ≥ 60 (almost every internal residue receives at
//     least one direction's contribution)
//   - atoms_accumulated > 200 (~7 cap atoms per residue × 60 residues)
//   - residual_vec_prev / residual_vec_next populated on most BB
//     atoms whose contributing direction hit
//   - SER residues' i-1 / i+1 directions, if SER itself is i±1, would
//     hit orca_input_orientation rows; assert at least one
//     ResidueMatch has frame_type=orca_input_orientation

#include "TestEnvironment.h"
#include <gtest/gtest.h>

#include "AminoAcidType.h"
#include "ConformationAtom.h"
#include "EnrichmentResult.h"
#include "GeometryResult.h"
#include "PdbFileReader.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "RuntimeEnvironment.h"
#include "Session.h"
#include "SpatialIndexResult.h"
#include "TripeptideDftTable.h"
#include "TripeptideNeighborShieldingResult.h"

#include <cmath>
#include <filesystem>
#include <string>

using namespace nmr;
namespace fs = std::filesystem;


namespace {
constexpr const char* kLarsen1UbqPm6 =
    "/mnt/expansion/larsen_archive/structures/1UBQ_pm6dh3plus.pdb";
}


class TripeptideNeighborShieldingTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (RuntimeEnvironment::TensorCs15Dsn().empty()) {
            GTEST_SKIP() << "tensorcs15 DSN not configured";
        }
        if (!fs::exists(kLarsen1UbqPm6)) {
            GTEST_SKIP() << "Larsen 1UBQ PM6-D3H+ PDB not at "
                         << kLarsen1UbqPm6;
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


TEST_F(TripeptideNeighborShieldingTest, RunsOn1UbqPm6) {
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

    auto tn = TripeptideNeighborShieldingResult::Compute(
        conf, *session.TripeptideDftTablePtr());
    ASSERT_NE(tn, nullptr);

    EXPECT_GE(tn->ResiduesWithAnyNeighbor(), 60)
        << "expected most internal residues to get >=1 neighbor "
           "contribution; got " << tn->ResiduesWithAnyNeighbor();
    EXPECT_GT(tn->AtomsAccumulated(), 200)
        << "expected >200 per-atom Δσ accumulations; got "
        << tn->AtomsAccumulated();

    // Frame_type discriminator: SER residues as i-1 or i+1 neighbors
    // should produce ResidueMatch entries with frame_type ==
    // "orca_input_orientation" (the project SER PBE regen).
    int n_pbe_dir = 0;
    int n_opbe_dir = 0;
    for (const auto& m : tn->ResidueMatches()) {
        for (const auto& ft : {m.prev_frame_type, m.next_frame_type}) {
            if (ft == "orca_input_orientation") ++n_pbe_dir;
            else if (ft == "gaussian_standard_orientation") ++n_opbe_dir;
        }
    }
    EXPECT_GT(n_pbe_dir, 0)
        << "at least one SER-side neighbor lookup should hit ORCA PBE; "
        << "got " << n_pbe_dir << " PBE vs " << n_opbe_dir << " OPBE";
    EXPECT_GT(n_opbe_dir, n_pbe_dir)
        << "OPBE should dominate (only SER is PBE); "
        << "got " << n_opbe_dir << " OPBE, " << n_pbe_dir << " PBE";

    // Residual sanity: per-direction residual_vecs populated on
    // matched atoms. The cap-side Kabsch aligns flanking-ALA N/CA/C
    // onto residue i's N/CA/C; the AAA reference is also at standard
    // angles, so AXA and AAA residuals should both be small.
    int n_atoms_with_neighbor = 0;
    int n_extreme_prev = 0, n_extreme_next = 0;
    double max_prev = 0.0, max_next = 0.0;
    for (size_t i = 0; i < conf.AtomCount(); ++i) {
        const auto& ca = conf.AtomAt(i);
        if (!ca.tripeptide_neighbor_has_match) continue;
        ++n_atoms_with_neighbor;
        const double pr = ca.tripeptide_neighbor_residual_vec_prev.norm();
        const double nr = ca.tripeptide_neighbor_residual_vec_next.norm();
        if (pr > max_prev) max_prev = pr;
        if (nr > max_next) max_next = nr;
        if (pr > 3.0) ++n_extreme_prev;
        if (nr > 3.0) ++n_extreme_next;
    }
    EXPECT_GT(n_atoms_with_neighbor, 200);
    EXPECT_LT(n_extreme_prev, 20)
        << "more than 20 atoms have prev-direction residual > 3 Å "
        << "(max=" << max_prev << " Å) — likely Kabsch issue";
    EXPECT_LT(n_extreme_next, 20)
        << "more than 20 atoms have next-direction residual > 3 Å "
        << "(max=" << max_next << " Å) — likely Kabsch issue";

    // NPY emission.
    const std::string out_dir = "/tmp/tripeptide_neighbor_smoke_out";
    fs::create_directories(out_dir);
    int n_npy = tn->WriteFeatures(conf, out_dir);
    EXPECT_EQ(n_npy, 3);
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_neighbor_shielding.npy"));
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_neighbor_residual_vec_prev.npy"));
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_neighbor_residual_vec_next.npy"));
}

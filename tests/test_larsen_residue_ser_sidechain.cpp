// SER OG ↔ backbone O regression test for LarsenResidue perception.
//
// The 2026-05-10 pre-perception heuristic mis-bound SER's sidechain OG
// (Locant::Gamma) to the protein's backbone O slot because the
// canonical Gaussian SER ordering interleaves the BB carbonyl atoms
// inside the sidechain run. Perception derives identity from the bond
// graph rather than position, so the typed identity for OG carries
// (Element::O, Locant::Gamma) and for BB O carries (Element::O,
// BackboneRole::CarbonylOxygen) — the two cannot collide.
//
// This test queries a SER row from the tensorcs15 replica, runs
// perception, and asserts the per-atom typed identities are correct.
//
// Skipped if the tensorcs15 DSN is not configured.

#include "TestEnvironment.h"
#include <gtest/gtest.h>

#include "AminoAcidType.h"
#include "LarsenResidue.h"
#include "RuntimeEnvironment.h"
#include "SemanticEnums.h"
#include "Session.h"
#include "TripeptideDftTable.h"

#include <string>

using namespace nmr;


class LarsenResidueSerSidechainTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (RuntimeEnvironment::TensorCs15Dsn().empty()) {
            GTEST_SKIP() << "tensorcs15 DSN not configured";
        }
        ASSERT_EQ(session.LoadTripeptideDftTable(), kOk);
        ASSERT_TRUE(session.HasTripeptideDftTable());
    }

    Session session;
};


TEST_F(LarsenResidueSerSidechainTest, OgAndBackboneOHaveDistinctIdentities) {
    const TripeptideDftTable* table = session.TripeptideDftTablePtr();
    ASSERT_NE(table, nullptr);

    // ASA SER (ORCA PBE rows). Default chi1=-60 hits the standard grid.
    TripeptideDftRecord rec = table->QueryNearest(
        'S', /*phi=*/-60.0, /*psi=*/120.0,
        /*chi1=*/-60.0, /*chi2=*/0.0, /*chi3=*/0.0, /*chi4=*/0.0,
        /*n_chi_axes=*/1);
    ASSERT_TRUE(rec.IsHit()) << "no SER pose at phi=-60,psi=120,chi1=-60";
    ASSERT_TRUE(rec.larsen.has_value());

    const LarsenResidue& cen = rec.larsen->central;
    ASSERT_EQ(cen.residue, AminoAcid::SER);
    ASSERT_EQ(cen.atoms.size(), 11u);  // BB6 + CB+HB2+HB3+OG+HG = 11

    // Find the SER OG and the backbone O by typed identity.
    // OG: Element::O, Locant::Gamma, BackboneRole::None.
    // O:  Element::O, BackboneRole::CarbonylOxygen.
    int og_idx = -1;
    int bb_o_idx = -1;
    for (std::size_t i = 0; i < cen.atoms.size(); ++i) {
        const auto& id = cen.atoms[i].identity;
        if (id.element != Element::O) continue;
        if (id.backbone_role == BackboneRole::CarbonylOxygen) {
            bb_o_idx = static_cast<int>(i);
        } else if (id.locant == Locant::Gamma) {
            og_idx = static_cast<int>(i);
        }
    }
    ASSERT_GE(og_idx,   0) << "no OG (Locant::Gamma) found in SER central";
    ASSERT_GE(bb_o_idx, 0) << "no BB O (CarbonylOxygen) found in SER central";

    EXPECT_NE(og_idx, bb_o_idx)
        << "OG and BB O collapsed to the same atom — perception bug "
           "(must carry distinct typed identities)";

    // Slot cache should pick up BB O via its typed BackboneRole, not OG.
    ASSERT_GE(cen.O_idx, 0);
    EXPECT_EQ(cen.O_idx, bb_o_idx)
        << "central.O_idx pointed at sidechain OG instead of BB O";

    // Final invariant: LookupByIdentity for the BB-O identity returns
    // the BB-O index, not OG's index.
    AtomMechanicalIdentity bb_o_id;
    bb_o_id.element       = Element::O;
    bb_o_id.backbone_role = BackboneRole::CarbonylOxygen;
    EXPECT_EQ(cen.LookupByIdentity(bb_o_id), bb_o_idx);

    AtomMechanicalIdentity og_id;
    og_id.element = Element::O;
    og_id.locant  = Locant::Gamma;
    EXPECT_EQ(cen.LookupByIdentity(og_id), og_idx);
}

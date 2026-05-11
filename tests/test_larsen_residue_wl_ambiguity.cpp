// K=3 Weisfeiler-Lehman ambiguity regression test for LarsenResidue.
//
// canonical_assignment_ambiguous is the per-PerAtom flag that distinguishes
//   (a) singleton WL signature classes (chemistry-deterministic — strict
//       identity match required downstream), from
//   (b) multi-atom WL signature classes (graph-automorphic — relaxed
//       identity match + nearest-spatial tiebreak downstream).
//
// The previous always-relaxed path (round-3 fix for the PHE/TYR CD/CE
// scramble) dropped BranchAddress + DiastereotopicIndex for every atom,
// including chemistry-distinct branches like ILE CG1/CG2 — silently
// allowing nearest-spatial to swap them under non-canonical chi
// orientations. This test locks in the corrected dispatch by asserting:
//
//   - ILE CG1 and CG2 perceive with canonical_assignment_ambiguous=false.
//     K=3 WL distinguishes them at K=1 (CG1's neighbour set includes
//     the chain-extension CD1; CG2's is a terminal methyl).
//   - PHE CD1/CD2/CE1/CE2 perceive with canonical_assignment_ambiguous=true.
//     They are graph-automorphic about the para axis; no K splits them.
//   - PHE CZ (para) perceives with canonical_assignment_ambiguous=false
//     (singleton — only one para carbon).
//
// Skipped if the tensorcs15 DSN is not configured (test machine doesn't
// have the local DB).

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


class LarsenResidueWlAmbiguityTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (RuntimeEnvironment::TensorCs15Dsn().empty()) {
            GTEST_SKIP() << "tensorcs15 DSN not configured";
        }
        ASSERT_EQ(session.LoadTripeptideDftTable(), kOk)
            << session.LastError();
        ASSERT_TRUE(session.HasTripeptideDftTable());
    }

    Session session;
};


// ILE CG1 (Locant::Gamma, BranchAddress.outer=1; methylene C extended
// to CD1) and CG2 (Locant::Gamma, BranchAddress.outer=2; terminal methyl)
// must perceive as chemistry-distinct singletons. K=3 WL distinguishes
// them at K=1 — different neighbour-element multisets.
TEST_F(LarsenResidueWlAmbiguityTest, IleCg1AndCg2AreChemistryDistinct) {
    const TripeptideDftTable* table = session.TripeptideDftTablePtr();
    ASSERT_NE(table, nullptr);

    TripeptideDftRecord rec = table->QueryNearest(
        'I', /*phi=*/-180.0, /*psi=*/-180.0,
        /*chi1=*/0.0, /*chi2=*/0.0, /*chi3=*/0.0, /*chi4=*/0.0,
        /*n_chi_axes=*/0);
    ASSERT_TRUE(rec.IsHit()) << "no ILE row at phi=-180,psi=-180";
    ASSERT_TRUE(rec.larsen.has_value());

    const LarsenResidue& cen = rec.larsen->central;
    ASSERT_EQ(cen.residue, AminoAcid::ILE);

    int cg1_idx = -1, cg2_idx = -1;
    for (std::size_t i = 0; i < cen.atoms.size(); ++i) {
        const auto& id = cen.atoms[i].identity;
        if (id.element != Element::C)        continue;
        if (id.locant  != Locant::Gamma)     continue;
        if (id.backbone_role != BackboneRole::None) continue;
        if (id.branch.outer == 1) cg1_idx = static_cast<int>(i);
        if (id.branch.outer == 2) cg2_idx = static_cast<int>(i);
    }
    ASSERT_GE(cg1_idx, 0) << "ILE CG1 not found by typed identity";
    ASSERT_GE(cg2_idx, 0) << "ILE CG2 not found by typed identity";

    EXPECT_FALSE(cen.atoms[cg1_idx].canonical_assignment_ambiguous)
        << "ILE CG1 must perceive as a singleton WL class — chemistry "
           "extends to CD1, K=1 WL splits it from CG2's terminal methyl";
    EXPECT_FALSE(cen.atoms[cg2_idx].canonical_assignment_ambiguous)
        << "ILE CG2 must perceive as a singleton WL class — terminal "
           "methyl, K=1 WL splits it from CG1's methylene+CD1 chain";
}


// PHE CD1/CD2/CE1/CE2 are graph-automorphic about the ring's para axis;
// no K rounds of WL split them. They must perceive ambiguous=true so
// AssembleCentralTyped relaxes BranchAddress and resolves the within-
// pair assignment by nearest-spatial.
//
// PHE CZ (para carbon, singleton) must perceive ambiguous=false.
TEST_F(LarsenResidueWlAmbiguityTest, PheCdAndCeAreGraphAutomorphic) {
    const TripeptideDftTable* table = session.TripeptideDftTablePtr();
    ASSERT_NE(table, nullptr);

    TripeptideDftRecord rec = table->QueryNearest(
        'F', /*phi=*/-180.0, /*psi=*/-180.0,
        /*chi1=*/0.0, /*chi2=*/0.0, /*chi3=*/0.0, /*chi4=*/0.0,
        /*n_chi_axes=*/0);
    ASSERT_TRUE(rec.IsHit()) << "no PHE row at phi=-180,psi=-180";
    ASSERT_TRUE(rec.larsen.has_value());

    const LarsenResidue& cen = rec.larsen->central;
    ASSERT_EQ(cen.residue, AminoAcid::PHE);

    int cd1_idx = -1, cd2_idx = -1, ce1_idx = -1, ce2_idx = -1, cz_idx = -1;
    for (std::size_t i = 0; i < cen.atoms.size(); ++i) {
        const auto& id = cen.atoms[i].identity;
        if (id.element != Element::C)              continue;
        if (id.backbone_role != BackboneRole::None) continue;
        if (id.locant == Locant::Delta && id.branch.outer == 1)
            cd1_idx = static_cast<int>(i);
        if (id.locant == Locant::Delta && id.branch.outer == 2)
            cd2_idx = static_cast<int>(i);
        if (id.locant == Locant::Epsilon && id.branch.outer == 1)
            ce1_idx = static_cast<int>(i);
        if (id.locant == Locant::Epsilon && id.branch.outer == 2)
            ce2_idx = static_cast<int>(i);
        if (id.locant == Locant::Zeta) cz_idx = static_cast<int>(i);
    }
    ASSERT_GE(cd1_idx, 0) << "PHE CD1 not found";
    ASSERT_GE(cd2_idx, 0) << "PHE CD2 not found";
    ASSERT_GE(ce1_idx, 0) << "PHE CE1 not found";
    ASSERT_GE(ce2_idx, 0) << "PHE CE2 not found";
    ASSERT_GE(cz_idx,  0) << "PHE CZ not found";

    EXPECT_TRUE(cen.atoms[cd1_idx].canonical_assignment_ambiguous)
        << "PHE CD1 must perceive as ambiguous — graph-automorphic with CD2";
    EXPECT_TRUE(cen.atoms[cd2_idx].canonical_assignment_ambiguous)
        << "PHE CD2 must perceive as ambiguous — graph-automorphic with CD1";
    EXPECT_TRUE(cen.atoms[ce1_idx].canonical_assignment_ambiguous)
        << "PHE CE1 must perceive as ambiguous — graph-automorphic with CE2";
    EXPECT_TRUE(cen.atoms[ce2_idx].canonical_assignment_ambiguous)
        << "PHE CE2 must perceive as ambiguous — graph-automorphic with CE1";

    EXPECT_FALSE(cen.atoms[cz_idx].canonical_assignment_ambiguous)
        << "PHE CZ must perceive as singleton — unique para carbon, no "
           "automorphic sibling";
}


// `QueryNearest` adds `ORDER BY calc_id ASC LIMIT 1` to the chi-fallback
// SQL. Without the ORDER BY, the row returned depends on the planner's
// emission order — i.e., session-arbitrary. This test verifies that
// when the chi-fallback drops to a depth that admits multiple rows
// (n_chi=0 against a 4-chi residue like ARG), the chosen calc_id is
// stable across consecutive queries.
TEST_F(LarsenResidueWlAmbiguityTest, ChiFallbackIsDeterministic) {
    const TripeptideDftTable* table = session.TripeptideDftTablePtr();
    ASSERT_NE(table, nullptr);

    // ARG at (phi=-180, psi=-180) with n_chi_axes=0 → WHERE has only
    // tripeptide + phi + psi terms, so every ARA row at that grid
    // point matches; the DB has many (chi1, chi2, chi3, chi4) combos
    // at that phi/psi point.
    TripeptideDftRecord rec1 = table->QueryNearest(
        'R', /*phi=*/-180.0, /*psi=*/-180.0,
        /*chi1=*/0.0, /*chi2=*/0.0, /*chi3=*/0.0, /*chi4=*/0.0,
        /*n_chi_axes=*/0);
    TripeptideDftRecord rec2 = table->QueryNearest(
        'R', /*phi=*/-180.0, /*psi=*/-180.0,
        /*chi1=*/0.0, /*chi2=*/0.0, /*chi3=*/0.0, /*chi4=*/0.0,
        /*n_chi_axes=*/0);

    ASSERT_TRUE(rec1.IsHit()) << "no ARA row at phi=-180,psi=-180";
    ASSERT_TRUE(rec2.IsHit()) << "no ARA row at phi=-180,psi=-180 on re-query";
    EXPECT_EQ(rec1.calc_id, rec2.calc_id)
        << "chi-fallback returned different rows across two consecutive "
           "calls at the same grid point — `ORDER BY calc_id ASC` was "
           "expected to make this deterministic. Got calc_id=" <<
           rec1.calc_id << " then " << rec2.calc_id;
    // Also sanity check that ORDER BY calc_id ASC is the lowest:
    // re-query under the same conditions and assert the chosen
    // calc_id is the minimum across the visible set. We can't easily
    // enumerate the full match set without raw SQL access, so the
    // determinism check above is the load-bearing assertion. The
    // lowest-calc_id property is a side-effect.
}

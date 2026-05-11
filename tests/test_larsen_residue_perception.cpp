// Perception parity test for LarsenResidue / PerceiveLarsenTripeptide.
//
// For every (tripeptide, frame_type) combination present in the
// tensorcs15 DB, fetch one row and verify that:
//
//   1. PerceiveLarsenTripeptide returns a non-empty optional.
//   2. The result has exactly 5 pieces (ACE, NCapAla, Central, CCapAla, NME)
//      with all required slots populated for each piece's Kind.
//   3. Per piece, total atom count matches canonical chemistry.
//   4. The central piece's BB N/CA/C are 3 distinct local atoms with
//      element N/C/C and the canonical N-CA, CA-C bond distances.
//   5. The central piece's residue field matches the queried letter.
//
// This is a DB-only test — it does not require any PDB file. The Python
// POC at scripts/perceive_larsen_tripeptide.py covers the same ground
// at the algorithm level; this test locks the C++ port behavior in CI.
//
// Skipped if the tensorcs15 DSN is not configured (test machine doesn't
// have the local DB).

#include "TestEnvironment.h"
#include <gtest/gtest.h>

#include "AminoAcidType.h"
#include "LarsenResidue.h"
#include "OperationLog.h"
#include "RuntimeEnvironment.h"
#include "Session.h"
#include "TripeptideDftTable.h"
#include "Types.h"

#include <cmath>
#include <string>
#include <vector>

using namespace nmr;


namespace {

// The (tripeptide, frame_type, central_letter) test matrix. One row
// per DB (tripeptide, frame_type) combination. We probe each combo
// at (phi=-180, psi=-180) with n_chi_axes=0 — verified via psql to
// exist for all 20 (tripeptide, frame_type) pairs (COUNT(DISTINCT
// tripeptide) = 20 at that grid point). QueryNearest does an exact
// SQL match rounded to grid, so picking a point that always exists
// keeps the test deterministic regardless of which chi values are
// physically meaningful for each residue.
struct Combo {
    const char* tripeptide;
    const char* frame_type;
    char        letter;
};

constexpr Combo kCombos[] = {
    {"AAA", "gaussian_standard_orientation", 'A'},
    {"ACA", "gaussian_standard_orientation", 'C'},
    {"ADA", "gaussian_standard_orientation", 'D'},
    {"AEA", "gaussian_standard_orientation", 'E'},
    {"AFA", "gaussian_standard_orientation", 'F'},
    {"AGA", "gaussian_standard_orientation", 'G'},
    {"AHA", "gaussian_standard_orientation", 'H'},
    {"AIA", "gaussian_standard_orientation", 'I'},
    {"AKA", "gaussian_standard_orientation", 'K'},
    {"ALA", "gaussian_standard_orientation", 'L'},
    {"AMA", "gaussian_standard_orientation", 'M'},
    {"ANA", "gaussian_standard_orientation", 'N'},
    {"APA", "gaussian_standard_orientation", 'P'},
    {"AQA", "gaussian_standard_orientation", 'Q'},
    {"ARA", "gaussian_standard_orientation", 'R'},
    {"ASA", "orca_input_orientation",        'S'},
    {"ATA", "gaussian_standard_orientation", 'T'},
    {"AVA", "gaussian_standard_orientation", 'V'},
    {"AWA", "gaussian_standard_orientation", 'W'},
    {"AYA", "gaussian_standard_orientation", 'Y'},
};

constexpr int kNCombos = sizeof(kCombos) / sizeof(kCombos[0]);

AminoAcid LetterToAa(char letter) {
    for (const auto& t : AllAminoAcidTypes()) {
        if (t.one_letter_code == letter) return t.index;
    }
    return AminoAcid::Unknown;
}

}  // namespace


class LarsenResiduePerceptionTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (RuntimeEnvironment::TensorCs15Dsn().empty()) {
            GTEST_SKIP() << "tensorcs15 DSN not configured "
                            "(set [databases].tensorcs15 in "
                            "~/.nmr_tools.toml)";
        }
        ASSERT_EQ(session.LoadTripeptideDftTable(), kOk)
            << session.LastError();
        ASSERT_TRUE(session.HasTripeptideDftTable());
    }

    Session session;
};


// Iterate every (tripeptide, frame_type) combination and verify
// perception succeeds with all canonical atoms assigned.
TEST_F(LarsenResiduePerceptionTest, AllCombinationsPerceiveCleanly) {
    const TripeptideDftTable* table = session.TripeptideDftTablePtr();
    ASSERT_NE(table, nullptr);

    int n_ok = 0;
    int n_skipped_db_miss = 0;
    for (const Combo& c : kCombos) {
        const AminoAcid aa = LetterToAa(c.letter);
        ASSERT_NE(aa, AminoAcid::Unknown) << "unknown letter " << c.letter;

        // phi=-180, psi=-180, n_chi_axes=0 → drop all chi from
        // WHERE clause. Verified via psql that all 20 (tripeptide,
        // frame_type) pairs have at least one row at this grid point.
        TripeptideDftRecord rec = table->QueryNearest(
            c.letter, /*phi=*/-180.0, /*psi=*/-180.0,
            /*chi1=*/0.0, /*chi2=*/0.0, /*chi3=*/0.0, /*chi4=*/0.0,
            /*n_chi_axes=*/0);
        if (!rec.IsHit()) {
            ++n_skipped_db_miss;
            ADD_FAILURE() << c.tripeptide
                << " (frame_type=" << c.frame_type
                << "): no row at phi=-180,psi=-180 — DB schema drift?";
            continue;
        }

        EXPECT_EQ(rec.frame_type, c.frame_type)
            << c.tripeptide << ": frame_type expected " << c.frame_type
            << " got " << rec.frame_type;

        // QueryNearest perceives at fetch time; verify the result is
        // populated.
        ASSERT_TRUE(rec.larsen.has_value())
            << c.tripeptide << " (calc_id=" << rec.calc_id
            << "): perception returned nullopt";

        const LarsenTripeptide& trip = *rec.larsen;

        // (a) 5 pieces with required slots.
        EXPECT_EQ(trip.ace.kind,     LarsenResidue::Kind::AceCap);
        EXPECT_EQ(trip.n_cap.kind,   LarsenResidue::Kind::NCapAla);
        EXPECT_EQ(trip.central.kind, LarsenResidue::Kind::Central);
        EXPECT_EQ(trip.c_cap.kind,   LarsenResidue::Kind::CCapAla);
        EXPECT_EQ(trip.nme.kind,     LarsenResidue::Kind::NmeCap);

        EXPECT_TRUE(trip.ace.HasAllRequiredSlots())   << c.tripeptide;
        EXPECT_TRUE(trip.n_cap.HasAllRequiredSlots()) << c.tripeptide;
        EXPECT_TRUE(trip.central.HasAllRequiredSlots()) << c.tripeptide;
        EXPECT_TRUE(trip.c_cap.HasAllRequiredSlots()) << c.tripeptide;
        EXPECT_TRUE(trip.nme.HasAllRequiredSlots())   << c.tripeptide;

        // (b) Total atom count matches canonical AAA construct shape:
        //     ACE 6 + N-cap ALA 10 + X + C-cap ALA 10 + NME 6
        // X depends on the central residue and (for HIS) variant.
        EXPECT_EQ(trip.ace.atoms.size(),   6u) << c.tripeptide;
        EXPECT_EQ(trip.n_cap.atoms.size(), 10u) << c.tripeptide;
        EXPECT_EQ(trip.c_cap.atoms.size(), 10u) << c.tripeptide;
        EXPECT_EQ(trip.nme.atoms.size(),   6u) << c.tripeptide;
        EXPECT_GT(trip.central.atoms.size(), 0u) << c.tripeptide;

        // (c) Central residue assignment matches expected.
        EXPECT_EQ(trip.central.residue, aa) << c.tripeptide;

        // (d) Central BB N/CA/C are 3 distinct elements/atoms.
        const auto& cen = trip.central;
        ASSERT_GE(cen.N_idx,  0) << c.tripeptide;
        ASSERT_GE(cen.CA_idx, 0) << c.tripeptide;
        ASSERT_GE(cen.C_idx,  0) << c.tripeptide;
        ASSERT_GE(cen.O_idx,  0) << c.tripeptide;
        EXPECT_NE(cen.N_idx,  cen.CA_idx);
        EXPECT_NE(cen.CA_idx, cen.C_idx);
        EXPECT_NE(cen.N_idx,  cen.C_idx);
        EXPECT_EQ(cen.atoms[cen.N_idx].element,  Element::N);
        EXPECT_EQ(cen.atoms[cen.CA_idx].element, Element::C);
        EXPECT_EQ(cen.atoms[cen.C_idx].element,  Element::C);
        EXPECT_EQ(cen.atoms[cen.O_idx].element,  Element::O);

        // (e) Canonical N-CA and CA-C bond distances are in
        //     the standard peptide range (1.3–1.6 Å).
        const Vec3 N  = cen.atoms[cen.N_idx].position;
        const Vec3 CA = cen.atoms[cen.CA_idx].position;
        const Vec3 C  = cen.atoms[cen.C_idx].position;
        const double d_NCA = (CA - N).norm();
        const double d_CAC = (C - CA).norm();
        EXPECT_GT(d_NCA, 1.3) << c.tripeptide << " N-CA = " << d_NCA;
        EXPECT_LT(d_NCA, 1.6) << c.tripeptide << " N-CA = " << d_NCA;
        EXPECT_GT(d_CAC, 1.3) << c.tripeptide << " CA-C = " << d_CAC;
        EXPECT_LT(d_CAC, 1.6) << c.tripeptide << " CA-C = " << d_CAC;

        ++n_ok;
    }

    EXPECT_EQ(n_skipped_db_miss, 0)
        << "DB-miss fixtures: " << n_skipped_db_miss
        << " — update kCombos default poses to match the available grid";
    EXPECT_EQ(n_ok, kNCombos)
        << "passed: " << n_ok << "/" << kNCombos;
}

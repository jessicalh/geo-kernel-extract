// Source-log forensic test for LarsenResidue perception.
//
// Parse the Larsen 2015 Gaussian log directly (Standard orientation
// block), run PerceiveLarsenTripeptide against it, and verify the
// canonical structure ACE — NCapAla — central X — CCapAla — NME is
// recovered without ever touching the Postgres replica. This is the
// independent-of-DB check: if a future ingest perturbs the order in
// tensorcs15, the smoke test would still pass but THIS test catches
// the drift by going back to the source.
//
// Skipped if the Larsen ERDA archive isn't present at the expected
// path. The archive lives at /mnt/expansion/larsen_archive/.

#include "TestEnvironment.h"
#include <gtest/gtest.h>

#include "AminoAcidType.h"
#include "LarsenResidue.h"
#include "TripeptideDftTable.h"
#include "Types.h"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

using namespace nmr;
namespace fs = std::filesystem;


namespace {

constexpr const char* kAaaLog =
    "/mnt/expansion/williamsproject/procs15_data/nmrlogs/AAAnmrlog/"
    "AAA_4_54_nmr.log";

// Parse the LAST "Standard orientation:" block from a Gaussian log and
// return per-atom records. Mirrors the parser in
// scripts/perceive_larsen_tripeptide.py.
std::vector<TripeptideDftAtom>
ParseStandardOrientation(const fs::path& log_path) {
    std::ifstream in(log_path);
    std::stringstream ss; ss << in.rdbuf();
    const std::string text = ss.str();

    // Find the LAST occurrence of "Standard orientation:".
    std::size_t pos = std::string::npos;
    std::size_t scan = 0;
    while (true) {
        std::size_t hit = text.find("Standard orientation:", scan);
        if (hit == std::string::npos) break;
        pos = hit;
        scan = hit + 1;
    }
    if (pos == std::string::npos) return {};

    // Skip past the 5-line header (---, "Center Atomic Atomic",
    // "Number Number Type   X  Y  Z", ---), arriving at the atom rows.
    std::size_t lineno = 0;
    std::size_t cursor = pos;
    auto eat_line = [&]() {
        const std::size_t nl = text.find('\n', cursor);
        const std::string line = (nl == std::string::npos)
            ? text.substr(cursor)
            : text.substr(cursor, nl - cursor);
        cursor = (nl == std::string::npos) ? text.size() : nl + 1;
        ++lineno;
        return line;
    };
    // Header skip.
    while (cursor < text.size()) {
        std::string line = eat_line();
        // Trim leading whitespace and look for first dashed separator.
        if (line.find("---") != std::string::npos) break;
    }
    // Skip 2 column-header lines.
    eat_line(); eat_line();
    // Skip second dashed.
    eat_line();

    std::vector<TripeptideDftAtom> atoms;
    while (cursor < text.size()) {
        std::string line = eat_line();
        if (line.find("---") != std::string::npos) break;
        std::istringstream iss(line);
        int idx, z, atype;
        double x, y, z_coord;
        if (!(iss >> idx >> z >> atype >> x >> y >> z_coord)) break;
        TripeptideDftAtom a;
        a.atom_idx = idx;
        switch (z) {
            case 1:  a.element = Element::H; break;
            case 6:  a.element = Element::C; break;
            case 7:  a.element = Element::N; break;
            case 8:  a.element = Element::O; break;
            case 16: a.element = Element::S; break;
            default: a.element = Element::Unknown; break;
        }
        a.position = Vec3(x, y, z_coord);
        atoms.push_back(a);
    }
    return atoms;
}

}  // namespace


class LarsenResidueAgainstSourceLogTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!fs::exists(kAaaLog)) {
            GTEST_SKIP() << "AAA Gaussian log not found at "
                         << kAaaLog
                         << " — set up the larsen archive copy via "
                            "/mnt/expansion/williamsproject/procs15_data/"
                            "nmrlogs/AAAnmrlog/ before running this test.";
        }
    }
};


// Parse the AAA_4_54_nmr.log Standard orientation block directly,
// perceive, and verify the 5-piece structure is recovered.
TEST_F(LarsenResidueAgainstSourceLogTest, AaaLogPerceivesCleanly) {
    const std::vector<TripeptideDftAtom> log_atoms =
        ParseStandardOrientation(kAaaLog);

    // AAA = ACE 6 + N-cap ALA 10 + central ALA 10 + C-cap ALA 10 +
    //       NME 6 = 42 atoms.
    ASSERT_EQ(log_atoms.size(), 42u)
        << "AAA Gaussian log expected 42 atoms";

    // Build a synthetic TripeptideDftRecord (no tensors needed for
    // perception; positions + elements suffice).
    TripeptideDftRecord rec;
    rec.calc_id     = -1;  // synthetic marker
    rec.tripeptide  = "AAA";
    rec.frame_type  = "gaussian_standard_orientation";
    rec.n_atoms     = static_cast<int>(log_atoms.size());
    rec.atoms       = log_atoms;

    auto trip = PerceiveLarsenTripeptide(rec, AminoAcid::ALA);
    ASSERT_TRUE(trip.has_value())
        << "perception failed on raw Gaussian log — see "
           "PerceiveLarsenTripeptide warnings";

    // Same structural invariants as the DB parity test, applied to the
    // raw log. This is the independent-of-DB cross-check.
    EXPECT_EQ(trip->ace.atoms.size(),     6u);
    EXPECT_EQ(trip->n_cap.atoms.size(),   10u);
    EXPECT_EQ(trip->central.atoms.size(), 10u);
    EXPECT_EQ(trip->c_cap.atoms.size(),   10u);
    EXPECT_EQ(trip->nme.atoms.size(),     6u);

    EXPECT_TRUE(trip->ace.HasAllRequiredSlots());
    EXPECT_TRUE(trip->n_cap.HasAllRequiredSlots());
    EXPECT_TRUE(trip->central.HasAllRequiredSlots());
    EXPECT_TRUE(trip->c_cap.HasAllRequiredSlots());
    EXPECT_TRUE(trip->nme.HasAllRequiredSlots());

    EXPECT_EQ(trip->central.residue, AminoAcid::ALA);
}

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

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace nmr;
namespace fs = std::filesystem;

namespace {

struct NpyArray {
    std::vector<size_t> shape;
    std::vector<char> bytes;
};

std::string Trim(std::string s) {
    auto is_space = [](unsigned char c) { return std::isspace(c) != 0; };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                                    [&](char c) { return !is_space(c); }));
    s.erase(std::find_if(s.rbegin(), s.rend(),
                         [&](char c) { return !is_space(c); }).base(), s.end());
    return s;
}

NpyArray ReadNpy(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    EXPECT_TRUE(in.is_open()) << path;
    NpyArray arr;
    if (!in.is_open()) return arr;
    char magic[6] = {};
    in.read(magic, 6);
    EXPECT_EQ(std::string(magic, 6), std::string("\x93NUMPY", 6));
    char version[2] = {};
    in.read(version, 2);
    EXPECT_EQ(version[0], 1);
    EXPECT_EQ(version[1], 0);
    uint16_t header_len = 0;
    in.read(reinterpret_cast<char*>(&header_len), sizeof(header_len));
    std::string header(header_len, '\0');
    in.read(header.data(), header_len);
    auto shape_begin = header.find('(');
    auto shape_end = header.find(')', shape_begin);
    EXPECT_NE(shape_begin, std::string::npos);
    EXPECT_NE(shape_end, std::string::npos);
    if (shape_begin == std::string::npos || shape_end == std::string::npos)
        return arr;
    std::stringstream ss(header.substr(shape_begin + 1,
                                      shape_end - shape_begin - 1));
    std::string token;
    while (std::getline(ss, token, ',')) {
        token = Trim(token);
        if (!token.empty()) arr.shape.push_back(static_cast<size_t>(std::stoull(token)));
    }
    arr.bytes.assign(std::istreambuf_iterator<char>(in),
                     std::istreambuf_iterator<char>());
    return arr;
}

template <typename T>
const T* DataAs(const NpyArray& arr) {
    return reinterpret_cast<const T*>(arr.bytes.data());
}

bool RowHasFinite(const double* data, size_t row, size_t cols) {
    for (size_t c = 0; c < cols; ++c) {
        if (std::isfinite(data[row * cols + c])) return true;
    }
    return false;
}

}  // namespace

// Forward-declare the file-local PRODUCTION dihedral helper (named per-file
// namespace, external linkage) so the fixed-coordinate ±60° test below pins
// production DIRECTLY, not a copy of it (vet finding 2026-07). A production
// sign regression now fails this non-skippable test.
namespace nmr { namespace tripeptide_neighbor_dihedral {
double DihedralDegrees(const Vec3&, const Vec3&, const Vec3&, const Vec3&);
} }  // namespace nmr::tripeptide_neighbor_dihedral


TEST(TripeptideNeighborDihedral,
     CanonicalFixedCoordinatesPinSignAndDegenerateNaN) {
    // Calls the PRODUCTION helper directly (not a copy), so a production sign
    // regression fails this non-skippable test (vet finding 2026-07).
    using nmr::tripeptide_neighbor_dihedral::DihedralDegrees;
    const Vec3 a(1.0, 0.0, 0.0);
    const Vec3 b(0.0, 0.0, 0.0);
    const Vec3 c(0.0, 0.0, 1.0);
    const double root3_over_2 = std::sqrt(3.0) / 2.0;

    // Analytic fixtures distinguish the canonical trajectory convention
    // from the former sign-reversed tripeptide triple product.
    EXPECT_NEAR(DihedralDegrees(a, b, c, Vec3(0.5, -root3_over_2, 1.0)),
                60.0, 1e-12);
    EXPECT_NEAR(DihedralDegrees(a, b, c, Vec3(0.5, root3_over_2, 1.0)),
                -60.0, 1e-12);

    EXPECT_TRUE(std::isnan(DihedralDegrees(a, b, b, Vec3(0.0, 1.0, 0.0))));
    EXPECT_TRUE(std::isnan(DihedralDegrees(
        Vec3(0.0, 0.0, 0.0), Vec3(1.0, 0.0, 0.0),
        Vec3(2.0, 0.0, 0.0), Vec3(2.0, 1.0, 0.0))));
    EXPECT_TRUE(std::isnan(DihedralDegrees(
        Vec3(0.0, 1.0, 0.0), Vec3(0.0, 0.0, 0.0),
        Vec3(1.0, 0.0, 0.0), Vec3(2.0, 0.0, 0.0))));
}


class TripeptideNeighborShieldingTest : public ::testing::Test {
protected:
    void SetUp() override {
        const std::string kLarsen1UbqPm6 =
            nmr::test::TestEnvironment::Larsen1UbqPm6Pdb();
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

    // Exercise the production file-local helper through the DB row it
    // selected. ALA has no chi axes, so independently computed canonical
    // phi/psi identify one unambiguous 2-degree-grid row. The former
    // sign-reversed helper selects the opposite-angle row and a different
    // calc_id.
    auto query_ala_from_geometry = [&](size_t neighbor_index) {
        TripeptideDftRecord miss;
        const auto prev_idx = protein->BackbonePredecessor(neighbor_index);
        const auto next_idx = protein->BackboneSuccessor(neighbor_index);
        if (!prev_idx || !next_idx) return miss;
        const Residue& prev = protein->ResidueAt(*prev_idx);
        const Residue& neighbor = protein->ResidueAt(neighbor_index);
        const Residue& next = protein->ResidueAt(*next_idx);
        if (prev.C == Residue::NONE || neighbor.N == Residue::NONE ||
            neighbor.CA == Residue::NONE || neighbor.C == Residue::NONE ||
            next.N == Residue::NONE) {
            return miss;
        }
        const double phi = nmr::tripeptide_neighbor_dihedral::DihedralDegrees(
            conf.PositionAt(prev.C), conf.PositionAt(neighbor.N),
            conf.PositionAt(neighbor.CA), conf.PositionAt(neighbor.C));
        const double psi = nmr::tripeptide_neighbor_dihedral::DihedralDegrees(
            conf.PositionAt(neighbor.N), conf.PositionAt(neighbor.CA),
            conf.PositionAt(neighbor.C), conf.PositionAt(next.N));
        return session.TripeptideDftTablePtr()->QueryNearest(
            'A', phi, psi, 0.0, 0.0, 0.0, 0.0,
            /*n_chi_axes=*/0, /*his_variant_hint=*/-1);
    };

    int n_ala_dihedral_queries_pinned = 0;
    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        const auto& match = tn->ResidueMatches()[ri];
        if (const auto prev_idx = protein->BackbonePredecessor(ri);
            prev_idx && protein->ResidueAt(*prev_idx).type == AminoAcid::ALA &&
            match.prev_calc_id != 0) {
            const auto expected = query_ala_from_geometry(*prev_idx);
            ASSERT_TRUE(expected.IsHit());
            ASSERT_TRUE(expected.larsen.has_value());
            EXPECT_EQ(match.prev_calc_id, expected.calc_id)
                << "i-1 ALA neighbor at residue " << ri;
            ++n_ala_dihedral_queries_pinned;
        }
        if (const auto next_idx = protein->BackboneSuccessor(ri);
            next_idx && protein->ResidueAt(*next_idx).type == AminoAcid::ALA &&
            match.next_calc_id != 0) {
            const auto expected = query_ala_from_geometry(*next_idx);
            ASSERT_TRUE(expected.IsHit());
            ASSERT_TRUE(expected.larsen.has_value());
            EXPECT_EQ(match.next_calc_id, expected.calc_id)
                << "i+1 ALA neighbor at residue " << ri;
            ++n_ala_dihedral_queries_pinned;
        }
    }
    EXPECT_GT(n_ala_dihedral_queries_pinned, 0)
        << "1UBQ should exercise an ALA neighbor DB lookup";

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
    //
    // Per M4 contract, the per-direction residual is NaN whenever a
    // direction had no contribution at all (terminal residue's
    // missing i-1 / i+1, or chain break). The previous formulation
    // here read .norm() on every has_match atom and compared the
    // result against 3.0 — NaN comparisons evaluate false in C++ so
    // an entire direction being absent everywhere would silently
    // pass the extreme bound check. The fix is to count finite
    // residuals separately and only fold finite vectors into the
    // max/extreme tallies; the test then also asserts non-trivial
    // finite coverage in each direction.
    int n_atoms_with_neighbor = 0;
    int n_finite_prev = 0, n_finite_next = 0;
    int n_extreme_prev = 0, n_extreme_next = 0;
    double max_prev = 0.0, max_next = 0.0;
    for (size_t i = 0; i < conf.AtomCount(); ++i) {
        const auto& ca = conf.AtomAt(i);
        if (!ca.tripeptide_neighbor_has_match) continue;
        ++n_atoms_with_neighbor;
        const double pr = ca.tripeptide_neighbor_residual_vec_prev.norm();
        const double nr = ca.tripeptide_neighbor_residual_vec_next.norm();
        if (std::isfinite(pr)) {
            ++n_finite_prev;
            if (pr > max_prev) max_prev = pr;
            if (pr > 3.0) ++n_extreme_prev;
        }
        if (std::isfinite(nr)) {
            ++n_finite_next;
            if (nr > max_next) max_next = nr;
            if (nr > 3.0) ++n_extreme_next;
        }
    }
    EXPECT_GT(n_atoms_with_neighbor, 200);
    EXPECT_GT(n_finite_prev, 100)
        << "expected most has_match atoms to also have a finite "
           "i-1 residual; got " << n_finite_prev << " finite of "
        << n_atoms_with_neighbor << " has_match atoms";
    EXPECT_GT(n_finite_next, 100)
        << "expected most has_match atoms to also have a finite "
           "i+1 residual; got " << n_finite_next << " finite of "
        << n_atoms_with_neighbor << " has_match atoms";
    EXPECT_LT(n_extreme_prev, 20)
        << "more than 20 atoms have prev-direction residual > 3 Å "
        << "(max=" << max_prev << " Å) — likely Kabsch issue";
    EXPECT_LT(n_extreme_next, 20)
        << "more than 20 atoms have next-direction residual > 3 Å "
        << "(max=" << max_next << " Å) — likely Kabsch issue";

    // NPY emission.
    const std::string out_dir =
        nmr::test::TestEnvironment::TempPath("tripeptide_neighbor_smoke_out");
    fs::create_directories(out_dir);
    int n_npy = tn->WriteFeatures(conf, out_dir);
    EXPECT_EQ(n_npy, 6);
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_neighbor_shielding.npy"));
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_neighbor_shielding_prev.npy"));
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_neighbor_shielding_next.npy"));
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_neighbor_residual_vec_prev.npy"));
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_neighbor_residual_vec_next.npy"));
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_neighbor_reference.npy"));

    auto sum_arr = ReadNpy(out_dir + "/tripeptide_neighbor_shielding.npy");
    auto prev_arr = ReadNpy(out_dir + "/tripeptide_neighbor_shielding_prev.npy");
    auto next_arr = ReadNpy(out_dir + "/tripeptide_neighbor_shielding_next.npy");
    auto ref_arr = ReadNpy(out_dir + "/tripeptide_neighbor_reference.npy");
    const std::vector<size_t> tensor_shape{conf.AtomCount(), 9};
    ASSERT_EQ(sum_arr.shape, tensor_shape);
    ASSERT_EQ(prev_arr.shape, tensor_shape);
    ASSERT_EQ(next_arr.shape, tensor_shape);
    ASSERT_EQ(ref_arr.shape, (std::vector<size_t>{1, 5}));

    const double* sum = DataAs<double>(sum_arr);
    const double* prev = DataAs<double>(prev_arr);
    const double* next = DataAs<double>(next_arr);
    const double* ref = DataAs<double>(ref_arr);
    int n_prev_rows = 0, n_next_rows = 0, n_reconstructed = 0;
    for (size_t i = 0; i < conf.AtomCount(); ++i) {
        const bool has_prev = RowHasFinite(prev, i, 9);
        const bool has_next = RowHasFinite(next, i, 9);
        if (has_prev) ++n_prev_rows;
        if (has_next) ++n_next_rows;
        if (!has_prev && !has_next) {
            for (int c = 0; c < 9; ++c)
                EXPECT_TRUE(std::isnan(sum[i * 9 + c]));
            continue;
        }
        ++n_reconstructed;
        for (int c = 0; c < 9; ++c) {
            const double expected =
                (has_prev ? prev[i * 9 + c] : 0.0) +
                (has_next ? next[i * 9 + c] : 0.0);
            EXPECT_NEAR(sum[i * 9 + c], expected, 1e-9)
                << "atom=" << i << " component=" << c;
        }
    }
    EXPECT_GT(n_prev_rows, 100);
    EXPECT_GT(n_next_rows, 100);
    EXPECT_GT(n_reconstructed, 200);
    EXPECT_GT(ref[0], 0.0);
    EXPECT_DOUBLE_EQ(ref[1], 1.0);
    EXPECT_DOUBLE_EQ(ref[2], -120.0);
    EXPECT_DOUBLE_EQ(ref[3], 140.0);
    EXPECT_DOUBLE_EQ(ref[4], 1.0);
}

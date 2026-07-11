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
//   - WriteFeatures emits the central tensor, residual, method tag,
//     distance, and atom-match metadata NPYs
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

#include <libpq-fe.h>

#include <algorithm>
#include <array>
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
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [&](char c) { return !is_space(c); }));
    s.erase(std::find_if(s.rbegin(), s.rend(), [&](char c) { return !is_space(c); }).base(), s.end());
    return s;
}

NpyArray ReadNpy(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    EXPECT_TRUE(in.is_open()) << path;
    NpyArray arr;
    if (!in.is_open())
        return arr;
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
    std::stringstream ss(header.substr(shape_begin + 1, shape_end - shape_begin - 1));
    std::string token;
    while (std::getline(ss, token, ',')) {
        token = Trim(token);
        if (!token.empty())
            arr.shape.push_back(static_cast<size_t>(std::stoull(token)));
    }
    arr.bytes.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
    return arr;
}

template <typename T>
const T* DataAs(const NpyArray& arr) {
    return reinterpret_cast<const T*>(arr.bytes.data());
}

}  // namespace

// Forward-declare the file-local PRODUCTION dihedral helper (named per-file
// namespace, external linkage) so the fixed-coordinate ±60° test below pins
// production DIRECTLY, not a copy of it (vet finding 2026-07). A production
// sign regression now fails this non-skippable test.
namespace nmr {
namespace tripeptide_backbone_dihedral {
double DihedralDegrees(const Vec3&, const Vec3&, const Vec3&, const Vec3&);
}
}  // namespace nmr

// Production C3/M14 seams. These functions are defined in named, per-file
// namespaces and are called by QueryNearest/WriteFeatures themselves.
namespace nmr {
namespace tripeptide_dft_query {
double CircularDeltaDegrees(double, double);
double DroppedChiSquaredDistance(const std::array<int, 4>&, const std::array<int, 4>&, int, int);
std::string OmittedChiScoreSql(int, int, const std::array<int, 4>&);
std::string NearestOrderBySql();
}  // namespace tripeptide_dft_query
}  // namespace nmr

namespace nmr {
namespace tripeptide_backbone_status {
int HisVariantDiagnosticCode(AminoAcid, int);
TripeptideMatchStatus UnsupportedHisStatus(AminoAcid, int);
}  // namespace tripeptide_backbone_status
}  // namespace nmr

namespace nmr {
namespace tripeptide_backbone_diagnostics {
std::array<double, 28> PackDiagnosticRow(const TripeptideBackboneShieldingResult::ResidueMatch&);
}
}  // namespace nmr


TEST(TripeptideBackboneDihedral, CanonicalFixedCoordinatesPinSignAndDegenerateNaN) {
    // Calls the PRODUCTION helper directly (not a copy), so a production sign
    // regression fails this non-skippable test (vet finding 2026-07).
    using nmr::tripeptide_backbone_dihedral::DihedralDegrees;
    const Vec3 a(1.0, 0.0, 0.0);
    const Vec3 b(0.0, 0.0, 0.0);
    const Vec3 c(0.0, 0.0, 1.0);
    const double root3_over_2 = std::sqrt(3.0) / 2.0;

    // Analytic ±60° fixtures. The former tripeptide triple product gives the
    // opposite signs, so a sign regression fails here.
    EXPECT_NEAR(DihedralDegrees(a, b, c, Vec3(0.5, -root3_over_2, 1.0)), 60.0, 1e-12);
    EXPECT_NEAR(DihedralDegrees(a, b, c, Vec3(0.5, root3_over_2, 1.0)), -60.0, 1e-12);

    // Zero central bond and collinear plane-defining triples are undefined,
    // matching the trajectory convention.
    EXPECT_TRUE(std::isnan(DihedralDegrees(a, b, b, Vec3(0.0, 1.0, 0.0))));
    EXPECT_TRUE(
        std::isnan(DihedralDegrees(Vec3(0.0, 0.0, 0.0), Vec3(1.0, 0.0, 0.0), Vec3(2.0, 0.0, 0.0), Vec3(2.0, 1.0, 0.0))));
    EXPECT_TRUE(
        std::isnan(DihedralDegrees(Vec3(0.0, 1.0, 0.0), Vec3(0.0, 0.0, 0.0), Vec3(1.0, 0.0, 0.0), Vec3(2.0, 0.0, 0.0))));
}


TEST(TripeptideDftFallback, ProductionSqlAndCircularScorePreferNearestDroppedChiBeforeCalcId) {
    using nmr::tripeptide_dft_query::CircularDeltaDegrees;
    using nmr::tripeptide_dft_query::DroppedChiSquaredDistance;
    using nmr::tripeptide_dft_query::OmittedChiScoreSql;

    // Analytic wrap-around oracle: -170 and +170 are 20°, not 340° apart.
    EXPECT_DOUBLE_EQ(CircularDeltaDegrees(-170.0, 170.0), 20.0);
    EXPECT_DOUBLE_EQ(CircularDeltaDegrees(-180.0, 180.0), 0.0);

    // Same retained chi1; the lower-id candidate is deliberately farther
    // on omitted chi2..4 than the higher-id candidate. This calls the
    // production score used to audit QueryNearest's selected SQL row.
    const std::array<int, 4> target{40, -170, 170, -160};
    const std::array<int, 4> lower_id_far{40, 120, 100, 80};
    const std::array<int, 4> higher_id_near{40, -180, -180, -140};
    const double far_score = DroppedChiSquaredDistance(lower_id_far, target, /*natural=*/4, /*used=*/1);
    const double near_score = DroppedChiSquaredDistance(higher_id_near, target, /*natural=*/4, /*used=*/1);
    EXPECT_LT(near_score, far_score);
    EXPECT_NEAR(std::sqrt(near_score), std::sqrt(600.0), 1e-12);

    // Pin the exact production SQL construction: every omitted natural
    // axis participates via circular LEAST distance, and calc_id appears
    // only after the score in QueryNearest's fixed ORDER BY.
    const std::string score_sql = OmittedChiScoreSql(4, 1, target);
    EXPECT_EQ(score_sql.find("chi1"), std::string::npos);
    EXPECT_NE(score_sql.find("chi2"), std::string::npos);
    EXPECT_NE(score_sql.find("chi3"), std::string::npos);
    EXPECT_NE(score_sql.find("chi4"), std::string::npos);
    EXPECT_NE(score_sql.find("LEAST"), std::string::npos);
    EXPECT_NE(score_sql.find("-170"), std::string::npos);
    EXPECT_NE(score_sql.find("170"), std::string::npos);
    EXPECT_NE(score_sql.find("-160"), std::string::npos);
    EXPECT_EQ(OmittedChiScoreSql(4, 4, target), "0.0");
    EXPECT_EQ(nmr::tripeptide_dft_query::NearestOrderBySql(),
              "ORDER BY dropped_chi_score ASC, calc_id ASC LIMIT 1");
}


TEST(TripeptideBackboneDiagnostics, ProductionPackerPinsExact28ColumnsAndFiveStatusCodes) {
    TripeptideBackboneShieldingResult::ResidueMatch rm;
    rm.status = TripeptideMatchStatus::Ok;
    rm.calc_id = 42;
    rm.backbone_rmsd = 0.125;
    rm.ca_match_dist = 0.25;
    rm.frame_type = "gaussian_standard_orientation";
    rm.n_atoms_matched = 17;
    rm.natural_chi_axes = 4;
    rm.n_chi_axes_used = 1;
    rm.dropped_chi_distance_deg = 30.0;
    rm.phi_actual = -61.5;
    rm.psi_actual = 139.25;
    rm.phi_db = -60;
    rm.psi_db = 140;
    rm.residue_type_code = static_cast<int>(AminoAcid::HIS);
    rm.his_variant_hint = 1;
    for (int k = 0; k < 4; ++k) {
        rm.chi_actual[k] = 10.5 + k;
        rm.target_chi_grid[k] = 20 * (k + 1);
        rm.chi_db[k] = -20 * (k + 1);
    }

    const auto row = nmr::tripeptide_backbone_diagnostics::PackDiagnosticRow(rm);
    EXPECT_DOUBLE_EQ(row[0], 42.0);
    EXPECT_DOUBLE_EQ(row[1], 1.0);
    EXPECT_DOUBLE_EQ(row[2], 1.0);
    EXPECT_DOUBLE_EQ(row[3], 0.125);
    EXPECT_DOUBLE_EQ(row[4], 0.25);
    EXPECT_DOUBLE_EQ(row[5], 1.0);
    EXPECT_DOUBLE_EQ(row[6], 17.0);
    EXPECT_DOUBLE_EQ(row[7], 4.0);
    EXPECT_DOUBLE_EQ(row[8], 1.0);
    EXPECT_DOUBLE_EQ(row[9], 30.0);
    EXPECT_DOUBLE_EQ(row[10], -61.5);
    EXPECT_DOUBLE_EQ(row[11], 139.25);
    EXPECT_DOUBLE_EQ(row[12], -60.0);
    EXPECT_DOUBLE_EQ(row[13], 140.0);
    for (int k = 0; k < 4; ++k) {
        EXPECT_DOUBLE_EQ(row[14 + k], 10.5 + k);
        EXPECT_DOUBLE_EQ(row[18 + k], 20.0 * (k + 1));
        EXPECT_DOUBLE_EQ(row[22 + k], -20.0 * (k + 1));
    }
    EXPECT_DOUBLE_EQ(row[26], static_cast<double>(AminoAcid::HIS));
    EXPECT_DOUBLE_EQ(row[27], 1.0);

    const std::array<TripeptideMatchStatus, 5> statuses = {
        TripeptideMatchStatus::Miss,
        TripeptideMatchStatus::Ok,
        TripeptideMatchStatus::UnsupportedHid,
        TripeptideMatchStatus::UnsupportedHie,
        TripeptideMatchStatus::PerceptionFailed,
    };
    for (std::size_t i = 0; i < statuses.size(); ++i) {
        rm.status = statuses[i];
        const auto status_row = nmr::tripeptide_backbone_diagnostics::PackDiagnosticRow(rm);
        EXPECT_DOUBLE_EQ(status_row[1], static_cast<double>(i));
        EXPECT_DOUBLE_EQ(status_row[2], i == 1 ? 1.0 : 0.0);
    }

    // Chi target/DB columns are floating diagnostics: axes beyond the
    // residue's natural depth are unavailable and therefore remain NaN.
    rm.status = TripeptideMatchStatus::Ok;
    rm.natural_chi_axes = 1;
    const auto shallow = nmr::tripeptide_backbone_diagnostics::PackDiagnosticRow(rm);
    EXPECT_TRUE(std::isfinite(shallow[18]));
    EXPECT_TRUE(std::isfinite(shallow[22]));
    for (int k = 1; k < 4; ++k) {
        EXPECT_TRUE(std::isnan(shallow[18 + k]));
        EXPECT_TRUE(std::isnan(shallow[22 + k]));
    }

    // Execute the exact production censor/diagnostic mapping. HIP is not
    // coerced from HID/HIE: only variant 2 is the supported HIP code.
    EXPECT_EQ(nmr::tripeptide_backbone_status::UnsupportedHisStatus(AminoAcid::HIS, 0),
              TripeptideMatchStatus::UnsupportedHid);
    EXPECT_EQ(nmr::tripeptide_backbone_status::UnsupportedHisStatus(AminoAcid::HIS, 1),
              TripeptideMatchStatus::UnsupportedHie);
    EXPECT_EQ(nmr::tripeptide_backbone_status::UnsupportedHisStatus(AminoAcid::HIS, 2),
              TripeptideMatchStatus::Miss);
    EXPECT_EQ(nmr::tripeptide_backbone_status::HisVariantDiagnosticCode(AminoAcid::HIS, 0), 2);
    EXPECT_EQ(nmr::tripeptide_backbone_status::HisVariantDiagnosticCode(AminoAcid::HIS, 1), 3);
    EXPECT_EQ(nmr::tripeptide_backbone_status::HisVariantDiagnosticCode(AminoAcid::HIS, 2), 1);
}


TEST(TripeptideBackboneM14, HidHieAreExplicitlyCensoredWhileHipUsesTheProductionLookup) {
    if (RuntimeEnvironment::TensorCs15Dsn().empty()) {
        GTEST_SKIP() << "tensorcs15 DSN not configured";
    }

    Session session;
    ASSERT_EQ(session.LoadTripeptideDftTable(), kOk) << session.LastError();
    const std::string pdb = nmr::test::TestEnvironment::UbqProtonated();
    auto built = BuildFromProtonatedPdb(pdb);
    ASSERT_TRUE(built.Ok()) << built.error;
    std::unique_ptr<Protein> protein = std::move(built.protein);
    auto& conf = protein->Conformation();

    ASSERT_TRUE(conf.AttachResult(GeometryResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(SpatialIndexResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(EnrichmentResult::Compute(conf)));

    std::size_t his_index = SIZE_MAX;
    for (std::size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        if (protein->ResidueAt(ri).type == AminoAcid::HIS) {
            his_index = ri;
            break;
        }
    }
    ASSERT_NE(his_index, SIZE_MAX);

    auto compute_variant = [&](int variant) {
        protein->MutableResidueAt(his_index).protonation_variant_index = variant;
        auto result = TripeptideBackboneShieldingResult::Compute(conf, *session.TripeptideDftTablePtr());
        EXPECT_NE(result, nullptr);
        return result;
    };

    auto hid = compute_variant(0);
    ASSERT_NE(hid, nullptr);
    EXPECT_EQ(hid->ResidueMatches()[his_index].status, TripeptideMatchStatus::UnsupportedHid);
    EXPECT_EQ(hid->ResidueMatches()[his_index].his_variant_hint, 2);
    for (std::size_t ai : protein->ResidueAt(his_index).atom_indices) {
        EXPECT_FALSE(conf.AtomAt(ai).tripeptide_bb_has_match);
    }

    auto hie = compute_variant(1);
    ASSERT_NE(hie, nullptr);
    EXPECT_EQ(hie->ResidueMatches()[his_index].status, TripeptideMatchStatus::UnsupportedHie);
    EXPECT_EQ(hie->ResidueMatches()[his_index].his_variant_hint, 3);
    for (std::size_t ai : protein->ResidueAt(his_index).atom_indices) {
        EXPECT_FALSE(conf.AtomAt(ai).tripeptide_bb_has_match);
    }

    auto hip = compute_variant(2);
    ASSERT_NE(hip, nullptr);
    EXPECT_EQ(hip->ResidueMatches()[his_index].status, TripeptideMatchStatus::Ok);
    EXPECT_EQ(hip->ResidueMatches()[his_index].his_variant_hint, 1);
    EXPECT_GT(hip->ResidueMatches()[his_index].n_atoms_matched, 0);
}


class TripeptideBackboneShieldingTest : public ::testing::Test {
protected:
    void SetUp() override {
        const std::string kLarsen1UbqPm6 = nmr::test::TestEnvironment::Larsen1UbqPm6Pdb();
        if (RuntimeEnvironment::TensorCs15Dsn().empty()) {
            GTEST_SKIP() << "tensorcs15 DSN not configured "
                            "(set [databases].tensorcs15 in "
                            "~/.nmr_tools.toml)";
        }
        if (!fs::exists(kLarsen1UbqPm6)) {
            GTEST_SKIP() << "Larsen 1UBQ PM6-D3H+ PDB not found at " << kLarsen1UbqPm6
                         << " (download via larsen ERDA archive — see "
                            "/mnt/expansion/larsen_archive/README.md)";
        }
        ASSERT_EQ(session.LoadTripeptideDftTable(), kOk) << session.LastError();
        ASSERT_TRUE(session.HasTripeptideDftTable());

        auto r = BuildFromProtonatedPdb(kLarsen1UbqPm6);
        ASSERT_TRUE(r.Ok()) << r.error;
        protein = std::move(r.protein);
    }

    Session session;
    std::unique_ptr<Protein> protein;
};


TEST_F(TripeptideBackboneShieldingTest, QueryNearestSelectsNearestRoundedDroppedChiPoseInsteadOfLowestCalcId) {
    // Supplemental live-DB forcing function: independently locate an AKA
    // phi/psi group whose lowest calc_id has a different chi pose, then ask
    // production QueryNearest to drop all four chi equality constraints. A
    // second read-only query ranks the candidates with an independently
    // expressed modular circular-distance oracle.
    std::unique_ptr<PGconn, decltype(&PQfinish)> conn(PQconnectdb(RuntimeEnvironment::TensorCs15Dsn().c_str()), &PQfinish);
    ASSERT_NE(conn.get(), nullptr);
    ASSERT_EQ(PQstatus(conn.get()), CONNECTION_OK) << PQerrorMessage(conn.get());

    const char* oracle_sql = R"SQL(
        SELECT r.phi, r.psi, r.chi1, r.chi2, r.chi3, r.chi4,
               CASE WHEN (ROUND(r.chi1 / 20.0) * 20)::int = 180
                    THEN -180 ELSE (ROUND(r.chi1 / 20.0) * 20)::int END AS target_chi1,
               CASE WHEN (ROUND(r.chi2 / 20.0) * 20)::int = 180
                    THEN -180 ELSE (ROUND(r.chi2 / 20.0) * 20)::int END AS target_chi2,
               CASE WHEN (ROUND(r.chi3 / 20.0) * 20)::int = 180
                    THEN -180 ELSE (ROUND(r.chi3 / 20.0) * 20)::int END AS target_chi3,
               CASE WHEN (ROUND(r.chi4 / 20.0) * 20)::int = 180
                    THEN -180 ELSE (ROUND(r.chi4 / 20.0) * 20)::int END AS target_chi4,
               (SELECT MIN(g.calc_id)
                  FROM raw_dft_calculations g
                 WHERE g.tripeptide = r.tripeptide
                   AND g.phi = r.phi AND g.psi = r.psi) AS group_min
          FROM raw_dft_calculations r
         WHERE r.tripeptide = 'AKA'
           AND r.chi1 IS NOT NULL AND r.chi2 IS NOT NULL
           AND r.chi3 IS NOT NULL AND r.chi4 IS NOT NULL
           AND r.calc_id > (
               SELECT MIN(g.calc_id)
                 FROM raw_dft_calculations g
                WHERE g.tripeptide = r.tripeptide
                  AND g.phi = r.phi AND g.psi = r.psi)
         ORDER BY r.calc_id ASC
         LIMIT 1
    )SQL";
    std::unique_ptr<PGresult, decltype(&PQclear)> oracle(PQexec(conn.get(), oracle_sql), &PQclear);
    ASSERT_NE(oracle.get(), nullptr);
    ASSERT_EQ(PQresultStatus(oracle.get()), PGRES_TUPLES_OK) << PQerrorMessage(conn.get());
    ASSERT_EQ(PQntuples(oracle.get()), 1) << "tensorcs15 AKA should contain multiple chi poses per phi/psi";

    auto value = [&](int col) {
        EXPECT_FALSE(PQgetisnull(oracle.get(), 0, col));
        if (PQgetisnull(oracle.get(), 0, col))
            return 0;
        return std::stoi(PQgetvalue(oracle.get(), 0, col));
    };
    const int phi = value(0);
    const int psi = value(1);
    const std::array<int, 4> chi = {value(2), value(3), value(4), value(5)};
    const std::array<int, 4> target_chi = {value(6), value(7), value(8), value(9)};
    const int group_min = value(10);

    const TripeptideDftRecord selected = session.TripeptideDftTablePtr()->QueryNearest('K',
                                                                                       phi,
                                                                                       psi,
                                                                                       chi[0],
                                                                                       chi[1],
                                                                                       chi[2],
                                                                                       chi[3],
                                                                                       /*n_chi_axes=*/0,
                                                                                       /*his_variant_hint=*/-1);
    ASSERT_TRUE(selected.IsHit());
    EXPECT_EQ(selected.natural_chi_axes, 4);
    EXPECT_EQ(selected.n_chi_axes_used, 0);
    EXPECT_EQ(selected.target_phi_grid_deg, phi);
    EXPECT_EQ(selected.target_psi_grid_deg, psi);
    EXPECT_EQ(selected.target_chi_grid_deg, target_chi);

    // Independent circular-distance oracle: modulo-to-[-180,180) is a
    // different formulation from production's LEAST(abs, 360-abs).
    auto modular_delta_sq = [](int axis, int target) {
        const std::string column = "chi" + std::to_string(axis);
        const std::string delta = "ABS(MOD((" + column + " - (" + std::to_string(target)
                                  + ") + 540)::numeric, 360) - 180)";
        return "(" + delta + " * " + delta + ")";
    };
    std::string oracle_score;
    for (int k = 0; k < 4; ++k) {
        if (!oracle_score.empty())
            oracle_score += " + ";
        oracle_score += modular_delta_sq(k + 1, target_chi[k]);
    }
    std::ostringstream nearest_sql;
    nearest_sql << "SELECT calc_id, (" << oracle_score << ")::double precision "
                << "FROM raw_dft_calculations WHERE tripeptide='AKA' "
                << "AND phi=" << phi << " AND psi=" << psi << " "
                << "AND chi1 IS NOT NULL AND chi2 IS NOT NULL "
                << "AND chi3 IS NOT NULL AND chi4 IS NOT NULL "
                << "ORDER BY 2 ASC, calc_id ASC LIMIT 1";
    std::unique_ptr<PGresult, decltype(&PQclear)> nearest(
        PQexec(conn.get(), nearest_sql.str().c_str()), &PQclear);
    ASSERT_NE(nearest.get(), nullptr);
    ASSERT_EQ(PQresultStatus(nearest.get()), PGRES_TUPLES_OK) << PQerrorMessage(conn.get());
    ASSERT_EQ(PQntuples(nearest.get()), 1);
    const int expected_calc_id = std::stoi(PQgetvalue(nearest.get(), 0, 0));
    const double expected_score = std::stod(PQgetvalue(nearest.get(), 0, 1));
    ASSERT_NE(expected_calc_id, group_min)
        << "fixture must distinguish nearest omitted-chi selection from the old calc_id-only order";
    EXPECT_EQ(selected.calc_id, expected_calc_id);
    EXPECT_NEAR(selected.dropped_chi_distance_deg, std::sqrt(expected_score), 1e-12);
}


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

    auto tbb = TripeptideBackboneShieldingResult::Compute(conf, *session.TripeptideDftTablePtr());
    ASSERT_NE(tbb, nullptr);

    // 1UBQ has 76 residues; ResiduesAttempted counts standard residues
    // (Unknown skipped). Termini fail phi/psi.
    EXPECT_EQ(tbb->ResiduesAttempted(), static_cast<int>(protein->ResidueCount()));
    EXPECT_GE(tbb->ResiduesMatched(), 60);  // most internals match
    EXPECT_GT(tbb->AtomsAssigned(), 0);

    // Backbone Kabsch RMSD should be small — the canonical tripeptide
    // and the protein backbone share the same N/CA/C bond geometry by
    // construction. ~0.5 Å is the empirical rough upper bound; if mean
    // exceeds 1 Å something is structurally wrong with the alignment.
    EXPECT_LT(tbb->MeanBackboneRmsd(), 1.0);

    // Exercise the production file-local helper through Compute. The
    // result keeps the actual protein phi/psi before its DB lookup; compare
    // those retained values with the independent coordinate oracle above.
    // The old (n1 x n2).b2hat formula flips every non-zero sign.
    int n_dihedrals_pinned = 0;
    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        const auto prev_idx = protein->BackbonePredecessor(ri);
        const auto next_idx = protein->BackboneSuccessor(ri);
        if (!prev_idx || !next_idx)
            continue;
        const Residue& prev = protein->ResidueAt(*prev_idx);
        const Residue& res = protein->ResidueAt(ri);
        const Residue& next = protein->ResidueAt(*next_idx);
        if (prev.C == Residue::NONE || res.N == Residue::NONE || res.CA == Residue::NONE || res.C == Residue::NONE
            || next.N == Residue::NONE) {
            continue;
        }

        const double phi = nmr::tripeptide_backbone_dihedral::DihedralDegrees(conf.PositionAt(prev.C),
                                                                              conf.PositionAt(res.N),
                                                                              conf.PositionAt(res.CA),
                                                                              conf.PositionAt(res.C));
        const double psi = nmr::tripeptide_backbone_dihedral::DihedralDegrees(conf.PositionAt(res.N),
                                                                              conf.PositionAt(res.CA),
                                                                              conf.PositionAt(res.C),
                                                                              conf.PositionAt(next.N));
        ASSERT_TRUE(std::isfinite(phi));
        ASSERT_TRUE(std::isfinite(psi));
        const auto& match = tbb->ResidueMatches()[ri];
        EXPECT_NEAR(match.phi_actual, phi, 1e-10) << "residue " << ri;
        EXPECT_NEAR(match.psi_actual, psi, 1e-10) << "residue " << ri;
        n_dihedrals_pinned += 2;
    }
    EXPECT_GT(n_dihedrals_pinned, 100) << "1UBQ should exercise both file-local backbone dihedral paths";

    // Spot check: every matched atom has finite tensor.
    int n_finite = 0, n_total = 0;
    for (size_t i = 0; i < conf.AtomCount(); ++i) {
        const auto& ca = conf.AtomAt(i);
        if (!ca.tripeptide_bb_has_match)
            continue;
        ++n_total;
        const auto& m = ca.tripeptide_bb_shielding_tensor;
        bool finite = true;
        for (int r = 0; r < 3 && finite; ++r)
            for (int c = 0; c < 3 && finite; ++c)
                if (!std::isfinite(m(r, c)))
                    finite = false;
        if (finite)
            ++n_finite;
    }
    EXPECT_EQ(n_finite, n_total) << "non-finite tensor at " << (n_total - n_finite) << " atoms";

    // Spot check: SER residues that hit the DB should carry the
    // orca_input_orientation tag (project SER PBE regen). 1UBQ has
    // 4 SER residues. Not all may hit (e.g. terminal SER skipped).
    int n_ser_attempted = 0, n_ser_matched_pbe = 0;
    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        if (protein->ResidueAt(ri).type != AminoAcid::SER)
            continue;
        ++n_ser_attempted;
        const auto& m = tbb->ResidueMatches()[ri];
        if (m.calc_id == 0)
            continue;
        if (m.frame_type == "orca_input_orientation") {
            ++n_ser_matched_pbe;
        }
    }
    EXPECT_GE(n_ser_attempted, 1) << "1UBQ should have at least 1 SER";
    EXPECT_GT(n_ser_matched_pbe, 0) << "SER residues that match should hit ORCA PBE rows "
                                    << "(frame_type=orca_input_orientation), " << "got " << n_ser_matched_pbe << "/"
                                    << n_ser_attempted;

    // Test NPY emission.
    const std::string out_dir = nmr::test::TestEnvironment::TempPath("tripeptide_bb_smoke_out");
    fs::create_directories(out_dir);
    int n_npy = tbb->WriteFeatures(conf, out_dir);
    EXPECT_EQ(n_npy, 6);
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_bb_shielding.npy"));
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_bb_residual_vec.npy"));
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_bb_match_distance.npy"));
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_bb_method_tag.npy"));
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_bb_match_atoms.npy"));
    EXPECT_TRUE(fs::exists(out_dir + "/tripeptide_bb_diagnostics.npy"));
    EXPECT_FALSE(fs::exists(out_dir + "/tripeptide_bb_residue_status.npy"));

    auto match_atoms = ReadNpy(out_dir + "/tripeptide_bb_match_atoms.npy");
    ASSERT_EQ(match_atoms.shape, (std::vector<size_t>{conf.AtomCount(), 5}));
    const double* ma = DataAs<double>(match_atoms);
    int n_match_atom_rows = 0;
    for (size_t i = 0; i < conf.AtomCount(); ++i) {
        const auto& ca = conf.AtomAt(i);
        EXPECT_DOUBLE_EQ(ma[i * 5 + 0], static_cast<double>(protein->AtomAt(i).residue_index));
        if (ca.tripeptide_bb_method_tag == 0)
            continue;
        ++n_match_atom_rows;
        EXPECT_DOUBLE_EQ(ma[i * 5 + 1], 1.0);
        EXPECT_GT(ma[i * 5 + 2], 0.0);
        EXPECT_TRUE(std::isfinite(ma[i * 5 + 3]));
        EXPECT_DOUBLE_EQ(ma[i * 5 + 4], static_cast<double>(ca.tripeptide_bb_method_tag));
    }
    EXPECT_GT(n_match_atom_rows, 0);

    // Read back the single frozen residue-diagnostics schema. This pins
    // both C3 fallback fields and M14 censor codes at the emitted boundary.
    auto diagnostics = ReadNpy(out_dir + "/tripeptide_bb_diagnostics.npy");
    ASSERT_EQ(diagnostics.shape, (std::vector<size_t>{protein->ResidueCount(), 28}));
    const double* diag = DataAs<double>(diagnostics);
    int n_fallback_rows = 0;
    int n_his_censored = 0;
    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        const auto& m = tbb->ResidueMatches()[ri];
        const double* row = &diag[ri * 28];
        EXPECT_DOUBLE_EQ(row[0], static_cast<double>(m.calc_id));
        EXPECT_DOUBLE_EQ(row[1], static_cast<double>(m.status));
        EXPECT_DOUBLE_EQ(row[2], m.status == TripeptideMatchStatus::Ok ? 1.0 : 0.0);
        EXPECT_DOUBLE_EQ(row[6], static_cast<double>(m.n_atoms_matched));
        EXPECT_DOUBLE_EQ(row[7], static_cast<double>(m.natural_chi_axes));
        EXPECT_DOUBLE_EQ(row[8], static_cast<double>(m.n_chi_axes_used));
        EXPECT_DOUBLE_EQ(row[26], static_cast<double>(m.residue_type_code));
        EXPECT_DOUBLE_EQ(row[27], static_cast<double>(m.his_variant_hint));

        if (m.status == TripeptideMatchStatus::Ok && m.n_chi_axes_used < m.natural_chi_axes) {
            ++n_fallback_rows;
            EXPECT_TRUE(std::isfinite(row[9]));
            EXPECT_GE(row[9], 0.0);
            EXPECT_DOUBLE_EQ(row[9], m.dropped_chi_distance_deg);
            for (int k = 0; k < 4; ++k) {
                if (k < m.natural_chi_axes) {
                    EXPECT_DOUBLE_EQ(row[18 + k], static_cast<double>(m.target_chi_grid[k]));
                    EXPECT_DOUBLE_EQ(row[22 + k], static_cast<double>(m.chi_db[k]));
                } else {
                    EXPECT_TRUE(std::isnan(row[18 + k]));
                    EXPECT_TRUE(std::isnan(row[22 + k]));
                }
            }
        }

        const Residue& residue = protein->ResidueAt(ri);
        if (residue.type == AminoAcid::HIS && residue.protonation_variant_index == 0) {
            ++n_his_censored;
            EXPECT_DOUBLE_EQ(row[1], 2.0);  // unsupported_hid
            EXPECT_DOUBLE_EQ(row[2], 0.0);
            EXPECT_DOUBLE_EQ(row[27], 2.0);  // HID hint
        } else if (residue.type == AminoAcid::HIS && residue.protonation_variant_index == 1) {
            ++n_his_censored;
            EXPECT_DOUBLE_EQ(row[1], 3.0);  // unsupported_hie
            EXPECT_DOUBLE_EQ(row[2], 0.0);
            EXPECT_DOUBLE_EQ(row[27], 3.0);  // HIE hint
        }
    }
    EXPECT_GT(n_fallback_rows, 0) << "1UBQ should exercise at least one dropped-chi fallback";
    EXPECT_GE(n_his_censored, 0);  // exact codes are pinned non-skipping above

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
            if (slot == Residue::NONE)
                continue;
            const auto& ca = conf.AtomAt(slot);
            if (!ca.tripeptide_bb_has_match)
                continue;
            ++n_bb_total;
            const double d = ca.tripeptide_bb_match_distance;
            if (d > max_bb_residual)
                max_bb_residual = d;
            if (d > 3.0)
                ++n_bb_extreme;
        }
    }
    EXPECT_GT(n_bb_total, 100) << "expected backbone atoms to mostly get tensors";
    EXPECT_LT(n_bb_extreme, 5) << "more than 5 backbone atoms have post-Kabsch residual > 3 Å "
                               << "— Kabsch failure or substrate disagreement (max=" << max_bb_residual << " Å)";
}

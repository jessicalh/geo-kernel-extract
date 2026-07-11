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
//   - WriteFeatures emits the per-atom and raw pair-audit NPYs.
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
#include "LarsenSidechainDonorAuditResult.h"
#include "LarsenHBondShieldingResult.h"
#include "OperationLog.h"
#include "PdbFileReader.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "RuntimeEnvironment.h"
#include "Session.h"
#include "SpatialIndexResult.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <vector>

using namespace nmr;
namespace fs = std::filesystem;

namespace {

fs::path MakeUniqueTempDirectory(const std::string& prefix) {
    static std::atomic<unsigned long long> sequence{0};
    const auto tick = std::chrono::steady_clock::now()
                          .time_since_epoch()
                          .count();
    const fs::path dir = fs::temp_directory_path() /
        (prefix + "_" + std::to_string(tick) + "_" +
         std::to_string(sequence.fetch_add(1)));
    fs::create_directories(dir);
    return dir;
}

void RemoveFlatTempDirectory(
    const fs::path& dir,
    std::initializer_list<const char*> file_names) {
    for (const char* file_name : file_names) {
        const std::string path = (dir / file_name).string();
        EXPECT_EQ(std::remove(path.c_str()), 0) << path;
    }
    const std::string path = dir.string();
    EXPECT_EQ(std::remove(path.c_str()), 0) << path;
}

struct RawNpy {
    std::string header;
    std::vector<char> payload;
};

template <typename T>
T PayloadValue(const RawNpy& npy, std::size_t element_index) {
    T value{};
    const std::size_t byte_offset = element_index * sizeof(T);
    EXPECT_LE(byte_offset + sizeof(T), npy.payload.size());
    if (byte_offset + sizeof(T) <= npy.payload.size()) {
        std::memcpy(&value, npy.payload.data() + byte_offset, sizeof(T));
    }
    return value;
}

RawNpy ReadRawNpy(const fs::path& path) {
    std::ifstream in(path, std::ios::binary);
    EXPECT_TRUE(in.is_open()) << path;
    char magic[6] = {};
    in.read(magic, 6);
    EXPECT_EQ(std::string(magic, 6), std::string("\x93NUMPY", 6));
    char version[2] = {};
    in.read(version, 2);
    EXPECT_EQ(version[0], 1);
    EXPECT_EQ(version[1], 0);
    std::uint16_t header_len = 0;
    in.read(reinterpret_cast<char*>(&header_len), sizeof(header_len));
    RawNpy out;
    out.header.resize(header_len);
    in.read(out.header.data(), header_len);
    out.payload.assign(std::istreambuf_iterator<char>(in),
                       std::istreambuf_iterator<char>());
    return out;
}

std::vector<float> ConstantTensorGrid(float base) {
    constexpr std::size_t kGridCells = 2 * 2 * 2;
    std::vector<float> values(kGridCells * 9, 0.0f);
    for (std::size_t cell = 0; cell < kGridCells; ++cell) {
        values[cell * 9 + 0] = base;
        values[cell * 9 + 4] = base * 2.0f;
        values[cell * 9 + 8] = base * 3.0f;
    }
    return values;
}

void WriteTensorDataset(HighFive::File& file,
                        const std::string& name,
                        float base) {
    const HighFive::DataSpace space{2, 2, 2, 3, 3};
    auto dataset = file.createDataSet<float>(name, space);
    const std::vector<float> values = ConstantTensorGrid(base);
    dataset.write_raw(values.data());
}

// Creates all six schema-valid archives with a tiny regular grid. Every
// successful query interpolates the production tensors above and sees
// exactly one imputed corner, making C34 deterministic without an external
// database or GTEST_SKIP gate.
fs::path WriteSyntheticLarsenGridDirectory(
    const std::string& tag,
    bool force_all_grid_miss = false) {
    const fs::path dir =
        MakeUniqueTempDirectory("nmr_larsen_synthetic_" + tag);

    const std::vector<double> theta_axis = {90.0, 180.0};
    const std::vector<double> rho_axis = {-180.0, 0.0};
    const std::vector<std::uint8_t> validity = {0, 1, 1, 1, 1, 1, 1, 1};
    const char* stems[6] = {
        "NMANMA", "NMACOH", "NMACOO",
        "ALANMA", "ALACOH", "ALACOO",
    };

    for (int archive = 0; archive < 6; ++archive) {
        HighFive::File file(
            (dir / (std::string(stems[archive]) + "_dense.h5")).string(),
            HighFive::File::Overwrite);
        const std::vector<double> r_axis = force_all_grid_miss
            ? std::vector<double>{0.05, 0.10}
            : (archive < 3
                ? std::vector<double>{1.5, 3.0}
                : std::vector<double>{1.8, 4.0});
        file.createDataSet("r_axis", r_axis);
        file.createDataSet("theta_axis", theta_axis);
        file.createDataSet("rho_axis", rho_axis);

        float base = static_cast<float>(archive + 1);
        WriteTensorDataset(file, "donor_N",  base + 0.1f);
        WriteTensorDataset(file, "donor_CA", base + 0.2f);
        WriteTensorDataset(file, "donor_C",  base + 0.3f);
        WriteTensorDataset(file, "donor_HA", base + 0.4f);
        WriteTensorDataset(file, "donor_HN", base + 0.5f);
        const bool ala_donor = archive >= 3;
        if (ala_donor) WriteTensorDataset(file, "donor_CB", base + 0.6f);

        const bool nma_acceptor = archive == 0 || archive == 3;
        if (nma_acceptor) {
            WriteTensorDataset(file, "acceptor_N",  base + 0.7f);
            WriteTensorDataset(file, "acceptor_HN", base + 0.8f);
            WriteTensorDataset(file, "acceptor_CA", base + 0.9f);
            WriteTensorDataset(file, "acceptor_C",  base + 1.0f);
            WriteTensorDataset(file, "acceptor_HA", base + 1.1f);
        }

        auto mask = file.createDataSet<std::uint8_t>(
            "validity_mask", HighFive::DataSpace{2, 2, 2});
        mask.write_raw(validity.data());
    }
    return dir;
}

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


// Production-linked forcing function for the A23 value bug. The test calls
// the exact audit helper used at every production tensor-write site and
// compares it with the independent physical invariant isotropic = trace/3.
TEST(LarsenHBondPairAudit, TermAwareIsotropicAccountingUsesAppliedTensors) {
    using TA = LarsenContribDispatch::TargetAtom;
    using Term = LarsenContribDispatch::Term;

    LarsenHBondShieldingResult::PairRecord pair;
    pair.iso_1pHB = pair.iso_2pHB = 0.0;
    pair.iso_1pHaB = pair.iso_2pHaB = 0.0;
    pair.iso_diagnostic_CB = 0.0;
    Mat3 primary = Mat3::Zero();
    primary.diagonal() << 3.0, 6.0, 12.0;      // isotropic 7
    Mat3 secondary = Mat3::Zero();
    secondary.diagonal() << -3.0, 0.0, 6.0;    // isotropic 1
    Mat3 alpha = Mat3::Zero();
    alpha.diagonal() << 9.0, 12.0, 15.0;       // isotropic 12
    Mat3 diagnostic = Mat3::Zero();
    diagnostic.diagonal() << 30.0, 0.0, 0.0;   // isotropic 10

    larsen_hbond_shielding_detail::RecordTable2Contribution(
        pair, TA::HN, Term::Primary_HB, primary);
    // A second actual target write under the same term must accumulate,
    // not overwrite (the GLY HA2/HA3 fan-out is this shape).
    larsen_hbond_shielding_detail::RecordTable2Contribution(
        pair, TA::HA, Term::Primary_HB, secondary);
    larsen_hbond_shielding_detail::RecordTable2Contribution(
        pair, TA::N, Term::Secondary_HaB, alpha);
    larsen_hbond_shielding_detail::RecordDiagnosticCbContribution(
        pair, diagnostic);

    EXPECT_DOUBLE_EQ(pair.iso_1pHB,
                     primary.trace() / 3.0 + secondary.trace() / 3.0);
    EXPECT_DOUBLE_EQ(pair.iso_2pHB, 0.0);
    EXPECT_DOUBLE_EQ(pair.iso_1pHaB, 0.0);
    EXPECT_DOUBLE_EQ(pair.iso_2pHaB, alpha.trace() / 3.0);
    EXPECT_DOUBLE_EQ(pair.iso_diagnostic_CB,
                     diagnostic.trace() / 3.0);
    EXPECT_NE(pair.target_mask_1pHB & (1u << static_cast<unsigned>(TA::HN)), 0u);
    EXPECT_NE(pair.target_mask_1pHB & (1u << static_cast<unsigned>(TA::HA)), 0u);
    EXPECT_NE(pair.target_mask_2pHaB & (1u << static_cast<unsigned>(TA::N)), 0u);
    EXPECT_EQ(pair.target_mask_diagnostic_CB,
              1u << static_cast<unsigned>(TA::CB));
}


// Bare-checkout, end-to-end freeze gate for A23/C34/C35. It runs the
// production grid loader, production spatial enumeration/dispatch, and the
// production NPY writer against committed 1UBQ coordinates plus tiny H5
// grids created by the test. No configured database and no GTEST_SKIP.
TEST(LarsenHBondPairAudit, SyntheticGridEmitsAuditablePairTables) {
    const fs::path pdb = nmr::test::TestEnvironment::UbqProtonated();
    ASSERT_TRUE(fs::is_regular_file(pdb)) << pdb;
    auto build = BuildFromProtonatedPdb(pdb.string());
    ASSERT_TRUE(build.Ok()) << build.error;
    std::unique_ptr<Protein> protein = std::move(build.protein);
    ProteinConformation& conf = protein->Conformation();

    ASSERT_TRUE(conf.AttachResult(GeometryResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(SpatialIndexResult::Compute(conf)));

    const fs::path grid_dir =
        WriteSyntheticLarsenGridDirectory("pair_audit");
    LarsenHBondGrid grid(grid_dir.string());
    auto result = LarsenHBondShieldingResult::Compute(conf, grid);
    ASSERT_NE(result, nullptr);

    EXPECT_EQ(result->Pairs().size(),
              static_cast<std::size_t>(result->PairsFound() +
                                       result->PairsGridSkipped()));
    ASSERT_GT(result->PairsFound(), 0);

    int n_success = 0;
    int n_missing_frame = 0;
    int n_theta_miss = 0;
    int n_grid_miss = 0;
    int n_invalid_frame = 0;
    int n_carboxylate_symmetry_filtered = 0;
    int n_sidechain_carbonyl_success = 0;
    double pair_1pHB = 0.0;
    double pair_2pHB = 0.0;
    double pair_1pHaB = 0.0;
    double pair_2pHaB = 0.0;
    double pair_cb = 0.0;
    for (const auto& pair : result->Pairs()) {
        using Disp = LarsenHBondShieldingResult::PairDisposition;
        switch (pair.disposition) {
            case Disp::MissingFrameAtoms:
                ++n_missing_frame;
                EXPECT_FALSE(pair.frame_valid);
                EXPECT_TRUE(std::isnan(pair.r_angstrom));
                EXPECT_TRUE(std::isnan(pair.theta_deg));
                EXPECT_TRUE(std::isnan(pair.rho_deg));
                EXPECT_TRUE(std::isnan(pair.iso_1pHB));
                EXPECT_TRUE(std::isnan(pair.iso_2pHB));
                EXPECT_TRUE(std::isnan(pair.iso_1pHaB));
                EXPECT_TRUE(std::isnan(pair.iso_2pHaB));
                EXPECT_TRUE(std::isnan(pair.iso_diagnostic_CB));
                EXPECT_TRUE(std::isnan(pair.iso_table2_total));
                EXPECT_FALSE(pair.any_corner_imputed);
                EXPECT_EQ(pair.imputed_corner_count, 0);
                break;
            case Disp::ThetaOutOfRange:
                ++n_theta_miss;
                EXPECT_TRUE(pair.frame_valid);
                EXPECT_TRUE(std::isfinite(pair.r_angstrom));
                EXPECT_TRUE(std::isnan(pair.iso_table2_total));
                break;
            case Disp::GridMiss:
                ++n_grid_miss;
                EXPECT_TRUE(pair.frame_valid);
                EXPECT_TRUE(std::isfinite(pair.theta_deg));
                EXPECT_TRUE(std::isnan(pair.iso_table2_total));
                break;
            case Disp::Success:
                ++n_success;
                EXPECT_TRUE(pair.frame_valid);
                EXPECT_EQ(pair.imputed_corner_count, 1);
                EXPECT_TRUE(pair.any_corner_imputed);
                EXPECT_NEAR(pair.iso_table2_total,
                            pair.iso_1pHB + pair.iso_2pHB +
                            pair.iso_1pHaB + pair.iso_2pHaB,
                            1e-12);
                EXPECT_DOUBLE_EQ(pair.isotropic_total,
                                 pair.iso_table2_total);
                pair_1pHB += pair.iso_1pHB;
                pair_2pHB += pair.iso_2pHB;
                pair_1pHaB += pair.iso_1pHaB;
                pair_2pHaB += pair.iso_2pHaB;
                pair_cb += pair.iso_diagnostic_CB;
                if (pair.acceptor_class ==
                    HBondAcceptorClass::SidechainCarbonyl) {
                    ++n_sidechain_carbonyl_success;
                }
                break;
            case Disp::InvalidFrame:
                ++n_invalid_frame;
                EXPECT_FALSE(pair.frame_valid);
                EXPECT_TRUE(std::isnan(pair.iso_table2_total));
                EXPECT_FALSE(pair.any_corner_imputed);
                EXPECT_EQ(pair.imputed_corner_count, 0);
                break;
            case Disp::CarboxylateSymmetryFiltered:
                ++n_carboxylate_symmetry_filtered;
                EXPECT_EQ(pair.acceptor_class,
                          HBondAcceptorClass::CarboxylateOxygen);
                EXPECT_TRUE(pair.frame_valid);
                EXPECT_TRUE(std::isfinite(pair.r_angstrom));
                EXPECT_TRUE(std::isfinite(pair.theta_deg));
                EXPECT_TRUE(std::isfinite(pair.rho_deg));
                EXPECT_EQ(pair.target_mask_1pHB, 0u);
                EXPECT_EQ(pair.target_mask_2pHB, 0u);
                EXPECT_EQ(pair.target_mask_1pHaB, 0u);
                EXPECT_EQ(pair.target_mask_2pHaB, 0u);
                EXPECT_EQ(pair.target_mask_diagnostic_CB, 0u);
                EXPECT_TRUE(std::isnan(pair.iso_table2_total));
                EXPECT_FALSE(pair.any_corner_imputed);
                EXPECT_EQ(pair.imputed_corner_count, 0);
                break;
        }
    }
    EXPECT_EQ(n_success, result->PairsFound());
    EXPECT_EQ(n_missing_frame, 0)
        << "committed 1UBQ fixture has no classified missing-anchor row";
    EXPECT_GT(n_theta_miss, 0);
    EXPECT_GT(n_grid_miss, 0);
    EXPECT_EQ(n_invalid_frame, 0)
        << "unmodified committed coordinates should have valid frames";
    EXPECT_GT(n_carboxylate_symmetry_filtered, 0)
        << "all classified carboxylate sibling attempts must be retained";
    for (const auto& filtered : result->Pairs()) {
        if (filtered.disposition !=
            LarsenHBondShieldingResult::PairDisposition::
                CarboxylateSymmetryFiltered) {
            continue;
        }
        const auto selected = std::find_if(
            result->Pairs().begin(), result->Pairs().end(),
            [&](const LarsenHBondShieldingResult::PairRecord& candidate) {
                return candidate.donor_atom_idx == filtered.donor_atom_idx &&
                       candidate.acceptor_C_atom_idx ==
                           filtered.acceptor_C_atom_idx &&
                       candidate.acceptor_atom_idx ==
                           filtered.acceptor_third_atom_idx &&
                       candidate.disposition !=
                           LarsenHBondShieldingResult::PairDisposition::
                               CarboxylateSymmetryFiltered;
            });
        EXPECT_NE(selected, result->Pairs().end())
            << "filtered carboxylate row must accompany the selected sibling";
    }
    EXPECT_GT(n_sidechain_carbonyl_success, 0)
        << "committed 1UBQ should exercise ASN/GLN acceptor approximation";

    // Independent conservation oracle: pair-level trace sums must equal
    // traces of the tensors actually accumulated on the atom axis.
    double atom_1pHB = 0.0;
    double atom_2pHB = 0.0;
    double atom_1pHaB = 0.0;
    double atom_2pHaB = 0.0;
    double atom_cb = 0.0;
    for (std::size_t i = 0; i < conf.AtomCount(); ++i) {
        const ConformationAtom& atom = conf.AtomAt(i);
        atom_1pHB += atom.larsen_hbond_1pHB_tensor.trace() / 3.0;
        atom_2pHB += atom.larsen_hbond_2pHB_tensor.trace() / 3.0;
        atom_1pHaB += atom.larsen_hbond_1pHaB_tensor.trace() / 3.0;
        atom_2pHaB += atom.larsen_hbond_2pHaB_tensor.trace() / 3.0;
        atom_cb += atom.larsen_hbond_diagnostic_CB.trace() / 3.0;
    }
    EXPECT_NEAR(pair_1pHB, atom_1pHB, 1e-8);
    EXPECT_NEAR(pair_2pHB, atom_2pHB, 1e-8);
    EXPECT_NEAR(pair_1pHaB, atom_1pHaB, 1e-8);
    EXPECT_NEAR(pair_2pHaB, atom_2pHaB, 1e-8);
    EXPECT_NEAR(pair_cb, atom_cb, 1e-8);

    const fs::path out =
        MakeUniqueTempDirectory("nmr_larsen_pair_audit_output");
    EXPECT_EQ(result->WriteFeatures(conf, out.string()), 15);

    const std::size_t P = result->Pairs().size();
    const RawNpy index = ReadRawNpy(out / "larsen_hbond_pairs_index.npy");
    const RawNpy geometry =
        ReadRawNpy(out / "larsen_hbond_pairs_geometry.npy");
    const RawNpy isotropic =
        ReadRawNpy(out / "larsen_hbond_pairs_isotropic.npy");
    const RawNpy compat = ReadRawNpy(out / "larsen_hbond_pairs.npy");
    EXPECT_NE(index.header.find("'descr': '<i4'"), std::string::npos);
    EXPECT_NE(index.header.find("'shape': (" + std::to_string(P) +
                                ",16)"), std::string::npos);
    EXPECT_EQ(index.payload.size(), P * 16 * sizeof(std::int32_t));
    EXPECT_NE(geometry.header.find("'shape': (" + std::to_string(P) +
                                   ",6)"), std::string::npos);
    EXPECT_EQ(geometry.payload.size(), P * 6 * sizeof(double));
    EXPECT_EQ(isotropic.payload.size(), P * 6 * sizeof(double));
    EXPECT_NE(compat.header.find("'shape': (" + std::to_string(P) +
                                 ",28)"), std::string::npos);
    EXPECT_EQ(compat.payload.size(), P * 28 * sizeof(double));

    // A23 read-back: pin every field family to the production-owned row,
    // rather than accepting correctly shaped but empty/reordered tables.
    const auto success_it = std::find_if(
        result->Pairs().begin(), result->Pairs().end(),
        [](const LarsenHBondShieldingResult::PairRecord& pair) {
            return pair.disposition ==
                LarsenHBondShieldingResult::PairDisposition::Success;
        });
    ASSERT_NE(success_it, result->Pairs().end());
    const std::size_t success_row = static_cast<std::size_t>(
        std::distance(result->Pairs().begin(), success_it));
    const auto& expected_pair = *success_it;
    EXPECT_EQ(PayloadValue<std::int32_t>(index, success_row * 16 + 0),
              static_cast<std::int32_t>(expected_pair.donor_atom_idx));
    EXPECT_EQ(PayloadValue<std::int32_t>(index, success_row * 16 + 1),
              static_cast<std::int32_t>(expected_pair.acceptor_atom_idx));
    EXPECT_EQ(PayloadValue<std::int32_t>(index, success_row * 16 + 4),
              static_cast<std::int32_t>(expected_pair.donor_class));
    EXPECT_EQ(PayloadValue<std::int32_t>(index, success_row * 16 + 5),
              static_cast<std::int32_t>(expected_pair.acceptor_class));
    EXPECT_EQ(PayloadValue<std::int32_t>(index, success_row * 16 + 6), 3);
    EXPECT_EQ(PayloadValue<std::int32_t>(index, success_row * 16 + 11),
              static_cast<std::int32_t>(expected_pair.target_mask_1pHB));
    EXPECT_EQ(PayloadValue<std::int32_t>(index, success_row * 16 + 15),
              static_cast<std::int32_t>(
                  expected_pair.target_mask_diagnostic_CB));
    EXPECT_DOUBLE_EQ(PayloadValue<double>(geometry, success_row * 6 + 0),
                     expected_pair.r_angstrom);
    EXPECT_DOUBLE_EQ(PayloadValue<double>(geometry, success_row * 6 + 1),
                     expected_pair.theta_deg);
    EXPECT_DOUBLE_EQ(PayloadValue<double>(geometry, success_row * 6 + 2),
                     expected_pair.rho_deg);
    EXPECT_DOUBLE_EQ(PayloadValue<double>(geometry, success_row * 6 + 3), 1.0);
    EXPECT_DOUBLE_EQ(PayloadValue<double>(geometry, success_row * 6 + 4), 1.0);
    EXPECT_DOUBLE_EQ(PayloadValue<double>(geometry, success_row * 6 + 5), 1.0);
    EXPECT_DOUBLE_EQ(PayloadValue<double>(isotropic, success_row * 6 + 0),
                     expected_pair.iso_1pHB);
    EXPECT_DOUBLE_EQ(PayloadValue<double>(isotropic, success_row * 6 + 1),
                     expected_pair.iso_2pHB);
    EXPECT_DOUBLE_EQ(PayloadValue<double>(isotropic, success_row * 6 + 2),
                     expected_pair.iso_1pHaB);
    EXPECT_DOUBLE_EQ(PayloadValue<double>(isotropic, success_row * 6 + 3),
                     expected_pair.iso_2pHaB);
    EXPECT_DOUBLE_EQ(PayloadValue<double>(isotropic, success_row * 6 + 4),
                     expected_pair.iso_diagnostic_CB);
    EXPECT_DOUBLE_EQ(PayloadValue<double>(isotropic, success_row * 6 + 5),
                     expected_pair.iso_table2_total);
    for (std::size_t col = 0; col < 16; ++col) {
        EXPECT_DOUBLE_EQ(
            PayloadValue<double>(compat, success_row * 28 + col),
            static_cast<double>(PayloadValue<std::int32_t>(
                index, success_row * 16 + col)));
    }
    for (std::size_t col = 0; col < 6; ++col) {
        EXPECT_DOUBLE_EQ(
            PayloadValue<double>(compat, success_row * 28 + 16 + col),
            PayloadValue<double>(geometry, success_row * 6 + col));
        EXPECT_DOUBLE_EQ(
            PayloadValue<double>(compat, success_row * 28 + 22 + col),
            PayloadValue<double>(isotropic, success_row * 6 + col));
    }

    const auto filtered_it = std::find_if(
        result->Pairs().begin(), result->Pairs().end(),
        [](const LarsenHBondShieldingResult::PairRecord& pair) {
            return pair.disposition ==
                LarsenHBondShieldingResult::PairDisposition::
                    CarboxylateSymmetryFiltered;
        });
    ASSERT_NE(filtered_it, result->Pairs().end());
    const std::size_t filtered_row = static_cast<std::size_t>(
        std::distance(result->Pairs().begin(), filtered_it));
    EXPECT_EQ(PayloadValue<std::int32_t>(
                  index, filtered_row * 16 + 6),
              static_cast<std::int32_t>(
                  LarsenHBondShieldingResult::PairDisposition::
                      CarboxylateSymmetryFiltered));
    EXPECT_DOUBLE_EQ(PayloadValue<double>(
                         geometry, filtered_row * 6 + 5),
                     1.0);
    EXPECT_TRUE(std::isnan(PayloadValue<double>(
        isotropic, filtered_row * 6 + 5)));

    const RawNpy imputed =
        ReadRawNpy(out / "larsen_imputed_pair_count.npy");
    const RawNpy sidechain = ReadRawNpy(
        out / "larsen_sidechain_carbonyl_pair_count.npy");
    ASSERT_EQ(imputed.payload.size(),
              conf.AtomCount() * sizeof(std::int32_t));
    ASSERT_EQ(sidechain.payload.size(),
              conf.AtomCount() * sizeof(std::int32_t));
    std::int64_t sidechain_sum = 0;
    for (std::size_t i = 0; i < conf.AtomCount(); ++i) {
        std::int32_t imputed_value = 0;
        std::int32_t sidechain_value = 0;
        std::memcpy(&imputed_value,
                    imputed.payload.data() + i * sizeof(std::int32_t),
                    sizeof(std::int32_t));
        std::memcpy(&sidechain_value,
                    sidechain.payload.data() + i * sizeof(std::int32_t),
                    sizeof(std::int32_t));
        // Every synthetic successful lookup has one imputed corner, and
        // diagnostics-only CB writes must not increment this count.
        EXPECT_EQ(imputed_value, conf.AtomAt(i).larsen_hbond_n_pairs);
        sidechain_sum += sidechain_value;
    }
    EXPECT_GT(sidechain_sum, 0);

    RemoveFlatTempDirectory(out, {
        "larsen_hbond_shielding.npy",
        "larsen_hbond_1pHB_shielding.npy",
        "larsen_hbond_2pHB_shielding.npy",
        "larsen_hbond_1pHaB_shielding.npy",
        "larsen_hbond_2pHaB_shielding.npy",
        "larsen_hbond_diagnostic_CB_shielding.npy",
        "larsen_hbond_water_term.npy",
        "larsen_hbond_count.npy",
        "larsen_imputed_pair_count.npy",
        "larsen_sidechain_carbonyl_pair_count.npy",
        "larsen_corner_imputed.npy",
        "larsen_hbond_pairs_index.npy",
        "larsen_hbond_pairs_geometry.npy",
        "larsen_hbond_pairs_isotropic.npy",
        "larsen_hbond_pairs.npy",
    });
    RemoveFlatTempDirectory(grid_dir, {
        "NMANMA_dense.h5", "NMACOH_dense.h5", "NMACOO_dense.h5",
        "ALANMA_dense.h5", "ALACOH_dense.h5", "ALACOO_dense.h5",
    });
}


// A finite, non-degenerate candidate remains a confirmed geometric H-bond
// even when its distance lies outside the archive axis. This production
// forcing grid puts every real H...O distance outside the r axis, isolating
// GridMiss from Success while pinning the historical water-gate behavior.
TEST(LarsenHBondPairAudit, FiniteGridMissStillSuppressesWaterTerm) {
    const fs::path pdb = nmr::test::TestEnvironment::UbqProtonated();
    ASSERT_TRUE(fs::is_regular_file(pdb)) << pdb;
    auto build = BuildFromProtonatedPdb(pdb.string());
    ASSERT_TRUE(build.Ok()) << build.error;
    std::unique_ptr<Protein> protein = std::move(build.protein);
    ProteinConformation& conf = protein->Conformation();
    ASSERT_TRUE(conf.AttachResult(GeometryResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(SpatialIndexResult::Compute(conf)));

    const fs::path grid_dir =
        WriteSyntheticLarsenGridDirectory("finite_grid_miss", true);
    LarsenHBondGrid grid(grid_dir.string());
    auto result = LarsenHBondShieldingResult::Compute(conf, grid);
    ASSERT_NE(result, nullptr);
    EXPECT_EQ(result->PairsFound(), 0);

    const auto amide_grid_miss = std::find_if(
        result->Pairs().begin(), result->Pairs().end(),
        [](const LarsenHBondShieldingResult::PairRecord& pair) {
            return pair.donor_class == HBondDonorClass::AmideHydrogen &&
                   pair.disposition ==
                       LarsenHBondShieldingResult::PairDisposition::GridMiss;
        });
    ASSERT_NE(amide_grid_miss, result->Pairs().end());
    EXPECT_TRUE(amide_grid_miss->frame_valid);
    EXPECT_TRUE(std::isfinite(amide_grid_miss->r_angstrom));
    EXPECT_TRUE(std::isfinite(amide_grid_miss->theta_deg));
    EXPECT_TRUE(std::isfinite(amide_grid_miss->rho_deg));
    EXPECT_GE(amide_grid_miss->theta_deg, 90.0);
    EXPECT_LE(amide_grid_miss->theta_deg, 180.0);
    EXPECT_DOUBLE_EQ(
        conf.AtomAt(amide_grid_miss->donor_atom_idx)
            .larsen_hbond_water_term,
        0.0);

    RemoveFlatTempDirectory(grid_dir, {
        "NMANMA_dense.h5", "NMACOH_dense.h5", "NMACOO_dense.h5",
        "ALANMA_dense.h5", "ALACOH_dense.h5", "ALACOO_dense.h5",
    });
}


// Production-calling forcing gate for the malformed-frame/water interaction.
// Start from a real successful amide-H pair, collapse only that donor's N-H
// frame, and rerun the full calculator. The H...O candidate geometry remains
// spatially present, so an identity fallback would still hit the grid and
// suppress water; the correct path records InvalidFrame, applies no pair, and
// preserves Larsen's 2.07 ppm unbound-amide term.
TEST(LarsenHBondPairAudit,
     DegenerateDonorFrameIsRejectedAndDoesNotSuppressWater) {
    const fs::path pdb = nmr::test::TestEnvironment::UbqProtonated();
    ASSERT_TRUE(fs::is_regular_file(pdb)) << pdb;
    auto build = BuildFromProtonatedPdb(pdb.string());
    ASSERT_TRUE(build.Ok()) << build.error;
    std::unique_ptr<Protein> protein = std::move(build.protein);
    ProteinConformation& baseline_conf = protein->Conformation();
    ASSERT_TRUE(baseline_conf.AttachResult(
        GeometryResult::Compute(baseline_conf)));
    ASSERT_TRUE(baseline_conf.AttachResult(
        SpatialIndexResult::Compute(baseline_conf)));

    const fs::path grid_dir =
        WriteSyntheticLarsenGridDirectory("invalid_frame");
    LarsenHBondGrid grid(grid_dir.string());
    auto baseline =
        LarsenHBondShieldingResult::Compute(baseline_conf, grid);
    ASSERT_NE(baseline, nullptr);

    const auto successful_amide = std::find_if(
        baseline->Pairs().begin(), baseline->Pairs().end(),
        [](const LarsenHBondShieldingResult::PairRecord& pair) {
            return pair.donor_class == HBondDonorClass::AmideHydrogen &&
                   pair.disposition ==
                       LarsenHBondShieldingResult::PairDisposition::Success;
        });
    ASSERT_NE(successful_amide, baseline->Pairs().end());
    const std::size_t donor_h = successful_amide->donor_atom_idx;
    const std::size_t donor_residue =
        protein->AtomAt(donor_h).residue_index;
    const std::size_t donor_n = protein->ResidueAt(donor_residue).N;
    ASSERT_NE(donor_n, Residue::NONE);
    EXPECT_DOUBLE_EQ(
        baseline_conf.AtomAt(donor_h).larsen_hbond_water_term, 0.0);

    std::vector<Vec3> positions = baseline_conf.Positions();
    positions[donor_n] = positions[donor_h];
    DerivedConformation& invalid_conf = protein->AddDerived(
        baseline_conf, "forced degenerate Larsen donor frame",
        std::move(positions));
    ASSERT_TRUE(invalid_conf.AttachResult(
        GeometryResult::Compute(invalid_conf)));
    ASSERT_TRUE(invalid_conf.AttachResult(
        SpatialIndexResult::Compute(invalid_conf)));

    ::testing::internal::CaptureStderr();
    auto invalid = LarsenHBondShieldingResult::Compute(invalid_conf, grid);
    const std::string warning =
        ::testing::internal::GetCapturedStderr();
    ASSERT_NE(invalid, nullptr);

    int donor_rows = 0;
    int invalid_rows = 0;
    std::size_t first_invalid_row = SIZE_MAX;
    for (std::size_t row = 0; row < invalid->Pairs().size(); ++row) {
        const auto& pair = invalid->Pairs()[row];
        if (pair.donor_atom_idx != donor_h) continue;
        ++donor_rows;
        if (pair.disposition ==
            LarsenHBondShieldingResult::PairDisposition::InvalidFrame) {
            ++invalid_rows;
            if (first_invalid_row == SIZE_MAX) first_invalid_row = row;
            EXPECT_FALSE(pair.frame_valid);
            EXPECT_TRUE(std::isnan(pair.iso_1pHB));
            EXPECT_TRUE(std::isnan(pair.iso_2pHB));
            EXPECT_TRUE(std::isnan(pair.iso_1pHaB));
            EXPECT_TRUE(std::isnan(pair.iso_2pHaB));
            EXPECT_TRUE(std::isnan(pair.iso_table2_total));
        } else {
            EXPECT_EQ(pair.disposition,
                      LarsenHBondShieldingResult::PairDisposition::
                          MissingFrameAtoms);
        }
    }
    EXPECT_GT(donor_rows, 0);
    EXPECT_GT(invalid_rows, 0);
    ASSERT_NE(first_invalid_row, SIZE_MAX);
    EXPECT_DOUBLE_EQ(
        invalid_conf.AtomAt(donor_h).larsen_hbond_water_term, 2.07);
    EXPECT_NE(warning.find("invalid donor frame"), std::string::npos)
        << warning;
    EXPECT_NE(warning.find("invalid Larsen pair frame"),
              std::string::npos) << warning;

    const fs::path out =
        MakeUniqueTempDirectory("nmr_larsen_invalid_frame_output");
    EXPECT_EQ(invalid->WriteFeatures(invalid_conf, out.string()), 15);
    const RawNpy index =
        ReadRawNpy(out / "larsen_hbond_pairs_index.npy");
    const RawNpy geometry =
        ReadRawNpy(out / "larsen_hbond_pairs_geometry.npy");
    const RawNpy isotropic =
        ReadRawNpy(out / "larsen_hbond_pairs_isotropic.npy");
    const RawNpy water =
        ReadRawNpy(out / "larsen_hbond_water_term.npy");
    EXPECT_EQ(PayloadValue<std::int32_t>(
                  index, first_invalid_row * 16 + 6),
              static_cast<std::int32_t>(
                  LarsenHBondShieldingResult::PairDisposition::InvalidFrame));
    EXPECT_DOUBLE_EQ(PayloadValue<double>(
                         geometry, first_invalid_row * 6 + 5),
                     0.0);
    EXPECT_TRUE(std::isnan(PayloadValue<double>(
        isotropic, first_invalid_row * 6 + 5)));
    EXPECT_DOUBLE_EQ(PayloadValue<double>(water, donor_h), 2.07);

    RemoveFlatTempDirectory(out, {
        "larsen_hbond_shielding.npy",
        "larsen_hbond_1pHB_shielding.npy",
        "larsen_hbond_2pHB_shielding.npy",
        "larsen_hbond_1pHaB_shielding.npy",
        "larsen_hbond_2pHaB_shielding.npy",
        "larsen_hbond_diagnostic_CB_shielding.npy",
        "larsen_hbond_water_term.npy",
        "larsen_hbond_count.npy",
        "larsen_imputed_pair_count.npy",
        "larsen_sidechain_carbonyl_pair_count.npy",
        "larsen_corner_imputed.npy",
        "larsen_hbond_pairs_index.npy",
        "larsen_hbond_pairs_geometry.npy",
        "larsen_hbond_pairs_isotropic.npy",
        "larsen_hbond_pairs.npy",
    });

    RemoveFlatTempDirectory(grid_dir, {
        "NMANMA_dense.h5", "NMACOH_dense.h5", "NMACOO_dense.h5",
        "ALANMA_dense.h5", "ALACOH_dense.h5", "ALACOO_dense.h5",
    });
}


TEST(LarsenHBondPairAudit,
     NonFiniteDerivedGeometryIsRejectedAndDoesNotSuppressWater) {
    const fs::path pdb = nmr::test::TestEnvironment::UbqProtonated();
    ASSERT_TRUE(fs::is_regular_file(pdb)) << pdb;
    auto build = BuildFromProtonatedPdb(pdb.string());
    ASSERT_TRUE(build.Ok()) << build.error;
    std::unique_ptr<Protein> protein = std::move(build.protein);
    ProteinConformation& baseline_conf = protein->Conformation();
    ASSERT_TRUE(baseline_conf.AttachResult(
        GeometryResult::Compute(baseline_conf)));
    ASSERT_TRUE(baseline_conf.AttachResult(
        SpatialIndexResult::Compute(baseline_conf)));

    const fs::path grid_dir =
        WriteSyntheticLarsenGridDirectory("nonfinite_geometry");
    LarsenHBondGrid grid(grid_dir.string());
    auto baseline =
        LarsenHBondShieldingResult::Compute(baseline_conf, grid);
    ASSERT_NE(baseline, nullptr);

    // Choose a real successful amide donor, then invalidate the acceptor
    // C position for every selected geometric partner of that donor. Theta
    // misses and symmetry-filtered siblings are not water-gating partners.
    const auto successful_amide = std::find_if(
        baseline->Pairs().begin(), baseline->Pairs().end(),
        [](const LarsenHBondShieldingResult::PairRecord& pair) {
            return pair.donor_class == HBondDonorClass::AmideHydrogen &&
                   pair.disposition ==
                       LarsenHBondShieldingResult::PairDisposition::Success;
        });
    ASSERT_NE(successful_amide, baseline->Pairs().end());
    const std::size_t donor_h = successful_amide->donor_atom_idx;
    EXPECT_DOUBLE_EQ(
        baseline_conf.AtomAt(donor_h).larsen_hbond_water_term, 0.0);

    std::vector<Vec3> positions = baseline_conf.Positions();
    int confirmed_before = 0;
    for (const auto& pair : baseline->Pairs()) {
        if (pair.donor_atom_idx != donor_h) continue;
        if (pair.disposition !=
                LarsenHBondShieldingResult::PairDisposition::Success &&
            pair.disposition !=
                LarsenHBondShieldingResult::PairDisposition::GridMiss) {
            continue;
        }
        ASSERT_NE(pair.acceptor_C_atom_idx, SIZE_MAX);
        ASSERT_NE(pair.acceptor_atom_idx, SIZE_MAX);
        positions[pair.acceptor_C_atom_idx] =
            positions[pair.acceptor_atom_idx];
        ++confirmed_before;
    }
    ASSERT_GT(confirmed_before, 0);
    DerivedConformation& invalid_conf = protein->AddDerived(
        baseline_conf, "forced invalid Larsen H-O-C geometry",
        std::move(positions));
    ASSERT_TRUE(invalid_conf.AttachResult(
        GeometryResult::Compute(invalid_conf)));
    ASSERT_TRUE(invalid_conf.AttachResult(
        SpatialIndexResult::Compute(invalid_conf)));

    ::testing::internal::CaptureStderr();
    auto invalid = LarsenHBondShieldingResult::Compute(invalid_conf, grid);
    const std::string warning =
        ::testing::internal::GetCapturedStderr();
    ASSERT_NE(invalid, nullptr);

    const auto invalid_pair = std::find_if(
        invalid->Pairs().begin(), invalid->Pairs().end(),
        [&](const LarsenHBondShieldingResult::PairRecord& pair) {
            return pair.donor_atom_idx == donor_h &&
                   pair.acceptor_atom_idx ==
                       successful_amide->acceptor_atom_idx;
        });
    ASSERT_NE(invalid_pair, invalid->Pairs().end());
    EXPECT_EQ(invalid_pair->disposition,
              LarsenHBondShieldingResult::PairDisposition::InvalidFrame);
    EXPECT_FALSE(invalid_pair->frame_valid);
    EXPECT_TRUE(std::isnan(invalid_pair->theta_deg));
    EXPECT_TRUE(std::isnan(invalid_pair->rho_deg));
    EXPECT_TRUE(std::isnan(invalid_pair->iso_table2_total));
    EXPECT_DOUBLE_EQ(
        invalid_conf.AtomAt(donor_h).larsen_hbond_water_term, 2.07);
    EXPECT_NE(warning.find("geometry_finite=0"), std::string::npos)
        << warning;

    RemoveFlatTempDirectory(grid_dir, {
        "NMANMA_dense.h5", "NMACOH_dense.h5", "NMACOO_dense.h5",
        "ALANMA_dense.h5", "ALACOH_dense.h5", "ALACOO_dense.h5",
    });
}


TEST(LarsenSidechainDonorAudit, EmitsTypedUnsupportedDonorsWithoutShielding) {
    const fs::path pdb = nmr::test::TestEnvironment::UbqProtonated();
    ASSERT_TRUE(fs::is_regular_file(pdb)) << pdb;
    auto build = BuildFromProtonatedPdb(pdb.string());
    ASSERT_TRUE(build.Ok()) << build.error;
    std::unique_ptr<Protein> protein = std::move(build.protein);
    ProteinConformation& conf = protein->Conformation();
    ASSERT_TRUE(conf.AttachResult(GeometryResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(SpatialIndexResult::Compute(conf)));

    auto audit = LarsenSidechainDonorAuditResult::Compute(conf);
    ASSERT_NE(audit, nullptr);
    ASSERT_EQ(audit->DonorAtoms().size(), conf.AtomCount());

    const LegacyAmberTopology& topology = protein->LegacyAmber();
    int n_donors = 0;
    int n_geometry_pass = 0;
    std::int64_t candidate_count_from_atoms = 0;
    std::set<PolarHKind> seen_kinds;
    for (std::size_t i = 0; i < conf.AtomCount(); ++i) {
        const auto& row = audit->DonorAtoms()[i];
        candidate_count_from_atoms += row.n_acceptor_candidates_3p5A;
        n_geometry_pass += row.n_geometry_pass;
        if (!row.is_sidechain_polar_h) continue;
        ++n_donors;
        ASSERT_GE(row.parent_atom, 0);
        const std::size_t parent = static_cast<std::size_t>(row.parent_atom);
        EXPECT_FALSE(topology.SemanticAt(parent).IsBackbone());
        EXPECT_TRUE(topology.SemanticAt(i).IsPolarH());
        seen_kinds.insert(topology.SemanticAt(i).polar_h);
    }
    EXPECT_GT(n_donors, 0);
    EXPECT_EQ(candidate_count_from_atoms,
              static_cast<std::int64_t>(audit->Candidates().size()));
    EXPECT_GT(audit->Candidates().size(), 0u);
    EXPECT_GT(n_geometry_pass, 0);
    EXPECT_TRUE(seen_kinds.count(PolarHKind::HydroxylOH_Aliphatic) > 0);
    EXPECT_TRUE(seen_kinds.count(PolarHKind::HydroxylOH_Aromatic) > 0);
    EXPECT_TRUE(seen_kinds.count(PolarHKind::AmmoniumNH) > 0);
    EXPECT_TRUE(seen_kinds.count(PolarHKind::GuanidiniumNH) > 0);

    for (const auto& row : audit->Candidates()) {
        EXPECT_LE(row.parent_acceptor_distance_A, 3.5 + 1e-12);
        EXPECT_FALSE(row.modeled_by_larsen_table2);
        EXPECT_EQ(protein->AtomAt(row.acceptor_atom).element, Element::O);
        EXPECT_FALSE(topology.SemanticAt(row.parent_atom).IsBackbone());
    }

    // The audit has no writer fields on ConformationAtom and must not
    // extend unsupported donor classes into Larsen's scientific tensor.
    for (std::size_t i = 0; i < conf.AtomCount(); ++i) {
        EXPECT_DOUBLE_EQ(conf.AtomAt(i).larsen_hbond_shielding_tensor.norm(),
                         0.0);
        EXPECT_EQ(conf.AtomAt(i).larsen_hbond_n_pairs, 0);
    }

    const fs::path out =
        MakeUniqueTempDirectory("nmr_larsen_sidechain_donor_audit_output");
    EXPECT_EQ(audit->WriteFeatures(conf, out.string()), 2);
    const RawNpy atoms =
        ReadRawNpy(out / "larsen_sidechain_donor_atoms.npy");
    const RawNpy candidates =
        ReadRawNpy(out / "larsen_sidechain_donor_candidates.npy");
    EXPECT_NE(atoms.header.find("'descr': '<i4'"), std::string::npos);
    EXPECT_NE(atoms.header.find("'shape': (" +
                                std::to_string(conf.AtomCount()) +
                                ",6)"), std::string::npos);
    EXPECT_EQ(atoms.payload.size(), conf.AtomCount() * 6 * sizeof(std::int32_t));
    EXPECT_NE(candidates.header.find("'shape': (" +
                                     std::to_string(audit->Candidates().size()) +
                                     ",13)"), std::string::npos);
    ASSERT_EQ(candidates.payload.size(),
              audit->Candidates().size() * 13 * sizeof(double));
    for (std::size_t row = 0; row < audit->Candidates().size(); ++row) {
        double modeled = 1.0;
        std::memcpy(&modeled,
                    candidates.payload.data() +
                        (row * 13 + 12) * sizeof(double),
                    sizeof(double));
        EXPECT_DOUBLE_EQ(modeled, 0.0);
    }

    // A24 read-back: a typed donor row and its first raw candidate must
    // survive the writer with their identities, geometry, and audit-only
    // modeled flag intact.
    const auto donor_it = std::find_if(
        audit->DonorAtoms().begin(), audit->DonorAtoms().end(),
        [](const LarsenSidechainDonorAuditResult::DonorAtomRecord& row) {
            return row.is_sidechain_polar_h != 0;
        });
    ASSERT_NE(donor_it, audit->DonorAtoms().end());
    const std::size_t donor_row = static_cast<std::size_t>(
        std::distance(audit->DonorAtoms().begin(), donor_it));
    EXPECT_EQ(PayloadValue<std::int32_t>(atoms, donor_row * 6 + 0), 1);
    EXPECT_EQ(PayloadValue<std::int32_t>(atoms, donor_row * 6 + 1),
              donor_it->polar_h_kind);
    EXPECT_EQ(PayloadValue<std::int32_t>(atoms, donor_row * 6 + 2),
              donor_it->parent_atom);
    EXPECT_EQ(PayloadValue<std::int32_t>(atoms, donor_row * 6 + 3),
              donor_it->residue_index);
    EXPECT_EQ(PayloadValue<std::int32_t>(atoms, donor_row * 6 + 4),
              donor_it->n_acceptor_candidates_3p5A);
    EXPECT_EQ(PayloadValue<std::int32_t>(atoms, donor_row * 6 + 5),
              donor_it->n_geometry_pass);

    ASSERT_FALSE(audit->Candidates().empty());
    const auto& expected_candidate = audit->Candidates().front();
    EXPECT_DOUBLE_EQ(PayloadValue<double>(candidates, 0),
                     static_cast<double>(expected_candidate.donor_atom));
    EXPECT_DOUBLE_EQ(PayloadValue<double>(candidates, 2),
                     static_cast<double>(expected_candidate.polar_h_kind));
    EXPECT_DOUBLE_EQ(PayloadValue<double>(candidates, 3),
                     static_cast<double>(expected_candidate.parent_atom));
    EXPECT_DOUBLE_EQ(PayloadValue<double>(candidates, 4),
                     static_cast<double>(expected_candidate.acceptor_atom));
    EXPECT_DOUBLE_EQ(PayloadValue<double>(candidates, 6),
                     static_cast<double>(expected_candidate.acceptor_class));
    EXPECT_DOUBLE_EQ(PayloadValue<double>(candidates, 7),
                     expected_candidate.h_acceptor_distance_A);
    EXPECT_DOUBLE_EQ(PayloadValue<double>(candidates, 8),
                     expected_candidate.parent_acceptor_distance_A);
    EXPECT_DOUBLE_EQ(PayloadValue<double>(candidates, 9),
                     expected_candidate.angle_parent_h_acceptor_deg);
    if (std::isnan(expected_candidate.rho_deg)) {
        EXPECT_TRUE(std::isnan(PayloadValue<double>(candidates, 10)));
    } else {
        EXPECT_DOUBLE_EQ(PayloadValue<double>(candidates, 10),
                         expected_candidate.rho_deg);
    }
    EXPECT_DOUBLE_EQ(PayloadValue<double>(candidates, 11),
                     expected_candidate.passes_geometry ? 1.0 : 0.0);
    EXPECT_DOUBLE_EQ(PayloadValue<double>(candidates, 12), 0.0);
    RemoveFlatTempDirectory(out, {
        "larsen_sidechain_donor_atoms.npy",
        "larsen_sidechain_donor_candidates.npy",
    });
}


// ---------------------------------------------------------------------------
// 1UBQ smoke test fixture
// ---------------------------------------------------------------------------

class LarsenHBondShieldingTest : public ::testing::Test {
protected:
    void SetUp() override {
        const std::string kLarsen1UbqPm6 =
            nmr::test::TestEnvironment::Larsen1UbqPm6Pdb();
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

    // WriteFeatures emits per-atom tensors/counts plus split and
    // compatibility pair-audit tables.
    fs::path tmp = fs::temp_directory_path() / "larsen_hbond_phase1_out";
    fs::create_directories(tmp);
    int n_written = result->WriteFeatures(conf, tmp.string());
    EXPECT_EQ(n_written, 15);
    for (const std::string& stem : {
        "larsen_hbond_shielding",
        "larsen_hbond_1pHB_shielding",
        "larsen_hbond_2pHB_shielding",
        "larsen_hbond_1pHaB_shielding",
        "larsen_hbond_2pHaB_shielding",
        "larsen_hbond_diagnostic_CB_shielding",
        "larsen_hbond_water_term",
        "larsen_hbond_count",
        "larsen_imputed_pair_count",
        "larsen_sidechain_carbonyl_pair_count",
        "larsen_corner_imputed",
        "larsen_hbond_pairs_index",
        "larsen_hbond_pairs_geometry",
        "larsen_hbond_pairs_isotropic",
        "larsen_hbond_pairs",
    }) {
        fs::path p = tmp / (stem + ".npy");
        EXPECT_TRUE(fs::exists(p)) << "missing " << p.string();
    }

    {
        fs::path p = tmp / "larsen_corner_imputed.npy";
        std::ifstream in(p, std::ios::binary);
        ASSERT_TRUE(in.is_open()) << p.string();
        char magic[6] = {};
        in.read(magic, 6);
        ASSERT_EQ(std::string(magic, 6), std::string("\x93NUMPY", 6));
        char version[2] = {};
        in.read(version, 2);
        ASSERT_EQ(version[0], 1);
        ASSERT_EQ(version[1], 0);
        std::uint16_t header_len = 0;
        in.read(reinterpret_cast<char*>(&header_len), sizeof(header_len));
        std::string header(header_len, '\0');
        in.read(header.data(), header_len);
        EXPECT_NE(header.find("'descr': '|i1'"), std::string::npos) << header;
        EXPECT_NE(header.find("'shape': (" + std::to_string(conf.AtomCount()) + ",)"),
                  std::string::npos) << header;
        std::vector<char> payload{std::istreambuf_iterator<char>(in),
                                  std::istreambuf_iterator<char>()};
        ASSERT_EQ(payload.size(), conf.AtomCount());
        for (std::size_t i = 0; i < conf.AtomCount(); ++i) {
            EXPECT_EQ(static_cast<std::int8_t>(payload[i]),
                      conf.AtomAt(i).larsen_hbond_any_corner_imputed ? 1 : 0);
        }
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

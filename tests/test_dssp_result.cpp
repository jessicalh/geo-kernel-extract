#include "TestEnvironment.h"
#include "DirectionalTestHelpers.h"
#include <gtest/gtest.h>
#include "Atom.h"
#include "PdbFileReader.h"
#include "GeometryResult.h"
#include "DsspResult.h"
#include "Dssp8TimeSeriesTrajectoryResult.h"
#include "PhysicalConstants.h"
#include "NpyWriter.h"
#include "Protein.h"
#include "Residue.h"
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <utility>
#include <unistd.h>
#include <vector>
#include <cstdint>
#include <cmath>

using namespace nmr;
namespace fs = std::filesystem;

namespace {

struct NpyArray {
    std::string descr;
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
    const std::string descr_key = "'descr': '";
    auto descr_pos = header.find(descr_key);
    EXPECT_NE(descr_pos, std::string::npos);
    if (descr_pos == std::string::npos) return arr;
    descr_pos += descr_key.size();
    auto descr_end = header.find('\'', descr_pos);
    EXPECT_NE(descr_end, std::string::npos);
    if (descr_end == std::string::npos) return arr;
    arr.descr = header.substr(descr_pos, descr_end - descr_pos);
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

void ExpectShapeAndDescr(const NpyArray& arr,
                         const std::vector<size_t>& shape,
                         const std::string& descr) {
    EXPECT_EQ(arr.shape, shape);
    EXPECT_EQ(arr.descr, descr);
}

void RemoveDsspOutput(const fs::path& output_dir) {
    std::error_code ec;
    for (const char* filename : {
             "dssp_observed.npy", "dssp_backbone.npy", "dssp_ss8.npy",
             "dssp_ppii.npy", "dssp_hbond_energy.npy", "dssp_chi.npy",
             "dssp_torsion_angle.npy", "dssp_torsion_sin.npy",
             "dssp_torsion_cos.npy", "dssp_torsion_valid.npy",
             "dssp_hbond_partner_residue_index.npy"}) {
        fs::remove(output_dir / filename, ec);
    }
    fs::remove(output_dir, ec);
}

// Two residues carrying the B15 analytic chi fixtures. Residue 0 is +60°;
// residue 1 is the y-reflected -60° case, translated away so topology
// perception cannot connect the fixtures. Empty atom names intentionally
// use Protein's supported substrate-free calculator-fixture path.
std::unique_ptr<Protein> BuildFixedChiSignFixture() {
    auto protein = std::make_unique<Protein>();
    std::vector<Vec3> positions;

    auto add_residue = [&](double x_offset, double d_y) {
        Residue residue;
        residue.type = AminoAcid::Unknown;
        residue.chain_id = "A";
        residue.sequence_number =
            static_cast<int>(protein->ResidueCount()) + 1;
        const size_t ri = protein->AddResidue(std::move(residue));

        const std::vector<Vec3> coords = {
            Vec3(x_offset + 1.0, 0.0, 0.0),
            Vec3(x_offset + 0.0, 0.0, 0.0),
            Vec3(x_offset + 0.0, 0.0, 1.0),
            Vec3(x_offset + 0.5, d_y, 1.0),
        };
        size_t chi_atoms[4] = {};
        for (size_t k = 0; k < coords.size(); ++k) {
            auto atom = Atom::Create(Element::C);
            atom->residue_index = ri;
            chi_atoms[k] = protein->AddAtom(std::move(atom));
            protein->MutableResidueAt(ri).atom_indices.push_back(chi_atoms[k]);
            positions.push_back(coords[k]);
        }
        for (size_t k = 0; k < 4; ++k) {
            protein->MutableResidueAt(ri).chi[0].a[k] = chi_atoms[k];
        }
    };

    const double root3_over_2 = std::sqrt(3.0) / 2.0;
    add_residue(/*x_offset=*/0.0, /*d_y=*/-root3_over_2);  // +60°
    add_residue(/*x_offset=*/10.0, /*d_y=*/root3_over_2);  // -60°

    protein->FinalizeConstruction(positions);
    protein->AddConformation(std::move(positions), "fixed chi signs");
    return protein;
}

}  // namespace

TEST(DsspResultHelpers, PpiiPredicateAndSs8Compatibility) {
    EXPECT_TRUE(DsspIsPpii('P'));
    EXPECT_FALSE(DsspIsPpii('C'));
    EXPECT_EQ(Dssp8TimeSeriesTrajectoryResult::Ss8Code('P'), 7u);
}

TEST(DsspResultHelpers,
     ChiNpyMatchesCanonicalFixedCoordinateSignsAndShape) {
    auto protein = BuildFixedChiSignFixture();
    auto dssp = DsspResult::CreateForTesting(
        std::vector<DsspResidue>(protein->ResidueCount()));
    const auto& conf = protein->Conformation();

    const fs::path output_dir = fs::temp_directory_path() /
        ("dssp_chi_sign_" + std::to_string(::getpid()));
    fs::create_directories(output_dir);
    ASSERT_EQ(dssp->WriteFeatures(conf, output_dir.string()), 11);

    const auto chi = ReadNpy(output_dir / "dssp_chi.npy");
    ExpectShapeAndDescr(chi, {conf.AtomCount(), 12}, "<f8");
    ASSERT_EQ(chi.bytes.size(), conf.AtomCount() * 12 * sizeof(double));
    const double* values = DataAs<double>(chi);

    const double expected_cos = 0.5;
    const double expected_sin = std::sqrt(3.0) / 2.0;
    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        const double sign = (ri == 0) ? 1.0 : -1.0;
        for (size_t ai : protein->ResidueAt(ri).atom_indices) {
            EXPECT_NEAR(values[ai * 12 + 0], expected_cos, 1e-12);
            EXPECT_NEAR(values[ai * 12 + 1], sign * expected_sin, 1e-12);
            EXPECT_DOUBLE_EQ(values[ai * 12 + 2], 1.0);
            for (int k = 1; k < 4; ++k) {
                EXPECT_DOUBLE_EQ(values[ai * 12 + k * 3 + 0], 0.0);
                EXPECT_DOUBLE_EQ(values[ai * 12 + k * 3 + 1], 0.0);
                EXPECT_DOUBLE_EQ(values[ai * 12 + k * 3 + 2], 0.0);
            }
        }
    }

    // The analytic coordinates are pinned to ±60°. This independently
    // checks DSSP's unchanged cos/sin/exists writer and catches either a
    // sign or column regression.
    std::error_code ec;
    for (const char* f : {"dssp_observed.npy", "dssp_backbone.npy",
                          "dssp_ss8.npy", "dssp_ppii.npy",
                          "dssp_hbond_energy.npy", "dssp_chi.npy"}) {
        fs::remove(output_dir / f, ec);
    }
    fs::remove(output_dir, ec);
}


class DsspResultTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!std::filesystem::exists(nmr::test::TestEnvironment::UbqProtonated())) {
            GTEST_SKIP() << "1UBQ PDB not found";
        }
        auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
        if (!r.Ok()) GTEST_SKIP() << r.error;
        protein = std::move(r.protein);
    }
    std::unique_ptr<Protein> protein;
};

TEST_F(DsspResultTest, ComputeSucceeds) {
    auto& conf = protein->Conformation();
    auto dssp = DsspResult::Compute(conf);
    ASSERT_NE(dssp, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(dssp)));
}

TEST_F(DsspResultTest,
       DirectionalSignedGeometryAndPpiiRerunProductionLibdsspOnO3) {
    using nmr::test::directional::Positions;
    using nmr::test::directional::SeededTransform;

    auto& original = protein->Conformation();
    auto original_result = DsspResult::Compute(original);
    ASSERT_NE(original_result, nullptr);
    ASSERT_TRUE(original.AttachResult(std::move(original_result)));
    const fs::path root = fs::temp_directory_path() /
        ("dssp_directional_" + std::to_string(::getpid()));
    RemoveDsspOutput(root / "original");
    RemoveDsspOutput(root / "proper");
    RemoveDsspOutput(root / "improper");
    std::error_code ec;
    fs::remove(root, ec);
    const fs::path original_dir = root / "original";
    ASSERT_TRUE(fs::create_directories(original_dir));
    ASSERT_EQ(original.Result<DsspResult>().WriteFeatures(
                  original, original_dir.string()), 11);
    const NpyArray original_observed =
        ReadNpy(original_dir / "dssp_observed.npy");
    const NpyArray original_ppii = ReadNpy(original_dir / "dssp_ppii.npy");
    const NpyArray original_torsion =
        ReadNpy(original_dir / "dssp_torsion_angle.npy");
    const NpyArray original_valid =
        ReadNpy(original_dir / "dssp_torsion_valid.npy");
    const NpyArray original_backbone =
        ReadNpy(original_dir / "dssp_backbone.npy");
    const NpyArray original_ss8 = ReadNpy(original_dir / "dssp_ss8.npy");
    const NpyArray original_chi = ReadNpy(original_dir / "dssp_chi.npy");
    const NpyArray original_sin =
        ReadNpy(original_dir / "dssp_torsion_sin.npy");
    const NpyArray original_cos =
        ReadNpy(original_dir / "dssp_torsion_cos.npy");
    const NpyArray original_hbond_energy =
        ReadNpy(original_dir / "dssp_hbond_energy.npy");
    const NpyArray original_hbond_partner = ReadNpy(
        original_dir / "dssp_hbond_partner_residue_index.npy");
    ExpectShapeAndDescr(original_observed, {original.AtomCount()}, "|i1");
    ExpectShapeAndDescr(original_ppii, {original.AtomCount()}, "|i1");
    ExpectShapeAndDescr(
        original_torsion, {protein->ResidueCount(), 6}, "<f8");
    ExpectShapeAndDescr(
        original_valid, {protein->ResidueCount(), 6}, "|u1");
    ExpectShapeAndDescr(
        original_backbone, {original.AtomCount(), 5}, "<f8");
    ExpectShapeAndDescr(original_ss8, {original.AtomCount(), 8}, "<f8");
    ExpectShapeAndDescr(original_chi, {original.AtomCount(), 12}, "<f8");
    ExpectShapeAndDescr(
        original_sin, {protein->ResidueCount(), 6}, "<f8");
    ExpectShapeAndDescr(
        original_cos, {protein->ResidueCount(), 6}, "<f8");
    ExpectShapeAndDescr(
        original_hbond_energy, {original.AtomCount(), 4}, "<f8");
    ExpectShapeAndDescr(
        original_hbond_partner, {original.AtomCount(), 4}, "<i4");

    const int8_t* observed_0 = DataAs<int8_t>(original_observed);
    const int8_t* ppii_0 = DataAs<int8_t>(original_ppii);
    const double* torsion_0 = DataAs<double>(original_torsion);
    const std::uint8_t* valid_0 = DataAs<std::uint8_t>(original_valid);
    const double* backbone_0 = DataAs<double>(original_backbone);
    const double* ss8_0 = DataAs<double>(original_ss8);
    const double* chi_0 = DataAs<double>(original_chi);
    const double* sin_0 = DataAs<double>(original_sin);
    const double* cos_0 = DataAs<double>(original_cos);
    const double* hbond_energy_0 = DataAs<double>(original_hbond_energy);
    const std::int32_t* hbond_partner_0 =
        DataAs<std::int32_t>(original_hbond_partner);
    std::size_t original_ppii_count = 0;
    for (std::size_t atom = 0; atom < original.AtomCount(); ++atom)
        original_ppii_count += ppii_0[atom] == 1 ? 1u : 0u;
    ASSERT_GT(original_ppii_count, 0u)
        << "1UBQ no longer exercises libdssp's PPII class";

    // DsspResult reruns libdssp through its PDB interchange path.  The rigidly
    // transformed coordinates are rounded to 0.001 A there, and libdssp's
    // solvent-accessible-surface calculation amplifies that quantization.
    // For seed 0xD55A2026 the recorded proper-rotation maximum is 3.755501 A^2;
    // 4.0 A^2 is the explicit serialization/discretization bound.  This is
    // separate from the signed torsion and exact categorical assertions.
    constexpr double kLibdsspSasaAbsToleranceA2 = 4.0;
    double max_o3_sasa_error_A2 = 0.0;

    for (const bool improper : {false, true}) {
        const auto transform =
            SeededTransform(0xD55A2026ULL, improper);
        auto& moved = protein->AddConformation(
            Positions(transform, original.Positions()),
            improper ? "DSSP improper rerun" : "DSSP proper rerun");
        auto moved_result = DsspResult::Compute(moved);
        ASSERT_NE(moved_result, nullptr);
        ASSERT_TRUE(moved.AttachResult(std::move(moved_result)));
        const fs::path moved_dir = root /
            (improper ? "improper" : "proper");
        ASSERT_TRUE(fs::create_directories(moved_dir));
        ASSERT_EQ(moved.Result<DsspResult>().WriteFeatures(
                      moved, moved_dir.string()), 11);
        const NpyArray moved_observed =
            ReadNpy(moved_dir / "dssp_observed.npy");
        const NpyArray moved_ppii = ReadNpy(moved_dir / "dssp_ppii.npy");
        const NpyArray moved_torsion =
            ReadNpy(moved_dir / "dssp_torsion_angle.npy");
        const NpyArray moved_valid =
            ReadNpy(moved_dir / "dssp_torsion_valid.npy");
        const NpyArray moved_backbone =
            ReadNpy(moved_dir / "dssp_backbone.npy");
        const NpyArray moved_ss8 = ReadNpy(moved_dir / "dssp_ss8.npy");
        const NpyArray moved_chi = ReadNpy(moved_dir / "dssp_chi.npy");
        const NpyArray moved_sin =
            ReadNpy(moved_dir / "dssp_torsion_sin.npy");
        const NpyArray moved_cos =
            ReadNpy(moved_dir / "dssp_torsion_cos.npy");
        const NpyArray moved_hbond_energy =
            ReadNpy(moved_dir / "dssp_hbond_energy.npy");
        const NpyArray moved_hbond_partner = ReadNpy(
            moved_dir / "dssp_hbond_partner_residue_index.npy");
        ASSERT_EQ(moved_observed.shape, original_observed.shape);
        ASSERT_EQ(moved_ppii.shape, original_ppii.shape);
        ASSERT_EQ(moved_torsion.shape, original_torsion.shape);
        ASSERT_EQ(moved_valid.shape, original_valid.shape);
        ASSERT_EQ(moved_backbone.shape, original_backbone.shape);
        ASSERT_EQ(moved_ss8.shape, original_ss8.shape);
        ASSERT_EQ(moved_chi.shape, original_chi.shape);
        ASSERT_EQ(moved_sin.shape, original_sin.shape);
        ASSERT_EQ(moved_cos.shape, original_cos.shape);
        ASSERT_EQ(moved_hbond_energy.shape, original_hbond_energy.shape);
        ASSERT_EQ(moved_hbond_partner.shape, original_hbond_partner.shape);
        const int8_t* observed_1 = DataAs<int8_t>(moved_observed);
        const int8_t* ppii_1 = DataAs<int8_t>(moved_ppii);
        const double* torsion_1 = DataAs<double>(moved_torsion);
        const std::uint8_t* valid_1 = DataAs<std::uint8_t>(moved_valid);
        const double* backbone_1 = DataAs<double>(moved_backbone);
        const double* ss8_1 = DataAs<double>(moved_ss8);
        const double* chi_1 = DataAs<double>(moved_chi);
        const double* sin_1 = DataAs<double>(moved_sin);
        const double* cos_1 = DataAs<double>(moved_cos);
        const double* hbond_energy_1 = DataAs<double>(moved_hbond_energy);
        const std::int32_t* hbond_partner_1 =
            DataAs<std::int32_t>(moved_hbond_partner);

        std::size_t ppii_changes = 0;
        for (std::size_t atom = 0; atom < original.AtomCount(); ++atom) {
            EXPECT_EQ(observed_1[atom], observed_0[atom])
                << "dssp_observed.npy atom=" << atom;
            if (improper) {
                ppii_changes += ppii_1[atom] != ppii_0[atom] ? 1u : 0u;
            } else {
                EXPECT_EQ(ppii_1[atom], ppii_0[atom])
                    << "dssp_ppii.npy atom=" << atom;
            }
            double source_ss8_sum = 0.0;
            double moved_ss8_sum = 0.0;
            for (std::size_t column = 0; column < 8; ++column) {
                const double source_value = ss8_0[atom * 8 + column];
                const double moved_value = ss8_1[atom * 8 + column];
                EXPECT_TRUE(source_value == 0.0 || source_value == 1.0);
                EXPECT_TRUE(moved_value == 0.0 || moved_value == 1.0);
                source_ss8_sum += source_value;
                moved_ss8_sum += moved_value;
                // The ordinary eight-state surface collapses libdssp's
                // reflection-sensitive P extension into coil.  Its emitted
                // one-hot category is therefore exactly invariant for both
                // proper and improper rigid transforms.
                EXPECT_DOUBLE_EQ(moved_value, source_value)
                    << "dssp_ss8.npy atom=" << atom
                    << " column=" << column;
            }
            EXPECT_TRUE(source_ss8_sum == 0.0 || source_ss8_sum == 1.0);
            EXPECT_TRUE(moved_ss8_sum == 0.0 || moved_ss8_sum == 1.0);
            // dssp_backbone.npy columns 2:5 are the serialized scalar SASA
            // and categorical helix/sheet flags.  SASA is invariant up to
            // the recorded PDB/libdssp quantization envelope; the flags are
            // exact after the PPII-to-coil collapse under all O(3).
            if (std::isnan(backbone_0[atom * 5 + 2])) {
                EXPECT_TRUE(std::isnan(backbone_1[atom * 5 + 2]));
            } else {
                max_o3_sasa_error_A2 = std::max(
                    max_o3_sasa_error_A2,
                    std::abs(backbone_1[atom * 5 + 2] -
                             backbone_0[atom * 5 + 2]));
                EXPECT_NEAR(backbone_1[atom * 5 + 2],
                            backbone_0[atom * 5 + 2],
                            kLibdsspSasaAbsToleranceA2);
            }
            EXPECT_DOUBLE_EQ(backbone_1[atom * 5 + 3],
                             backbone_0[atom * 5 + 3]);
            EXPECT_DOUBLE_EQ(backbone_1[atom * 5 + 4],
                             backbone_0[atom * 5 + 4]);

            // Static partner identities and their energy slots are checked
            // together.  Their physical law is invariant, while the
            // production libdssp boundary has the same explicitly recorded
            // 0.001-A quantization envelope as the trajectory owner.
            for (std::size_t slot = 0; slot < 4; ++slot) {
                const std::size_t offset = atom * 4 + slot;
                EXPECT_EQ(hbond_partner_1[offset], hbond_partner_0[offset])
                    << "dssp_hbond_partner_residue_index.npy atom=" << atom
                    << " slot=" << slot;
                if (std::isnan(hbond_energy_0[offset])) {
                    EXPECT_TRUE(std::isnan(hbond_energy_1[offset]))
                        << "dssp_hbond_energy.npy atom=" << atom
                        << " slot=" << slot;
                } else {
                    EXPECT_NEAR(hbond_energy_1[offset],
                                hbond_energy_0[offset], 5.0e-3)
                        << "dssp_hbond_energy.npy atom=" << atom
                        << " slot=" << slot;
                    if (hbond_partner_0[offset] < 0)
                        EXPECT_DOUBLE_EQ(hbond_energy_0[offset], 0.0);
                }
            }
        }
        if (improper) {
            EXPECT_GT(ppii_changes, 0u)
                << "dssp_ppii.npy was incorrectly treated as reflection-even";
        }

        std::size_t signed_angles_checked = 0;
        for (std::size_t i = 0;
             i < protein->ResidueCount() * 6; ++i) {
            EXPECT_EQ(valid_1[i], valid_0[i]);
            if (!valid_0[i]) continue;
            const double expected = improper ? -torsion_0[i] : torsion_0[i];
            EXPECT_NEAR(
                nmr::test::directional::CircularDifference(
                    torsion_1[i], expected),
                0.0, 5.0e-2)
                << "dssp_torsion_angle.npy flat_component=" << i;
            ++signed_angles_checked;

            const double expected_sin = improper ? -sin_0[i] : sin_0[i];
            EXPECT_NEAR(sin_1[i], expected_sin, 5.0e-2)
                << "dssp_torsion_sin.npy flat_component=" << i;
            EXPECT_NEAR(cos_1[i], cos_0[i], 5.0e-2)
                << "dssp_torsion_cos.npy flat_component=" << i;
        }
        EXPECT_GT(signed_angles_checked, 100u);

        // The atom-broadcast DSSP backbone table owns the same signed phi/psi
        // information in columns 0:2.  Check its serialized payload directly;
        // columns 2:5 were checked independently above because their scalar
        // and collapsed-category laws differ from these signed angles.
        for (std::size_t atom = 0; atom < original.AtomCount(); ++atom) {
            for (std::size_t column = 0; column < 2; ++column) {
                const double source = backbone_0[atom * 5 + column];
                const double emitted = backbone_1[atom * 5 + column];
                if (std::isnan(source)) {
                    EXPECT_TRUE(std::isnan(emitted));
                } else {
                    EXPECT_NEAR(
                        nmr::test::directional::CircularDifference(
                            emitted, improper ? -source : source),
                        0.0, 5.0e-2)
                        << "dssp_backbone.npy atom=" << atom
                        << " column=" << column;
                }
            }
        }

        // dssp_chi.npy is computed by the owning writer from the transformed
        // conformation coordinates: cos/exists are invariant and sin is a
        // pseudoscalar.  Missing chis retain the exact zero/zero/zero row.
        for (std::size_t atom = 0; atom < original.AtomCount(); ++atom) {
            for (std::size_t chi = 0; chi < 4; ++chi) {
                const std::size_t base = atom * 12 + chi * 3;
                EXPECT_NEAR(chi_1[base], chi_0[base], 2.0e-11)
                    << "dssp_chi.npy cos atom=" << atom << " chi=" << chi;
                EXPECT_NEAR(chi_1[base + 1],
                            improper ? -chi_0[base + 1] : chi_0[base + 1],
                            2.0e-11)
                    << "dssp_chi.npy sin atom=" << atom << " chi=" << chi;
                EXPECT_DOUBLE_EQ(chi_1[base + 2], chi_0[base + 2])
                    << "dssp_chi.npy exists atom=" << atom << " chi=" << chi;
            }
        }
    }

    std::cerr << "DSSP O(3)-rerun serialized SASA max_abs_A2="
              << max_o3_sasa_error_A2
              << " tolerance_A2=" << kLibdsspSasaAbsToleranceA2 << '\n';

    RemoveDsspOutput(root / "original");
    RemoveDsspOutput(root / "proper");
    RemoveDsspOutput(root / "improper");
    fs::remove(root, ec);
}

TEST_F(DsspResultTest, HasSecondaryStructure) {
    auto& conf = protein->Conformation();
    auto dssp = DsspResult::Compute(conf);
    ASSERT_NE(dssp, nullptr);
    conf.AttachResult(std::move(dssp));

    const auto& d = conf.Result<DsspResult>();

    // Check that we have per-residue data
    EXPECT_EQ(d.AllResidues().size(), protein->ResidueCount());
}

TEST_F(DsspResultTest, HelixResiduesDetected) {
    auto& conf = protein->Conformation();
    auto dssp = DsspResult::Compute(conf);
    ASSERT_NE(dssp, nullptr);
    conf.AttachResult(std::move(dssp));

    const auto& d = conf.Result<DsspResult>();

    // 1UBQ has an alpha helix from residue 23 to 34 (PDB numbering).
    // Find the residue indices corresponding to these sequence numbers.
    int helix_count = 0;
    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        int seq = protein->ResidueAt(ri).sequence_number;
        if (seq >= 23 && seq <= 34) {
            char ss = d.SecondaryStructure(ri);
            // H = alpha helix, G = 3-10 helix, I = pi helix
            if (ss == 'H' || ss == 'G' || ss == 'I')
                helix_count++;
        }
    }
    // At least some residues in this range should be helical
    EXPECT_GT(helix_count, 5) << "Expected helix residues in range 23-34";
}

TEST_F(DsspResultTest, SheetResiduesDetected) {
    auto& conf = protein->Conformation();
    auto dssp = DsspResult::Compute(conf);
    ASSERT_NE(dssp, nullptr);
    conf.AttachResult(std::move(dssp));

    const auto& d = conf.Result<DsspResult>();

    // 1UBQ has beta strands. Check for at least some E/B assignments.
    int sheet_count = 0;
    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        char ss = d.SecondaryStructure(ri);
        if (ss == 'E' || ss == 'B')
            sheet_count++;
    }
    EXPECT_GT(sheet_count, 5) << "Expected some beta strand residues";
}

TEST_F(DsspResultTest, PhiPsiInRange) {
    auto& conf = protein->Conformation();
    auto dssp = DsspResult::Compute(conf);
    ASSERT_NE(dssp, nullptr);
    conf.AttachResult(std::move(dssp));

    const auto& d = conf.Result<DsspResult>();

    // Phi and Psi should be in [-pi, pi] for internal residues.
    // Terminal residues may have special values (360 from DSSP = 2*pi).
    int valid_count = 0;
    for (size_t ri = 1; ri < protein->ResidueCount() - 1; ++ri) {
        double phi = d.Phi(ri);
        double psi = d.Psi(ri);
        // Allow some tolerance for DSSP conversion: [-pi-0.1, pi+0.1]
        // DSSP returns 360 for undefined angles, which converts to ~6.28 rad
        if (std::abs(phi) <= PI + 0.1 && std::abs(psi) <= PI + 0.1)
            valid_count++;
    }
    // Most internal residues should have valid phi/psi
    EXPECT_GT(valid_count, static_cast<int>(protein->ResidueCount()) / 2);
}

TEST_F(DsspResultTest, SASANonNegative) {
    auto& conf = protein->Conformation();
    auto dssp = DsspResult::Compute(conf);
    ASSERT_NE(dssp, nullptr);
    conf.AttachResult(std::move(dssp));

    const auto& d = conf.Result<DsspResult>();

    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        EXPECT_GE(d.SASA(ri), 0.0) << "Negative SASA at residue " << ri;
    }
}

TEST_F(DsspResultTest, SomeSASANonZero) {
    auto& conf = protein->Conformation();
    auto dssp = DsspResult::Compute(conf);
    ASSERT_NE(dssp, nullptr);
    conf.AttachResult(std::move(dssp));

    const auto& d = conf.Result<DsspResult>();

    int nonzero = 0;
    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        if (d.SASA(ri) > 0.0) nonzero++;
    }
    // Surface residues should have nonzero SASA
    EXPECT_GT(nonzero, 10);
}

TEST_F(DsspResultTest, dssp_unobserved_residue_npy) {
    auto& conf = protein->Conformation();
    auto dssp = DsspResult::Compute(conf);
    ASSERT_NE(dssp, nullptr);

    auto& residues = const_cast<std::vector<DsspResidue>&>(dssp->AllResidues());
    ASSERT_GE(residues.size(), 2u);
    ASSERT_FALSE(protein->ResidueAt(0).atom_indices.empty());
    ASSERT_FALSE(protein->ResidueAt(1).atom_indices.empty());

    residues[0].observed = true;
    residues[0].secondary_structure = 'H';
    residues[0].phi = 0.4;
    residues[0].psi = -0.7;
    residues[0].sasa = 12.5;
    residues[0].acceptors[0].energy = -1.1;
    residues[0].acceptors[1].energy = -2.2;
    residues[0].donors[0].energy = -3.3;
    residues[0].donors[1].energy = -4.4;

    residues[1].observed = false;
    residues[1].secondary_structure = 'C';
    residues[1].phi = 1.2;
    residues[1].psi = 2.3;
    residues[1].sasa = 99.0;
    residues[1].acceptors[0].energy = -9.0;
    residues[1].acceptors[1].energy = -8.0;
    residues[1].donors[0].energy = -7.0;
    residues[1].donors[1].energy = -6.0;

    const fs::path output_dir = fs::temp_directory_path() /
        ("dssp_unobserved_" + std::to_string(::getpid()));
    fs::create_directories(output_dir);
    ASSERT_EQ(dssp->WriteFeatures(conf, output_dir.string()), 11);

    const size_t N = conf.AtomCount();
    const size_t obs_ai = protein->ResidueAt(0).atom_indices.front();
    const size_t miss_ai = protein->ResidueAt(1).atom_indices.front();

    auto observed = ReadNpy(output_dir / "dssp_observed.npy");
    auto backbone = ReadNpy(output_dir / "dssp_backbone.npy");
    auto ss8 = ReadNpy(output_dir / "dssp_ss8.npy");
    auto ppii = ReadNpy(output_dir / "dssp_ppii.npy");
    auto hb = ReadNpy(output_dir / "dssp_hbond_energy.npy");
    auto chi = ReadNpy(output_dir / "dssp_chi.npy");

    ExpectShapeAndDescr(observed, {N}, "|i1");
    ExpectShapeAndDescr(backbone, {N, 5}, "<f8");
    ExpectShapeAndDescr(ss8, {N, 8}, "<f8");
    ExpectShapeAndDescr(ppii, {N}, "|i1");
    ExpectShapeAndDescr(hb, {N, 4}, "<f8");
    ExpectShapeAndDescr(chi, {N, 12}, "<f8");

    const int8_t* obs = DataAs<int8_t>(observed);
    const double* bb = DataAs<double>(backbone);
    const double* ss = DataAs<double>(ss8);
    const int8_t* pp = DataAs<int8_t>(ppii);
    const double* hbe = DataAs<double>(hb);
    const double* ch = DataAs<double>(chi);

    EXPECT_EQ(obs[obs_ai], 1);
    EXPECT_DOUBLE_EQ(bb[obs_ai * 5 + 0], -0.4);
    EXPECT_DOUBLE_EQ(bb[obs_ai * 5 + 1], 0.7);
    EXPECT_DOUBLE_EQ(bb[obs_ai * 5 + 2], 12.5);
    EXPECT_DOUBLE_EQ(bb[obs_ai * 5 + 3], 1.0);
    EXPECT_DOUBLE_EQ(bb[obs_ai * 5 + 4], 0.0);
    EXPECT_DOUBLE_EQ(ss[obs_ai * 8 + 0], 1.0);
    for (int c = 1; c < 8; ++c)
        EXPECT_DOUBLE_EQ(ss[obs_ai * 8 + c], 0.0);
    EXPECT_EQ(pp[obs_ai], 0);
    EXPECT_DOUBLE_EQ(hbe[obs_ai * 4 + 0], -1.1);
    EXPECT_DOUBLE_EQ(hbe[obs_ai * 4 + 1], -2.2);
    EXPECT_DOUBLE_EQ(hbe[obs_ai * 4 + 2], -3.3);
    EXPECT_DOUBLE_EQ(hbe[obs_ai * 4 + 3], -4.4);

    EXPECT_EQ(obs[miss_ai], 0);
    EXPECT_TRUE(std::isnan(bb[miss_ai * 5 + 0]));
    EXPECT_TRUE(std::isnan(bb[miss_ai * 5 + 1]));
    EXPECT_TRUE(std::isnan(bb[miss_ai * 5 + 2]));
    EXPECT_DOUBLE_EQ(bb[miss_ai * 5 + 3], 0.0);
    EXPECT_DOUBLE_EQ(bb[miss_ai * 5 + 4], 0.0);
    for (int c = 0; c < 8; ++c)
        EXPECT_DOUBLE_EQ(ss[miss_ai * 8 + c], 0.0);
    EXPECT_EQ(pp[miss_ai], -1);
    for (int c = 0; c < 4; ++c)
        EXPECT_TRUE(std::isnan(hbe[miss_ai * 4 + c]));
    for (int c = 0; c < 12; ++c)
        EXPECT_TRUE(std::isfinite(ch[miss_ai * 12 + c]));
}

TEST_F(DsspResultTest, dssp_ppii_npy_preserves_ss8_coil_compatibility) {
    auto& conf = protein->Conformation();
    std::vector<DsspResidue> residues(protein->ResidueCount());
    ASSERT_GE(residues.size(), 2u);
    ASSERT_FALSE(protein->ResidueAt(0).atom_indices.empty());
    ASSERT_FALSE(protein->ResidueAt(1).atom_indices.empty());

    residues[0].observed = true;
    residues[0].secondary_structure = 'P';
    residues[1].observed = false;
    residues[1].secondary_structure = 'C';

    auto dssp = DsspResult::CreateForTesting(std::move(residues));
    const fs::path output_dir = fs::temp_directory_path() /
        ("dssp_ppii_" + std::to_string(::getpid()));
    fs::create_directories(output_dir);
    ASSERT_EQ(dssp->WriteFeatures(conf, output_dir.string()), 11);

    const size_t N = conf.AtomCount();
    const size_t p_ai = protein->ResidueAt(0).atom_indices.front();
    const size_t miss_ai = protein->ResidueAt(1).atom_indices.front();

    auto ppii = ReadNpy(output_dir / "dssp_ppii.npy");
    auto ss8 = ReadNpy(output_dir / "dssp_ss8.npy");
    ExpectShapeAndDescr(ppii, {N}, "|i1");
    ExpectShapeAndDescr(ss8, {N, 8}, "<f8");

    const int8_t* pp = DataAs<int8_t>(ppii);
    const double* ss = DataAs<double>(ss8);
    EXPECT_EQ(pp[p_ai], 1);
    EXPECT_EQ(pp[miss_ai], -1);
    for (int c = 0; c < 8; ++c) {
        EXPECT_DOUBLE_EQ(ss[p_ai * 8 + c], c == 7 ? 1.0 : 0.0);
        EXPECT_DOUBLE_EQ(ss[miss_ai * 8 + c], 0.0);
    }
}


// Project-formula consistency check for the C2/B6 sign conversion.
//
// dssp_backbone.npy columns 0,1 store project-convention phi/psi, produced
// by negating libdssp's plain-IUPAC values. This test recomputes the same
// project-signed atan2 dihedral directly from the 1UBQ coordinates and
// checks the NPY writer's conversion. It is a consistency test, not the
// independent convention authority; the direct-libdssp Ramachandran test
// supplies that external pin.
TEST_F(DsspResultTest, DsspBackboneNpyPhiPsiMatchesProjectConventionGeometry) {
    auto& conf = protein->Conformation();
    auto dssp = DsspResult::Compute(conf);
    ASSERT_NE(dssp, nullptr);

    const fs::path output_dir = fs::temp_directory_path() /
        ("dssp_backbone_project_" + std::to_string(::getpid()));
    fs::create_directories(output_dir);
    ASSERT_EQ(dssp->WriteFeatures(conf, output_dir.string()), 11);

    auto backbone = ReadNpy(output_dir / "dssp_backbone.npy");
    auto observed = ReadNpy(output_dir / "dssp_observed.npy");
    const size_t N = conf.AtomCount();
    ExpectShapeAndDescr(backbone, {N, 5}, "<f8");
    ExpectShapeAndDescr(observed, {N}, "|i1");
    const double* bb = DataAs<double>(backbone);
    const int8_t* obs = DataAs<int8_t>(observed);

    // Project dihedral: D(p1,p2,p3,p4) = atan2((n1 x b2hat).n2, n1.n2).
    auto dihedral = [](const Vec3& p1, const Vec3& p2,
                       const Vec3& p3, const Vec3& p4) -> double {
        const Vec3 b1 = p2 - p1;
        const Vec3 b2 = p3 - p2;
        const Vec3 b3 = p4 - p3;
        const double b2n = b2.norm();
        const Vec3 n1 = b1.cross(b2);
        const Vec3 n2 = b2.cross(b3);
        const Vec3 m1 = n1.cross(b2 / b2n);
        return std::atan2(m1.dot(n2), n1.dot(n2));
    };

    constexpr double kTol = 0.05;       // >> libdssp drift, << sign flip
    constexpr double kMinAngle = 0.3;   // only assert where a flip is stark
    int phi_checked = 0, psi_checked = 0;

    for (size_t ri = 1; ri + 1 < protein->ResidueCount(); ++ri) {
        const Residue& prev = protein->ResidueAt(ri - 1);
        const Residue& res  = protein->ResidueAt(ri);
        const Residue& next = protein->ResidueAt(ri + 1);

        // Same-chain adjacency; backbone atoms present.
        if (res.chain_id != prev.chain_id || res.chain_id != next.chain_id)
            continue;
        if (prev.C == Residue::NONE || res.N == Residue::NONE ||
            res.CA == Residue::NONE || res.C == Residue::NONE ||
            next.N == Residue::NONE)
            continue;

        const size_t row = res.CA;  // dssp_backbone rows are per-atom
        if (row >= N || obs[row] != 1) continue;
        const double phi_npy = bb[row * 5 + 0];
        const double psi_npy = bb[row * 5 + 1];
        if (!std::isfinite(phi_npy) || !std::isfinite(psi_npy)) continue;

        const double phi_geo = dihedral(
            conf.PositionAt(prev.C), conf.PositionAt(res.N),
            conf.PositionAt(res.CA), conf.PositionAt(res.C));
        const double psi_geo = dihedral(
            conf.PositionAt(res.N), conf.PositionAt(res.CA),
            conf.PositionAt(res.C), conf.PositionAt(next.N));

        if (std::isfinite(phi_geo) && std::abs(phi_geo) > kMinAngle) {
            EXPECT_NEAR(phi_npy, phi_geo, kTol)
                << "dssp_backbone phi (col 0) is not project-convention "
                   "geometry at residue " << ri
                << " (a sign flip would land near " << -phi_geo << ")";
            ++phi_checked;
        }
        if (std::isfinite(psi_geo) && std::abs(psi_geo) > kMinAngle) {
            EXPECT_NEAR(psi_npy, psi_geo, kTol)
                << "dssp_backbone psi (col 1) is not project-convention "
                   "geometry at residue " << ri
                << " (a sign flip would land near " << -psi_geo << ")";
            ++psi_checked;
        }
    }

    // Guard against a silent skip-all: 1UBQ is a complete 76-residue
    // crystal structure; the great majority of internal residues qualify.
    EXPECT_GT(phi_checked, 30)
        << "phi pin did not exercise — fixture/observed-mask regression?";
    EXPECT_GT(psi_checked, 30)
        << "psi pin did not exercise — fixture/observed-mask regression?";

    // NB: std::filesystem::remove_all resolves into libtorch's bundled
    // libstdc++ and faults in this process; the test tree removes files
    // individually instead (cf. test_coulomb_result / test_mopac_result).
    std::error_code ec;
    for (const char* f : {"dssp_observed.npy", "dssp_backbone.npy",
                          "dssp_ss8.npy", "dssp_ppii.npy",
                          "dssp_hbond_energy.npy", "dssp_chi.npy"}) {
        fs::remove(output_dir / f, ec);
    }
    fs::remove(output_dir, ec);
}

#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include "PdbFileReader.h"
#include "GeometryResult.h"
#include "DsspResult.h"
#include "PhysicalConstants.h"
#include "NpyWriter.h"
#include "Residue.h"
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
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

}  // namespace


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
    ASSERT_EQ(dssp->WriteFeatures(conf, output_dir.string()), 5);

    const size_t N = conf.AtomCount();
    const size_t obs_ai = protein->ResidueAt(0).atom_indices.front();
    const size_t miss_ai = protein->ResidueAt(1).atom_indices.front();

    auto observed = ReadNpy(output_dir / "dssp_observed.npy");
    auto backbone = ReadNpy(output_dir / "dssp_backbone.npy");
    auto ss8 = ReadNpy(output_dir / "dssp_ss8.npy");
    auto hb = ReadNpy(output_dir / "dssp_hbond_energy.npy");
    auto chi = ReadNpy(output_dir / "dssp_chi.npy");

    ExpectShapeAndDescr(observed, {N}, "|i1");
    ExpectShapeAndDescr(backbone, {N, 5}, "<f8");
    ExpectShapeAndDescr(ss8, {N, 8}, "<f8");
    ExpectShapeAndDescr(hb, {N, 4}, "<f8");
    ExpectShapeAndDescr(chi, {N, 12}, "<f8");

    const int8_t* obs = DataAs<int8_t>(observed);
    const double* bb = DataAs<double>(backbone);
    const double* ss = DataAs<double>(ss8);
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
    for (int c = 0; c < 4; ++c)
        EXPECT_TRUE(std::isnan(hbe[miss_ai * 4 + c]));
    for (int c = 0; c < 12; ++c)
        EXPECT_TRUE(std::isfinite(ch[miss_ai * 12 + c]));
}

#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include "PdbFileReader.h"
#include "GeometryResult.h"
#include "DsspResult.h"
#include "Dssp8TimeSeriesTrajectoryResult.h"
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

}  // namespace

TEST(DsspResultHelpers, PpiiPredicateAndSs8Compatibility) {
    EXPECT_TRUE(DsspIsPpii('P'));
    EXPECT_FALSE(DsspIsPpii('C'));
    EXPECT_EQ(Dssp8TimeSeriesTrajectoryResult::Ss8Code('P'), 7u);
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
    ASSERT_EQ(dssp->WriteFeatures(conf, output_dir.string()), 6);

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
    ASSERT_EQ(dssp->WriteFeatures(conf, output_dir.string()), 6);

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


// Independent forcing function for the C2/B6 sign convention.
//
// dssp_backbone.npy columns 0,1 store IUPAC phi/psi, produced by the
// writer as -libdssp (the value-move under test: the old writer stored
// raw libdssp = negated IUPAC). The ledger's cited pin,
// DihedralTimeSeries CrossResultConsistencyDsspPlanarGeometry, only runs
// with the external gitignored MD trajectory fixture and SKIPS wherever
// that fixture is absent. This test recomputes IUPAC phi/psi from the
// 1UBQ backbone coordinates with an independent atan2 dihedral — the
// standard math definition, evaluated straight from coordinates, NOT
// through the libdssp/libcifpp code path — and asserts the written NPY
// matches. The tolerance absorbs libdssp's algorithmic drift on the same
// coordinates (~1e-3 rad) while a sign flip (O(pi)) fails loudly. Always
// runs: the fixture is the in-tree 1UBQ PDB.
TEST_F(DsspResultTest, DsspBackboneNpyPhiPsiPinnedToIndependentIupacGeometry) {
    auto& conf = protein->Conformation();
    auto dssp = DsspResult::Compute(conf);
    ASSERT_NE(dssp, nullptr);

    const fs::path output_dir = fs::temp_directory_path() /
        ("dssp_backbone_iupac_" + std::to_string(::getpid()));
    fs::create_directories(output_dir);
    ASSERT_EQ(dssp->WriteFeatures(conf, output_dir.string()), 6);

    auto backbone = ReadNpy(output_dir / "dssp_backbone.npy");
    auto observed = ReadNpy(output_dir / "dssp_observed.npy");
    const size_t N = conf.AtomCount();
    ExpectShapeAndDescr(backbone, {N, 5}, "<f8");
    ExpectShapeAndDescr(observed, {N}, "|i1");
    const double* bb = DataAs<double>(backbone);
    const int8_t* obs = DataAs<int8_t>(observed);

    // Independent IUPAC dihedral: D(p1,p2,p3,p4) = atan2((n1 x b2hat).n2,
    // n1.n2). Standard definition; deliberately re-derived here so the
    // pin does not borrow the production result path.
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
                << "dssp_backbone phi (col 0) is not independent IUPAC "
                   "geometry at residue " << ri
                << " (a sign flip would land near " << -phi_geo << ")";
            ++phi_checked;
        }
        if (std::isfinite(psi_geo) && std::abs(psi_geo) > kMinAngle) {
            EXPECT_NEAR(psi_npy, psi_geo, kTol)
                << "dssp_backbone psi (col 1) is not independent IUPAC "
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

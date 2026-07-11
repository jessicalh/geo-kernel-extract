#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iterator>
#include <cstdio>
#include <sstream>

#include "McConnellResult.h"
#include "MopacResult.h"
#include "SidechainCarbonylAnisotropyResult.h"
#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "Protein.h"
#include "PdbFileReader.h"
#include "OrcaRunLoader.h"
#include "PhysicalConstants.h"

#include <filesystem>
#include <Eigen/Geometry>
#include <unistd.h>
namespace fs = std::filesystem;
using namespace nmr;

namespace {

double MaxAbs(const Mat3& m) {
    double v = 0.0;
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c)
            v = std::max(v, std::abs(m(r, c)));
    return v;
}

bool SphericalAllZero(const SphericalTensor& st) {
    if (st.T0 != 0.0) return false;
    for (double v : st.T1) if (v != 0.0) return false;
    for (double v : st.T2) if (v != 0.0) return false;
    return true;
}

std::string ReadText(const fs::path& path) {
    std::ifstream in(path);
    return std::string(std::istreambuf_iterator<char>(in),
                       std::istreambuf_iterator<char>());
}

std::unique_ptr<Protein> BuildSyntheticPeptideCOProtein(
        bool degenerate = false,
        const Mat3& rotation = Mat3::Identity()) {
    auto protein = std::make_unique<Protein>();

    auto atom_n = Atom::Create(Element::N);
    auto atom_c = Atom::Create(Element::C);
    auto atom_o = Atom::Create(Element::O);
    auto atom_h = Atom::Create(Element::H);
    atom_n->pdb_atom_name = "N";
    atom_c->pdb_atom_name = "C";
    atom_o->pdb_atom_name = "O";
    atom_h->pdb_atom_name = "H";
    atom_n->residue_index = 0;
    atom_c->residue_index = 0;
    atom_o->residue_index = 0;
    atom_h->residue_index = 0;

    protein->AddAtom(std::move(atom_n));  // 0: amide N
    protein->AddAtom(std::move(atom_c));  // 1: carbonyl C
    protein->AddAtom(std::move(atom_o));  // 2: carbonyl O
    protein->AddAtom(std::move(atom_h));  // 3: target probe

    Residue res;
    res.type = AminoAcid::ALA;
    res.sequence_number = 1;
    res.chain_id = "A";
    res.atom_indices = {0, 1, 2, 3};
    res.N = 0;
    res.C = 1;
    res.O = 2;
    res.H = 3;
    protein->AddResidue(res);

    std::vector<Vec3> positions = {
        degenerate ? Vec3(-1.45, 0.0, 0.0) : Vec3(0.0, 1.45, 0.0),
        Vec3(0.0, 0.0, 0.0),
        Vec3(1.2, 0.0, 0.0),
        Vec3(0.6, 3.0, 0.0)
    };
    for (Vec3& pos : positions) pos = rotation * pos;

    protein->FinalizeConstruction(positions);
    protein->AddCrystalConformation(positions, 0, 0, 0, "synthetic_peptide_co");
    return protein;
}

std::unique_ptr<Protein> BuildSyntheticXHCategoryProtein() {
    auto protein = std::make_unique<Protein>();
    for (int ri = 0; ri < 3; ++ri) {
        Residue res;
        res.type = AminoAcid::Unknown;
        res.sequence_number = ri + 1;
        res.chain_id = "A";
        protein->AddResidue(std::move(res));
    }

    std::vector<Vec3> positions;
    auto add = [&](size_t ri, Element element, Vec3 pos) {
        auto atom = Atom::Create(element);
        atom->residue_index = ri;
        const size_t ai = protein->AddAtom(std::move(atom));
        protein->MutableResidueAt(ri).atom_indices.push_back(ai);
        positions.push_back(pos);
        return ai;
    };

    const size_t backbone_n = add(0, Element::N, Vec3(0.0, 0.0, 0.0));
    add(0, Element::H, Vec3(1.0, 0.0, 0.0));
    add(1, Element::C, Vec3(10.0, 0.0, 0.0));
    add(1, Element::H, Vec3(11.0, 0.0, 0.0));
    add(2, Element::S, Vec3(20.0, 0.0, 0.0));
    add(2, Element::H, Vec3(21.3, 0.0, 0.0));
    protein->MutableResidueAt(0).N = backbone_n;

    protein->FinalizeConstruction(positions);
    return protein;
}

std::unique_ptr<Protein> BuildSyntheticSidechainCOProtein() {
    auto protein = std::make_unique<Protein>();
    std::vector<Vec3> positions;

    auto add_residue = [&](AminoAcid type, int sequence_number) {
        Residue residue;
        residue.type = type;
        residue.sequence_number = sequence_number;
        residue.chain_id = "A";
        return protein->AddResidue(std::move(residue));
    };
    auto add_atom = [&](std::size_t residue_index, const char* name,
                        Element element, const Vec3& position) {
        auto atom = Atom::Create(element);
        atom->pdb_atom_name = name;
        atom->residue_index = residue_index;
        const std::size_t atom_index = protein->AddAtom(std::move(atom));
        protein->MutableResidueAt(residue_index).atom_indices.push_back(
            atom_index);
        positions.push_back(position);
        return atom_index;
    };

    // ASN primary amide: one typed SidechainCO source.  The amide N is the
    // second same-planar-group atom that fixes the local plane.
    const std::size_t asn = add_residue(AminoAcid::ASN, 1);
    add_atom(asn, "CG", Element::C, Vec3(0.0, 0.0, 0.0));
    add_atom(asn, "OD1", Element::O, Vec3(1.23, 0.0, 0.0));
    add_atom(asn, "ND2", Element::N, Vec3(-0.60, 1.05, 0.0));

    // ASP carboxylate: two typed SidechainCO sources, each using the other
    // oxygen as its deterministic in-plane reference.
    const std::size_t asp = add_residue(AminoAcid::ASP, 2);
    add_atom(asp, "CG", Element::C, Vec3(5.0, 0.0, 0.0));
    add_atom(asp, "OD1", Element::O, Vec3(6.25, 0.0, 0.0));
    add_atom(asp, "OD2", Element::O, Vec3(4.375, 1.08253175, 0.0));

    // A real peptide C=O in the same synthetic topology is a negative
    // control: source enumeration must exclude it by typed BondCategory.
    const std::size_t ala = add_residue(AminoAcid::ALA, 3);
    const std::size_t ala_n =
        add_atom(ala, "N", Element::N, Vec3(12.0, 1.35, 0.0));
    const std::size_t ala_c =
        add_atom(ala, "C", Element::C, Vec3(12.0, 0.0, 0.0));
    const std::size_t ala_o =
        add_atom(ala, "O", Element::O, Vec3(13.20, 0.0, 0.0));
    protein->MutableResidueAt(ala).N = ala_n;
    protein->MutableResidueAt(ala).C = ala_c;
    protein->MutableResidueAt(ala).O = ala_o;

    protein->FinalizeConstruction(positions);
    protein->AddCrystalConformation(
        positions, 0.0, 0.0, 0.0, "synthetic_sidechain_co");
    return protein;
}

template<typename T>
struct NpyArray {
    std::vector<std::size_t> shape;
    std::vector<T> values;
};

std::string TrimNpyToken(std::string token) {
    auto is_space = [](unsigned char c) { return std::isspace(c) != 0; };
    token.erase(token.begin(),
                std::find_if(token.begin(), token.end(),
                             [&](char c) { return !is_space(c); }));
    token.erase(std::find_if(token.rbegin(), token.rend(),
                             [&](char c) { return !is_space(c); }).base(),
                token.end());
    return token;
}

template<typename T>
NpyArray<T> ReadNpy(const fs::path& path, const char* dtype) {
    std::ifstream in(path, std::ios::binary);
    EXPECT_TRUE(in.is_open()) << path;
    NpyArray<T> result;
    if (!in.is_open()) return result;

    char magic[6] = {};
    in.read(magic, 6);
    EXPECT_EQ(std::string(magic, 6), std::string("\x93NUMPY", 6));
    char version[2] = {};
    in.read(version, 2);
    EXPECT_EQ(version[0], 1);
    EXPECT_EQ(version[1], 0);
    std::uint16_t header_length = 0;
    in.read(reinterpret_cast<char*>(&header_length), sizeof(header_length));
    std::string header(header_length, '\0');
    in.read(header.data(), header_length);
    EXPECT_NE(header.find(std::string("'descr': '") + dtype + "'"),
              std::string::npos);

    const std::size_t shape_begin = header.find('(');
    const std::size_t shape_end = header.find(')', shape_begin);
    EXPECT_NE(shape_begin, std::string::npos);
    EXPECT_NE(shape_end, std::string::npos);
    if (shape_begin == std::string::npos || shape_end == std::string::npos)
        return result;

    std::stringstream shape_stream(
        header.substr(shape_begin + 1, shape_end - shape_begin - 1));
    std::string token;
    std::size_t value_count = 1;
    while (std::getline(shape_stream, token, ',')) {
        token = TrimNpyToken(std::move(token));
        if (token.empty()) continue;
        const std::size_t extent = static_cast<std::size_t>(
            std::stoull(token));
        result.shape.push_back(extent);
        value_count *= extent;
    }
    if (result.shape.empty()) value_count = 0;
    result.values.resize(value_count);
    if (value_count > 0) {
        in.read(reinterpret_cast<char*>(result.values.data()),
                static_cast<std::streamsize>(value_count * sizeof(T)));
    }
    return result;
}

fs::path SidechainCOTempDir(const char* stem) {
    const fs::path dir = fs::temp_directory_path() /
        (std::string(stem) + "_" + std::to_string(::getpid()));
    fs::create_directories(dir);
    return dir;
}

void RemoveSidechainCOOutputs(const fs::path& dir) {
    constexpr std::array<const char*, 6> files = {
        "sidechain_co_source_bonds.npy",
        "sidechain_co_frame.npy",
        "sidechain_co_frame_quality.npy",
        "sidechain_co_fixed_T2.npy",
        "sidechain_co_bo_T2.npy",
        "sidechain_co_scalar_audit.npy",
    };
    for (const char* file : files) {
        std::remove((dir / file).string().c_str());
    }
    ::rmdir(dir.string().c_str());
}

size_t FirstBondWithCategory(const Protein& protein, BondCategory category) {
    for (size_t bi = 0; bi < protein.BondCount(); ++bi) {
        if (protein.BondAt(bi).category == category)
            return bi;
    }
    return SIZE_MAX;
}

size_t FirstBondBetween(const Protein& protein, size_t a, size_t b) {
    for (size_t bi = 0; bi < protein.BondCount(); ++bi) {
        const Bond& bond = protein.BondAt(bi);
        if ((bond.atom_index_a == a && bond.atom_index_b == b) ||
            (bond.atom_index_a == b && bond.atom_index_b == a))
            return bi;
    }
    return SIZE_MAX;
}

}  // namespace



// ============================================================================
// Analytical test: known geometry, hand-calculable result
// ============================================================================

TEST(McConnellAnalytical, DipolarKernelAtKnownGeometry) {
    // Atom at (3, 0, 0). Bond midpoint at origin. Bond direction (0, 0, 1).
    //
    // d = (3, 0, 0), r = 3, d_hat = (1, 0, 0)
    // b_hat = (0, 0, 1)
    // cos_theta = d_hat . b_hat = 0
    //
    // Unit axial source shape:
    //   Qhat = u u^T - I/3 = diag(-1/3, -1/3, 2/3)
    // Clean McConnell response:
    //   A = D(r) Qhat
    // PCS scalar branch:
    //   T0(A) = trace(A)/3 = n^T Qhat n / r^3 = (-1/3)/27 = -1/81
    //
    // Symmetric dipolar kernel K:
    //   K_11 = (3*1*1 - 1) / 27 =  2/27 =  0.07407
    //   K_22 = (3*0*0 - 1) / 27 = -1/27 = -0.03704
    //   K_33 = (3*0*0 - 1) / 27 = -1/27 = -0.03704
    //   Off-diagonal = 0 (d_hat only has x component)
    //   Trace = 0 (traceless)
    //
    // Response A = D Qhat:
    //   diag(-2/81, 1/81, -2/81), so trace(A)/3 = -1/81.

    // Build a minimal protein with one bond
    auto protein = std::make_unique<Protein>();

    auto atom_a = Atom::Create(Element::C);
    auto atom_b = Atom::Create(Element::O);
    auto atom_field = Atom::Create(Element::H);
    atom_a->residue_index = 0;
    atom_b->residue_index = 0;
    atom_field->residue_index = 0;

    protein->AddAtom(std::move(atom_a));  // 0: bond endpoint C at origin
    protein->AddAtom(std::move(atom_b));  // 1: bond endpoint O at (0,0,1)
    protein->AddAtom(std::move(atom_field));  // 2: field point H at (3,0,0)

    Residue res;
    res.type = AminoAcid::ALA;
    res.sequence_number = 1;
    res.chain_id = "A";
    res.atom_indices = {0, 1, 2};
    protein->AddResidue(res);

    // Positions: C at (-0.5,0,0), O at (0.5,0,0) gives midpoint at origin,
    // but we want bond direction (0,0,1). So:
    // C at (0,0,-0.5), O at (0,0,0.5), midpoint at origin, direction = (0,0,1)
    std::vector<Vec3> positions = {
        Vec3(0.0, 0.0, -0.5),   // C
        Vec3(0.0, 0.0,  0.5),   // O
        Vec3(3.0, 0.0,  0.0)    // H (field point)
    };

    protein->FinalizeConstruction(positions);

    // Manually set bond category for testing
    // (FinalizeConstruction may not classify correctly for this toy geometry)

    auto& conf = protein->AddCrystalConformation(positions, 0, 0, 0, "test");

    conf.AttachResult(GeometryResult::Compute(conf));
    conf.AttachResult(SpatialIndexResult::Compute(conf));

    auto mc = McConnellResult::Compute(conf);
    ASSERT_NE(mc, nullptr);
    conf.AttachResult(std::move(mc));

    // Check the field point atom (index 2)
    const auto& ca = conf.AtomAt(2);

    // Should have at least one bond neighbour
    ASSERT_GT(ca.bond_neighbours.size(), 0u)
        << "Field atom should see at least one bond";

    // Find the bond to the C-O pair
    const BondNeighbourhood& bn = ca.bond_neighbours[0];

    // Distance should be 3.0 A (from midpoint at origin to (3,0,0))
    EXPECT_NEAR(bn.distance_to_midpoint, 3.0, 0.01)
        << "Distance to bond midpoint should be 3.0 A";

    // PCS scalar: n^T Qhat n / r^3 = -1/81.
    EXPECT_NEAR(bn.mcconnell_scalar, -1.0/81.0, 1e-6)
        << "McConnell PCS scalar should be -1/81";

    // Dipolar kernel K should be traceless
    double trace_K = bn.dipolar_tensor.trace();
    EXPECT_NEAR(trace_K, 0.0, 1e-10) << "Dipolar kernel must be traceless";

    // K diagonal: (2/27, -1/27, -1/27)
    EXPECT_NEAR(bn.dipolar_tensor(0,0),  2.0/27.0, 1e-6);
    EXPECT_NEAR(bn.dipolar_tensor(1,1), -1.0/27.0, 1e-6);
    EXPECT_NEAR(bn.dipolar_tensor(2,2), -1.0/27.0, 1e-6);

    // K off-diagonal: all zero (d_hat = (1,0,0))
    EXPECT_NEAR(bn.dipolar_tensor(0,1), 0.0, 1e-10);
    EXPECT_NEAR(bn.dipolar_tensor(0,2), 0.0, 1e-10);
    EXPECT_NEAR(bn.dipolar_tensor(1,2), 0.0, 1e-10);

    EXPECT_NEAR(ca.mc_shielding_contribution.T0, -1.0/81.0, 1e-6)
        << "T0 of D(r)Qhat should equal the PCS scalar";

    std::cout << "  Analytical test passed:\n"
              << "    PCS scalar = " << bn.mcconnell_scalar
              << " (expected -1/81 = " << -1.0/81.0 << ")\n"
              << "    T0 = " << ca.mc_shielding_contribution.T0 << "\n"
              << "    K trace = " << trace_K << "\n"
              << "    K diag = (" << bn.dipolar_tensor(0,0) << ", "
              << bn.dipolar_tensor(1,1) << ", "
              << bn.dipolar_tensor(2,2) << ")\n";
}


// ============================================================================
// Full protein test
// ============================================================================

class McConnellProteinTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!fs::exists(nmr::test::TestEnvironment::UbqProtonated())) {
            GTEST_SKIP() << "1UBQ.pdb not found";
        }
        auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
        if (!r.Ok()) GTEST_SKIP() << "Failed to load 1UBQ";
        protein = std::move(r.protein);

        auto& conf = protein->Conformation();
        conf.AttachResult(GeometryResult::Compute(conf));
        conf.AttachResult(SpatialIndexResult::Compute(conf));
    }

    std::unique_ptr<Protein> protein;
};


TEST_F(McConnellProteinTest, ComputeAndAttach) {
    auto& conf = protein->Conformation();
    auto mc = McConnellResult::Compute(conf);
    ASSERT_NE(mc, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(mc)));
    ASSERT_TRUE(conf.HasResult<McConnellResult>());
}


TEST_F(McConnellProteinTest, EveryAtomHasBondNeighbours) {
    auto& conf = protein->Conformation();
    conf.AttachResult(McConnellResult::Compute(conf));

    int with_neighbours = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        if (!conf.AtomAt(ai).bond_neighbours.empty())
            with_neighbours++;
    }

    // Every atom in a protein should see nearby bonds
    EXPECT_GT(with_neighbours, static_cast<int>(conf.AtomCount()) - 5)
        << "Almost every atom should have bond neighbours";

    std::cout << "  Atoms with bond neighbours: " << with_neighbours
              << " / " << conf.AtomCount() << "\n";
}


TEST_F(McConnellProteinTest, DipolarKernelsAreTraceless) {
    auto& conf = protein->Conformation();
    conf.AttachResult(McConnellResult::Compute(conf));

    double max_trace = 0.0;
    int checked = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& bn : conf.AtomAt(ai).bond_neighbours) {
            double trace = bn.dipolar_tensor.trace();
            max_trace = std::max(max_trace, std::abs(trace));
            checked++;
        }
    }

    EXPECT_LT(max_trace, 1e-10)
        << "All dipolar kernels K must be traceless";

    std::cout << "  Checked " << checked << " kernels, max |trace| = "
              << max_trace << "\n";
}


TEST_F(McConnellProteinTest, McConnellScalarNonZero) {
    auto& conf = protein->Conformation();
    conf.AttachResult(McConnellResult::Compute(conf));

    // Most atoms should have non-zero CO sum
    int with_co = 0;
    double max_co = 0.0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        double co = std::abs(conf.AtomAt(ai).mcconnell_co_sum);
        if (co > 1e-6) with_co++;
        max_co = std::max(max_co, co);
    }

    EXPECT_GT(with_co, static_cast<int>(conf.AtomCount() / 2))
        << "Most atoms should have non-zero CO McConnell sum";

    std::cout << "  Atoms with |CO sum| > 1e-6: " << with_co
              << " / " << conf.AtomCount()
              << ", max |CO sum| = " << max_co << "\n";
}


TEST_F(McConnellProteinTest, ShieldingContributionHasT0AndT2) {
    auto& conf = protein->Conformation();
    conf.AttachResult(McConnellResult::Compute(conf));

    int nonzero_t0 = 0, nonzero_t2 = 0;
    double max_t0 = 0, max_t2 = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& sc = conf.AtomAt(ai).mc_shielding_contribution;
        if (std::abs(sc.T0) > 1e-8) nonzero_t0++;
        double t2mag = sc.T2Magnitude();
        if (t2mag > 1e-8) nonzero_t2++;
        max_t0 = std::max(max_t0, std::abs(sc.T0));
        max_t2 = std::max(max_t2, t2mag);
    }

    EXPECT_GT(nonzero_t0, 0) << "Full McConnell tensor must have non-zero T0";
    EXPECT_GT(nonzero_t2, 0) << "Full McConnell tensor must have non-zero T2";

    std::cout << "  T0 nonzero: " << nonzero_t0 << ", max |T0| = " << max_t0 << "\n";
    std::cout << "  T2 nonzero: " << nonzero_t2 << ", max |T2| = " << max_t2 << "\n";
}


TEST_F(McConnellProteinTest, PCSScalarMatchesTraceBranchForSingleBond) {
    auto& conf = protein->Conformation();
    conf.AttachResult(McConnellResult::Compute(conf));

    int checked = 0;
    double max_diff = 0.0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& bn : conf.AtomAt(ai).bond_neighbours) {
            Vec3 atom_pos = conf.PositionAt(ai);
            Vec3 midpoint = conf.bond_midpoints[bn.bond_index];
            Vec3 axis = conf.bond_directions[bn.bond_index];
            auto kernel = McConnellResult::ComputePairKernel(
                atom_pos, midpoint, axis);
            double r = kernel.distance;
            if (r < MIN_DISTANCE) continue;

            double r3 = r * r * r;
            double pcs = kernel.direction.dot(
                kernel.source_shape * kernel.direction) / r3;
            double t0 = SphericalTensor::Decompose(kernel.response).T0;
            max_diff = std::max(max_diff, std::abs(t0 - bn.mcconnell_scalar));
            max_diff = std::max(max_diff, std::abs(t0 - pcs));
            checked++;
        }
    }

    EXPECT_LT(max_diff, 1e-10)
        << "T0(DQhat) must equal the PCS scalar n^T Qhat n/r^3";

    std::cout << "  Checked " << checked
              << " pairs, max PCS scalar branch diff = "
              << max_diff << "\n";
}


TEST(McConnellImplementationChecks, RotationEquivariance) {
    Vec3 source_center(0.3, -0.7, 1.1);
    Vec3 source_axis(0.4, 0.8, -0.2);
    Vec3 target(3.2, 1.4, -0.6);

    const auto k = McConnellResult::ComputePairKernel(
        target, source_center, source_axis);
    const SphericalTensor st = SphericalTensor::Decompose(k.response);

    const Mat3 R = Eigen::AngleAxisd(
        0.73, Vec3(0.2, -0.4, 0.9).normalized()).toRotationMatrix();
    const auto kr = McConnellResult::ComputePairKernel(
        R * target, R * source_center, R * source_axis);
    const SphericalTensor sr = SphericalTensor::Decompose(kr.response);

    EXPECT_NEAR(sr.T0, st.T0, 1e-12);

    Vec3 t1(st.T1[0], st.T1[1], st.T1[2]);
    Vec3 t1_rot(sr.T1[0], sr.T1[1], sr.T1[2]);
    EXPECT_LT((t1_rot - R * t1).norm(), 1e-12);

    SphericalTensor st_t2 = st;
    st_t2.T0 = 0.0;
    st_t2.T1 = {0.0, 0.0, 0.0};
    SphericalTensor sr_t2 = sr;
    sr_t2.T0 = 0.0;
    sr_t2.T1 = {0.0, 0.0, 0.0};
    const Mat3 expected_t2 = R * st_t2.Reconstruct() * R.transpose();
    EXPECT_LT(MaxAbs(sr_t2.Reconstruct() - expected_t2), 1e-12);

    EXPECT_LT(MaxAbs(kr.response - R * k.response * R.transpose()), 1e-12);
}


TEST(McConnellImplementationChecks, PairPCSScalarIdentity) {
    const auto k = McConnellResult::ComputePairKernel(
        Vec3(2.0, -1.0, 3.0),
        Vec3(-0.5, 0.25, 0.75),
        Vec3(0.3, 0.7, 0.2));
    ASSERT_GT(k.distance, 0.0);
    const double r3 = k.distance * k.distance * k.distance;
    const double pcs = k.direction.dot(k.source_shape * k.direction) / r3;
    const double t0 = SphericalTensor::Decompose(k.response).T0;
    EXPECT_NEAR(t0, k.response.trace() / 3.0, 1e-15);
    EXPECT_NEAR(t0, pcs, 1e-15);
}


TEST(McConnellImplementationChecks, PeptideCORhombicSourceShapePinnedAnalytic) {
    auto protein = BuildSyntheticPeptideCOProtein();
    auto& conf = protein->Conformation();
    conf.AttachResult(GeometryResult::Compute(conf));
    const size_t co_bond = FirstBondWithCategory(*protein, BondCategory::PeptideCO);
    ASSERT_NE(co_bond, SIZE_MAX);

    const auto k = McConnellResult::ComputePeptideCORhombicPairKernel(
        conf, co_bond, conf.PositionAt(3));

    // C->O = +x, plane normal = +z, e_in = +y.  From
    // (-5.4, +4.0, -14) the mean-subtracted principal values are
    // (-4/15, 137/15, -133/15); dividing by axial scale 137/10 gives
    // (out, para, in) = (-8/411, 2/3, -266/411).
    Mat3 expected = Mat3::Zero();
    expected(0, 0) = 2.0 / 3.0;
    expected(1, 1) = -266.0 / 411.0;
    expected(2, 2) = -8.0 / 411.0;

    EXPECT_LT(MaxAbs(k.source_shape - expected), 1e-12);
    EXPECT_NEAR(k.source_shape.trace(), 0.0, 1e-15);
}


TEST(McConnellImplementationChecks, PeptideCORhombicSourceShapeRotationEquivariance) {
    const Mat3 R = Eigen::AngleAxisd(
        0.91, Vec3(0.3, -0.7, 0.2).normalized()).toRotationMatrix();

    auto protein = BuildSyntheticPeptideCOProtein();
    auto rotated_protein = BuildSyntheticPeptideCOProtein(false, R);
    auto& conf = protein->Conformation();
    auto& rotated_conf = rotated_protein->Conformation();
    conf.AttachResult(GeometryResult::Compute(conf));
    rotated_conf.AttachResult(GeometryResult::Compute(rotated_conf));
    const size_t co_bond = FirstBondWithCategory(*protein, BondCategory::PeptideCO);
    const size_t rotated_co_bond =
        FirstBondWithCategory(*rotated_protein, BondCategory::PeptideCO);
    ASSERT_NE(co_bond, SIZE_MAX);
    ASSERT_NE(rotated_co_bond, SIZE_MAX);

    const auto k = McConnellResult::ComputePeptideCORhombicPairKernel(
        conf, co_bond, conf.PositionAt(3));
    const auto kr = McConnellResult::ComputePeptideCORhombicPairKernel(
        rotated_conf, rotated_co_bond, rotated_conf.PositionAt(3));

    EXPECT_LT(MaxAbs(kr.source_shape -
                     R * k.source_shape * R.transpose()), 1e-12);
    EXPECT_LT(MaxAbs(kr.response -
                     R * k.response * R.transpose()), 1e-12);
}


TEST(McConnellImplementationChecks, PeptideCORhombicContributionNonZeroAndTraceless) {
    auto protein = BuildSyntheticPeptideCOProtein();
    auto& conf = protein->Conformation();
    conf.AttachResult(GeometryResult::Compute(conf));
    conf.AttachResult(SpatialIndexResult::Compute(conf));
    conf.AttachResult(McConnellResult::Compute(conf));

    const size_t co_bond = FirstBondWithCategory(*protein, BondCategory::PeptideCO);
    ASSERT_NE(co_bond, SIZE_MAX);
    const auto k = McConnellResult::ComputePeptideCORhombicPairKernel(
        conf, co_bond, conf.PositionAt(3));
    EXPECT_NEAR(k.source_shape.trace(), 0.0, 1e-15);

    // The emitted rhombic array is the additive delta relative to the
    // legacy axial PeptideCO shape.  Here n = e_in = +y, r = 3 A, so
    // D = diag(-1/27, 2/27, -1/27) and
    // Q_delta = diag(0, -43/137, +43/137).
    Mat3 expected_delta = Mat3::Zero();
    expected_delta(1, 1) = -86.0 / 3699.0;
    expected_delta(2, 2) = -43.0 / 3699.0;

    const Mat3 actual =
        conf.AtomAt(3).mcconnell_peptide_co_rhombic.Reconstruct();
    EXPECT_GT(MaxAbs(actual), 1e-6);
    EXPECT_LT(MaxAbs(actual - expected_delta), 1e-12);

    // C33: the canonical fixed PeptideCO channel now contains the full
    // rhombic response.  The separately emitted audit remains exactly the
    // unweighted rhombic-minus-axial delta.
    const auto axial = McConnellResult::ComputePairKernel(
        conf.PositionAt(3), conf.bond_midpoints[co_bond],
        conf.bond_directions[co_bond]);
    const auto rhombic = McConnellResult::ComputePeptideCORhombicPairKernel(
        conf, co_bond, conf.PositionAt(3));
    const size_t cat = static_cast<size_t>(McConnellSourceCategory::PeptideCO);
    const Mat3 fixed = conf.AtomAt(3)
        .mcconnell_source_tensors[cat]
                                   [static_cast<size_t>(McConnellChannel::Fixed)]
        .Reconstruct();
    EXPECT_LT(MaxAbs(fixed - rhombic.response), 1e-12);
    EXPECT_LT(MaxAbs(fixed - (axial.response + actual)), 1e-12);

    // The same production selector owns the BO channel.  Pin a non-unit BO
    // so an axial-only regression cannot hide behind fixed-channel coverage.
    const auto channels = mcconnell_result_detail::SelectChannelResponses(
        McConnellSourceCategory::PeptideCO,
        axial.response, rhombic.response, 1.75);
    EXPECT_LT(MaxAbs(channels.fixed - rhombic.response), 1e-12);
    EXPECT_LT(MaxAbs(channels.bond_order - 1.75 * rhombic.response), 1e-12);
    EXPECT_LT(MaxAbs(channels.rhombic_audit - actual), 1e-12);
}


TEST(McConnellImplementationChecks, XHBondsUseDedicatedProductionCategories) {
    auto protein = BuildSyntheticXHCategoryProtein();
    ASSERT_EQ(protein->BondCount(), 3u);

    std::array<McConnellSourceCategory, 3> categories{};
    for (size_t bi = 0; bi < protein->BondCount(); ++bi) {
        ASSERT_TRUE(mcconnell_result_detail::IsXHBond(
            *protein, protein->BondAt(bi)));
        categories[bi] = mcconnell_result_detail::SourceCategory(
            *protein, protein->BondAt(bi));
    }

    // Bonds are emitted in atom-pair order by CovalentTopology::Resolve.
    EXPECT_EQ(categories[0], McConnellSourceCategory::BackboneXH);
    EXPECT_EQ(categories[1], McConnellSourceCategory::SidechainXH);
    EXPECT_EQ(categories[2], McConnellSourceCategory::SH);
    EXPECT_STREQ(McConnellSourceCategoryStem(categories[0]), "backbone_xh");
    EXPECT_STREQ(McConnellSourceCategoryStem(categories[1]), "sidechain_xh");
    EXPECT_STREQ(McConnellSourceCategoryStem(categories[2]), "s_h");
}


TEST(SidechainCarbonylAnisotropyProduction,
     TypedSourcesFramesAndNoMopacReadBack) {
    using namespace sidechain_carbonyl_anisotropy_detail;

    auto protein = BuildSyntheticSidechainCOProtein();
    auto& conf = protein->Conformation();
    ASSERT_TRUE(conf.AttachResult(GeometryResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(SpatialIndexResult::Compute(conf)));

    std::vector<SourceBond> expected_sources;
    std::size_t peptide_co_count = 0;
    for (std::size_t bond_index = 0;
         bond_index < protein->BondCount(); ++bond_index) {
        const Bond& bond = protein->BondAt(bond_index);
        if (bond.category == BondCategory::PeptideCO) ++peptide_co_count;
        if (bond.category == BondCategory::SidechainCO) {
            expected_sources.push_back(
                ClassifySourceBond(*protein, bond_index));
        }
    }
    ASSERT_EQ(peptide_co_count, 1u);
    ASSERT_EQ(expected_sources.size(), 3u);

    std::size_t amide_sources = 0;
    std::size_t carboxylate_sources = 0;
    for (const SourceBond& source : expected_sources) {
        EXPECT_TRUE(source.source_valid);
        EXPECT_EQ(source.bond_category, BondCategory::SidechainCO);
        EXPECT_NE(source.plane_reference_atom, SIZE_MAX);
        if (source.oxygen_semantic_class ==
            OxygenSemanticClass::SidechainAmide) {
            ++amide_sources;
            EXPECT_EQ(source.planar_group_kind,
                      PlanarGroupKind::SidechainAmide);
        } else if (source.oxygen_semantic_class ==
                   OxygenSemanticClass::SidechainCarboxylate) {
            ++carboxylate_sources;
            EXPECT_EQ(source.planar_group_kind,
                      PlanarGroupKind::Carboxylate);
        }

        const SourceFrame frame = BuildSourceFrame(conf, source);
        EXPECT_TRUE(frame.frame_valid);
        EXPECT_TRUE(frame.origin.allFinite());
        EXPECT_NEAR(frame.x_axis.norm(), 1.0, 1e-14);
        EXPECT_NEAR(frame.y_axis.norm(), 1.0, 1e-14);
        EXPECT_NEAR(frame.z_axis.norm(), 1.0, 1e-14);
        EXPECT_NEAR(frame.x_axis.dot(frame.y_axis), 0.0, 1e-14);
        EXPECT_NEAR(frame.x_axis.dot(frame.z_axis), 0.0, 1e-14);
        EXPECT_NEAR(frame.y_axis.dot(frame.z_axis), 0.0, 1e-14);
        EXPECT_NEAR(frame.z_axis.cross(frame.x_axis).dot(frame.y_axis),
                    1.0, 1e-14);
        const Vec3 expected_x =
            (conf.PositionAt(source.oxygen_atom) -
             conf.PositionAt(source.carbon_atom)).normalized();
        EXPECT_NEAR(frame.x_axis.dot(expected_x), 1.0, 1e-14);
        EXPECT_LT(frame.orthogonality_error, 1e-14);
        EXPECT_GT(frame.normal_norm, 0.0);
    }
    EXPECT_EQ(amide_sources, 1u);
    EXPECT_EQ(carboxylate_sources, 2u);

    ASSERT_TRUE(conf.AttachResult(McConnellResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(
        SidechainCarbonylAnisotropyResult::Compute(conf)));

    const fs::path output_dir =
        SidechainCOTempDir("sidechain_co_no_mopac");
    ASSERT_EQ(conf.Result<SidechainCarbonylAnisotropyResult>()
                  .WriteFeatures(conf, output_dir.string()),
              6);

    const auto source_rows = ReadNpy<std::int32_t>(
        output_dir / "sidechain_co_source_bonds.npy", "<i4");
    const auto frame_rows = ReadNpy<double>(
        output_dir / "sidechain_co_frame.npy", "<f8");
    const auto quality_rows = ReadNpy<double>(
        output_dir / "sidechain_co_frame_quality.npy", "<f8");
    const auto fixed_rows = ReadNpy<double>(
        output_dir / "sidechain_co_fixed_T2.npy", "<f8");
    const auto bo_rows = ReadNpy<double>(
        output_dir / "sidechain_co_bo_T2.npy", "<f8");
    const auto audit_rows = ReadNpy<double>(
        output_dir / "sidechain_co_scalar_audit.npy", "<f8");

    const std::size_t Q = expected_sources.size();
    const std::size_t N = conf.AtomCount();
    ASSERT_EQ(source_rows.shape, (std::vector<std::size_t>{Q, 8}));
    ASSERT_EQ(frame_rows.shape, (std::vector<std::size_t>{Q, 12}));
    ASSERT_EQ(quality_rows.shape, (std::vector<std::size_t>{Q, 4}));
    ASSERT_EQ(fixed_rows.shape, (std::vector<std::size_t>{N, 9}));
    ASSERT_EQ(bo_rows.shape, (std::vector<std::size_t>{N, 9}));
    ASSERT_EQ(audit_rows.shape, (std::vector<std::size_t>{N, 4}));

    for (std::size_t row = 0; row < Q; ++row) {
        const SourceBond& source = expected_sources[row];
        EXPECT_EQ(source_rows.values[row * 8 + 0],
                  static_cast<std::int32_t>(source.bond_index));
        EXPECT_EQ(source_rows.values[row * 8 + 1],
                  static_cast<std::int32_t>(source.carbon_atom));
        EXPECT_EQ(source_rows.values[row * 8 + 2],
                  static_cast<std::int32_t>(source.oxygen_atom));
        EXPECT_EQ(source_rows.values[row * 8 + 3],
                  static_cast<std::int32_t>(source.residue_index));
        EXPECT_EQ(source_rows.values[row * 8 + 4],
                  static_cast<std::int32_t>(BondCategory::SidechainCO));
        EXPECT_NE(source_rows.values[row * 8 + 4],
                  static_cast<std::int32_t>(BondCategory::PeptideCO));
        EXPECT_EQ(source_rows.values[row * 8 + 5],
                  static_cast<std::int32_t>(source.planar_group_kind));
        EXPECT_EQ(source_rows.values[row * 8 + 6],
                  static_cast<std::int32_t>(
                      source.oxygen_semantic_class));
        EXPECT_EQ(source_rows.values[row * 8 + 7], 1);

        const SourceFrame production_frame =
            BuildSourceFrame(conf, source);
        const std::array<double, 12> packed_frame = {
            production_frame.origin.x(), production_frame.origin.y(),
            production_frame.origin.z(), production_frame.x_axis.x(),
            production_frame.x_axis.y(), production_frame.x_axis.z(),
            production_frame.y_axis.x(), production_frame.y_axis.y(),
            production_frame.y_axis.z(), production_frame.z_axis.x(),
            production_frame.z_axis.y(), production_frame.z_axis.z(),
        };
        for (std::size_t component = 0; component < 12; ++component) {
            EXPECT_DOUBLE_EQ(frame_rows.values[row * 12 + component],
                             packed_frame[component]);
        }
        EXPECT_DOUBLE_EQ(quality_rows.values[row * 4 + 0],
                         production_frame.bond_length);
        EXPECT_DOUBLE_EQ(quality_rows.values[row * 4 + 1],
                         production_frame.orthogonality_error);
        EXPECT_DOUBLE_EQ(quality_rows.values[row * 4 + 2],
                         production_frame.normal_norm);
        EXPECT_DOUBLE_EQ(quality_rows.values[row * 4 + 3], 1.0);
    }

    bool any_nonzero_fixed = false;
    const std::size_t sidechain_category = static_cast<std::size_t>(
        McConnellSourceCategory::SidechainCO);
    const std::size_t fixed_channel = static_cast<std::size_t>(
        McConnellChannel::Fixed);
    for (std::size_t atom_index = 0; atom_index < N; ++atom_index) {
        std::array<double, 9> expected_fixed{};
        const SphericalTensor& tensor = conf.AtomAt(atom_index)
            .mcconnell_source_tensors[sidechain_category][fixed_channel];
        tensor.PackFull9(expected_fixed.data());
        for (std::size_t component = 0; component < 9; ++component) {
            EXPECT_DOUBLE_EQ(fixed_rows.values[atom_index * 9 + component],
                             expected_fixed[component]);
            EXPECT_TRUE(std::isnan(
                bo_rows.values[atom_index * 9 + component]));
            any_nonzero_fixed = any_nonzero_fixed ||
                std::abs(expected_fixed[component]) > 0.0;
        }

        std::size_t expected_count = 0;
        double expected_nearest = std::numeric_limits<double>::infinity();
        for (const BondNeighbourhood& neighbour :
             conf.AtomAt(atom_index).bond_neighbours) {
            if (neighbour.bond_category != BondCategory::SidechainCO)
                continue;
            ++expected_count;
            expected_nearest = std::min(
                expected_nearest, neighbour.distance_to_midpoint);
        }
        EXPECT_DOUBLE_EQ(audit_rows.values[atom_index * 4 + 0],
                         tensor.T2Magnitude());
        EXPECT_TRUE(std::isnan(audit_rows.values[atom_index * 4 + 1]));
        EXPECT_DOUBLE_EQ(audit_rows.values[atom_index * 4 + 2],
                         static_cast<double>(expected_count));
        if (expected_count > 0) {
            EXPECT_DOUBLE_EQ(audit_rows.values[atom_index * 4 + 3],
                             expected_nearest);
        } else {
            EXPECT_TRUE(std::isnan(
                audit_rows.values[atom_index * 4 + 3]));
        }
    }
    EXPECT_TRUE(any_nonzero_fixed);

    RemoveSidechainCOOutputs(output_dir);
}


TEST(SidechainCarbonylAnisotropyProduction,
     MopacPresenceMakesBondOrderRowsFinite) {
    auto protein = BuildSyntheticSidechainCOProtein();
    auto& conf = protein->Conformation();
    ASSERT_TRUE(conf.AttachResult(GeometryResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(SpatialIndexResult::Compute(conf)));

    // A MopacResult is attached through the explicit test hook.  Its empty
    // topology-order table represents finite zero Wiberg weights; importantly,
    // the production presence branch must emit zeros rather than absence NaNs.
    conf.ForceAttachResultForTesting(std::make_unique<MopacResult>());
    ASSERT_TRUE(conf.AttachResult(McConnellResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(
        SidechainCarbonylAnisotropyResult::Compute(conf)));

    const fs::path output_dir =
        SidechainCOTempDir("sidechain_co_with_mopac");
    ASSERT_EQ(conf.Result<SidechainCarbonylAnisotropyResult>()
                  .WriteFeatures(conf, output_dir.string()),
              6);
    const auto bo_rows = ReadNpy<double>(
        output_dir / "sidechain_co_bo_T2.npy", "<f8");
    const auto audit_rows = ReadNpy<double>(
        output_dir / "sidechain_co_scalar_audit.npy", "<f8");
    ASSERT_EQ(bo_rows.shape,
              (std::vector<std::size_t>{conf.AtomCount(), 9}));
    ASSERT_EQ(audit_rows.shape,
              (std::vector<std::size_t>{conf.AtomCount(), 4}));
    for (double value : bo_rows.values) {
        EXPECT_TRUE(std::isfinite(value));
        EXPECT_DOUBLE_EQ(value, 0.0);
    }
    for (std::size_t atom_index = 0;
         atom_index < conf.AtomCount(); ++atom_index) {
        EXPECT_TRUE(std::isfinite(audit_rows.values[atom_index * 4 + 1]));
        EXPECT_DOUBLE_EQ(audit_rows.values[atom_index * 4 + 1], 0.0);
    }

    // Pin the upstream production weighting used by the copied BO channel
    // with a non-unit value; no test-side McConnell formula is re-derived.
    const Mat3 axial = (Mat3() <<
        1.0, 0.2, 0.0,
        0.0, -0.5, 0.1,
        0.0, 0.0, -0.5).finished();
    const auto weighted = mcconnell_result_detail::SelectChannelResponses(
        McConnellSourceCategory::SidechainCO, axial, Mat3::Zero(), 1.75);
    EXPECT_LT(MaxAbs(weighted.bond_order - 1.75 * axial), 1e-14);

    RemoveSidechainCOOutputs(output_dir);
}


TEST(McConnellImplementationChecks, PeptideCORhombicFallsBackToAxial) {
    {
        auto protein = BuildSyntheticPeptideCOProtein();
        auto& conf = protein->Conformation();
        conf.AttachResult(GeometryResult::Compute(conf));
        const size_t nc_bond = FirstBondBetween(*protein, 0, 1);
        ASSERT_NE(nc_bond, SIZE_MAX);
        ASSERT_NE(protein->BondAt(nc_bond).category, BondCategory::PeptideCO);

        const auto axial = McConnellResult::ComputePairKernel(
            conf.PositionAt(3),
            conf.bond_midpoints[nc_bond],
            conf.bond_directions[nc_bond]);
        const auto rhombic = McConnellResult::ComputePeptideCORhombicPairKernel(
            conf, nc_bond, conf.PositionAt(3));

        EXPECT_LT(MaxAbs(rhombic.source_shape - axial.source_shape), 1e-12);
        EXPECT_LT(MaxAbs(rhombic.response - axial.response), 1e-12);
    }

    {
        auto protein = BuildSyntheticPeptideCOProtein(true);
        auto& conf = protein->Conformation();
        conf.AttachResult(GeometryResult::Compute(conf));
        const size_t co_bond = FirstBondWithCategory(*protein, BondCategory::PeptideCO);
        ASSERT_NE(co_bond, SIZE_MAX);

        const auto axial = McConnellResult::ComputePairKernel(
            conf.PositionAt(3),
            conf.bond_midpoints[co_bond],
            conf.bond_directions[co_bond]);
        const auto rhombic = McConnellResult::ComputePeptideCORhombicPairKernel(
            conf, co_bond, conf.PositionAt(3));

        EXPECT_LT(MaxAbs(rhombic.source_shape - axial.source_shape), 1e-12);
        EXPECT_LT(MaxAbs(rhombic.response - axial.response), 1e-12);

        conf.AttachResult(SpatialIndexResult::Compute(conf));
        conf.AttachResult(McConnellResult::Compute(conf));
        EXPECT_LT(MaxAbs(conf.AtomAt(3)
                             .mcconnell_peptide_co_rhombic.Reconstruct()),
                  1e-12);
    }
}


TEST_F(McConnellProteinTest, AromaticCategoryIsExactlyZero) {
    const Protein& p = protein->Conformation().ProteinRef();
    int aromatic_bonds = 0;
    for (size_t bi = 0; bi < p.BondCount(); ++bi)
        if (p.BondAt(bi).category == BondCategory::Aromatic)
            ++aromatic_bonds;
    ASSERT_GT(aromatic_bonds, 0);

    auto& conf = protein->Conformation();
    conf.AttachResult(McConnellResult::Compute(conf));
    const size_t cat = static_cast<size_t>(
        McConnellSourceCategory::AromaticZeroed);
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        EXPECT_TRUE(SphericalAllZero(
            conf.AtomAt(ai).mcconnell_source_tensors[cat][0]));
        EXPECT_TRUE(SphericalAllZero(
            conf.AtomAt(ai).mcconnell_source_tensors[cat][1]));
    }
}


TEST_F(McConnellProteinTest, MissingMopacZerosOnlyBoChannel) {
    auto& conf = protein->Conformation();
    ASSERT_FALSE(conf.HasResult<MopacResult>());
    conf.AttachResult(McConnellResult::Compute(conf));

    double max_fixed_t2 = 0.0;
    double max_bo_abs = 0.0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& ca = conf.AtomAt(ai);
        for (size_t c = 0; c < kMcConnellSourceCategoryCount; ++c) {
            max_fixed_t2 = std::max(
                max_fixed_t2,
                ca.mcconnell_source_tensors[c][0].T2Magnitude());
            const auto& bo = ca.mcconnell_source_tensors[c][1];
            max_bo_abs = std::max(max_bo_abs, std::abs(bo.T0));
            for (double v : bo.T1) max_bo_abs = std::max(max_bo_abs, std::abs(v));
            for (double v : bo.T2) max_bo_abs = std::max(max_bo_abs, std::abs(v));
        }
    }
    EXPECT_GT(max_fixed_t2, 0.0);
    EXPECT_EQ(max_bo_abs, 0.0);
}


TEST_F(McConnellProteinTest, NearFieldAuditCountsAreReported) {
    auto& conf = protein->Conformation();
    conf.AttachResult(McConnellResult::Compute(conf));

    int accepted = 0;
    int rejected = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        accepted += conf.AtomAt(ai).mcconnell_near_field_accepted_lt3A;
        rejected += conf.AtomAt(ai).mcconnell_near_field_rejected_lt3A;
    }
    EXPECT_GE(accepted, 0);
    EXPECT_GE(rejected, 0);
    std::cout << "  McConnell near-field audit (<3 A): accepted="
              << accepted << " rejected=" << rejected << "\n";
}


TEST_F(McConnellProteinTest, WriteFeaturesEmitsTwentySixArraysAndManifest) {
    auto& conf = protein->Conformation();
    conf.AttachResult(McConnellResult::Compute(conf));

    const fs::path out_dir = fs::temp_directory_path() /
        ("mcconnell_features_" + std::to_string(::getpid()));
    fs::create_directories(out_dir);
    const auto& mc = conf.Result<McConnellResult>();
    EXPECT_EQ(mc.WriteFeatures(conf, out_dir.string()), 26);

    for (size_t c = 0; c < kMcConnellSourceCategoryCount; ++c) {
        const auto cat = static_cast<McConnellSourceCategory>(c);
        for (size_t h = 0; h < kMcConnellChannelCount; ++h) {
            const auto channel = static_cast<McConnellChannel>(h);
            const fs::path p = out_dir / (
                std::string("mc_") + McConnellSourceCategoryStem(cat) +
                "_" + McConnellChannelStem(channel) + ".npy");
            EXPECT_TRUE(fs::exists(p)) << p;
        }
    }
    EXPECT_TRUE(fs::exists(out_dir / "mc_peptide_co_rhombic.npy"));
    EXPECT_TRUE(fs::exists(out_dir / "mc_nearest_co_dir.npy"));
    EXPECT_TRUE(fs::exists(out_dir / "mc_nearest_co_midpoint.npy"));
    EXPECT_TRUE(fs::exists(out_dir / "mc_nearest_co_T2.npy"));
    EXPECT_TRUE(fs::exists(out_dir / "mc_nearest_cn_T2.npy"));

    const fs::path manifest = out_dir / "extraction_manifest.json";
    ASSERT_TRUE(fs::exists(manifest));
    const std::string text = ReadText(manifest);
    EXPECT_NE(text.find("\"source_model\": \"unit susceptibility shape; axial scale learned; peptide C=O fixed/BO responses use the full pinned rhombic shape\""),
              std::string::npos);
    EXPECT_NE(text.find("\"bo_source\": \"MOPAC Wiberg bond order\""),
              std::string::npos);
    EXPECT_NE(text.find("\"aromatic_zeroed_when_ring_active\": true"),
              std::string::npos);
    EXPECT_NE(text.find("\"irrep_layout\": \"0e,1e_x,1e_y,1e_z,2e_m-2..+2\""),
              std::string::npos);
    EXPECT_NE(text.find("\"units\": \"Angstrom^-3\""),
              std::string::npos);
    EXPECT_NE(text.find("\"rhombic_status\": \"peptide_co_pinned_present\""),
              std::string::npos);
    EXPECT_NE(text.find("\"rhombic_array\": \"mc_peptide_co_rhombic.npy\""),
              std::string::npos);
    EXPECT_NE(text.find("\"source\": \"Hooper & Kaiser 1965 Table III, EF-corrected acetamide A, Abraham-anchored sign\""),
              std::string::npos);
    EXPECT_TRUE(fs::exists(out_dir / "mc_nearfield_counts.npy"));

    for (size_t c = 0; c < kMcConnellSourceCategoryCount; ++c) {
        const auto cat = static_cast<McConnellSourceCategory>(c);
        for (size_t h = 0; h < kMcConnellChannelCount; ++h) {
            const auto channel = static_cast<McConnellChannel>(h);
            const fs::path p = out_dir / (
                std::string("mc_") + McConnellSourceCategoryStem(cat) +
                "_" + McConnellChannelStem(channel) + ".npy");
            std::remove(p.string().c_str());
        }
    }
    std::remove((out_dir / "mc_peptide_co_rhombic.npy").string().c_str());
    std::remove((out_dir / "mc_nearfield_counts.npy").string().c_str());
    std::remove(manifest.string().c_str());
    ::rmdir(out_dir.string().c_str());
}


TEST_F(McConnellProteinTest, NearestCOTrackingWorks) {
    auto& conf = protein->Conformation();
    conf.AttachResult(McConnellResult::Compute(conf));

    int has_co = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        double d = conf.AtomAt(ai).nearest_CO_dist;
        if (d < NO_DATA_SENTINEL) {
            has_co++;
            EXPECT_GT(d, 0.0);
            EXPECT_LT(d, MCCONNELL_CUTOFF_A);
        }
    }

    EXPECT_GT(has_co, static_cast<int>(conf.AtomCount() / 2))
        << "Most atoms should have a nearest CO bond within range";

    std::cout << "  Atoms with nearest CO: " << has_co
              << " / " << conf.AtomCount() << "\n";
}


// ============================================================================
// ORCA protein test (protonated, more atoms)
// ============================================================================

TEST(McConnellOrcaTest, RunOnProtonatedProtein) {
    OrcaRunFiles files;
    files.pdb_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.pdb";
    files.xyz_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.xyz";
    files.prmtop_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.prmtop";

    if (!fs::exists(files.xyz_path) || !fs::exists(files.prmtop_path))
        GTEST_SKIP() << "ORCA test data not found";

    auto load = BuildFromOrca(files);
    ASSERT_TRUE(load.Ok()) << load.error;

    auto& conf = load.protein->Conformation();
    conf.AttachResult(GeometryResult::Compute(conf));
    conf.AttachResult(SpatialIndexResult::Compute(conf));
    auto mc = McConnellResult::Compute(conf);
    ASSERT_NE(mc, nullptr);
    conf.AttachResult(std::move(mc));

    // Summary: bond neighbour counts, T0/T2 ranges
    double min_t0 = 1e30, max_t0 = -1e30;
    double max_t2 = 0;
    int total_bn = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& ca = conf.AtomAt(ai);
        total_bn += static_cast<int>(ca.bond_neighbours.size());
        double t0 = ca.mc_shielding_contribution.T0;
        min_t0 = std::min(min_t0, t0);
        max_t0 = std::max(max_t0, t0);
        max_t2 = std::max(max_t2, ca.mc_shielding_contribution.T2Magnitude());
    }

    std::cout << "  ORCA protein McConnell summary:\n"
              << "    atoms=" << conf.AtomCount()
              << " bonds=" << load.protein->BondCount()
              << " total_bond_neighbours=" << total_bn << "\n"
              << "    T0 range: [" << min_t0 << ", " << max_t0 << "]\n"
              << "    max |T2| = " << max_t2 << "\n";

    EXPECT_GT(total_bn, 0);
    EXPECT_GT(max_t0 - min_t0, 0.001) << "Should have a range of T0 values";
    EXPECT_GT(max_t2, 0.001) << "T2 should be non-zero";
}

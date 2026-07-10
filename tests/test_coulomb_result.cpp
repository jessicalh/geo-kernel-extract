#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>

#include "CoulombResult.h"
#include "MopacCoulombResult.h"
#include "ChargeAssignmentResult.h"
#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "ApbsFieldResult.h"
#include "Protein.h"
#include "PdbFileReader.h"
#include "OrcaRunLoader.h"
#include "ChargeSource.h"
#include "PhysicalConstants.h"

#include <filesystem>
namespace fs = std::filesystem;
using namespace nmr;

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
    EXPECT_EQ(std::string(magic, 6), std::string("\x93NUMPY", 6)) << path;

    char version[2] = {};
    in.read(version, 2);
    EXPECT_EQ(version[0], 1) << path;
    EXPECT_EQ(version[1], 0) << path;

    uint16_t header_len = 0;
    in.read(reinterpret_cast<char*>(&header_len), sizeof(header_len));
    std::string header(header_len, '\0');
    in.read(header.data(), header_len);

    const std::string descr_key = "'descr': '";
    auto descr_pos = header.find(descr_key);
    if (descr_pos == std::string::npos) {
        ADD_FAILURE() << header;
        return arr;
    }
    descr_pos += descr_key.size();
    auto descr_end = header.find('\'', descr_pos);
    if (descr_end == std::string::npos) {
        ADD_FAILURE() << header;
        return arr;
    }
    arr.descr = header.substr(descr_pos, descr_end - descr_pos);

    auto shape_key = header.find("'shape': (");
    if (shape_key == std::string::npos) {
        ADD_FAILURE() << header;
        return arr;
    }
    auto shape_begin = header.find('(', shape_key);
    auto shape_end = header.find(')', shape_begin);
    if (shape_begin == std::string::npos || shape_end == std::string::npos) {
        ADD_FAILURE() << header;
        return arr;
    }
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

const double* Doubles(const NpyArray& arr) {
    return reinterpret_cast<const double*>(arr.bytes.data());
}

void RemoveCoulombFeatureDir(const fs::path& out_dir) {
    for (const char* name : {
            "coulomb_efg.npy",
            "coulomb_efg_t2.npy",
            "coulomb_E.npy",
            "coulomb_E_backbone.npy",
            "coulomb_E_sidechain.npy",
            "coulomb_E_aromatic.npy",
            "coulomb_efg_backbone.npy",
            "coulomb_efg_sidechain.npy",
            "coulomb_efg_aromatic.npy",
            "coulomb_scalars.npy",
            "coulomb_aromatic_E_proj.npy",
            "coulomb_aromatic_n_src.npy",
            "coulomb_E_solvent.npy",
            "coulomb_efg_solvent.npy",
        }) {
        std::error_code ec;
        fs::remove(out_dir / name, ec);
    }
    std::error_code ec;
    fs::remove(out_dir, ec);
}

}  // namespace



// ============================================================================
// Analytical test: two point charges, hand-calculable E-field and EFG
// ============================================================================

TEST(CoulombAnalytical, TwoChargesKnownGeometry) {
    // Atom 0: charge +0.5 e at (0, 0, 0)
    // Atom 1: charge -0.3 e at (0, 0, 0)  -- dummy, both at same position
    // Actually: set up atoms with known positions and charges.
    //
    // Field point (atom 2) at (3, 0, 0), charge 0 (observer).
    // Source atom 0 at origin, charge +0.5 e.
    // Source atom 1 at (0, 3, 0), charge -0.3 e.
    //
    // E at atom 2 from atom 0:
    //   r = (3,0,0) - (0,0,0) = (3,0,0), |r| = 3
    //   E_0 = ke * 0.5 * (3,0,0) / 27 = ke * (0.5/27) * (3,0,0)
    //       = ke * (1/54, 0, 0) * 3 = ke * (1/18, 0, 0)
    //   Wait: E = ke * q * r / r^3 = ke * 0.5 * (3,0,0) / 27
    //       = ke * (1.5/27, 0, 0) = ke * (1/18, 0, 0)
    //
    // E at atom 2 from atom 1:
    //   r = (3,0,0) - (0,3,0) = (3,-3,0), |r| = 3*sqrt(2)
    //   r^3 = 27 * 2*sqrt(2) = 54*sqrt(2)
    //   E_1 = ke * (-0.3) * (3,-3,0) / (54*sqrt(2))
    //       = ke * (-0.3/54/sqrt(2)) * (3, -3, 0)
    //
    // Let's just verify the properties rather than exact values:
    // 1. EFG is traceless
    // 2. E is non-zero
    // 3. The total has the right order of magnitude

    auto protein = std::make_unique<Protein>();

    auto a0 = Atom::Create(Element::C);
    auto a1 = Atom::Create(Element::N);
    auto a2 = Atom::Create(Element::H);
    auto a3 = Atom::Create(Element::O);
    auto a4 = Atom::Create(Element::O);
    a0->residue_index = 0;
    a1->residue_index = 0;
    a2->residue_index = 0;
    a3->residue_index = 0;
    a4->residue_index = 0;
    a2->parent_atom_index = 0;

    protein->AddAtom(std::move(a0));
    protein->AddAtom(std::move(a1));
    protein->AddAtom(std::move(a2));
    protein->AddAtom(std::move(a3));
    protein->AddAtom(std::move(a4));

    Residue res;
    res.type = AminoAcid::ALA;
    res.sequence_number = 1;
    res.chain_id = "A";
    res.atom_indices = {0, 1, 2, 3, 4};
    protein->AddResidue(res);

    std::vector<Vec3> positions = {
        Vec3(0.0, 0.0, 0.0),   // atom 0: source, q = +0.5
        Vec3(0.0, 3.0, 0.0),   // atom 1: source, q = -0.3
        Vec3(3.0, 0.0, 0.0),   // atom 2: observer, q = 0
        Vec3(-1.2, 0.0, 0.0),  // atom 3: bonded zero-charge neighbour
        Vec3(0.0, 1.2, 0.0)    // atom 4: bonded zero-charge neighbour
    };

    protein->FinalizeConstruction(positions);
    // Synthetic topology has no semantic parent inference; set the explicit
    // H parent after finalization, which is the authority that may rebuild
    // parent indices.
    protein->MutableAtomAt(2).parent_atom_index = 0;
    // C11 forcing: atom 0 is a heavy target with multiple bonds. Reverse its
    // topology order so any "first bond" projection would change, while the
    // correct H-only contract remains NaN and order-independent.
    ASSERT_GE(protein->AtomAt(0).bond_indices.size(), 2u);
    std::reverse(protein->MutableAtomAt(0).bond_indices.begin(),
                 protein->MutableAtomAt(0).bond_indices.end());
    auto& conf = protein->AddCrystalConformation(positions, 0, 0, 0, "test");

    // Attach dependencies
    conf.AttachResult(GeometryResult::Compute(conf));

    // Typed dependency marker; source charges below are this analytical
    // fixture's independent authority, not a force-field table lookup.
    conf.ForceAttachResultForTesting(
        std::make_unique<ChargeAssignmentResult>());

    conf.MutableAtomAt(0).partial_charge = 0.5;
    conf.MutableAtomAt(1).partial_charge = -0.3;
    conf.MutableAtomAt(2).partial_charge = 0.0;
    conf.MutableAtomAt(3).partial_charge = 0.0;
    conf.MutableAtomAt(4).partial_charge = 0.0;
    conf.MutableAtomAt(0).mopac_charge = 0.5;
    conf.MutableAtomAt(1).mopac_charge = -0.3;
    conf.MutableAtomAt(2).mopac_charge = 0.0;
    conf.MutableAtomAt(3).mopac_charge = 0.0;
    conf.MutableAtomAt(4).mopac_charge = 0.0;

    conf.AttachResult(SpatialIndexResult::Compute(conf));

    auto coulomb = CoulombResult::Compute(conf);
    ASSERT_NE(coulomb, nullptr);
    conf.AttachResult(std::move(coulomb));

    // Check observer atom (index 2)
    const auto& ca = conf.AtomAt(2);

    // E-field should be non-zero
    EXPECT_GT(ca.coulomb_E_magnitude, 0.0)
        << "Observer should see non-zero E-field";

    // EFG should be traceless
    double trace = ca.coulomb_EFG_total.trace();
    EXPECT_NEAR(trace, 0.0, 1e-10)
        << "EFG must be traceless (Gauss's law)";

    // Verify analytically: E at atom 2 from atom 0 only
    // r = (3,0,0), |r| = 3, r^3 = 27
    // E = ke * 0.5 * (3,0,0) / 27 = ke * (1/18, 0, 0)
    double E_x_from_0 = COULOMB_KE * 0.5 * 3.0 / 27.0;  // ke/18
    double E_y_from_0 = 0.0;

    // E from atom 1: r = (3,-3,0), |r| = 3*sqrt(2), r^3 = 54*sqrt(2)
    double r1 = 3.0 * std::sqrt(2.0);
    double r1_3 = r1 * r1 * r1;
    double E_x_from_1 = COULOMB_KE * (-0.3) * 3.0 / r1_3;
    double E_y_from_1 = COULOMB_KE * (-0.3) * (-3.0) / r1_3;

    double E_x_expected = E_x_from_0 + E_x_from_1;
    double E_y_expected = E_y_from_0 + E_y_from_1;

    EXPECT_NEAR(ca.coulomb_E_total(0), E_x_expected, 1e-8)
        << "E_x should match analytical value";
    EXPECT_NEAR(ca.coulomb_E_total(1), E_y_expected, 1e-8)
        << "E_y should match analytical value";
    EXPECT_NEAR(ca.coulomb_E_total(2), 0.0, 1e-10)
        << "E_z should be zero (all atoms in xy plane)";

    // C11 source-of-truth check: projections have chemical meaning only for
    // H with a parent. Heavy atoms are NaN irrespective of their topology;
    // the H projection is parent-to-H dot E (the +x direction here).
    EXPECT_TRUE(std::isnan(conf.AtomAt(0).coulomb_E_bond_proj));
    EXPECT_TRUE(std::isnan(conf.AtomAt(1).coulomb_E_bond_proj));
    EXPECT_TRUE(std::isnan(conf.AtomAt(0).aromatic_E_bond_proj));
    EXPECT_TRUE(std::isnan(conf.AtomAt(1).aromatic_E_bond_proj));
    EXPECT_NEAR(ca.coulomb_E_bond_proj, ca.coulomb_E_total.x(), 1e-12);
    EXPECT_NEAR(ca.aromatic_E_bond_proj, ca.coulomb_E_aromatic.x(), 1e-12);

    // The MOPAC calculator has the same C11 H-only projection contract.
    // Its dependency is deliberately bypassed here because the charge values
    // above are the independent test fixture source of truth.
    auto mopac_coulomb = MopacCoulombResult::Compute(conf);
    ASSERT_NE(mopac_coulomb, nullptr);
    EXPECT_TRUE(std::isnan(conf.AtomAt(0).mopac_coulomb_E_bond_proj));
    EXPECT_TRUE(std::isnan(conf.AtomAt(1).mopac_coulomb_E_bond_proj));
    EXPECT_NEAR(conf.AtomAt(2).mopac_coulomb_E_bond_proj,
                conf.AtomAt(2).mopac_coulomb_E_total.x(), 1e-12);
    EXPECT_DOUBLE_EQ(
        conf.AtomAt(2).mopac_coulomb_aromatic_E_magnitude,
        conf.AtomAt(2).mopac_coulomb_E_aromatic.norm());

    // Verify EFG diagonal from atom 0 alone (d along x):
    // V_xx = ke * q * (3*1 - 1)/r^3 = ke * 0.5 * 2/27
    // V_yy = ke * q * (3*0 - 1)/r^3 = ke * 0.5 * (-1)/27
    // V_zz = ke * q * (3*0 - 1)/r^3 = ke * 0.5 * (-1)/27
    // Sum(V_diag) = ke * 0.5 * (2 - 1 - 1)/27 = 0  ✓

    std::cout << "  Analytical Coulomb test:\n"
              << "    E = (" << ca.coulomb_E_total(0) << ", "
              << ca.coulomb_E_total(1) << ", "
              << ca.coulomb_E_total(2) << ") V/A\n"
              << "    |E| = " << ca.coulomb_E_magnitude << " V/A\n"
              << "    EFG trace = " << trace << "\n"
              << "    EFG T2 magnitude = "
              << ca.coulomb_EFG_total_spherical.T2Magnitude() << "\n";
}


// ============================================================================
// Full protein test
// ============================================================================

class CoulombProteinTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!fs::exists(nmr::test::TestEnvironment::UbqProtonated()) || !fs::exists(nmr::test::TestEnvironment::Ff14sbParams())) {
            GTEST_SKIP() << "1UBQ.pdb or ff14sb_params.dat not found";
        }
        auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
        if (!r.Ok()) GTEST_SKIP() << "Failed to load 1UBQ";
        protein = std::move(r.protein);

        auto& conf = protein->Conformation();
        conf.AttachResult(GeometryResult::Compute(conf));
        conf.AttachResult(ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams()));
        conf.AttachResult(SpatialIndexResult::Compute(conf));
    }

    std::unique_ptr<Protein> protein;
};


TEST_F(CoulombProteinTest, ComputeAndAttach) {
    auto& conf = protein->Conformation();
    auto coulomb = CoulombResult::Compute(conf);
    ASSERT_NE(coulomb, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(coulomb)));
    ASSERT_TRUE(conf.HasResult<CoulombResult>());
}


TEST_F(CoulombProteinTest, EFieldNonZeroForAllAtoms) {
    auto& conf = protein->Conformation();
    conf.AttachResult(CoulombResult::Compute(conf));

    int nonzero = 0;
    double sum_mag = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        double mag = conf.AtomAt(ai).coulomb_E_magnitude;
        if (mag > 1e-10) nonzero++;
        sum_mag += mag;
    }

    double mean = sum_mag / static_cast<double>(conf.AtomCount());
    EXPECT_EQ(nonzero, static_cast<int>(conf.AtomCount()))
        << "Every atom should have non-zero E-field";

    // Mean E should be on the order of 0.1-10 V/A for a protein
    EXPECT_GT(mean, 0.01) << "Mean |E| too small";
    EXPECT_LT(mean, 100.0) << "Mean |E| too large";

    std::cout << "  Non-zero E: " << nonzero << " / " << conf.AtomCount()
              << ", mean |E| = " << mean << " V/A\n";
}


TEST_F(CoulombProteinTest, EFGIsTraceless) {
    auto& conf = protein->Conformation();
    conf.AttachResult(CoulombResult::Compute(conf));

    double max_trace = 0.0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        double trace = std::abs(conf.AtomAt(ai).coulomb_EFG_total.trace());
        max_trace = std::max(max_trace, trace);
    }

    EXPECT_LT(max_trace, 1e-8)
        << "All EFG tensors must be traceless (Gauss's law)";

    std::cout << "  Max |EFG trace| = " << max_trace << "\n";
}


TEST_F(CoulombProteinTest, DecompositionSumsToTotal) {
    auto& conf = protein->Conformation();
    conf.AttachResult(CoulombResult::Compute(conf));

    double max_E_diff = 0.0;
    double max_EFG_diff = 0.0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& ca = conf.AtomAt(ai);

        // E_total = E_backbone + E_sidechain + E_aromatic
        Vec3 E_sum = ca.coulomb_E_backbone + ca.coulomb_E_sidechain
                   + ca.coulomb_E_aromatic;
        double E_diff = (ca.coulomb_E_total - E_sum).norm();
        max_E_diff = std::max(max_E_diff, E_diff);

        // EFG_total = EFG_backbone + EFG_sidechain + EFG_aromatic
        Mat3 EFG_sum = ca.coulomb_EFG_backbone + ca.coulomb_EFG_sidechain
                     + ca.coulomb_EFG_aromatic;
        double EFG_diff = (ca.coulomb_EFG_total - EFG_sum).norm();
        max_EFG_diff = std::max(max_EFG_diff, EFG_diff);
    }

    EXPECT_LT(max_E_diff, 1e-8)
        << "E decomposition must sum to total";
    EXPECT_LT(max_EFG_diff, 1e-8)
        << "EFG decomposition must sum to total";

    std::cout << "  Max |E_total - (bb + sc + arom)| = " << max_E_diff << "\n";
    std::cout << "  Max |EFG_total - (bb + sc + arom)| = " << max_EFG_diff << "\n";
}


TEST_F(CoulombProteinTest, BackboneFractionIsReasonable) {
    auto& conf = protein->Conformation();
    conf.AttachResult(CoulombResult::Compute(conf));

    double sum_proj = 0;
    int count = 0;
    int positive = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        double proj = conf.AtomAt(ai).coulomb_E_backbone_frac;
        // Projection of E_backbone along E_total direction.
        // Positive = backbone aligned with total. Negative = opposed.
        // Bounded by |E_backbone|.
        sum_proj += proj;
        if (proj > 0) positive++;
        count++;
    }

    double mean_proj = sum_proj / count;
    // Most atoms should have positive projection (backbone dominates)
    EXPECT_GT(positive, count / 2)
        << "Majority of atoms should have backbone aligned with total E";
    EXPECT_GT(mean_proj, 0.0)
        << "Mean backbone projection should be positive";

    std::cout << "  Mean backbone E projection = " << mean_proj << " V/A\n"
              << "  Positive projection: " << positive << " / " << count << "\n";
}


TEST_F(CoulombProteinTest, AromaticEIsStrongerNearAromaticRings) {
    auto& conf = protein->Conformation();
    conf.AttachResult(CoulombResult::Compute(conf));

    ASSERT_GT(conf.ring_geometries.size(), 0u)
        << "1UBQ should have aromatic rings";

    constexpr double kNearRingCenter = 6.0;
    constexpr double kFarRingCenter = 12.0;

    double sum_near = 0.0;
    double sum_far = 0.0;
    int near_count = 0;
    int far_count = 0;
    double max_mag_diff = 0.0;

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        double nearest_ring_center = std::numeric_limits<double>::infinity();
        for (const auto& geom : conf.ring_geometries) {
            double d = (conf.PositionAt(ai) - geom.center).norm();
            nearest_ring_center = std::min(nearest_ring_center, d);
        }

        const auto& ca = conf.AtomAt(ai);
        max_mag_diff = std::max(
            max_mag_diff,
            std::abs(ca.aromatic_E_magnitude - ca.coulomb_E_aromatic.norm()));

        if (nearest_ring_center <= kNearRingCenter) {
            sum_near += ca.aromatic_E_magnitude;
            near_count++;
        } else if (nearest_ring_center >= kFarRingCenter) {
            sum_far += ca.aromatic_E_magnitude;
            far_count++;
        }
    }

    ASSERT_GT(near_count, 0);
    ASSERT_GT(far_count, 0);

    double mean_near = sum_near / near_count;
    double mean_far = sum_far / far_count;

    EXPECT_LT(max_mag_diff, 1e-12)
        << "aromatic_E_magnitude should match |coulomb_E_aromatic|";
    EXPECT_GT(mean_near, mean_far)
        << "Aromatic E-field should be stronger near aromatic rings";

    std::cout << "  Aromatic E near rings = " << mean_near << " V/A over "
              << near_count << " atoms; far = " << mean_far << " V/A over "
              << far_count << " atoms\n";
}


TEST_F(CoulombProteinTest, T2IsNonZero) {
    auto& conf = protein->Conformation();
    conf.AttachResult(CoulombResult::Compute(conf));

    int nonzero_t2 = 0;
    double max_t2 = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        double t2 = conf.AtomAt(ai).coulomb_EFG_total_spherical.T2Magnitude();
        if (t2 > 1e-8) nonzero_t2++;
        max_t2 = std::max(max_t2, t2);
    }

    EXPECT_EQ(nonzero_t2, static_cast<int>(conf.AtomCount()))
        << "Every atom should have non-zero EFG T2";
    EXPECT_GT(max_t2, 0.01) << "Max T2 should be appreciable";

    std::cout << "  T2 nonzero: " << nonzero_t2
              << ", max |T2| = " << max_t2 << "\n";
}

TEST_F(CoulombProteinTest, WriteFeaturesKeepsFull9AndAddsT2Companion) {
    auto& conf = protein->Conformation();
    conf.AttachResult(CoulombResult::Compute(conf));
    const auto& result = conf.Result<CoulombResult>();

    const fs::path out_dir = fs::temp_directory_path() /
        ("coulomb_write_features_" + std::to_string(::getpid()));
    RemoveCoulombFeatureDir(out_dir);
    fs::create_directories(out_dir);

    const int written = result.WriteFeatures(conf, out_dir.string());
    EXPECT_EQ(written, 12);

    auto full = ReadNpy(out_dir / "coulomb_efg.npy");
    auto t2 = ReadNpy(out_dir / "coulomb_efg_t2.npy");
    const size_t N = conf.AtomCount();
    ASSERT_EQ(full.shape, (std::vector<size_t>{N, 9}));
    ASSERT_EQ(t2.shape, (std::vector<size_t>{N, 5}));
    EXPECT_EQ(full.descr, "<f8");
    EXPECT_EQ(t2.descr, "<f8");

    const double* full_data = Doubles(full);
    const double* t2_data = Doubles(t2);
    for (size_t i = 0; i < N; ++i) {
        for (size_t c = 0; c < 4; ++c) {
            EXPECT_NEAR(full_data[i*9 + c], 0.0, 1e-10)
                << "Coulomb EFG full9 T0/T1 structural zero at atom "
                << i << " col " << c;
        }
        for (size_t k = 0; k < 5; ++k) {
            EXPECT_DOUBLE_EQ(t2_data[i*5 + k], full_data[i*9 + 4 + k])
                << "coulomb_efg_t2 must mirror coulomb_efg columns 4:9";
        }
    }

    RemoveCoulombFeatureDir(out_dir);
}


// ============================================================================
// ORCA protein test (protonated, with prmtop charges)
// ============================================================================

TEST(CoulombOrcaTest, RunOnProtonatedProtein) {
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

    PrmtopChargeSource charge_source(files.prmtop_path);
    conf.AttachResult(ChargeAssignmentResult::Compute(conf, charge_source));
    conf.AttachResult(SpatialIndexResult::Compute(conf));

    auto coulomb = CoulombResult::Compute(conf);
    ASSERT_NE(coulomb, nullptr);
    conf.AttachResult(std::move(coulomb));

    // Summary
    double min_E = 1e30, max_E = 0;
    double max_t2 = 0;
    double max_trace = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& ca = conf.AtomAt(ai);
        double E_mag = ca.coulomb_E_magnitude;
        min_E = std::min(min_E, E_mag);
        max_E = std::max(max_E, E_mag);
        max_t2 = std::max(max_t2, ca.coulomb_EFG_total_spherical.T2Magnitude());
        max_trace = std::max(max_trace, std::abs(ca.coulomb_EFG_total.trace()));
    }

    std::cout << "  ORCA protein Coulomb summary:\n"
              << "    atoms=" << conf.AtomCount() << "\n"
              << "    |E| range: [" << min_E << ", " << max_E << "] V/A\n"
              << "    max |T2| = " << max_t2 << " V/A^2\n"
              << "    max |EFG trace| = " << max_trace << "\n";

    // Sanity checks
    EXPECT_GT(max_E, 0.1) << "Should have appreciable E-fields";
    EXPECT_LT(max_trace, 1e-8) << "EFG must be traceless";
    EXPECT_GT(max_t2, 0.01) << "T2 should be non-zero";
}


// ============================================================================
// APBS comparison: the Coulomb solvent surface is an exact read-back alias
// of the independently computed canonical APBS reaction field.
// ============================================================================

TEST(CoulombApbsComparison, SolventContributionIsReasonable) {
    if (!fs::exists(nmr::test::TestEnvironment::UbqProtonated()) || !fs::exists(nmr::test::TestEnvironment::Ff14sbParams()))
        GTEST_SKIP() << "Test data not found";

    auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
    if (!r.Ok()) GTEST_SKIP() << r.error;

    auto& conf = r.protein->Conformation();
    conf.AttachResult(GeometryResult::Compute(conf));
    conf.AttachResult(ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams()));
    conf.AttachResult(SpatialIndexResult::Compute(conf));

    // APBS first (needs charges)
    auto apbs = ApbsFieldResult::Compute(conf);
    if (!apbs) GTEST_SKIP() << "APBS failed";
    conf.AttachResult(std::move(apbs));

    // Then Coulomb (will copy the already-defined APBS reaction field)
    conf.AttachResult(CoulombResult::Compute(conf));

    // A11 frozen-surface forcing: cross the two independent writers and pin
    // the actual NPY descriptors, shapes, and payload bytes.  An in-memory
    // equality alone would not catch one alias being emitted from a stale
    // total-PB field.
    const fs::path out_dir = fs::temp_directory_path() /
        ("coulomb_apbs_alias_" + std::to_string(::getpid()));
    fs::create_directories(out_dir);
    EXPECT_EQ(conf.Result<ApbsFieldResult>().WriteFeatures(
                  conf, out_dir.string()), 8);
    EXPECT_EQ(conf.Result<CoulombResult>().WriteFeatures(
                  conf, out_dir.string()), 14);

    const auto apbs_E = ReadNpy(out_dir / "apbs_E.npy");
    const auto coulomb_E = ReadNpy(out_dir / "coulomb_E_solvent.npy");
    EXPECT_EQ(apbs_E.descr, "<f8");
    EXPECT_EQ(coulomb_E.descr, apbs_E.descr);
    EXPECT_EQ(apbs_E.shape,
              (std::vector<size_t>{conf.AtomCount(), 3u}));
    EXPECT_EQ(coulomb_E.shape, apbs_E.shape);
    EXPECT_EQ(coulomb_E.bytes, apbs_E.bytes);

    const auto apbs_efg = ReadNpy(out_dir / "apbs_efg.npy");
    const auto coulomb_efg =
        ReadNpy(out_dir / "coulomb_efg_solvent.npy");
    EXPECT_EQ(apbs_efg.descr, "<f8");
    EXPECT_EQ(coulomb_efg.descr, apbs_efg.descr);
    EXPECT_EQ(apbs_efg.shape,
              (std::vector<size_t>{conf.AtomCount(), 5u}));
    EXPECT_EQ(coulomb_efg.shape, apbs_efg.shape);
    EXPECT_EQ(coulomb_efg.bytes, apbs_efg.bytes);

    int has_solvent = 0;
    double mean_ratio = 0;
    int count = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& ca = conf.AtomAt(ai);
        for (int d = 0; d < 3; ++d) {
            EXPECT_DOUBLE_EQ(ca.coulomb_E_solvent(d), ca.apbs_efield(d));
        }
        for (int r_i = 0; r_i < 3; ++r_i) {
            for (int c_i = 0; c_i < 3; ++c_i) {
                EXPECT_DOUBLE_EQ(ca.coulomb_EFG_solvent(r_i, c_i),
                                 ca.apbs_efg(r_i, c_i));
            }
        }
        EXPECT_DOUBLE_EQ(ca.coulomb_EFG_solvent_spherical.T0,
                         ca.apbs_efg_spherical.T0);
        for (int component = 0; component < 3; ++component) {
            EXPECT_DOUBLE_EQ(
                ca.coulomb_EFG_solvent_spherical.T1[component],
                ca.apbs_efg_spherical.T1[component]);
        }
        for (int component = 0; component < 5; ++component) {
            EXPECT_DOUBLE_EQ(
                ca.coulomb_EFG_solvent_spherical.T2[component],
                ca.apbs_efg_spherical.T2[component]);
        }
        double solv_mag = ca.coulomb_E_solvent.norm();
        double vac_mag = ca.coulomb_E_magnitude;
        if (solv_mag > 1e-10) has_solvent++;
        if (vac_mag > 1e-6) {
            mean_ratio += solv_mag / vac_mag;
            count++;
        }
    }
    if (count > 0) mean_ratio /= count;

    std::cout << "  APBS vs vacuum Coulomb:\n"
              << "    Atoms with non-zero solvent E: " << has_solvent
              << " / " << conf.AtomCount() << "\n"
              << "    Mean |E_solvent|/|E_vacuum| = " << mean_ratio << "\n";

    // Solvation should modify the E-field (screening), so solvent
    // contribution should be non-zero for most atoms
    EXPECT_GT(has_solvent, static_cast<int>(conf.AtomCount() / 2))
        << "Most atoms should have non-zero solvent contribution";
    for (const char* name : {
            "apbs_E.npy", "apbs_efg.npy", "apbs_phi.npy",
            "apbs_E_clamp_mask.npy", "apbs_E_clamp_scale.npy",
            "apbs_nonfinite_sanitizer_mask.npy",
            "apbs_E_total_diagnostic.npy",
            "apbs_efg_total_diagnostic.npy"}) {
        std::error_code ec;
        fs::remove(out_dir / name, ec);
    }
    RemoveCoulombFeatureDir(out_dir);
}

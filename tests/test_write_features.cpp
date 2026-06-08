#include "TestEnvironment.h"
//
// test_write_features.cpp
//
// Runs the full pipeline on P84477 and writes all features to NPY files.
// This is the baseline for the topology refactor: after refactoring,
// the output directory should be byte-identical.
//

#include <gtest/gtest.h>
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>

#include "ChargeAssignmentResult.h"
#include "ConformationResult.h"
#include "CoulombResult.h"
#include "EnrichmentResult.h"
#include "HBondResult.h"
#include "LarsenHBondShieldingResult.h"
#include "McConnellResult.h"
#include "MopacResult.h"
#include "OperationRunner.h"
#include "OrcaRunLoader.h"
#include "PdbFileReader.h"
#include "ProteinConformation.h"
#include "SpatialIndexResult.h"
#include "ChargeSource.h"
#include "OperationLog.h"

namespace fs = std::filesystem;
using namespace nmr;

namespace {

struct NpyArray {
    std::string descr;
    std::vector<size_t> shape;
    std::vector<char> bytes;
};

std::string ReadTextFile(const fs::path& path) {
    std::ifstream in(path);
    return std::string(std::istreambuf_iterator<char>(in),
                       std::istreambuf_iterator<char>());
}

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
        ADD_FAILURE() << "NPY header missing descr: " << header;
        return arr;
    }
    descr_pos += descr_key.size();
    auto descr_end = header.find('\'', descr_pos);
    if (descr_end == std::string::npos) {
        ADD_FAILURE() << "NPY header malformed descr: " << header;
        return arr;
    }
    arr.descr = header.substr(descr_pos, descr_end - descr_pos);

    auto shape_key = header.find("'shape': (");
    if (shape_key == std::string::npos) {
        ADD_FAILURE() << "NPY header missing shape: " << header;
        return arr;
    }
    auto shape_begin = header.find('(', shape_key);
    auto shape_end = header.find(')', shape_begin);
    if (shape_begin == std::string::npos || shape_end == std::string::npos) {
        ADD_FAILURE() << "NPY header malformed shape: " << header;
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

std::string g_readback_context;

void SetReadbackContext(std::string context) {
    g_readback_context = std::move(context);
}

void ExpectDoubleEqOrNan(double actual, double expected) {
    if (::testing::Test::HasFailure()) return;
    if (std::isnan(expected)) {
        if (!std::isnan(actual)) {
            ADD_FAILURE() << g_readback_context
                          << " actual=" << actual << " expected=NaN";
        }
    } else if (actual != expected) {
        ADD_FAILURE() << g_readback_context
                      << " actual=" << actual << " expected=" << expected;
    }
}

void PackMat3RowMajor(const Mat3& m, double* out) {
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            out[r*3 + c] = m(r, c);
        }
    }
}

void ExpectCatalogAndSdkPresence(const fs::path& repo_root,
                                 const std::vector<std::string>& stems) {
    const std::string catalog = ReadTextFile(repo_root / "python/nmr_extract/_catalog.py");
    const std::string protein = ReadTextFile(repo_root / "python/nmr_extract/_protein.py");
    ASSERT_FALSE(catalog.empty());
    ASSERT_FALSE(protein.empty());
    for (const std::string& stem : stems) {
        EXPECT_NE(catalog.find("ArraySpec(\"" + stem + "\""), std::string::npos)
            << stem << " missing from _catalog.py";
        EXPECT_NE(protein.find("\"" + stem + "\""), std::string::npos)
            << stem << " missing from _protein.py loader exposure";
    }
}

void AssertCompleteEmitReadback(ProteinConformation& conf,
                                const RunResult& result,
                                const fs::path& out_dir) {
    const std::vector<std::string> new_stems = {
        "enrichment_role",
        "enrichment_hybridisation",
        "enrichment_flags",
        "ff_partial_charge",
        "ff_pb_radius",
        "coulomb_aromatic_E_proj",
        "coulomb_aromatic_n_src",
        "hbond_nearest_dir",
        "hbond_nearest_tensor",
        "hbond_nearest_spherical",
        "hbond_flags",
        "mc_nearest_co_dir",
        "mc_nearest_co_midpoint",
        "mc_nearest_co_T2",
        "mc_nearest_cn_T2",
        "larsen_corner_imputed",
        "mopac_bond_neighbors",
        "ring_direction_to_center",
        "disp_per_ring_tensor",
        "disp_per_ring_spherical",
        "piquad_quad_scalar",
    };

#ifdef NMR_TEST_DATA_DIR
    const fs::path repo_root = fs::path(NMR_TEST_DATA_DIR).parent_path().parent_path();
    ExpectCatalogAndSdkPresence(repo_root, new_stems);
#endif

    fs::create_directories(out_dir);
    for (const auto& stem : new_stems) {
        fs::remove(out_dir / (stem + ".npy"));
    }

    int arrays = ConformationResult::WriteAllFeatures(conf, out_dir.string());
    EXPECT_GT(arrays, 25) << "Expected 25+ arrays from calculators + identity";

    std::cout << "\n  WriteFeatures readback written to " << out_dir << "\n"
              << "  Arrays: " << arrays << "\n"
              << "  Atoms: " << conf.AtomCount() << "\n"
              << "  Results attached: " << result.attached.size() << "\n";

    for (const auto& entry : fs::directory_iterator(out_dir)) {
        if (!entry.is_regular_file()) continue;
        auto sz = fs::file_size(entry);
        std::cout << "    " << entry.path().filename().string()
                  << " (" << sz << " bytes)\n";
    }

    const size_t N = conf.AtomCount();

    if (conf.HasResult<EnrichmentResult>()) {
        auto role = ReadNpy(out_dir / "enrichment_role.npy");
        ExpectShapeAndDescr(role, {N}, "<i4");
        const int32_t* v = DataAs<int32_t>(role);
        for (size_t i = 0; i < N; ++i)
            EXPECT_EQ(v[i], static_cast<int32_t>(conf.AtomAt(i).role));

        auto hybrid = ReadNpy(out_dir / "enrichment_hybridisation.npy");
        ExpectShapeAndDescr(hybrid, {N}, "<i4");
        v = DataAs<int32_t>(hybrid);
        for (size_t i = 0; i < N; ++i)
            EXPECT_EQ(v[i], static_cast<int32_t>(conf.AtomAt(i).hybridisation));

        auto flags = ReadNpy(out_dir / "enrichment_flags.npy");
        ExpectShapeAndDescr(flags, {N, 8}, "|i1");
        const int8_t* f = DataAs<int8_t>(flags);
        for (size_t i = 0; i < N; ++i) {
            const auto& ca = conf.AtomAt(i);
            EXPECT_EQ(f[i*8 + 0], ca.is_backbone ? 1 : 0);
            EXPECT_EQ(f[i*8 + 1], ca.is_amide_H ? 1 : 0);
            EXPECT_EQ(f[i*8 + 2], ca.is_alpha_H ? 1 : 0);
            EXPECT_EQ(f[i*8 + 3], ca.is_methyl ? 1 : 0);
            EXPECT_EQ(f[i*8 + 4], ca.is_aromatic_H ? 1 : 0);
            EXPECT_EQ(f[i*8 + 5], ca.is_hbond_donor ? 1 : 0);
            EXPECT_EQ(f[i*8 + 6], ca.is_hbond_acceptor ? 1 : 0);
            EXPECT_EQ(f[i*8 + 7], ca.is_on_aromatic_residue ? 1 : 0);
        }
    }

    if (conf.HasResult<ChargeAssignmentResult>()) {
        SetReadbackContext("ff_partial_charge");
        auto charges = ReadNpy(out_dir / "ff_partial_charge.npy");
        ExpectShapeAndDescr(charges, {N}, "<f8");
        const double* q = DataAs<double>(charges);
        for (size_t i = 0; i < N; ++i)
            ExpectDoubleEqOrNan(q[i], conf.AtomAt(i).partial_charge);

        SetReadbackContext("ff_pb_radius");
        auto radii = ReadNpy(out_dir / "ff_pb_radius.npy");
        ExpectShapeAndDescr(radii, {N}, "<f8");
        const double* r = DataAs<double>(radii);
        for (size_t i = 0; i < N; ++i)
            ExpectDoubleEqOrNan(r[i], conf.AtomAt(i).pb_radius);
    }

    {
        size_t rows = 0;
        for (size_t i = 0; i < N; ++i) rows += conf.AtomAt(i).ring_neighbours.size();
        if (rows > 0) {
            SetReadbackContext("ring sparse arrays");
            auto dir = ReadNpy(out_dir / "ring_direction_to_center.npy");
            auto disp_tensor = ReadNpy(out_dir / "disp_per_ring_tensor.npy");
            auto disp_spherical = ReadNpy(out_dir / "disp_per_ring_spherical.npy");
            auto quad_scalar = ReadNpy(out_dir / "piquad_quad_scalar.npy");
            ExpectShapeAndDescr(dir, {rows, 3}, "<f8");
            ExpectShapeAndDescr(disp_tensor, {rows, 9}, "<f8");
            ExpectShapeAndDescr(disp_spherical, {rows, 9}, "<f8");
            ExpectShapeAndDescr(quad_scalar, {rows}, "<f8");

            const double* dv = DataAs<double>(dir);
            const double* mt = DataAs<double>(disp_tensor);
            const double* st = DataAs<double>(disp_spherical);
            const double* qs = DataAs<double>(quad_scalar);
            size_t row = 0;
            for (size_t i = 0; i < N; ++i) {
                for (const auto& rn : conf.AtomAt(i).ring_neighbours) {
                    ExpectDoubleEqOrNan(dv[row*3 + 0], rn.direction_to_center.x());
                    ExpectDoubleEqOrNan(dv[row*3 + 1], rn.direction_to_center.y());
                    ExpectDoubleEqOrNan(dv[row*3 + 2], rn.direction_to_center.z());
                    double packed_mat[9];
                    PackMat3RowMajor(rn.disp_tensor, packed_mat);
                    double packed_sph[9];
                    rn.disp_spherical.PackFull9(packed_sph);
                    for (int c = 0; c < 9; ++c) {
                        ExpectDoubleEqOrNan(mt[row*9 + c], packed_mat[c]);
                        ExpectDoubleEqOrNan(st[row*9 + c], packed_sph[c]);
                    }
                    ExpectDoubleEqOrNan(qs[row], rn.quad_scalar);
                    ++row;
                }
            }
            EXPECT_EQ(row, rows);
        }
    }

    if (conf.HasResult<CoulombResult>()) {
        SetReadbackContext("coulomb_aromatic_E_proj");
        auto proj = ReadNpy(out_dir / "coulomb_aromatic_E_proj.npy");
        ExpectShapeAndDescr(proj, {N}, "<f8");
        const double* p = DataAs<double>(proj);
        for (size_t i = 0; i < N; ++i)
            ExpectDoubleEqOrNan(p[i], conf.AtomAt(i).aromatic_E_bond_proj);

        auto count = ReadNpy(out_dir / "coulomb_aromatic_n_src.npy");
        ExpectShapeAndDescr(count, {N}, "<i4");
        const int32_t* c = DataAs<int32_t>(count);
        for (size_t i = 0; i < N; ++i)
            EXPECT_EQ(c[i], static_cast<int32_t>(conf.AtomAt(i).aromatic_n_sidechain_atoms));
    }

    if (conf.HasResult<HBondResult>()) {
        SetReadbackContext("hbond_nearest arrays");
        auto dir = ReadNpy(out_dir / "hbond_nearest_dir.npy");
        auto tensor = ReadNpy(out_dir / "hbond_nearest_tensor.npy");
        auto spherical = ReadNpy(out_dir / "hbond_nearest_spherical.npy");
        auto flags = ReadNpy(out_dir / "hbond_flags.npy");
        ExpectShapeAndDescr(dir, {N, 3}, "<f8");
        ExpectShapeAndDescr(tensor, {N, 9}, "<f8");
        ExpectShapeAndDescr(spherical, {N, 9}, "<f8");
        ExpectShapeAndDescr(flags, {N, 3}, "|i1");
        const double* dv = DataAs<double>(dir);
        const double* tv = DataAs<double>(tensor);
        const double* sv = DataAs<double>(spherical);
        const int8_t* fv = DataAs<int8_t>(flags);
        for (size_t i = 0; i < N; ++i) {
            const auto& ca = conf.AtomAt(i);
            ExpectDoubleEqOrNan(dv[i*3 + 0], ca.hbond_nearest_dir.x());
            ExpectDoubleEqOrNan(dv[i*3 + 1], ca.hbond_nearest_dir.y());
            ExpectDoubleEqOrNan(dv[i*3 + 2], ca.hbond_nearest_dir.z());
            double packed[9];
            PackMat3RowMajor(ca.hbond_nearest_tensor, packed);
            for (int c = 0; c < 9; ++c)
                ExpectDoubleEqOrNan(tv[i*9 + c], packed[c]);
            ca.hbond_nearest_spherical.PackFull9(packed);
            for (int c = 0; c < 9; ++c)
                ExpectDoubleEqOrNan(sv[i*9 + c], packed[c]);
            EXPECT_EQ(fv[i*3 + 0], ca.hbond_is_backbone ? 1 : 0);
            EXPECT_EQ(fv[i*3 + 1], ca.hbond_is_donor ? 1 : 0);
            EXPECT_EQ(fv[i*3 + 2], ca.hbond_is_acceptor ? 1 : 0);
        }
    }

    if (conf.HasResult<McConnellResult>()) {
        SetReadbackContext("mc_nearest arrays");
        auto co_dir = ReadNpy(out_dir / "mc_nearest_co_dir.npy");
        auto co_mid = ReadNpy(out_dir / "mc_nearest_co_midpoint.npy");
        auto co_t2 = ReadNpy(out_dir / "mc_nearest_co_T2.npy");
        auto cn_t2 = ReadNpy(out_dir / "mc_nearest_cn_T2.npy");
        ExpectShapeAndDescr(co_dir, {N, 3}, "<f8");
        ExpectShapeAndDescr(co_mid, {N, 3}, "<f8");
        ExpectShapeAndDescr(co_t2, {N, 9}, "<f8");
        ExpectShapeAndDescr(cn_t2, {N, 9}, "<f8");
        const double* cod = DataAs<double>(co_dir);
        const double* com = DataAs<double>(co_mid);
        const double* cot = DataAs<double>(co_t2);
        const double* cnt = DataAs<double>(cn_t2);
        for (size_t i = 0; i < N; ++i) {
            const auto& ca = conf.AtomAt(i);
            ExpectDoubleEqOrNan(cod[i*3 + 0], ca.dir_nearest_CO.x());
            ExpectDoubleEqOrNan(cod[i*3 + 1], ca.dir_nearest_CO.y());
            ExpectDoubleEqOrNan(cod[i*3 + 2], ca.dir_nearest_CO.z());
            ExpectDoubleEqOrNan(com[i*3 + 0], ca.nearest_CO_midpoint.x());
            ExpectDoubleEqOrNan(com[i*3 + 1], ca.nearest_CO_midpoint.y());
            ExpectDoubleEqOrNan(com[i*3 + 2], ca.nearest_CO_midpoint.z());
            double co_pack[9], cn_pack[9];
            ca.T2_CO_nearest.PackFull9(co_pack);
            ca.T2_CN_nearest.PackFull9(cn_pack);
            for (int c = 0; c < 9; ++c) {
                ExpectDoubleEqOrNan(cot[i*9 + c], co_pack[c]);
                ExpectDoubleEqOrNan(cnt[i*9 + c], cn_pack[c]);
            }
        }

    }

    if (conf.HasResult<MopacResult>()) {
        SetReadbackContext("mopac_bond_neighbors");
        size_t rows = 0;
        for (size_t i = 0; i < N; ++i) rows += conf.AtomAt(i).mopac_bond_neighbours.size();
        auto arr = ReadNpy(out_dir / "mopac_bond_neighbors.npy");
        ExpectShapeAndDescr(arr, {rows, 4}, "<f8");
        const double* d = DataAs<double>(arr);
        size_t row = 0;
        for (size_t i = 0; i < N; ++i) {
            for (const auto& nb : conf.AtomAt(i).mopac_bond_neighbours) {
                EXPECT_DOUBLE_EQ(d[row*4 + 0], static_cast<double>(i));
                EXPECT_DOUBLE_EQ(d[row*4 + 1], static_cast<double>(nb.other_atom));
                ExpectDoubleEqOrNan(d[row*4 + 2], nb.wiberg_order);
                const double topo = nb.topology_bond_index == SIZE_MAX
                    ? -1.0
                    : static_cast<double>(nb.topology_bond_index);
                EXPECT_DOUBLE_EQ(d[row*4 + 3], topo);
                ++row;
            }
        }
        EXPECT_EQ(row, rows);
    }

    if (conf.HasResult<LarsenHBondShieldingResult>()) {
        auto arr = ReadNpy(out_dir / "larsen_corner_imputed.npy");
        ExpectShapeAndDescr(arr, {N}, "|i1");
        const int8_t* d = DataAs<int8_t>(arr);
        for (size_t i = 0; i < N; ++i)
            EXPECT_EQ(d[i], conf.AtomAt(i).larsen_hbond_any_corner_imputed ? 1 : 0);
    }
}

TEST(WriteFeatures, UbqProtonatedReadback) {
    if (!fs::exists(nmr::test::TestEnvironment::UbqProtonated())) {
        GTEST_SKIP() << "1UBQ not found";
    }

    auto build = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(build.Ok()) << build.error;
    auto& conf = build.protein->Conformation();

    RunOptions opts;
    opts.charge_source = build.charges.get();
    opts.net_charge = build.net_charge;
    opts.skip_apbs = true;

    auto result = OperationRunner::Run(conf, opts);
    ASSERT_TRUE(result.Ok()) << result.error;
    ASSERT_TRUE(conf.HasResult<MopacResult>())
        << "MOPAC must be attached to read back mopac_bond_neighbors";
    ASSERT_TRUE(conf.HasResult<HBondResult>())
        << "HBond must be attached to read back hbond_nearest_*";

    const fs::path out_dir = fs::temp_directory_path() /
        ("write_features_ubq_readback_" + std::to_string(::getpid()));
    fs::remove_all(out_dir);
    AssertCompleteEmitReadback(conf, result, out_dir);
}

}  // namespace



TEST(WriteFeatures, A0A7C5FAR6OrcaReadback) {
    std::string dir = std::string(nmr::test::TestEnvironment::OrcaDir());

    // Load the AMBER-prepared ORCA pose (prmtop-derived charges) with the full pipeline.
    OrcaRunFiles files;
    files.pdb_path = dir + "A0A7C5FAR6_WT.pdb";
    files.xyz_path = dir + "A0A7C5FAR6_WT.xyz";
    files.prmtop_path = dir + "A0A7C5FAR6_WT.prmtop";
    if (!fs::exists(files.xyz_path) || !fs::exists(files.prmtop_path))
        GTEST_SKIP() << "ORCA test data not found";

    auto load = BuildFromOrca(files);
    ASSERT_TRUE(load.Ok()) << load.error;
    auto& conf = load.protein->Conformation();

    RunOptions opts;
    opts.charge_source = load.charges.get();
    opts.net_charge = load.net_charge;

    // Find NMR output
    for (const auto& entry : fs::directory_iterator(dir)) {
        std::string name = entry.path().filename().string();
        if (name.find("A0A7C5FAR6_WT") == 0 && name.find("_nmr.out") != std::string::npos) {
            opts.orca_nmr_path = entry.path().string();
            break;
        }
    }

    auto result = OperationRunner::Run(conf, opts);
    ASSERT_TRUE(result.Ok()) << result.error;

    const std::vector<std::string> new_stems = {
        "enrichment_role",
        "enrichment_hybridisation",
        "enrichment_flags",
        "ff_partial_charge",
        "ff_pb_radius",
        "coulomb_aromatic_E_proj",
        "coulomb_aromatic_n_src",
        "hbond_nearest_dir",
        "hbond_nearest_tensor",
        "hbond_nearest_spherical",
        "hbond_flags",
        "mc_nearest_co_dir",
        "mc_nearest_co_midpoint",
        "mc_nearest_co_T2",
        "mc_nearest_cn_T2",
        "larsen_corner_imputed",
        "mopac_bond_neighbors",
        "ring_direction_to_center",
        "disp_per_ring_tensor",
        "disp_per_ring_spherical",
        "piquad_quad_scalar",
    };

#ifdef NMR_TEST_DATA_DIR
    const fs::path repo_root = fs::path(NMR_TEST_DATA_DIR).parent_path().parent_path();
    ExpectCatalogAndSdkPresence(repo_root, new_stems);
#endif

    const fs::path out_dir = nmr::test::TestEnvironment::BaselineFeatures();
    fs::create_directories(out_dir);
    for (const auto& stem : new_stems) {
        fs::remove(out_dir / (stem + ".npy"));
    }

    // Write all features
    int arrays = ConformationResult::WriteAllFeatures(conf, out_dir.string());
    EXPECT_GT(arrays, 25) << "Expected 25+ arrays from 8 calculators + identity";

    std::cout << "\n  WriteFeatures baseline written to " << nmr::test::TestEnvironment::BaselineFeatures() << "\n"
              << "  Arrays: " << arrays << "\n"
              << "  Atoms: " << conf.AtomCount() << "\n"
              << "  Results attached: " << result.attached.size() << "\n";

    // List what we wrote
    for (const auto& entry : fs::directory_iterator(nmr::test::TestEnvironment::BaselineFeatures())) {
        if (!entry.is_regular_file()) continue;
        auto sz = fs::file_size(entry);
        std::cout << "    " << entry.path().filename().string()
                  << " (" << sz << " bytes)\n";
    }

    const size_t N = conf.AtomCount();

    if (conf.HasResult<EnrichmentResult>()) {
        auto role = ReadNpy(out_dir / "enrichment_role.npy");
        ExpectShapeAndDescr(role, {N}, "<i4");
        const int32_t* v = DataAs<int32_t>(role);
        for (size_t i = 0; i < N; ++i)
            EXPECT_EQ(v[i], static_cast<int32_t>(conf.AtomAt(i).role));

        auto hybrid = ReadNpy(out_dir / "enrichment_hybridisation.npy");
        ExpectShapeAndDescr(hybrid, {N}, "<i4");
        v = DataAs<int32_t>(hybrid);
        for (size_t i = 0; i < N; ++i)
            EXPECT_EQ(v[i], static_cast<int32_t>(conf.AtomAt(i).hybridisation));

        auto flags = ReadNpy(out_dir / "enrichment_flags.npy");
        ExpectShapeAndDescr(flags, {N, 8}, "|i1");
        const int8_t* f = DataAs<int8_t>(flags);
        for (size_t i = 0; i < N; ++i) {
            const auto& ca = conf.AtomAt(i);
            EXPECT_EQ(f[i*8 + 0], ca.is_backbone ? 1 : 0);
            EXPECT_EQ(f[i*8 + 1], ca.is_amide_H ? 1 : 0);
            EXPECT_EQ(f[i*8 + 2], ca.is_alpha_H ? 1 : 0);
            EXPECT_EQ(f[i*8 + 3], ca.is_methyl ? 1 : 0);
            EXPECT_EQ(f[i*8 + 4], ca.is_aromatic_H ? 1 : 0);
            EXPECT_EQ(f[i*8 + 5], ca.is_hbond_donor ? 1 : 0);
            EXPECT_EQ(f[i*8 + 6], ca.is_hbond_acceptor ? 1 : 0);
            EXPECT_EQ(f[i*8 + 7], ca.is_on_aromatic_residue ? 1 : 0);
        }
    }

    if (conf.HasResult<ChargeAssignmentResult>()) {
        auto charges = ReadNpy(out_dir / "ff_partial_charge.npy");
        ExpectShapeAndDescr(charges, {N}, "<f8");
        const double* q = DataAs<double>(charges);
        for (size_t i = 0; i < N; ++i)
            ExpectDoubleEqOrNan(q[i], conf.AtomAt(i).partial_charge);

        auto radii = ReadNpy(out_dir / "ff_pb_radius.npy");
        ExpectShapeAndDescr(radii, {N}, "<f8");
        const double* r = DataAs<double>(radii);
        for (size_t i = 0; i < N; ++i)
            ExpectDoubleEqOrNan(r[i], conf.AtomAt(i).pb_radius);
    }

    {
        size_t rows = 0;
        for (size_t i = 0; i < N; ++i) rows += conf.AtomAt(i).ring_neighbours.size();
        if (rows > 0) {
            auto dir = ReadNpy(out_dir / "ring_direction_to_center.npy");
            auto disp_tensor = ReadNpy(out_dir / "disp_per_ring_tensor.npy");
            auto disp_spherical = ReadNpy(out_dir / "disp_per_ring_spherical.npy");
            auto quad_scalar = ReadNpy(out_dir / "piquad_quad_scalar.npy");
            ExpectShapeAndDescr(dir, {rows, 3}, "<f8");
            ExpectShapeAndDescr(disp_tensor, {rows, 9}, "<f8");
            ExpectShapeAndDescr(disp_spherical, {rows, 9}, "<f8");
            ExpectShapeAndDescr(quad_scalar, {rows}, "<f8");

            const double* dv = DataAs<double>(dir);
            const double* mt = DataAs<double>(disp_tensor);
            const double* st = DataAs<double>(disp_spherical);
            const double* qs = DataAs<double>(quad_scalar);
            size_t row = 0;
            for (size_t i = 0; i < N; ++i) {
                for (const auto& rn : conf.AtomAt(i).ring_neighbours) {
                    ExpectDoubleEqOrNan(dv[row*3 + 0], rn.direction_to_center.x());
                    ExpectDoubleEqOrNan(dv[row*3 + 1], rn.direction_to_center.y());
                    ExpectDoubleEqOrNan(dv[row*3 + 2], rn.direction_to_center.z());
                    double packed_mat[9];
                    PackMat3RowMajor(rn.disp_tensor, packed_mat);
                    double packed_sph[9];
                    rn.disp_spherical.PackFull9(packed_sph);
                    for (int c = 0; c < 9; ++c) {
                        ExpectDoubleEqOrNan(mt[row*9 + c], packed_mat[c]);
                        ExpectDoubleEqOrNan(st[row*9 + c], packed_sph[c]);
                    }
                    ExpectDoubleEqOrNan(qs[row], rn.quad_scalar);
                    ++row;
                }
            }
            EXPECT_EQ(row, rows);
        }
    }

    if (conf.HasResult<CoulombResult>()) {
        auto proj = ReadNpy(out_dir / "coulomb_aromatic_E_proj.npy");
        ExpectShapeAndDescr(proj, {N}, "<f8");
        const double* p = DataAs<double>(proj);
        for (size_t i = 0; i < N; ++i)
            ExpectDoubleEqOrNan(p[i], conf.AtomAt(i).aromatic_E_bond_proj);

        auto count = ReadNpy(out_dir / "coulomb_aromatic_n_src.npy");
        ExpectShapeAndDescr(count, {N}, "<i4");
        const int32_t* c = DataAs<int32_t>(count);
        for (size_t i = 0; i < N; ++i)
            EXPECT_EQ(c[i], static_cast<int32_t>(conf.AtomAt(i).aromatic_n_sidechain_atoms));
    }

    if (conf.HasResult<HBondResult>()) {
        auto dir = ReadNpy(out_dir / "hbond_nearest_dir.npy");
        auto tensor = ReadNpy(out_dir / "hbond_nearest_tensor.npy");
        auto spherical = ReadNpy(out_dir / "hbond_nearest_spherical.npy");
        auto flags = ReadNpy(out_dir / "hbond_flags.npy");
        ExpectShapeAndDescr(dir, {N, 3}, "<f8");
        ExpectShapeAndDescr(tensor, {N, 9}, "<f8");
        ExpectShapeAndDescr(spherical, {N, 9}, "<f8");
        ExpectShapeAndDescr(flags, {N, 3}, "|i1");
        const double* dv = DataAs<double>(dir);
        const double* tv = DataAs<double>(tensor);
        const double* sv = DataAs<double>(spherical);
        const int8_t* fv = DataAs<int8_t>(flags);
        for (size_t i = 0; i < N; ++i) {
            const auto& ca = conf.AtomAt(i);
            ExpectDoubleEqOrNan(dv[i*3 + 0], ca.hbond_nearest_dir.x());
            ExpectDoubleEqOrNan(dv[i*3 + 1], ca.hbond_nearest_dir.y());
            ExpectDoubleEqOrNan(dv[i*3 + 2], ca.hbond_nearest_dir.z());
            double packed[9];
            PackMat3RowMajor(ca.hbond_nearest_tensor, packed);
            for (int c = 0; c < 9; ++c)
                ExpectDoubleEqOrNan(tv[i*9 + c], packed[c]);
            ca.hbond_nearest_spherical.PackFull9(packed);
            for (int c = 0; c < 9; ++c)
                ExpectDoubleEqOrNan(sv[i*9 + c], packed[c]);
            EXPECT_EQ(fv[i*3 + 0], ca.hbond_is_backbone ? 1 : 0);
            EXPECT_EQ(fv[i*3 + 1], ca.hbond_is_donor ? 1 : 0);
            EXPECT_EQ(fv[i*3 + 2], ca.hbond_is_acceptor ? 1 : 0);
        }
    }

    if (conf.HasResult<McConnellResult>()) {
        auto co_dir = ReadNpy(out_dir / "mc_nearest_co_dir.npy");
        auto co_mid = ReadNpy(out_dir / "mc_nearest_co_midpoint.npy");
        auto co_t2 = ReadNpy(out_dir / "mc_nearest_co_T2.npy");
        auto cn_t2 = ReadNpy(out_dir / "mc_nearest_cn_T2.npy");
        ExpectShapeAndDescr(co_dir, {N, 3}, "<f8");
        ExpectShapeAndDescr(co_mid, {N, 3}, "<f8");
        ExpectShapeAndDescr(co_t2, {N, 9}, "<f8");
        ExpectShapeAndDescr(cn_t2, {N, 9}, "<f8");
        const double* cod = DataAs<double>(co_dir);
        const double* com = DataAs<double>(co_mid);
        const double* cot = DataAs<double>(co_t2);
        const double* cnt = DataAs<double>(cn_t2);
        for (size_t i = 0; i < N; ++i) {
            const auto& ca = conf.AtomAt(i);
            ExpectDoubleEqOrNan(cod[i*3 + 0], ca.dir_nearest_CO.x());
            ExpectDoubleEqOrNan(cod[i*3 + 1], ca.dir_nearest_CO.y());
            ExpectDoubleEqOrNan(cod[i*3 + 2], ca.dir_nearest_CO.z());
            ExpectDoubleEqOrNan(com[i*3 + 0], ca.nearest_CO_midpoint.x());
            ExpectDoubleEqOrNan(com[i*3 + 1], ca.nearest_CO_midpoint.y());
            ExpectDoubleEqOrNan(com[i*3 + 2], ca.nearest_CO_midpoint.z());
            double co_pack[9], cn_pack[9];
            ca.T2_CO_nearest.PackFull9(co_pack);
            ca.T2_CN_nearest.PackFull9(cn_pack);
            for (int c = 0; c < 9; ++c) {
                ExpectDoubleEqOrNan(cot[i*9 + c], co_pack[c]);
                ExpectDoubleEqOrNan(cnt[i*9 + c], cn_pack[c]);
            }
        }

    }

    if (conf.HasResult<MopacResult>()) {
        size_t rows = 0;
        for (size_t i = 0; i < N; ++i) rows += conf.AtomAt(i).mopac_bond_neighbours.size();
        auto arr = ReadNpy(out_dir / "mopac_bond_neighbors.npy");
        ExpectShapeAndDescr(arr, {rows, 4}, "<f8");
        const double* d = DataAs<double>(arr);
        size_t row = 0;
        for (size_t i = 0; i < N; ++i) {
            for (const auto& nb : conf.AtomAt(i).mopac_bond_neighbours) {
                EXPECT_DOUBLE_EQ(d[row*4 + 0], static_cast<double>(i));
                EXPECT_DOUBLE_EQ(d[row*4 + 1], static_cast<double>(nb.other_atom));
                ExpectDoubleEqOrNan(d[row*4 + 2], nb.wiberg_order);
                const double topo = nb.topology_bond_index == SIZE_MAX
                    ? -1.0
                    : static_cast<double>(nb.topology_bond_index);
                EXPECT_DOUBLE_EQ(d[row*4 + 3], topo);
                ++row;
            }
        }
        EXPECT_EQ(row, rows);
    }

    if (conf.HasResult<LarsenHBondShieldingResult>()) {
        auto arr = ReadNpy(out_dir / "larsen_corner_imputed.npy");
        ExpectShapeAndDescr(arr, {N}, "|i1");
        const int8_t* d = DataAs<int8_t>(arr);
        for (size_t i = 0; i < N; ++i)
            EXPECT_EQ(d[i], conf.AtomAt(i).larsen_hbond_any_corner_imputed ? 1 : 0);
    }
}

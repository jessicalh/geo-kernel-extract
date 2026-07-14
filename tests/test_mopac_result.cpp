#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include "PdbFileReader.h"
#include "MopacResult.h"
#include "RuntimeEnvironment.h"
#include "DirectionalTestHelpers.h"
#include <filesystem>
#include <cmath>
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

using namespace nmr;

#ifndef NMR_TEST_PYTHON_EXECUTABLE
#error "NMR_TEST_PYTHON_EXECUTABLE must be defined"
#endif
#ifndef NMR_NPY_ALLOW_PICKLE_FALSE_SCRIPT
#error "NMR_NPY_ALLOW_PICKLE_FALSE_SCRIPT must be defined"
#endif

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

NpyArray ReadNpy(const std::filesystem::path& path) {
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

void ExpectDoubleEqOrNan(double actual, double expected) {
    if (std::isnan(expected)) {
        EXPECT_TRUE(std::isnan(actual));
    } else {
        EXPECT_DOUBLE_EQ(actual, expected);
    }
}

void CleanupMopacWriteFeaturesDir(const std::filesystem::path& out_dir) {
    static const char* files[] = {
        "mopac_charges.npy",
        "mopac_scalars.npy",
        "mopac_bond_orders.npy",
        "mopac_bond_neighbors.npy",
        "mopac_global.npy",
        "mopac_atom_populations.npy",
        "mopac_atomic_orbital_populations.npy",
        "mopac_atomic_orbital_population_totals.npy",
        "mopac_bond_valencies.npy",
        "mopac_bond_orders_unique.npy",
        "mopac_topology_bond_orders_full.npy",
    };
    std::error_code ec;
    for (const char* file : files) {
        std::filesystem::remove(out_dir / file, ec);
        ec.clear();
    }
    std::filesystem::remove(out_dir, ec);
}

}  // namespace


class MopacResultTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!std::filesystem::exists(nmr::test::TestEnvironment::UbqProtonated())) {
            GTEST_SKIP() << "1UBQ PDB not found";
        }
        auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
        if (!r.Ok()) GTEST_SKIP() << r.error;
        protein = std::move(r.protein);

        RuntimeEnvironment::Load();
    }
    std::unique_ptr<Protein> protein;
};


TEST_F(MopacResultTest, ComputeOnFullProtein) {
    auto& conf = protein->Conformation();

    auto result = MopacResult::Compute(conf, 0);
    ASSERT_NE(result, nullptr) << "MOPAC computation failed";

    // Verify charges are non-zero and reasonable
    bool any_nonzero = false;
    for (size_t i = 0; i < conf.AtomCount(); ++i) {
        double q = result->ChargeAt(i);
        EXPECT_FALSE(std::isnan(q)) << "NaN charge at atom " << i;
        EXPECT_GT(q, -3.0) << "Charge too negative at atom " << i;
        EXPECT_LT(q, 3.0) << "Charge too positive at atom " << i;
        if (std::abs(q) > 1e-6) any_nonzero = true;
    }
    EXPECT_TRUE(any_nonzero) << "All MOPAC charges are zero";

    // Heat of formation should be a large negative number for a protein
    EXPECT_LT(result->HeatOfFormation(), 0.0)
        << "Heat of formation should be negative for a protein";
}


TEST_F(MopacResultTest, BondOrdersReasonable) {
    auto& conf = protein->Conformation();

    auto result = MopacResult::Compute(conf, 0);
    ASSERT_NE(result, nullptr);

    const auto& orders = result->AllBondOrders();
    EXPECT_GT(orders.size(), 0u) << "No bond orders reported";

    for (const auto& bo : orders) {
        EXPECT_GE(bo.wiberg_order, 0.01)
            << "Bond order below threshold between "
            << bo.atom_a << " and " << bo.atom_b;
        EXPECT_LE(bo.wiberg_order, 4.0)
            << "Bond order too high between "
            << bo.atom_a << " and " << bo.atom_b;
    }
}


TEST_F(MopacResultTest, ChargesStoredOnConformationAtom) {
    auto& conf = protein->Conformation();

    auto result = MopacResult::Compute(conf, 0);
    ASSERT_NE(result, nullptr);

    for (size_t i = 0; i < conf.AtomCount(); ++i) {
        double ca_charge = conf.AtomAt(i).mopac_charge;
        double result_charge = result->ChargeAt(i);
        EXPECT_DOUBLE_EQ(ca_charge, result_charge)
            << "Charge mismatch at atom " << i;
    }
}


TEST_F(MopacResultTest, OrbitalPopulationsPresent) {
    auto& conf = protein->Conformation();

    auto result = MopacResult::Compute(conf, 0);
    ASSERT_NE(result, nullptr);

    // Heavy atoms should have non-zero s and p populations
    bool any_spop = false;
    bool any_ppop = false;
    for (size_t i = 0; i < conf.AtomCount(); ++i) {
        if (conf.AtomAt(i).mopac_s_pop > 0.1) any_spop = true;
        if (conf.AtomAt(i).mopac_p_pop > 0.1) any_ppop = true;
    }
    EXPECT_TRUE(any_spop) << "No atom has significant s-orbital population";
    EXPECT_TRUE(any_ppop) << "No atom has significant p-orbital population";
}

TEST_F(MopacResultTest, AtomicOrbitalPopulationTotalsWrittenFromSyntheticRow) {
    auto& conf = protein->Conformation();

    MopacResult result;
    auto& record = const_cast<MopacRunRecord&>(result.RunRecord());
    MopacAtomicOrbitalPopulation row;
    row.atom_index = 0;
    row.element = "C";
    row.populations = {
        1.0, 0.1, 0.2, 0.3, 0.01, 0.02, 0.03, 0.04, 0.05
    };
    record.atomic_orbital_populations.push_back(row);

    const auto out_dir = std::filesystem::temp_directory_path() /
        ("mopac_synthetic_ao_population_totals_" + std::to_string(::getpid()));
    CleanupMopacWriteFeaturesDir(out_dir);
    std::filesystem::create_directories(out_dir);

    result.WriteFeatures(conf, out_dir.string());
    ASSERT_EQ(nmr::test::directional::RunNumpyAllowPickleFalse(
                  NMR_TEST_PYTHON_EXECUTABLE,
                  NMR_NPY_ALLOW_PICKLE_FALSE_SCRIPT,
                  {out_dir / "mopac_global.npy",
                   out_dir / "mopac_atom_populations.npy",
                   out_dir / "mopac_atomic_orbital_populations.npy"}),
              0);

    auto raw = ReadNpy(out_dir / "mopac_atomic_orbital_populations.npy");
    auto totals = ReadNpy(out_dir / "mopac_atomic_orbital_population_totals.npy");
    ASSERT_EQ(raw.shape, (std::vector<size_t>{1, 9}));
    ASSERT_EQ(totals.shape, (std::vector<size_t>{1, 3}));

    const double* raw_data = Doubles(raw);
    const double* totals_data = Doubles(totals);
    const double expected_raw[] = {
        1.0, 0.1, 0.2, 0.3, 0.01, 0.02, 0.03, 0.04, 0.05
    };
    for (size_t c = 0; c < 9; ++c) {
        EXPECT_DOUBLE_EQ(raw_data[c], expected_raw[c]);
    }
    EXPECT_DOUBLE_EQ(totals_data[0], 1.0);
    EXPECT_DOUBLE_EQ(totals_data[1], 0.6);
    EXPECT_DOUBLE_EQ(totals_data[2], 0.15);

    CleanupMopacWriteFeaturesDir(out_dir);
}

TEST_F(MopacResultTest, AtomicOrbitalPopulationTotalsWrittenFromRawRows) {
    auto& conf = protein->Conformation();

    auto result = MopacResult::Compute(conf, 0);
    ASSERT_NE(result, nullptr);

    const auto& rows = result->AtomicOrbitalPopulations();
    ASSERT_GT(rows.size(), 0u) << "No printed AO population rows parsed";

    const auto out_dir = std::filesystem::temp_directory_path() /
        ("mopac_ao_population_totals_" + std::to_string(::getpid()));
    CleanupMopacWriteFeaturesDir(out_dir);
    std::filesystem::create_directories(out_dir);

    result->WriteFeatures(conf, out_dir.string());

    auto raw = ReadNpy(out_dir / "mopac_atomic_orbital_populations.npy");
    auto totals = ReadNpy(out_dir / "mopac_atomic_orbital_population_totals.npy");
    ASSERT_EQ(raw.shape, (std::vector<size_t>{rows.size(), 9}));
    ASSERT_EQ(totals.shape, (std::vector<size_t>{rows.size(), 3}));

    const double* raw_data = Doubles(raw);
    const double* totals_data = Doubles(totals);
    for (size_t i = 0; i < rows.size(); ++i) {
        for (size_t c = 0; c < 9; ++c) {
            ExpectDoubleEqOrNan(raw_data[i*9 + c], rows[i].populations[c]);
        }
        ExpectDoubleEqOrNan(totals_data[i*3 + 0], raw_data[i*9 + 0]);
        ExpectDoubleEqOrNan(totals_data[i*3 + 1],
                            raw_data[i*9 + 1] + raw_data[i*9 + 2] + raw_data[i*9 + 3]);
        ExpectDoubleEqOrNan(totals_data[i*3 + 2],
                            raw_data[i*9 + 4] + raw_data[i*9 + 5] + raw_data[i*9 + 6] +
                            raw_data[i*9 + 7] + raw_data[i*9 + 8]);
    }

    CleanupMopacWriteFeaturesDir(out_dir);
}


TEST_F(MopacResultTest, BondOrderAtomPairLookup) {
    auto& conf = protein->Conformation();

    auto result = MopacResult::Compute(conf, 0);
    ASSERT_NE(result, nullptr);

    // Pick a covalent bond from the topology and verify its MOPAC bond order
    if (protein->BondCount() > 0) {
        const Bond& bond = protein->BondAt(0);
        double topo_order = result->TopologyBondOrder(0);
        double pair_order = result->BondOrder(bond.atom_index_a, bond.atom_index_b);

        // Both paths should agree
        EXPECT_DOUBLE_EQ(topo_order, pair_order)
            << "Topology vs pair lookup mismatch for bond 0";

        // Covalent bonds should have meaningful bond orders
        // (most covalent bonds > 0.5, but some may be weak)
        if (topo_order > 0.0) {
            EXPECT_GT(topo_order, 0.1)
                << "Covalent bond has suspiciously low MOPAC bond order";
        }
    }

    // BondOrder for a non-existent pair should return 0
    EXPECT_DOUBLE_EQ(result->BondOrder(0, conf.AtomCount() + 100), 0.0);
}


TEST_F(MopacResultTest, ValencyReasonable) {
    auto& conf = protein->Conformation();

    auto result = MopacResult::Compute(conf, 0);
    ASSERT_NE(result, nullptr);

    for (size_t i = 0; i < conf.AtomCount(); ++i) {
        double val = result->ValencyAt(i);
        EXPECT_GE(val, 0.0) << "Negative valency at atom " << i;
        // Carbon typically ~4, oxygen ~2, hydrogen ~1, nitrogen ~3
        EXPECT_LT(val, 6.0) << "Valency too high at atom " << i
            << " (" << val << ")";
    }
}


TEST_F(MopacResultTest, MopacBondNeighboursSorted) {
    auto& conf = protein->Conformation();

    auto result = MopacResult::Compute(conf, 0);
    ASSERT_NE(result, nullptr);

    // Verify that bond neighbours are sorted descending by wiberg_order
    for (size_t i = 0; i < conf.AtomCount(); ++i) {
        const auto& nbs = conf.AtomAt(i).mopac_bond_neighbours;
        for (size_t j = 1; j < nbs.size(); ++j) {
            EXPECT_GE(nbs[j-1].wiberg_order, nbs[j].wiberg_order)
                << "Bond neighbours not sorted at atom " << i
                << " indices " << (j-1) << "," << j;
        }
    }
}

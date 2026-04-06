#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include "PdbFileReader.h"
#include "MopacResult.h"
#include "RuntimeEnvironment.h"
#include <filesystem>
#include <cmath>

using namespace nmr;


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

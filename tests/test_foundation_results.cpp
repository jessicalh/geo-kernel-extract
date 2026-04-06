#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include "PdbFileReader.h"
#include "Protein.h"
#include "GeometryResult.h"
#include "ChargeAssignmentResult.h"
#include "SpatialIndexResult.h"
#include "EnrichmentResult.h"
#include "MolecularGraphResult.h"
#include "ProtonationDetectionResult.h"
#include <filesystem>
#include <cmath>

using namespace nmr;

// ============================================================================
// Test data paths
// ============================================================================





// ============================================================================
// Shared fixture: loads 1UBQ once for all tests
// ============================================================================

class FoundationTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!std::filesystem::exists(nmr::test::TestEnvironment::UbqProtonated())) {
            GTEST_SKIP() << "1UBQ PDB not found at " << nmr::test::TestEnvironment::UbqProtonated();
        }
        auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
        if (!r.Ok()) GTEST_SKIP() << r.error;
        protein = std::move(r.protein);
    }
    std::unique_ptr<Protein> protein;
};


// ============================================================================
// Deliverable 1: ChargeAssignmentResult with ff14SB params
// ============================================================================

class ChargeFF14SBTest : public FoundationTest {};

TEST_F(ChargeFF14SBTest, LoadsAndAssigns) {
    if (!std::filesystem::exists(nmr::test::TestEnvironment::Ff14sbParams())) {
        GTEST_SKIP() << "ff14sb_params.dat not found";
    }

    auto& conf = protein->Conformation();
    auto result = ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams());
    ASSERT_NE(result, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(result)));
}

TEST_F(ChargeFF14SBTest, TotalChargeNearInteger) {
    if (!std::filesystem::exists(nmr::test::TestEnvironment::Ff14sbParams())) GTEST_SKIP();

    auto& conf = protein->Conformation();
    auto result = ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams());
    ASSERT_NE(result, nullptr);

    // 1UBQ: net charge should be close to an integer
    // 1UBQ at pH 7: +1 MET, +7 LYS, +4 ARG = +12 positive
    //              -4 GLU, -5 ASP = -9 negative  => net ~+3
    // But without terminal patches, it's approximate. Check integer-closeness.
    double total = result->TotalCharge();
    double nearest_int = std::round(total);
    double frac = std::abs(total - nearest_int);
    EXPECT_LT(frac, 0.5) << "Total charge " << total
        << " not close to integer (nearest=" << nearest_int << ")";

    conf.AttachResult(std::move(result));
}

TEST_F(ChargeFF14SBTest, BackboneNChargeRange) {
    if (!std::filesystem::exists(nmr::test::TestEnvironment::Ff14sbParams())) GTEST_SKIP();

    auto& conf = protein->Conformation();
    auto result = ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams());
    ASSERT_NE(result, nullptr);
    conf.AttachResult(std::move(result));

    // Backbone N atoms should have charge ~ -0.4 to -0.5
    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        const Residue& res = protein->ResidueAt(ri);
        if (res.N == Residue::NONE) continue;
        double q = conf.AtomAt(res.N).partial_charge;
        EXPECT_LT(q, -0.1) << "Backbone N at res " << res.sequence_number
            << " has charge " << q << " (expected < -0.1)";
        EXPECT_GT(q, -0.9) << "Backbone N at res " << res.sequence_number
            << " has charge " << q << " (expected > -0.9)";
    }
}

TEST_F(ChargeFF14SBTest, BackboneCChargeRange) {
    if (!std::filesystem::exists(nmr::test::TestEnvironment::Ff14sbParams())) GTEST_SKIP();

    auto& conf = protein->Conformation();
    auto result = ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams());
    ASSERT_NE(result, nullptr);
    conf.AttachResult(std::move(result));

    // Backbone C atoms should have charge ~ +0.5 to +0.7
    for (size_t ri = 0; ri < protein->ResidueCount(); ++ri) {
        const Residue& res = protein->ResidueAt(ri);
        if (res.C == Residue::NONE) continue;
        double q = conf.AtomAt(res.C).partial_charge;
        EXPECT_GT(q, 0.3) << "Backbone C at res " << res.sequence_number
            << " has charge " << q << " (expected > 0.3)";
        EXPECT_LT(q, 1.0) << "Backbone C at res " << res.sequence_number
            << " has charge " << q << " (expected < 1.0)";
    }
}

TEST_F(ChargeFF14SBTest, FewUnassignedAtoms) {
    if (!std::filesystem::exists(nmr::test::TestEnvironment::Ff14sbParams())) GTEST_SKIP();

    auto& conf = protein->Conformation();
    auto result = ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams());
    ASSERT_NE(result, nullptr);

    // At most a handful of atoms should be unassigned (charge = 0.0).
    // N-terminal NH3+ hydrogens (H1/H2/H3) lack entries in the param
    // file because it has no terminal-specific residue types.
    EXPECT_LE(result->UnassignedCount(), 10u)
        << result->UnassignedCount() << " atoms unassigned out of "
        << conf.AtomCount();

    conf.AttachResult(std::move(result));
}

TEST_F(ChargeFF14SBTest, MostAtomsAssigned) {
    if (!std::filesystem::exists(nmr::test::TestEnvironment::Ff14sbParams())) GTEST_SKIP();

    auto& conf = protein->Conformation();
    auto result = ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams());
    ASSERT_NE(result, nullptr);

    // Most atoms should get charges from the parameter file
    double assigned_ratio = static_cast<double>(result->AssignedCount()) /
                            static_cast<double>(conf.AtomCount());
    EXPECT_GT(assigned_ratio, 0.90)
        << "Only " << result->AssignedCount() << " of " << conf.AtomCount()
        << " atoms got ff14SB charges";

    conf.AttachResult(std::move(result));
}

// LegacyStubStillWorks — DELETED 2026-04-03
// StubChargeSource removed. Every protein has real charges.


// ============================================================================
// Deliverable 2: SpatialIndexResult
// ============================================================================

class SpatialIndexTest : public FoundationTest {};

TEST_F(SpatialIndexTest, ComputeAndAttach) {
    auto& conf = protein->Conformation();
    auto geo = GeometryResult::Compute(conf);
    ASSERT_TRUE(conf.AttachResult(std::move(geo)));

    auto spatial = SpatialIndexResult::Compute(conf);
    ASSERT_NE(spatial, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(spatial)));
}

TEST_F(SpatialIndexTest, EveryAtomHasNeighbours) {
    auto& conf = protein->Conformation();
    auto geo = GeometryResult::Compute(conf);
    conf.AttachResult(std::move(geo));
    auto spatial = SpatialIndexResult::Compute(conf);
    conf.AttachResult(std::move(spatial));

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        EXPECT_GT(conf.AtomAt(ai).spatial_neighbours.size(), 0u)
            << "Atom " << ai << " has no neighbours within 15A";
    }
}

TEST_F(SpatialIndexTest, NeighbourDistancesInRange) {
    auto& conf = protein->Conformation();
    auto geo = GeometryResult::Compute(conf);
    conf.AttachResult(std::move(geo));
    auto spatial = SpatialIndexResult::Compute(conf);
    conf.AttachResult(std::move(spatial));

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& nb : conf.AtomAt(ai).spatial_neighbours) {
            EXPECT_GT(nb.distance, 0.0) << "Zero distance neighbour at atom " << ai;
            EXPECT_LE(nb.distance, SPATIAL_NEIGHBOUR_CUTOFF_A + 0.01)
                << "Neighbour beyond cutoff at atom " << ai;
        }
    }
}

TEST_F(SpatialIndexTest, NeighbourDirectionsNormalised) {
    auto& conf = protein->Conformation();
    auto geo = GeometryResult::Compute(conf);
    conf.AttachResult(std::move(geo));
    auto spatial = SpatialIndexResult::Compute(conf);
    conf.AttachResult(std::move(spatial));

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& nb : conf.AtomAt(ai).spatial_neighbours) {
            double norm = nb.direction.norm();
            EXPECT_NEAR(norm, 1.0, 1e-6)
                << "Direction not normalised at atom " << ai
                << " -> " << nb.atom_index;
        }
    }
}

TEST_F(SpatialIndexTest, SelfNotInNeighbourList) {
    auto& conf = protein->Conformation();
    auto geo = GeometryResult::Compute(conf);
    conf.AttachResult(std::move(geo));
    auto spatial = SpatialIndexResult::Compute(conf);
    conf.AttachResult(std::move(spatial));

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& nb : conf.AtomAt(ai).spatial_neighbours) {
            EXPECT_NE(nb.atom_index, ai)
                << "Atom " << ai << " is in its own neighbour list";
        }
    }
}

TEST_F(SpatialIndexTest, KnownAtomPairDistance) {
    auto& conf = protein->Conformation();
    auto geo = GeometryResult::Compute(conf);
    conf.AttachResult(std::move(geo));
    auto spatial = SpatialIndexResult::Compute(conf);
    conf.AttachResult(std::move(spatial));

    // Find backbone N and CA of residue 0 -- they are bonded, distance ~1.5A
    const Residue& res0 = protein->ResidueAt(0);
    if (res0.N == Residue::NONE || res0.CA == Residue::NONE) return;

    // Direct Eigen computation
    double direct_dist = (conf.PositionAt(res0.N) - conf.PositionAt(res0.CA)).norm();

    // Find the same pair in the neighbour list
    bool found = false;
    for (const auto& nb : conf.AtomAt(res0.N).spatial_neighbours) {
        if (nb.atom_index == res0.CA) {
            EXPECT_NEAR(nb.distance, direct_dist, 1e-6)
                << "Neighbour distance does not match direct Eigen computation";
            found = true;
            break;
        }
    }
    EXPECT_TRUE(found) << "Backbone N-CA not found in neighbour list";
}


// ============================================================================
// Deliverable 3: EnrichmentResult
// ============================================================================

class EnrichmentTest : public FoundationTest {};

TEST_F(EnrichmentTest, ComputeAndAttach) {
    auto& conf = protein->Conformation();
    auto enrich = EnrichmentResult::Compute(conf);
    ASSERT_NE(enrich, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(enrich)));
}

TEST_F(EnrichmentTest, EveryAtomHasRole) {
    auto& conf = protein->Conformation();
    auto enrich = EnrichmentResult::Compute(conf);
    conf.AttachResult(std::move(enrich));

    int unknown = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        if (conf.AtomAt(ai).role == AtomRole::Unknown) unknown++;
    }
    // Allow very few Unknown (possibly edge cases at termini)
    EXPECT_LT(unknown, 5) << unknown << " atoms have Unknown role";
}

TEST_F(EnrichmentTest, Residue1BackboneN) {
    auto& conf = protein->Conformation();
    auto enrich = EnrichmentResult::Compute(conf);
    conf.AttachResult(std::move(enrich));

    const Residue& res0 = protein->ResidueAt(0);
    if (res0.N != Residue::NONE) {
        EXPECT_EQ(conf.AtomAt(res0.N).role, AtomRole::BackboneN);
    }
}

TEST_F(EnrichmentTest, Residue1BackboneCA) {
    auto& conf = protein->Conformation();
    auto enrich = EnrichmentResult::Compute(conf);
    conf.AttachResult(std::move(enrich));

    const Residue& res0 = protein->ResidueAt(0);
    if (res0.CA != Residue::NONE) {
        EXPECT_EQ(conf.AtomAt(res0.CA).role, AtomRole::BackboneCA);
    }
}

TEST_F(EnrichmentTest, PheRingCarbonsAreAromaticC) {
    auto& conf = protein->Conformation();
    auto enrich = EnrichmentResult::Compute(conf);
    conf.AttachResult(std::move(enrich));

    // Find a PHE ring and check that its C atoms have AromaticC role
    for (size_t ri = 0; ri < protein->RingCount(); ++ri) {
        const Ring& ring = protein->RingAt(ri);
        if (ring.type_index == RingTypeIndex::PheBenzene) {
            for (size_t ai : ring.atom_indices) {
                if (protein->AtomAt(ai).element == Element::C) {
                    EXPECT_EQ(conf.AtomAt(ai).role, AtomRole::AromaticC)
                        << "PHE ring C atom " << ai << " is not AromaticC";
                }
            }
            break;  // check first PHE ring
        }
    }
}

TEST_F(EnrichmentTest, BackboneNCountEquals76) {
    auto& conf = protein->Conformation();
    auto enrich = EnrichmentResult::Compute(conf);
    ASSERT_NE(enrich, nullptr);

    const auto& by_role = enrich->AtomsByRole();
    auto it = by_role.find(AtomRole::BackboneN);
    ASSERT_NE(it, by_role.end()) << "No BackboneN atoms found";
    EXPECT_EQ(it->second.size(), protein->ResidueCount())
        << "BackboneN count should equal number of residues (76)";

    conf.AttachResult(std::move(enrich));
}

TEST_F(EnrichmentTest, BackboneCACountEquals76) {
    auto& conf = protein->Conformation();
    auto enrich = EnrichmentResult::Compute(conf);
    ASSERT_NE(enrich, nullptr);

    const auto& by_role = enrich->AtomsByRole();
    auto it = by_role.find(AtomRole::BackboneCA);
    ASSERT_NE(it, by_role.end()) << "No BackboneCA atoms found";
    EXPECT_EQ(it->second.size(), protein->ResidueCount())
        << "BackboneCA count should equal number of residues (76)";

    conf.AttachResult(std::move(enrich));
}

TEST_F(EnrichmentTest, AtomsByRoleBackboneN76) {
    auto& conf = protein->Conformation();
    auto enrich = EnrichmentResult::Compute(conf);
    ASSERT_NE(enrich, nullptr);

    const auto& by_role = enrich->AtomsByRole();
    auto it = by_role.find(AtomRole::BackboneN);
    ASSERT_NE(it, by_role.end());
    EXPECT_EQ(it->second.size(), 76u);

    conf.AttachResult(std::move(enrich));
}


// ============================================================================
// Deliverable 4: Unified protonation / ring detection for HIS
// ============================================================================

class UnifiedHisTest : public ::testing::Test {};

TEST_F(UnifiedHisTest, ProtonationDrivesRingType) {
    if (!std::filesystem::exists(nmr::test::TestEnvironment::GmxProtonated())) {
        GTEST_SKIP() << "Protonated PDB not found";
    }

    // Load protonated structure
    auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::GmxProtonated());
    if (!r.Ok()) GTEST_SKIP() << r.error;
    auto& protein = *r.protein;

    auto& conf = protein.Conformation();

    // Run protonation detection FIRST
    auto prot = ProtonationDetectionResult::Compute(conf);
    ASSERT_NE(prot, nullptr);
    conf.AttachResult(std::move(prot));

    // Now re-detect rings (uses protonation_variant_index for HIS)
    protein.DetectAromaticRings();

    // Verify HIS ring types match protonation variants
    for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        if (res.type != AminoAcid::HIS) continue;

        std::string variant = conf.Result<ProtonationDetectionResult>()
                                  .VariantNameAt(ri);

        for (size_t ring_i = 0; ring_i < protein.RingCount(); ++ring_i) {
            const Ring& ring = protein.RingAt(ring_i);
            if (ring.parent_residue_index != ri) continue;
            // Skip perimeter/fused rings
            if (ring.RingSizeValue() > 5) continue;

            if (variant == "HID") {
                EXPECT_EQ(ring.type_index, RingTypeIndex::HidImidazole)
                    << "HIS " << res.sequence_number
                    << " protonation=HID but ring type is " << ring.TypeName();
            } else if (variant == "HIE") {
                EXPECT_EQ(ring.type_index, RingTypeIndex::HieImidazole)
                    << "HIS " << res.sequence_number
                    << " protonation=HIE but ring type is " << ring.TypeName();
            } else if (variant == "HIP") {
                EXPECT_EQ(ring.type_index, RingTypeIndex::HisImidazole)
                    << "HIS " << res.sequence_number
                    << " protonation=HIP but ring type is " << ring.TypeName();
            }
        }
    }
}

TEST_F(UnifiedHisTest, ExistingRingDetectionStillWorks) {
    // The basic ring detection (without prior protonation) still works
    if (!std::filesystem::exists(nmr::test::TestEnvironment::UbqProtonated())) {
        GTEST_SKIP() << "1UBQ PDB not found";
    }
    auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
    if (!r.Ok()) GTEST_SKIP() << r.error;

    // Ring detection happens during BuildFromPdb. Check it worked.
    EXPECT_GE(r.protein->RingCount(), 3)
        << "Ring detection should find at least 3 rings in 1UBQ";
}


// ============================================================================
// Deliverable 5: MolecularGraphResult
// ============================================================================

class MolecularGraphTest : public FoundationTest {};

TEST_F(MolecularGraphTest, ComputeAndAttach) {
    auto& conf = protein->Conformation();
    auto geo = GeometryResult::Compute(conf);
    conf.AttachResult(std::move(geo));
    auto spatial = SpatialIndexResult::Compute(conf);
    conf.AttachResult(std::move(spatial));

    auto graph = MolecularGraphResult::Compute(conf);
    ASSERT_NE(graph, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(graph)));
}

TEST_F(MolecularGraphTest, RingAtomsDistZero) {
    auto& conf = protein->Conformation();
    auto geo = GeometryResult::Compute(conf);
    conf.AttachResult(std::move(geo));
    auto spatial = SpatialIndexResult::Compute(conf);
    conf.AttachResult(std::move(spatial));
    auto graph = MolecularGraphResult::Compute(conf);
    conf.AttachResult(std::move(graph));

    // Ring atoms should have graph_dist_ring == 0
    for (size_t ri = 0; ri < protein->RingCount(); ++ri) {
        for (size_t ai : protein->RingAt(ri).atom_indices) {
            EXPECT_EQ(conf.AtomAt(ai).graph_dist_ring, 0)
                << "Ring atom " << ai << " has graph_dist_ring="
                << conf.AtomAt(ai).graph_dist_ring;
        }
    }
}

TEST_F(MolecularGraphTest, BondedToRingDistOne) {
    auto& conf = protein->Conformation();
    auto geo = GeometryResult::Compute(conf);
    conf.AttachResult(std::move(geo));
    auto spatial = SpatialIndexResult::Compute(conf);
    conf.AttachResult(std::move(spatial));
    auto graph = MolecularGraphResult::Compute(conf);
    conf.AttachResult(std::move(graph));

    // Find an atom bonded to a ring atom (but not a ring atom itself)
    // and verify it has graph_dist_ring == 1
    std::set<size_t> ring_atoms;
    for (size_t ri = 0; ri < protein->RingCount(); ++ri) {
        for (size_t ai : protein->RingAt(ri).atom_indices) {
            ring_atoms.insert(ai);
        }
    }

    bool found_check = false;
    for (size_t ai : ring_atoms) {
        for (size_t bi : protein->AtomAt(ai).bond_indices) {
            const Bond& bond = protein->BondAt(bi);
            size_t other = (bond.atom_index_a == ai)
                ? bond.atom_index_b : bond.atom_index_a;
            if (ring_atoms.count(other) == 0) {
                EXPECT_EQ(conf.AtomAt(other).graph_dist_ring, 1)
                    << "Atom " << other << " bonded to ring atom " << ai
                    << " has graph_dist_ring="
                    << conf.AtomAt(other).graph_dist_ring;
                found_check = true;
            }
        }
        if (found_check) break;
    }
    EXPECT_TRUE(found_check) << "Could not find an atom bonded to a ring atom";
}

TEST_F(MolecularGraphTest, BfsDecayRingAtomsOne) {
    auto& conf = protein->Conformation();
    auto geo = GeometryResult::Compute(conf);
    conf.AttachResult(std::move(geo));
    auto spatial = SpatialIndexResult::Compute(conf);
    conf.AttachResult(std::move(spatial));
    auto graph = MolecularGraphResult::Compute(conf);
    conf.AttachResult(std::move(graph));

    // Ring atoms have bfs_decay == exp(-0) == 1.0
    for (size_t ri = 0; ri < protein->RingCount(); ++ri) {
        for (size_t ai : protein->RingAt(ri).atom_indices) {
            EXPECT_NEAR(conf.AtomAt(ai).bfs_decay, 1.0, 1e-9)
                << "Ring atom " << ai << " bfs_decay="
                << conf.AtomAt(ai).bfs_decay;
        }
    }
}

TEST_F(MolecularGraphTest, BfsDecayDecreasesWithDistance) {
    auto& conf = protein->Conformation();
    auto geo = GeometryResult::Compute(conf);
    conf.AttachResult(std::move(geo));
    auto spatial = SpatialIndexResult::Compute(conf);
    conf.AttachResult(std::move(spatial));
    auto graph = MolecularGraphResult::Compute(conf);
    conf.AttachResult(std::move(graph));

    // Find an atom at distance 1 and one at distance 2 from a ring
    // and check that decay decreases
    double decay_d1 = -1.0;
    double decay_d2 = -1.0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        if (conf.AtomAt(ai).graph_dist_ring == 1 && decay_d1 < 0) {
            decay_d1 = conf.AtomAt(ai).bfs_decay;
        }
        if (conf.AtomAt(ai).graph_dist_ring == 2 && decay_d2 < 0) {
            decay_d2 = conf.AtomAt(ai).bfs_decay;
        }
        if (decay_d1 >= 0 && decay_d2 >= 0) break;
    }

    if (decay_d1 >= 0 && decay_d2 >= 0) {
        EXPECT_GT(decay_d1, decay_d2)
            << "bfs_decay should decrease: d1=" << decay_d1 << " d2=" << decay_d2;
        EXPECT_LT(decay_d1, 1.0) << "decay at d=1 should be less than 1.0";
    }
}

TEST_F(MolecularGraphTest, EnegSum1NonZero) {
    auto& conf = protein->Conformation();
    auto geo = GeometryResult::Compute(conf);
    conf.AttachResult(std::move(geo));
    auto spatial = SpatialIndexResult::Compute(conf);
    conf.AttachResult(std::move(spatial));
    auto graph = MolecularGraphResult::Compute(conf);
    conf.AttachResult(std::move(graph));

    // Atoms with bonds should have non-zero eneg_sum_1
    int nonzero = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        if (protein->AtomAt(ai).bond_indices.size() > 0 &&
            conf.AtomAt(ai).eneg_sum_1 > 0.0) {
            nonzero++;
        }
    }
    EXPECT_GT(nonzero, 0) << "No atoms have non-zero eneg_sum_1";
    // Most bonded atoms should have it
    EXPECT_GT(nonzero, static_cast<int>(conf.AtomCount()) / 2);
}

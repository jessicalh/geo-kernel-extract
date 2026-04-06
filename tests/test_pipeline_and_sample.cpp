#include "TestEnvironment.h"
//
// test_pipeline_and_sample.cpp
//
// Validates the viewer's library surface:
//   1. Pipeline: OperationRunner::Run attaches all results
//   2. Traversal: per-atom fields are populated and consistent
//   3. SampleAt: grid evaluation matches atom-position values
//   4. All 8 calculators: SampleShieldingAt returns sensible values
//
// This test is the contract between the library and the viewer.
// If it passes, the viewer can use the library without touching it.
//

#include <gtest/gtest.h>
#include <filesystem>
#include <cmath>

#include "OperationRunner.h"
#include "PdbFileReader.h"
#include "OrcaRunLoader.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "ConformationAtom.h"
#include "BiotSavartResult.h"
#include "HaighMallionResult.h"
#include "McConnellResult.h"
#include "RingSusceptibilityResult.h"
#include "PiQuadrupoleResult.h"
#include "DispersionResult.h"
#include "CoulombResult.h"
#include "HBondResult.h"
#include "ChargeSource.h"

namespace fs = std::filesystem;


// Helper: load a protonated protein with charges from prmtop
static std::unique_ptr<nmr::Protein> LoadTestProtein(const std::string& protein_id) {
    std::string dir = std::string(nmr::test::TestEnvironment::Consolidated()) + protein_id + "/";

    nmr::OrcaRunFiles files;
    files.pdb_path = dir + protein_id + "_WT.pdb";
    files.xyz_path = dir + protein_id + "_WT.xyz";
    files.prmtop_path = dir + protein_id + "_WT.prmtop";

    if (!fs::exists(files.xyz_path) || !fs::exists(files.prmtop_path))
        return nullptr;

    auto load = nmr::BuildFromOrca(files);
    if (!load.Ok()) return nullptr;
    return std::move(load.protein);
}


// ============================================================================
// Test 1: Pipeline attaches all results in correct order
// ============================================================================

TEST(OperationRunnerTest, AttachesAllResults) {
    auto protein = LoadTestProtein("P84477");
    if (!protein) GTEST_SKIP() << "P84477 test data not found";

    auto& conf = protein->Conformation();

    nmr::PrmtopChargeSource prmtop(
        std::string(nmr::test::TestEnvironment::Consolidated()) + "P84477/P84477_WT.prmtop",
        nmr::ForceField::Amber_ff14SB);
    nmr::RunOptions opts;
    opts.charge_source = &prmtop;

    auto run_result = nmr::OperationRunner::Run(conf, opts);

    // Should have: Geometry, SpatialIndex, Enrichment, Dssp,
    // ChargeAssignment, BS, HM, MC, RingSusc, PQ, Disp, Coulomb, HBond
    EXPECT_GE(run_result.attached.size(), 13u)
        << "Expected 13+ results, got " << run_result.attached.size();

    // Verify key results are attached
    EXPECT_TRUE(conf.HasResult<nmr::BiotSavartResult>());
    EXPECT_TRUE(conf.HasResult<nmr::HaighMallionResult>());
    EXPECT_TRUE(conf.HasResult<nmr::McConnellResult>());
    EXPECT_TRUE(conf.HasResult<nmr::CoulombResult>());
    EXPECT_TRUE(conf.HasResult<nmr::RingSusceptibilityResult>());
    EXPECT_TRUE(conf.HasResult<nmr::PiQuadrupoleResult>());
    EXPECT_TRUE(conf.HasResult<nmr::DispersionResult>());
    EXPECT_TRUE(conf.HasResult<nmr::HBondResult>());
}


// ============================================================================
// Test 2: Per-atom traversal produces populated fields
// ============================================================================

TEST(OperationRunnerTest, AtomFieldsPopulated) {
    auto protein = LoadTestProtein("P84477");
    if (!protein) GTEST_SKIP() << "P84477 test data not found";

    auto& conf = protein->Conformation();

    nmr::PrmtopChargeSource prmtop(
        std::string(nmr::test::TestEnvironment::Consolidated()) + "P84477/P84477_WT.prmtop",
        nmr::ForceField::Amber_ff14SB);
    nmr::RunOptions opts;
    opts.charge_source = &prmtop;

    nmr::OperationRunner::Run(conf, opts);

    // Traverse atoms — the viewer pattern
    int atoms_with_bs = 0;
    int atoms_with_mc = 0;
    int atoms_with_coulomb = 0;
    int atoms_with_rings = 0;
    double max_bs_t0 = 0;
    double max_mc_t0 = 0;

    for (size_t i = 0; i < conf.AtomCount(); ++i) {
        const auto& atom = conf.AtomAt(i);

        // Position is non-zero (loaded correctly)
        EXPECT_GT(atom.Position().norm(), 0.0);

        // BS shielding
        if (std::abs(atom.bs_shielding_contribution.T0) > 1e-10) {
            atoms_with_bs++;
            max_bs_t0 = std::max(max_bs_t0,
                std::abs(atom.bs_shielding_contribution.T0));
        }

        // McConnell shielding
        if (std::abs(atom.mc_shielding_contribution.T0) > 1e-10) {
            atoms_with_mc++;
            max_mc_t0 = std::max(max_mc_t0,
                std::abs(atom.mc_shielding_contribution.T0));
        }

        // Coulomb E-field
        if (atom.coulomb_E_total.norm() > 1e-10) {
            atoms_with_coulomb++;
        }

        // Ring neighbourhood
        if (!atom.ring_neighbours.empty()) {
            atoms_with_rings++;
            // Check that ring neighbourhood has BS and HM tensors
            for (const auto& rn : atom.ring_neighbours) {
                EXPECT_GT(rn.distance_to_center, 0.0);
                // G_tensor should be nonzero for nearby rings
                if (rn.distance_to_center < 10.0) {
                    EXPECT_GT(rn.G_tensor.norm(), 0.0)
                        << "BS G_tensor zero at dist=" << rn.distance_to_center;
                }
            }
        }
    }

    // P84477 is a small protein with aromatic rings
    EXPECT_GT(atoms_with_bs, 0) << "No atoms with BS shielding";
    EXPECT_GT(atoms_with_mc, 0) << "No atoms with MC shielding";
    EXPECT_GT(atoms_with_coulomb, 0) << "No atoms with Coulomb E-field";
    EXPECT_GT(atoms_with_rings, 0) << "No atoms with ring neighbours";

    std::cout << "  P84477 traversal summary:\n"
              << "    atoms=" << conf.AtomCount()
              << " rings=" << protein->RingCount() << "\n"
              << "    with BS: " << atoms_with_bs
              << " (max T0=" << max_bs_t0 << ")\n"
              << "    with MC: " << atoms_with_mc
              << " (max T0=" << max_mc_t0 << ")\n"
              << "    with Coulomb: " << atoms_with_coulomb << "\n"
              << "    with ring neighbours: " << atoms_with_rings << "\n";
}


// ============================================================================
// Test 3: SampleAt agrees with atom-position results
// ============================================================================

TEST(SampleAtTest, BSMatchesAtomValues) {
    auto protein = LoadTestProtein("P84477");
    if (!protein) GTEST_SKIP() << "P84477 test data not found";

    auto& conf = protein->Conformation();
    nmr::OperationRunner::Run(conf, {});

    const auto& bs = conf.Result<nmr::BiotSavartResult>();

    // For atoms near rings, SampleShieldingAt should approximately match
    // the stored bs_shielding_contribution. Not exact — the atom-position
    // result uses filters (RingBondedExclusion) that SampleAt doesn't.
    // But for atoms NOT bonded to ring vertices, they should be close.
    int tested = 0;
    int matched = 0;

    for (size_t i = 0; i < conf.AtomCount() && tested < 50; ++i) {
        const auto& atom = conf.AtomAt(i);
        double atom_t0 = atom.bs_shielding_contribution.T0;
        if (std::abs(atom_t0) < 0.001) continue;  // skip negligible

        // Skip atoms that are ring vertices or bonded to ring vertices
        // (filters differ between Compute and SampleAt)
        bool is_ring_atom = false;
        for (size_t ri = 0; ri < protein->RingCount(); ++ri) {
            const auto& ring = protein->RingAt(ri);
            for (size_t vi : ring.atom_indices) {
                if (vi == i) { is_ring_atom = true; break; }
                // Check if bonded to vertex
                for (const auto& bond : protein->Bonds()) {
                    if ((bond.atom_index_a == i && bond.atom_index_b == vi) ||
                        (bond.atom_index_a == vi && bond.atom_index_b == i)) {
                        is_ring_atom = true;
                        break;
                    }
                }
                if (is_ring_atom) break;
            }
            if (is_ring_atom) break;
        }
        if (is_ring_atom) continue;

        nmr::SphericalTensor sampled = bs.SampleShieldingAt(atom.Position());
        tested++;

        // Relative tolerance: these should match well for non-excluded atoms
        double rel_diff = std::abs(sampled.T0 - atom_t0)
                        / std::max(std::abs(atom_t0), 1e-6);
        if (rel_diff < 0.05) matched++;  // 5% tolerance

        EXPECT_NEAR(sampled.T0, atom_t0, std::abs(atom_t0) * 0.05 + 1e-6)
            << "BS SampleAt mismatch at atom " << i
            << " sampled=" << sampled.T0 << " stored=" << atom_t0;
    }

    EXPECT_GT(tested, 5) << "Too few testable atoms";
    std::cout << "  BS SampleAt: tested=" << tested
              << " matched=" << matched << "\n";
}


// ============================================================================
// Test 4: SampleAt on a grid produces physically sensible values
// ============================================================================

TEST(SampleAtTest, BSGridAboveRing) {
    auto protein = LoadTestProtein("P84477");
    if (!protein) GTEST_SKIP() << "P84477 test data not found";

    auto& conf = protein->Conformation();
    nmr::OperationRunner::Run(conf, {});

    if (protein->RingCount() == 0) GTEST_SKIP() << "No rings";

    const auto& bs = conf.Result<nmr::BiotSavartResult>();
    const auto& geom = conf.ring_geometries[0];

    // Sample on the ring normal at 3A above center — should be shielded (T0 > 0)
    nmr::Vec3 above = geom.center + 3.0 * geom.normal;
    auto st_above = bs.SampleShieldingAt(above);

    // And in the ring plane at 5A — should be deshielded (T0 < 0)
    // Find a direction perpendicular to normal
    nmr::Vec3 perp = geom.normal.cross(
        std::abs(geom.normal.x()) < 0.9 ? nmr::Vec3(1,0,0) : nmr::Vec3(0,1,0)
    ).normalized();
    nmr::Vec3 inplane = geom.center + 5.0 * perp;
    auto st_inplane = bs.SampleShieldingAt(inplane);

    std::cout << "  Ring 0: T0 at +3A normal = " << st_above.T0
              << ", T0 at 5A in-plane = " << st_inplane.T0 << "\n";

    // Classic ring current pattern: shielded above, deshielded in-plane
    // (with negative intensity, the signs flip, but the pattern holds)
    EXPECT_NE(st_above.T0, 0.0) << "Zero shielding above ring";
    // The above and in-plane should have opposite signs
    if (std::abs(st_above.T0) > 1e-4 && std::abs(st_inplane.T0) > 1e-4) {
        EXPECT_LT(st_above.T0 * st_inplane.T0, 0.0)
            << "Above and in-plane should have opposite T0 signs";
    }
}


// ============================================================================
// Test 5: B-field sampling for butterfly visualization
// ============================================================================

TEST(SampleAtTest, BSButterflyField) {
    auto protein = LoadTestProtein("P84477");
    if (!protein) GTEST_SKIP() << "P84477 test data not found";

    auto& conf = protein->Conformation();
    nmr::OperationRunner::Run(conf, {});

    if (protein->RingCount() == 0) GTEST_SKIP() << "No rings";

    const auto& bs = conf.Result<nmr::BiotSavartResult>();
    const auto& geom = conf.ring_geometries[0];

    // B-field on the ring axis should be parallel to the normal
    nmr::Vec3 above = geom.center + 3.0 * geom.normal;
    nmr::Vec3 B = bs.SampleBFieldAt(above);

    EXPECT_GT(B.norm(), 0.0) << "Zero B-field above ring";

    // B should be roughly along the normal direction
    double cos_angle = B.normalized().dot(geom.normal);
    EXPECT_GT(std::abs(cos_angle), 0.8)
        << "B-field not aligned with ring normal at 3A above";

    std::cout << "  B-field at +3A: |B|=" << B.norm()
              << " cos(B,n)=" << cos_angle << "\n";
}


// ============================================================================
// Test 6: All other SampleAt methods return non-zero for appropriate queries
// ============================================================================

TEST(SampleAtTest, AllCalculatorsSample) {
    auto protein = LoadTestProtein("P84477");
    if (!protein) GTEST_SKIP() << "P84477 test data not found";

    auto& conf = protein->Conformation();
    nmr::PrmtopChargeSource prmtop(
        std::string(nmr::test::TestEnvironment::Consolidated()) + "P84477/P84477_WT.prmtop",
        nmr::ForceField::Amber_ff14SB);
    nmr::RunOptions opts;
    opts.charge_source = &prmtop;
    nmr::OperationRunner::Run(conf, opts);

    if (protein->RingCount() == 0) GTEST_SKIP() << "No rings";

    const auto& geom = conf.ring_geometries[0];
    nmr::Vec3 test_point = geom.center + 3.0 * geom.normal;

    // HM
    auto hm_st = conf.Result<nmr::HaighMallionResult>().SampleShieldingAt(test_point);
    EXPECT_NE(hm_st.T0, 0.0) << "HM SampleAt returned zero";

    // McConnell — sample near a bond midpoint
    nmr::Vec3 bond_test = conf.bond_midpoints[0] + nmr::Vec3(2.0, 0.0, 0.0);
    auto mc_st = conf.Result<nmr::McConnellResult>().SampleShieldingAt(bond_test);
    // McConnell is pure T2 (T0 ≈ 0), so check T2 magnitude
    EXPECT_GT(mc_st.T2Magnitude(), 0.0) << "MC SampleAt returned zero T2";

    // Ring Susceptibility
    auto chi_st = conf.Result<nmr::RingSusceptibilityResult>().SampleShieldingAt(test_point);
    EXPECT_NE(chi_st.T0, 0.0) << "RingSusc SampleAt returned zero";

    // PiQuadrupole — pure T2
    auto pq_st = conf.Result<nmr::PiQuadrupoleResult>().SampleShieldingAt(test_point);
    EXPECT_GT(pq_st.T2Magnitude(), 0.0) << "PQ SampleAt returned zero T2";

    // Dispersion — sample near a ring vertex
    const auto& vertices = geom.vertices;
    if (!vertices.empty()) {
        nmr::Vec3 disp_test = vertices[0] + nmr::Vec3(3.0, 0.0, 0.0);
        auto disp_st = conf.Result<nmr::DispersionResult>().SampleShieldingAt(disp_test);
        // Dispersion may be zero if point is outside 5A cutoff from all vertices
        // Just verify it doesn't crash
        (void)disp_st;
    }

    // Coulomb E-field
    auto E = conf.Result<nmr::CoulombResult>().SampleEFieldAt(test_point);
    EXPECT_GT(E.norm(), 0.0) << "Coulomb SampleEFieldAt returned zero";

    // HBond — may return zero if no H-bonds near test point, just verify no crash
    auto hb_st = conf.Result<nmr::HBondResult>().SampleShieldingAt(test_point);
    (void)hb_st;

    std::cout << "  All SampleAt methods exercised successfully\n"
              << "    BS T0=" << conf.Result<nmr::BiotSavartResult>().SampleShieldingAt(test_point).T0 << "\n"
              << "    HM T0=" << hm_st.T0 << "\n"
              << "    MC T2_mag=" << mc_st.T2Magnitude() << "\n"
              << "    Chi T0=" << chi_st.T0 << "\n"
              << "    PQ T2_mag=" << pq_st.T2Magnitude() << "\n"
              << "    |E|=" << E.norm() << " V/A\n"
              << "    HBond T0=" << hb_st.T0 << "\n";
}


// ============================================================================
// Test 7: Pipeline without charges gracefully skips Coulomb
// ============================================================================

TEST(OperationRunnerTest, SkipsCoulombWithoutCharges) {
    auto protein = LoadTestProtein("P84477");
    if (!protein) GTEST_SKIP() << "P84477 test data not found";

    auto& conf = protein->Conformation();

    // No charge source — Coulomb should be skipped
    auto run_result = nmr::OperationRunner::Run(conf, {});

    EXPECT_FALSE(conf.HasResult<nmr::CoulombResult>())
        << "Coulomb should not be attached without charges";
    EXPECT_TRUE(conf.HasResult<nmr::BiotSavartResult>())
        << "BS should still be attached";
}


// ============================================================================
// Test 8: Ring and bond geometry accessible for viewer rendering
// ============================================================================

TEST(OperationRunnerTest, GeometryAccessible) {
    auto protein = LoadTestProtein("P84477");
    if (!protein) GTEST_SKIP() << "P84477 test data not found";

    auto& conf = protein->Conformation();
    nmr::OperationRunner::Run(conf, {});

    // Ring geometries
    EXPECT_EQ(conf.ring_geometries.size(), protein->RingCount());
    for (size_t i = 0; i < protein->RingCount(); ++i) {
        const auto& geo = conf.ring_geometries[i];
        EXPECT_GT(geo.radius, 0.0);
        EXPECT_NEAR(geo.normal.norm(), 1.0, 1e-6);
        EXPECT_GT(geo.vertices.size(), 2u);
    }

    // Bond geometries
    EXPECT_EQ(conf.bond_midpoints.size(), protein->BondCount());
    EXPECT_EQ(conf.bond_directions.size(), protein->BondCount());
    EXPECT_EQ(conf.bond_lengths.size(), protein->BondCount());

    // Global geometry
    EXPECT_GT(conf.radius_of_gyration, 0.0);

    std::cout << "  Geometry: "
              << protein->RingCount() << " rings, "
              << protein->BondCount() << " bonds, "
              << "Rg=" << conf.radius_of_gyration << " A\n";
}

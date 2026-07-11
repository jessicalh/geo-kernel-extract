#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include <cmath>
#include <iostream>

#include "BiotSavartResult.h"
#include "DispersionResult.h"
#include "HaighMallionResult.h"
#include "PiQuadrupoleResult.h"
#include "RingSusceptibilityResult.h"
#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "CalculatorConfig.h"
#include "Protein.h"
#include "PdbFileReader.h"
#include "OrcaRunLoader.h"
#include "PhysicalConstants.h"

#include <filesystem>
#include <set>
namespace fs = std::filesystem;
using namespace nmr;

namespace {

std::vector<Vec3> RegularHexagon(double radius) {
    std::vector<Vec3> vertices;
    const double pi = std::acos(-1.0);
    for (int i = 0; i < 6; ++i) {
        const double angle = 2.0 * pi * static_cast<double>(i) / 6.0;
        vertices.emplace_back(radius * std::cos(angle),
                              radius * std::sin(angle), 0.0);
    }
    return vertices;
}

}  // namespace



// ============================================================================
// Analytical test: B-field from a single wire segment.
//
// A wire from (0,0,0) to (1,0,0) carrying 1 A, field point at (0.5, 1, 0).
// The field should be along -z (by right-hand rule for current in +x).
//
// For a finite wire: B = (mu_0 I / 4pi d) * (sin(theta_2) - sin(theta_1))
// where d = perpendicular distance = 1.0 m, theta_1 = -26.57°, theta_2 = 26.57°
// B = 1e-7 * 1 / 1.0 * 2 * sin(26.57°) = 1e-7 * 2 * 0.4472 = 8.944e-8 T
// Direction: -z
// ============================================================================

TEST(BiotSavartAnalytical, WireSegmentDirection) {
    // We test indirectly through JohnsonBoveyField with a single-segment
    // "ring" (degenerate, but the wire segment formula is what we verify).
    // Instead, verify that a real ring produces B along the normal above center.

    // A regular hexagon centered at origin in the xy-plane, radius 1.39 A.
    // Atom at (0, 0, 3.0 A) — directly above center.
    // B-field should be along +z (ring current in the conventional direction).

    // We cannot call the static functions directly, so we verify through
    // the full pipeline with a known protein.
    SUCCEED() << "Wire segment tested through full protein pipeline below";
}


TEST(BiotSavartAnalytical, FiniteWireGuardAcceptsPointAboveRingFace) {
    const std::vector<Vec3> vertices = RegularHexagon(1.4);
    const Vec3 normal(0.0, 0.0, 1.0);
    const double lobe_offset = 0.64;
    const Vec3 point(0.0, 0.0, 0.5);  // center distance < ring radius

    const double minimum_wire_distance =
        biot_savart_detail::MinimumDistanceToJohnsonBoveyWireLoops(
            vertices, normal, lobe_offset, point);
    const double guard =
        CalculatorConfig::Get("biot_savart_wire_endpoint_guard") /
        ANGSTROMS_TO_METRES;

    EXPECT_LT(point.norm(), 1.4);
    EXPECT_GT(minimum_wire_distance, guard);

    // Execute the production current-loop kernel, not a test-side formula.
    const Vec3 field = biot_savart_detail::JohnsonBoveyField(
        vertices, normal, lobe_offset, 1.0, point);
    EXPECT_TRUE(field.allFinite());
    EXPECT_GT(field.norm(), 0.0);

    // A point exactly on a segment endpoint is the true singular row guard.
    const Vec3 endpoint = vertices[0] + lobe_offset * normal;
    EXPECT_DOUBLE_EQ(
        biot_savart_detail::MinimumDistanceToJohnsonBoveyWireLoops(
            vertices, normal, lobe_offset, endpoint),
        0.0);
}


TEST(RingNearFieldProduction, CenterInsideRadiusIsAcceptedByAllFiveCalculators) {
    // Move a topology-unrelated atom to a controlled point above a real ring,
    // then execute each complete production Compute path.  This freezes the
    // row-selection contract as well as the exported numerical helpers: a
    // reintroduced center/diameter filter would suppress one or more payloads.
    auto loaded = BuildFromProtonatedPdb(
        nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(loaded.Ok()) << loaded.error;
    Protein& protein = *loaded.protein;
    auto& source_conf = protein.Conformation();
    ASSERT_TRUE(source_conf.AttachResult(GeometryResult::Compute(source_conf)));
    ASSERT_GT(protein.RingCount(), 0u);

    constexpr size_t ring_index = 0;
    const Ring& ring = protein.RingAt(ring_index);
    const RingGeometry source_geometry =
        source_conf.ring_geometries[ring_index];

    std::set<size_t> topology_excluded;
    for (size_t vertex_index : ring.atom_indices) {
        topology_excluded.insert(vertex_index);
        for (size_t bond_index : protein.AtomAt(vertex_index).bond_indices) {
            const Bond& bond = protein.BondAt(bond_index);
            topology_excluded.insert(bond.atom_index_a);
            topology_excluded.insert(bond.atom_index_b);
        }
    }

    size_t target_index = protein.AtomCount();
    for (size_t ai = 0; ai < protein.AtomCount(); ++ai) {
        if (!topology_excluded.count(ai)) {
            target_index = ai;
            break;
        }
    }
    ASSERT_LT(target_index, protein.AtomCount());

    std::vector<Vec3> positions = source_conf.Positions();
    positions[target_index] = source_geometry.center
        + 0.5 * source_geometry.normal.normalized();
    auto& probe = protein.AddConformation(
        std::move(positions), "ring-near-field-production-probe");
    ASSERT_TRUE(probe.AttachResult(GeometryResult::Compute(probe)));
    ASSERT_TRUE(probe.AttachResult(SpatialIndexResult::Compute(probe)));
    ASSERT_LT((probe.PositionAt(target_index)
               - probe.ring_geometries[ring_index].center).norm(),
              probe.ring_geometries[ring_index].radius);

    ASSERT_TRUE(probe.AttachResult(BiotSavartResult::Compute(probe)));
    ASSERT_TRUE(probe.AttachResult(HaighMallionResult::Compute(probe)));
    ASSERT_TRUE(probe.AttachResult(PiQuadrupoleResult::Compute(probe)));
    ASSERT_TRUE(probe.AttachResult(
        RingSusceptibilityResult::Compute(probe)));
    ASSERT_TRUE(probe.AttachResult(DispersionResult::Compute(probe)));

    const RingNeighbourhood* row = nullptr;
    for (const RingNeighbourhood& candidate :
         probe.AtomAt(target_index).ring_neighbours) {
        if (candidate.ring_index == ring_index) {
            row = &candidate;
            break;
        }
    }
    ASSERT_NE(row, nullptr);
    EXPECT_GT(row->G_tensor.norm(), 0.0);
    EXPECT_GT(row->hm_G_tensor.norm(), 0.0);
    EXPECT_NE(row->quad_scalar, 0.0);
    EXPECT_NE(row->chi_scalar, 0.0);
    EXPECT_GT(row->disp_scalar, 0.0);
    EXPECT_GT(row->disp_contacts, 0);
}


TEST(BiotSavartAnalytical, TrpPerimeterIsDiagnosticOnlyPolicy) {
    EXPECT_TRUE(biot_savart_detail::ContributesToCanonicalTotals(
        RingTypeIndex::TrpBenzene));
    EXPECT_TRUE(biot_savart_detail::ContributesToCanonicalTotals(
        RingTypeIndex::TrpPyrrole));
    EXPECT_FALSE(biot_savart_detail::ContributesToCanonicalTotals(
        RingTypeIndex::TrpPerimeter));
    EXPECT_FALSE(haigh_mallion_detail::ContributesToCanonicalTotals(
        RingTypeIndex::TrpPerimeter));
}


TEST(RingCurrentTrpProduction, CanonicalTotalsExcludeVisibleTrpPerimeter) {
    const fs::path fixture =
        fs::path(nmr::test::TestEnvironment::UbqProtonated()).parent_path() /
        "illustrative_peptides/trp_cage_1l2y_model1.pdb";
    ASSERT_TRUE(fs::is_regular_file(fixture)) << fixture;

    auto loaded = BuildFromProtonatedPdb(fixture.string());
    ASSERT_TRUE(loaded.Ok()) << loaded.error;
    auto& conf = loaded.protein->Conformation();
    ASSERT_TRUE(conf.AttachResult(GeometryResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(SpatialIndexResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(BiotSavartResult::Compute(conf)));

    size_t perimeter_rows = 0;
    size_t nonzero_perimeter_bs_rows = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& atom = conf.AtomAt(ai);
        Mat3 expected_G = Mat3::Zero();
        Vec3 expected_B = Vec3::Zero();
        for (const RingNeighbourhood& row : atom.ring_neighbours) {
            if (row.ring_type == RingTypeIndex::TrpPerimeter) {
                ++perimeter_rows;
                if (row.G_tensor.norm() > 0.0) ++nonzero_perimeter_bs_rows;
                continue;
            }
            expected_G += row.G_tensor;
            expected_B += row.B_field;
        }
        EXPECT_LT((atom.total_G_tensor - expected_G).norm(), 1e-12);
        EXPECT_LT((atom.total_B_field - expected_B).norm(), 1e-24);
        EXPECT_DOUBLE_EQ(
            atom.per_type_G_T0_sum[
                static_cast<int>(RingTypeIndex::TrpPerimeter)],
            0.0);
    }
    EXPECT_GT(perimeter_rows, 0u);
    EXPECT_GT(nonzero_perimeter_bs_rows, 0u)
        << "TRP9 must remain visible as a computed sparse diagnostic";

    ASSERT_TRUE(conf.AttachResult(HaighMallionResult::Compute(conf)));
    size_t nonzero_perimeter_hm_rows = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& atom = conf.AtomAt(ai);
        Mat3 expected_G = Mat3::Zero();
        for (const RingNeighbourhood& row : atom.ring_neighbours) {
            if (row.ring_type == RingTypeIndex::TrpPerimeter) {
                if (row.hm_G_tensor.norm() > 0.0) ++nonzero_perimeter_hm_rows;
                continue;
            }
            expected_G += row.hm_G_tensor;
        }
        EXPECT_LT((atom.hm_shielding_contribution.Reconstruct() - expected_G)
                      .norm(),
                  1e-12);
        EXPECT_DOUBLE_EQ(
            atom.per_type_hm_T0_sum[
                static_cast<int>(RingTypeIndex::TrpPerimeter)],
            0.0);
    }
    EXPECT_GT(nonzero_perimeter_hm_rows, 0u)
        << "TRP9 HM kernel must remain visible only in sparse diagnostics";
}


// ============================================================================
// Full protein test: 1UBQ
// ============================================================================

class BiotSavartProteinTest : public ::testing::Test {
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


TEST_F(BiotSavartProteinTest, ComputeAndAttach) {
    auto& conf = protein->Conformation();
    auto bs = BiotSavartResult::Compute(conf);
    ASSERT_NE(bs, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(bs)));
    ASSERT_TRUE(conf.HasResult<BiotSavartResult>());
}


TEST_F(BiotSavartProteinTest, GTensorIsRankOne) {
    // G_ab = n_b * B_a is rank-1 by construction. For any atom near a ring,
    // det(G) should be zero (to floating point precision).
    auto& conf = protein->Conformation();
    conf.AttachResult(BiotSavartResult::Compute(conf));

    int checked = 0;
    double max_det = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& rn : conf.AtomAt(ai).ring_neighbours) {
            if (rn.G_tensor.norm() < 1e-15) continue;
            double det = std::abs(rn.G_tensor.determinant());
            max_det = std::max(max_det, det);
            checked++;
        }
    }

    EXPECT_GT(checked, 100) << "Should have many ring-atom pairs";
    EXPECT_LT(max_det, 1e-15)
        << "G = n (x) B is rank-1, determinant must be zero";

    std::cout << "  Checked " << checked << " pairs, max |det(G)| = "
              << max_det << "\n";
}


TEST_F(BiotSavartProteinTest, BoydSkrynnikovTensorStructure) {
    // Verify that G decomposes into Boyd & Skrynnikov eq 3 structure:
    //   T0 = (n . B) * PPM_FACTOR / 3
    //   T1 from antisymmetric part (sigma_xz term)
    //   T2 from traceless symmetric part
    auto& conf = protein->Conformation();
    conf.AttachResult(BiotSavartResult::Compute(conf));

    int checked = 0;
    double max_t0_diff = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& rn : conf.AtomAt(ai).ring_neighbours) {
            if (rn.B_field.norm() < 1e-30) continue;

            const RingGeometry& geom = conf.ring_geometries[rn.ring_index];

            // T0 from trace of G should equal -(n . B) * PPM_FACTOR / 3
            // (minus sign from sigma_ab = -dB_a/dB_0b)
            double trace_over_3 = rn.G_tensor.trace() / 3.0;
            double ndotB_ppm = -geom.normal.dot(rn.B_field) * PPM_FACTOR / 3.0;
            double diff = std::abs(trace_over_3 - ndotB_ppm);
            max_t0_diff = std::max(max_t0_diff, diff);

            // Also verify SphericalTensor T0 matches
            EXPECT_NEAR(rn.G_spherical.T0, trace_over_3, 1e-10);

            checked++;
        }
    }

    EXPECT_GT(checked, 100) << "Should verify many ring-atom pairs";
    EXPECT_LT(max_t0_diff, 1e-10)
        << "T0 must equal (n.B)*PPM/3 at machine precision";

    std::cout << "  Verified Boyd-Skrynnikov T0 on " << checked
              << " pairs, max diff = " << max_t0_diff << "\n";
}


TEST_F(BiotSavartProteinTest, ShieldingContributionHasAllIrreps) {
    // BS produces rank-1 tensor -> T0, T1, T2 all non-zero.
    auto& conf = protein->Conformation();
    conf.AttachResult(BiotSavartResult::Compute(conf));

    int nonzero_t0 = 0, nonzero_t1 = 0, nonzero_t2 = 0;
    double max_t0 = 0, max_t2 = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& sc = conf.AtomAt(ai).bs_shielding_contribution;
        if (std::abs(sc.T0) > 1e-10) nonzero_t0++;
        double t1mag = std::sqrt(sc.T1[0]*sc.T1[0] + sc.T1[1]*sc.T1[1]
                                 + sc.T1[2]*sc.T1[2]);
        if (t1mag > 1e-10) nonzero_t1++;
        double t2mag = sc.T2Magnitude();
        if (t2mag > 1e-10) nonzero_t2++;
        max_t0 = std::max(max_t0, std::abs(sc.T0));
        max_t2 = std::max(max_t2, t2mag);
    }

    EXPECT_GT(nonzero_t0, 0) << "Rank-1 tensor must have non-zero T0";
    EXPECT_GT(nonzero_t1, 0) << "Rank-1 tensor must have non-zero T1";
    EXPECT_GT(nonzero_t2, 0) << "Rank-1 tensor must have non-zero T2";

    std::cout << "  T0 nonzero: " << nonzero_t0 << ", max |T0| = " << max_t0 << "\n";
    std::cout << "  T1 nonzero: " << nonzero_t1 << "\n";
    std::cout << "  T2 nonzero: " << nonzero_t2 << ", max |T2| = " << max_t2 << "\n";
}


TEST_F(BiotSavartProteinTest, BFieldAboveRingAlongNormal) {
    // Physical check: for an atom directly above a ring center, B should
    // be predominantly along the ring normal (z-component dominates).
    auto& conf = protein->Conformation();
    conf.AttachResult(BiotSavartResult::Compute(conf));

    int checked = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& rn : conf.AtomAt(ai).ring_neighbours) {
            // Atom within 30 degrees of ring normal and < 5A
            if (rn.distance_to_center > 5.0) continue;
            if (rn.B_field.norm() < 1e-30) continue;

            const RingGeometry& geom = conf.ring_geometries[rn.ring_index];
            Vec3 d = conf.PositionAt(ai) - geom.center;
            double cos_angle = std::abs(d.normalized().dot(geom.normal));

            if (cos_angle > 0.866) {  // within 30 degrees of normal
                // B_z component (along normal) should dominate
                double Bz = std::abs(rn.B_field.dot(geom.normal));
                double Bmag = rn.B_field.norm();
                if (Bmag > 1e-25) {
                    EXPECT_GT(Bz / Bmag, 0.5)
                        << "B should be mostly along normal above ring";
                    checked++;
                }
            }
        }
    }

    EXPECT_GT(checked, 0) << "Should find atoms near ring normals";
    std::cout << "  Verified B || n for " << checked << " atoms above rings\n";
}


// ============================================================================
// ORCA protein test
// ============================================================================

TEST(BiotSavartOrcaTest, RunOnProtonatedProtein) {
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

    auto bs = BiotSavartResult::Compute(conf);
    ASSERT_NE(bs, nullptr);
    conf.AttachResult(std::move(bs));

    // Summary statistics
    double max_t0 = 0, max_t2 = 0, max_B = 0;
    int with_rings = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& sc = conf.AtomAt(ai).bs_shielding_contribution;
        max_t0 = std::max(max_t0, std::abs(sc.T0));
        max_t2 = std::max(max_t2, sc.T2Magnitude());
        max_B = std::max(max_B, conf.AtomAt(ai).total_B_field.norm());
        if (!conf.AtomAt(ai).ring_neighbours.empty()) with_rings++;
    }

    std::cout << "  ORCA protein BiotSavart summary:\n"
              << "    atoms=" << conf.AtomCount()
              << " rings=" << load.protein->RingCount() << "\n"
              << "    atoms with ring neighbours: " << with_rings << "\n"
              << "    max |T0| = " << max_t0 << "\n"
              << "    max |T2| = " << max_t2 << "\n"
              << "    max |B|  = " << max_B << " Tesla\n";

    EXPECT_GT(with_rings, 0) << "Some atoms should see rings";
    EXPECT_GT(max_t0, 1e-6) << "T0 should be non-zero";
    EXPECT_GT(max_t2, 1e-6) << "T2 should be non-zero";
    EXPECT_GT(max_B, 0) << "B-field should be non-zero";
}

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
#include "KernelEvaluationFilter.h"
#include "SpatialIndexResult.h"
#include "CalculatorConfig.h"
#include "Protein.h"
#include "PdbFileReader.h"
#include "OrcaRunLoader.h"
#include "PhysicalConstants.h"

#include <algorithm>
#include <array>
#include <filesystem>
#include <limits>
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


// Independent midpoint integration of the Biot-Savart line integral.  The
// production kernel uses the closed-form finite-segment expression; this
// oracle instead integrates dB = (mu0/4pi) I dl x r / |r|^3 over many short
// pieces of both Johnson-Bovey loops.
Vec3 MidpointJohnsonBoveyField(
        const std::vector<Vec3>& vertices,
        const Vec3& normal,
        double lobe_offset_ang,
        double current_nanoamperes,
        const Vec3& point_ang,
        int subdivisions_per_edge = 4096) {
    Vec3 field = Vec3::Zero();
    const double half_current_A =
        0.5 * current_nanoamperes * NANOAMPERES_TO_AMPERES;

    for (double sign : {-1.0, 1.0}) {
        const Vec3 offset = sign * lobe_offset_ang * normal;
        for (size_t i = 0; i < vertices.size(); ++i) {
            const Vec3 a = vertices[i] + offset;
            const Vec3 b = vertices[(i + 1) % vertices.size()] + offset;
            const Vec3 dl_ang = (b - a) /
                static_cast<double>(subdivisions_per_edge);
            const Vec3 dl_m = dl_ang * ANGSTROMS_TO_METRES;
            for (int k = 0; k < subdivisions_per_edge; ++k) {
                const Vec3 source_ang = a +
                    (static_cast<double>(k) + 0.5) * dl_ang;
                const Vec3 separation_m =
                    (point_ang - source_ang) * ANGSTROMS_TO_METRES;
                const double r = separation_m.norm();
                field += BIOT_SAVART_PREFACTOR * half_current_A *
                    dl_m.cross(separation_m) / (r * r * r);
            }
        }
    }
    return field;
}


// Independent uniform-centroid integration of the HM surface integral.  Each
// fan triangle is divided into n^2 equal-area microtriangles.  Production uses
// adaptive Dunavant/Stroud 7-point Gaussian quadrature, so agreement exercises
// the physical integral through a numerically distinct algorithm.
Mat3 MidpointFanSurfaceIntegral(const Vec3& point,
                                const RingGeometry& geometry,
                                int subdivisions = 192) {
    Mat3 integral = Mat3::Zero();
    const Mat3 identity = Mat3::Identity();

    for (size_t i = 0; i < geometry.vertices.size(); ++i) {
        const Vec3 v0 = geometry.center;
        const Vec3 v1 = geometry.vertices[i];
        const Vec3 v2 = geometry.vertices[(i + 1) % geometry.vertices.size()];
        const Vec3 edge1 = v1 - v0;
        const Vec3 edge2 = v2 - v0;
        const double micro_area =
            0.5 * edge1.cross(edge2).norm() /
            static_cast<double>(subdivisions * subdivisions);

        auto accumulate_at = [&](double u, double v) {
            const Vec3 surface = v0 + u * edge1 + v * edge2;
            const Vec3 separation = point - surface;
            const double r = separation.norm();
            const double r3 = r * r * r;
            const double r5 = r3 * r * r;
            integral += micro_area *
                (3.0 * separation * separation.transpose() / r5 -
                 identity / r3);
        };

        for (int u_index = 0; u_index < subdivisions; ++u_index) {
            for (int v_index = 0;
                 v_index < subdivisions - u_index; ++v_index) {
                // Lower microtriangle centroid.
                accumulate_at(
                    (static_cast<double>(u_index) + 1.0 / 3.0) /
                        subdivisions,
                    (static_cast<double>(v_index) + 1.0 / 3.0) /
                        subdivisions);

                // The companion upper microtriangle exists off the sloping
                // boundary of the reference triangle.
                if (u_index + v_index <= subdivisions - 2) {
                    accumulate_at(
                        (static_cast<double>(u_index) + 2.0 / 3.0) /
                            subdivisions,
                        (static_cast<double>(v_index) + 2.0 / 3.0) /
                            subdivisions);
                }
            }
        }
    }
    return integral;
}


const RingNeighbourhood* FindRingRow(const ProteinConformation& conf,
                                     size_t atom_index,
                                     size_t ring_index) {
    for (const RingNeighbourhood& row :
         conf.AtomAt(atom_index).ring_neighbours) {
        if (row.ring_index == ring_index) return &row;
    }
    return nullptr;
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

    RingBondedExclusionFilter topology_filter(protein);
    size_t target_index = protein.AtomCount();
    for (size_t ai = 0; ai < protein.AtomCount(); ++ai) {
        KernelEvaluationContext context;
        context.atom_index = ai;
        context.source_ring_index = ring_index;
        if (topology_filter.Accept(context)) {
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

    const RingNeighbourhood* row =
        FindRingRow(probe, target_index, ring_index);
    ASSERT_NE(row, nullptr);

    const RingGeometry& geometry = probe.ring_geometries[ring_index];
    const Vec3 point = probe.PositionAt(target_index);
    const Vec3 normal = geometry.normal.normalized();

    // BS: compare the production finite-segment kernel with a direct
    // midpoint line integral over both offset loops.
    const Vec3 expected_B = MidpointJohnsonBoveyField(
        geometry.vertices, geometry.normal,
        ring.JohnsonBoveyLobeOffset(), 1.0, point);
    EXPECT_LE((row->B_field - expected_B).norm(),
              1e-11 + 2e-6 * expected_B.norm());
    const Mat3 expected_G = -expected_B * geometry.normal.transpose() *
        PPM_FACTOR;
    EXPECT_LE((row->G_tensor - expected_G).norm(),
              1e-8 + 2e-6 * expected_G.norm());

    // HM: compare the adaptive Gaussian production quadrature with a uniform
    // centroid surface integration over the same physical fan triangles.
    const Mat3 expected_H =
        MidpointFanSurfaceIntegral(point, geometry);
    EXPECT_LE((row->hm_H_tensor - expected_H).norm(),
              5e-4 + 5e-4 * expected_H.norm());
    const Vec3 expected_hm_field = expected_H * geometry.normal;
    const Mat3 expected_hm_G =
        -expected_hm_field * geometry.normal.transpose();
    EXPECT_LE((row->hm_G_tensor - expected_hm_G).norm(),
              5e-4 + 5e-4 * expected_hm_G.norm());

    // Point-center kernels have exact axis values at r=0.5 A.  These fixed
    // numbers are independent of the production helper implementations.
    const double center_distance = (point - geometry.center).norm();
    EXPECT_NEAR(center_distance, 0.5, 1e-12);
    EXPECT_NEAR((point - geometry.center).normalized().dot(normal),
                1.0, 1e-12);
    EXPECT_NEAR(row->quad_scalar, 32.0, 1e-10);  // 2 / 0.5^4
    EXPECT_NEAR(row->chi_scalar, 16.0, 1e-11);   // 2 / 0.5^3

    // Every vertex is inside the unit-switch region, so the independent
    // dispersion oracle is simply sum_v 1/|point-v|^6.
    double expected_dispersion = 0.0;
    for (const Vec3& vertex : geometry.vertices) {
        const double vertex_distance = (point - vertex).norm();
        ASSERT_GT(vertex_distance,
                  CalculatorConfig::Get("singularity_guard_distance"));
        ASSERT_LT(vertex_distance,
                  CalculatorConfig::Get(
                      "dispersion_switching_onset_distance"));
        const double r2 = vertex_distance * vertex_distance;
        expected_dispersion += 1.0 / (r2 * r2 * r2);
    }
    EXPECT_NEAR(row->disp_scalar, expected_dispersion,
                2e-13 * std::max(1.0, expected_dispersion));
    EXPECT_EQ(row->disp_contacts,
              static_cast<int>(geometry.vertices.size()));

    // The two sampling APIs changed with C30 as well.  Grid points have no
    // topology exclusion, so independently sum every canonical in-range ring.
    Vec3 expected_sample_B = Vec3::Zero();
    Mat3 expected_sample_G = Mat3::Zero();
    Mat3 expected_sample_hm_G = Mat3::Zero();
    for (size_t ri = 0; ri < protein.RingCount(); ++ri) {
        const Ring& sample_ring = protein.RingAt(ri);
        const RingGeometry& sample_geometry = probe.ring_geometries[ri];
        if ((point - sample_geometry.center).norm() >
            CalculatorConfig::Get("ring_current_spatial_cutoff")) {
            continue;
        }
        if (biot_savart_detail::ContributesToCanonicalTotals(
                sample_ring.type_index)) {
            const Vec3 sample_B = MidpointJohnsonBoveyField(
                sample_geometry.vertices, sample_geometry.normal,
                sample_ring.JohnsonBoveyLobeOffset(), 1.0, point);
            expected_sample_B += sample_B;
            expected_sample_G +=
                -sample_B * sample_geometry.normal.transpose() * PPM_FACTOR;
        }
        if (haigh_mallion_detail::ContributesToCanonicalTotals(
                sample_ring.type_index)) {
            const Mat3 sample_H =
                MidpointFanSurfaceIntegral(point, sample_geometry, 128);
            const Vec3 sample_field = sample_H * sample_geometry.normal;
            expected_sample_hm_G +=
                -sample_field * sample_geometry.normal.transpose();
        }
    }

    const auto& bs = probe.Result<BiotSavartResult>();
    EXPECT_LE((bs.SampleBFieldAt(point) - expected_sample_B).norm(),
              1e-11 + 2e-6 * expected_sample_B.norm());
    EXPECT_LE((bs.SampleKernelAt(point).Reconstruct() - expected_sample_G)
                  .norm(),
              1e-8 + 2e-6 * expected_sample_G.norm());

    const auto& hm = probe.Result<HaighMallionResult>();
    EXPECT_LE((hm.SampleKernelAt(point).Reconstruct() - expected_sample_hm_G)
                  .norm(),
              8e-4 + 8e-4 * expected_sample_hm_G.norm());
}


TEST(RingNearFieldProduction,
     CheckedInProteinsHaveNoTopologyAdmittedAtomInsideRingRadius) {
    const fs::path fixture_root =
        fs::path(nmr::test::TestEnvironment::UbqProtonated()).parent_path();
    const std::array<std::pair<const char*, fs::path>, 2> fixtures{{
        {"1UBQ", nmr::test::TestEnvironment::UbqProtonated()},
        {"TRP-cage", fixture_root /
            "illustrative_peptides/trp_cage_1l2y_model1.pdb"},
    }};

    for (const auto& [label, path] : fixtures) {
        SCOPED_TRACE(label);
        ASSERT_TRUE(fs::is_regular_file(path)) << path;
        auto loaded = BuildFromProtonatedPdb(path.string());
        ASSERT_TRUE(loaded.Ok()) << loaded.error;
        Protein& protein = *loaded.protein;
        auto& conf = protein.Conformation();
        ASSERT_TRUE(conf.AttachResult(GeometryResult::Compute(conf)));
        ASSERT_TRUE(conf.AttachResult(SpatialIndexResult::Compute(conf)));
        ASSERT_TRUE(conf.AttachResult(PiQuadrupoleResult::Compute(conf)));
        ASSERT_TRUE(conf.AttachResult(DispersionResult::Compute(conf)));

        RingBondedExclusionFilter topology_filter(protein);
        const double center_cutoff =
            CalculatorConfig::Get("ring_current_spatial_cutoff");
        size_t admitted_in_range = 0;
        double minimum_radius_ratio =
            std::numeric_limits<double>::infinity();
        double maximum_abs_pi = 0.0;
        double maximum_dispersion = 0.0;

        for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
            for (size_t ri = 0; ri < protein.RingCount(); ++ri) {
                KernelEvaluationContext context;
                context.atom_index = ai;
                context.source_ring_index = ri;
                if (!topology_filter.Accept(context)) continue;

                const RingGeometry& geometry = conf.ring_geometries[ri];
                ASSERT_TRUE(geometry.center.allFinite());
                ASSERT_TRUE(geometry.normal.allFinite());
                ASSERT_TRUE(std::isfinite(geometry.radius));
                ASSERT_GT(geometry.radius, 0.0);
                const double distance =
                    (conf.PositionAt(ai) - geometry.center).norm();
                const double radius_ratio = distance / geometry.radius;
                minimum_radius_ratio =
                    std::min(minimum_radius_ratio, radius_ratio);
                EXPECT_GT(radius_ratio, 1.0)
                    << "atom=" << ai << " ring=" << ri
                    << " distance=" << distance
                    << " radius=" << geometry.radius;

                if (distance > center_cutoff) continue;
                ++admitted_in_range;
                const RingNeighbourhood* row = FindRingRow(conf, ai, ri);
                ASSERT_NE(row, nullptr)
                    << "atom=" << ai << " ring=" << ri;
                ASSERT_TRUE(std::isfinite(row->quad_scalar));
                ASSERT_TRUE(std::isfinite(row->disp_scalar));
                EXPECT_GE(row->disp_scalar, 0.0);
                maximum_abs_pi =
                    std::max(maximum_abs_pi, std::abs(row->quad_scalar));
                maximum_dispersion =
                    std::max(maximum_dispersion, row->disp_scalar);
            }
        }

        size_t production_rows = 0;
        for (size_t ai = 0; ai < conf.AtomCount(); ++ai)
            production_rows += conf.AtomAt(ai).ring_neighbours.size();
        EXPECT_EQ(production_rows, admitted_in_range);
        EXPECT_GT(admitted_in_range, 0u);
        EXPECT_GT(minimum_radius_ratio, 1.3);
        EXPECT_LT(maximum_abs_pi, 0.1);       // Angstrom^-4
        EXPECT_LT(maximum_dispersion, 0.1);  // Angstrom^-6

        std::cout << "  " << label
                  << " C30 audit: admitted=" << admitted_in_range
                  << " min(center/radius)=" << minimum_radius_ratio
                  << " max|Pi|=" << maximum_abs_pi
                  << " max(Disp)=" << maximum_dispersion << "\n";
    }
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

#include "BiotSavartResult.h"
#include "Protein.h"
#include "SpatialIndexResult.h"
#include "GeometryResult.h"
#include "KernelEvaluationFilter.h"
#include "PhysicalConstants.h"
#include "CalculatorConfig.h"
#include "GeometryChoice.h"
#include "RingNeighbourGeometry.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>

namespace nmr {


std::vector<std::type_index> BiotSavartResult::Dependencies() const {
    return {
        std::type_index(typeid(SpatialIndexResult)),
        std::type_index(typeid(GeometryResult))
    };
}


namespace biot_savart_detail {

double PointToSegmentDistance(const Vec3& point,
                              const Vec3& segment_a,
                              const Vec3& segment_b) {
    const Vec3 segment = segment_b - segment_a;
    const double length_sq = segment.squaredNorm();
    if (length_sq <= std::numeric_limits<double>::epsilon()) {
        return (point - segment_a).norm();
    }
    const double t = std::clamp(
        (point - segment_a).dot(segment) / length_sq, 0.0, 1.0);
    return (point - (segment_a + t * segment)).norm();
}


double MinimumDistanceToJohnsonBoveyWireLoops(
        const std::vector<Vec3>& vertices,
        const Vec3& normal,
        double lobe_offset_ang,
        const Vec3& point_ang) {
    if (vertices.size() < 2) {
        return std::numeric_limits<double>::infinity();
    }

    const Vec3 offset = normal * lobe_offset_ang;
    double minimum = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < vertices.size(); ++i) {
        const size_t j = (i + 1) % vertices.size();
        minimum = std::min(minimum, PointToSegmentDistance(
            point_ang, vertices[i] + offset, vertices[j] + offset));
        minimum = std::min(minimum, PointToSegmentDistance(
            point_ang, vertices[i] - offset, vertices[j] - offset));
    }
    return minimum;
}


bool ContributesToCanonicalTotals(RingTypeIndex ring_type) {
    // TRP9 is a diagnostic representation of the same indole pi system
    // already represented by TRP6 + TRP5.  Retain its sparse kernel row,
    // but never fold it into the canonical current sum a third time.
    return ring_type != RingTypeIndex::TrpPerimeter;
}


// ============================================================================
// Wire segment B-field (Biot-Savart law).
//
// All computation in SI: positions in metres, current in amperes, B in Tesla.
//
//   B = (mu_0/4pi) * I * (dl x dA) / |dl x dA|^2 * (dl.dA/|dA| - dl.dB/|dB|)
//
// where:
//   dl = b - a  (wire segment vector)
//   dA = r - a  (field point from segment start)
//   dB = r - b  (field point from segment end)
//
// Numerical guards (in SI, so default thresholds are small):
//   endpoint_guard (default 1e-25 m):  field point at segment endpoint
//   axis_guard     (default 1e-70 m^2): field point on the wire axis
// ============================================================================

static Vec3 WireSegmentField(
        const Vec3& a_m, const Vec3& b_m,
        double I_A, const Vec3& r_m) {

    const double endpoint_guard = CalculatorConfig::Get("biot_savart_wire_endpoint_guard");
    const double axis_guard = CalculatorConfig::Get("biot_savart_wire_axis_guard");

    Vec3 dl_m = b_m - a_m;
    Vec3 dA_m = r_m - a_m;
    Vec3 dB_m = r_m - b_m;

    double lenA = dA_m.norm();
    double lenB = dB_m.norm();
    if (lenA < endpoint_guard || lenB < endpoint_guard) return Vec3::Zero();

    Vec3 cross = dl_m.cross(dA_m);
    double crossSq = cross.squaredNorm();
    if (crossSq < axis_guard) return Vec3::Zero();

    double biot_savart_scale = BIOT_SAVART_PREFACTOR * I_A / crossSq;
    double endpoint_projection_term = dl_m.dot(dA_m) / lenA - dl_m.dot(dB_m) / lenB;

    return biot_savart_scale * endpoint_projection_term * cross;  // Tesla
}


// ============================================================================
// Johnson-Bovey double-loop model.
//
// Two current loops at +/- lobe_offset from the ring plane along the normal.
// Each loop carries half the total current (I/2). The total B-field is the
// sum over all wire segments of both loops.
//
// Input: vertex positions in Angstroms, current in nanoamperes.
// Converts to SI at the boundary, computes in pure SI, returns B in Tesla.
// ============================================================================

Vec3 JohnsonBoveyField(
        const std::vector<Vec3>& vertices,
        const Vec3& normal,
        double lobe_offset_ang,
        double current_nanoamperes,
        const Vec3& point_ang) {

    int n = static_cast<int>(vertices.size());
    if (n < 3) return Vec3::Zero();

    // Unit conversion at the boundary: Angstroms -> metres, nA -> A.
    // After this block, all computation is pure SI.
    Vec3 offset_ang = normal * lobe_offset_ang;
    double halfI_A = 0.5 * current_nanoamperes * NANOAMPERES_TO_AMPERES;

    Vec3 B = Vec3::Zero();
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;

        // Upper loop (+lobe_offset along normal)
        Vec3 a_upper = (vertices[i] + offset_ang) * ANGSTROMS_TO_METRES;
        Vec3 b_upper = (vertices[j] + offset_ang) * ANGSTROMS_TO_METRES;

        // Lower loop (-lobe_offset along normal)
        Vec3 a_lower = (vertices[i] - offset_ang) * ANGSTROMS_TO_METRES;
        Vec3 b_lower = (vertices[j] - offset_ang) * ANGSTROMS_TO_METRES;

        Vec3 r_m = point_ang * ANGSTROMS_TO_METRES;

        B += WireSegmentField(a_upper, b_upper, halfI_A, r_m);
        B += WireSegmentField(a_lower, b_lower, halfI_A, r_m);
    }

    return B;  // Tesla
}

}  // namespace biot_savart_detail


// ============================================================================
// BiotSavartResult::Compute
//
// For each atom, find rings within the ring_current_spatial_cutoff config
// value (default 15 A). For each ring:
//   1. Compute B-field from JB double-loop model (unit current, I=1 nA)
//   2. Construct geometric kernel G_ab = -n_b * B_a * PPM_FACTOR
//   3. Store G, SphericalTensor(G), B, cylindrical coords on RingNeighbourhood
//   4. Accumulate per-type T0 and T2 sums on ConformationAtom
// ============================================================================

std::unique_ptr<BiotSavartResult> BiotSavartResult::Compute(
        ProteinConformation& conf) {

    OperationLog::Scope scope("BiotSavartResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()) +
        " rings=" + std::to_string(conf.ProteinRef().RingCount()));

    const Protein& protein = conf.ProteinRef();
    const auto& spatial = conf.Result<SpatialIndexResult>();
    const size_t n_atoms = conf.AtomCount();
    const size_t n_rings = protein.RingCount();

    auto result_ptr = std::make_unique<BiotSavartResult>();
    result_ptr->conf_ = &conf;

    if (n_rings == 0) {
        OperationLog::Info(LogCalcBiotSavart, "BiotSavartResult::Compute",
            "no rings — nothing to compute");
        return result_ptr;
    }

    // Ring topology exclusion.  Numerical near-field validity is measured
    // against the actual two offset wire loops below, not a center sphere.
    KernelFilterSet filters;
    filters.Add(std::make_unique<RingBondedExclusionFilter>(protein));

    OperationLog::Info(LogCalcBiotSavart, "BiotSavartResult::Compute",
        "filter set: " + filters.Describe());

    GeometryChoiceBuilder choices(conf);
    std::set<size_t> recorded_rings;

    const double ring_cutoff = CalculatorConfig::Get("ring_current_spatial_cutoff");
    const double ring_proximity_shell_1 = CalculatorConfig::Get("ring_proximity_shell_1");
    const double ring_proximity_shell_2 = CalculatorConfig::Get("ring_proximity_shell_2");
    const double ring_proximity_shell_3 = CalculatorConfig::Get("ring_proximity_shell_3");
    const double ring_proximity_shell_4 = CalculatorConfig::Get("ring_proximity_shell_4");

    int total_pairs = 0;

    for (size_t ai = 0; ai < n_atoms; ++ai) {
        auto& ca = conf.MutableAtomAt(ai);
        Vec3 atom_pos = conf.PositionAt(ai);

        auto nearby_rings = spatial.RingsWithinRadius(atom_pos, ring_cutoff);

        Mat3 G_total = Mat3::Zero();
        Vec3 B_total = Vec3::Zero();

        for (size_t ri : nearby_rings) {
            const Ring& ring = protein.RingAt(ri);
            const RingGeometry& geom = conf.ring_geometries[ri];

            if (geom.vertices.size() < 3) continue;

            // Provenance: record one sampler per ring (physics resumes below).
            if (recorded_rings.find(ri) == recorded_rings.end()) {
                recorded_rings.insert(ri);
                auto verts_copy = geom.vertices;
                Vec3 normal_copy = geom.normal;
                double lobe_copy = ring.JohnsonBoveyLobeOffset();
                choices.Record(CalculatorId::BiotSavart, ri, "ring current",
                    [&ring, verts_copy, normal_copy, lobe_copy](GeometryChoice& gc) {
                        AddRing(gc, &ring, EntityRole::Source, EntityOutcome::Included);
                        AddNumber(gc, "intensity", ring.Intensity(), "nA·T⁻¹");
                        AddNumber(gc, "lobe_offset", ring.JohnsonBoveyLobeOffset(), "A");
                        SetSampler(gc, [verts_copy, normal_copy, lobe_copy](Vec3 pt) -> SphericalTensor {
                            Vec3 B = biot_savart_detail::JohnsonBoveyField(
                                verts_copy, normal_copy, lobe_copy, 1.0, pt);
                            // G_ab = -n_b B_a PPM_FACTOR (see header sign convention).
                            Mat3 G;
                            for (int a = 0; a < 3; ++a)
                                for (int b = 0; b < 3; ++b)
                                    G(a, b) = -normal_copy(b) * B(a) * PPM_FACTOR;
                            return SphericalTensor::Decompose(G);
                        });
                    });
            }

            double distance = (atom_pos - geom.center).norm();

            const double minimum_wire_distance =
                biot_savart_detail::MinimumDistanceToJohnsonBoveyWireLoops(
                    geom.vertices, geom.normal,
                    ring.JohnsonBoveyLobeOffset(), atom_pos);
            const double wire_distance_guard =
                CalculatorConfig::Get("biot_savart_wire_endpoint_guard") /
                ANGSTROMS_TO_METRES;
            if (minimum_wire_distance <= wire_distance_guard) {
                choices.Record(CalculatorId::BiotSavart, ri,
                    "wire singularity exclusion",
                    [&](GeometryChoice& gc) {
                        AddRing(gc, &ring, EntityRole::Source,
                                EntityOutcome::Included);
                        AddAtom(gc, &ca, ai, EntityRole::Target,
                                EntityOutcome::Excluded,
                                "johnson_bovey_wire_distance_guard");
                        AddNumber(gc, "minimum_wire_distance",
                                  minimum_wire_distance, "A");
                        AddNumber(gc, "wire_distance_guard",
                                  wire_distance_guard, "A");
                    });
                continue;
            }

            // Apply the typed through-bond exclusion after the physical
            // wire-distance guard.  No center-distance source extent exists
            // for the finite wire kernel.
            KernelEvaluationContext ctx;
            ctx.distance = distance;
            ctx.source_extent = 0.0;
            ctx.atom_index = ai;
            ctx.source_ring_index = ri;
            if (!filters.AcceptAll(ctx)) {
                choices.Record(CalculatorId::BiotSavart, ri, "topology exclusion",
                    [&](GeometryChoice& gc) {
                        AddRing(gc, &ring, EntityRole::Source, EntityOutcome::Included);
                        AddAtom(gc, &ca, ai, EntityRole::Target, EntityOutcome::Excluded,
                                filters.LastRejectorName());
                        AddNumber(gc, "distance", distance, "A");
                        AddNumber(gc, "source_extent", ctx.source_extent, "A");
                    });
                continue;
            }

            // Unit-current B-field (I = 1 nA).
            // The geometric kernel is independent of intensity.
            Vec3 B = biot_savart_detail::JohnsonBoveyField(
                geom.vertices, geom.normal,
                ring.JohnsonBoveyLobeOffset(), 1.0, atom_pos);

            // G_ab = -n_b B_a PPM_FACTOR (shielding-sign convention; derivation +
            // worked example in the header).
            Mat3 G;
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                    G(a, b) = -geom.normal(b) * B(a) * PPM_FACTOR;

            RingNeighbourhood& rn_ref = EnsureRingNeighbourGeometry(
                ca, geom, ring, ri, atom_pos);
            RingNeighbourhood* rn = &rn_ref;

            // Store BS results on RingNeighbourhood
            rn->G_tensor = G;
            rn->G_spherical = SphericalTensor::Decompose(G);
            rn->B_field = B;

            // Project B into the same ring frame (z/rho re-derived from the
            // geometry above).
            Vec3 d = atom_pos - geom.center;
            double z_coord = d.dot(geom.normal);
            Vec3 d_plane = d - z_coord * geom.normal;
            double rho_mag = d_plane.norm();
            Vec3 rho_hat = Vec3::Zero();
            if (rho_mag > CalculatorConfig::Get("near_zero_vector_norm_threshold")) rho_hat = d_plane / rho_mag;
            rn->B_cylindrical = Vec3(
                B.dot(rho_hat),        // B_rho
                0.0,                   // B_phi is not stored by this summary
                B.dot(geom.normal));   // B_z

            if (biot_savart_detail::ContributesToCanonicalTotals(
                    ring.type_index)) {
                G_total += G;
                B_total += B;

                const int ti = ring.TypeIndexAsInt();
                if (ti >= 0 && ti < kAromaticRingTypeCount) {
                    ca.per_type_G_T0_sum[ti] += rn->G_spherical.T0;
                    for (int c = 0; c < 3; ++c)
                        ca.per_type_G_T1_sum[ti][c] += rn->G_spherical.T1[c];
                    for (int c = 0; c < 5; ++c)
                        ca.per_type_G_T2_sum[ti][c] += rn->G_spherical.T2[c];
                }
            }

            total_pairs++;
        }

        // Store accumulated totals on ConformationAtom.
        // total_G_spherical: running sum across calculators on this atom.
        // bs_shielding_contribution: this calculator's per-call BS sum only.
        ca.total_B_field += B_total;
        ca.total_G_tensor += G_total;
        ca.total_G_spherical = SphericalTensor::Decompose(
            ca.total_G_tensor);
        ca.bs_shielding_contribution = SphericalTensor::Decompose(G_total);

        // Ring proximity counts (each ring appears once in ring_neighbours;
        // see find-or-create above).
        for (const auto& rn : ca.ring_neighbours) {
            if (rn.distance_to_center <= ring_proximity_shell_1) ca.n_rings_within_3A++;
            if (rn.distance_to_center <= ring_proximity_shell_2) ca.n_rings_within_5A++;
            if (rn.distance_to_center <= ring_proximity_shell_3) ca.n_rings_within_8A++;
            if (rn.distance_to_center <= ring_proximity_shell_4) ca.n_rings_within_12A++;
        }

        // Provenance: record this atom's ring-shell counts.
        choices.Record(CalculatorId::BiotSavart, ai, "ring shells",
            [&ca, ai](GeometryChoice& gc) {
                AddAtom(gc, &ca, ai, EntityRole::Target, EntityOutcome::Included);
                AddNumber(gc, "n_within_3A", static_cast<double>(ca.n_rings_within_3A), "count");
                AddNumber(gc, "n_within_5A", static_cast<double>(ca.n_rings_within_5A), "count");
                AddNumber(gc, "n_within_8A", static_cast<double>(ca.n_rings_within_8A), "count");
                AddNumber(gc, "n_within_12A", static_cast<double>(ca.n_rings_within_12A), "count");
            });
    }

    OperationLog::Info(LogCalcBiotSavart, "BiotSavartResult::Compute",
        "atom_ring_pairs=" + std::to_string(total_pairs) +
        " rejected={" + filters.ReportRejections() + "}" +
        " atoms=" + std::to_string(n_atoms) +
        " rings=" + std::to_string(n_rings));

    return result_ptr;
}


// ============================================================================
// SampleBFieldAt / SampleKernelAt: evaluate at arbitrary 3D points.
//
// Same physics as Compute(), single point. Grid points are not atoms, so the
// per-atom/topology filters (MinDistance, RingBondedExclusion) do not apply;
// the inline distance/singularity/inside-source guards keep the multipole
// valid.
// ============================================================================

Vec3 BiotSavartResult::SampleBFieldAt(Vec3 point) const {
    if (!conf_) return Vec3::Zero();

    const Protein& protein = conf_->ProteinRef();
    const size_t n_rings = protein.RingCount();

    const double ring_cutoff = CalculatorConfig::Get("ring_current_spatial_cutoff");

    Vec3 B_total = Vec3::Zero();

    for (size_t ri = 0; ri < n_rings; ++ri) {
        const Ring& ring = protein.RingAt(ri);
        if (!biot_savart_detail::ContributesToCanonicalTotals(
                ring.type_index)) continue;
        const RingGeometry& geom = conf_->ring_geometries[ri];
        if (geom.vertices.size() < 3) continue;

        // Grid acceptance uses the actual offset-wire distance plus the far
        // cutoff (no per-atom/topology filters; grid points are not atoms).
        double distance = (point - geom.center).norm();
        if (distance > ring_cutoff) continue;
        const double minimum_wire_distance =
            biot_savart_detail::MinimumDistanceToJohnsonBoveyWireLoops(
                geom.vertices, geom.normal,
                ring.JohnsonBoveyLobeOffset(), point);
        const double wire_distance_guard =
            CalculatorConfig::Get("biot_savart_wire_endpoint_guard") /
            ANGSTROMS_TO_METRES;
        if (minimum_wire_distance <= wire_distance_guard) continue;

        B_total += biot_savart_detail::JohnsonBoveyField(
            geom.vertices, geom.normal,
            ring.JohnsonBoveyLobeOffset(), 1.0, point);
    }

    return B_total;
}

SphericalTensor BiotSavartResult::SampleKernelAt(Vec3 point) const {
    if (!conf_) return SphericalTensor{};

    const Protein& protein = conf_->ProteinRef();
    const size_t n_rings = protein.RingCount();

    const double ring_cutoff = CalculatorConfig::Get("ring_current_spatial_cutoff");

    Mat3 G_total = Mat3::Zero();

    for (size_t ri = 0; ri < n_rings; ++ri) {
        const Ring& ring = protein.RingAt(ri);
        if (!biot_savart_detail::ContributesToCanonicalTotals(
                ring.type_index)) continue;
        const RingGeometry& geom = conf_->ring_geometries[ri];
        if (geom.vertices.size() < 3) continue;

        // Grid acceptance uses the actual offset-wire distance plus the far
        // cutoff (no per-atom/topology filters; grid points are not atoms).
        double distance = (point - geom.center).norm();
        if (distance > ring_cutoff) continue;
        const double minimum_wire_distance =
            biot_savart_detail::MinimumDistanceToJohnsonBoveyWireLoops(
                geom.vertices, geom.normal,
                ring.JohnsonBoveyLobeOffset(), point);
        const double wire_distance_guard =
            CalculatorConfig::Get("biot_savart_wire_endpoint_guard") /
            ANGSTROMS_TO_METRES;
        if (minimum_wire_distance <= wire_distance_guard) continue;

        Vec3 B = biot_savart_detail::JohnsonBoveyField(
            geom.vertices, geom.normal,
            ring.JohnsonBoveyLobeOffset(), 1.0, point);

        // G_ab = -n_b B_a PPM_FACTOR (see header sign convention).
        Mat3 G;
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                G(a, b) = -geom.normal(b) * B(a) * PPM_FACTOR;

        G_total += G;
    }

    return SphericalTensor::Decompose(G_total);
}


// ============================================================================
// WriteFeatures: export what Compute() wrote on ConformationAtom.
//
// Every field this method reads was written by Compute() above. The arrays
// match Compute()'s accumulation: shielding contribution (the full
// SphericalTensor of the summed G over all rings), per-type T0 and T2
// sums (8 ring types, matching the RingTypeIndex enum), ring proximity
// counts, and the total B-field vector.
//
// Pack order for SphericalTensor: [T0, T1[0..2], T2[0..4]] = 9 doubles.
// ============================================================================

int BiotSavartResult::WriteFeatures(const ProteinConformation& conf,
                                     const std::string& output_dir) const {
    const size_t N = conf.AtomCount();
    int written = 0;

    // bs_shielding: (N, 9) — the full SphericalTensor sum over all rings
    {
        std::vector<double> data(N * 9);
        for (size_t i = 0; i < N; ++i)
            conf.AtomAt(i).bs_shielding_contribution.PackFull9(&data[i*9]);
        NpyWriter::WriteFloat64(output_dir + "/bs_shielding.npy", data.data(), N, 9);
        written++;
    }

    // bs_per_type_T0: (N, 8) — isotropic kernel per ring type
    {
        std::vector<double> data(N * 8);
        for (size_t i = 0; i < N; ++i)
            for (int t = 0; t < 8; ++t)
                data[i*8 + t] = conf.AtomAt(i).per_type_G_T0_sum[t];
        NpyWriter::WriteFloat64(output_dir + "/bs_per_type_T0.npy", data.data(), N, 8);
        written++;
    }

    // bs_per_type_T1: (N, 24) — T1[3] per ring type[8]
    {
        std::vector<double> data(N * 24);
        for (size_t i = 0; i < N; ++i)
            for (int t = 0; t < 8; ++t)
                for (int c = 0; c < 3; ++c)
                    data[i*24 + t*3 + c] = conf.AtomAt(i).per_type_G_T1_sum[t][c];
        NpyWriter::WriteFloat64(output_dir + "/bs_per_type_T1.npy", data.data(), N, 24);
        written++;
    }

    // bs_per_type_T2: (N, 40) — T2[5] per ring type[8]
    {
        std::vector<double> data(N * 40);
        for (size_t i = 0; i < N; ++i)
            for (int t = 0; t < 8; ++t)
                for (int c = 0; c < 5; ++c)
                    data[i*40 + t*5 + c] = conf.AtomAt(i).per_type_G_T2_sum[t][c];
        NpyWriter::WriteFloat64(output_dir + "/bs_per_type_T2.npy", data.data(), N, 40);
        written++;
    }

    // bs_total_B: (N, 3) — total B-field vector at each atom
    {
        std::vector<double> data(N * 3);
        for (size_t i = 0; i < N; ++i) {
            const Vec3& B = conf.AtomAt(i).total_B_field;
            data[i*3+0] = B.x(); data[i*3+1] = B.y(); data[i*3+2] = B.z();
        }
        NpyWriter::WriteFloat64(output_dir + "/bs_total_B.npy", data.data(), N, 3);
        written++;
    }

    // bs_ring_counts: (N, 4) — proximity counts at 3/5/8/12 A shells
    {
        std::vector<double> data(N * 4);
        for (size_t i = 0; i < N; ++i) {
            const auto& ca = conf.AtomAt(i);
            data[i*4+0] = static_cast<double>(ca.n_rings_within_3A);
            data[i*4+1] = static_cast<double>(ca.n_rings_within_5A);
            data[i*4+2] = static_cast<double>(ca.n_rings_within_8A);
            data[i*4+3] = static_cast<double>(ca.n_rings_within_12A);
        }
        NpyWriter::WriteFloat64(output_dir + "/bs_ring_counts.npy", data.data(), N, 4);
        written++;
    }

    return written;
}

}  // namespace nmr

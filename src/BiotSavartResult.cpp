#include "BiotSavartResult.h"
#include "Protein.h"
#include "SpatialIndexResult.h"
#include "GeometryResult.h"
#include "KernelEvaluationFilter.h"
#include "PhysicalConstants.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <cmath>

namespace nmr {


std::vector<std::type_index> BiotSavartResult::Dependencies() const {
    return {
        std::type_index(typeid(SpatialIndexResult)),
        std::type_index(typeid(GeometryResult))
    };
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
// Numerical guards (in SI, so thresholds are small):
//   lenA, lenB < 1e-25 m:  field point at segment endpoint
//   crossSq < 1e-70 m^2:   field point on the wire axis
// ============================================================================

static Vec3 WireSegmentField(
        const Vec3& a_m, const Vec3& b_m,
        double I_A, const Vec3& r_m) {

    Vec3 dl_m = b_m - a_m;
    Vec3 dA_m = r_m - a_m;
    Vec3 dB_m = r_m - b_m;

    double lenA = dA_m.norm();
    double lenB = dB_m.norm();
    if (lenA < 1e-25 || lenB < 1e-25) return Vec3::Zero();

    Vec3 cross = dl_m.cross(dA_m);
    double crossSq = cross.squaredNorm();
    if (crossSq < 1e-70) return Vec3::Zero();

    // mu_0/(4*pi) = 1e-7 T*m/A (exact in SI)
    constexpr double mu0_4pi = 1e-7;
    double factor = mu0_4pi * I_A / crossSq;
    double dotTerm = dl_m.dot(dA_m) / lenA - dl_m.dot(dB_m) / lenB;

    return factor * dotTerm * cross;  // Tesla
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

static Vec3 JohnsonBoveyField(
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

        // Upper loop (z = +d)
        Vec3 a_upper = (vertices[i] + offset_ang) * ANGSTROMS_TO_METRES;
        Vec3 b_upper = (vertices[j] + offset_ang) * ANGSTROMS_TO_METRES;

        // Lower loop (z = -d)
        Vec3 a_lower = (vertices[i] - offset_ang) * ANGSTROMS_TO_METRES;
        Vec3 b_lower = (vertices[j] - offset_ang) * ANGSTROMS_TO_METRES;

        Vec3 r_m = point_ang * ANGSTROMS_TO_METRES;

        B += WireSegmentField(a_upper, b_upper, halfI_A, r_m);
        B += WireSegmentField(a_lower, b_lower, halfI_A, r_m);
    }

    return B;  // Tesla
}


// ============================================================================
// BiotSavartResult::Compute
//
// For each atom, find rings within RING_CALC_CUTOFF (15A). For each ring:
//   1. Compute B-field from JB double-loop model (unit current, I=1 nA)
//   2. Construct geometric kernel G_ab = n_b * B_a * PPM_FACTOR
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

    // Filter set: DipolarNearFieldFilter with source_extent = ring diameter,
    // plus RingBondedExclusionFilter for topological exclusion of ring
    // vertices and their bonded neighbours.
    KernelFilterSet filters;
    filters.Add(std::make_unique<DipolarNearFieldFilter>());
    filters.Add(std::make_unique<RingBondedExclusionFilter>(protein));

    OperationLog::Info(LogCalcBiotSavart, "BiotSavartResult::Compute",
        "filter set: " + filters.Describe());

    int total_pairs = 0;

    for (size_t ai = 0; ai < n_atoms; ++ai) {
        auto& ca = conf.MutableAtomAt(ai);
        Vec3 atom_pos = conf.PositionAt(ai);

        auto nearby_rings = spatial.RingsWithinRadius(atom_pos, RING_CALC_CUTOFF);

        Mat3 G_total = Mat3::Zero();
        Vec3 B_total = Vec3::Zero();

        for (size_t ri : nearby_rings) {
            const Ring& ring = protein.RingAt(ri);
            const RingGeometry& geom = conf.ring_geometries[ri];

            if (geom.vertices.size() < 3) continue;

            double distance = (atom_pos - geom.center).norm();
            if (distance < MIN_DISTANCE) continue;

            // Apply filter: source extent = ring diameter (2 * radius)
            KernelEvaluationContext ctx;
            ctx.distance = distance;
            ctx.source_extent = 2.0 * geom.radius;
            ctx.atom_index = ai;
            ctx.source_ring_index = ri;
            if (!filters.AcceptAll(ctx)) continue;

            // B-field from JB model with unit current (1.0 nA).
            // The geometric kernel is independent of intensity.
            Vec3 B = JohnsonBoveyField(
                geom.vertices, geom.normal,
                ring.JBLobeOffset(), 1.0, atom_pos);

            // Geometric kernel: G_ab = -n_b * B_a * PPM_FACTOR
            // The minus sign comes from the shielding tensor definition:
            //   sigma_ab = -dB_a^sec / dB_{0,b}
            // With this convention, sigma = I * G gives the correct sign
            // using literature I values (I < 0 for diamagnetic rings).
            // Verified: I=-12, G_T0=-0.116 at (0,0,3A) above PHE
            //   -> sigma = (-12)(-0.116) = +1.40 ppm (shielded). Correct.
            Mat3 G;
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                    G(a, b) = -geom.normal(b) * B(a) * PPM_FACTOR;

            // Find or create RingNeighbourhood for this ring
            RingNeighbourhood* rn = nullptr;
            for (auto& existing : ca.ring_neighbours) {
                if (existing.ring_index == ri) {
                    rn = &existing;
                    break;
                }
            }
            if (!rn) {
                RingNeighbourhood new_rn;
                new_rn.ring_index = ri;
                new_rn.ring_type = ring.type_index;
                new_rn.distance_to_center = distance;

                Vec3 d = atom_pos - geom.center;
                new_rn.direction_to_center = d.normalized();

                // Cylindrical coordinates in ring frame
                double z = d.dot(geom.normal);
                Vec3 d_plane = d - z * geom.normal;
                double rho = d_plane.norm();
                double theta = std::atan2(d_plane.norm(), std::abs(z));
                new_rn.z = z;
                new_rn.rho = rho;
                new_rn.theta = theta;

                ca.ring_neighbours.push_back(new_rn);
                rn = &ca.ring_neighbours.back();
            }

            // Store BS results on RingNeighbourhood
            rn->G_tensor = G;
            rn->G_spherical = SphericalTensor::Decompose(G);
            rn->B_field = B;

            // B-field in cylindrical coordinates (ring frame)
            Vec3 d = atom_pos - geom.center;
            double z_coord = d.dot(geom.normal);
            Vec3 d_plane = d - z_coord * geom.normal;
            double rho_mag = d_plane.norm();
            Vec3 rho_hat = Vec3::Zero();
            if (rho_mag > NEAR_ZERO_NORM) rho_hat = d_plane / rho_mag;
            rn->B_cylindrical = Vec3(
                B.dot(rho_hat),        // B_rho
                0.0,                   // B_phi (zero by axial symmetry choice)
                B.dot(geom.normal));   // B_z

            // Accumulate totals
            G_total += G;
            B_total += B;

            // Per-type T0 and T2 sums
            int ti = ring.TypeIndexAsInt();
            if (ti >= 0 && ti < 8) {
                ca.per_type_G_T0_sum[ti] += rn->G_spherical.T0;
                for (int c = 0; c < 5; ++c)
                    ca.per_type_G_T2_sum[ti][c] += rn->G_spherical.T2[c];
            }

            total_pairs++;
        }

        // Store accumulated totals on ConformationAtom
        ca.total_B_field += B_total;
        ca.total_G_tensor += G_total;
        ca.total_G_spherical = SphericalTensor::Decompose(
            ca.total_G_tensor);
        ca.bs_shielding_contribution = SphericalTensor::Decompose(G_total);

        // Ring proximity counts
        for (const auto& rn : ca.ring_neighbours) {
            if (rn.distance_to_center <= RING_COUNT_SHELL_1) ca.n_rings_within_3A++;
            if (rn.distance_to_center <= RING_COUNT_SHELL_2) ca.n_rings_within_5A++;
            if (rn.distance_to_center <= RING_COUNT_SHELL_3) ca.n_rings_within_8A++;
            if (rn.distance_to_center <= RING_COUNT_SHELL_4) ca.n_rings_within_12A++;
        }
    }

    OperationLog::Info(LogCalcBiotSavart, "BiotSavartResult::Compute",
        "atom_ring_pairs=" + std::to_string(total_pairs) +
        " rejected={" + filters.ReportRejections() + "}" +
        " atoms=" + std::to_string(n_atoms) +
        " rings=" + std::to_string(n_rings));

    return result_ptr;
}


// ============================================================================
// SampleBFieldAt / SampleShieldingAt: evaluate at arbitrary 3D points.
//
// Same physics as Compute(), but for a single point. No filter set
// (grid points are not atoms — no self-exclusion or bonded-exclusion
// applies). DipolarNearFieldFilter still applied for physical validity.
// ============================================================================

Vec3 BiotSavartResult::SampleBFieldAt(Vec3 point) const {
    if (!conf_) return Vec3::Zero();

    const Protein& protein = conf_->ProteinRef();
    const size_t n_rings = protein.RingCount();

    Vec3 B_total = Vec3::Zero();

    for (size_t ri = 0; ri < n_rings; ++ri) {
        const Ring& ring = protein.RingAt(ri);
        const RingGeometry& geom = conf_->ring_geometries[ri];
        if (geom.vertices.size() < 3) continue;

        double distance = (point - geom.center).norm();
        if (distance < MIN_DISTANCE) continue;

        // DipolarNearFieldFilter: multipole invalid inside source
        if (distance < geom.radius) continue;

        // Distance cutoff
        if (distance > RING_CALC_CUTOFF) continue;

        B_total += JohnsonBoveyField(
            geom.vertices, geom.normal,
            ring.JBLobeOffset(), 1.0, point);
    }

    return B_total;
}

SphericalTensor BiotSavartResult::SampleShieldingAt(Vec3 point) const {
    if (!conf_) return SphericalTensor{};

    const Protein& protein = conf_->ProteinRef();
    const size_t n_rings = protein.RingCount();

    Mat3 G_total = Mat3::Zero();

    for (size_t ri = 0; ri < n_rings; ++ri) {
        const Ring& ring = protein.RingAt(ri);
        const RingGeometry& geom = conf_->ring_geometries[ri];
        if (geom.vertices.size() < 3) continue;

        double distance = (point - geom.center).norm();
        if (distance < MIN_DISTANCE) continue;
        if (distance < geom.radius) continue;
        if (distance > RING_CALC_CUTOFF) continue;

        Vec3 B = JohnsonBoveyField(
            geom.vertices, geom.normal,
            ring.JBLobeOffset(), 1.0, point);

        // G_ab = -n_b * B_a * PPM_FACTOR
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

static void PackST(const SphericalTensor& st, double* out) {
    out[0] = st.T0;
    for (int i = 0; i < 3; ++i) out[1+i] = st.T1[i];
    for (int i = 0; i < 5; ++i) out[4+i] = st.T2[i];
}

int BiotSavartResult::WriteFeatures(const ProteinConformation& conf,
                                     const std::string& output_dir) const {
    const size_t N = conf.AtomCount();
    int written = 0;

    // bs_shielding: (N, 9) — the full SphericalTensor sum over all rings
    {
        std::vector<double> data(N * 9);
        for (size_t i = 0; i < N; ++i)
            PackST(conf.AtomAt(i).bs_shielding_contribution, &data[i*9]);
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

#include "DispersionResult.h"
#include "Protein.h"
#include "SpatialIndexResult.h"
#include "GeometryResult.h"
#include "KernelEvaluationFilter.h"
#include "PhysicalConstants.h"
#include "CalculatorConfig.h"
#include "GeometryChoice.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <cmath>
#include <set>

namespace nmr {


std::vector<std::type_index> DispersionResult::Dependencies() const {
    return {
        std::type_index(typeid(SpatialIndexResult)),
        std::type_index(typeid(GeometryResult))
    };
}


// ============================================================================
// Smooth switching function for dispersion cutoff (CHARMM functional form).
//
// S(r) = 1                                                for r <= R_switch
// S(r) = (Rc²-r²)²(Rc²+2r²-3Rs²) / (Rc²-Rs²)³         for R_switch < r < R_cut
// S(r) = 0                                                for r >= R_cut
//
// C¹ continuous at both boundaries. S(R_switch) = 1, S(R_cut) = 0,
// S'(R_switch) = 0, S'(R_cut) = 0.
//
// Physics reason for smooth taper: the 1/r^6 interaction does not
// physically stop at any distance. The switching function tapers it
// to zero over a finite range so that features vary smoothly with
// atomic position. Essential when the same kernel is evaluated on
// 100+ MD frames where atom positions fluctuate by ~0.5-1A.
//
// Reference: Brooks et al., J. Comput. Chem. 4, 187 (1983) — the
// CHARMM switching function for non-bonded interactions.
//
// R_switch = 4.3 A: onset of taper. Below this, full strength.
// R_cut = 5.0 A: zero beyond this. At R_cut, 1/r^6 = 6.4e-5 A^-6,
//   which is 0.03% of the value at 2A (the typical nearest non-bonded
//   contact). Truncation error from stopping here is < 0.1%.
// ============================================================================

// Dispersion range limits (Angstroms).
// R_CUT: at 5A, C6/r^6 = C6/15625. For the total sum over ~6 vertices
// within range, the contribution beyond 5A is < 0.1% of the total.
// This is a numerical precision choice: we are truncating a convergent
// sum, not asserting that dispersion stops at 5A.
//
// R_SWITCH: onset of the smooth taper. Set at 4.3A (0.7A before R_CUT)
// to give a gentle taper over the range where 1/r^6 is already small.
// The taper width (0.7A) is comparable to atomic position fluctuations
// in MD (~0.5A RMS), ensuring no atom ever jumps across the entire
// taper in one frame.
static double DispSwitchingFunction(double r) {
    if (r <= CalculatorConfig::Get("dispersion_switching_onset_distance")) return 1.0;
    if (r >= CalculatorConfig::Get("dispersion_vertex_distance_cutoff")) return 0.0;

    double rc2 = CalculatorConfig::Get("dispersion_vertex_distance_cutoff") * CalculatorConfig::Get("dispersion_vertex_distance_cutoff");
    double rs2 = CalculatorConfig::Get("dispersion_switching_onset_distance") * CalculatorConfig::Get("dispersion_switching_onset_distance");
    double r2 = r * r;
    double num = (rc2 - r2) * (rc2 - r2) * (rc2 + 2.0 * r2 - 3.0 * rs2);
    double den = (rc2 - rs2) * (rc2 - rs2) * (rc2 - rs2);
    return num / den;
}


// ============================================================================
// London dispersion kernel from one ring vertex at one atom.
//
// Per vertex, with unit C6 = 1:
//
//   K_ab = S(r) * (3 d_a d_b / r^8 - delta_ab / r^6)    (Angstrom^-6)
//   scalar = S(r) / r^6                                    (Angstrom^-6)
//
// where d = r_atom - r_vertex, r = |d|, S(r) is the switching function.
//
// The tensor is traceless per vertex:
//   Tr(K) = S(r) * (3|d|^2/r^8 - 3/r^6) = S(r) * 0 = 0.
//
// Vertex exclusion: atoms covalently bonded to a ring vertex are
// excluded because the through-space 1/r^6 kernel does not model
// through-bond electronic coupling. This is checked via the protein's
// bond connectivity, not a distance heuristic.
// ============================================================================

struct DispVertexResult {
    Mat3 K = Mat3::Zero();
    double scalar = 0.0;
    bool valid = false;
};


static DispVertexResult ComputeDispVertex(
        const Vec3& atom_pos,
        const Vec3& vertex_pos,
        double r) {

    DispVertexResult result;

    if (r < CalculatorConfig::Get("singularity_guard_distance") || r > CalculatorConfig::Get("dispersion_vertex_distance_cutoff")) return result;

    double S = DispSwitchingFunction(r);
    if (S < CalculatorConfig::Get("dispersion_switching_noise_floor")) return result;  // below switching threshold

    Vec3 d = atom_pos - vertex_pos;
    double r2 = r * r;
    double r6 = r2 * r2 * r2;
    double r8 = r6 * r2;

    result.scalar = S / r6;

    // K_ab = S(r) * (3 d_a d_b / r^8 - delta_ab / r^6)
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            result.K(a, b) = S * (3.0 * d(a) * d(b) / r8
                                - (a == b ? 1.0 : 0.0) / r6);

    result.valid = true;
    return result;
}


// ============================================================================
// Build the set of atoms bonded to any vertex of a ring.
// Used to exclude through-bond pairs from the through-space 1/r^6 kernel.
// ============================================================================

static std::set<size_t> BondedToVertices(
        const Ring& ring, const Protein& protein) {
    std::set<size_t> bonded;
    for (size_t vi : ring.atom_indices) {
        bonded.insert(vi);  // the vertex itself
        const auto& atom = protein.AtomAt(vi);
        for (size_t bi : atom.bond_indices) {
            const auto& bond = protein.BondAt(bi);
            bonded.insert(bond.atom_index_a);
            bonded.insert(bond.atom_index_b);
        }
    }
    return bonded;
}


// ============================================================================
// DispersionResult::Compute
//
// For each atom, find rings within RING_CALC_CUTOFF (15A coarse spatial
// query). Apply DipolarNearFieldFilter at ring level (same as all other
// ring calculators). For each passing ring, sum the dispersion kernel
// over vertices, excluding atoms bonded to any vertex (through-bond
// exclusion). Smooth switching function tapers the 1/r^6 contribution
// to zero at DISP_VERTEX_R_CUT.
// ============================================================================

std::unique_ptr<DispersionResult> DispersionResult::Compute(
        ProteinConformation& conf) {

    OperationLog::Scope scope("DispersionResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()) +
        " rings=" + std::to_string(conf.ProteinRef().RingCount()));

    const Protein& protein = conf.ProteinRef();
    const auto& spatial = conf.Result<SpatialIndexResult>();
    const size_t n_atoms = conf.AtomCount();
    const size_t n_rings = protein.RingCount();

    auto result_ptr = std::make_unique<DispersionResult>();
    result_ptr->conf_ = &conf;

    if (n_rings == 0) {
        OperationLog::Info(LogCalcOther, "DispersionResult::Compute",
            "no rings — nothing to compute");
        return result_ptr;
    }

    // Ring-level filter: DipolarNearFieldFilter with source_extent = ring
    // diameter. Same physics as all other ring calculators: the per-vertex
    // summation is a discrete approximation that breaks down when the
    // field point is inside the ring.
    KernelFilterSet filters;
    filters.Add(std::make_unique<DipolarNearFieldFilter>());

    OperationLog::Info(LogCalcOther, "DispersionResult::Compute",
        "filter set: " + filters.Describe() +
        " | vertex range: [MIN_DISTANCE=" + std::to_string(CalculatorConfig::Get("singularity_guard_distance")) +
        ", R_CUT=" + std::to_string(CalculatorConfig::Get("dispersion_vertex_distance_cutoff")) +
        "] A, switch onset=" + std::to_string(CalculatorConfig::Get("dispersion_switching_onset_distance")) + " A" +
        " | through-bond vertex exclusion: yes");

    GeometryChoiceBuilder choices(conf);
    std::set<size_t> recorded_rings;

    // Pre-build bonded-to-vertex sets for each ring (once, not per atom).
    std::vector<std::set<size_t>> ring_bonded(n_rings);
    for (size_t ri = 0; ri < n_rings; ++ri)
        ring_bonded[ri] = BondedToVertices(protein.RingAt(ri), protein);

    int total_pairs = 0;
    int total_contacts = 0;
    int bonded_exclusions = 0;

    for (size_t ai = 0; ai < n_atoms; ++ai) {
        auto& ca = conf.MutableAtomAt(ai);
        Vec3 atom_pos = conf.PositionAt(ai);

        auto nearby_rings = spatial.RingsWithinRadius(atom_pos, CalculatorConfig::Get("ring_current_spatial_cutoff"));

        Mat3 disp_total = Mat3::Zero();

        for (size_t ri : nearby_rings) {
            const Ring& ring = protein.RingAt(ri);
            const RingGeometry& geom = conf.ring_geometries[ri];

            if (geom.vertices.empty()) continue;

            double dist_to_center = (atom_pos - geom.center).norm();

            // Ring-level filter
            KernelEvaluationContext ctx;
            ctx.distance = dist_to_center;
            ctx.source_extent = 2.0 * geom.radius;  // ring diameter (A)
            ctx.atom_index = ai;
            if (!filters.AcceptAll(ctx)) {
                // ---- GeometryChoice: near-field exclusion ----
                choices.Record(CalculatorId::Dispersion, ri, "near-field exclusion",
                    [&](GeometryChoice& gc) {
                        AddRing(gc, &ring, EntityRole::Source, EntityOutcome::Included);
                        AddAtom(gc, &ca, ai, EntityRole::Target, EntityOutcome::Excluded,
                                filters.LastRejectorName());
                        AddNumber(gc, "distance", dist_to_center, "A");
                        AddNumber(gc, "source_extent", ctx.source_extent, "A");
                    });
                continue;
            }

            // Through-bond exclusion: skip this ring entirely if the
            // field atom is bonded to any vertex (it's part of or
            // immediately adjacent to the ring).
            if (ring_bonded[ri].count(ai)) {
                // ---- GeometryChoice: through-bond exclusion ----
                choices.Record(CalculatorId::Dispersion, ri, "through-bond exclusion",
                    [&](GeometryChoice& gc) {
                        AddRing(gc, &ring, EntityRole::Source, EntityOutcome::Included);
                        AddAtom(gc, &ca, ai, EntityRole::Target, EntityOutcome::Excluded,
                                "ring_bonded");
                        AddNumber(gc, "distance", dist_to_center, "A");
                    });
                bonded_exclusions++;
                continue;
            }

            // ---- GeometryChoice: dispersion taper ----
            if (recorded_rings.insert(ri).second) {
                choices.Record(CalculatorId::Dispersion, ri, "dispersion taper",
                    [&ring](GeometryChoice& gc) {
                        AddRing(gc, &ring, EntityRole::Source, EntityOutcome::Included);
                        AddNumber(gc, "switch_onset", CalculatorConfig::Get("dispersion_switching_onset_distance"), "A");
                        AddNumber(gc, "cutoff", CalculatorConfig::Get("dispersion_vertex_distance_cutoff"), "A");
                    });
            }

            // Sum dispersion kernel over ring vertices
            Mat3 K_ring = Mat3::Zero();
            double s_ring = 0.0;
            int contacts = 0;

            for (size_t vi = 0; vi < ring.atom_indices.size(); ++vi) {
                Vec3 vpos = geom.vertices[vi];
                double r = (atom_pos - vpos).norm();

                DispVertexResult vr = ComputeDispVertex(atom_pos, vpos, r);
                if (!vr.valid) {
                    // ---- GeometryChoice: switching function noise floor ----
                    // Only fires when r is in the taper range but S < 1e-15
                    if (r > CalculatorConfig::Get("dispersion_switching_onset_distance") && r < CalculatorConfig::Get("dispersion_vertex_distance_cutoff")) {
                        choices.Record(CalculatorId::Dispersion, ri, "switching noise floor",
                            [&ring, &ca, ai, r](GeometryChoice& gc) {
                                AddRing(gc, &ring, EntityRole::Source, EntityOutcome::Included);
                                AddAtom(gc, &ca, ai, EntityRole::Target, EntityOutcome::Excluded,
                                        "switching_noise_floor");
                                AddNumber(gc, "vertex_distance", r, "A");
                            });
                    }
                    continue;
                }

                K_ring += vr.K;
                s_ring += vr.scalar;
                contacts++;
            }

            if (contacts == 0) continue;

            // Find or create RingNeighbourhood
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
                new_rn.distance_to_center = dist_to_center;
                Vec3 d = atom_pos - geom.center;
                if (d.norm() > CalculatorConfig::Get("near_zero_vector_norm_threshold"))
                    new_rn.direction_to_center = d.normalized();

                double z = d.dot(geom.normal);
                Vec3 d_plane = d - z * geom.normal;
                new_rn.z = z;
                new_rn.rho = d_plane.norm();
                new_rn.theta = std::atan2(d_plane.norm(), std::abs(z));

                ca.ring_neighbours.push_back(new_rn);
                rn = &ca.ring_neighbours.back();
            }

            // Store dispersion results on RingNeighbourhood
            rn->disp_tensor = K_ring;
            rn->disp_spherical = SphericalTensor::Decompose(K_ring);
            rn->disp_scalar = s_ring;
            rn->disp_contacts = contacts;

            // Per-type accumulation
            int ti = ring.TypeIndexAsInt();
            if (ti >= 0 && ti < 8) {
                ca.per_type_disp_scalar_sum[ti] += s_ring;
                for (int c = 0; c < 5; ++c)
                    ca.per_type_disp_T2_sum[ti][c] += rn->disp_spherical.T2[c];
            }

            disp_total += K_ring;
            total_contacts += contacts;
            total_pairs++;
        }

        ca.disp_shielding_contribution = SphericalTensor::Decompose(disp_total);
    }

    OperationLog::Info(LogCalcOther, "DispersionResult::Compute",
        "atom_ring_pairs=" + std::to_string(total_pairs) +
        " vertex_contacts=" + std::to_string(total_contacts) +
        " bonded_exclusions=" + std::to_string(bonded_exclusions) +
        " rejected={" + filters.ReportRejections() + "}" +
        " atoms=" + std::to_string(n_atoms) +
        " rings=" + std::to_string(n_rings));

    return result_ptr;
}


SphericalTensor DispersionResult::SampleShieldingAt(Vec3 point) const {
    if (!conf_) return SphericalTensor{};

    const Protein& protein = conf_->ProteinRef();
    Mat3 K_total = Mat3::Zero();

    for (size_t ri = 0; ri < protein.RingCount(); ++ri) {
        const RingGeometry& geom = conf_->ring_geometries[ri];

        // Ring-level distance check
        double ring_dist = (point - geom.center).norm();
        if (ring_dist < CalculatorConfig::Get("singularity_guard_distance")) continue;
        if (ring_dist < geom.radius) continue;
        if (ring_dist > CalculatorConfig::Get("ring_current_spatial_cutoff")) continue;

        // Sum over ring vertices
        for (const auto& vertex : geom.vertices) {
            double r = (point - vertex).norm();
            if (r < CalculatorConfig::Get("singularity_guard_distance") || r > CalculatorConfig::Get("dispersion_vertex_distance_cutoff")) continue;

            auto vr = ComputeDispVertex(point, vertex, r);
            if (vr.valid) K_total += vr.K;
        }
    }

    return SphericalTensor::Decompose(K_total);
}


static void PackST_D(const SphericalTensor& st, double* out) {
    out[0] = st.T0;
    for (int i = 0; i < 3; ++i) out[1+i] = st.T1[i];
    for (int i = 0; i < 5; ++i) out[4+i] = st.T2[i];
}

int DispersionResult::WriteFeatures(const ProteinConformation& conf,
                                     const std::string& output_dir) const {
    const size_t N = conf.AtomCount();

    std::vector<double> shielding(N * 9);
    std::vector<double> per_type_T0(N * 8);
    std::vector<double> per_type_T2(N * 40);

    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        PackST_D(ca.disp_shielding_contribution, &shielding[i*9]);
        for (int t = 0; t < 8; ++t) {
            per_type_T0[i*8 + t] = ca.per_type_disp_scalar_sum[t];
            for (int c = 0; c < 5; ++c)
                per_type_T2[i*40 + t*5 + c] = ca.per_type_disp_T2_sum[t][c];
        }
    }

    NpyWriter::WriteFloat64(output_dir + "/disp_shielding.npy", shielding.data(), N, 9);
    NpyWriter::WriteFloat64(output_dir + "/disp_per_type_T0.npy", per_type_T0.data(), N, 8);
    NpyWriter::WriteFloat64(output_dir + "/disp_per_type_T2.npy", per_type_T2.data(), N, 40);
    return 3;
}

}  // namespace nmr

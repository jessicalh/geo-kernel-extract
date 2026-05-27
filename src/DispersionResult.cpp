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

// Number of T2 (rank-2) spherical-tensor components per ring type.
static constexpr int kT2Components = 5;


std::vector<std::type_index> DispersionResult::Dependencies() const {
    return {
        std::type_index(typeid(SpatialIndexResult)),
        std::type_index(typeid(GeometryResult))
    };
}


// ============================================================================
// Smooth switching function for the dispersion cutoff (CHARMM functional form).
//
// S(r) = 1                                            for r <= R_switch
// S(r) = (Rc²-r²)²(Rc²+2r²-3Rs²) / (Rc²-Rs²)³         for R_switch < r < R_cut
// S(r) = 0                                            for r >= R_cut
//
// C¹ continuous at both boundaries: S(R_switch)=1, S(R_cut)=0,
// S'(R_switch)=S'(R_cut)=0. Reference: Brooks et al., J. Comput. Chem. 4,
// 187 (1983) — the CHARMM switching function for non-bonded interactions.
//
// Why a smooth taper: the 1/r^6 interaction does not physically stop at any
// distance, so a hard cutoff would make features jump as atoms cross it. The
// taper truncates a convergent sum smoothly instead. The onset (R_switch,
// default 4.3 A) and cutoff (R_cut, default 5.0 A) are read from
// CalculatorConfig; the defaults below are illustrative. The 0.7 A taper
// width is the same order as MD position fluctuations (~0.5 A RMS), reducing
// cutoff jumps, and at R_cut the raw 1/r^6 term is below 0.5% of a typical
// 2 A contact before the switch takes it to zero.
// ============================================================================

static double DispersionSwitchingFunction(double r) {
    const double r_switch = CalculatorConfig::Get("dispersion_switching_onset_distance");
    const double r_cut    = CalculatorConfig::Get("dispersion_vertex_distance_cutoff");
    if (r <= r_switch) return 1.0;
    if (r >= r_cut)    return 0.0;

    const double switch_sq  = r_switch * r_switch;
    const double cutoff_sq  = r_cut * r_cut;
    const double r2         = r * r;
    const double numerator   = (cutoff_sq - r2) * (cutoff_sq - r2) * (cutoff_sq + 2.0 * r2 - 3.0 * switch_sq);
    const double denominator = (cutoff_sq - switch_sq) * (cutoff_sq - switch_sq) * (cutoff_sq - switch_sq);
    return numerator / denominator;
}


// ============================================================================
// per-vertex kernel: London dispersion from one ring vertex at one atom.
//
// Per vertex, with unit C6 = 1:
//
//   K_ab   = S(r) * (3 d_a d_b / r^8 - delta_ab / r^6)   (Angstrom^-6)
//   scalar = S(r) / r^6                                  (Angstrom^-6)
//
// where d = r_atom - r_vertex, r = |d|, and S(r) is the switching function.
// The tensor is traceless per vertex:
//   Tr(K) = S(r) * (3|d|^2/r^8 - 3/r^6) = S(r) * 0 = 0.
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

    // singularity guard: skip coincident / near-coincident points
    if (r < CalculatorConfig::Get("singularity_guard_distance")) return result;
    // outer cutoff: kernel is zero beyond the switching range
    if (r > CalculatorConfig::Get("dispersion_vertex_distance_cutoff")) return result;

    double S = DispersionSwitchingFunction(r);
    if (S < CalculatorConfig::Get("dispersion_switching_noise_floor")) return result;  // below noise floor

    Vec3 d_av = atom_pos - vertex_pos;   // d = r_atom − r_vertex (header symbol)
    double r2 = r * r;
    double r6 = r2 * r2 * r2;
    double r8 = r6 * r2;

    result.scalar = S / r6;

    // K_ab = S(r) * (3 d_a d_b / r^8 - delta_ab / r^6)
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            result.K(a, b) = S * (3.0 * d_av(a) * d_av(b) / r8
                                - (a == b ? 1.0 : 0.0) / r6);

    result.valid = true;
    return result;
}


// ============================================================================
// Build the set of ring vertices and atoms bonded to any vertex.
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
// Dataflow: for each atom, query nearby rings, reject near-field and
// through-bond rings, sum the per-vertex kernel over the survivors, and
// store the result (per ring-neighbourhood, per ring type, and as the
// atom total).
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

    // Ring near-field filter (DipolarNearFieldFilter, source_extent = ring
    // diameter): the discrete vertex sum is invalid when the field point is
    // inside the ring.
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
            // reject near-field atom
            if (!filters.AcceptAll(ctx)) {
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

            // reject bonded atom: skip this ring entirely if the field atom is
            // bonded to any vertex (part of, or immediately adjacent to, the
            // ring) — the through-space 1/r^6 kernel does not model through-bond
            // coupling.
            if (ring_bonded[ri].count(ai)) {
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

            // record taper parameters
            if (recorded_rings.insert(ri).second) {
                choices.Record(CalculatorId::Dispersion, ri, "dispersion taper",
                    [&ring](GeometryChoice& gc) {
                        AddRing(gc, &ring, EntityRole::Source, EntityOutcome::Included);
                        AddNumber(gc, "switch_onset", CalculatorConfig::Get("dispersion_switching_onset_distance"), "A");
                        AddNumber(gc, "cutoff", CalculatorConfig::Get("dispersion_vertex_distance_cutoff"), "A");
                    });
            }

            // --- vertex kernel sum ---
            Mat3 K_ring = Mat3::Zero();
            double s_ring = 0.0;
            int contacts = 0;

            for (size_t vi = 0; vi < ring.atom_indices.size(); ++vi) {
                Vec3 vpos = geom.vertices[vi];
                double r = (atom_pos - vpos).norm();

                DispVertexResult vertex_result = ComputeDispVertex(atom_pos, vpos, r);
                if (!vertex_result.valid) {
                    // switching floor: vertex tapered below noise floor
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

                K_ring += vertex_result.K;
                s_ring += vertex_result.scalar;
                contacts++;
            }

            if (contacts == 0) continue;

            // locate or create the ring-neighbourhood record
            RingNeighbourhood* ring_neighbour = nullptr;
            for (auto& existing : ca.ring_neighbours) {
                if (existing.ring_index == ri) {
                    ring_neighbour = &existing;
                    break;
                }
            }
            if (!ring_neighbour) {
                RingNeighbourhood new_rn;
                new_rn.ring_index = ri;
                new_rn.ring_type = ring.type_index;
                new_rn.distance_to_center = dist_to_center;
                Vec3 center_to_atom = atom_pos - geom.center;
                if (center_to_atom.norm() > CalculatorConfig::Get("near_zero_vector_norm_threshold"))
                    new_rn.direction_to_center = center_to_atom.normalized();

                // ring-frame coordinates (z along normal, rho in plane, theta polar)
                double z = center_to_atom.dot(geom.normal);
                Vec3 d_plane = center_to_atom - z * geom.normal;
                new_rn.z = z;
                new_rn.rho = d_plane.norm();
                new_rn.theta = std::atan2(d_plane.norm(), std::abs(z));

                ca.ring_neighbours.push_back(new_rn);
                ring_neighbour = &ca.ring_neighbours.back();
            }

            // store kernel on the ring-neighbourhood record (spherical decomposition)
            ring_neighbour->disp_tensor = K_ring;
            ring_neighbour->disp_spherical = SphericalTensor::Decompose(K_ring);
            ring_neighbour->disp_scalar = s_ring;
            ring_neighbour->disp_contacts = contacts;

            // Per-type accumulation (Pro pyrrolidine index 8 is non-aromatic, excluded here)
            int ti = ring.TypeIndexAsInt();
            if (ti >= 0 && ti < kAromaticRingTypeCount) {
                ca.per_type_disp_scalar_sum[ti] += s_ring;
                for (int c = 0; c < kT2Components; ++c)
                    ca.per_type_disp_T2_sum[ti][c] += ring_neighbour->disp_spherical.T2[c];
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


SphericalTensor DispersionResult::SampleKernelAt(Vec3 point) const {
    if (!conf_) return SphericalTensor{};

    const Protein& protein = conf_->ProteinRef();
    Mat3 K_total = Mat3::Zero();

    for (size_t ri = 0; ri < protein.RingCount(); ++ri) {
        const RingGeometry& geom = conf_->ring_geometries[ri];

        // grid-sampling guard: skip points inside the ring-radius sphere and
        // out-of-range points (raw field sample, no provenance recording)
        double ring_dist = (point - geom.center).norm();
        if (ring_dist < CalculatorConfig::Get("singularity_guard_distance")) continue;
        if (ring_dist < geom.radius) continue;
        if (ring_dist > CalculatorConfig::Get("ring_current_spatial_cutoff")) continue;

        // Sum over ring vertices
        for (const auto& vertex : geom.vertices) {
            double r = (point - vertex).norm();
            if (r < CalculatorConfig::Get("singularity_guard_distance") || r > CalculatorConfig::Get("dispersion_vertex_distance_cutoff")) continue;

            auto vertex_result = ComputeDispVertex(point, vertex, r);
            if (vertex_result.valid) K_total += vertex_result.K;
        }
    }

    return SphericalTensor::Decompose(K_total);
}


int DispersionResult::WriteFeatures(const ProteinConformation& conf,
                                     const std::string& output_dir) const {
    const size_t N = conf.AtomCount();

    // per_type_T2 layout: 8 ring types × 5 T2 components = 40
    constexpr int kPerTypeT2Width = kAromaticRingTypeCount * kT2Components;  // 40

    std::vector<double> shielding(N * 9);
    std::vector<double> per_type_T0(N * kAromaticRingTypeCount);
    std::vector<double> per_type_T2(N * kPerTypeT2Width);

    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        ca.disp_shielding_contribution.PackFull9(&shielding[i*9]);
        for (int t = 0; t < kAromaticRingTypeCount; ++t) {
            per_type_T0[i*kAromaticRingTypeCount + t] = ca.per_type_disp_scalar_sum[t];
            for (int c = 0; c < kT2Components; ++c)
                per_type_T2[i*kPerTypeT2Width + t*kT2Components + c] = ca.per_type_disp_T2_sum[t][c];
        }
    }

    NpyWriter::WriteFloat64(output_dir + "/disp_shielding.npy", shielding.data(), N, 9);
    NpyWriter::WriteFloat64(output_dir + "/disp_per_type_T0.npy", per_type_T0.data(), N, kAromaticRingTypeCount);
    NpyWriter::WriteFloat64(output_dir + "/disp_per_type_T2.npy", per_type_T2.data(), N, kPerTypeT2Width);
    return 3;
}

}  // namespace nmr

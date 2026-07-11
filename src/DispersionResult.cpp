#include "DispersionResult.h"
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
// Per-vertex clean London dispersion scalar from one ring vertex at one atom.
//
// Per vertex, with unit C6 = 1:
//
//   scalar = S(r) / r^6                                  (Angstrom^-6)
//
// where d = r_atom - r_vertex, r = |d|, and S(r) is the switching function.
// The former unit-C6 rank-2 tensor is intentionally no longer computed
// or emitted.
// ============================================================================

namespace dispersion_detail {

VertexResult ComputeVertex(double r) {

    VertexResult result;
    result.distance = r;

    // singularity guard: skip coincident / near-coincident points
    if (r < CalculatorConfig::Get("singularity_guard_distance")) return result;
    // outer cutoff: kernel is zero beyond the switching range
    if (r > CalculatorConfig::Get("dispersion_vertex_distance_cutoff")) return result;

    double S = DispersionSwitchingFunction(r);
    if (S < CalculatorConfig::Get("dispersion_switching_noise_floor")) return result;  // below noise floor

    double r2 = r * r;
    double r6 = r2 * r2 * r2;

    result.scalar = S / r6;

    result.valid = true;
    return result;
}


std::vector<VertexResult> ComputeRingVertices(
        const Vec3& atom_pos,
        const RingGeometry& geom) {
    std::vector<VertexResult> results;
    results.reserve(geom.vertices.size());
    for (const Vec3& vertex : geom.vertices) {
        results.push_back(ComputeVertex((atom_pos - vertex).norm()));
    }
    return results;
}

}  // namespace dispersion_detail


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

    OperationLog::Info(LogCalcOther, "DispersionResult::Compute",
        "vertex range: [MIN_DISTANCE=" + std::to_string(CalculatorConfig::Get("singularity_guard_distance")) +
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

        for (size_t ri : nearby_rings) {
            const Ring& ring = protein.RingAt(ri);
            const RingGeometry& geom = conf.ring_geometries[ri];

            if (geom.vertices.empty()) continue;

            double dist_to_center = (atom_pos - geom.center).norm();

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
            double s_ring = 0.0;
            int contacts = 0;

            const std::vector<dispersion_detail::VertexResult>
                vertex_results =
                    dispersion_detail::ComputeRingVertices(atom_pos, geom);
            for (const auto& vertex_result : vertex_results) {
                const double r = vertex_result.distance;
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

                s_ring += vertex_result.scalar;
                contacts++;
            }

            if (contacts == 0) continue;

            RingNeighbourhood& ring_neighbour_ref =
                EnsureRingNeighbourGeometry(ca, geom, ring, ri, atom_pos);
            RingNeighbourhood* ring_neighbour = &ring_neighbour_ref;

            // Store only the clean scalar/contact rescue.
            ring_neighbour->disp_scalar = s_ring;
            ring_neighbour->disp_contacts = contacts;

            // Per-type accumulation (Pro pyrrolidine index 8 is non-aromatic, excluded here)
            int ti = ring.TypeIndexAsInt();
            if (ti >= 0 && ti < kAromaticRingTypeCount) {
                ca.per_type_disp_scalar_sum[ti] += s_ring;
            }

            total_contacts += contacts;
            total_pairs++;
        }
    }

    OperationLog::Info(LogCalcOther, "DispersionResult::Compute",
        "atom_ring_pairs=" + std::to_string(total_pairs) +
        " vertex_contacts=" + std::to_string(total_contacts) +
        " bonded_exclusions=" + std::to_string(bonded_exclusions) +
        " atoms=" + std::to_string(n_atoms) +
        " rings=" + std::to_string(n_rings));

    return result_ptr;
}


int DispersionResult::WriteFeatures(const ProteinConformation& conf,
                                     const std::string& output_dir) const {
    const size_t N = conf.AtomCount();

    std::vector<double> per_type_T0(N * kAromaticRingTypeCount);

    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        for (int t = 0; t < kAromaticRingTypeCount; ++t) {
            per_type_T0[i*kAromaticRingTypeCount + t] = ca.per_type_disp_scalar_sum[t];
        }
    }

    NpyWriter::WriteFloat64(
        output_dir + "/disp_per_type_T0.npy",
        per_type_T0.data(), N, kAromaticRingTypeCount);
    NpyWriter::WriteFloat64(
        output_dir + "/aromatic_r6_proximity_per_type_T0.npy",
        per_type_T0.data(), N, kAromaticRingTypeCount);
    return 2;
}

}  // namespace nmr

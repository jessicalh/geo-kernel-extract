#include "CoulombResult.h"
#include "Protein.h"
#include "ChargeAssignmentResult.h"
#include "SpatialIndexResult.h"
#include "ApbsFieldResult.h"
#include "KernelEvaluationFilter.h"
#include "PhysicalConstants.h"
#include "CalculatorConfig.h"
#include "GeometryChoice.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <cmath>
#include <vector>

namespace nmr {


std::vector<std::type_index> CoulombResult::Dependencies() const {
    return {
        std::type_index(typeid(ChargeAssignmentResult)),
        std::type_index(typeid(SpatialIndexResult))
    };
}


// ============================================================================
// CoulombResult::Compute
//
// E_a(i) = ke * sum_{j!=i} q_j * (r_i - r_j)_a / |r_i - r_j|^3
//                                   (j within coulomb_efield_cutoff of i)
//
// V_ab(i) = ke * sum_{j!=i} q_j * [3 (r_i-r_j)_a (r_i-r_j)_b / |r_i-r_j|^5
//                                   - delta_ab / |r_i-r_j|^3]
//                                   (j within coulomb_efield_cutoff of i)
//
// ke = 14.3996 V*A  (Coulomb's constant in {e, A, eV} units)
//
// Decomposed by source atom classification:
//   backbone:  N, CA, C, O, H, HA, CB (from residue backbone cache)
//   aromatic:  atoms that are members of any Ring
//   sidechain: everything else
// ============================================================================

std::unique_ptr<CoulombResult> CoulombResult::Compute(
        ProteinConformation& conf) {

    OperationLog::Scope scope("CoulombResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()));

    const Protein& protein = conf.ProteinRef();
    const size_t n_atoms = conf.AtomCount();

    auto result_ptr = std::make_unique<CoulombResult>();
    result_ptr->conf_ = &conf;

    // atom classes (backbone / aromatic) from Residue backbone cache and
    // Ring atom membership; no EnrichmentResult dependency.
    std::vector<bool> is_backbone(n_atoms, false);
    std::vector<bool> is_aromatic_atom(n_atoms, false);

    for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        auto mark_bb = [&](size_t idx) {
            if (idx != Residue::NONE && idx < n_atoms) is_backbone[idx] = true;
        };
        mark_bb(res.N);
        mark_bb(res.CA);
        mark_bb(res.C);
        mark_bb(res.O);
        mark_bb(res.H);
        mark_bb(res.HA);
        mark_bb(res.CB);
    }

    for (size_t ri = 0; ri < protein.RingCount(); ++ri) {
        for (size_t ai : protein.RingAt(ri).atom_indices) {
            if (ai < n_atoms) is_aromatic_atom[ai] = true;
        }
    }

    // primary bond direction (for E_bond_proj). H: parent heavy atom -> H.
    // Heavy: first bond (arbitrary but consistent).
    std::vector<Vec3> primary_bond_dir(n_atoms, Vec3::Zero());
    for (size_t ai = 0; ai < n_atoms; ++ai) {
        const Atom& atom = protein.AtomAt(ai);
        if (atom.element == Element::H && atom.parent_atom_index != SIZE_MAX) {
            Vec3 parent_to_hydrogen = conf.PositionAt(ai) - conf.PositionAt(atom.parent_atom_index);
            double len = parent_to_hydrogen.norm();
            if (len > CalculatorConfig::Get("near_zero_vector_norm_threshold")) primary_bond_dir[ai] = parent_to_hydrogen / len;
        } else if (!atom.bond_indices.empty()) {
            const Bond& b = protein.BondAt(atom.bond_indices[0]);
            size_t other = (b.atom_index_a == ai) ? b.atom_index_b : b.atom_index_a;
            Vec3 atom_to_bond_neighbor = conf.PositionAt(other) - conf.PositionAt(ai);
            double len = atom_to_bond_neighbor.norm();
            if (len > CalculatorConfig::Get("near_zero_vector_norm_threshold")) primary_bond_dir[ai] = atom_to_bond_neighbor / len;
        }
    }

    // ------------------------------------------------------------------
    // Coulomb sum with spatial cutoff
    // ------------------------------------------------------------------

    // Filter set: MinDistanceFilter plus SelfSourceFilter (field undefined at
    // source itself).
    // Coulomb is a point-source sum — no DipolarNearFieldFilter needed.
    KernelFilterSet filters;
    filters.Add(std::make_unique<MinDistanceFilter>());
    filters.Add(std::make_unique<SelfSourceFilter>());

    GeometryChoiceBuilder choices(conf);

    const auto& spatial = conf.Result<SpatialIndexResult>();
    const double coulomb_cutoff = CalculatorConfig::Get("coulomb_efield_cutoff");

    // diagnostic count: feeds the summary log line, not the field sum
    int aromatic_source_count = 0;
    for (size_t j = 0; j < n_atoms; ++j)
        if (is_aromatic_atom[j]) aromatic_source_count++;

    for (size_t i = 0; i < n_atoms; ++i) {
        Vec3 pos_i = conf.PositionAt(i);

        Vec3 E_total = Vec3::Zero();
        Vec3 E_backbone = Vec3::Zero();
        Vec3 E_sidechain = Vec3::Zero();
        Vec3 E_aromatic = Vec3::Zero();

        Mat3 EFG_total = Mat3::Zero();
        Mat3 EFG_backbone = Mat3::Zero();
        Mat3 EFG_sidechain = Mat3::Zero();
        Mat3 EFG_aromatic = Mat3::Zero();

        // diagnostic count: feeds one stored scalar, not the field sum
        int n_sidechain_aromatic_sources = 0;

        auto neighbours = spatial.AtomsWithinRadius(pos_i, coulomb_cutoff);
        // legacy radius-search diagnostic: n_atoms - 1 minus the raw neighbour
        // match count; NOT the count dropped by filters / charge floor
        // (those thin the summed set further below)
        int sources_outside_radius = static_cast<int>(n_atoms) - 1
                                    - static_cast<int>(neighbours.size());

        // --- per-source field + EFG sum ---
        for (size_t j : neighbours) {
            // self-exclusion
            KernelEvaluationContext ctx;
            ctx.atom_index = i;
            ctx.source_atom_a = j;
            ctx.distance = (pos_i - conf.PositionAt(j)).norm();
            if (!filters.AcceptAll(ctx)) continue;

            double q_j = conf.AtomAt(j).partial_charge;
            if (std::abs(q_j) < CalculatorConfig::Get("coulomb_charge_noise_floor")) continue;

            Vec3 r_vec = pos_i - conf.PositionAt(j);
            double r_mag = r_vec.norm();

            double r3 = r_mag * r_mag * r_mag;
            double r5 = r3 * r_mag * r_mag;

            // E_a = q_j * r_a / r^3
            Vec3 E_j = q_j * r_vec / r3;

            // V_ab = q_j * (3 r_a r_b / r^5 - delta_ab / r^3)
            Mat3 V_j = q_j * (3.0 * r_vec * r_vec.transpose() / r5
                              - Mat3::Identity() / r3);

            E_total += E_j;
            EFG_total += V_j;

            // Classify source atom
            if (is_aromatic_atom[j]) {
                E_aromatic += E_j;
                EFG_aromatic += V_j;
                // aromatic-and-sidechain sources contributing to this target
                if (!is_backbone[j]) n_sidechain_aromatic_sources++;
            } else if (is_backbone[j]) {
                E_backbone += E_j;
                EFG_backbone += V_j;
            } else {
                E_sidechain += E_j;
                EFG_sidechain += V_j;
            }
        }

        // --- units, tracelessness, sanitise, clamp ---
        // Apply Coulomb constant: convert from e/A^2 to V/A
        E_total     *= COULOMB_KE;
        E_backbone  *= COULOMB_KE;
        E_sidechain *= COULOMB_KE;
        E_aromatic  *= COULOMB_KE;
        EFG_total     *= COULOMB_KE;
        EFG_backbone  *= COULOMB_KE;
        EFG_sidechain *= COULOMB_KE;
        EFG_aromatic  *= COULOMB_KE;

        // Traceless projection on all EFG matrices.
        // Each individual term is traceless (Gauss's law), but floating-point
        // accumulation breaks this. Project to enforce the physics.
        auto project_traceless = [](Mat3& m) {
            m -= (m.trace() / 3.0) * Mat3::Identity();
        };
        project_traceless(EFG_total);
        project_traceless(EFG_backbone);
        project_traceless(EFG_sidechain);  // kept for symmetry; not stored or written to NPY
        project_traceless(EFG_aromatic);

        // Sanitise NaN/Inf
        auto sanitise_vec = [](Vec3& v) {
            for (int d = 0; d < 3; ++d)
                if (std::isnan(v(d)) || std::isinf(v(d))) { v = Vec3::Zero(); return; }
        };
        auto sanitise_mat = [](Mat3& m) {
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                    if (std::isnan(m(a,b)) || std::isinf(m(a,b))) m(a,b) = 0.0;
        };
        sanitise_vec(E_total);
        sanitise_vec(E_backbone);
        sanitise_vec(E_sidechain);
        sanitise_vec(E_aromatic);
        sanitise_mat(EFG_total);
        sanitise_mat(EFG_backbone);
        sanitise_mat(EFG_sidechain);
        sanitise_mat(EFG_aromatic);

        // Clamp extreme E-field magnitudes. Scope: E vectors only — the EFG
        // matrices are guarded separately (project_traceless + sanitise_mat +
        // the 0.1 A MinDistanceFilter); no consumer couples E and EFG.
        double E_mag = E_total.norm();
        if (E_mag > CalculatorConfig::Get("efield_magnitude_sanity_clamp")) {
            double scale = CalculatorConfig::Get("efield_magnitude_sanity_clamp") / E_mag;

            // ---- GeometryChoice: E-field clamp ----
            choices.Record(CalculatorId::Coulomb, i, "E-field clamp",
                [&conf, i, E_mag, scale](GeometryChoice& gc) {
                    AddAtom(gc, &conf.AtomAt(i), i, EntityRole::Target, EntityOutcome::Triggered);
                    AddNumber(gc, "actual_E_magnitude", E_mag, "V/A");
                    AddNumber(gc, "scale_factor", scale, "");
                });

            E_total     *= scale;
            E_backbone  *= scale;
            E_sidechain *= scale;
            E_aromatic  *= scale;
        }

        // ---- GeometryChoice: cutoff summary ----
        if (sources_outside_radius > 0) {
            int sources_within = static_cast<int>(neighbours.size());
            choices.Record(CalculatorId::Coulomb, i, "coulomb cutoff",
                [&conf, i, sources_within, sources_outside_radius, coulomb_cutoff](GeometryChoice& gc) {
                    AddAtom(gc, &conf.AtomAt(i), i, EntityRole::Target, EntityOutcome::Included);
                    AddNumber(gc, "sources_within_cutoff", static_cast<double>(sources_within), "count");
                    // recorded key kept as "sources_beyond_cutoff" (audit contract);
                    // value is the legacy n_atoms-1-minus-match-count diagnostic
                    AddNumber(gc, "sources_beyond_cutoff", static_cast<double>(sources_outside_radius), "count");
                    AddNumber(gc, "cutoff_distance", coulomb_cutoff, "A");
                });
        }

        // ===== store: physics outputs, then derived scalars =====
        auto& ca = conf.MutableAtomAt(i);

        ca.coulomb_E_total     = E_total;
        ca.coulomb_E_backbone  = E_backbone;
        ca.coulomb_E_sidechain = E_sidechain;
        ca.coulomb_E_aromatic  = E_aromatic;

        ca.coulomb_EFG_total   = EFG_total;
        ca.coulomb_EFG_total_spherical = SphericalTensor::Decompose(EFG_total);

        ca.coulomb_EFG_backbone = EFG_backbone;
        ca.coulomb_EFG_backbone_spherical = SphericalTensor::Decompose(EFG_backbone);

        ca.coulomb_EFG_aromatic = EFG_aromatic;
        ca.coulomb_EFG_aromatic_spherical = SphericalTensor::Decompose(EFG_aromatic);

        // bond projection
        ca.coulomb_E_magnitude = E_total.norm();
        // E projected along primary bond direction (for Buckingham E_z)
        ca.coulomb_E_bond_proj = E_total.dot(primary_bond_dir[i]);

        // backbone alignment
        // Backbone projection: component of E_backbone along E_total direction.
        // Positive = backbone field aligned with total; negative = opposed.
        // Bounded by |E_backbone|. Stable near cancellation (unlike |bb|/|total|).
        // This is a signed projection (V/A), not a ratio — field name is historical.
        if (ca.coulomb_E_magnitude > CalculatorConfig::Get("near_zero_vector_norm_threshold")) {
            Vec3 E_hat = E_total / ca.coulomb_E_magnitude;
            ca.coulomb_E_backbone_frac = E_backbone.dot(E_hat);
        } else {
            ca.coulomb_E_backbone_frac = 0.0;
        }

        // aromatic scalars
        ca.aromatic_E_magnitude = E_aromatic.norm();
        ca.aromatic_E_bond_proj = E_aromatic.dot(primary_bond_dir[i]);
        ca.aromatic_n_sidechain_atoms = n_sidechain_aromatic_sources;

        // Solvent contribution: APBS (solvated) minus vacuum Coulomb.
        // Only meaningful if ApbsFieldResult is present.
        // Both are in V/A, so this is a proper subtraction.
        if (conf.HasResult<ApbsFieldResult>()) {
            ca.coulomb_E_solvent = ca.apbs_efield - E_total;
            ca.coulomb_EFG_solvent = ca.apbs_efg - EFG_total;
        }

        // T2 only; Buckingham T0 not included here (see CoulombResult.h).
        ca.coulomb_shielding_contribution = SphericalTensor::Decompose(EFG_total);
    }

    OperationLog::Info(LogCalcOther, "CoulombResult::Compute",
        "atoms=" + std::to_string(n_atoms) +
        " aromatic_sources=" + std::to_string(aromatic_source_count) +
        " rejected={" + filters.ReportRejections() + "}");

    return result_ptr;
}


// ============================================================================
// Query methods
// ============================================================================

Vec3 CoulombResult::EFieldAt(size_t atom_index) const {
    return conf_->AtomAt(atom_index).coulomb_E_total;
}

Mat3 CoulombResult::EFGAt(size_t atom_index) const {
    return conf_->AtomAt(atom_index).coulomb_EFG_total;
}

SphericalTensor CoulombResult::EFGSphericalAt(size_t atom_index) const {
    return conf_->AtomAt(atom_index).coulomb_EFG_total_spherical;
}


Vec3 CoulombResult::SampleEFieldAt(Vec3 point) const {
    if (!conf_) return Vec3::Zero();

    Vec3 E = Vec3::Zero();

    for (size_t j = 0; j < conf_->AtomCount(); ++j) {
        double q = conf_->AtomAt(j).partial_charge;
        if (std::abs(q) < CalculatorConfig::Get("coulomb_charge_noise_floor")) continue;

        Vec3 d = point - conf_->PositionAt(j);
        double r = d.norm();
        if (r < CalculatorConfig::Get("singularity_guard_distance")) continue;

        double r3 = r * r * r;
        E += q * d / r3;
    }

    return E * COULOMB_KE;  // V/A
}


// ============================================================================
// WriteFeatures: coulomb_shielding (9), E-field (3), EFG decompositions,
// scalar features (magnitude, bond projection, backbone fraction).
// ============================================================================

int CoulombResult::WriteFeatures(const ProteinConformation& conf,
                                  const std::string& output_dir) const {
    const size_t N = conf.AtomCount();

    std::vector<double> shielding(N * 9);
    std::vector<double> efield(N * 3);
    // EFG: T2 only (5 components). Coulomb EFG built from q·(3r⊗r/r⁵−I/r³) —
    // symmetric per charge → T0+T1 structural zeros.
    std::vector<double> efg_bb(N * 5);
    std::vector<double> efg_aro(N * 5);
    std::vector<double> scalars(N * 4);

    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        ca.coulomb_shielding_contribution.PackFull9(&shielding[i*9]);

        efield[i*3+0] = ca.coulomb_E_total.x();
        efield[i*3+1] = ca.coulomb_E_total.y();
        efield[i*3+2] = ca.coulomb_E_total.z();

        ca.coulomb_EFG_backbone_spherical.PackT2(&efg_bb[i*5]);
        ca.coulomb_EFG_aromatic_spherical.PackT2(&efg_aro[i*5]);

        scalars[i*4+0] = ca.coulomb_E_magnitude;
        scalars[i*4+1] = ca.coulomb_E_bond_proj;
        scalars[i*4+2] = ca.coulomb_E_backbone_frac;
        scalars[i*4+3] = ca.aromatic_E_magnitude;
    }

    NpyWriter::WriteFloat64(output_dir + "/coulomb_shielding.npy", shielding.data(), N, 9);
    NpyWriter::WriteFloat64(output_dir + "/coulomb_E.npy", efield.data(), N, 3);
    NpyWriter::WriteFloat64(output_dir + "/coulomb_efg_backbone.npy", efg_bb.data(), N, 5);
    NpyWriter::WriteFloat64(output_dir + "/coulomb_efg_aromatic.npy", efg_aro.data(), N, 5);
    NpyWriter::WriteFloat64(output_dir + "/coulomb_scalars.npy", scalars.data(), N, 4);
    return 5;
}

}  // namespace nmr

#include "CoulombResult.h"
#include "Protein.h"
#include "ChargeAssignmentResult.h"
#include "SpatialIndexResult.h"
#include "ApbsFieldResult.h"
#include "KernelEvaluationFilter.h"
#include "PhysicalConstants.h"
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
//
// V_ab(i) = ke * sum_{j!=i} q_j * [3 (r_i-r_j)_a (r_i-r_j)_b / |r_i-r_j|^5
//                                   - delta_ab / |r_i-r_j|^3]
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

    // ------------------------------------------------------------------
    // Pre-build atom classification vectors from topology (no EnrichmentResult
    // dependency — uses Residue backbone cache and Ring atom membership).
    // ------------------------------------------------------------------

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

    // ------------------------------------------------------------------
    // For E_bond_proj: find each atom's primary bond direction.
    // H atoms: direction from parent heavy atom to H.
    // Heavy atoms: direction of first bond (arbitrary but consistent).
    // ------------------------------------------------------------------

    std::vector<Vec3> primary_bond_dir(n_atoms, Vec3::Zero());
    for (size_t ai = 0; ai < n_atoms; ++ai) {
        const Atom& atom = protein.AtomAt(ai);
        if (atom.element == Element::H && atom.parent_atom_index != SIZE_MAX) {
            // H atom: bond direction from parent to H
            Vec3 d = conf.PositionAt(ai) - conf.PositionAt(atom.parent_atom_index);
            double len = d.norm();
            if (len > NEAR_ZERO_NORM) primary_bond_dir[ai] = d / len;
        } else if (!atom.bond_indices.empty()) {
            // Heavy atom: first bond direction
            const Bond& b = protein.BondAt(atom.bond_indices[0]);
            size_t other = (b.atom_index_a == ai) ? b.atom_index_b : b.atom_index_a;
            Vec3 d = conf.PositionAt(other) - conf.PositionAt(ai);
            double len = d.norm();
            if (len > NEAR_ZERO_NORM) primary_bond_dir[ai] = d / len;
        }
    }

    // ------------------------------------------------------------------
    // Main N^2 Coulomb sum
    // ------------------------------------------------------------------

    // Filter set: SelfSourceFilter (field undefined at source itself).
    // Coulomb is a point-source sum — no DipolarNearFieldFilter needed.
    KernelFilterSet filters;
    filters.Add(std::make_unique<MinDistanceFilter>());
    filters.Add(std::make_unique<SelfSourceFilter>());

    GeometryChoiceBuilder choices(conf);

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

        int n_sidechain_aromatic_sources = 0;

        for (size_t j = 0; j < n_atoms; ++j) {
            // Self-exclusion via filter framework (not inline check)
            KernelEvaluationContext ctx;
            ctx.atom_index = i;
            ctx.source_atom_a = j;
            ctx.distance = (pos_i - conf.PositionAt(j)).norm();
            if (!filters.AcceptAll(ctx)) continue;

            double q_j = conf.AtomAt(j).partial_charge;
            if (std::abs(q_j) < 1e-15) continue;

            Vec3 r = pos_i - conf.PositionAt(j);
            double r_mag = r.norm();

            double r3 = r_mag * r_mag * r_mag;
            double r5 = r3 * r_mag * r_mag;

            // E_a = q_j * r_a / r^3
            Vec3 E_j = q_j * r / r3;

            // V_ab = q_j * (3 r_a r_b / r^5 - delta_ab / r^3)
            Mat3 V_j = q_j * (3.0 * r * r.transpose() / r5
                              - Mat3::Identity() / r3);

            E_total += E_j;
            EFG_total += V_j;

            // Classify source atom
            if (is_aromatic_atom[j]) {
                E_aromatic += E_j;
                EFG_aromatic += V_j;
                // Count aromatic sidechain sources near this atom
                // (atoms on aromatic residues that are sidechain)
                if (!is_backbone[j]) n_sidechain_aromatic_sources++;
            } else if (is_backbone[j]) {
                E_backbone += E_j;
                EFG_backbone += V_j;
            } else {
                E_sidechain += E_j;
                EFG_sidechain += V_j;
            }
        }

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
        project_traceless(EFG_sidechain);
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

        // Clamp extreme E-field magnitudes
        double E_mag = E_total.norm();
        if (E_mag > APBS_SANITY_LIMIT) {
            double scale = APBS_SANITY_LIMIT / E_mag;

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

        // ------------------------------------------------------------------
        // Store on ConformationAtom
        // ------------------------------------------------------------------
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

        // Derived scalars
        ca.coulomb_E_magnitude = E_total.norm();

        // E projected along primary bond direction (for Buckingham E_z)
        ca.coulomb_E_bond_proj = E_total.dot(primary_bond_dir[i]);

        // Backbone projection: component of E_backbone along E_total direction.
        // Positive = backbone field aligned with total; negative = opposed.
        // Bounded by |E_backbone|. Stable near cancellation (unlike |bb|/|total|).
        if (ca.coulomb_E_magnitude > NEAR_ZERO_NORM) {
            Vec3 E_hat = E_total / ca.coulomb_E_magnitude;
            ca.coulomb_E_backbone_frac = E_backbone.dot(E_hat);
        } else {
            ca.coulomb_E_backbone_frac = 0.0;
        }

        // Aromatic E-field derived scalars
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

        // Shielding contribution: SphericalTensor of the total EFG.
        //
        // WARNING: This is ONLY the T2 geometric kernel. T0 is zero here
        // (EFG is traceless by Gauss's law). The actual T0 Coulomb shielding
        // comes from Buckingham's A*E_z + B*E_z^2, which requires element-
        // dependent parameters and a bond direction — it is NOT a pure
        // geometric kernel like McConnell's M/r^3.
        //
        // Unlike McConnell (where χ·K tensor contraction produces an
        // asymmetric tensor with non-zero T0+T1+T2 from geometry alone),
        // Coulomb has two SEPARATE geometric kernels:
        //   E_a  (rank-1) → T0 shielding via Buckingham A,B parameters
        //   V_ab (rank-2, symmetric, traceless) → T2 shielding via γ
        // There is no single "full tensor" that unifies them.
        //
        // This field stores Decompose(EFG) only. When Buckingham parameters
        // are available (from ParameterCorrectionResult), the full shielding
        // contribution including T0 from E should be computed and stored.
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

    const Protein& protein = conf_->ProteinRef();
    Vec3 E = Vec3::Zero();

    for (size_t j = 0; j < conf_->AtomCount(); ++j) {
        double q = conf_->AtomAt(j).partial_charge;
        if (std::abs(q) < 1e-15) continue;

        Vec3 d = point - conf_->PositionAt(j);
        double r = d.norm();
        if (r < MIN_DISTANCE) continue;

        double r3 = r * r * r;
        E += q * d / r3;
    }

    return E * COULOMB_KE;  // V/A
}


// ============================================================================
// WriteFeatures: coulomb_shielding (9), E-field (3), EFG decompositions,
// scalar features (magnitude, bond projection, backbone fraction).
// ============================================================================

static void PackST_C(const SphericalTensor& st, double* out) {
    out[0] = st.T0;
    for (int i = 0; i < 3; ++i) out[1+i] = st.T1[i];
    for (int i = 0; i < 5; ++i) out[4+i] = st.T2[i];
}

int CoulombResult::WriteFeatures(const ProteinConformation& conf,
                                  const std::string& output_dir) const {
    const size_t N = conf.AtomCount();

    std::vector<double> shielding(N * 9);
    std::vector<double> efield(N * 3);
    std::vector<double> efg_bb(N * 9);
    std::vector<double> efg_aro(N * 9);
    std::vector<double> scalars(N * 4);

    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        PackST_C(ca.coulomb_shielding_contribution, &shielding[i*9]);

        efield[i*3+0] = ca.coulomb_E_total.x();
        efield[i*3+1] = ca.coulomb_E_total.y();
        efield[i*3+2] = ca.coulomb_E_total.z();

        PackST_C(ca.coulomb_EFG_backbone_spherical, &efg_bb[i*9]);
        PackST_C(ca.coulomb_EFG_aromatic_spherical, &efg_aro[i*9]);

        scalars[i*4+0] = ca.coulomb_E_magnitude;
        scalars[i*4+1] = ca.coulomb_E_bond_proj;
        scalars[i*4+2] = ca.coulomb_E_backbone_frac;
        scalars[i*4+3] = ca.aromatic_E_magnitude;
    }

    NpyWriter::WriteFloat64(output_dir + "/coulomb_shielding.npy", shielding.data(), N, 9);
    NpyWriter::WriteFloat64(output_dir + "/coulomb_E.npy", efield.data(), N, 3);
    NpyWriter::WriteFloat64(output_dir + "/coulomb_efg_backbone.npy", efg_bb.data(), N, 9);
    NpyWriter::WriteFloat64(output_dir + "/coulomb_efg_aromatic.npy", efg_aro.data(), N, 9);
    NpyWriter::WriteFloat64(output_dir + "/coulomb_scalars.npy", scalars.data(), N, 4);
    return 5;
}

}  // namespace nmr

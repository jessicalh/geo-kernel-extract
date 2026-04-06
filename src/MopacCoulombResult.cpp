#include "MopacCoulombResult.h"
#include "Protein.h"
#include "MopacResult.h"
#include "SpatialIndexResult.h"
#include "KernelEvaluationFilter.h"
#include "PhysicalConstants.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <cmath>
#include <vector>

namespace nmr {


std::vector<std::type_index> MopacCoulombResult::Dependencies() const {
    return {
        std::type_index(typeid(MopacResult)),
        std::type_index(typeid(SpatialIndexResult))
    };
}


// ============================================================================
// MopacCoulombResult::Compute
//
// Same kernel as CoulombResult — same dipolar EFG, same Coulomb constant,
// same decomposition — but reading mopac_charge (PM7 Mulliken, conformation-
// dependent) instead of partial_charge (ff14SB, fixed per atom type).
//
// E_a(i) = ke * sum_{j!=i} q_mopac_j * (r_i - r_j)_a / |r_i - r_j|^3
// V_ab(i) = ke * sum_{j!=i} q_mopac_j * [3(r_i-r_j)_a(r_i-r_j)_b/|r_i-r_j|^5
//                                         - delta_ab / |r_i-r_j|^3]
// ============================================================================

std::unique_ptr<MopacCoulombResult> MopacCoulombResult::Compute(
        ProteinConformation& conf) {

    OperationLog::Scope scope("MopacCoulombResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()));

    const Protein& protein = conf.ProteinRef();
    const size_t n_atoms = conf.AtomCount();

    auto result_ptr = std::make_unique<MopacCoulombResult>();
    result_ptr->conf_ = &conf;

    // ------------------------------------------------------------------
    // Atom classification: backbone, aromatic, sidechain.
    // Same topology walk as CoulombResult — from Residue backbone cache
    // and Ring atom membership. No EnrichmentResult dependency.
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
    // Primary bond direction for E_bond_proj (same as CoulombResult).
    // ------------------------------------------------------------------

    std::vector<Vec3> primary_bond_dir(n_atoms, Vec3::Zero());
    for (size_t ai = 0; ai < n_atoms; ++ai) {
        const Atom& atom = protein.AtomAt(ai);
        if (atom.element == Element::H && atom.parent_atom_index != SIZE_MAX) {
            Vec3 d = conf.PositionAt(ai) - conf.PositionAt(atom.parent_atom_index);
            double len = d.norm();
            if (len > NEAR_ZERO_NORM) primary_bond_dir[ai] = d / len;
        } else if (!atom.bond_indices.empty()) {
            const Bond& b = protein.BondAt(atom.bond_indices[0]);
            size_t other = (b.atom_index_a == ai) ? b.atom_index_b : b.atom_index_a;
            Vec3 d = conf.PositionAt(other) - conf.PositionAt(ai);
            double len = d.norm();
            if (len > NEAR_ZERO_NORM) primary_bond_dir[ai] = d / len;
        }
    }

    // ------------------------------------------------------------------
    // Main N^2 Coulomb sum with MOPAC charges
    // ------------------------------------------------------------------

    KernelFilterSet filters;
    filters.Add(std::make_unique<SelfSourceFilter>());

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

        for (size_t j = 0; j < n_atoms; ++j) {
            KernelEvaluationContext ctx;
            ctx.atom_index = i;
            ctx.source_atom_a = j;
            ctx.distance = (pos_i - conf.PositionAt(j)).norm();
            if (!filters.AcceptAll(ctx)) continue;

            // MOPAC QM charge instead of ff14SB fixed charge
            double q_j = conf.AtomAt(j).mopac_charge;
            if (std::abs(q_j) < 1e-15) continue;

            Vec3 r = pos_i - conf.PositionAt(j);
            double r_mag = r.norm();

            if (r_mag < MIN_DISTANCE) continue;

            double r3 = r_mag * r_mag * r_mag;
            double r5 = r3 * r_mag * r_mag;

            // E_a = q_j * r_a / r^3
            Vec3 E_j = q_j * r / r3;

            // V_ab = q_j * (3 r_a r_b / r^5 - delta_ab / r^3)
            Mat3 V_j = q_j * (3.0 * r * r.transpose() / r5
                              - Mat3::Identity() / r3);

            E_total += E_j;
            EFG_total += V_j;

            if (is_aromatic_atom[j]) {
                E_aromatic += E_j;
                EFG_aromatic += V_j;
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

        // Traceless projection: each term is traceless by Gauss's law,
        // but floating-point accumulation breaks this.
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
            E_total     *= scale;
            E_backbone  *= scale;
            E_sidechain *= scale;
            E_aromatic  *= scale;
        }

        // ------------------------------------------------------------------
        // Store on ConformationAtom
        // ------------------------------------------------------------------
        auto& ca = conf.MutableAtomAt(i);

        ca.mopac_coulomb_E_total     = E_total;
        ca.mopac_coulomb_E_backbone  = E_backbone;
        ca.mopac_coulomb_E_sidechain = E_sidechain;
        ca.mopac_coulomb_E_aromatic  = E_aromatic;

        ca.mopac_coulomb_EFG_total   = EFG_total;
        ca.mopac_coulomb_EFG_total_spherical = SphericalTensor::Decompose(EFG_total);

        ca.mopac_coulomb_EFG_backbone = EFG_backbone;
        ca.mopac_coulomb_EFG_backbone_spherical = SphericalTensor::Decompose(EFG_backbone);

        ca.mopac_coulomb_EFG_aromatic = EFG_aromatic;
        ca.mopac_coulomb_EFG_aromatic_spherical = SphericalTensor::Decompose(EFG_aromatic);

        ca.mopac_coulomb_E_magnitude = E_total.norm();

        ca.mopac_coulomb_E_bond_proj = E_total.dot(primary_bond_dir[i]);

        if (ca.mopac_coulomb_E_magnitude > NEAR_ZERO_NORM) {
            Vec3 E_hat = E_total / ca.mopac_coulomb_E_magnitude;
            ca.mopac_coulomb_E_backbone_frac = E_backbone.dot(E_hat);
        } else {
            ca.mopac_coulomb_E_backbone_frac = 0.0;
        }

        // Shielding contribution: the total EFG SphericalTensor.
        // Pure T2 (EFG is traceless). gamma converts this to shielding.
        ca.mopac_coulomb_shielding_contribution =
            SphericalTensor::Decompose(EFG_total);
    }

    OperationLog::Info(LogCalcOther, "MopacCoulombResult::Compute",
        "atoms=" + std::to_string(n_atoms) +
        " rejected={" + filters.ReportRejections() + "}");

    return result_ptr;
}


// ============================================================================
// Query methods
// ============================================================================

Vec3 MopacCoulombResult::EFieldAt(size_t atom_index) const {
    return conf_->AtomAt(atom_index).mopac_coulomb_E_total;
}

Mat3 MopacCoulombResult::EFGAt(size_t atom_index) const {
    return conf_->AtomAt(atom_index).mopac_coulomb_EFG_total;
}

SphericalTensor MopacCoulombResult::EFGSphericalAt(size_t atom_index) const {
    return conf_->AtomAt(atom_index).mopac_coulomb_EFG_total_spherical;
}


// ============================================================================
// WriteFeatures: mopac_coulomb_shielding (9), E-field (3),
// EFG decompositions, scalar features.
// ============================================================================

static void PackST_MCC(const SphericalTensor& st, double* out) {
    out[0] = st.T0;
    for (int i = 0; i < 3; ++i) out[1+i] = st.T1[i];
    for (int i = 0; i < 5; ++i) out[4+i] = st.T2[i];
}

int MopacCoulombResult::WriteFeatures(const ProteinConformation& conf,
                                       const std::string& output_dir) const {
    const size_t N = conf.AtomCount();

    std::vector<double> shielding(N * 9);
    std::vector<double> efield(N * 3);
    std::vector<double> efg_bb(N * 9);
    std::vector<double> efg_aro(N * 9);
    std::vector<double> scalars(N * 4);

    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        PackST_MCC(ca.mopac_coulomb_shielding_contribution, &shielding[i*9]);

        efield[i*3+0] = ca.mopac_coulomb_E_total.x();
        efield[i*3+1] = ca.mopac_coulomb_E_total.y();
        efield[i*3+2] = ca.mopac_coulomb_E_total.z();

        PackST_MCC(ca.mopac_coulomb_EFG_backbone_spherical, &efg_bb[i*9]);
        PackST_MCC(ca.mopac_coulomb_EFG_aromatic_spherical, &efg_aro[i*9]);

        scalars[i*4+0] = ca.mopac_coulomb_E_magnitude;
        scalars[i*4+1] = ca.mopac_coulomb_E_bond_proj;
        scalars[i*4+2] = ca.mopac_coulomb_E_backbone_frac;
        scalars[i*4+3] = ca.mopac_coulomb_E_aromatic.norm();
    }

    NpyWriter::WriteFloat64(output_dir + "/mopac_coulomb_shielding.npy",
                            shielding.data(), N, 9);
    NpyWriter::WriteFloat64(output_dir + "/mopac_coulomb_E.npy",
                            efield.data(), N, 3);
    NpyWriter::WriteFloat64(output_dir + "/mopac_coulomb_efg_backbone.npy",
                            efg_bb.data(), N, 9);
    NpyWriter::WriteFloat64(output_dir + "/mopac_coulomb_efg_aromatic.npy",
                            efg_aro.data(), N, 9);
    NpyWriter::WriteFloat64(output_dir + "/mopac_coulomb_scalars.npy",
                            scalars.data(), N, 4);
    return 5;
}

}  // namespace nmr

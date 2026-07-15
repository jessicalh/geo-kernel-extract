#include "MopacCoulombResult.h"
#include "Protein.h"
#include "MopacResult.h"
#include "SpatialIndexResult.h"
#include "KernelEvaluationFilter.h"
#include "PhysicalConstants.h"
#include "CalculatorConfig.h"
#include "GeometryChoice.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <cmath>
#include <limits>
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
// same decomposition — but reading mopac_charge (the legacy F15.6
// projection of the PM7 Coulson charge, conformation-dependent) instead of
// force-field partial_charge.
//
// E_a(i) = ke * sum_{j!=i} q_mopac_j * (r_i - r_j)_a / |r_i - r_j|^3
// V_ab(i) = ke * sum_{j!=i} q_mopac_j * [3(r_i-r_j)_a(r_i-r_j)_b/|r_i-r_j|^5
//                                         - delta_ab / |r_i-r_j|^3]
// ============================================================================

std::unique_ptr<MopacCoulombResult> MopacCoulombResult::Compute(
        ProteinConformation& conf) {

    // Near-verbatim copy of CoulombResult::Compute. Differs only in:
    //   - charge source: mopac_charge (legacy F15.6 PM7 Coulson projection)
    //                    vs force-field partial_charge
    //   - source set: all pairs (full N^2) vs AtomsWithinRadius cutoff
    //   - no APBS solvent subtraction
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
    std::vector<bool> has_primary_bond_dir(n_atoms, false);
    for (size_t ai = 0; ai < n_atoms; ++ai) {
        const Atom& atom = protein.AtomAt(ai);
        if (atom.element == Element::H &&
            atom.parent_atom_index != SIZE_MAX &&
            atom.parent_atom_index < n_atoms) {
            Vec3 d = conf.PositionAt(ai) - conf.PositionAt(atom.parent_atom_index);
            double len = d.norm();
            if (len > CalculatorConfig::Get("near_zero_vector_norm_threshold")) {
                primary_bond_dir[ai] = d / len;
                has_primary_bond_dir[ai] = true;
            }
        }
    }

    // ------------------------------------------------------------------
    // Main N^2 Coulomb sum with MOPAC charges
    // ------------------------------------------------------------------

    KernelFilterSet filters;
    filters.Add(std::make_unique<MinDistanceFilter>());
    filters.Add(std::make_unique<SelfSourceFilter>());

    GeometryChoiceBuilder choices(conf);

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

        int charge_floor_skipped = 0;
        int n_sidechain_aromatic_sources = 0;

        for (size_t j = 0; j < n_atoms; ++j) {
            KernelEvaluationContext ctx;
            ctx.atom_index = i;
            ctx.source_atom_a = j;
            ctx.distance = (pos_i - conf.PositionAt(j)).norm();
            if (!filters.AcceptAll(ctx)) continue;

            // MOPAC QM charge instead of force-field charge
            double q_j = conf.AtomAt(j).mopac_charge;
            if (std::abs(q_j) < CalculatorConfig::Get("coulomb_charge_noise_floor")) { charge_floor_skipped++; continue; }

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

            if (is_aromatic_atom[j]) {
                E_aromatic += E_j;
                EFG_aromatic += V_j;
                if (!is_backbone[j]) ++n_sidechain_aromatic_sources;
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
        const double E_mag = E_total.norm();
        std::uint8_t clamp_mask = 0;
        double clamp_scale = 1.0;
        if (E_mag > CalculatorConfig::Get("efield_magnitude_sanity_clamp")) {
            clamp_mask = 1;
            clamp_scale =
                CalculatorConfig::Get("efield_magnitude_sanity_clamp") / E_mag;

            // ---- GeometryChoice: E-field clamp ----
            choices.Record(CalculatorId::MopacCoulomb, i, "mopac E-field clamp",
                [&conf, i, E_mag, clamp_scale](GeometryChoice& gc) {
                    AddAtom(gc, &conf.AtomAt(i), i, EntityRole::Target, EntityOutcome::Triggered);
                    AddNumber(gc, "actual_E_magnitude", E_mag, "V/A");
                    AddNumber(gc, "scale_factor", clamp_scale, "");
                });

            E_total     *= clamp_scale;
            E_backbone  *= clamp_scale;
            E_sidechain *= clamp_scale;
            E_aromatic  *= clamp_scale;
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

        ca.mopac_coulomb_EFG_sidechain = EFG_sidechain;
        ca.mopac_coulomb_EFG_sidechain_spherical = SphericalTensor::Decompose(EFG_sidechain);

        ca.mopac_coulomb_EFG_aromatic = EFG_aromatic;
        ca.mopac_coulomb_EFG_aromatic_spherical = SphericalTensor::Decompose(EFG_aromatic);

        ca.mopac_coulomb_E_magnitude = E_total.norm();
        ca.mopac_coulomb_aromatic_E_magnitude = E_aromatic.norm();
        ca.mopac_coulomb_aromatic_E_bond_proj = has_primary_bond_dir[i]
            ? E_aromatic.dot(primary_bond_dir[i])
            : std::numeric_limits<double>::quiet_NaN();
        ca.mopac_coulomb_aromatic_n_sidechain_atoms =
            n_sidechain_aromatic_sources;
        ca.mopac_coulomb_efield_clamp_mask = clamp_mask;
        ca.mopac_coulomb_efield_clamp_scale = clamp_scale;

        ca.mopac_coulomb_E_bond_proj = has_primary_bond_dir[i]
            ? E_total.dot(primary_bond_dir[i])
            : std::numeric_limits<double>::quiet_NaN();

        if (ca.mopac_coulomb_E_magnitude > CalculatorConfig::Get("near_zero_vector_norm_threshold")) {
            Vec3 E_hat = E_total / ca.mopac_coulomb_E_magnitude;
            ca.mopac_coulomb_E_backbone_frac = E_backbone.dot(E_hat);
        } else {
            ca.mopac_coulomb_E_backbone_frac = 0.0;
        }

        // Shielding contribution: the total EFG SphericalTensor.
        // Pure T2 (EFG is traceless). gamma converts this to shielding.
        ca.mopac_coulomb_shielding_contribution =
            SphericalTensor::Decompose(EFG_total);

        // ---- GeometryChoice: charge noise floor count ----
        if (charge_floor_skipped > 0) {
            choices.Record(CalculatorId::MopacCoulomb, i, "mopac charge floor",
                [&ca, i, charge_floor_skipped](GeometryChoice& gc) {
                    AddAtom(gc, &ca, i, EntityRole::Target, EntityOutcome::Included);
                    AddNumber(gc, "zero_charge_skipped", static_cast<double>(charge_floor_skipped), "count");
                });
        }
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
// WriteFeatures: mopac_coulomb_efg (9), E-field (3),
// EFG decompositions, scalar features, aromatic projection/source count,
// and clamp provenance.
// ============================================================================

int MopacCoulombResult::WriteFeatures(const ProteinConformation& conf,
                                       const std::string& output_dir) const {
    const size_t N = conf.AtomCount();

    std::vector<double> efg_total(N * 9);
    std::vector<double> efield(N * 3);
    std::vector<double> efield_bb(N * 3);
    std::vector<double> efield_sc(N * 3);
    std::vector<double> efield_aro(N * 3);
    // EFG schema rev 2026-05-18: T2 only (5 components). Same symmetric
    // outer-product physics as Coulomb → T0+T1 structural zeros.
    std::vector<double> efg_bb(N * 5);
    std::vector<double> efg_sc(N * 5);
    std::vector<double> efg_aro(N * 5);
    std::vector<double> scalars(N * 4);
    std::vector<double> aromatic_E_proj(N);
    std::vector<int32_t> aromatic_n_src(N);
    std::vector<std::uint8_t> clamp_mask(N);
    std::vector<double> clamp_scale(N);

    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        ca.mopac_coulomb_shielding_contribution.PackFull9(&efg_total[i*9]);

        efield[i*3+0] = ca.mopac_coulomb_E_total.x();
        efield[i*3+1] = ca.mopac_coulomb_E_total.y();
        efield[i*3+2] = ca.mopac_coulomb_E_total.z();

        efield_bb[i*3+0] = ca.mopac_coulomb_E_backbone.x();
        efield_bb[i*3+1] = ca.mopac_coulomb_E_backbone.y();
        efield_bb[i*3+2] = ca.mopac_coulomb_E_backbone.z();

        efield_sc[i*3+0] = ca.mopac_coulomb_E_sidechain.x();
        efield_sc[i*3+1] = ca.mopac_coulomb_E_sidechain.y();
        efield_sc[i*3+2] = ca.mopac_coulomb_E_sidechain.z();

        efield_aro[i*3+0] = ca.mopac_coulomb_E_aromatic.x();
        efield_aro[i*3+1] = ca.mopac_coulomb_E_aromatic.y();
        efield_aro[i*3+2] = ca.mopac_coulomb_E_aromatic.z();

        ca.mopac_coulomb_EFG_backbone_spherical.PackT2(&efg_bb[i*5]);
        ca.mopac_coulomb_EFG_sidechain_spherical.PackT2(&efg_sc[i*5]);
        ca.mopac_coulomb_EFG_aromatic_spherical.PackT2(&efg_aro[i*5]);

        scalars[i*4+0] = ca.mopac_coulomb_E_magnitude;
        scalars[i*4+1] = ca.mopac_coulomb_E_bond_proj;
        scalars[i*4+2] = ca.mopac_coulomb_E_backbone_frac;
        scalars[i*4+3] = ca.mopac_coulomb_aromatic_E_magnitude;
        aromatic_E_proj[i] = ca.mopac_coulomb_aromatic_E_bond_proj;
        aromatic_n_src[i] = static_cast<int32_t>(
            ca.mopac_coulomb_aromatic_n_sidechain_atoms);
        clamp_mask[i] = ca.mopac_coulomb_efield_clamp_mask;
        clamp_scale[i] = ca.mopac_coulomb_efield_clamp_scale;
    }

    int files_written = 0;
    auto record_write = [&](bool success, const std::string& filename) {
        if (success) {
            ++files_written;
            return;
        }
        OperationLog::Error(
            "MopacCoulombResult::WriteFeatures",
            "failed to write " + output_dir + "/" + filename);
    };

    record_write(
        NpyWriter::WriteFloat64(output_dir + "/mopac_coulomb_efg.npy",
                                efg_total.data(), N, 9),
        "mopac_coulomb_efg.npy");
    record_write(
        NpyWriter::WriteFloat64(output_dir + "/mopac_coulomb_E.npy",
                                efield.data(), N, 3),
        "mopac_coulomb_E.npy");
    record_write(
        NpyWriter::WriteFloat64(output_dir + "/mopac_coulomb_E_backbone.npy",
                                efield_bb.data(), N, 3),
        "mopac_coulomb_E_backbone.npy");
    record_write(
        NpyWriter::WriteFloat64(output_dir + "/mopac_coulomb_E_sidechain.npy",
                                efield_sc.data(), N, 3),
        "mopac_coulomb_E_sidechain.npy");
    record_write(
        NpyWriter::WriteFloat64(output_dir + "/mopac_coulomb_E_aromatic.npy",
                                efield_aro.data(), N, 3),
        "mopac_coulomb_E_aromatic.npy");
    record_write(
        NpyWriter::WriteFloat64(output_dir + "/mopac_coulomb_efg_backbone.npy",
                                efg_bb.data(), N, 5),
        "mopac_coulomb_efg_backbone.npy");
    record_write(
        NpyWriter::WriteFloat64(output_dir + "/mopac_coulomb_efg_sidechain.npy",
                                efg_sc.data(), N, 5),
        "mopac_coulomb_efg_sidechain.npy");
    record_write(
        NpyWriter::WriteFloat64(output_dir + "/mopac_coulomb_efg_aromatic.npy",
                                efg_aro.data(), N, 5),
        "mopac_coulomb_efg_aromatic.npy");
    record_write(
        NpyWriter::WriteFloat64(output_dir + "/mopac_coulomb_scalars.npy",
                                scalars.data(), N, 4),
        "mopac_coulomb_scalars.npy");
    record_write(
        NpyWriter::WriteFloat64(
            output_dir + "/mopac_coulomb_aromatic_E_proj.npy",
            aromatic_E_proj.data(), N),
        "mopac_coulomb_aromatic_E_proj.npy");
    record_write(
        NpyWriter::WriteInt32(
            output_dir + "/mopac_coulomb_aromatic_n_src.npy",
            aromatic_n_src.data(), N),
        "mopac_coulomb_aromatic_n_src.npy");
    record_write(
        NpyWriter::WriteUInt8(
            output_dir + "/mopac_coulomb_E_clamp_mask.npy",
            clamp_mask.data(), N),
        "mopac_coulomb_E_clamp_mask.npy");
    record_write(
        NpyWriter::WriteFloat64(
            output_dir + "/mopac_coulomb_E_clamp_scale.npy",
            clamp_scale.data(), N),
        "mopac_coulomb_E_clamp_scale.npy");
    return files_written;
}

}  // namespace nmr

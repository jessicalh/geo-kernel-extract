#include "EeqCoulombResult.h"

#include "CalculatorConfig.h"
#include "EeqResult.h"
#include "GeometryChoice.h"
#include "KernelEvaluationFilter.h"
#include "NpyWriter.h"
#include "OperationLog.h"
#include "PhysicalConstants.h"
#include "Protein.h"
#include "SpatialIndexResult.h"

#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace nmr {

std::vector<std::type_index> EeqCoulombResult::Dependencies() const {
    return {
        std::type_index(typeid(EeqResult)),
        std::type_index(typeid(SpatialIndexResult))
    };
}

std::unique_ptr<EeqCoulombResult> EeqCoulombResult::Compute(
        ProteinConformation& conf) {
    OperationLog::Scope scope("EeqCoulombResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()));

    const Protein& protein = conf.ProteinRef();
    const size_t n_atoms = conf.AtomCount();
    auto result = std::make_unique<EeqCoulombResult>();
    result->conf_ = &conf;

    // EEQ-Coulomb is a new calculator, so it does not inherit the legacy
    // FF-Coulomb nonfinite-to-zero sanitizer. Validate every source before
    // writing any owned field: an invalid charge/coordinate makes the result
    // absent and is logged by the normal calculator failure path.
    for (size_t i = 0; i < n_atoms; ++i) {
        if (!conf.PositionAt(i).allFinite()) {
            OperationLog::Error("EeqCoulombResult::Compute",
                "non-finite position at atom " + std::to_string(i) +
                "; result not attached");
            return nullptr;
        }
        if (!std::isfinite(conf.AtomAt(i).eeq_charge)) {
            OperationLog::Error("EeqCoulombResult::Compute",
                "non-finite EEQ source charge at atom " +
                std::to_string(i) + "; result not attached");
            return nullptr;
        }
    }

    std::vector<bool> is_backbone(n_atoms, false);
    std::vector<bool> is_aromatic(n_atoms, false);
    for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
        const Residue& residue = protein.ResidueAt(ri);
        auto mark = [&](size_t index) {
            if (index != Residue::NONE && index < n_atoms)
                is_backbone[index] = true;
        };
        mark(residue.N);
        mark(residue.CA);
        mark(residue.C);
        mark(residue.O);
        mark(residue.H);
        mark(residue.HA);
        mark(residue.CB);
    }
    for (size_t ri = 0; ri < protein.RingCount(); ++ri) {
        for (size_t ai : protein.RingAt(ri).atom_indices) {
            if (ai < n_atoms) is_aromatic[ai] = true;
        }
    }

    // H-only parent-to-H projection direction, matching the corrected FF
    // Coulomb scalar contract. Unsupported targets retain a NaN sentinel.
    std::vector<Vec3> primary_bond_dir(n_atoms, Vec3::Zero());
    std::vector<bool> has_primary_bond_dir(n_atoms, false);
    for (size_t ai = 0; ai < n_atoms; ++ai) {
        const Atom& atom = protein.AtomAt(ai);
        if (atom.element != Element::H ||
            atom.parent_atom_index == SIZE_MAX ||
            atom.parent_atom_index >= n_atoms) {
            continue;
        }
        const Vec3 d = conf.PositionAt(ai) -
                       conf.PositionAt(atom.parent_atom_index);
        const double length = d.norm();
        if (length > CalculatorConfig::Get("near_zero_vector_norm_threshold")) {
            primary_bond_dir[ai] = d / length;
            has_primary_bond_dir[ai] = true;
        }
    }

    KernelFilterSet filters;
    filters.Add(std::make_unique<MinDistanceFilter>());
    filters.Add(std::make_unique<SelfSourceFilter>());
    GeometryChoiceBuilder choices(conf);
    const auto& spatial = conf.Result<SpatialIndexResult>();
    const double cutoff = CalculatorConfig::Get("coulomb_efield_cutoff");

    int aromatic_source_count = 0;
    for (bool value : is_aromatic)
        if (value) ++aromatic_source_count;

    for (size_t i = 0; i < n_atoms; ++i) {
        const Vec3 pos_i = conf.PositionAt(i);
        Vec3 E_total = Vec3::Zero();
        Vec3 E_backbone = Vec3::Zero();
        Vec3 E_sidechain = Vec3::Zero();
        Vec3 E_aromatic = Vec3::Zero();
        Mat3 EFG_total = Mat3::Zero();
        Mat3 EFG_backbone = Mat3::Zero();
        Mat3 EFG_sidechain = Mat3::Zero();
        Mat3 EFG_aromatic = Mat3::Zero();
        int aromatic_sidechain_sources = 0;

        const auto neighbours = spatial.AtomsWithinRadius(pos_i, cutoff);
        // AtomsWithinRadius includes the target itself.  Cutoff provenance
        // describes candidate *sources*, so exclude self from the within
        // count and compute the complementary beyond count over N-1 sources.
        // This is intentionally independent of the later filters/charge
        // floor: those thin the summed set, not the geometric radius set.
        int sources_within_radius = 0;
        for (size_t j : neighbours) {
            if (j != i) ++sources_within_radius;
        }
        const int sources_outside_radius =
            static_cast<int>(n_atoms) - 1 - sources_within_radius;

        for (size_t j : neighbours) {
            KernelEvaluationContext ctx;
            ctx.atom_index = i;
            ctx.source_atom_a = j;
            ctx.distance = (pos_i - conf.PositionAt(j)).norm();
            if (!filters.AcceptAll(ctx)) continue;

            const double q_j = conf.AtomAt(j).eeq_charge;
            if (std::abs(q_j) <
                CalculatorConfig::Get("coulomb_charge_noise_floor")) {
                continue;
            }

            const Vec3 r = pos_i - conf.PositionAt(j);
            const double r_mag = r.norm();
            const double r3 = r_mag * r_mag * r_mag;
            const double r5 = r3 * r_mag * r_mag;
            const Vec3 E_j = q_j * r / r3;
            const Mat3 V_j = q_j *
                (3.0 * r * r.transpose() / r5 - Mat3::Identity() / r3);

            E_total += E_j;
            EFG_total += V_j;
            if (is_aromatic[j]) {
                E_aromatic += E_j;
                EFG_aromatic += V_j;
                if (!is_backbone[j]) ++aromatic_sidechain_sources;
            } else if (is_backbone[j]) {
                E_backbone += E_j;
                EFG_backbone += V_j;
            } else {
                E_sidechain += E_j;
                EFG_sidechain += V_j;
            }
        }

        E_total *= COULOMB_KE;
        E_backbone *= COULOMB_KE;
        E_sidechain *= COULOMB_KE;
        E_aromatic *= COULOMB_KE;
        EFG_total *= COULOMB_KE;
        EFG_backbone *= COULOMB_KE;
        EFG_sidechain *= COULOMB_KE;
        EFG_aromatic *= COULOMB_KE;

        auto project_traceless = [](Mat3& matrix) {
            matrix -= (matrix.trace() / 3.0) * Mat3::Identity();
        };
        project_traceless(EFG_total);
        project_traceless(EFG_backbone);
        project_traceless(EFG_sidechain);
        project_traceless(EFG_aromatic);

        if (!E_total.allFinite() || !E_backbone.allFinite() ||
            !E_sidechain.allFinite() || !E_aromatic.allFinite() ||
            !EFG_total.allFinite() || !EFG_backbone.allFinite() ||
            !EFG_sidechain.allFinite() || !EFG_aromatic.allFinite()) {
            OperationLog::Error("EeqCoulombResult::Compute",
                "non-finite Coulomb E/EFG at target atom " +
                std::to_string(i) + "; result not attached");
            return nullptr;
        }

        const double E_mag = E_total.norm();
        const double clamp =
            CalculatorConfig::Get("efield_magnitude_sanity_clamp");
        if (E_mag > clamp) {
            const double scale = clamp / E_mag;
            choices.Record(CalculatorId::EEQCoulomb, i, "EEQ E-field clamp",
                [&conf, i, E_mag, scale](GeometryChoice& gc) {
                    AddAtom(gc, &conf.AtomAt(i), i,
                            EntityRole::Target, EntityOutcome::Triggered);
                    AddNumber(gc, "actual_E_magnitude", E_mag, "V/A");
                    AddNumber(gc, "scale_factor", scale, "");
                });
            E_total *= scale;
            E_backbone *= scale;
            E_sidechain *= scale;
            E_aromatic *= scale;
        }

        if (sources_outside_radius > 0) {
            choices.Record(CalculatorId::EEQCoulomb, i, "EEQ Coulomb cutoff",
                [&conf, i, sources_within_radius, sources_outside_radius,
                 cutoff]
                (GeometryChoice& gc) {
                    AddAtom(gc, &conf.AtomAt(i), i,
                            EntityRole::Target, EntityOutcome::Included);
                    AddNumber(gc, "sources_within_cutoff",
                              static_cast<double>(sources_within_radius),
                              "count");
                    AddNumber(gc, "sources_beyond_cutoff",
                              static_cast<double>(sources_outside_radius),
                              "count");
                    AddNumber(gc, "cutoff_distance", cutoff, "A");
                });
        }

        auto& atom = conf.MutableAtomAt(i);
        atom.eeq_coulomb_E_total = E_total;
        atom.eeq_coulomb_E_backbone = E_backbone;
        atom.eeq_coulomb_E_sidechain = E_sidechain;
        atom.eeq_coulomb_E_aromatic = E_aromatic;
        atom.eeq_coulomb_EFG_total = EFG_total;
        atom.eeq_coulomb_EFG_total_spherical =
            SphericalTensor::Decompose(EFG_total);
        atom.eeq_coulomb_EFG_backbone = EFG_backbone;
        atom.eeq_coulomb_EFG_backbone_spherical =
            SphericalTensor::Decompose(EFG_backbone);
        atom.eeq_coulomb_EFG_sidechain = EFG_sidechain;
        atom.eeq_coulomb_EFG_sidechain_spherical =
            SphericalTensor::Decompose(EFG_sidechain);
        atom.eeq_coulomb_EFG_aromatic = EFG_aromatic;
        atom.eeq_coulomb_EFG_aromatic_spherical =
            SphericalTensor::Decompose(EFG_aromatic);

        atom.eeq_coulomb_E_magnitude = E_total.norm();
        atom.eeq_coulomb_E_bond_proj = has_primary_bond_dir[i]
            ? E_total.dot(primary_bond_dir[i])
            : std::numeric_limits<double>::quiet_NaN();
        if (atom.eeq_coulomb_E_magnitude >
            CalculatorConfig::Get("near_zero_vector_norm_threshold")) {
            atom.eeq_coulomb_E_backbone_frac = E_backbone.dot(
                E_total / atom.eeq_coulomb_E_magnitude);
        } else {
            atom.eeq_coulomb_E_backbone_frac = 0.0;
        }
        atom.eeq_coulomb_aromatic_E_magnitude = E_aromatic.norm();
        atom.eeq_coulomb_aromatic_E_bond_proj = has_primary_bond_dir[i]
            ? E_aromatic.dot(primary_bond_dir[i])
            : std::numeric_limits<double>::quiet_NaN();
        atom.eeq_coulomb_aromatic_n_sidechain_atoms =
            aromatic_sidechain_sources;
        atom.eeq_coulomb_shielding_contribution =
            SphericalTensor::Decompose(EFG_total);
    }

    OperationLog::Info(LogCalcOther, "EeqCoulombResult::Compute",
        "atoms=" + std::to_string(n_atoms) +
        " aromatic_sources=" + std::to_string(aromatic_source_count) +
        " rejected={" + filters.ReportRejections() + "}");
    return result;
}

Vec3 EeqCoulombResult::EFieldAt(size_t atom_index) const {
    return conf_->AtomAt(atom_index).eeq_coulomb_E_total;
}

Mat3 EeqCoulombResult::EFGAt(size_t atom_index) const {
    return conf_->AtomAt(atom_index).eeq_coulomb_EFG_total;
}

SphericalTensor EeqCoulombResult::EFGSphericalAt(size_t atom_index) const {
    return conf_->AtomAt(atom_index).eeq_coulomb_EFG_total_spherical;
}

int EeqCoulombResult::WriteFeatures(
        const ProteinConformation& conf,
        const std::string& output_dir) const {
    const size_t N = conf.AtomCount();
    std::vector<double> efg_total(N * 9);
    std::vector<double> E(N * 3);
    std::vector<double> E_backbone(N * 3);
    std::vector<double> E_sidechain(N * 3);
    std::vector<double> E_aromatic(N * 3);
    std::vector<double> efg_backbone(N * 5);
    std::vector<double> efg_sidechain(N * 5);
    std::vector<double> efg_aromatic(N * 5);
    std::vector<double> scalars(N * 4);
    std::vector<double> aromatic_projection(N);
    std::vector<int32_t> aromatic_source_count(N);

    for (size_t i = 0; i < N; ++i) {
        const auto& atom = conf.AtomAt(i);
        atom.eeq_coulomb_shielding_contribution.PackFull9(
            &efg_total[i*9]);
        for (int d = 0; d < 3; ++d) {
            E[i*3+d] = atom.eeq_coulomb_E_total(d);
            E_backbone[i*3+d] = atom.eeq_coulomb_E_backbone(d);
            E_sidechain[i*3+d] = atom.eeq_coulomb_E_sidechain(d);
            E_aromatic[i*3+d] = atom.eeq_coulomb_E_aromatic(d);
        }
        atom.eeq_coulomb_EFG_backbone_spherical.PackT2(
            &efg_backbone[i*5]);
        atom.eeq_coulomb_EFG_sidechain_spherical.PackT2(
            &efg_sidechain[i*5]);
        atom.eeq_coulomb_EFG_aromatic_spherical.PackT2(
            &efg_aromatic[i*5]);
        scalars[i*4+0] = atom.eeq_coulomb_E_magnitude;
        scalars[i*4+1] = atom.eeq_coulomb_E_bond_proj;
        scalars[i*4+2] = atom.eeq_coulomb_E_backbone_frac;
        scalars[i*4+3] = atom.eeq_coulomb_aromatic_E_magnitude;
        aromatic_projection[i] = atom.eeq_coulomb_aromatic_E_bond_proj;
        aromatic_source_count[i] = static_cast<int32_t>(
            atom.eeq_coulomb_aromatic_n_sidechain_atoms);
    }

    int files_written = 0;
    auto record_write = [&](bool success, const char* filename) {
        if (success) {
            ++files_written;
        } else {
            OperationLog::Error(
                "EeqCoulombResult::WriteFeatures",
                "failed to write " + output_dir + "/" + filename);
        }
    };

    record_write(NpyWriter::WriteFloat64(
                     output_dir + "/eeq_coulomb_efg.npy",
                     efg_total.data(), N, 9),
                 "eeq_coulomb_efg.npy");
    record_write(NpyWriter::WriteFloat64(
                     output_dir + "/eeq_coulomb_E.npy",
                     E.data(), N, 3),
                 "eeq_coulomb_E.npy");
    record_write(NpyWriter::WriteFloat64(
                     output_dir + "/eeq_coulomb_E_backbone.npy",
                     E_backbone.data(), N, 3),
                 "eeq_coulomb_E_backbone.npy");
    record_write(NpyWriter::WriteFloat64(
                     output_dir + "/eeq_coulomb_E_sidechain.npy",
                     E_sidechain.data(), N, 3),
                 "eeq_coulomb_E_sidechain.npy");
    record_write(NpyWriter::WriteFloat64(
                     output_dir + "/eeq_coulomb_E_aromatic.npy",
                     E_aromatic.data(), N, 3),
                 "eeq_coulomb_E_aromatic.npy");
    record_write(NpyWriter::WriteFloat64(
                     output_dir + "/eeq_coulomb_efg_backbone.npy",
                     efg_backbone.data(), N, 5),
                 "eeq_coulomb_efg_backbone.npy");
    record_write(NpyWriter::WriteFloat64(
                     output_dir + "/eeq_coulomb_efg_sidechain.npy",
                     efg_sidechain.data(), N, 5),
                 "eeq_coulomb_efg_sidechain.npy");
    record_write(NpyWriter::WriteFloat64(
                     output_dir + "/eeq_coulomb_efg_aromatic.npy",
                     efg_aromatic.data(), N, 5),
                 "eeq_coulomb_efg_aromatic.npy");
    record_write(NpyWriter::WriteFloat64(
                     output_dir + "/eeq_coulomb_scalars.npy",
                     scalars.data(), N, 4),
                 "eeq_coulomb_scalars.npy");
    record_write(NpyWriter::WriteFloat64(
                     output_dir + "/eeq_coulomb_aromatic_E_proj.npy",
                     aromatic_projection.data(), N),
                 "eeq_coulomb_aromatic_E_proj.npy");
    record_write(NpyWriter::WriteInt32(
                     output_dir + "/eeq_coulomb_aromatic_n_src.npy",
                     aromatic_source_count.data(), N),
                 "eeq_coulomb_aromatic_n_src.npy");
    return files_written;
}

}  // namespace nmr

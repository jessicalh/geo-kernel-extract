#include "LocalBackboneGeometryTrajectoryResult.h"

#include "LocalBackboneGeometryResult.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "TrajectoryProtein.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <algorithm>
#include <limits>
#include <string>

namespace nmr {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

bool HasAtom(std::size_t atom_index) {
    return atom_index != Residue::NONE;
}

}  // namespace


std::unique_ptr<LocalBackboneGeometryTrajectoryResult>
LocalBackboneGeometryTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto result = std::make_unique<LocalBackboneGeometryTrajectoryResult>();
    const Protein& protein = tp.ProteinRef();
    const std::size_t R = protein.ResidueCount();

    result->tau_n_ca_c_.assign(R, {});
    result->angle_n_ca_cb_.assign(R, {});
    result->angle_cb_ca_c_.assign(R, {});
    result->angle_cprev_n_ca_.assign(R, {});
    result->angle_ca_c_nnext_.assign(R, {});
    result->cb_deviation_.assign(R, {});
    result->cb_local_vector_.assign(R, {});

    result->has_n_ca_c_.assign(R, 0u);
    result->has_cb_.assign(R, 0u);
    result->has_prev_c_.assign(R, 0u);
    result->has_next_n_.assign(R, 0u);
    result->is_glycine_.assign(R, 0u);
    result->is_proline_.assign(R, 0u);

    for (std::size_t ri = 0; ri < R; ++ri) {
        const Residue& residue = protein.ResidueAt(ri);
        result->has_n_ca_c_[ri] =
            HasAtom(residue.N) && HasAtom(residue.CA) && HasAtom(residue.C) ? 1u : 0u;
        result->has_cb_[ri] = HasAtom(residue.CB) ? 1u : 0u;
        result->is_glycine_[ri] = residue.type == AminoAcid::GLY ? 1u : 0u;
        result->is_proline_[ri] = residue.type == AminoAcid::PRO ? 1u : 0u;

        // The graph query is the single adjacency authority. Its successful
        // predecessor/successor contract guarantees the corresponding C/N
        // backbone slot exists on the returned residue.
        result->has_prev_c_[ri] = protein.BackbonePredecessor(ri).has_value() ? 1u : 0u;
        result->has_next_n_[ri] = protein.BackboneSuccessor(ri).has_value() ? 1u : 0u;
    }

    return result;
}


void LocalBackboneGeometryTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp;
    (void)traj;
    const Protein& protein = conf.ProteinRef();
    const std::size_t R = protein.ResidueCount();

    for (std::size_t ri = 0; ri < R; ++ri) {
        const auto m = local_backbone_geometry::MeasureResidue(conf, ri);
        tau_n_ca_c_[ri].push_back(m.tau_n_ca_c);
        angle_n_ca_cb_[ri].push_back(m.angle_n_ca_cb);
        angle_cb_ca_c_[ri].push_back(m.angle_cb_ca_c);
        angle_cprev_n_ca_[ri].push_back(m.angle_cprev_n_ca);
        angle_ca_c_nnext_[ri].push_back(m.angle_ca_c_nnext);
        cb_deviation_[ri].push_back(m.cb_deviation);
        cb_local_vector_[ri].push_back(m.cb_residual_vector);
    }

    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    source_attached_per_frame_.push_back(1u);
    ++n_frames_;
}


void LocalBackboneGeometryTrajectoryResult::Finalize(
        TrajectoryProtein& tp,
        Trajectory& traj) {
    (void)tp;
    (void)traj;
    finalized_ = true;
    OperationLog::Info(LogCalcOther,
                       "LocalBackboneGeometryTrajectoryResult::Finalize",
                       "finalized " + std::to_string(tau_n_ca_c_.size())
                           + " residues across " + std::to_string(n_frames_) + " dispatched frames");
}


void LocalBackboneGeometryTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    (void)tp;
    const std::size_t R = tau_n_ca_c_.size();
    const std::size_t T = n_frames_;
    const std::size_t source_attached_count = static_cast<std::size_t>(
        std::count(source_attached_per_frame_.begin(), source_attached_per_frame_.end(), std::uint8_t{1}));

    auto group = file.createGroup("/trajectory/local_backbone_geometry");
    group.createAttribute("result_name", Name());
    group.createAttribute("n_residues", R);
    group.createAttribute("n_frames", T);
    group.createAttribute("source_attached_count", source_attached_count);
    group.createAttribute("finalized", finalized_);
    group.createAttribute("angle_units", std::string("radians"));
    group.createAttribute("length_units", std::string("Angstrom"));
    group.createAttribute("residue_axis", std::string("protein_residue_index"));
    group.createAttribute("source_attached_policy", std::string(
        "always_attached -- positions and invariant typed topology are present for every dispatched frame"));
    group.createAttribute("adjacency_policy", std::string(
        "Protein::BackbonePredecessor/BackboneSuccessor peptide-bond graph; never residue-index arithmetic"));
    group.createAttribute("ideal_cb_formula", std::string(
        "n=N-CA; c=C-CA; ideal_CB=CA-0.58273431*cross(c,n)+0.56802827*n-0.54067466*c"));
    group.createAttribute("cb_deviation_definition", std::string(
        "norm(observed_CB-ideal_CB); NaN if N,CA,C,CB is absent or the N-CA-C frame is degenerate"));
    group.createAttribute("cb_local_vector_definition", std::string(
        "observed_CB-ideal_CB in global Cartesian xyz Angstrom; norm equals cb_deviation; identical NaN gate"));

    group.createDataSet("frame_indices", frame_indices_)
         .createAttribute("units", std::string("frame_index"));
    group.createDataSet("frame_times", frame_times_)
         .createAttribute("units", std::string("ps"));
    group.createDataSet("source_attached_per_frame", source_attached_per_frame_)
         .createAttribute("units", std::string("dimensionless"));

    auto emit_f64_2d = [&](const std::string& name,
                           const std::vector<std::vector<double>>& values,
                           const std::string& units) {
        std::vector<double> flat(R * T, kNaN);
        for (std::size_t ri = 0; ri < R && ri < values.size(); ++ri) {
            for (std::size_t ti = 0; ti < T && ti < values[ri].size(); ++ti)
                flat[ri * T + ti] = values[ri][ti];
        }
        HighFive::DataSpace space(std::vector<std::size_t>{R, T});
        auto dataset = group.createDataSet<double>(name, space);
        dataset.write_raw(flat.data());
        dataset.createAttribute("units", units);
    };

    emit_f64_2d("tau_N_CA_C", tau_n_ca_c_, "radians");
    emit_f64_2d("angle_N_CA_CB", angle_n_ca_cb_, "radians");
    emit_f64_2d("angle_CB_CA_C", angle_cb_ca_c_, "radians");
    emit_f64_2d("angle_Cprev_N_CA", angle_cprev_n_ca_, "radians");
    emit_f64_2d("angle_CA_C_Nnext", angle_ca_c_nnext_, "radians");
    emit_f64_2d("cb_deviation", cb_deviation_, "Angstrom");
    {
        auto dataset = group.getDataSet("cb_deviation");
        dataset.createAttribute("coordinate_frame",
                                std::string("intrinsic_chiral_lookup"));
        dataset.createAttribute("transformation", std::string(
            "rotation-invariant under proper rotations; ideal-L-CB construction "
            "is chirality-conditioned and has no improper-transform contract"));
        dataset.createAttribute("parity", std::string("mixed"));
        dataset.createAttribute("validity", std::string(
            "NaN if N,CA,C,CB is absent or the N-CA-C frame is degenerate"));
    }

    {
        std::vector<double> flat(R * T * 3, kNaN);
        for (std::size_t ri = 0; ri < R && ri < cb_local_vector_.size(); ++ri) {
            for (std::size_t ti = 0; ti < T && ti < cb_local_vector_[ri].size(); ++ti) {
                const Vec3& value = cb_local_vector_[ri][ti];
                const std::size_t base = (ri * T + ti) * 3;
                flat[base + 0] = value.x();
                flat[base + 1] = value.y();
                flat[base + 2] = value.z();
            }
        }
        HighFive::DataSpace space(std::vector<std::size_t>{R, T, std::size_t{3}});
        auto dataset = group.createDataSet<double>("cb_local_vector", space);
        dataset.write_raw(flat.data());
        dataset.createAttribute("units", std::string("Angstrom"));
        dataset.createAttribute("axis_3", std::string("global_cartesian_x_y_z"));
        dataset.createAttribute("coordinate_frame",
                                std::string("conformation_cartesian_xyz"));
        dataset.createAttribute("transformation", std::string(
            "polar displacement under proper rotations: v'=R v; the ideal-L-CB "
            "construction mixes bond vectors with a cross product and therefore "
            "has no single improper-transform parity"));
        dataset.createAttribute("parity", std::string("mixed"));
        dataset.createAttribute("validity", std::string(
            "NaN if N,CA,C,CB is absent or the N-CA-C frame is degenerate"));
    }

    auto emit_mask = [&](const std::string& name,
                         const std::vector<std::uint8_t>& values,
                         const std::string& note) {
        auto dataset = group.createDataSet(name, values);
        dataset.createAttribute("units", std::string("dimensionless"));
        dataset.createAttribute("note", note);
    };
    emit_mask("has_N_CA_C", has_n_ca_c_, "1 when typed Residue N, CA, and C slots all exist");
    emit_mask("has_CB", has_cb_, "1 when the typed Residue CB slot exists; GLY normally has 0");
    emit_mask("has_prev_C", has_prev_c_,
              "1 when Protein::BackbonePredecessor finds a peptide-bond predecessor carrying C");
    emit_mask("has_next_N", has_next_n_,
              "1 when Protein::BackboneSuccessor finds a peptide-bond successor carrying N");
    emit_mask("is_glycine", is_glycine_, "1 when residue type is GLY");
    emit_mask("is_proline", is_proline_, "1 when residue type is PRO");

    OperationLog::Info(LogCalcOther,
                       "LocalBackboneGeometryTrajectoryResult::WriteH5Group",
                       "wrote /trajectory/local_backbone_geometry with " + std::to_string(R)
                           + " residues x " + std::to_string(T) + " dispatched frames");
}

}  // namespace nmr

#include "LocalBackboneGeometryTrajectoryResult.h"

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
#include <cmath>
#include <limits>
#include <string>

namespace nmr {

// Production geometry kernels have external linkage in this per-file named
// namespace so fixture-independent forcing tests execute the exact functions
// used by Compute. They remain local to this result rather than becoming a
// shared geometry utility/header family.
namespace local_backbone_geometry {

namespace {

constexpr double kNearZero = 1e-12;

Vec3 NanVec3() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    return Vec3(nan, nan, nan);
}

bool Finite(const Vec3& v) {
    return v.allFinite();
}

}  // namespace

double BondAngleRadians(const Vec3& endpoint_a,
                        const Vec3& vertex,
                        const Vec3& endpoint_b) {
    const Vec3 u = endpoint_a - vertex;
    const Vec3 v = endpoint_b - vertex;
    if (!Finite(u) || !Finite(v))
        return std::numeric_limits<double>::quiet_NaN();
    const double un = u.norm();
    const double vn = v.norm();
    if (un <= kNearZero || vn <= kNearZero)
        return std::numeric_limits<double>::quiet_NaN();
    const double cosine = std::clamp(u.dot(v) / (un * vn), -1.0, 1.0);
    return std::acos(cosine);
}

Vec3 IdealCbPosition(const Vec3& n_pos,
                     const Vec3& ca_pos,
                     const Vec3& c_pos) {
    if (!Finite(n_pos) || !Finite(ca_pos) || !Finite(c_pos))
        return NanVec3();

    const Vec3 n = n_pos - ca_pos;
    const Vec3 c = c_pos - ca_pos;
    const double nn = n.norm();
    const double cn = c.norm();
    if (nn <= kNearZero || cn <= kNearZero)
        return NanVec3();

    const Vec3 cross_cn = c.cross(n);
    // Relative collinearity gate: a zero/near-zero N-CA-C plane cannot
    // define the chiral ideal-CB construction.
    if (cross_cn.norm() <= kNearZero * nn * cn)
        return NanVec3();

    return ca_pos - 0.58273431 * cross_cn
                  + 0.56802827 * n
                  - 0.54067466 * c;
}

Vec3 CbDeviationVector(const Vec3& n_pos,
                       const Vec3& ca_pos,
                       const Vec3& c_pos,
                       const Vec3& observed_cb) {
    if (!Finite(observed_cb))
        return NanVec3();
    const Vec3 ideal = IdealCbPosition(n_pos, ca_pos, c_pos);
    if (!Finite(ideal))
        return NanVec3();
    return observed_cb - ideal;
}

double CbDeviation(const Vec3& n_pos,
                   const Vec3& ca_pos,
                   const Vec3& c_pos,
                   const Vec3& observed_cb) {
    const Vec3 residual = CbDeviationVector(n_pos, ca_pos, c_pos, observed_cb);
    return Finite(residual)
        ? residual.norm()
        : std::numeric_limits<double>::quiet_NaN();
}

}  // namespace local_backbone_geometry

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
        const Residue& residue = protein.ResidueAt(ri);
        const auto prev_idx = protein.BackbonePredecessor(ri);
        const auto next_idx = protein.BackboneSuccessor(ri);

        double tau = kNaN;
        double angle_n_ca_cb = kNaN;
        double angle_cb_ca_c = kNaN;
        double angle_cprev_n_ca = kNaN;
        double angle_ca_c_nnext = kNaN;
        double cb_deviation = kNaN;
        Vec3 cb_local_vector(kNaN, kNaN, kNaN);

        if (HasAtom(residue.N) && HasAtom(residue.CA) && HasAtom(residue.C)) {
            const Vec3& n = conf.PositionAt(residue.N);
            const Vec3& ca = conf.PositionAt(residue.CA);
            const Vec3& c = conf.PositionAt(residue.C);
            tau = local_backbone_geometry::BondAngleRadians(n, ca, c);

            if (HasAtom(residue.CB)) {
                const Vec3& cb = conf.PositionAt(residue.CB);
                angle_n_ca_cb = local_backbone_geometry::BondAngleRadians(n, ca, cb);
                angle_cb_ca_c = local_backbone_geometry::BondAngleRadians(cb, ca, c);
                cb_local_vector = local_backbone_geometry::CbDeviationVector(n, ca, c, cb);
                if (cb_local_vector.allFinite())
                    cb_deviation = cb_local_vector.norm();
            }
        } else {
            // The two CB valence angles need only their own three atoms, so
            // preserve them even when the other backbone endpoint is absent.
            if (HasAtom(residue.N) && HasAtom(residue.CA) && HasAtom(residue.CB)) {
                angle_n_ca_cb = local_backbone_geometry::BondAngleRadians(
                    conf.PositionAt(residue.N), conf.PositionAt(residue.CA), conf.PositionAt(residue.CB));
            }
            if (HasAtom(residue.CB) && HasAtom(residue.CA) && HasAtom(residue.C)) {
                angle_cb_ca_c = local_backbone_geometry::BondAngleRadians(
                    conf.PositionAt(residue.CB), conf.PositionAt(residue.CA), conf.PositionAt(residue.C));
            }
        }

        if (prev_idx && HasAtom(residue.N) && HasAtom(residue.CA)) {
            const Residue& prev = protein.ResidueAt(*prev_idx);
            angle_cprev_n_ca = local_backbone_geometry::BondAngleRadians(
                conf.PositionAt(prev.C), conf.PositionAt(residue.N), conf.PositionAt(residue.CA));
        }
        if (next_idx && HasAtom(residue.CA) && HasAtom(residue.C)) {
            const Residue& next = protein.ResidueAt(*next_idx);
            angle_ca_c_nnext = local_backbone_geometry::BondAngleRadians(
                conf.PositionAt(residue.CA), conf.PositionAt(residue.C), conf.PositionAt(next.N));
        }

        tau_n_ca_c_[ri].push_back(tau);
        angle_n_ca_cb_[ri].push_back(angle_n_ca_cb);
        angle_cb_ca_c_[ri].push_back(angle_cb_ca_c);
        angle_cprev_n_ca_[ri].push_back(angle_cprev_n_ca);
        angle_ca_c_nnext_[ri].push_back(angle_ca_c_nnext);
        cb_deviation_[ri].push_back(cb_deviation);
        cb_local_vector_[ri].push_back(cb_local_vector);
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

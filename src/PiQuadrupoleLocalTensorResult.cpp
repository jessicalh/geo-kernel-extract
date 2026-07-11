#include "PiQuadrupoleLocalTensorResult.h"

#include "BiotSavartResult.h"
#include "CalculatorConfig.h"
#include "DispersionResult.h"
#include "GeometryResult.h"
#include "HaighMallionResult.h"
#include "NpyWriter.h"
#include "OperationLog.h"
#include "PiQuadrupoleResult.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "RingSusceptibilityResult.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace nmr {

namespace pi_quadrupole_local_tensor_detail {

LocalFrame BuildLocalFrame(const RingGeometry& geom) {
    LocalFrame frame;
    const double guard =
        CalculatorConfig::Get("near_zero_vector_norm_threshold");

    const double normal_norm = geom.normal.norm();
    if (!std::isfinite(normal_norm) || normal_norm <= guard ||
        geom.vertices.empty()) {
        return frame;
    }
    frame.z_axis = geom.normal / normal_norm;

    const Vec3 reference = geom.vertices[0] - geom.center;
    const Vec3 reference_plane =
        reference - reference.dot(frame.z_axis) * frame.z_axis;
    const double reference_norm = reference_plane.norm();
    if (!std::isfinite(reference_norm) || reference_norm <= guard) {
        return frame;
    }
    frame.x_axis = reference_plane / reference_norm;

    const Vec3 y = frame.z_axis.cross(frame.x_axis);
    const double y_norm = y.norm();
    if (!std::isfinite(y_norm) || y_norm <= guard) return frame;
    frame.y_axis = y / y_norm;
    frame.valid = true;
    return frame;
}


TensorEvaluation ComputeLocalTensor(const Vec3& d_local) {
    TensorEvaluation evaluation;
    const double r = d_local.norm();
    if (!std::isfinite(r) ||
        r <= CalculatorConfig::Get("singularity_guard_distance")) {
        return evaluation;
    }

    const Vec3 n_local(0.0, 0.0, 1.0);
    const double r2 = r * r;
    const double r5 = r2 * r2 * r;
    const double r7 = r5 * r2;
    const double r9 = r7 * r2;
    const double dn = d_local.dot(n_local);
    const double dn2 = dn * dn;
    const double diagonal = 3.0 / r5 - 15.0 * dn2 / r7;

    Mat3 G;
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            G(a, b) =
                105.0 * dn2 * d_local(a) * d_local(b) / r9
                - 30.0 * dn *
                    (n_local(a) * d_local(b) +
                     n_local(b) * d_local(a)) / r7
                - 15.0 * d_local(a) * d_local(b) / r7
                + 6.0 * n_local(a) * n_local(b) / r5
                + (a == b ? diagonal : 0.0);
        }
    }

    evaluation.tensor = G;
    evaluation.valid = true;
    return evaluation;
}


static void SetMissingPayload(RingNeighbourhood& neighbour) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    neighbour.piquad_local_tensor = Mat3::Constant(nan);
    neighbour.piquad_local_spherical.T0 = nan;
    neighbour.piquad_local_spherical.T1.fill(nan);
    neighbour.piquad_local_spherical.T2.fill(nan);
    neighbour.piquad_local_x_axis = Vec3::Constant(nan);
    neighbour.piquad_local_y_axis = Vec3::Constant(nan);
    neighbour.piquad_local_z_axis = Vec3::Constant(nan);
    neighbour.piquad_local_distance = nan;
    neighbour.piquad_local_cos_theta = nan;
    neighbour.piquad_local_evaluated = true;
    neighbour.piquad_local_valid = false;
}

}  // namespace pi_quadrupole_local_tensor_detail


std::vector<std::type_index>
PiQuadrupoleLocalTensorResult::Dependencies() const {
    return {
        std::type_index(typeid(GeometryResult)),
        std::type_index(typeid(BiotSavartResult)),
        std::type_index(typeid(HaighMallionResult)),
        std::type_index(typeid(RingSusceptibilityResult)),
        std::type_index(typeid(PiQuadrupoleResult)),
        std::type_index(typeid(DispersionResult)),
    };
}


std::unique_ptr<PiQuadrupoleLocalTensorResult>
PiQuadrupoleLocalTensorResult::Compute(ProteinConformation& conf) {
    OperationLog::Scope scope("PiQuadrupoleLocalTensorResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()));

    const Protein& protein = conf.ProteinRef();
    const double norm_guard =
        CalculatorConfig::Get("near_zero_vector_norm_threshold");
    size_t evaluated_rows = 0;
    size_t invalid_frames = 0;
    size_t singular_tensors = 0;

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const Vec3 atom_pos = conf.PositionAt(ai);
        for (RingNeighbourhood& neighbour :
             conf.MutableAtomAt(ai).ring_neighbours) {
            pi_quadrupole_local_tensor_detail::SetMissingPayload(neighbour);
            ++evaluated_rows;

            if (neighbour.ring_index >= protein.RingCount() ||
                neighbour.ring_index >= conf.ring_geometries.size()) {
                const std::string message =
                    "ring-neighbour row has out-of-range ring_index=" +
                    std::to_string(neighbour.ring_index) +
                    " for atom_index=" + std::to_string(ai) +
                    " (RingCount=" + std::to_string(protein.RingCount()) +
                    ", ring_geometries=" +
                    std::to_string(conf.ring_geometries.size()) + ")";
                OperationLog::Error(
                    "PiQuadrupoleLocalTensorResult::Compute", message);
                throw std::runtime_error(message);
            }

            const RingGeometry& geom =
                conf.ring_geometries[neighbour.ring_index];
            const Vec3 d_global = atom_pos - geom.center;
            const double distance = d_global.norm();
            neighbour.piquad_local_distance = distance;

            const double normal_norm = geom.normal.norm();
            if (std::isfinite(distance) && distance > norm_guard &&
                std::isfinite(normal_norm) && normal_norm > norm_guard) {
                neighbour.piquad_local_cos_theta =
                    d_global.dot(geom.normal / normal_norm) / distance;
            }

            const auto frame =
                pi_quadrupole_local_tensor_detail::BuildLocalFrame(geom);
            if (!frame.valid) {
                ++invalid_frames;
                continue;
            }

            const Vec3 d_local(
                d_global.dot(frame.x_axis),
                d_global.dot(frame.y_axis),
                d_global.dot(frame.z_axis));
            const auto tensor =
                pi_quadrupole_local_tensor_detail::ComputeLocalTensor(
                    d_local);
            if (!tensor.valid) {
                ++singular_tensors;
                continue;
            }

            neighbour.piquad_local_tensor = tensor.tensor;
            neighbour.piquad_local_spherical =
                SphericalTensor::Decompose(tensor.tensor);
            neighbour.piquad_local_x_axis = frame.x_axis;
            neighbour.piquad_local_y_axis = frame.y_axis;
            neighbour.piquad_local_z_axis = frame.z_axis;
            neighbour.piquad_local_valid = true;
        }
    }

    if (invalid_frames > 0) {
        OperationLog::Warn(
            "PiQuadrupoleLocalTensorResult::Compute",
            "invalid local ring frames=" + std::to_string(invalid_frames) +
            "; rows retained with NaN payload and frame_valid=0");
    }
    if (singular_tensors > 0) {
        OperationLog::Warn(
            "PiQuadrupoleLocalTensorResult::Compute",
            "singular/nonfinite local tensors=" +
            std::to_string(singular_tensors) +
            "; rows retained with NaN payload and frame_valid=0");
    }
    OperationLog::Info(
        LogCalcOther, "PiQuadrupoleLocalTensorResult::Compute",
        "evaluated_rows=" + std::to_string(evaluated_rows));

    return std::make_unique<PiQuadrupoleLocalTensorResult>();
}


int PiQuadrupoleLocalTensorResult::WriteFeatures(
        const ProteinConformation& conf,
        const std::string& output_dir) const {
    const size_t N = conf.AtomCount();
    size_t P = 0;
    for (size_t ai = 0; ai < N; ++ai) {
        P += conf.AtomAt(ai).ring_neighbours.size();
    }

    std::vector<double> full_tensor(P * 9);
    std::vector<double> t2_tensor(P * 5);
    std::vector<double> local_frame(P * 9);
    std::vector<double> local_geometry(P * 8);

    size_t row = 0;
    for (size_t ai = 0; ai < N; ++ai) {
        for (const RingNeighbourhood& neighbour :
             conf.AtomAt(ai).ring_neighbours) {
            if (!neighbour.piquad_local_evaluated) {
                const std::string message =
                    "unevaluated ring-neighbour row at atom_index=" +
                    std::to_string(ai) + " ring_index=" +
                    std::to_string(neighbour.ring_index) +
                    "; a row was appended after Compute or the result was "
                    "not computed";
                OperationLog::Error(
                    "PiQuadrupoleLocalTensorResult::WriteFeatures", message);
                throw std::runtime_error(message);
            }

            neighbour.piquad_local_spherical.PackFull9(
                &full_tensor[row * 9]);
            neighbour.piquad_local_spherical.PackT2(
                &t2_tensor[row * 5]);

            double* frame_out = &local_frame[row * 9];
            frame_out[0] = neighbour.piquad_local_x_axis.x();
            frame_out[1] = neighbour.piquad_local_x_axis.y();
            frame_out[2] = neighbour.piquad_local_x_axis.z();
            frame_out[3] = neighbour.piquad_local_y_axis.x();
            frame_out[4] = neighbour.piquad_local_y_axis.y();
            frame_out[5] = neighbour.piquad_local_y_axis.z();
            frame_out[6] = neighbour.piquad_local_z_axis.x();
            frame_out[7] = neighbour.piquad_local_z_axis.y();
            frame_out[8] = neighbour.piquad_local_z_axis.z();

            double* geometry = &local_geometry[row * 8];
            geometry[0] = static_cast<double>(ai);
            geometry[1] = static_cast<double>(neighbour.ring_index);
            geometry[2] = static_cast<double>(neighbour.ring_type);
            geometry[3] = neighbour.piquad_local_distance;
            geometry[4] = neighbour.piquad_local_cos_theta;
            geometry[5] = neighbour.quad_scalar;
            geometry[6] = neighbour.piquad_local_valid ? 1.0 : 0.0;
            geometry[7] = 1.0;  // Protein::RingAt is aromatic-only.
            ++row;
        }
    }

    NpyWriter::WriteFloat64(
        output_dir + "/piquad_local_tensor.npy",
        full_tensor.data(), P, 9);
    NpyWriter::WriteFloat64(
        output_dir + "/piquad_local_T2.npy",
        t2_tensor.data(), P, 5);
    NpyWriter::WriteFloat64(
        output_dir + "/piquad_local_frame.npy",
        local_frame.data(), P, 9);
    NpyWriter::WriteFloat64(
        output_dir + "/piquad_local_geometry.npy",
        local_geometry.data(), P, 8);
    return 4;
}

}  // namespace nmr

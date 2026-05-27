#include "AIMNet2ChargeResponseGradientTimeSeriesTrajectoryResult.h"

#include "AIMNet2ChargeResponseGradientResult.h"
#include "ConformationAtom.h"
#include "OperationLog.h"
#include "ProteinConformation.h"
#include "TrajectoryProtein.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <algorithm>
#include <cstddef>
#include <limits>
#include <string>

namespace nmr {

std::unique_ptr<AIMNet2ChargeResponseGradientTimeSeriesTrajectoryResult>
AIMNet2ChargeResponseGradientTimeSeriesTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<AIMNet2ChargeResponseGradientTimeSeriesTrajectoryResult>();
    r->per_atom_vector_.assign(tp.AtomCount(), {});
    r->per_atom_scalar_.assign(tp.AtomCount(), {});
    return r;
}

void AIMNet2ChargeResponseGradientTimeSeriesTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;
    const std::size_t N = per_atom_vector_.size();
    // Source-attached gate. Always-attached policy means the source
    // should be present every frame in the production config;
    // RunConfiguration requires AIMNet2ChargeResponseGradientResult there.
    // A custom config that omits
    // the Require would land NaN-fill rows here; the F1 non-finite
    // gradient guard at the source calculator can also trip this gate
    // by returning nullptr on degenerate-backward frames (codex F1
    // 2026-05-20).
    const bool source_present = conf.HasResult<AIMNet2ChargeResponseGradientResult>();
    if (!source_present) {
        OperationLog::Warn(
            "AIMNet2ChargeResponseGradientTimeSeriesTrajectoryResult::Compute",
            "AIMNet2ChargeResponseGradientResult not attached at frame " +
            std::to_string(frame_idx) +
            " — NaN-fill emitted + mask=0. Always-attached policy "
            "requires AIMNet2ChargeResponseGradientResult; the run's "
            "PerFrameExtractionSet must RequireConformationResult"
            "(AIMNet2ChargeResponseGradientResult).");
    }
    // Per `feedback_capture_at_the_boundary` "absent, not faked":
    // NaN-fill on source-absent frames (NOT zero); mask records absence.
    const double nan_d = std::numeric_limits<double>::quiet_NaN();
    const Vec3 nan_v(nan_d, nan_d, nan_d);
    for (std::size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        per_atom_vector_[i].push_back(
            source_present ? ca.aimnet2_charge_response_gradient_vector : nan_v);
        per_atom_scalar_[i].push_back(
            source_present ? ca.aimnet2_charge_response_gradient_scalar : nan_d);
    }
    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    source_attached_per_frame_.push_back(source_present ? 1u : 0u);
    ++n_frames_;
}

void AIMNet2ChargeResponseGradientTimeSeriesTrajectoryResult::Finalize(
        TrajectoryProtein& tp, Trajectory& traj) {
    (void)tp; (void)traj;
    finalized_ = true;
    OperationLog::Info(LogCalcOther,
        "AIMNet2ChargeResponseGradientTimeSeriesTrajectoryResult::Finalize",
        "finalized across " + std::to_string(n_frames_) +
        " frames, " + std::to_string(per_atom_vector_.size()) + " atoms");
}

void AIMNet2ChargeResponseGradientTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp, HighFive::File& file) const {
    const std::size_t N = per_atom_vector_.size();
    const std::size_t T = n_frames_;

    auto grp = file.createGroup("/trajectory/aimnet2_charge_response_gradient_time_series");

    grp.createAttribute("result_name",            Name());
    grp.createAttribute("n_atoms",                N);
    grp.createAttribute("n_frames",               T);
    grp.createAttribute("finalized",              finalized_);
    grp.createAttribute("units_vector",           std::string("e^2/Å"));
    grp.createAttribute("units_scalar",           std::string("e^2/Å"));
    // Vec3 metadata follows the existing TR convention: layout +
    // normalization + parity emitted as separate attrs (codex review
    // 2026-05-20). Charge-response gradient vector is a Cartesian-ordered
    // Vec3 with odd parity (gradient of a scalar w.r.t. atomic coordinates).
    grp.createAttribute("irrep_layout_vector",    std::string("x,y,z"));
    grp.createAttribute("normalization_vector",   std::string("cartesian"));
    grp.createAttribute("parity_vector",          std::string("1o"));
    grp.createAttribute("irrep_layout_scalar",    std::string("T0"));
    grp.createAttribute("parity_scalar",          std::string("0e"));
    grp.createAttribute("source",                 std::string(
        "AIMNet2ChargeResponseGradientResult.{aimnet2_charge_response_gradient_vector (Vec3), "
        "aimnet2_charge_response_gradient_scalar (double, L2 norm of vector)}. "
        "Gradient of L = sum_j q_j^2 (AIMNet2 Hirshfeld charges, e^2) "
        "with respect to atomic coordinates (e^2/Å). Both emitted "
        "for downstream convenience per feedback_methods_accumulate."));
    grp.createAttribute("source_attached_policy", std::string("always_attached"));

    // Flat (N, T, 3) double for vector + (N, T) double for scalar.
    // Chunk (N, frame_chunk=min(T, 64), 3) / (N, frame_chunk) per movie-target.
    const std::size_t frame_chunk = std::min<std::size_t>(T, 64);

    {
        std::vector<double> flat(N * T * 3);
        for (std::size_t i = 0; i < N; ++i) {
            const auto& atom_frames = per_atom_vector_[i];
            for (std::size_t f = 0; f < T; ++f) {
                const auto& v = atom_frames[f];
                const std::size_t base = (i * T + f) * 3;
                flat[base + 0] = v.x();
                flat[base + 1] = v.y();
                flat[base + 2] = v.z();
            }
        }
        HighFive::DataSpace space({N, T, std::size_t(3)});
        HighFive::DataSetCreateProps props;
        props.add(HighFive::Chunking(std::vector<hsize_t>{
            static_cast<hsize_t>(N),
            static_cast<hsize_t>(std::max<std::size_t>(frame_chunk, 1u)),
            static_cast<hsize_t>(3)
        }));
        auto ds = grp.createDataSet<double>("charge_response_gradient_vector", space, props);
        ds.write_raw(flat.data());
    }

    {
        std::vector<double> flat(N * T);
        for (std::size_t i = 0; i < N; ++i) {
            const auto& atom_frames = per_atom_scalar_[i];
            for (std::size_t f = 0; f < T; ++f) {
                flat[i * T + f] = atom_frames[f];
            }
        }
        HighFive::DataSpace space({N, T});
        HighFive::DataSetCreateProps props;
        props.add(HighFive::Chunking(std::vector<hsize_t>{
            static_cast<hsize_t>(N),
            static_cast<hsize_t>(std::max<std::size_t>(frame_chunk, 1u))
        }));
        auto ds = grp.createDataSet<double>("charge_response_gradient_scalar", space, props);
        ds.write_raw(flat.data());
    }

    grp.createDataSet("frame_indices", frame_indices_);
    grp.createDataSet("frame_times",   frame_times_)
       .createAttribute("units", std::string("ps"));

    // Canonical SDK contract: source_attached_per_frame uint8 per-frame
    // mask. Compute's HasResult gate fills the mask correctly; normally
    // all-1 under always-attached policy.
    grp.createDataSet("source_attached_per_frame", source_attached_per_frame_);

    OperationLog::Info(LogCalcOther,
        "AIMNet2ChargeResponseGradientTimeSeriesTrajectoryResult::WriteH5Group",
        "wrote /trajectory/aimnet2_charge_response_gradient_time_series with " +
        std::to_string(N) + " atoms × " + std::to_string(T) + " frames "
        "(Vec3 + scalar, chunk " + std::to_string(frame_chunk) + " frames)");
}

}  // namespace nmr

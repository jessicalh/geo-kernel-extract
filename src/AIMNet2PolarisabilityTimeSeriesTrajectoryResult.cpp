#include "AIMNet2PolarisabilityTimeSeriesTrajectoryResult.h"

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
#include <string>

namespace nmr {

std::unique_ptr<AIMNet2PolarisabilityTimeSeriesTrajectoryResult>
AIMNet2PolarisabilityTimeSeriesTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<AIMNet2PolarisabilityTimeSeriesTrajectoryResult>();
    r->per_atom_vector_.assign(tp.AtomCount(), {});
    r->per_atom_scalar_.assign(tp.AtomCount(), {});
    return r;
}

void AIMNet2PolarisabilityTimeSeriesTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;
    const std::size_t N = per_atom_vector_.size();
    for (std::size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        per_atom_vector_[i].push_back(ca.aimnet2_polarisability_vector);
        per_atom_scalar_[i].push_back(ca.aimnet2_polarisability_scalar);
    }
    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    ++n_frames_;
}

void AIMNet2PolarisabilityTimeSeriesTrajectoryResult::Finalize(
        TrajectoryProtein& tp, Trajectory& traj) {
    (void)tp; (void)traj;
    finalized_ = true;
    OperationLog::Info(LogCalcOther,
        "AIMNet2PolarisabilityTimeSeriesTrajectoryResult::Finalize",
        "finalized across " + std::to_string(n_frames_) +
        " frames, " + std::to_string(per_atom_vector_.size()) + " atoms");
}

void AIMNet2PolarisabilityTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp, HighFive::File& file) const {
    const std::size_t N = per_atom_vector_.size();
    const std::size_t T = n_frames_;

    auto grp = file.createGroup("/trajectory/aimnet2_polarisability_time_series");

    grp.createAttribute("result_name",            Name());
    grp.createAttribute("n_atoms",                N);
    grp.createAttribute("n_frames",               T);
    grp.createAttribute("finalized",              finalized_);
    grp.createAttribute("units_vector",           std::string("e^2/Angstrom"));
    grp.createAttribute("units_scalar",           std::string("e^2/Angstrom"));
    grp.createAttribute("irrep_layout_vector",    std::string("1o"));
    grp.createAttribute("parity_vector",          std::string("1o"));
    grp.createAttribute("irrep_layout_scalar",    std::string("T0"));
    grp.createAttribute("parity_scalar",          std::string("0e"));
    grp.createAttribute("source",                 std::string(
        "AIMNet2PolarisabilityResult.{aimnet2_polarisability_vector (Vec3), "
        "aimnet2_polarisability_scalar (double, L2 norm of vector)}. "
        "Gradient of L = sum_j q_j^2 (AIMNet2 Hirshfeld charges, e^2) "
        "with respect to atomic coordinates (e^2/Angstrom). Both emitted "
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
        auto ds = grp.createDataSet<double>("polarisability_vector", space, props);
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
        auto ds = grp.createDataSet<double>("polarisability_scalar", space, props);
        ds.write_raw(flat.data());
    }

    grp.createDataSet("frame_indices", frame_indices_);
    grp.createDataSet("frame_times",   frame_times_)
       .createAttribute("units", std::string("ps"));

    // Canonical SDK contract: source_attached_per_frame all-1 mask for
    // always_attached TRs (OBJECT_MODEL.md "Conditional-attach TR
    // discipline" subsection). Trivially all-1 since AIMNet2PolarisabilityResult
    // is RequireConformationResult'd at the PerFrameExtractionSet config.
    std::vector<std::uint8_t> all_attached(T, 1u);
    grp.createDataSet("source_attached_per_frame", all_attached);

    OperationLog::Info(LogCalcOther,
        "AIMNet2PolarisabilityTimeSeriesTrajectoryResult::WriteH5Group",
        "wrote /trajectory/aimnet2_polarisability_time_series with " +
        std::to_string(N) + " atoms × " + std::to_string(T) + " frames "
        "(Vec3 + scalar, chunk " + std::to_string(frame_chunk) + " frames)");
}

}  // namespace nmr

#include "AIMNet2EmbeddingTimeSeriesTrajectoryResult.h"

#include "ConformationAtom.h"
#include "OperationLog.h"
#include "ProteinConformation.h"
#include "TrajectoryProtein.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <highfive/H5PropertyList.hpp>

#include <algorithm>
#include <cstddef>
#include <string>

namespace nmr {

std::unique_ptr<AIMNet2EmbeddingTimeSeriesTrajectoryResult>
AIMNet2EmbeddingTimeSeriesTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<AIMNet2EmbeddingTimeSeriesTrajectoryResult>();
    r->per_atom_embedding_.assign(tp.AtomCount(), {});
    return r;
}

void AIMNet2EmbeddingTimeSeriesTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;
    const std::size_t N = per_atom_embedding_.size();
    for (std::size_t i = 0; i < N; ++i) {
        per_atom_embedding_[i].push_back(conf.AtomAt(i).aimnet2_aim);
    }
    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    ++n_frames_;
}

void AIMNet2EmbeddingTimeSeriesTrajectoryResult::Finalize(
        TrajectoryProtein& tp, Trajectory& traj) {
    (void)tp; (void)traj;
    finalized_ = true;
    OperationLog::Info(LogCalcOther,
        "AIMNet2EmbeddingTimeSeriesTrajectoryResult::Finalize",
        "finalized across " + std::to_string(n_frames_) +
        " frames, " + std::to_string(per_atom_embedding_.size()) + " atoms");
}

void AIMNet2EmbeddingTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp, HighFive::File& file) const {
    const std::size_t N = per_atom_embedding_.size();
    const std::size_t T = n_frames_;
    constexpr std::size_t D = AIMNET2_AIM_DIMS;

    auto grp = file.createGroup("/trajectory/aimnet2_embedding_time_series");

    grp.createAttribute("result_name",            Name());
    grp.createAttribute("n_atoms",                N);
    grp.createAttribute("n_frames",               T);
    grp.createAttribute("embedding_dim",          D);
    grp.createAttribute("finalized",              finalized_);
    grp.createAttribute("units",                  std::string("dimensionless"));
    // No spherical-tensor structure for a 256-dim AIMNet2 feature vector;
    // explicit attrs so SDK introspection can dispatch (vs missing-attr).
    grp.createAttribute("irrep_layout",           std::string("feature_vector"));
    grp.createAttribute("parity",                 std::string("0e"));
    grp.createAttribute("source",                 std::string(
        "AIMNet2Result.aimnet2_aim (AIMNET2_AIM_DIMS=256, "
        "std::array<float, 256> per ConformationAtom). Always-attached "
        "source: AIMNet2Result is required by PerFrameExtractionSet; if "
        "the model is missing OperationRunner aborts before any TR fires."));
    grp.createAttribute("source_attached_policy", std::string("always_attached"));
    grp.createAttribute("optional_large",         true);

    // (N, T, D) float32. Chunk (1, frame_chunk=min(T, 64), D) per atom
    // so per-atom writes match the chunk row. Per-atom write_raw from
    // the existing per_atom_embedding_[i] contiguous storage (a
    // std::vector<std::array<float, D>>) avoids the (N*T*D) flat
    // allocation that would double peak memory on a fleet-scale run
    // (3.76 GB for 1P9J × 750 frames × 256 dims uncompressed; doubling
    // to 7.5 GB at flatten time would OOM larger fleet proteins).
    const std::size_t frame_chunk = std::min<std::size_t>(T, 64);
    HighFive::DataSpace space({N, T, D});
    HighFive::DataSetCreateProps props;
    props.add(HighFive::Chunking(std::vector<hsize_t>{
        static_cast<hsize_t>(1),
        static_cast<hsize_t>(std::max<std::size_t>(frame_chunk, 1u)),
        static_cast<hsize_t>(D)
    }));
    auto ds = grp.createDataSet<float>("embedding", space, props);

    // Per-atom hyperslab writes. std::vector<std::array<float, D>> is
    // contiguous (T arrays × D floats), so the per-atom slice can be
    // passed to write_raw via reinterpret_cast on the array's begin().
    for (std::size_t i = 0; i < N; ++i) {
        const auto& atom_frames = per_atom_embedding_[i];
        if (atom_frames.size() != T) continue;
        const std::vector<std::size_t> offset = {i, std::size_t(0), std::size_t(0)};
        const std::vector<std::size_t> count  = {std::size_t(1), T, D};
        ds.select(offset, count).write_raw(
            reinterpret_cast<const float*>(atom_frames.data()));
    }

    grp.createDataSet("frame_indices", frame_indices_);
    grp.createDataSet("frame_times",   frame_times_)
       .createAttribute("units", std::string("ps"));

    // Canonical SDK contract: source_attached_per_frame all-1 mask for
    // always_attached TRs (OBJECT_MODEL.md "Conditional-attach TR
    // discipline" subsection). Trivially all-1 since AIMNet2Result is
    // RequireConformationResult'd at the PerFrameExtractionSet config.
    std::vector<std::uint8_t> all_attached(T, 1u);
    grp.createDataSet("source_attached_per_frame", all_attached);

    OperationLog::Info(LogCalcOther,
        "AIMNet2EmbeddingTimeSeriesTrajectoryResult::WriteH5Group",
        "wrote /trajectory/aimnet2_embedding_time_series with " +
        std::to_string(N) + " atoms × " + std::to_string(T) + " frames × " +
        std::to_string(D) + " dims (float32, chunk " +
        std::to_string(frame_chunk) + " frames)");
}

}  // namespace nmr

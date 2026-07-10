#include "AIMNet2AimProjectionWelfordTrajectoryResult.h"

#include "AIMNet2Result.h"
#include "ConformationAtom.h"
#include "OperationLog.h"
#include "ProteinConformation.h"
#include "Trajectory.h"
#include "TrajectoryAtom.h"
#include "TrajectoryMoments.h"
#include "TrajectoryProtein.h"
#include "generated/AIMNet2AimProjection.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace nmr {

std::unique_ptr<AIMNet2AimProjectionWelfordTrajectoryResult>
AIMNet2AimProjectionWelfordTrajectoryResult::Create(
        const TrajectoryProtein& tp) {
    (void)tp;
    return std::make_unique<
        AIMNet2AimProjectionWelfordTrajectoryResult>();
}

void AIMNet2AimProjectionWelfordTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)traj;

    bool source_usable = conf.HasResult<AIMNet2Result>();
    if (!source_usable) {
        OperationLog::Warn(
            "AIMNet2AimProjectionWelfordTrajectoryResult::Compute",
            "AIMNet2Result not attached at frame " +
            std::to_string(frame_idx) +
            " — Welford update skipped; source mask=0 recorded");
    } else if (conf.AtomCount() != tp.AtomCount()) {
        OperationLog::Error(
            "AIMNet2AimProjectionWelfordTrajectoryResult::Compute",
            "conformation/trajectory atom-count mismatch at frame " +
            std::to_string(frame_idx) +
            " — Welford update skipped; source mask=0 recorded");
        source_usable = false;
    }

    // AIMNet2Result is a hard-failure calculator and normally guarantees
    // finite stored projections. Keep a defensive sentinel gate here so a
    // malformed test/custom source cannot NaN-poison Welford state. The bad
    // frame is logged and masked; no value is silently replaced with zero.
    if (source_usable) {
        for (std::size_t i = 0; i < conf.AtomCount() && source_usable; ++i) {
            const auto& projection =
                conf.AtomAt(i).aimnet2_aim_projection;
            for (std::size_t k = 0; k < kAimnet2AimProjectionDims; ++k) {
                if (!std::isfinite(projection[k])) {
                    OperationLog::Error(
                        "AIMNet2AimProjectionWelfordTrajectoryResult::Compute",
                        "non-finite stored projection at frame " +
                        std::to_string(frame_idx) + ", atom " +
                        std::to_string(i) + ", component " +
                        std::to_string(k) +
                        " — whole-frame Welford update skipped; source "
                        "mask=0 recorded");
                    source_usable = false;
                    break;
                }
            }
        }
    }

    if (source_usable) {
        for (std::size_t i = 0; i < tp.AtomCount(); ++i) {
            const auto& projection =
                conf.AtomAt(i).aimnet2_aim_projection;
            auto& state = tp.MutableAtomAt(i)
                .aimnet2_aim_projection_welford;
            const std::size_t n_new = state.n_frames + 1;
            for (std::size_t k = 0; k < kAimnet2AimProjectionDims; ++k) {
                // Pure read-back: projection[k] is already owned/stored by
                // AIMNet2Result::Compute. No raw aim or basis access here.
                WelfordUpdate(state.projection[k],
                              static_cast<double>(projection[k]),
                              n_new,
                              frame_idx);
            }
            state.n_frames = n_new;
        }
        ++source_attached_count_;
    }

    frame_indices_.push_back(static_cast<std::uint64_t>(frame_idx));
    frame_times_.push_back(time_ps);
    source_attached_per_frame_.push_back(source_usable ? 1u : 0u);
    ++n_frames_;
}

void AIMNet2AimProjectionWelfordTrajectoryResult::Finalize(
        TrajectoryProtein& tp,
        Trajectory& traj) {
    (void)traj;

    if (source_attached_count_ > 0) {
        for (std::size_t i = 0; i < tp.AtomCount(); ++i) {
            auto& state = tp.MutableAtomAt(i)
                .aimnet2_aim_projection_welford;
            for (std::size_t k = 0; k < kAimnet2AimProjectionDims; ++k) {
                WelfordFinalize(state.projection[k], state.n_frames);
            }
        }
    }
    finalized_ = true;
    OperationLog::Info(
        LogCalcOther,
        "AIMNet2AimProjectionWelfordTrajectoryResult::Finalize",
        "finalized across " + std::to_string(n_frames_) +
        " frames; source attached on " +
        std::to_string(source_attached_count_) + " frames");
}

void AIMNet2AimProjectionWelfordTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    const std::size_t N = tp.AtomCount();
    constexpr std::size_t K = kAimnet2AimProjectionDims;

    auto group = file.createGroup(
        "/trajectory/aimnet2_aim_projection_welford");
    group.createAttribute("result_name", Name());
    group.createAttribute("n_atoms", N);
    group.createAttribute("n_frames", n_frames_);
    group.createAttribute("source_attached_count", source_attached_count_);
    group.createAttribute("projection_dim", static_cast<std::uint64_t>(K));
    group.createAttribute("basis_id",
        std::string(kAimnet2AimProjectionBasisId));
    group.createAttribute("units", std::string("dimensionless"));
    group.createAttribute("irrep_layout", std::string("feature_vector"));
    group.createAttribute("parity", std::string("0e"));
    group.createAttribute("source", std::string(
        "AIMNet2Result.aimnet2_aim projected by committed element basis"));
    group.createAttribute("source_attached_policy", std::string(
        "always_attached_with_skip_on_absent_source"));
    group.createAttribute("finalized", finalized_);

    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> mean(N*K, nan);
    std::vector<double> m2(N*K, nan);
    std::vector<double> stddev(N*K, nan);
    std::vector<double> min_value(
        N*K, std::numeric_limits<double>::infinity());
    std::vector<double> max_value(
        N*K, -std::numeric_limits<double>::infinity());
    std::vector<std::uint64_t> min_frame(N*K, 0u);
    std::vector<std::uint64_t> max_frame(N*K, 0u);
    std::vector<std::uint64_t> n_per_atom(N, 0u);

    if (source_attached_count_ > 0) {
        for (std::size_t i = 0; i < N; ++i) {
            const auto& state = tp.AtomAt(i)
                .aimnet2_aim_projection_welford;
            n_per_atom[i] = static_cast<std::uint64_t>(state.n_frames);
            for (std::size_t k = 0; k < K; ++k) {
                const std::size_t offset = i*K + k;
                const auto& moments = state.projection[k];
                mean[offset] = moments.mean;
                m2[offset] = moments.m2;
                stddev[offset] = moments.std;
                min_value[offset] = moments.min;
                max_value[offset] = moments.max;
                min_frame[offset] =
                    static_cast<std::uint64_t>(moments.min_frame);
                max_frame[offset] =
                    static_cast<std::uint64_t>(moments.max_frame);
            }
        }
    }

    const HighFive::DataSpace component_space({N, K});
    auto write_double_2d = [&](const std::string& name,
                               const std::vector<double>& values,
                               const std::string& units) {
        auto dataset = group.createDataSet<double>(name, component_space);
        dataset.write_raw(values.data());
        dataset.createAttribute("units", units);
    };
    auto write_u64_2d = [&](const std::string& name,
                            const std::vector<std::uint64_t>& values) {
        auto dataset = group.createDataSet<std::uint64_t>(
            name, component_space);
        dataset.write_raw(values.data());
        dataset.createAttribute("units", std::string("frame_index"));
    };

    write_double_2d("projection_mean", mean, "dimensionless");
    write_double_2d("projection_m2", m2, "dimensionless_squared");
    write_double_2d("projection_std", stddev, "dimensionless");
    write_double_2d("projection_min", min_value, "dimensionless");
    write_double_2d("projection_max", max_value, "dimensionless");
    write_u64_2d("projection_min_frame", min_frame);
    write_u64_2d("projection_max_frame", max_frame);

    group.createDataSet("n_frames_per_atom", n_per_atom)
        .createAttribute("units", std::string("frame_count"));
    group.createDataSet("frame_indices", frame_indices_)
        .createAttribute("units", std::string("frame_index"));
    group.createDataSet("frame_times", frame_times_)
        .createAttribute("units", std::string("ps"));
    group.createDataSet("source_attached_per_frame",
                        source_attached_per_frame_);

    OperationLog::Info(
        LogCalcOther,
        "AIMNet2AimProjectionWelfordTrajectoryResult::WriteH5Group",
        "wrote /trajectory/aimnet2_aim_projection_welford for " +
        std::to_string(N) + " atoms x " + std::to_string(K) +
        " components");
}

}  // namespace nmr

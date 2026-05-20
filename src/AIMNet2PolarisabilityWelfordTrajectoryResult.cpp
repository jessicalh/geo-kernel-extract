#include "AIMNet2PolarisabilityWelfordTrajectoryResult.h"

#include "AIMNet2PolarisabilityResult.h"
#include "ConformationAtom.h"
#include "OperationLog.h"
#include "ProteinConformation.h"
#include "TrajectoryAtom.h"
#include "TrajectoryMoments.h"  // WelfordUpdate
#include "TrajectoryProtein.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <cstddef>
#include <string>

namespace nmr {

std::unique_ptr<AIMNet2PolarisabilityWelfordTrajectoryResult>
AIMNet2PolarisabilityWelfordTrajectoryResult::Create(
        const TrajectoryProtein& tp) {
    (void)tp;
    return std::make_unique<AIMNet2PolarisabilityWelfordTrajectoryResult>();
}

void AIMNet2PolarisabilityWelfordTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)traj;

    const bool source_present =
        conf.HasResult<AIMNet2PolarisabilityResult>();
    if (!source_present) {
        OperationLog::Warn(
            "AIMNet2PolarisabilityWelfordTrajectoryResult::Compute",
            "AIMNet2PolarisabilityResult not attached at frame " +
            std::to_string(frame_idx) +
            " — Welford accumulation skipped; mask=0 recorded.");
    } else {
        // Update per-atom Welford on the four channels (3 Vec3
        // components + scalar) via the canonical WelfordUpdate free
        // function. No delta variants in this minimum-viable v0.
        const std::size_t N = conf.AtomCount();
        for (std::size_t i = 0; i < N; ++i) {
            auto& ws = tp.MutableAtomAt(i).aimnet2_polarisability_welford;
            const auto& ca = conf.AtomAt(i);
            const Vec3&  v = ca.aimnet2_polarisability_vector;
            const double s = ca.aimnet2_polarisability_scalar;
            const std::size_t n_new = ws.n_frames + 1;
            WelfordUpdate(ws.polarisability_vector[0], v.x(), n_new, frame_idx);
            WelfordUpdate(ws.polarisability_vector[1], v.y(), n_new, frame_idx);
            WelfordUpdate(ws.polarisability_vector[2], v.z(), n_new, frame_idx);
            WelfordUpdate(ws.polarisability_scalar,    s,      n_new, frame_idx);
            ws.n_frames = n_new;
        }
        ++source_attached_count_;
    }

    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    source_attached_per_frame_.push_back(source_present ? 1u : 0u);
    ++n_frames_;
}

void AIMNet2PolarisabilityWelfordTrajectoryResult::Finalize(
        TrajectoryProtein& tp, Trajectory& traj) {
    (void)tp; (void)traj;
    finalized_ = true;
    OperationLog::Info(LogCalcOther,
        "AIMNet2PolarisabilityWelfordTrajectoryResult::Finalize",
        "finalized across " + std::to_string(n_frames_) +
        " frames; source attached on " +
        std::to_string(source_attached_count_) + " frames");
}

void AIMNet2PolarisabilityWelfordTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp, HighFive::File& file) const {
    const std::size_t N = tp.AtomCount();
    const std::size_t T = n_frames_;

    if (source_attached_count_ == 0) {
        OperationLog::Warn(
            "AIMNet2PolarisabilityWelfordTrajectoryResult::WriteH5Group",
            "source attached on 0/" + std::to_string(T) +
            " frames; skipping /trajectory/aimnet2_polarisability_welford/.");
        return;
    }

    auto grp = file.createGroup("/trajectory/aimnet2_polarisability_welford");

    grp.createAttribute("result_name",            Name());
    grp.createAttribute("n_atoms",                N);
    grp.createAttribute("n_frames",               T);
    grp.createAttribute("source_attached_count",  source_attached_count_);
    grp.createAttribute("finalized",              finalized_);
    grp.createAttribute("units_vector",           std::string("e^2/Angstrom"));
    grp.createAttribute("units_scalar",           std::string("e^2/Angstrom"));
    grp.createAttribute("irrep_layout_vector",    std::string("x,y,z"));
    grp.createAttribute("normalization_vector",   std::string("cartesian"));
    grp.createAttribute("parity_vector",          std::string("1o"));
    grp.createAttribute("irrep_layout_scalar",    std::string("T0"));
    grp.createAttribute("parity_scalar",          std::string("0e"));
    grp.createAttribute("source", std::string(
        "AIMNet2PolarisabilityResult.{aimnet2_polarisability_vector "
        "(Vec3), aimnet2_polarisability_scalar (double)}. Per-atom "
        "Welford rollup; minimum-viable v0 emits mean+M2+n only."));
    grp.createAttribute("source_attached_policy", std::string(
        "always_attached -- but Compute's HasResult<AIMNet2PolarisabilityResult>() "
        "gate skips Welford update + records mask=0 on absent frames "
        "(canonical \"absent, not faked\" gate)."));

    // (N, 3) Vec3-component channels + (N,) scalar channel. Allocate
    // flat buffers and copy from TrajectoryAtom Welford state. Sample
    // count `n` lives on the enclosing AIMNet2PolarisabilityWelfordState
    // (n_frames), not per-WelfordMoments — the per-channel counts are
    // identical because all four channels are updated together.
    std::vector<double>        vec_mean(N * 3), vec_m2(N * 3);
    std::vector<double>        scl_mean(N), scl_m2(N);
    std::vector<std::uint64_t> n_per_atom(N);
    for (std::size_t i = 0; i < N; ++i) {
        const auto& ws = tp.AtomAt(i).aimnet2_polarisability_welford;
        for (std::size_t c = 0; c < 3; ++c) {
            vec_mean[i * 3 + c] = ws.polarisability_vector[c].mean;
            vec_m2[i * 3 + c]   = ws.polarisability_vector[c].m2;
        }
        scl_mean[i]   = ws.polarisability_scalar.mean;
        scl_m2[i]     = ws.polarisability_scalar.m2;
        n_per_atom[i] = static_cast<std::uint64_t>(ws.n_frames);
    }

    {
        HighFive::DataSpace space({N, std::size_t(3)});
        grp.createDataSet<double>("vector_mean", space).write_raw(vec_mean.data());
        grp.createDataSet<double>("vector_m2",   space).write_raw(vec_m2.data());
    }
    {
        HighFive::DataSpace space({N});
        grp.createDataSet<double>("scalar_mean", space).write_raw(scl_mean.data());
        grp.createDataSet<double>("scalar_m2",   space).write_raw(scl_m2.data());
        grp.createDataSet<std::uint64_t>("n_per_atom", space).write_raw(n_per_atom.data());
    }

    grp.createDataSet("frame_indices",            frame_indices_);
    grp.createDataSet("frame_times",              frame_times_)
       .createAttribute("units", std::string("ps"));
    grp.createDataSet("source_attached_per_frame", source_attached_per_frame_);

    OperationLog::Info(LogCalcOther,
        "AIMNet2PolarisabilityWelfordTrajectoryResult::WriteH5Group",
        "wrote /trajectory/aimnet2_polarisability_welford with " +
        std::to_string(N) + " atoms (3-component Vec3 + scalar; " +
        std::to_string(source_attached_count_) + "/" + std::to_string(T) +
        " frames attached)");
}

}  // namespace nmr

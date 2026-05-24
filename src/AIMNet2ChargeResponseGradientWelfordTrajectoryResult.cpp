#include "AIMNet2ChargeResponseGradientWelfordTrajectoryResult.h"

#include "AIMNet2ChargeResponseGradientResult.h"
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

std::unique_ptr<AIMNet2ChargeResponseGradientWelfordTrajectoryResult>
AIMNet2ChargeResponseGradientWelfordTrajectoryResult::Create(
        const TrajectoryProtein& tp) {
    (void)tp;
    return std::make_unique<AIMNet2ChargeResponseGradientWelfordTrajectoryResult>();
}

void AIMNet2ChargeResponseGradientWelfordTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)traj;

    const bool source_present =
        conf.HasResult<AIMNet2ChargeResponseGradientResult>();
    if (!source_present) {
        OperationLog::Warn(
            "AIMNet2ChargeResponseGradientWelfordTrajectoryResult::Compute",
            "AIMNet2ChargeResponseGradientResult not attached at frame " +
            std::to_string(frame_idx) +
            " — Welford accumulation skipped; mask=0 recorded.");
    } else {
        // Update per-atom Welford on the four channels (3 Vec3
        // components + scalar) via the canonical WelfordUpdate free
        // function. No delta variants in this minimum-viable v0.
        const std::size_t N = conf.AtomCount();
        for (std::size_t i = 0; i < N; ++i) {
            auto& ws = tp.MutableAtomAt(i).aimnet2_charge_response_gradient_welford;
            const auto& ca = conf.AtomAt(i);
            const Vec3&  v = ca.aimnet2_charge_response_gradient_vector;
            const double s = ca.aimnet2_charge_response_gradient_scalar;
            const std::size_t n_new = ws.n_frames + 1;
            WelfordUpdate(ws.charge_response_gradient_vector[0], v.x(), n_new, frame_idx);
            WelfordUpdate(ws.charge_response_gradient_vector[1], v.y(), n_new, frame_idx);
            WelfordUpdate(ws.charge_response_gradient_vector[2], v.z(), n_new, frame_idx);
            WelfordUpdate(ws.charge_response_gradient_scalar,    s,      n_new, frame_idx);
            ws.n_frames = n_new;
        }
        ++source_attached_count_;
    }

    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    source_attached_per_frame_.push_back(source_present ? 1u : 0u);
    ++n_frames_;
}

void AIMNet2ChargeResponseGradientWelfordTrajectoryResult::Finalize(
        TrajectoryProtein& tp, Trajectory& traj) {
    (void)traj;
    // Per the canonical Welford TR pattern (HydrationGeometryWelford,
    // BsWelford, HmWelford, etc.): call WelfordFinalize per-atom-
    // per-channel to derive std and NaN-fill atoms with n_frames == 0
    // (uncomputable). This populates the std/min/max fields already
    // tracked by WelfordUpdate so downstream readers don't have to
    // recompute (math + science review HIGH 2026-05-20).
    // Only iterate when at least one frame attached the source —
    // otherwise the TrajectoryProtein's per-atom Welford slots may
    // not have been touched by Compute (test path) and MutableAtomAt
    // access is unsafe before the first Compute invocation populates
    // them. The WriteH5Group skip path also handles this case by
    // returning early when source_attached_count_ == 0.
    if (source_attached_count_ > 0) {
        const std::size_t N = tp.AtomCount();
        for (std::size_t i = 0; i < N; ++i) {
            auto& ws = tp.MutableAtomAt(i).aimnet2_charge_response_gradient_welford;
            WelfordFinalize(ws.charge_response_gradient_vector[0], ws.n_frames);
            WelfordFinalize(ws.charge_response_gradient_vector[1], ws.n_frames);
            WelfordFinalize(ws.charge_response_gradient_vector[2], ws.n_frames);
            WelfordFinalize(ws.charge_response_gradient_scalar,    ws.n_frames);
        }
    }
    finalized_ = true;
    OperationLog::Info(LogCalcOther,
        "AIMNet2ChargeResponseGradientWelfordTrajectoryResult::Finalize",
        "finalized across " + std::to_string(n_frames_) +
        " frames; source attached on " +
        std::to_string(source_attached_count_) + " frames");
}

void AIMNet2ChargeResponseGradientWelfordTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp, HighFive::File& file) const {
    const std::size_t N = tp.AtomCount();
    const std::size_t T = n_frames_;

    if (source_attached_count_ == 0) {
        OperationLog::Warn(
            "AIMNet2ChargeResponseGradientWelfordTrajectoryResult::WriteH5Group",
            "source attached on 0/" + std::to_string(T) +
            " frames; skipping /trajectory/aimnet2_charge_response_gradient_welford/.");
        return;
    }

    auto grp = file.createGroup("/trajectory/aimnet2_charge_response_gradient_welford");

    grp.createAttribute("result_name",            Name());
    grp.createAttribute("n_atoms",                N);
    grp.createAttribute("n_frames",               T);
    grp.createAttribute("source_attached_count",  source_attached_count_);
    grp.createAttribute("finalized",              finalized_);
    // Group-level units (base channel units). Per-dataset units below
    // distinguish base (mean/std/min/max) from squared (m2) and from
    // frame_index (min_frame/max_frame) — codex H3 2026-05-20.
    grp.createAttribute("units_vector",           std::string("e^2/Å"));
    grp.createAttribute("units_scalar",           std::string("e^2/Å"));
    grp.createAttribute("irrep_layout_vector",    std::string("x,y,z"));
    grp.createAttribute("normalization_vector",   std::string("cartesian"));
    grp.createAttribute("parity_vector",          std::string("1o"));
    grp.createAttribute("irrep_layout_scalar",    std::string("T0"));
    grp.createAttribute("parity_scalar",          std::string("0e"));
    grp.createAttribute("source", std::string(
        "AIMNet2ChargeResponseGradientResult.{aimnet2_charge_response_gradient_vector "
        "(Vec3), aimnet2_charge_response_gradient_scalar (double)}. Per-atom "
        "Welford rollup; emits mean + std + m2 + min/max + min_frame/"
        "max_frame + n_per_atom per channel (canonical "
        "WelfordFinalize-derived row matching sibling Welford TRs). "
        "scalar_mean = E[|v|] is NOT |E[v]| (= norm of vector_mean) — "
        "consumers must distinguish (sci review S10 2026-05-20)."));
    grp.createAttribute("source_attached_policy", std::string(
        "always_attached -- but Compute's HasResult<AIMNet2ChargeResponseGradientResult>() "
        "gate skips Welford update + records mask=0 on absent frames "
        "(canonical \"absent, not faked\" gate)."));
    grp.createAttribute("physics_description", std::string(
        "AIMNet2 charge-response gradient: dL/dr_i where L = sum_j q_j^2. "
        "NOT a Buckingham polarisability tensor (alpha_ab = d(mu_a)/d(E_b) is "
        "the conventional NMR-relevant alpha). The L = sum q^2 objective is "
        "a computationally cheap autograd-friendly summary of the per-atom "
        "charge sensitivity; physical-observable connection to NMR shielding "
        "is exploratory and calibration-ridge decides if signal carries beyond "
        "the AIMNet2 charge channel. Class renamed from project-shorthand "
        "'Polarisability' to ChargeResponseGradient 2026-05-20 (commit 58594f5)."));
    grp.createAttribute("noise_floor_caveat", std::string(
        "Per-atom std contains an autograd-noise contribution. AIMNet2 "
        "backward kernel uses non-deterministic cuBLAS matmul/scatter_add "
        "by default; gradient values can vary in the 6-10th significant "
        "figure on repeated runs of the same coords. For atoms in tightly-"
        "constrained regions the noise floor may match or exceed real "
        "physical variance. Math review M3 2026-05-20."));

    // Full canonical Welford row per sibling TR convention
    // (HydrationGeometryWelford, BsWelford, HmWelford, etc.):
    // mean + std + m2 + min + max + min_frame + max_frame per channel.
    // Sample count `n` shared across the 4 channels lives on the
    // enclosing AIMNet2ChargeResponseGradientWelfordState (n_frames).
    std::vector<double>        vec_mean(N * 3), vec_std(N * 3),
                               vec_m2(N * 3),   vec_min(N * 3), vec_max(N * 3);
    std::vector<std::uint64_t> vec_min_frame(N * 3), vec_max_frame(N * 3);
    std::vector<double>        scl_mean(N), scl_std(N), scl_m2(N),
                               scl_min(N),  scl_max(N);
    std::vector<std::uint64_t> scl_min_frame(N), scl_max_frame(N);
    std::vector<std::uint64_t> n_per_atom(N);
    for (std::size_t i = 0; i < N; ++i) {
        const auto& ws = tp.AtomAt(i).aimnet2_charge_response_gradient_welford;
        for (std::size_t c = 0; c < 3; ++c) {
            const auto& w = ws.charge_response_gradient_vector[c];
            vec_mean[i * 3 + c]      = w.mean;
            vec_std[i * 3 + c]       = w.std;
            vec_m2[i * 3 + c]        = w.m2;
            vec_min[i * 3 + c]       = w.min;
            vec_max[i * 3 + c]       = w.max;
            vec_min_frame[i * 3 + c] = static_cast<std::uint64_t>(w.min_frame);
            vec_max_frame[i * 3 + c] = static_cast<std::uint64_t>(w.max_frame);
        }
        const auto& s = ws.charge_response_gradient_scalar;
        scl_mean[i]      = s.mean;
        scl_std[i]       = s.std;
        scl_m2[i]        = s.m2;
        scl_min[i]       = s.min;
        scl_max[i]       = s.max;
        scl_min_frame[i] = static_cast<std::uint64_t>(s.min_frame);
        scl_max_frame[i] = static_cast<std::uint64_t>(s.max_frame);
        n_per_atom[i]    = static_cast<std::uint64_t>(ws.n_frames);
    }

    // Per-dataset unit attrs: mean/std/min/max in base units, m2 in
    // squared base units, *_frame in frame_index, n_per_atom in
    // frame_count. Codex H3 2026-05-20.
    const std::string base_v = "e^2/Å";
    const std::string sqr_v  = "(e^2/Å)^2";
    const std::string base_s = "e^2/Å";
    const std::string sqr_s  = "(e^2/Å)^2";
    const std::string fi     = "frame_index";
    const std::string fc     = "frame_count";

    auto write_double_attr = [](HighFive::DataSet& ds, const std::string& units) {
        ds.createAttribute("units", units);
    };
    {
        HighFive::DataSpace space({N, std::size_t(3)});
        auto d_vmn  = grp.createDataSet<double>("vector_mean",      space); d_vmn.write_raw(vec_mean.data());  write_double_attr(d_vmn, base_v);
        auto d_vsd  = grp.createDataSet<double>("vector_std",       space); d_vsd.write_raw(vec_std.data());   write_double_attr(d_vsd, base_v);
        auto d_vm2  = grp.createDataSet<double>("vector_m2",        space); d_vm2.write_raw(vec_m2.data());    write_double_attr(d_vm2, sqr_v);
        auto d_vmin = grp.createDataSet<double>("vector_min",       space); d_vmin.write_raw(vec_min.data());  write_double_attr(d_vmin, base_v);
        auto d_vmax = grp.createDataSet<double>("vector_max",       space); d_vmax.write_raw(vec_max.data());  write_double_attr(d_vmax, base_v);
        auto d_vmnf = grp.createDataSet<std::uint64_t>("vector_min_frame", space); d_vmnf.write_raw(vec_min_frame.data()); d_vmnf.createAttribute("units", fi);
        auto d_vmxf = grp.createDataSet<std::uint64_t>("vector_max_frame", space); d_vmxf.write_raw(vec_max_frame.data()); d_vmxf.createAttribute("units", fi);
    }
    {
        HighFive::DataSpace space({N});
        auto d_smn  = grp.createDataSet<double>("scalar_mean",      space); d_smn.write_raw(scl_mean.data());  write_double_attr(d_smn, base_s);
        auto d_ssd  = grp.createDataSet<double>("scalar_std",       space); d_ssd.write_raw(scl_std.data());   write_double_attr(d_ssd, base_s);
        auto d_sm2  = grp.createDataSet<double>("scalar_m2",        space); d_sm2.write_raw(scl_m2.data());    write_double_attr(d_sm2, sqr_s);
        auto d_smin = grp.createDataSet<double>("scalar_min",       space); d_smin.write_raw(scl_min.data());  write_double_attr(d_smin, base_s);
        auto d_smax = grp.createDataSet<double>("scalar_max",       space); d_smax.write_raw(scl_max.data());  write_double_attr(d_smax, base_s);
        auto d_smnf = grp.createDataSet<std::uint64_t>("scalar_min_frame", space); d_smnf.write_raw(scl_min_frame.data()); d_smnf.createAttribute("units", fi);
        auto d_smxf = grp.createDataSet<std::uint64_t>("scalar_max_frame", space); d_smxf.write_raw(scl_max_frame.data()); d_smxf.createAttribute("units", fi);
        auto d_npa  = grp.createDataSet<std::uint64_t>("n_per_atom",       space); d_npa.write_raw(n_per_atom.data());     d_npa.createAttribute("units", fc);
    }

    grp.createDataSet("frame_indices",            frame_indices_);
    grp.createDataSet("frame_times",              frame_times_)
       .createAttribute("units", std::string("ps"));
    grp.createDataSet("source_attached_per_frame", source_attached_per_frame_);

    OperationLog::Info(LogCalcOther,
        "AIMNet2ChargeResponseGradientWelfordTrajectoryResult::WriteH5Group",
        "wrote /trajectory/aimnet2_charge_response_gradient_welford with " +
        std::to_string(N) + " atoms (3-component Vec3 + scalar; " +
        std::to_string(source_attached_count_) + "/" + std::to_string(T) +
        " frames attached)");
}

}  // namespace nmr

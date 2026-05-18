#include "HydrationGeometryTimeSeriesTrajectoryResult.h"
#include "HydrationGeometryResult.h"
#include "ConformationAtom.h"
#include "ProteinConformation.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"
#include "OperationLog.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <cassert>
#include <limits>
#include <typeinfo>

namespace nmr {


std::vector<std::type_index>
HydrationGeometryTimeSeriesTrajectoryResult::Dependencies() const {
    return { std::type_index(typeid(HydrationGeometryResult)) };
}


std::unique_ptr<HydrationGeometryTimeSeriesTrajectoryResult>
HydrationGeometryTimeSeriesTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<HydrationGeometryTimeSeriesTrajectoryResult>();
    const std::size_t N = tp.AtomCount();
    r->dipole_vector_.assign(N, {});
    r->surface_normal_.assign(N, {});
    r->first_shell_count_.assign(N, {});
    r->half_shell_asymmetry_.assign(N, {});
    r->dipole_alignment_.assign(N, {});
    r->dipole_coherence_.assign(N, {});
    return r;
}


void HydrationGeometryTimeSeriesTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;
    const std::size_t N = conf.AtomCount();
    assert(N == tp.AtomCount());

    const bool source_attached = force_source_present_for_testing_
        || conf.HasResult<HydrationGeometryResult>();
    source_attached_per_frame_.push_back(source_attached ? 1u : 0u);

    for (std::size_t i = 0; i < N; ++i) {
        const auto& a = conf.AtomAt(i);
        if (source_attached) {
            dipole_vector_[i].push_back(a.water_dipole_vector);
            surface_normal_[i].push_back(a.water_surface_normal);
            first_shell_count_[i].push_back(
                static_cast<std::uint32_t>(a.sasa_first_shell_count));
            half_shell_asymmetry_[i].push_back(a.sasa_half_shell_asymmetry);
            dipole_alignment_[i].push_back(a.sasa_dipole_alignment);
            dipole_coherence_[i].push_back(a.sasa_dipole_coherence);
        } else {
            dipole_vector_[i].push_back(Vec3::Zero());
            surface_normal_[i].push_back(Vec3::Zero());
            first_shell_count_[i].push_back(0u);
            half_shell_asymmetry_[i].push_back(0.0);
            dipole_alignment_[i].push_back(0.0);
            dipole_coherence_[i].push_back(0.0);
        }
    }
    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    ++n_frames_;
}


void HydrationGeometryTimeSeriesTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                          Trajectory& traj) {
    (void)tp; (void)traj;
    finalized_ = true;
    OperationLog::Info(LogCalcOther,
        "HydrationGeometryTimeSeriesTrajectoryResult::Finalize",
        "captured (" + std::to_string(tp.AtomCount()) + " atoms × " +
        std::to_string(n_frames_) + " frames) hydration geometry");
}


void HydrationGeometryTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    const std::size_t N = tp.AtomCount();
    const std::size_t T = n_frames_;

    std::size_t source_attached_count = 0;
    for (auto v : source_attached_per_frame_) if (v) ++source_attached_count;
    if (source_attached_count == 0) {
        OperationLog::Warn(
            "HydrationGeometryTimeSeriesTrajectoryResult::WriteH5Group",
            "HydrationGeometryResult attached in 0/" +
            std::to_string(source_attached_per_frame_.size()) +
            " frames (no solvent loaded?); skipping "
            "/trajectory/hydration_geometry_time_series/ emission.");
        return;
    }

    auto grp = file.createGroup("/trajectory/hydration_geometry_time_series");
    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

    grp.createAttribute("result_name",                Name());
    grp.createAttribute("n_atoms",                    N);
    grp.createAttribute("n_frames",                   T);
    grp.createAttribute("source_attached_count",      source_attached_count);
    grp.createAttribute("finalized",                  finalized_);
    grp.createAttribute("dipole_vector_layout",       std::string("x,y,z"));
    grp.createAttribute("dipole_vector_parity",       std::string("1o"));
    // Units: Water::Dipole() (src/SolventEnvironment.h:32-38) returns
    // H_charge·(H_pos − O_pos) summed over both H atoms; charges in e,
    // positions in Å → result is e·Å. Convert to Debye via the standard
    // factor 1 e·Å = 4.80320 D if needed downstream. Per R5 codex
    // 2026-05-18: previous "Debye_unnormalised" label was wrong by that
    // factor; the values are in e·Å, raw sum (no normalisation to
    // unit dipole, no conversion to Debye).
    grp.createAttribute("dipole_vector_units",        std::string("e_Angstrom"));
    grp.createAttribute("surface_normal_layout",      std::string("x,y,z"));
    grp.createAttribute("surface_normal_parity",      std::string("1o"));
    // Reference-frame disambiguation: HydrationGeometry uses the SASA
    // outward surface normal as the reference vector for cos-alignment
    // and half-shell asymmetry; the sibling HydrationShell TR uses COM.
    // Same dataset name (`half_shell_asymmetry`) under different group
    // paths means the same physical question in different reference
    // frames — never aggregate cross-frame without checking this attr.
    grp.createAttribute("reference_frame",            std::string("SASA_normal"));
    grp.createAttribute("polarisation_signal_channels",
        std::string("dipole_alignment,dipole_coherence,half_shell_asymmetry"));

    auto emit_vec3 = [&](const std::string& name,
                         const std::vector<std::vector<Vec3>>& src,
                         const std::string& units) {
        std::vector<double> flat(N * T * 3);
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t t = 0; t < T; ++t) {
                const std::size_t base = (i * T + t) * 3;
                if (!source_attached_per_frame_[t]) {
                    flat[base+0] = kNaN; flat[base+1] = kNaN; flat[base+2] = kNaN;
                    continue;
                }
                const Vec3& v = src[i][t];
                flat[base+0] = v.x(); flat[base+1] = v.y(); flat[base+2] = v.z();
            }
        }
        HighFive::DataSpace space({N, T, std::size_t(3)});
        auto ds = grp.createDataSet<double>(name, space);
        ds.write_raw(flat.data());
        ds.createAttribute("units", units);
    };
    emit_vec3("dipole_vector",  dipole_vector_,  "e_Angstrom");
    emit_vec3("surface_normal", surface_normal_, "unit_vector");

    auto emit_scalar = [&](const std::string& name,
                           const std::vector<std::vector<double>>& src,
                           const std::string& units) {
        std::vector<double> flat(N * T);
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t t = 0; t < T; ++t) {
                flat[i*T+t] = source_attached_per_frame_[t]
                            ? src[i][t]
                            : kNaN;
            }
        }
        HighFive::DataSpace space({N, T});
        auto ds = grp.createDataSet<double>(name, space);
        ds.write_raw(flat.data());
        ds.createAttribute("units", units);
    };
    emit_scalar("half_shell_asymmetry", half_shell_asymmetry_, "fraction");
    emit_scalar("dipole_alignment",     dipole_alignment_,     "cos_angle");
    emit_scalar("dipole_coherence",     dipole_coherence_,     "order_parameter");

    // Shell count uint32; absent sentinel like WaterFieldTS.
    {
        constexpr std::uint32_t kAbsent = 0xFFFFFFFFu;
        std::vector<std::uint32_t> flat(N * T);
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t t = 0; t < T; ++t) {
                flat[i*T+t] = source_attached_per_frame_[t]
                            ? first_shell_count_[i][t]
                            : kAbsent;
            }
        }
        HighFive::DataSpace space({N, T});
        auto ds = grp.createDataSet<std::uint32_t>("first_shell_count", space);
        ds.write_raw(flat.data());
        ds.createAttribute("units", std::string("dimensionless"));
        ds.createAttribute("absent_sentinel", kAbsent);
    }

    auto fi_ds = grp.createDataSet("frame_indices", frame_indices_);
    fi_ds.createAttribute("units", std::string("frame_index"));
    auto ft_ds = grp.createDataSet("frame_times",   frame_times_);
    ft_ds.createAttribute("units", std::string("ps"));
    auto sa_ds = grp.createDataSet("source_attached_per_frame",
                                   source_attached_per_frame_);
    sa_ds.createAttribute("units", std::string("dimensionless"));
}

}  // namespace nmr

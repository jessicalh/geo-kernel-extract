#include "HydrationShellTimeSeriesTrajectoryResult.h"
#include "HydrationShellResult.h"
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
HydrationShellTimeSeriesTrajectoryResult::Dependencies() const {
    return { std::type_index(typeid(HydrationShellResult)) };
}


std::unique_ptr<HydrationShellTimeSeriesTrajectoryResult>
HydrationShellTimeSeriesTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<HydrationShellTimeSeriesTrajectoryResult>();
    const std::size_t N = tp.AtomCount();
    r->half_shell_asymmetry_.assign(N, {});
    r->mean_water_dipole_cos_.assign(N, {});
    r->nearest_ion_distance_.assign(N, {});
    r->nearest_ion_charge_.assign(N, {});
    return r;
}


void HydrationShellTimeSeriesTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;
    const std::size_t N = conf.AtomCount();
    assert(N == tp.AtomCount());

    const bool source_attached = force_source_present_for_testing_
        || conf.HasResult<HydrationShellResult>();
    source_attached_per_frame_.push_back(source_attached ? 1u : 0u);

    for (std::size_t i = 0; i < N; ++i) {
        const auto& a = conf.AtomAt(i);
        if (source_attached) {
            half_shell_asymmetry_[i].push_back(a.half_shell_asymmetry);
            mean_water_dipole_cos_[i].push_back(a.mean_water_dipole_cos);
            nearest_ion_distance_[i].push_back(a.nearest_ion_distance);
            nearest_ion_charge_[i].push_back(a.nearest_ion_charge);
        } else {
            half_shell_asymmetry_[i].push_back(0.0);
            mean_water_dipole_cos_[i].push_back(0.0);
            nearest_ion_distance_[i].push_back(0.0);
            nearest_ion_charge_[i].push_back(0.0);
        }
    }
    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    ++n_frames_;
}


void HydrationShellTimeSeriesTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                       Trajectory& traj) {
    (void)tp; (void)traj;
    finalized_ = true;
    OperationLog::Info(LogCalcOther,
        "HydrationShellTimeSeriesTrajectoryResult::Finalize",
        "captured (" + std::to_string(tp.AtomCount()) + " atoms × " +
        std::to_string(n_frames_) + " frames) hydration shell");
}


void HydrationShellTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    const std::size_t N = tp.AtomCount();
    const std::size_t T = n_frames_;

    std::size_t source_attached_count = 0;
    for (auto v : source_attached_per_frame_) if (v) ++source_attached_count;
    if (source_attached_count == 0) {
        OperationLog::Warn(
            "HydrationShellTimeSeriesTrajectoryResult::WriteH5Group",
            "HydrationShellResult attached in 0/" +
            std::to_string(source_attached_per_frame_.size()) +
            " frames; skipping /trajectory/hydration_shell_time_series/.");
        return;
    }

    auto grp = file.createGroup("/trajectory/hydration_shell_time_series");
    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

    grp.createAttribute("result_name",            Name());
    grp.createAttribute("n_atoms",                N);
    grp.createAttribute("n_frames",               T);
    grp.createAttribute("source_attached_count",  source_attached_count);
    grp.createAttribute("finalized",              finalized_);
    grp.createAttribute("reference_frame",
        std::string("COM"));  // COM, not SASA — distinguishes from HydrationGeometry
    grp.createAttribute("nearest_ion_distance_sentinel",
        std::string("+infinity when no ion within ion_cutoff (default 20 A)"));
    // Physically-zero vs no-ion ambiguity on nearest_ion_charge: source
    // emits 0.0 for both "no ion in cutoff" AND "nearest ion happens to
    // be neutral" (latter impossible for Na+/Cl- but data type permits).
    // Downstream that sees nearest_ion_charge == 0.0 MUST cross-check
    // nearest_ion_distance == +inf to distinguish the two cases.
    grp.createAttribute("nearest_ion_charge_zero_means_no_ion_in_cutoff",
        std::string("0.0 on nearest_ion_charge AND +inf on "
                    "nearest_ion_distance jointly indicate 'no ion in "
                    "ion_cutoff' — neither alone is unambiguous."));

    auto emit_scalar = [&](const std::string& name,
                           const std::vector<std::vector<double>>& src,
                           const std::string& units) {
        std::vector<double> flat(N * T);
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t t = 0; t < T; ++t) {
                flat[i*T+t] = source_attached_per_frame_[t] ? src[i][t] : kNaN;
            }
        }
        HighFive::DataSpace space({N, T});
        auto ds = grp.createDataSet<double>(name, space);
        ds.write_raw(flat.data());
        ds.createAttribute("units", units);
    };
    emit_scalar("half_shell_asymmetry",  half_shell_asymmetry_,  "fraction");
    emit_scalar("mean_water_dipole_cos", mean_water_dipole_cos_, "cos_angle");
    emit_scalar("nearest_ion_distance",  nearest_ion_distance_,  "Angstrom");
    emit_scalar("nearest_ion_charge",    nearest_ion_charge_,    "e");

    auto fi_ds = grp.createDataSet("frame_indices", frame_indices_);
    fi_ds.createAttribute("units", std::string("frame_index"));
    auto ft_ds = grp.createDataSet("frame_times",   frame_times_);
    ft_ds.createAttribute("units", std::string("ps"));
    auto sa_ds = grp.createDataSet("source_attached_per_frame",
                                   source_attached_per_frame_);
    sa_ds.createAttribute("units", std::string("dimensionless"));
}

}  // namespace nmr

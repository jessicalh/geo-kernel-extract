#include "WaterFieldTimeSeriesTrajectoryResult.h"
#include "WaterFieldResult.h"
#include "ProteinConformation.h"
#include "ConformationAtom.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"
#include "OperationLog.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <cassert>
#include <typeinfo>

namespace nmr {


std::vector<std::type_index>
WaterFieldTimeSeriesTrajectoryResult::Dependencies() const {
    return { std::type_index(typeid(WaterFieldResult)) };
}


std::unique_ptr<WaterFieldTimeSeriesTrajectoryResult>
WaterFieldTimeSeriesTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<WaterFieldTimeSeriesTrajectoryResult>();
    const std::size_t N = tp.AtomCount();
    r->efield_.assign(N, {});
    r->efield_first_.assign(N, {});
    r->efg_.assign(N, {});
    r->efg_first_.assign(N, {});
    r->n_first_.assign(N, {});
    r->n_second_.assign(N, {});
    return r;
}


// ── Compute ──────────────────────────────────────────────────────

void WaterFieldTimeSeriesTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;
    const std::size_t N = conf.AtomCount();
    // Mirror BondedEnergy::Compute discipline: the per-atom buffers were
    // sized to tp.AtomCount() at Create; the architectural invariant is
    // conf.AtomCount() == tp.AtomCount() each frame. If that ever breaks,
    // catch it loudly rather than silently producing ragged buffers.
    assert(N == tp.AtomCount());
    for (std::size_t i = 0; i < N; ++i) {
        const auto& a = conf.AtomAt(i);
        efield_[i].push_back(a.water_efield);
        efield_first_[i].push_back(a.water_efield_first);
        efg_[i].push_back(a.water_efg_spherical);
        efg_first_[i].push_back(a.water_efg_first_spherical);
        n_first_[i].push_back(static_cast<std::uint32_t>(a.water_n_first));
        n_second_[i].push_back(static_cast<std::uint32_t>(a.water_n_second));
    }
    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    ++n_frames_;
}


// ── Finalize ─────────────────────────────────────────────────────

void WaterFieldTimeSeriesTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                    Trajectory& traj) {
    (void)tp; (void)traj;
    finalized_ = true;
    OperationLog::Info(LogCalcOther,
        "WaterFieldTimeSeriesTrajectoryResult::Finalize",
        "captured (" + std::to_string(tp.AtomCount()) + " atoms × " +
        std::to_string(n_frames_) + " frames) water E-field + EFG");
}


// ── WriteH5Group ─────────────────────────────────────────────────

void WaterFieldTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    const std::size_t N = tp.AtomCount();
    const std::size_t T = n_frames_;
    auto grp = file.createGroup("/trajectory/water_field_time_series");

    grp.createAttribute("result_name",      Name());
    grp.createAttribute("n_atoms",          N);
    grp.createAttribute("n_frames",         T);
    grp.createAttribute("finalized",        finalized_);
    grp.createAttribute("efield_layout",    std::string("x,y,z"));
    grp.createAttribute("efield_units",     std::string("V/Angstrom"));
    // E-field is a polar (true) vector — parity-odd under inversion.
    // Matches ApbsEfieldTimeSeries / PositionsTimeSeries convention.
    grp.createAttribute("efield_parity",        std::string("1o"));
    grp.createAttribute("efield_normalization", std::string("cartesian"));

    grp.createAttribute("efg_irrep_layout",
        std::string("T0,T1_m-1,T1_m0,T1_m+1,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"));
    grp.createAttribute("efg_units",        std::string("V/Angstrom^2"));
    // EFG is symmetric-traceless rank-2 — all-parity-even in e3nn convention.
    // T0 is structurally zero (WaterFieldResult traceless-projects V before
    // Decompose); kept in the irrep layout for shape symmetry with BS shielding.
    grp.createAttribute("efg_parity",         std::string("0e+1e+2e"));
    grp.createAttribute("efg_normalization",  std::string("isometric_real_sph"));
    grp.createAttribute("efg_t0_structural_zero", true);

    grp.createAttribute("count_units",      std::string("dimensionless"));

    // Cutoff radii — efield uses a 15 Å sphere; n_first / n_second use
    // hydration-shell cutoffs at 3.5 / 5.5 Å (TIP3P standard). Downstream
    // must not assume same support across these channels.
    grp.createAttribute("efield_cutoff_A",   15.0);
    grp.createAttribute("n_first_cutoff_A",  3.5);
    grp.createAttribute("n_second_cutoff_A", 5.5);

    // E-field as (N, T, 3)
    auto emit_vec3 = [&](const std::string& name,
                         const std::vector<std::vector<Vec3>>& src,
                         const std::string& units) {
        std::vector<double> flat(N * T * 3);
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t t = 0; t < T; ++t) {
                const Vec3& v = src[i][t];
                const std::size_t base = (i * T + t) * 3;
                flat[base + 0] = v.x();
                flat[base + 1] = v.y();
                flat[base + 2] = v.z();
            }
        }
        HighFive::DataSpace space({N, T, std::size_t(3)});
        auto ds = grp.createDataSet<double>(name, space);
        ds.write_raw(flat.data());
        ds.createAttribute("units", units);
    };
    emit_vec3("efield",       efield_,       "V/Angstrom");
    emit_vec3("efield_first", efield_first_, "V/Angstrom");

    // EFG SphericalTensor as (N, T, 9): T0, T1[-1..+1], T2[-2..+2]
    auto emit_sph = [&](const std::string& name,
                        const std::vector<std::vector<SphericalTensor>>& src,
                        const std::string& units) {
        std::vector<double> flat(N * T * 9);
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t t = 0; t < T; ++t) {
                const SphericalTensor& st = src[i][t];
                const std::size_t base = (i * T + t) * 9;
                flat[base + 0] = st.T0;
                flat[base + 1] = st.T1[0];
                flat[base + 2] = st.T1[1];
                flat[base + 3] = st.T1[2];
                flat[base + 4] = st.T2[0];
                flat[base + 5] = st.T2[1];
                flat[base + 6] = st.T2[2];
                flat[base + 7] = st.T2[3];
                flat[base + 8] = st.T2[4];
            }
        }
        HighFive::DataSpace space({N, T, std::size_t(9)});
        auto ds = grp.createDataSet<double>(name, space);
        ds.write_raw(flat.data());
        ds.createAttribute("units", units);
    };
    emit_sph("efg",       efg_,       "V/Angstrom^2");
    emit_sph("efg_first", efg_first_, "V/Angstrom^2");

    // Shell counts as (N, T) uint32
    auto emit_count = [&](const std::string& name,
                          const std::vector<std::vector<std::uint32_t>>& src) {
        std::vector<std::uint32_t> flat(N * T);
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t t = 0; t < T; ++t) {
                flat[i * T + t] = src[i][t];
            }
        }
        HighFive::DataSpace space({N, T});
        auto ds = grp.createDataSet<std::uint32_t>(name, space);
        ds.write_raw(flat.data());
        ds.createAttribute("units", std::string("dimensionless"));
    };
    emit_count("n_first",  n_first_);
    emit_count("n_second", n_second_);

    auto fi_ds = grp.createDataSet("frame_indices", frame_indices_);
    fi_ds.createAttribute("units", std::string("frame_index"));
    auto ft_ds = grp.createDataSet("frame_times",   frame_times_);
    ft_ds.createAttribute("units", std::string("ps"));
}

}  // namespace nmr

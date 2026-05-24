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
#include <limits>
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
    assert(N == tp.AtomCount());

    // "Absent, not faked" — record whether WaterFieldResult attached
    // this frame. ConformationAtom defaults for water_efield etc. are
    // Vec3::Zero / default SphericalTensor / 0; reading them when the
    // source wasn't attached would emit a valid-looking all-zero
    // H5 group that mimics a solvated extraction.
    const bool source_attached = force_source_present_for_testing_
        || conf.HasResult<WaterFieldResult>();
    source_attached_per_frame_.push_back(source_attached ? 1u : 0u);

    for (std::size_t i = 0; i < N; ++i) {
        const auto& a = conf.AtomAt(i);
        if (source_attached) {
            efield_[i].push_back(a.water_efield);
            efield_first_[i].push_back(a.water_efield_first);
            efg_[i].push_back(a.water_efg_spherical);
            efg_first_[i].push_back(a.water_efg_first_spherical);
            n_first_[i].push_back(static_cast<std::uint32_t>(a.water_n_first));
            n_second_[i].push_back(static_cast<std::uint32_t>(a.water_n_second));
        } else {
            // Push zeros / defaults; NaN-fill at WriteH5Group via the mask.
            efield_[i].push_back(Vec3::Zero());
            efield_first_[i].push_back(Vec3::Zero());
            efg_[i].push_back(SphericalTensor{});
            efg_first_[i].push_back(SphericalTensor{});
            n_first_[i].push_back(0u);
            n_second_[i].push_back(0u);
        }
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

    // "Absent, not faked" — skip emission entirely when source never attached.
    std::size_t source_attached_count = 0;
    for (auto v : source_attached_per_frame_) if (v) ++source_attached_count;
    if (source_attached_count == 0) {
        OperationLog::Warn(
            "WaterFieldTimeSeriesTrajectoryResult::WriteH5Group",
            "WaterFieldResult attached in 0/" +
            std::to_string(source_attached_per_frame_.size()) +
            " frames (no solvent loaded? protein-only extraction?); skipping "
            "/trajectory/water_field_time_series/ emission.");
        return;
    }

    auto grp = file.createGroup("/trajectory/water_field_time_series");
    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

    grp.createAttribute("result_name",      Name());
    grp.createAttribute("n_atoms",          N);
    grp.createAttribute("n_frames",         T);
    grp.createAttribute("source_attached_count", source_attached_count);
    grp.createAttribute("finalized",        finalized_);
    grp.createAttribute("efield_layout",    std::string("x,y,z"));
    grp.createAttribute("efield_units",     std::string("V/Angstrom"));
    // E-field is a polar (true) vector — parity-odd under inversion.
    // Matches ApbsEfieldTimeSeries / PositionsTimeSeries convention.
    grp.createAttribute("efield_parity",        std::string("1o"));
    grp.createAttribute("efield_normalization", std::string("cartesian"));

    // Water EFG layout: T2-only, 5 real-spherical-tesseral components.
    // T0 = 0 (traceless projection in WaterFieldResult.cpp:147-150).
    // T1 = 0 (water EFG built from r⊗r outer products → symmetric →
    // antisymmetric pseudovector vanishes; src/WaterFieldResult.cpp:130,
    // src/Types.cpp:28). Both T0 and T1 are STRUCTURAL zeros, not
    // numerical zeros from a particular run — so they are NOT emitted
    // (no kernel pollution / no consumer ambiguity).
    grp.createAttribute("efg_irrep_layout",
        std::string("T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"));
    grp.createAttribute("efg_units",        std::string("V/Angstrom^2"));
    grp.createAttribute("efg_parity",         std::string("2e"));
    grp.createAttribute("efg_normalization",  std::string("isometric_real_sph"));
    grp.createAttribute("efg_t0_structural_zero", true);
    grp.createAttribute("efg_t1_structural_zero", true);

    grp.createAttribute("count_units",      std::string("dimensionless"));

    // Cutoff radii — efield uses a 15 Å sphere; n_first / n_second use
    // hydration-shell cutoffs at 3.5 / 5.5 Å (TIP3P standard). Downstream
    // must not assume same support across these channels.
    grp.createAttribute("efield_cutoff_A",   15.0);
    grp.createAttribute("n_first_cutoff_A",  3.5);
    grp.createAttribute("n_second_cutoff_A", 5.5);

    // E-field as (N, T, 3); NaN-filled on source-absent frames.
    auto emit_vec3 = [&](const std::string& name,
                         const std::vector<std::vector<Vec3>>& src,
                         const std::string& units) {
        std::vector<double> flat(N * T * 3);
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t t = 0; t < T; ++t) {
                const std::size_t base = (i * T + t) * 3;
                if (!source_attached_per_frame_[t]) {
                    flat[base + 0] = kNaN;
                    flat[base + 1] = kNaN;
                    flat[base + 2] = kNaN;
                    continue;
                }
                const Vec3& v = src[i][t];
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

    // EFG SphericalTensor flat as (N, T, 5) — T2 only. Water EFG is
    // symmetric (built from r⊗r outer products in WaterFieldResult.cpp:130)
    // and traceless-projected, so T0 = trace = 0 AND T1 = antisymmetric
    // pseudovector = 0 are both structurally zero. Only T2 carries
    // signal. Layout: m=-2,m=-1,m=0,m=+1,m=+2 (real-spherical-tesseral).
    auto emit_sph_t2 = [&](const std::string& name,
                           const std::vector<std::vector<SphericalTensor>>& src,
                           const std::string& units) {
        std::vector<double> flat(N * T * 5);
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t t = 0; t < T; ++t) {
                const std::size_t base = (i * T + t) * 5;
                if (!source_attached_per_frame_[t]) {
                    for (std::size_t k = 0; k < 5; ++k) flat[base + k] = kNaN;
                    continue;
                }
                const SphericalTensor& st = src[i][t];
                for (std::size_t k = 0; k < 5; ++k) flat[base + k] = st.T2[k];
            }
        }
        HighFive::DataSpace space({N, T, std::size_t(5)});
        auto ds = grp.createDataSet<double>(name, space);
        ds.write_raw(flat.data());
        ds.createAttribute("units", units);
    };
    emit_sph_t2("efg",       efg_,       "V/Angstrom^2");
    emit_sph_t2("efg_first", efg_first_, "V/Angstrom^2");

    // Shell counts as (N, T) uint32 — sentinel 0xFFFFFFFFu marks
    // source-absent frames (uint32 can't carry NaN; max-value sentinel
    // is the convention).
    auto emit_count = [&](const std::string& name,
                          const std::vector<std::vector<std::uint32_t>>& src) {
        std::vector<std::uint32_t> flat(N * T);
        constexpr std::uint32_t kAbsent = 0xFFFFFFFFu;
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t t = 0; t < T; ++t) {
                flat[i * T + t] = source_attached_per_frame_[t]
                                ? src[i][t]
                                : kAbsent;
            }
        }
        HighFive::DataSpace space({N, T});
        auto ds = grp.createDataSet<std::uint32_t>(name, space);
        ds.write_raw(flat.data());
        ds.createAttribute("units", std::string("dimensionless"));
        ds.createAttribute("absent_sentinel", kAbsent);
    };
    emit_count("n_first",  n_first_);
    emit_count("n_second", n_second_);

    auto fi_ds = grp.createDataSet("frame_indices", frame_indices_);
    fi_ds.createAttribute("units", std::string("frame_index"));
    auto ft_ds = grp.createDataSet("frame_times",   frame_times_);
    ft_ds.createAttribute("units", std::string("ps"));
    auto sa_ds = grp.createDataSet("source_attached_per_frame",
                                   source_attached_per_frame_);
    sa_ds.createAttribute("units", std::string("dimensionless"));
}

}  // namespace nmr

#include "BondedEnergyTimeSeriesTrajectoryResult.h"
#include "BondedEnergyResult.h"
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
BondedEnergyTimeSeriesTrajectoryResult::Dependencies() const {
    return { std::type_index(typeid(BondedEnergyResult)) };
}


std::unique_ptr<BondedEnergyTimeSeriesTrajectoryResult>
BondedEnergyTimeSeriesTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<BondedEnergyTimeSeriesTrajectoryResult>();
    const std::size_t N = tp.AtomCount();
    r->bond_.assign(N, {});
    r->angle_.assign(N, {});
    r->urey_bradley_.assign(N, {});
    r->proper_dih_.assign(N, {});
    r->improper_dih_.assign(N, {});
    r->cmap_.assign(N, {});
    r->total_.assign(N, {});
    return r;
}


// ── Compute ──────────────────────────────────────────────────────
//
// Append this frame's per-atom bonded breakdown to the growing buffers.
// Source is BondedEnergyResult on the conformation. The seven channels
// come directly from its public accessors; no derivation here.

void BondedEnergyTimeSeriesTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;
    const std::size_t N = conf.AtomCount();

    // "Absent, not faked" — record whether BondedEnergyResult attached
    // this frame. tp.HasBondedParams()==false silently skips the source
    // in OperationRunner; calling conf.Result<BondedEnergyResult>() in
    // that state hard-aborts.
    const bool source_attached = force_source_present_for_testing_
        || conf.HasResult<BondedEnergyResult>();
    source_attached_per_frame_.push_back(source_attached ? 1u : 0u);

    if (source_attached) {
        const auto& be = conf.Result<BondedEnergyResult>();
        assert(be.BondEnergy().size() == N);
        for (std::size_t i = 0; i < N; ++i) {
            bond_[i].push_back(be.BondEnergy()[i]);
            angle_[i].push_back(be.AngleEnergy()[i]);
            urey_bradley_[i].push_back(be.UBEnergy()[i]);
            proper_dih_[i].push_back(be.ProperDihEnergy()[i]);
            improper_dih_[i].push_back(be.ImproperDihEnergy()[i]);
            cmap_[i].push_back(be.CmapEnergy()[i]);
            total_[i].push_back(be.TotalBonded()[i]);
        }
    } else {
        // Push zero placeholders. NaN-filled at WriteH5Group time so
        // downstream distinguishes "no measurement" from "real zero."
        for (std::size_t i = 0; i < N; ++i) {
            bond_[i].push_back(0.0);
            angle_[i].push_back(0.0);
            urey_bradley_[i].push_back(0.0);
            proper_dih_[i].push_back(0.0);
            improper_dih_[i].push_back(0.0);
            cmap_[i].push_back(0.0);
            total_[i].push_back(0.0);
        }
    }
    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    ++n_frames_;
}


// ── Finalize ─────────────────────────────────────────────────────

void BondedEnergyTimeSeriesTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                      Trajectory& traj) {
    (void)tp; (void)traj;
    finalized_ = true;
    OperationLog::Info(LogCalcOther,
        "BondedEnergyTimeSeriesTrajectoryResult::Finalize",
        "captured (" + std::to_string(tp.AtomCount()) + " atoms × " +
        std::to_string(n_frames_) + " frames) bonded breakdown");
}


// ── WriteH5Group ─────────────────────────────────────────────────
//
// Seven (N, T) datasets, one per channel. Layout matches the source
// BondedEnergyResult's public accessor order:
//   bond / angle / urey_bradley / proper_dih / improper_dih / cmap / total
// CHARMM36m through GROMACS reports UB=0 / improper=0 / CMAP=0 on the
// 1P9J fleet fixture (verified 2026-05-18); the datasets are still
// emitted so downstream cannot trip on a "missing column" error and
// alternative force fields that DO populate those interactions hit
// the same schema.

void BondedEnergyTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    const std::size_t N = tp.AtomCount();
    const std::size_t T = n_frames_;

    // "Absent, not faked" — skip entirely when source never attached.
    std::size_t source_attached_count = 0;
    for (auto v : source_attached_per_frame_) if (v) ++source_attached_count;
    if (source_attached_count == 0) {
        OperationLog::Warn(
            "BondedEnergyTimeSeriesTrajectoryResult::WriteH5Group",
            "BondedEnergyResult attached in 0/" +
            std::to_string(source_attached_per_frame_.size()) +
            " frames (no .tpr bonded params?); skipping "
            "/trajectory/bonded_energy_time_series/ emission.");
        return;
    }

    auto grp = file.createGroup("/trajectory/bonded_energy_time_series");
    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_atoms",     N);
    grp.createAttribute("n_frames",    T);
    grp.createAttribute("source_attached_count", source_attached_count);
    grp.createAttribute("finalized",   finalized_);
    grp.createAttribute("units",       std::string("kJ/mol"));
    // The `evenly_among_participating_atoms` convention is one of
    // several valid attributions for bonded energies. Alternatives:
    // virial-style (force·r), central-atom (angle apex, dihedral i+j),
    // or atomic-energy-density methods. Calibration-absorbable for
    // kernel-style models. Per-atom-sum gives the PROTEIN-only
    // bonded contribution; it does NOT reproduce the whole-system .edr
    // term because FullSystemReader filters bonded interactions to
    // protein-slice atoms (src/FullSystemReader.cpp:329 — solvent
    // bonded interactions are dropped).
    grp.createAttribute("split_convention",
        std::string("evenly_among_participating_atoms"));
    grp.createAttribute("split_convention_note",
        std::string("one of several valid attributions; calibration-absorbable"));
    grp.createAttribute("system_scope",
        std::string("protein_slice_only; sum != .edr whole-system term"));

    auto emit_channel = [&](const std::string& name,
                            const std::vector<std::vector<double>>& src) {
        std::vector<double> flat(N * T);
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t t = 0; t < T; ++t) {
                flat[i * T + t] = source_attached_per_frame_[t]
                                ? src[i][t]
                                : kNaN;
            }
        }
        HighFive::DataSpace space({N, T});
        auto ds = grp.createDataSet<double>(name, space);
        ds.write_raw(flat.data());
        ds.createAttribute("units", std::string("kJ/mol"));
    };

    emit_channel("bond",          bond_);
    emit_channel("angle",         angle_);
    emit_channel("urey_bradley",  urey_bradley_);
    emit_channel("proper_dih",    proper_dih_);
    emit_channel("improper_dih",  improper_dih_);
    emit_channel("cmap_dih",      cmap_);
    emit_channel("total",         total_);

    auto fi_ds = grp.createDataSet("frame_indices", frame_indices_);
    fi_ds.createAttribute("units", std::string("frame_index"));
    auto ft_ds = grp.createDataSet("frame_times",   frame_times_);
    ft_ds.createAttribute("units", std::string("ps"));
    auto sa_ds = grp.createDataSet("source_attached_per_frame",
                                   source_attached_per_frame_);
    sa_ds.createAttribute("units", std::string("dimensionless"));
}

}  // namespace nmr

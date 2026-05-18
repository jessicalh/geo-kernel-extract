#include "GromacsEnergyTimeSeriesTrajectoryResult.h"
#include "ProteinConformation.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"
#include "OperationLog.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <typeinfo>

namespace nmr {


std::vector<std::type_index>
GromacsEnergyTimeSeriesTrajectoryResult::Dependencies() const {
    return { std::type_index(typeid(GromacsEnergyResult)) };
}


std::unique_ptr<GromacsEnergyTimeSeriesTrajectoryResult>
GromacsEnergyTimeSeriesTrajectoryResult::Create(const TrajectoryProtein& tp) {
    (void)tp;
    return std::make_unique<GromacsEnergyTimeSeriesTrajectoryResult>();
}


// ── Compute ──────────────────────────────────────────────────────
//
// Each frame, append a copy of the per-frame GromacsEnergyResult::Energy()
// snapshot. The data is whole-system, not per-atom — storage scales as
// O(T) not O(N·T).

void GromacsEnergyTimeSeriesTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;
    const auto& er = conf.Result<GromacsEnergyResult>();
    per_frame_energy_.push_back(er.Energy());
    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    ++n_frames_;
}


// ── Finalize ─────────────────────────────────────────────────────
//
// No per-atom buffer to adopt and no Welford to close out — system-scalar
// shape. Mark finalized so WriteH5Group emits the provenance flag.

void GromacsEnergyTimeSeriesTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                       Trajectory& traj) {
    (void)tp; (void)traj;
    finalized_ = true;

    OperationLog::Info(LogCalcOther,
        "GromacsEnergyTimeSeriesTrajectoryResult::Finalize",
        "captured " + std::to_string(n_frames_) +
        " per-frame GromacsEnergy snapshots");
}


// ── WriteH5Group ─────────────────────────────────────────────────
//
// /trajectory/gromacs_energy_time_series/ — one (T,) scalar dataset per
// logical channel; one (T, 9) dataset for each of the virial and pressure
// 3×3 tensors (laid out XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ matching GromacsEnergy
// struct layout).

void GromacsEnergyTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    (void)tp;
    const std::size_t T = per_frame_energy_.size();
    auto grp = file.createGroup("/trajectory/gromacs_energy_time_series");

    // Group provenance.
    grp.createAttribute("result_name",   Name());
    grp.createAttribute("n_frames",      T);
    grp.createAttribute("finalized",     finalized_);
    grp.createAttribute("units",         std::string("kJ/mol"));
    grp.createAttribute("tensor_layout",
        std::string("XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ"));

    // Helper: emit a single (T,) channel with its units attribute.
    auto emit_scalar = [&](const std::string& name,
                           const std::string& units,
                           double GromacsEnergy::* member) {
        std::vector<double> col(T);
        for (std::size_t f = 0; f < T; ++f)
            col[f] = per_frame_energy_[f].*member;
        auto ds = grp.createDataSet(name, col);
        ds.createAttribute("units", units);
    };

    // ── Electrostatic (kJ/mol) ────────────────────────────────────
    emit_scalar("coulomb_sr",    "kJ/mol", &GromacsEnergy::coulomb_sr);
    emit_scalar("coulomb_recip", "kJ/mol", &GromacsEnergy::coulomb_recip);
    emit_scalar("coulomb_14",    "kJ/mol", &GromacsEnergy::coulomb_14);

    // ── Bonded (kJ/mol) ───────────────────────────────────────────
    emit_scalar("bond",          "kJ/mol", &GromacsEnergy::bond);
    emit_scalar("angle",         "kJ/mol", &GromacsEnergy::angle);
    emit_scalar("urey_bradley",  "kJ/mol", &GromacsEnergy::urey_bradley);
    emit_scalar("proper_dih",    "kJ/mol", &GromacsEnergy::proper_dih);
    emit_scalar("improper_dih",  "kJ/mol", &GromacsEnergy::improper_dih);
    emit_scalar("cmap_dih",      "kJ/mol", &GromacsEnergy::cmap_dih);

    // ── Van der Waals (kJ/mol) ────────────────────────────────────
    emit_scalar("lj_sr",         "kJ/mol", &GromacsEnergy::lj_sr);
    emit_scalar("lj_14",         "kJ/mol", &GromacsEnergy::lj_14);
    emit_scalar("disper_corr",   "kJ/mol", &GromacsEnergy::disper_corr);

    // ── Thermodynamic state ───────────────────────────────────────
    emit_scalar("potential",     "kJ/mol", &GromacsEnergy::potential);
    emit_scalar("kinetic",       "kJ/mol", &GromacsEnergy::kinetic);
    emit_scalar("total_energy",  "kJ/mol", &GromacsEnergy::total_energy);
    emit_scalar("enthalpy",      "kJ/mol", &GromacsEnergy::enthalpy);
    emit_scalar("temperature",   "K",      &GromacsEnergy::temperature);
    emit_scalar("pressure",      "bar",    &GromacsEnergy::pressure);
    emit_scalar("volume",        "nm^3",   &GromacsEnergy::volume);
    emit_scalar("density",       "kg/m^3", &GromacsEnergy::density);

    // ── Box dimensions ────────────────────────────────────────────
    emit_scalar("box_x",         "nm",     &GromacsEnergy::box_x);
    emit_scalar("box_y",         "nm",     &GromacsEnergy::box_y);
    emit_scalar("box_z",         "nm",     &GromacsEnergy::box_z);

    // ── Per-group temperatures ────────────────────────────────────
    emit_scalar("T_protein",     "K",      &GromacsEnergy::T_protein);
    emit_scalar("T_non_protein", "K",      &GromacsEnergy::T_non_protein);

    // ── Virial / pressure tensors (T, 9) ──────────────────────────
    // Layout matches GromacsEnergy::vir / pres struct order:
    // XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ (3×3 row-major). The 3×3 is symmetric
    // for the virial but stored asymmetrically because the .edr file
    // can report off-diagonal terms in some integrators.
    auto emit_tensor = [&](const std::string& name,
                           const std::string& units,
                           double (GromacsEnergy::* arr)[9]) {
        std::vector<double> flat(T * 9);
        for (std::size_t f = 0; f < T; ++f) {
            for (std::size_t k = 0; k < 9; ++k) {
                flat[f * 9 + k] = (per_frame_energy_[f].*arr)[k];
            }
        }
        HighFive::DataSpace space({T, std::size_t(9)});
        auto ds = grp.createDataSet<double>(name, space);
        ds.write_raw(flat.data());
        ds.createAttribute("units", units);
        ds.createAttribute("tensor_layout",
            std::string("XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ"));
    };

    emit_tensor("virial",          "kJ/mol", &GromacsEnergy::vir);
    emit_tensor("pressure_tensor", "bar",    &GromacsEnergy::pres);

    // ── Frame indexing ────────────────────────────────────────────
    auto fi_ds = grp.createDataSet("frame_indices", frame_indices_);
    fi_ds.createAttribute("units", std::string("frame_index"));
    auto ft_ds = grp.createDataSet("frame_times",   frame_times_);
    ft_ds.createAttribute("units", std::string("ps"));
}

}  // namespace nmr

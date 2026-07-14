#include "RingPuckerTimeSeriesTrajectoryResult.h"

#include "LegacyAmberTopology.h"
#include "OperationLog.h"
#include "PlanarGeometryResult.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Ring.h"
#include "TrajectoryProtein.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <typeinfo>

namespace nmr {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

}  // anonymous namespace


std::unique_ptr<RingPuckerTimeSeriesTrajectoryResult>
RingPuckerTimeSeriesTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<RingPuckerTimeSeriesTrajectoryResult>();
    const Protein& protein = tp.ProteinRef();
    const LegacyAmberTopology& topo = protein.LegacyAmber();

    const std::size_t S = topo.SaturatedRingCount();
    const std::size_t A = topo.AromaticRingCount();

    r->pucker_Q_.assign(S, {});
    r->pucker_theta_.assign(S, {});
    r->aromatic_chi2_.assign(A, {});

    r->saturated_parent_residue_index_.assign(S, -1);
    for (std::size_t si = 0; si < S; ++si) {
        r->saturated_parent_residue_index_[si] =
            static_cast<std::int32_t>(topo.SaturatedRingAt(si).parent_residue_index);
    }
    r->aromatic_parent_residue_index_.assign(A, -1);
    for (std::size_t ai = 0; ai < A; ++ai) {
        r->aromatic_parent_residue_index_[ai] =
            static_cast<std::int32_t>(topo.AromaticRingAt(ai).parent_residue_index);
    }
    return r;
}


// ── Compute ──────────────────────────────────────────────────────────
//
// Per frame: if PlanarGeometryResult attached, copy per-ring pucker
// (Q, θ) and per-aromatic-ring χ₂. Production PerFrameExtractionSet
// hard-requires PlanarGeometryResult; the source-absent path is kept
// for custom/synthetic configurations.
//
// CROSS-RESULT READ (read-point): PG.PuckerQ() / PG.PuckerTheta() /
// PG.AromaticChi2() are length-stable per-frame arrays sized by the
// LegacyAmber substrate ring counts. We expect the array sizes to
// match SaturatedRingCount() / AromaticRingCount() unchanged across
// frames (topology is invariant).

void RingPuckerTimeSeriesTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;

    const bool source_attached = force_source_present_for_testing_
        || conf.HasResult<PlanarGeometryResult>();
    source_attached_per_frame_.push_back(source_attached ? 1u : 0u);

    const std::size_t S = pucker_Q_.size();
    const std::size_t A = aromatic_chi2_.size();

    if (source_attached) {
        const auto& pg = conf.Result<PlanarGeometryResult>();
        const auto& Q_arr = pg.PuckerQ();
        const auto& th_arr = pg.PuckerTheta();
        const auto& chi2_arr = pg.AromaticChi2();
        for (std::size_t si = 0; si < S; ++si) {
            pucker_Q_[si].push_back(si < Q_arr.size() ? Q_arr[si] : kNaN);
            pucker_theta_[si].push_back(si < th_arr.size() ? th_arr[si] : kNaN);
        }
        for (std::size_t ai = 0; ai < A; ++ai) {
            aromatic_chi2_[ai].push_back(
                ai < chi2_arr.size() ? chi2_arr[ai] : kNaN);
        }
    } else {
        for (std::size_t si = 0; si < S; ++si) {
            pucker_Q_[si].push_back(kNaN);
            pucker_theta_[si].push_back(kNaN);
        }
        for (std::size_t ai = 0; ai < A; ++ai) {
            aromatic_chi2_[ai].push_back(kNaN);
        }
    }

    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    ++n_frames_;
}


void RingPuckerTimeSeriesTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                    Trajectory& traj) {
    (void)tp; (void)traj;
    finalized_ = true;
    OperationLog::Info(LogCalcOther,
        "RingPuckerTimeSeriesTrajectoryResult::Finalize",
        "finalized across " + std::to_string(n_frames_) +
        " frames, " + std::to_string(pucker_Q_.size()) +
        " saturated rings, " + std::to_string(aromatic_chi2_.size()) +
        " aromatic rings.");
}


void RingPuckerTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    const std::size_t S = pucker_Q_.size();
    const std::size_t A = aromatic_chi2_.size();
    const std::size_t T = n_frames_;
    const std::size_t N = tp.AtomCount();

    std::size_t source_attached_count = 0;
    for (auto v : source_attached_per_frame_) if (v) ++source_attached_count;
    if (source_attached_count == 0) {
        OperationLog::Warn(
            "RingPuckerTimeSeriesTrajectoryResult::WriteH5Group",
            "PlanarGeometryResult was not attached in any of " +
            std::to_string(source_attached_per_frame_.size()) +
            " frames; skipping /trajectory/ring_pucker_time_series/ emission.");
        return;
    }

    auto grp = file.createGroup("/trajectory/ring_pucker_time_series");

    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_saturated_rings", S);
    grp.createAttribute("n_aromatic_rings",  A);
    grp.createAttribute("n_atoms",     N);
    grp.createAttribute("n_frames",    T);
    grp.createAttribute("finalized",   finalized_);
    grp.createAttribute("source_attached_count", source_attached_count);

    grp.createAttribute("pucker_convention", std::string(
        "Cremer-Pople 1975 J. Am. Chem. Soc. 97, 1354. 5-ring "
        "formulation per Eqs 11-16 (single (Q2, theta2) pair). "
        "Implementation lives in PlanarGeometryResult.cpp; sign "
        "convention via R'1 = sum z_j sin, R'2 = sum z_j cos, n = "
        "R'1 x R'2 (canonical Cremer-Pople 1975)."));
    grp.createAttribute("pucker_Q_units",       std::string("Angstrom"));
    grp.createAttribute("pucker_theta_units",   std::string("degrees"));
    grp.createAttribute("pucker_theta_range",   std::string("[0, 360)"));
    grp.createAttribute("pucker_5ring_endvtwist", std::string(
        "theta mod 72 deg gives envelope (E) vs twist (T) endo/exo "
        "classification; reference Cremer-Pople 1975."));
    grp.createAttribute("aromatic_chi2_units", std::string("radians"));
    grp.createAttribute("aromatic_chi2_convention", std::string(
        "IUPAC signed dihedral via atan2(y,x). Ca-Cb-Cg-Cd1 for "
        "PHE/TYR; Ca-Cb-Cg-Nd1 for HIS (Markley et al. 1998 IUPAC "
        "nomenclature, Pure Appl. Chem. 70:117 -- older sources use "
        "Ca-Cb-Cg-Cd2 which differs by ~120 deg since HIS is "
        "asymmetric); Ca-Cb-Cg-Cd1 for TRP. Matches "
        "DihedralTimeSeries.chi[1] for the parent aromatic residue. "
        "Per Akke & Weininger 2023 J. Phys. Chem. B 127, 591: chi2 "
        "IS the canonical ring-flip observable, but the per-frame "
        "value is INSTANTANEOUS -- ring-flip kinetics (slow ~10^1-10^2 "
        "s^-1 via CEST/CPMG; fast averaged) are NOT measurable from "
        "one frame."));
    grp.createAttribute("source", std::string("PlanarGeometryResult"));
    grp.createAttribute("source_attached_policy", std::string(
        "production PerFrameExtractionSet hard-requires "
        "PlanarGeometryResult; source_attached_per_frame is retained "
        "for custom/synthetic configurations. NaN "
        "values within attached frames indicate per-ring degenerate "
        "geometry (incomplete ring or collinear vertices)."));
    grp.createAttribute("saturated_ring_axis", std::string("saturated_ring_index"));
    grp.createAttribute("aromatic_ring_axis",  std::string("aromatic_ring_index"));

    // ── Per-frame metadata (T,) ──────────────────────────────────────
    grp.createDataSet("frame_indices", frame_indices_)
       .createAttribute("units", std::string("frame_index"));
    grp.createDataSet("frame_times",   frame_times_)
       .createAttribute("units", std::string("ps"));
    grp.createDataSet("source_attached_per_frame", source_attached_per_frame_)
       .createAttribute("units", std::string("dimensionless"));

    // ── Per-ring static (S,) and (A,) ────────────────────────────────
    grp.createDataSet("saturated_parent_residue_index",
                      saturated_parent_residue_index_)
       .createAttribute("units", std::string("residue_index"));
    grp.createDataSet("aromatic_parent_residue_index",
                      aromatic_parent_residue_index_)
       .createAttribute("units", std::string("residue_index"));

    // ── Per-saturated-ring per-frame (S, T) float64 ──────────────────
    auto emit_sat_2d = [&](const std::string& name,
                            const std::vector<std::vector<double>>& src,
                            const std::string& units) {
        std::vector<double> flat(S * T, kNaN);
        for (std::size_t si = 0; si < S; ++si) {
            const auto& row = src[si];
            for (std::size_t f = 0; f < T && f < row.size(); ++f) {
                flat[si * T + f] = row[f];
            }
        }
        const std::vector<std::size_t> dims = {S, T};
        HighFive::DataSpace space(dims);
        auto ds = grp.createDataSet<double>(name, space);
        ds.write_raw(flat.data());
        ds.createAttribute("units", units);
    };
    // Only emit saturated-ring datasets if S > 0 (HighFive doesn't like
    // zero-extent datasets).
    if (S > 0) {
        emit_sat_2d("pucker_Q",     pucker_Q_,     "Angstrom");
        emit_sat_2d("pucker_theta", pucker_theta_, "degrees");

        auto q_ds = grp.getDataSet("pucker_Q");
        q_ds.createAttribute("coordinate_frame",
            std::string("intrinsic_cremer_pople_5ring"));
        q_ds.createAttribute("irrep_layout", std::string("0e"));
        q_ds.createAttribute("parity", std::string("even"));
        q_ds.createAttribute("transformation", std::string(
            "exact O(3)- and translation-invariant pucker amplitude: Q'=Q"));
        q_ds.createAttribute("validity", std::string(
            "source_attached_per_frame gates source absence; NaN marks an "
            "incomplete ring or degenerate mean-plane construction"));

        auto theta_ds = grp.getDataSet("pucker_theta");
        theta_ds.createAttribute("coordinate_frame",
            std::string("intrinsic_oriented_cremer_pople_5ring_phase"));
        theta_ds.createAttribute("irrep_layout",
            std::string("none_periodic_phase_with_improper_offset"));
        theta_ds.createAttribute("parity", std::string("mixed"));
        theta_ds.createAttribute("transformation", std::string(
            "translation/proper-rotation invariant phase; under an improper "
            "orthogonal transform theta'=(theta+180 degrees) mod 360"));
        theta_ds.createAttribute("validity", std::string(
            "source_attached_per_frame gates source absence; NaN marks an "
            "incomplete/degenerate ring or pucker_Q < 1e-6 Angstrom; compare "
            "finite values with circular differences"));
    }

    // ── Per-aromatic-ring per-frame (A, T) float64 ───────────────────
    if (A > 0) {
        std::vector<double> flat(A * T, kNaN);
        for (std::size_t ai = 0; ai < A; ++ai) {
            const auto& row = aromatic_chi2_[ai];
            for (std::size_t f = 0; f < T && f < row.size(); ++f) {
                flat[ai * T + f] = row[f];
            }
        }
        const std::vector<std::size_t> dims = {A, T};
        HighFive::DataSpace space(dims);
        auto ds = grp.createDataSet<double>("aromatic_chi2", space);
        ds.write_raw(flat.data());
        ds.createAttribute("units", std::string("radians"));
        ds.createAttribute("coordinate_frame",
            std::string("intrinsic_signed_dihedral"));
        ds.createAttribute("irrep_layout", std::string("0o"));
        ds.createAttribute("parity", std::string("odd"));
        ds.createAttribute("transformation", std::string(
            "wrapped signed-dihedral pseudoscalar: chi2'=det(R) chi2 "
            "modulo 2pi; translation invariant"));
        ds.createAttribute("validity", std::string(
            "source_attached_per_frame gates source absence; NaN marks an "
            "unavailable or geometrically degenerate parent-residue chi2; "
            "compare finite values with circular differences"));
    }

    OperationLog::Info(LogCalcOther,
        "RingPuckerTimeSeriesTrajectoryResult::WriteH5Group",
        "wrote /trajectory/ring_pucker_time_series with " +
        std::to_string(S) + " saturated x " + std::to_string(A) +
        " aromatic x " + std::to_string(T) + " frames");
}


}  // namespace nmr

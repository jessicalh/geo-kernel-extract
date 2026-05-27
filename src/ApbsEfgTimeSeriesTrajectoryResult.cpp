#include "ApbsEfgTimeSeriesTrajectoryResult.h"
#include "ApbsFieldResult.h"
#include "TrajectoryProtein.h"
#include "ProteinConformation.h"
#include "ConformationAtom.h"
#include "OperationLog.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <limits>
#include <typeinfo>

namespace nmr {


std::unique_ptr<ApbsEfgTimeSeriesTrajectoryResult>
ApbsEfgTimeSeriesTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<ApbsEfgTimeSeriesTrajectoryResult>();
    r->per_atom_efg_.assign(tp.AtomCount(), std::vector<SphericalTensor>{});
    return r;
}


// ── Compute ──────────────────────────────────────────────────────
//
// Source-attached gate: in production ApbsFieldResult is
// RequireConformationResult'd in PerFrameExtractionSet, so
// HasResult<ApbsFieldResult>() is always true here. The gate is
// defensive — if a non-standard
// RunConfiguration omits the Require, this TR NaN-fills the absent
// frames and records mask=0 instead of capturing the zero-default
// apbs_efg_spherical (canonical "absent, not faked" — see
// feedback_capture_at_the_boundary).

void ApbsEfgTimeSeriesTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;

    const bool source_present = conf.HasResult<ApbsFieldResult>();
    if (!source_present) {
        OperationLog::Warn(
            "ApbsEfgTimeSeriesTrajectoryResult::Compute",
            "ApbsFieldResult not attached at frame " +
            std::to_string(frame_idx) +
            " — NaN-fill emitted + mask=0. Always-attached policy "
            "requires ApbsFieldResult; the run's PerFrameExtractionSet "
            "must RequireConformationResult(ApbsFieldResult).");
    }

    SphericalTensor nan_tensor{};
    if (!source_present) {
        const double nan_d = std::numeric_limits<double>::quiet_NaN();
        nan_tensor.T0 = nan_d;
        nan_tensor.T1 = {nan_d, nan_d, nan_d};
        nan_tensor.T2 = {nan_d, nan_d, nan_d, nan_d, nan_d};
    }

    const std::size_t N = conf.AtomCount();
    for (std::size_t i = 0; i < N; ++i) {
        per_atom_efg_[i].push_back(
            source_present ? conf.AtomAt(i).apbs_efg_spherical : nan_tensor);
    }
    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    source_attached_per_frame_.push_back(source_present ? 1u : 0u);
    ++n_frames_;
}


// ── Finalize ─────────────────────────────────────────────────────
//
// Transfer growing per-atom buffers into a contiguous atom-major
// DenseBuffer<SphericalTensor> owned by TrajectoryProtein.
//
// Idempotent via bounds-check (see feedback_bounds_check_over_state_flag):
// after the first Finalize swaps per_atom_efg_[i] to empty, src.size()
// is 0, the per-atom skip kicks in, atoms_written remains 0, we return
// without re-Adopting.

void ApbsEfgTimeSeriesTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                  Trajectory& traj) {
    (void)traj;
    const std::size_t N = tp.AtomCount();

    auto buffer = std::make_unique<DenseBuffer<SphericalTensor>>(N, n_frames_);
    std::size_t atoms_written = 0;
    for (std::size_t i = 0; i < N; ++i) {
        const auto& src = per_atom_efg_[i];
        if (src.size() != n_frames_) continue;
        for (std::size_t f = 0; f < n_frames_; ++f) {
            buffer->At(i, f) = src[f];
        }
        std::vector<SphericalTensor>().swap(per_atom_efg_[i]);
        ++atoms_written;
    }
    if (atoms_written == 0) return;

    tp.AdoptDenseBuffer<SphericalTensor>(
        std::move(buffer),
        std::type_index(typeid(ApbsEfgTimeSeriesTrajectoryResult)));

    finalized_ = true;

    OperationLog::Info(LogCalcOther,
        "ApbsEfgTimeSeriesTrajectoryResult::Finalize",
        "transferred (" + std::to_string(N) + " atoms x " +
        std::to_string(n_frames_) +
        " frames) APBS EFG T2 SphericalTensor to tp dense buffer");
}


// ── WriteH5Group ─────────────────────────────────────────────────
//
// Flat (N · T · 5) double array via explicit T2[k] component access.
// T2-only emission: APBS EFG is rank-2 symmetric-traceless after the
// symmetrization and trace projection applied at source, so T0 and T1
// are structurally zero.
//
// The 5-component trailing axis is the e3nn real-spherical
// (l=2, m=-2 ... m=+2) ordering matching SphericalTensor::T2 layout.
// The (irrep_layout, normalization, parity, units) attrs pin the
// convention for downstream Python consumers — parity "2e" because
// the gradient of a polar vector is an even-parity rank-2 tensor.

void ApbsEfgTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    auto* buffer = const_cast<TrajectoryProtein&>(tp)
        .GetDenseBuffer<SphericalTensor>(std::type_index(
            typeid(ApbsEfgTimeSeriesTrajectoryResult)));
    if (!buffer) {
        OperationLog::Warn(
            "ApbsEfgTimeSeriesTrajectoryResult::WriteH5Group",
            "no dense buffer present (Finalize not called?)");
        return;
    }

    const std::size_t N = buffer->AtomCount();
    const std::size_t T = buffer->StridePerAtom();

    auto grp = file.createGroup("/trajectory/apbs_efg_time_series");

    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_atoms",     N);
    grp.createAttribute("n_frames",    T);
    grp.createAttribute("finalized",   finalized_);

    grp.createAttribute("irrep_layout",
        std::string("T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"));
    grp.createAttribute("normalization", std::string("isometric_real_sph"));
    grp.createAttribute("parity",        std::string("2e"));
    grp.createAttribute("units",         std::string("V/Å^2"));
    grp.createAttribute("source", std::string(
        "ApbsFieldResult.apbs_efg_spherical (SphericalTensor, T2 "
        "components 0..4 only). APBS EFG = Hessian of φ from "
        "linearised Poisson-Boltzmann solve (traceless projection "
        "applied at source — ApbsFieldResult.cpp:262-272). T0+T1 "
        "structurally zero by Hessian symmetry + Poisson "
        "tracelessness; 5-component T2-only emission per the "
        "2026-05-18 EFG schema rev (task #166)."));
    grp.createAttribute("source_attached_policy", std::string(
        "always_attached -- ApbsFieldResult is RequireConformationResult'd "
        "in PerFrameExtractionSet (RunConfiguration.cpp:157). Compute's "
        "HasResult<ApbsFieldResult>() gate is defensive and emits "
        "NaN-fill + mask=0 on absent frames per canonical 'absent, "
        "not faked' (OBJECT_MODEL.md Conditional-attach TR discipline)."));

    // (N, T, 5) via explicit T2[k] component access. No reinterpret,
    // no struct-packing assumption. Atom-major:
    // [atom_0_frame_0_T2(-2..+2), atom_0_frame_1_T2(-2..+2), ...,
    //  atom_1_frame_0_T2(-2..+2), ...].
    std::vector<double> flat(N * T * 5);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t t = 0; t < T; ++t) {
            const SphericalTensor& st = buffer->At(i, t);
            const std::size_t base = (i * T + t) * 5;
            st.PackT2(&flat[base]);
        }
    }
    std::vector<std::size_t> dims = {N, T, std::size_t(5)};
    HighFive::DataSpace space(dims);
    auto ds = grp.createDataSet<double>("t2", space);
    ds.write_raw(flat.data());

    grp.createDataSet("frame_indices", frame_indices_);
    grp.createDataSet("frame_times",   frame_times_)
       .createAttribute("units", std::string("ps"));
    grp.createDataSet("source_attached_per_frame", source_attached_per_frame_);
}

}  // namespace nmr

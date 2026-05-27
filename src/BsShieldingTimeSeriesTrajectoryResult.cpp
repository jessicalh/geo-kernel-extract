#include "BsShieldingTimeSeriesTrajectoryResult.h"
#include "BiotSavartResult.h"
#include "TrajectoryProtein.h"
#include "ProteinConformation.h"
#include "ConformationAtom.h"
#include "OperationLog.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <typeinfo>

namespace nmr {


std::vector<std::type_index>
BsShieldingTimeSeriesTrajectoryResult::Dependencies() const {
    return { std::type_index(typeid(BiotSavartResult)) };
}


std::unique_ptr<BsShieldingTimeSeriesTrajectoryResult>
BsShieldingTimeSeriesTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<BsShieldingTimeSeriesTrajectoryResult>();
    r->per_atom_shielding_.assign(tp.AtomCount(),
                                  std::vector<SphericalTensor>{});
    return r;
}


// ── Compute ──────────────────────────────────────────────────────
//
// Append this frame's per-atom SphericalTensor to the growing buffers.
// Records frame_idx + time_ps so WriteH5Group emits the frame list
// that downstream consumers need to align columns with time.

void BsShieldingTimeSeriesTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;
    const std::size_t N = conf.AtomCount();
    for (std::size_t i = 0; i < N; ++i) {
        per_atom_shielding_[i].push_back(
            conf.AtomAt(i).bs_shielding_contribution);
    }
    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    ++n_frames_;
}


// ── Finalize ─────────────────────────────────────────────────────
//
// Transfer growing per-atom buffers into a contiguous atom-major
// DenseBuffer<SphericalTensor> owned by TrajectoryProtein.

void BsShieldingTimeSeriesTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                    Trajectory& traj) {
    (void)traj;
    // Idempotent: a second call would read the per-atom buffers already
    // swapped-empty below and index out of bounds (latent UB). Guard early.
    if (finalized_) return;
    const std::size_t N = tp.AtomCount();

    auto buffer =
        std::make_unique<DenseBuffer<SphericalTensor>>(N, n_frames_);
    for (std::size_t i = 0; i < N; ++i) {
        const auto& src = per_atom_shielding_[i];
        for (std::size_t f = 0; f < n_frames_; ++f) {
            buffer->At(i, f) = src[f];
        }
        std::vector<SphericalTensor>().swap(per_atom_shielding_[i]);
    }

    tp.AdoptDenseBuffer<SphericalTensor>(
        std::move(buffer),
        std::type_index(typeid(BsShieldingTimeSeriesTrajectoryResult)));

    finalized_ = true;

    OperationLog::Info(LogCalcOther,
        "BsShieldingTimeSeriesTrajectoryResult::Finalize",
        "transferred (" + std::to_string(N) + " atoms x " +
        std::to_string(n_frames_) +
        " frames) SphericalTensor time-series to tp dense buffer");
}


// ── WriteH5Group ─────────────────────────────────────────────────
//
// Flat (N · T · 9) double array via explicit component access on each
// SphericalTensor — .T0 / .T1[k] / .T2[k] — so the layout is
// independent of any assumption about struct packing. The 9-component
// trailing axis follows SphericalTensor::PackFull9:
//
//   index 0: T0   (l=0, m=0)
//   index 1: T1_x   (Cartesian Levi-Civita dual)
//   index 2: T1_y
//   index 3: T1_z
//   index 4: T2_{-2}
//   index 5: T2_{-1}
//   index 6: T2_{0}
//   index 7: T2_{+1}
//   index 8: T2_{+2}
//
// The payload order above is explicit; group attributes also record
// normalization, parity, units, and the emitted irrep_layout string.

void BsShieldingTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    auto* buffer = const_cast<TrajectoryProtein&>(tp)
        .GetDenseBuffer<SphericalTensor>(std::type_index(
            typeid(BsShieldingTimeSeriesTrajectoryResult)));
    if (!buffer) {
        OperationLog::Warn(
            "BsShieldingTimeSeriesTrajectoryResult::WriteH5Group",
            "no dense buffer present (Finalize not called?)");
        return;
    }

    const std::size_t N = buffer->AtomCount();
    const std::size_t T = buffer->StridePerAtom();

    auto grp = file.createGroup("/trajectory/bs_shielding_time_series");

    // Provenance header.
    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_atoms",     N);
    grp.createAttribute("n_frames",    T);
    grp.createAttribute("finalized",   finalized_);

    // Schema metadata.
    grp.createAttribute("irrep_layout",
        std::string("T0,T1_m-1,T1_m0,T1_m+1,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"));
    grp.createAttribute("normalization", std::string("isometric_real_sph"));
    grp.createAttribute("parity",        std::string("0e+1o+2e"));
    grp.createAttribute("units",         std::string("ppm"));

    // Flat (N, T, 9) via explicit component access. No reinterpret,
    // no struct-packing assumption. Atom-major: [atom_0_frame_0_T0,
    // ..._T1_x, ..., ..._T2_+2, atom_0_frame_1_T0, ...].
    std::vector<double> flat(N * T * 9);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t t = 0; t < T; ++t) {
            const SphericalTensor& st = buffer->At(i, t);
            const std::size_t base = (i * T + t) * 9;
            st.PackFull9(&flat[base]);
        }
    }

    std::vector<std::size_t> dims = {N, T, std::size_t(9)};
    HighFive::DataSpace space(dims);
    auto ds = grp.createDataSet<double>("xyz", space);
    ds.write_raw(flat.data());

    grp.createDataSet("frame_indices", frame_indices_);
    grp.createDataSet("frame_times",   frame_times_);
}

}  // namespace nmr

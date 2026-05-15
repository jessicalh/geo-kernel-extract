#include "TripeptideNeighborShieldingTimeSeriesTrajectoryResult.h"
#include "TripeptideNeighborShieldingResult.h"
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


std::unique_ptr<TripeptideNeighborShieldingTimeSeriesTrajectoryResult>
TripeptideNeighborShieldingTimeSeriesTrajectoryResult::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<TripeptideNeighborShieldingTimeSeriesTrajectoryResult>();
    r->per_atom_shielding_.assign(tp.AtomCount(),
                                  std::vector<SphericalTensor>{});
    return r;
}


// ── Compute ──────────────────────────────────────────────────────
//
// Append this frame's per-atom SphericalTensor to the growing buffers.
// Records frame_idx + time_ps so WriteH5Group emits the frame list
// that downstream consumers need to align columns with time.

void TripeptideNeighborShieldingTimeSeriesTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)tp; (void)traj;
    // Per-frame source-attached provenance — see analogous block in
    // TripeptideBackboneShieldingTimeSeriesTrajectoryResult::Compute.
    const bool source_attached = force_source_present_for_testing_
        || conf.HasResult<TripeptideNeighborShieldingResult>();
    source_present_per_frame_.push_back(source_attached ? 1u : 0u);

    const std::size_t N = conf.AtomCount();
    for (std::size_t i = 0; i < N; ++i) {
        per_atom_shielding_[i].push_back(
            conf.AtomAt(i).tripeptide_neighbor_shielding_spherical);
    }
    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    ++n_frames_;
}


// ── Finalize ─────────────────────────────────────────────────────
//
// Transfer growing per-atom buffers into a contiguous atom-major
// DenseBuffer<SphericalTensor> owned by TrajectoryProtein.

// Idempotent: the per-atom bounds check (src.size() != n_frames_)
// makes a second call short-circuit naturally — after the first
// Finalize swaps per_atom_shielding_[i] to empty, src.size() is 0,
// the bounds check skips that atom, the buffer comes out empty, and
// we don't AdoptDenseBuffer. The data flow makes idempotency fall
// out without a "have I run yet" state flag.
void TripeptideNeighborShieldingTimeSeriesTrajectoryResult::Finalize(TrajectoryProtein& tp,
                                                                    Trajectory& traj) {
    (void)traj;
    const std::size_t N = tp.AtomCount();

    auto buffer =
        std::make_unique<DenseBuffer<SphericalTensor>>(N, n_frames_);
    std::size_t atoms_written = 0;
    for (std::size_t i = 0; i < N; ++i) {
        const auto& src = per_atom_shielding_[i];
        if (src.size() != n_frames_) continue;  // skip mismatched / post-swap
        for (std::size_t f = 0; f < n_frames_; ++f) {
            buffer->At(i, f) = src[f];
        }
        std::vector<SphericalTensor>().swap(per_atom_shielding_[i]);
        ++atoms_written;
    }

    if (atoms_written == 0) return;  // idempotent no-op on second call

    tp.AdoptDenseBuffer<SphericalTensor>(
        std::move(buffer),
        std::type_index(typeid(TripeptideNeighborShieldingTimeSeriesTrajectoryResult)));

    finalized_ = true;

    OperationLog::Info(LogCalcOther,
        "TripeptideNeighborShieldingTimeSeriesTrajectoryResult::Finalize",
        "transferred (" + std::to_string(N) + " atoms x " +
        std::to_string(n_frames_) +
        " frames) SphericalTensor time-series to tp dense buffer");
}


// ── WriteH5Group ─────────────────────────────────────────────────
//
// Flat (N · T · 9) double array via explicit component access on each
// SphericalTensor — .T0 / .T1[k] / .T2[k] — so the layout is
// independent of any assumption about struct packing. The 9-component
// trailing axis is the lexicographic (l, m) ordering that e3nn, NequIP,
// MACE, and sphericart all expect:
//
//   index 0: T0   (l=0, m=0)
//   index 1: T1_{-1}   (l=1, m=-1)
//   index 2: T1_{0}
//   index 3: T1_{+1}
//   index 4: T2_{-2}
//   index 5: T2_{-1}
//   index 6: T2_{0}
//   index 7: T2_{+1}
//   index 8: T2_{+2}
//
// The irrep_layout, normalization, and parity attributes pin the
// semantics for downstream Python consumers. A future e3nn consumer
// reads these and constructs `Irreps("0e+1o+2e")` directly without
// guessing the convention.

void TripeptideNeighborShieldingTimeSeriesTrajectoryResult::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    auto* buffer = const_cast<TrajectoryProtein&>(tp)
        .GetDenseBuffer<SphericalTensor>(std::type_index(
            typeid(TripeptideNeighborShieldingTimeSeriesTrajectoryResult)));
    if (!buffer) {
        OperationLog::Warn(
            "TripeptideNeighborShieldingTimeSeriesTrajectoryResult::WriteH5Group",
            "no dense buffer present (Finalize not called?)");
        return;
    }

    // "Absent, not faked" — skip emission if source calc never attached.
    std::size_t source_present_count = 0;
    for (auto v : source_present_per_frame_)
        if (v) ++source_present_count;
    if (source_present_count == 0) {
        OperationLog::Warn(
            "TripeptideNeighborShieldingTimeSeriesTrajectoryResult::WriteH5Group",
            "TripeptideNeighborShieldingResult was not attached in any of "
            + std::to_string(source_present_per_frame_.size()) +
            " frames; skipping /trajectory/tripeptide_neighbor_shielding_time_series/ "
            "emission per 'absent, not faked' discipline.");
        return;
    }

    const std::size_t N = buffer->AtomCount();
    const std::size_t T = buffer->StridePerAtom();

    auto grp = file.createGroup("/trajectory/tripeptide_neighbor_shielding_time_series");

    // Provenance header.
    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_atoms",     N);
    grp.createAttribute("n_frames",    T);
    grp.createAttribute("finalized",   finalized_);

    // e3nn-consumable metadata.
    grp.createAttribute("irrep_layout",
        std::string("T0,T1_m-1,T1_m0,T1_m+1,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2"));
    grp.createAttribute("normalization", std::string("isometric_real_sph"));
    grp.createAttribute("parity",        std::string("0e+1o+2e"));
    grp.createAttribute("units",         std::string("ppm"));

    // Flat (N, T, 9) via explicit component access. NaN-fill rows
    // where source wasn't attached so readers can isfinite/isnan
    // discriminate "no measurement" from "real measurement = 0."
    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> flat(N * T * 9);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t t = 0; t < T; ++t) {
            const std::size_t base = (i * T + t) * 9;
            if (t >= source_present_per_frame_.size()
                || source_present_per_frame_[t] == 0) {
                for (std::size_t k = 0; k < 9; ++k) flat[base + k] = kNaN;
                continue;
            }
            const SphericalTensor& st = buffer->At(i, t);
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

    std::vector<std::size_t> dims = {N, T, std::size_t(9)};
    HighFive::DataSpace space(dims);
    auto ds = grp.createDataSet<double>("xyz", space);
    ds.write_raw(flat.data());

    grp.createDataSet("frame_indices", frame_indices_);
    grp.createDataSet("frame_times",   frame_times_);

    // Provenance mask.
    grp.createDataSet("source_attached_per_frame", source_present_per_frame_);
}

}  // namespace nmr

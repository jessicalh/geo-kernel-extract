#include "RingNeighbourhoodTrajectoryStats.h"

#include "CalculatorConfig.h"
#include "GeometryResult.h"
#include "OperationLog.h"
#include "PhysicalConstants.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Ring.h"
#include "SpatialIndexResult.h"
#include "TrajectoryProtein.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <typeinfo>

namespace nmr {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// Singularity threshold for the in-plane azimuth: atoms with rho
// below this OR with degenerate ring vertex-0 geometry return NaN
// for in_plane_angle. Matches the project-wide `MIN_DISTANCE` (0.1 Å)
// 1/r-singularity convention in PhysicalConstants.h; the
// pre-2026-05-21-review value of 1e-12 was 11 orders too tight and
// would have admitted geometrically-meaningless angles from near-axis
// atoms (ρ ≪ 0.1 Å).

}  // namespace

std::vector<std::type_index>
RingNeighbourhoodTrajectoryStats::Dependencies() const {
    return {
        std::type_index(typeid(SpatialIndexResult)),
        std::type_index(typeid(GeometryResult)),
    };
}

std::unique_ptr<RingNeighbourhoodTrajectoryStats>
RingNeighbourhoodTrajectoryStats::Create(const TrajectoryProtein& tp) {
    auto r = std::make_unique<RingNeighbourhoodTrajectoryStats>();
    r->n_atoms_ = tp.AtomCount();
    r->per_atom_ring_list_.assign(r->n_atoms_, {});
    r->per_atom_data_.assign(r->n_atoms_, {});
    return r;
}

void RingNeighbourhoodTrajectoryStats::InitStaticSnapshot_(
        const ProteinConformation& conf,
        const TrajectoryProtein& tp) {
    const auto& spatial = conf.Result<SpatialIndexResult>();
    const double cutoff = CalculatorConfig::Get("ring_current_spatial_cutoff");
    const std::size_t n_aromatic = tp.ProteinRef().RingCount();

    std::size_t r_max = 0;
    for (std::size_t ai = 0; ai < n_atoms_; ++ai) {
        const Vec3 atom_pos = conf.PositionAt(ai);
        auto nearby = spatial.RingsWithinRadius(atom_pos, cutoff);
        // Aromatic-only filter: the spatial index covers aromatic +
        // saturated rings (parallel to conf.ring_geometries); saturated
        // indices >= n_aromatic. ProPyrrolidine is excluded per
        // ring-current physics (no current, no contribution).
        for (std::size_t ri : nearby) {
            if (ri < n_aromatic) {
                per_atom_ring_list_[ai].push_back(ri);
            }
        }
        // Sort ascending by ring index for deterministic emission
        // (the spatial index returns by distance ordering which would
        // shuffle across runs if positions drift identically).
        std::sort(per_atom_ring_list_[ai].begin(),
                   per_atom_ring_list_[ai].end());
        r_max = std::max(r_max, per_atom_ring_list_[ai].size());
    }
    r_per_atom_max_ = r_max;

    OperationLog::Info(
        "RingNeighbourhoodTrajectoryStats::InitStaticSnapshot",
        "n_atoms=" + std::to_string(n_atoms_) +
        " n_aromatic_rings=" + std::to_string(n_aromatic) +
        " r_per_atom_max=" + std::to_string(r_per_atom_max_) +
        " cutoff=" + std::to_string(cutoff) + "A");
}

void RingNeighbourhoodTrajectoryStats::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)traj;

    if (n_frames_ == 0) {
        InitStaticSnapshot_(conf, tp);
    }

    const std::size_t slot_count = r_per_atom_max_ * kChannelCount;

    // r_per_atom_max_ == 0 means no aromatic-ring/atom pairs are in the
    // frame-0 cutoff set. Per-frame slab is zero-sized; we still
    // advance frame counters so the absence is uniform with other
    // TRs' source_attached_per_frame timeline.
    for (std::size_t ai = 0; ai < n_atoms_; ++ai) {
        if (slot_count == 0) continue;

        const Vec3 atom_pos = conf.PositionAt(ai);
        auto& dst = per_atom_data_[ai];
        const std::size_t base = dst.size();
        // NaN-pad the entire R_max * 4 slab; live slots overwritten
        // below. Trailing slots (r >= per_atom_ring_count[ai]) stay NaN.
        dst.resize(base + slot_count, kNaN);

        const auto& rings = per_atom_ring_list_[ai];
        for (std::size_t r_slot = 0; r_slot < rings.size(); ++r_slot) {
            const std::size_t ri = rings[r_slot];
            const RingGeometry& geom = conf.ring_geometries[ri];

            const Vec3 d = atom_pos - geom.center;
            const double dist = d.norm();
            const double z = d.dot(geom.normal);
            const Vec3 rho_vec = d - z * geom.normal;
            const double rho = rho_vec.norm();

            // in_plane_angle: azimuth in ring plane from vertex 0
            // toward atom, in [0, 2*pi). Right-handed about normal.
            // Vertex 0's in-plane direction defines the zero
            // reference; the atom's rho_vec is decomposed onto that
            // basis + the perpendicular `normal x v0_hat`.
            double in_plane_angle = kNaN;
            if (!geom.vertices.empty() && rho > MIN_DISTANCE) {
                const Vec3 v0_to_center = geom.vertices[0] - geom.center;
                const double v0_axial = v0_to_center.dot(geom.normal);
                const Vec3 v0_inplane = v0_to_center - v0_axial * geom.normal;
                const double v0_norm = v0_inplane.norm();
                if (v0_norm > MIN_DISTANCE) {
                    const Vec3 v0_hat = v0_inplane / v0_norm;
                    const Vec3 perp_hat = geom.normal.cross(v0_hat);
                    const double cos_phi = rho_vec.dot(v0_hat) / rho;
                    const double sin_phi = rho_vec.dot(perp_hat) / rho;
                    in_plane_angle = std::atan2(sin_phi, cos_phi);
                    // Wrap [-pi, pi] -> [0, 2*pi). The `< 0` branch maps
                    // negative angles forward; the second branch guards
                    // the FP edge where (-eps) + 2*pi rounds to exactly
                    // 2*pi (we want strict [0, 2*pi) range).
                    if (in_plane_angle < 0.0) in_plane_angle += 2.0 * M_PI;
                    if (in_plane_angle >= 2.0 * M_PI) in_plane_angle = 0.0;
                }
            }

            const std::size_t slot_offset = base + r_slot * kChannelCount;
            dst[slot_offset + 0] = dist;
            dst[slot_offset + 1] = rho;
            dst[slot_offset + 2] = z;
            dst[slot_offset + 3] = in_plane_angle;
        }
    }

    frame_indices_.push_back(frame_idx);
    frame_times_.push_back(time_ps);
    source_attached_per_frame_.push_back(1u);  // always_attached policy
    ++n_frames_;
}

void RingNeighbourhoodTrajectoryStats::Finalize(
        TrajectoryProtein& tp,
        Trajectory& traj) {
    (void)traj;

    // R_max == 0 case: nothing to transfer. Skip the DenseBuffer
    // adoption; WriteH5Group will skip the group entirely.
    if (r_per_atom_max_ == 0) {
        finalized_ = true;
        OperationLog::Info(
            "RingNeighbourhoodTrajectoryStats::Finalize",
            "r_per_atom_max=0 (no aromatic rings in range of any atom); "
            "skipping dense buffer transfer + H5 group emission");
        return;
    }

    // **Idempotency** (codex round 1 2026-05-21 HIGH finding,
    // `feedback_bounds_check_over_state_flag`): the per-atom buffers
    // are swapped EMPTY at the end of the first Finalize call (we
    // move them into the DenseBuffer's storage). A second Finalize
    // call would otherwise allocate a zero-sized DenseBuffer, fail
    // the size-mismatch check per atom, and STILL adopt the empty
    // buffer — overwriting the real data already adopted by tp.
    // Data-flow short-circuit: if any per-atom buffer is empty,
    // a previous Finalize call already ran. Bail.
    if (!per_atom_data_.empty() && per_atom_data_[0].empty()) {
        return;
    }

    const std::size_t stride = n_frames_ * r_per_atom_max_ * kChannelCount;
    auto buffer = std::make_unique<DenseBuffer<double>>(n_atoms_, stride);

    for (std::size_t ai = 0; ai < n_atoms_; ++ai) {
        const auto& src = per_atom_data_[ai];
        if (src.size() != stride) {
            OperationLog::Error(
                "RingNeighbourhoodTrajectoryStats::Finalize",
                "atom " + std::to_string(ai) + " buffer size " +
                std::to_string(src.size()) + " != expected " +
                std::to_string(stride));
            continue;
        }
        double* dst = buffer->AtomSlicePtr(ai);
        std::copy(src.begin(), src.end(), dst);
        std::vector<double>().swap(per_atom_data_[ai]);
    }

    tp.AdoptDenseBuffer<double>(
        std::move(buffer),
        std::type_index(typeid(RingNeighbourhoodTrajectoryStats)));

    finalized_ = true;

    OperationLog::Info(
        "RingNeighbourhoodTrajectoryStats::Finalize",
        "transferred (N=" + std::to_string(n_atoms_) +
        ", T=" + std::to_string(n_frames_) +
        ", R_max=" + std::to_string(r_per_atom_max_) +
        ", 4 channels) to tp dense buffer");
}

void RingNeighbourhoodTrajectoryStats::WriteH5Group(
        const TrajectoryProtein& tp,
        HighFive::File& file) const {
    // Skip the group entirely when there were no aromatic-ring/atom
    // pairs in this protein. "Absent, not faked" -- reader contract:
    // group absence = "no aromatic rings in cutoff of any atom."
    if (r_per_atom_max_ == 0) {
        OperationLog::Info(
            "RingNeighbourhoodTrajectoryStats::WriteH5Group",
            "r_per_atom_max=0; skipping /trajectory/"
            "ring_neighbourhood_trajectory_stats group");
        return;
    }

    auto* buffer = const_cast<TrajectoryProtein&>(tp)
        .GetDenseBuffer<double>(std::type_index(
            typeid(RingNeighbourhoodTrajectoryStats)));
    if (!buffer) {
        OperationLog::Warn(
            "RingNeighbourhoodTrajectoryStats::WriteH5Group",
            "no dense buffer present (Finalize not called?)");
        return;
    }

    const std::size_t N = buffer->AtomCount();
    const std::size_t T = n_frames_;
    const std::size_t R = r_per_atom_max_;

    auto grp = file.createGroup(
        "/trajectory/ring_neighbourhood_trajectory_stats");
    grp.createAttribute("result_name", Name());
    grp.createAttribute("n_atoms", N);
    grp.createAttribute("n_frames", T);
    grp.createAttribute("r_per_atom_max", R);
    grp.createAttribute("finalized", finalized_);
    grp.createAttribute("ring_current_spatial_cutoff_A",
        CalculatorConfig::Get("ring_current_spatial_cutoff"));
    grp.createAttribute("channel_layout",
        std::string("distance,rho,z,in_plane_angle"));
    grp.createAttribute("units",
        std::string("Angstrom,Angstrom,Angstrom,radians"));
    grp.createAttribute("in_plane_angle_range",
        std::string("[0, 2*pi); NaN when atom is on ring axis "
                    "(rho < MIN_DISTANCE = 0.1 Å, "
                    "PhysicalConstants.h:76) OR vertex 0 in-plane "
                    "direction is degenerate"));
    grp.createAttribute("z_sign_convention",
        std::string("z > 0 means atom is on same side as ring "
                    "normal vector"));
    grp.createAttribute("nan_semantics",
        std::string("r_slot >= per_atom_ring_count: unused slot "
                    "(ring_membership_per_atom value is -1 for those "
                    "slots). Within live slots, NaN in geometry "
                    "indicates the in_plane_angle singular case only."));
    grp.createAttribute("aromatic_only",
        std::string("true -- TR reads protein.RingAt(i) (aromatic). "
                    "ProPyrrolidine excluded; emit pucker geometry "
                    "via RingPuckerTimeSeries instead."));
    grp.createAttribute("source_attached_policy",
        std::string("always_attached -- positions present at "
                    "tp.Seed; (atom, ring) pair set frozen at first "
                    "Compute call (frame 0); per-frame geometry "
                    "computed fresh from conf.ring_geometries + "
                    "positions"));
    grp.createAttribute("static_snapshot_origin",
        std::string("frame_0 -- (atom, ring) cutoff-set captured "
                    "from frame-0 positions; subsequent frames "
                    "track geometry of those same pairs and may "
                    "drift past the 15A cutoff (consumer filters "
                    "via the distance channel)"));

    // Dynamic: geometry (N, T, R, 4) double
    std::vector<double> flat(N * T * R * kChannelCount);
    for (std::size_t i = 0; i < N; ++i) {
        const double* src = buffer->AtomSlicePtr(i);
        const std::size_t base = i * T * R * kChannelCount;
        std::copy(src, src + T * R * kChannelCount, flat.begin() + base);
    }
    std::vector<std::size_t> dims = {N, T, R, kChannelCount};
    HighFive::DataSpace space(dims);
    auto ds = grp.createDataSet<double>("geometry", space);
    ds.write_raw(flat.data());

    // Static: ring_membership_per_atom (N, R) int32, -1 sentinel.
    // This is the per-atom counterpart to TopologySidecar's per-vertex
    // ring_membership.npy (substrate truth); ring_membership_per_atom
    // is the per-atom cutoff-set snapshot at frame 0 (per-conformation
    // geometry, not substrate; emitted here rather than in
    // TopologySidecar per PATTERNS Lesson 1 "Protein is identity and
    // topology only" + OBJECT_MODEL spatial-neighbours-on-conformation).
    std::vector<std::int32_t> membership(N * R, -1);
    for (std::size_t i = 0; i < N; ++i) {
        const auto& rings = per_atom_ring_list_[i];
        for (std::size_t r = 0; r < rings.size() && r < R; ++r) {
            membership[i * R + r] = static_cast<std::int32_t>(rings[r]);
        }
    }
    std::vector<std::size_t> mem_dims = {N, R};
    HighFive::DataSpace mem_space(mem_dims);
    auto mem_ds = grp.createDataSet<std::int32_t>(
        "ring_membership_per_atom", mem_space);
    mem_ds.write_raw(membership.data());

    grp.createDataSet("frame_indices", frame_indices_);
    grp.createDataSet("frame_times", frame_times_);
    grp.createDataSet("source_attached_per_frame",
                       source_attached_per_frame_);
}

}  // namespace nmr

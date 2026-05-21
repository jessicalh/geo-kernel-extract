#include "DftPoseCoordinatorTrajectoryResult.h"

#include "ChiRotamerSelectionTrajectoryResult.h"
#include "OperationLog.h"
#include "RmsdSpikeSelectionTrajectoryResult.h"
#include "SelectionRecord.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"

#include <set>
#include <typeinfo>
#include <utility>

namespace nmr {

namespace {

// Resolve residue_index from a record's metadata. Returns -1 if the
// metadata field is absent (RmsdSpike: whole-protein event, no
// residue dimension) or unparseable. Narrow catch list: stoi throws
// invalid_argument on non-numeric, out_of_range on overflow; only
// those two cases legitimately map to "treat as whole-protein."
// Other exception types would be programming errors and should
// propagate.
int RecordResidueIndex(const SelectionRecord& rec) {
    auto it = rec.metadata.find("residue_index");
    if (it == rec.metadata.end()) return -1;
    try {
        return std::stoi(it->second);
    } catch (const std::invalid_argument&) {
        return -1;
    } catch (const std::out_of_range&) {
        return -1;
    }
}

}  // namespace

std::vector<std::type_index>
DftPoseCoordinatorTrajectoryResult::Dependencies() const {
    return {
        typeid(RmsdSpikeSelectionTrajectoryResult),
        typeid(ChiRotamerSelectionTrajectoryResult),
    };
}

std::unique_ptr<DftPoseCoordinatorTrajectoryResult>
DftPoseCoordinatorTrajectoryResult::Create(const TrajectoryProtein& tp) {
    (void)tp;
    return std::make_unique<DftPoseCoordinatorTrajectoryResult>();
}

// Per-frame is a no-op for the reducer — the SelectionBag is filled
// by upstream emitters during their own Compute. We synthesise the
// reduced set once at Finalize.
void DftPoseCoordinatorTrajectoryResult::Compute(
        const ProteinConformation& conf,
        TrajectoryProtein& tp,
        Trajectory& traj,
        std::size_t frame_idx,
        double time_ps) {
    (void)conf; (void)tp; (void)traj;
    (void)frame_idx; (void)time_ps;
}

void DftPoseCoordinatorTrajectoryResult::Finalize(
        TrajectoryProtein& tp,
        Trajectory& traj) {
    (void)tp;

    // Idempotency guard: a second Finalize call would recollect the
    // same upstream records and push duplicates into the SelectionBag
    // (codex round 2 2026-05-21 MEDIUM). Bag-shaped data flow makes
    // the data-flow short-circuit shared-bag-aware; a flag is the
    // right shape here.
    if (finalized_) return;

    // CROSS-RESULT READ: walk both upstream emitter kinds. Iterate in
    // each kind's push order (oldest first); the first record in each
    // (residue_index, ns_bucket) cell wins.
    //
    // **Iterator invalidation guard** (codex round 1 2026-05-21 HIGH
    // finding): `RecordBag::ByKind<T>()` returns pointers INTO the
    // bag's internal vector. If we pushed reduced records while
    // iterating, the vector could reallocate and invalidate the
    // remaining pointers. Two-phase approach: collect reduced
    // records into a local vector, THEN push them all after both
    // iterations complete.
    std::set<std::pair<int, std::size_t>> seen_keys;
    std::vector<SelectionRecord> reduced;

    auto collect = [&](const SelectionRecord* rec, const char* origin) {
        if (!rec) return;
        const int  ri      = RecordResidueIndex(*rec);
        // Bucket by ABSOLUTE elapsed time (codex round 2 HIGH): the
        // earlier `frame_idx / kNsBucketFrames` form was wrong at any
        // TRR cadence other than the 20 ps/frame default because
        // `frame_idx` counts both read AND skipped raw TRR frames.
        // `time_ps` is the authoritative trajectory clock.
        const auto bucket = static_cast<std::size_t>(
            rec->time_ps / kNsBucketPs);
        const auto key     = std::make_pair(ri, bucket);
        if (!seen_keys.insert(key).second) return;  // dup
        std::string reason = std::string("dft_pose_candidate_") +
            origin + "_frame_" + std::to_string(rec->frame_idx);

        SelectionRecord out(
            std::type_index(typeid(DftPoseCoordinatorTrajectoryResult)),
            rec->frame_idx, rec->time_ps,
            std::move(reason),
            rec->metadata);  // preserve upstream metadata
        out.metadata["upstream_kind"] = origin;
        out.metadata["ns_bucket"]     = std::to_string(bucket);
        reduced.push_back(std::move(out));
    };

    {
        const auto& bag = traj.Selections();
        for (const SelectionRecord* rec :
             bag.ByKind<RmsdSpikeSelectionTrajectoryResult>()) {
            collect(rec, "RmsdSpike");
        }
        for (const SelectionRecord* rec :
             bag.ByKind<ChiRotamerSelectionTrajectoryResult>()) {
            collect(rec, "ChiRotamer");
        }
    }  // bag pointers go out of scope; safe to push now.

    auto& mutable_bag = traj.MutableSelections();
    for (auto& out : reduced) {
        mutable_bag.Push(std::move(out));
        ++n_reduced_;
    }

    OperationLog::Info(
        "DftPoseCoordinatorTrajectoryResult::Finalize",
        "deduped " + std::to_string(n_reduced_) +
        " unique (residue, ns_bucket) candidates from RmsdSpike + "
        "ChiRotamer streams");
    finalized_ = true;
}

}  // namespace nmr

#pragma once
//
// RecordBag<Record>: typed container for a stream of records emitted
// during a trajectory run. Push during Compute, query during Finalize
// and at WriteH5.
//
// The same bag pattern serves two scopes:
//
//   Run scope:  Trajectory owns RecordBag<SelectionRecord> for
//               frame-level events (DFT pose candidates, RMSD spikes,
//               etc.). TrajectoryResults push directly to
//               traj.Selections() during their Compute; reducer
//               TrajectoryResults (e.g. DftPoseCoordinator) consume at
//               Finalize via ByKind<...>() queries and push their own
//               reduced set back.
//
//   Atom scope: TrajectoryAtom owns RecordBag<AtomEvent> for per-atom
//               events (rotamer transitions, hbond form/break, ring
//               flips). Same push/query grammar; different payload;
//               orthogonal axis (per-atom instead of per-run).
//
// The Record type must expose three top-level fields:
//
//     std::type_index kind       — discriminator queried against
//     std::size_t     frame_idx  — original frame index
//     double          time_ps    — simulation time at emission
//
// frame_idx and time_ps are explicit fields — not stashed in a
// metadata map — so window queries are direct reads, not string
// lookups.
//
// Query affordances (ByKind, ByKindSinceFrame, ByKindSinceTime,
// MostRecent, CountByKind) are sample shapes meant to document how the
// bag is used. Copy the pattern to add new queries; do not build
// parallel indexes alongside it.
//

#include <cstddef>
#include <typeindex>
#include <vector>

namespace nmr {

template <typename Record>
class RecordBag {
public:
    // ── Primitives ───────────────────────────────────────────────

    void Push(Record r) { records_.push_back(std::move(r)); }

    const std::vector<Record>& All() const { return records_; }

    std::size_t Count() const { return records_.size(); }

    // All distinct kinds present, in first-seen order.
    std::vector<std::type_index> Kinds() const {
        std::vector<std::type_index> out;
        for (const auto& r : records_) {
            bool seen = false;
            for (const auto& k : out) {
                if (k == r.kind) { seen = true; break; }
            }
            if (!seen) out.push_back(r.kind);
        }
        return out;
    }

    // ── Affordance queries ──────────────────────────────────────

    // Records of one kind in push order.
    std::vector<const Record*> ByKind(std::type_index kind) const {
        std::vector<const Record*> out;
        for (const auto& r : records_)
            if (r.kind == kind) out.push_back(&r);
        return out;
    }
    template <typename T>
    std::vector<const Record*> ByKind() const {
        return ByKind(std::type_index(typeid(T)));
    }

    // Records of one kind with frame_idx >= from_frame.
    std::vector<const Record*> ByKindSinceFrame(
            std::type_index kind, std::size_t from_frame) const {
        std::vector<const Record*> out;
        for (const auto& r : records_)
            if (r.kind == kind && r.frame_idx >= from_frame)
                out.push_back(&r);
        return out;
    }
    template <typename T>
    std::vector<const Record*> ByKindSinceFrame(std::size_t from_frame) const {
        return ByKindSinceFrame(std::type_index(typeid(T)), from_frame);
    }

    // Records of one kind with time_ps >= since_time_ps.
    std::vector<const Record*> ByKindSinceTime(
            std::type_index kind, double since_time_ps) const {
        std::vector<const Record*> out;
        for (const auto& r : records_)
            if (r.kind == kind && r.time_ps >= since_time_ps)
                out.push_back(&r);
        return out;
    }
    template <typename T>
    std::vector<const Record*> ByKindSinceTime(double since_time_ps) const {
        return ByKindSinceTime(std::type_index(typeid(T)), since_time_ps);
    }

    // Most recent record of kind by push order; nullptr if none.
    const Record* MostRecent(std::type_index kind) const {
        for (auto it = records_.rbegin(); it != records_.rend(); ++it)
            if (it->kind == kind) return &*it;
        return nullptr;
    }
    template <typename T>
    const Record* MostRecent() const {
        return MostRecent(std::type_index(typeid(T)));
    }

    // Count of records of one kind.
    std::size_t CountByKind(std::type_index kind) const {
        std::size_t n = 0;
        for (const auto& r : records_) if (r.kind == kind) ++n;
        return n;
    }
    template <typename T>
    std::size_t CountByKind() const {
        return CountByKind(std::type_index(typeid(T)));
    }

private:
    std::vector<Record> records_;
};

}  // namespace nmr

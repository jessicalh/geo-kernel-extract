#pragma once
//
// AtomEvent: one entry in a per-atom RecordBag on TrajectoryAtom.
//
// TrajectoryAtom is a bag of three pattern-shaped things:
//   (A) typed struct vectors for known-shape per-source data
//       (e.g. ring_neighbours on ConformationAtom — trajectory-scope
//        parallel to land as RingNeighbourhoodTrajectoryStats);
//   (B) typed accumulator fields with one writer each
//       (e.g. bs_t0_mean, bs_t0_std — what's on TrajectoryAtom today);
//   (C) the per-atom event bag (this record's home) — open-shape
//       events emitted by scan-mode detectors and lifetime/transition
//       accumulators.
//
// The per-atom bag is orthogonal to the run-scope selection bag on
// Trajectory: per-atom queries on one axis, per-run queries on the
// other, same RecordBag<T> verbs.
//
// Fields:
//   emitter    — typeid of the TrajectoryResult that pushed this event.
//   kind       — the emitter's own event-kind discriminator (typed).
//                One emitter typically produces several kinds (a chi
//                rotamer detector emits Chi1Transition, Chi2Transition,
//                Chi3Transition, Chi4Transition — one emitter, four
//                kinds). Queries on the bag discriminate on `kind`;
//                `emitter` is retained so consumers can filter by who
//                if they want.
//   frame_idx  — original frame index (explicit top-level field so
//                the bag can do windowed queries without string
//                lookups into metadata).
//   time_ps    — simulation time of the frame.
//   metadata   — free-form k/v for event-specific extras.
//

#include <cstddef>
#include <map>
#include <string>
#include <typeindex>
#include <utility>

namespace nmr {

struct AtomEvent {
    std::type_index emitter;
    std::type_index kind;
    std::size_t frame_idx = 0;
    double time_ps = 0.0;
    std::map<std::string, std::string> metadata;

    AtomEvent(std::type_index e,
              std::type_index k,
              std::size_t fi,
              double t,
              std::map<std::string, std::string> m = {})
        : emitter(e), kind(k), frame_idx(fi), time_ps(t),
          metadata(std::move(m)) {}
};

}  // namespace nmr

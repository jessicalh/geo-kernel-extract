#pragma once
//
// AtomEvent: one entry in a per-atom RecordBag<AtomEvent> on TrajectoryAtom —
// an open-shape event emitted by a scan-mode detector or a lifetime/
// transition accumulator.
//
// Fields:
//   emitter   — typeid of the TrajectoryResult that pushed the event.
//   kind      — the emitter's own event-kind discriminator. One emitter may
//               produce several kinds (a chi-rotamer detector emits
//               Chi1..Chi4Transition). Queries discriminate on kind; emitter
//               is kept so consumers can also filter by who.
//   frame_idx — original frame index (top-level so the bag can do windowed
//               queries without metadata string lookups).
//   time_ps   — simulation time of the frame.
//   metadata  — free-form k/v for event-specific extras.
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

#include "Relationship.h"

#include "../model/QtAtom.h"
#include "../model/QtProtein.h"

namespace h5reader::rediscover {

Stratum atomsWhere(std::function<bool(const model::QtAtom&)> pred) {
    // Curries the predicate; the body arrives at iteration. (atomsWhere of
    // SURFACE_DESIGN's verb table — generalises the procedural stratum loops.)
    return [pred = std::move(pred)](const Body& body) -> std::vector<std::size_t> {
        std::vector<std::size_t> out;
        if (!body.run.protein) return out;
        const model::QtProtein& p = *body.run.protein;
        for (std::size_t i = 0; i < p.atomCount(); ++i)
            if (pred(p.atom(i))) out.push_back(i);
        return out;
    };
}

SourceSelector slotsBackend() {
    // The frozen-membership ring backend: verbs::ringSlots, wrapped as RawSources.
    // No config to capture (the H5 supplies the cutoff/membership), so the
    // closure captures nothing but answers (atom, frame) -> [src].
    return [](const Body& body, std::size_t atom, std::size_t frame) -> std::vector<RawSource> {
        std::vector<RawSource> out;
        for (const verbs::RingSlot& s : verbs::ringSlots(body, atom, frame)) {
            RawSource r;
            r.kind = SourceKind::Ring;
            r.ring = s;
            out.push_back(r);
        }
        return out;
    };
}

SourceSelector nearBackend(CloudKind cloud, double cutoff) {
    // near(cloud, cutoff) returns the configured (atom, frame) -> [src] closure
    // (cloud + cutoff baked in) — the canonical curried verb of SURFACE_DESIGN.
    return [cloud, cutoff](const Body& body, std::size_t atom,
                           std::size_t frame) -> std::vector<RawSource> {
        std::vector<RawSource> out;
        for (const SourceRef& ref : verbs::near(body, cloud, atom, frame, cutoff)) {
            RawSource r;
            r.kind = cloud == CloudKind::RingCenters
                         ? SourceKind::Ring
                         : (cloud == CloudKind::ChargeSites ? SourceKind::Charge : SourceKind::Bond);
            r.ref = ref;
            out.push_back(r);
        }
        return out;
    };
}

}  // namespace h5reader::rediscover

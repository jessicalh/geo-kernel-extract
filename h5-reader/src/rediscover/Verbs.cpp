#include "Verbs.h"

#include "../io/QtTrajectoryH5.h"
#include "../model/Conformation.h"
#include "../model/QtAtom.h"
#include "../model/QtProtein.h"
#include "../model/QtRing.h"
#include "../model/QtSpecialBuffers.h"
#include "../model/QtTopology.h"

#include <cmath>

namespace h5reader::rediscover {
namespace verbs {

Vec3 pos(const Body& body, std::size_t atom, std::size_t frame) {
    return body.catalog.valueVec3(body, ArrayId::Positions, atom, frame);
}

AtomState at(const Body& body, std::size_t atom, std::size_t frame) {
    AtomState s;
    s.atom = atom;
    s.frame = frame;
    s.pos = pos(body, atom, frame);
    return s;
}

FrameWindow window(const Body& body, std::size_t atom, std::size_t centerFrame,
                   std::size_t before, std::size_t after) {
    return body.idx.temporal.range(atom, centerFrame, before, after);
}

std::vector<SourceRef> near(const Body& body, CloudKind cloud, std::size_t atom,
                            std::size_t frame, double cutoff) {
    return body.idx.spatial.near(cloud, frame, pos(body, atom, frame), cutoff);
}

std::vector<SourceRef> nearPoint(const Body& body, CloudKind cloud, const Vec3& query,
                                 std::size_t frame, double cutoff) {
    return body.idx.spatial.near(cloud, frame, query, cutoff);
}

double value(const Body& body, ArrayId id, std::size_t atom, std::size_t frame, int slot,
             int comp) {
    return body.catalog.value(body, id, atom, frame, slot, comp);
}

Vec3 valueVec3(const Body& body, ArrayId id, std::size_t atom, std::size_t frame) {
    return body.catalog.valueVec3(body, id, atom, frame);
}

model::SphericalTensor valueTensor(const Body& body, ArrayId id, std::size_t atom,
                                   std::size_t frame) {
    return body.catalog.valueTensor(body, id, atom, frame);
}

std::array<double, 5> valueT2(const Body& body, ArrayId id, std::size_t atom, std::size_t frame) {
    return body.catalog.valueT2(body, id, atom, frame);
}

bool present(const Body& body, ArrayId id, std::size_t atom, std::size_t frame) {
    return body.catalog.present(body, id, atom, frame);
}

std::vector<RingSlot> ringSlots(const Body& body, std::size_t atom, std::size_t frame) {
    std::vector<RingSlot> out;
    const io::QtTrajectoryH5* h5 = body.run.h5();
    const model::QtRingNeighbourhoodTimeSeries* rn = h5 ? h5->ringNeighbourhood() : nullptr;
    if (!rn) return out;
    for (std::size_t slot = 0; slot < rn->n_slots; ++slot) {
        const int32_t ringId = rn->ringIndexAt(atom, slot);
        if (ringId < 0) continue;  // unused slot
        const std::array<double, 4> ch = rn->at(atom, frame, slot);
        const double distance = ch[0];
        if (!std::isfinite(distance) || !(distance > 1e-6)) continue;
        out.push_back({ringId, distance, ch[1], ch[2], ch[3]});
    }
    return out;
}

const model::RingGeometry& ringGeom(const Body& body, std::size_t ring, std::size_t frame) {
    return body.idx.ringGeometry.at(ring, frame);
}

std::optional<int32_t> atomOf(const Body& body, const std::vector<int32_t>& scope,
                              const TypedAtomSelector& selector, QString* err_out) {
    return body.idx.typedAtoms.selectUnique(scope, selector, err_out);
}

std::vector<int32_t> selectAll(const Body& body, const std::vector<int32_t>& scope,
                               const TypedAtomSelector& selector) {
    return body.idx.typedAtoms.select(scope, selector);
}

namespace {
std::size_t heavyParent(const model::QtProtein& p, std::size_t atomIdx) {
    const model::QtAtom& a = p.atom(atomIdx);
    return a.parentAtomIndex >= 0 ? static_cast<std::size_t>(a.parentAtomIndex) : atomIdx;
}
}  // namespace

std::vector<int32_t> ringsOf(const Body& body, std::size_t atom) {
    std::vector<int32_t> out;
    if (!body.run.protein) return out;
    const model::QtProtein& p = *body.run.protein;
    const model::QtTopology& topo = p.topology();
    const std::size_t heavy = heavyParent(p, atom);
    for (int memb : topo.ringMembershipsForAtom(heavy)) {
        const model::QtRingMembership& m = topo.ringMembershipAt(static_cast<std::size_t>(memb));
        if (m.ringId < 0) continue;
        if (topo.ringAt(static_cast<std::size_t>(m.ringId)).IsAromatic()) out.push_back(m.ringId);
    }
    return out;
}

std::vector<int32_t> ownRingAtoms(const Body& body, std::size_t atom) {
    std::vector<int32_t> out;
    if (!body.run.protein) return out;
    const model::QtProtein& p = *body.run.protein;
    const model::QtTopology& topo = p.topology();
    for (int32_t ringId : ringsOf(body, atom)) {
        const model::QtRing& ring = topo.ringAt(static_cast<std::size_t>(ringId));
        for (int32_t ra : ring.atomIndices) out.push_back(ra);
    }
    return out;
}

int ownAromaticRing(const Body& body, std::size_t atom) {
    const std::vector<int32_t> rings = ringsOf(body, atom);
    return rings.empty() ? -1 : rings.front();
}

}  // namespace verbs
}  // namespace h5reader::rediscover

#include "RingGeometryCache.h"

#include "RunData.h"

#include "../model/ConformationGeometry.h"
#include "../model/QtProtein.h"

#include <stdexcept>

namespace h5reader::rediscover {

namespace {

model::RingGeometry canonicalGeometry(const model::Conformation& conf, std::size_t ring, std::size_t frame) {
    model::RingGeometry g = model::RingGeometryAt(conf, ring, frame);
    const std::vector<model::Vec3> verts = model::RingVertices(conf, ring, frame);
    if (verts.size() >= 3 && g.normal.norm() > 1e-12) {
        const model::Vec3 canon = (verts[1] - verts[0]).cross(verts[2] - verts[0]);
        if (canon.norm() > 1e-12 && g.normal.dot(canon) < 0.0) g.normal = -g.normal;
    }
    return g;
}

}  // namespace

RingGeometryCache::RingGeometryCache(const RunData& run) {
    if (!run.protein || !run.conformation) return;
    ringCount_ = run.protein->ringCount();
    frameCount_ = run.conformation->frameCount();
    geometries_.reserve(ringCount_ * frameCount_);
    for (std::size_t frame = 0; frame < frameCount_; ++frame)
        for (std::size_t ring = 0; ring < ringCount_; ++ring)
            geometries_.push_back(canonicalGeometry(*run.conformation, ring, frame));
}

const model::RingGeometry& RingGeometryCache::at(std::size_t ring, std::size_t frame) const {
    if (ring >= ringCount_ || frame >= frameCount_)
        throw std::out_of_range("RingGeometryCache::at out of range");
    return geometries_[frame * ringCount_ + ring];
}

}  // namespace h5reader::rediscover

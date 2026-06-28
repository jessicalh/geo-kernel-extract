#include "RingGeometryCache.h"

#include "RunData.h"

#include "../model/ConformationGeometry.h"
#include "../model/QtProtein.h"

#include <stdexcept>

namespace h5reader::rediscover {

RingGeometryCache::RingGeometryCache(const RunData& run) {
    if (!run.protein || !run.conformation) return;
    ringCount_ = run.protein->ringCount();
    frameCount_ = run.conformation->frameCount();
    geometries_.reserve(ringCount_ * frameCount_);
    for (std::size_t frame = 0; frame < frameCount_; ++frame)
        for (std::size_t ring = 0; ring < ringCount_; ++ring)
            geometries_.push_back(model::RingGeometryAt(*run.conformation, ring, frame));
}

const model::RingGeometry& RingGeometryCache::at(std::size_t ring, std::size_t frame) const {
    if (ring >= ringCount_ || frame >= frameCount_)
        throw std::out_of_range("RingGeometryCache::at out of range");
    return geometries_[frame * ringCount_ + ring];
}

}  // namespace h5reader::rediscover

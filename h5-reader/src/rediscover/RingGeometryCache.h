// RingGeometryCache — per-frame ring geometry with canonical normal flipping.

#pragma once

#include "../model/QtRing.h"

#include <cstddef>
#include <vector>

namespace h5reader::rediscover {

class RunData;

class RingGeometryCache {
public:
    RingGeometryCache() = default;
    explicit RingGeometryCache(const RunData& run);

    const model::RingGeometry& at(std::size_t ring, std::size_t frame) const;
    std::size_t ringCount() const { return ringCount_; }
    std::size_t frameCount() const { return frameCount_; }

private:
    std::size_t ringCount_ = 0;
    std::size_t frameCount_ = 0;
    std::vector<model::RingGeometry> geometries_;  // frame-major: frame*ringCount + ring
};

}  // namespace h5reader::rediscover

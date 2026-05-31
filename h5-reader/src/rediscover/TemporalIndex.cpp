#include "TemporalIndex.h"

#include <algorithm>

namespace h5reader::rediscover {

FrameWindow TemporalIndex::range(std::size_t atom, std::size_t centerFrame, std::size_t before,
                                 std::size_t after) const {
    FrameWindow w;
    w.atom = atom;
    w.center = centerFrame;
    if (frameCount_ == 0) return w;
    const std::size_t boundedCenter = std::min(centerFrame, frameCount_ - 1);
    w.begin = boundedCenter > before ? boundedCenter - before : 0;
    w.end = std::min(frameCount_, boundedCenter + after + 1);
    return w;
}

}  // namespace h5reader::rediscover

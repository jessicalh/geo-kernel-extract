// TemporalIndex — no-copy frame window arithmetic over the resident trajectory.

#pragma once

#include <cstddef>

namespace h5reader::rediscover {

struct FrameWindow {
    std::size_t atom = 0;
    std::size_t begin = 0;  // inclusive H5 row
    std::size_t end = 0;    // exclusive H5 row
    std::size_t center = 0;

    std::size_t size() const { return end > begin ? end - begin : 0; }
    bool contains(std::size_t frame) const { return frame >= begin && frame < end; }
};

class TemporalIndex {
public:
    explicit TemporalIndex(std::size_t frameCount = 0) : frameCount_(frameCount) {}

    FrameWindow range(std::size_t atom, std::size_t centerFrame, std::size_t before,
                      std::size_t after = 0) const;
    FrameWindow trailing(std::size_t atom, std::size_t centerFrame, std::size_t width) const {
        return range(atom, centerFrame, width, 0);
    }

    std::size_t frameCount() const { return frameCount_; }

private:
    std::size_t frameCount_ = 0;
};

}  // namespace h5reader::rediscover

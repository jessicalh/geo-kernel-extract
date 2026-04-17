// QtConformation — one trajectory: the H5 file and everything derived
// from its per-frame data.
//
// Holds the shared_ptr<const AnalysisFile> as the data backing for
// every QtFrame. Owns the QtFrame vector. Back-points at the
// QtProtein for topology lookups. Matches the north-star model:
// QtProtein = identity, QtConformation = trajectory (one H5),
// QtFrame = one sampled time point.

#pragma once

#include "QtFrame.h"
#include "analysis_file.h"

#include <cstddef>
#include <memory>
#include <vector>

namespace h5reader::model {

class QtProtein;

// Column offsets into ring_geometry/data's innermost dimension.
// The writer emits {center_x, center_y, center_z, normal_x, normal_y,
// normal_z, radius}; we parse the fields string array at load time to
// discover the actual order so this file does not silently mis-read if
// the extractor ever reorders. Invalid = -1, loader reports via
// ErrorBus and the overlay falls back to zero geometry.
struct RingGeometryLayout {
    int centerX = -1;
    int centerY = -1;
    int centerZ = -1;
    int normalX = -1;
    int normalY = -1;
    int normalZ = -1;
    int radius  = -1;

    bool IsValid() const {
        return centerX >= 0 && centerY >= 0 && centerZ >= 0
            && normalX >= 0 && normalY >= 0 && normalZ >= 0
            && radius  >= 0;
    }
};

class QtConformation {
public:
    QtConformation(const QtProtein*                 protein,
                   std::shared_ptr<const AnalysisFile> h5);

    ~QtConformation() = default;
    QtConformation(const QtConformation&)            = delete;
    QtConformation& operator=(const QtConformation&) = delete;

    // ----- Back-references -----
    const QtProtein*    protein() const { return protein_; }
    const AnalysisFile* h5()      const { return h5_.get(); }

    // ----- Trajectory parameters -----
    size_t frameCount()    const { return frames_.size(); }
    int    stride()        const { return static_cast<int>(h5_->meta.stride); }
    double startTimePicoseconds() const;
    double endTimePicoseconds()   const;

    // ----- Frame access -----
    const QtFrame& frame(size_t t) const { return frames_[t]; }
    const std::vector<QtFrame>& frames() const { return frames_; }

    // ----- Ring count (convenience, matches QtProtein::ringCount) -----
    size_t ringCount() const { return h5_->n_rings; }

    // Ring-geometry column layout. Decoded once at construction from
    // ring_geometry/fields. Invalid if the H5 lacks ring_geometry or
    // uses an unrecognised field name (the loader logs via ErrorBus).
    const RingGeometryLayout& ringGeometryLayout() const { return ringLayout_; }

private:
    const QtProtein*                    protein_;
    std::shared_ptr<const AnalysisFile> h5_;
    std::vector<QtFrame>                frames_;
    RingGeometryLayout                  ringLayout_;
};

}  // namespace h5reader::model

// Per-field presentation limits shared by the catalog, REST, and tests.

#pragma once

#include "DashboardSignal.h"

#include <QLatin1String>

namespace h5reader::model {

// Structural topology tables describe the molecule rather than a per-frame
// signal. Bond lengths remain displayable because they are scalar measurements.
// The AIMNet2 embedding is a model feature vector rather than a readable plot.
inline bool IsDashboardDisplayable(const SignalDescriptor& descriptor) {
    if (descriptor.valueShape == SignalValueShape::Embedding)
        return false;
    if (descriptor.family == QLatin1String("topology")
        && descriptor.valueShape == SignalValueShape::Category)
        return false;
    // Whole-trajectory rollups would be constant across a temporal strip and are
    // therefore not dashboard signals.
    if (descriptor.valueShape == SignalValueShape::RollupMoments)
        return false;
    return true;
}

}  // namespace h5reader::model

#pragma once

#include <QString>

#include <array>

namespace h5reader::model {

struct DisplayModeCapability {
    bool hasVisibleSurface = false;  // dialog: offer as a pickable display mode
    bool buildsPanelWidget = false;  // controller: build an AbstractStripPanel via setOwnedPanels
    bool emitsPanelRef = false;      // panel-model: track with a "panel" sentinel display ref
};

namespace detail {

struct DisplayModeCapabilityRow {
    const char* mode = "";
    DisplayModeCapability capability;
};

}  // namespace detail

DisplayModeCapability DisplayModeCapabilityFor(const QString& mode);

}  // namespace h5reader::model

#include "DisplayModeCapability.h"

#include "VisualizationRegistry.h"

namespace h5reader::model {

DisplayModeCapability DisplayModeCapabilityFor(const QString& mode) {
    return VisualizationRegistry::instance().capabilityForMode(mode);
}

}  // namespace h5reader::model

#include "VisualizationDefinition.h"

#include "TrajectoryFieldAvailability.h"

namespace h5reader::model {

bool VisualizationDescriptorOffersMode(const SignalDescriptor& descriptor,
                                       const QString& modeId) {
    return descriptor.temporalModes.contains(modeId) || descriptor.staticModes.contains(modeId);
}

bool VisualizationDescriptorDataAvailable(const VisualizationContext& ctx,
                                          const SignalDescriptor& descriptor) {
    if (!ctx.availability)
        return true;

    const TrajectoryFieldAvailabilityRecord* descriptorRecord =
        ctx.availability->recordForDescriptor(descriptor.id);
    if (!descriptorRecord)
        return true;

    if (!TrajectoryFieldAvailability::isVisibleState(descriptorRecord->state))
        return false;

    const TrajectoryFieldAvailabilityRecord* storageRecord =
        ctx.availability->recordForStoragePath(descriptor.storagePath);
    return !storageRecord || TrajectoryFieldAvailability::isVisibleState(storageRecord->state);
}

}  // namespace h5reader::model

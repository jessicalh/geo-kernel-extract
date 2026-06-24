#include "NewmanVisualization.h"

namespace h5reader::model {

namespace {
constexpr auto kMode = "static.newman";
}

QString NewmanVisualization::label() const {
    return QStringLiteral("Newman");
}

bool NewmanVisualization::supports(const SignalDescriptor& descriptor) const {
    if (!VisualizationDescriptorOffersMode(descriptor, QString::fromLatin1(kMode)))
        return false;
    // Backbone phi/psi only -- the residue-dihedral angle time series. (chi /
    // arbitrary 4-atom dihedrals are a follow-up; ComputeNewmanProjection is
    // residue-keyed.)
    return descriptor.storagePath == QStringLiteral("/trajectory/dihedral_time_series");
}

bool NewmanVisualization::isAvailable(const VisualizationContext& ctx,
                                      const SignalDescriptor& descriptor) const {
    return ctx.hasTrajectory && VisualizationDescriptorDataAvailable(ctx, descriptor);
}

DisplayModeCapability NewmanVisualization::capability() const {
    return DisplayModeCapability{true, true, true};
}

QStringList NewmanVisualization::legacyModeIds() const {
    return {QString::fromLatin1(kMode)};
}

}  // namespace h5reader::model

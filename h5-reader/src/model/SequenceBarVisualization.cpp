#include "SequenceBarVisualization.h"

namespace h5reader::model {

namespace {
constexpr auto kMode = "static.bar.sequence";
}

QString SequenceBarVisualization::label() const {
    return QStringLiteral("Bar (sequence)");
}

bool SequenceBarVisualization::supports(const SignalDescriptor& descriptor) const {
    if (!VisualizationDescriptorOffersMode(descriptor, QString::fromLatin1(kMode)))
        return false;
    return descriptor.storagePath == QStringLiteral("/trajectory/ired_order_parameters")
        || descriptor.storagePath == QStringLiteral("/trajectory/reorientational_dynamics")
        || descriptor.storagePath == QStringLiteral("/trajectory/dihedral_autocorrelation");
}

bool SequenceBarVisualization::isAvailable(const VisualizationContext& ctx,
                                           const SignalDescriptor& descriptor) const {
    return ctx.hasTrajectory && VisualizationDescriptorDataAvailable(ctx, descriptor);
}

DisplayModeCapability SequenceBarVisualization::capability() const {
    return DisplayModeCapability{true, true, true};
}

QStringList SequenceBarVisualization::legacyModeIds() const {
    return {QString::fromLatin1(kMode)};
}

}  // namespace h5reader::model

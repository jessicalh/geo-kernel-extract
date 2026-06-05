#include "LagCurveVisualization.h"

namespace h5reader::model {

namespace {
constexpr auto kMode = "static.curve.lag.animated";
}

QString LagCurveVisualization::label() const {
    return QStringLiteral("Curve (lag)");
}

bool LagCurveVisualization::supports(const SignalDescriptor& descriptor) const {
    if (!VisualizationDescriptorOffersMode(descriptor, QString::fromLatin1(kMode)))
        return false;
    return descriptor.storagePath == QStringLiteral("/trajectory/kernel_dynamics")
        || descriptor.storagePath == QStringLiteral("/trajectory/reorientational_dynamics")
        || descriptor.storagePath == QStringLiteral("/trajectory/dihedral_autocorrelation");
}

bool LagCurveVisualization::isAvailable(const VisualizationContext& ctx,
                                        const SignalDescriptor& descriptor) const {
    return ctx.hasTrajectory && VisualizationDescriptorDataAvailable(ctx, descriptor);
}

DisplayModeCapability LagCurveVisualization::capability() const {
    return DisplayModeCapability{true, true, true};
}

QStringList LagCurveVisualization::legacyModeIds() const {
    return {QString::fromLatin1(kMode)};
}

}  // namespace h5reader::model

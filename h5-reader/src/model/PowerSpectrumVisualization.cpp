#include "PowerSpectrumVisualization.h"

namespace h5reader::model {

namespace {
constexpr auto kMode = "static.spectrum.power";
}

QString PowerSpectrumVisualization::label() const {
    return QStringLiteral("Power spectrum");
}

bool PowerSpectrumVisualization::supports(const SignalDescriptor& descriptor) const {
    return VisualizationDescriptorOffersMode(descriptor, QString::fromLatin1(kMode))
        && descriptor.storagePath == QStringLiteral("/trajectory/kernel_dynamics");
}

bool PowerSpectrumVisualization::isAvailable(const VisualizationContext& ctx,
                                             const SignalDescriptor& descriptor) const {
    return ctx.hasTrajectory && VisualizationDescriptorDataAvailable(ctx, descriptor);
}

DisplayModeCapability PowerSpectrumVisualization::capability() const {
    return DisplayModeCapability{false, true, true};
}

QStringList PowerSpectrumVisualization::legacyModeIds() const {
    return {QString::fromLatin1(kMode)};
}

}  // namespace h5reader::model

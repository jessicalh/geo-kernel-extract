#include "FixedFrequencyVisualization.h"

namespace h5reader::model {

namespace {
constexpr auto kMode = "static.fixed_freq";
}

QString FixedFrequencyVisualization::label() const {
    return QStringLiteral("Fixed-freq J(w)");
}

bool FixedFrequencyVisualization::supports(const SignalDescriptor& descriptor) const {
    return VisualizationDescriptorOffersMode(descriptor, QString::fromLatin1(kMode))
        && descriptor.storagePath == QStringLiteral("/trajectory/reorientational_dynamics");
}

bool FixedFrequencyVisualization::isAvailable(const VisualizationContext& ctx,
                                              const SignalDescriptor& descriptor) const {
    return ctx.hasTrajectory && VisualizationDescriptorDataAvailable(ctx, descriptor);
}

DisplayModeCapability FixedFrequencyVisualization::capability() const {
    return DisplayModeCapability{true, true, true};
}

QStringList FixedFrequencyVisualization::legacyModeIds() const {
    return {QString::fromLatin1(kMode)};
}

}  // namespace h5reader::model

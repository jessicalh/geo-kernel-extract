#include "ChordCouplingVisualization.h"

namespace h5reader::model {

namespace {
constexpr auto kMode = "static.chord.coupling";
}

QString ChordCouplingVisualization::label() const {
    return QStringLiteral("Chord (coupling)");
}

bool ChordCouplingVisualization::supports(const SignalDescriptor& descriptor) const {
    return VisualizationDescriptorOffersMode(descriptor, QString::fromLatin1(kMode))
        && descriptor.storagePath == QStringLiteral("/trajectory/kernel_coherence");
}

bool ChordCouplingVisualization::isAvailable(const VisualizationContext& ctx,
                                             const SignalDescriptor& descriptor) const {
    return ctx.hasTrajectory && VisualizationDescriptorDataAvailable(ctx, descriptor);
}

DisplayModeCapability ChordCouplingVisualization::capability() const {
    return DisplayModeCapability{true, true, true};
}

QStringList ChordCouplingVisualization::legacyModeIds() const {
    return {QString::fromLatin1(kMode)};
}

}  // namespace h5reader::model

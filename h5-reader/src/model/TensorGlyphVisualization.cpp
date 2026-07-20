#include "TensorGlyphVisualization.h"

namespace h5reader::model {

namespace {
constexpr auto kMode = "static.tensor";
}

QString TensorGlyphVisualization::label() const {
    return QStringLiteral("Tensor glyph");
}

bool TensorGlyphVisualization::supports(const SignalDescriptor& descriptor) const {
    return descriptor.id == QStringLiteral("ml:experimental_shielding_t2")
        && descriptor.sourceKind == SignalSourceKind::ExperimentalShieldingMl
        && descriptor.valueShape == SignalValueShape::EfgT2
        && VisualizationDescriptorOffersMode(descriptor, QString::fromLatin1(kMode));
}

bool TensorGlyphVisualization::isAvailable(const VisualizationContext& ctx,
                                           const SignalDescriptor& descriptor) const {
    return ctx.hasSceneOverlay
        && ctx.tensorGlyphGestureEnabled
        && VisualizationDescriptorDataAvailable(ctx, descriptor);
}

DisplayModeCapability TensorGlyphVisualization::capability() const {
    return DisplayModeCapability{true, false, true};
}

QStringList TensorGlyphVisualization::legacyModeIds() const {
    return {QString::fromLatin1(kMode)};
}

}  // namespace h5reader::model

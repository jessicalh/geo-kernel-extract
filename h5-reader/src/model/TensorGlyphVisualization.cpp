#include "TensorGlyphVisualization.h"

namespace h5reader::model {

namespace {
constexpr auto kMode = "static.tensor";
}

QString TensorGlyphVisualization::label() const {
    return QStringLiteral("Tensor glyph");
}

bool TensorGlyphVisualization::supports(const SignalDescriptor& descriptor) const {
    return VisualizationDescriptorOffersMode(descriptor, QString::fromLatin1(kMode))
        && (descriptor.valueShape == SignalValueShape::Mat3PerRow
            || descriptor.valueShape == SignalValueShape::EfgT2
            || descriptor.valueShape == SignalValueShape::SphericalTensor
            || descriptor.valueShape == SignalValueShape::TensorComponents);
}

bool TensorGlyphVisualization::isAvailable(const VisualizationContext& ctx,
                                           const SignalDescriptor& descriptor) const {
    return ctx.hasSceneOverlay
        && ctx.tensorGlyphGestureEnabled
        && VisualizationDescriptorDataAvailable(ctx, descriptor);
}

DisplayModeCapability TensorGlyphVisualization::capability() const {
    return DisplayModeCapability{false, false, true};
}

QStringList TensorGlyphVisualization::legacyModeIds() const {
    return {QString::fromLatin1(kMode)};
}

}  // namespace h5reader::model

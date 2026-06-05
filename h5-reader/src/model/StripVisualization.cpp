#include "StripVisualization.h"

#include <QStringList>

namespace h5reader::model {

namespace {

bool isStripMode(const QString& modeId) {
    return modeId.startsWith(QStringLiteral("strip."));
}

bool isOfferableStripMode(const QString& modeId) {
    const QString lower = modeId.toLower();
    return lower.startsWith(QStringLiteral("strip."))
        && !lower.contains(QStringLiteral("spectrum"))
        && !lower.contains(QStringLiteral("fft"));
}

bool descriptorHasOfferableStripMode(const SignalDescriptor& descriptor) {
    for (const QString& mode : descriptor.temporalModes) {
        if (isOfferableStripMode(mode))
            return true;
    }
    return false;
}

std::optional<StripComponent> componentForChannel(const ChannelDescriptor& channel) {
    const QString id = channel.id.toLower();
    if (id == QStringLiteral("x"))
        return StripComponent::VectorX;
    if (id == QStringLiteral("y"))
        return StripComponent::VectorY;
    if (id == QStringLiteral("z"))
        return StripComponent::VectorZ;
    if (id == QStringLiteral("magnitude"))
        return StripComponent::VectorMagnitude;
    if (id == QStringLiteral("t0"))
        return StripComponent::TensorT0;
    if (id == QStringLiteral("t1"))
        return StripComponent::TensorT1;
    if (id == QStringLiteral("t2"))
        return StripComponent::TensorT2;
    if (id == QStringLiteral("component"))
        return StripComponent::TensorComponent;
    return std::nullopt;
}

bool modeSelectsVectorComponents(const QString& modeId) {
    return modeId == QStringLiteral("strip.vector.component");
}

}  // namespace

QString StripVisualization::label() const {
    return QStringLiteral("Strip");
}

bool StripVisualization::supports(const SignalDescriptor& descriptor) const {
    return descriptor.temporal && descriptorHasOfferableStripMode(descriptor);
}

bool StripVisualization::isAvailable(const VisualizationContext& ctx,
                                     const SignalDescriptor& descriptor) const {
    if (descriptor.temporal && !ctx.hasTrajectory)
        return false;
    return VisualizationDescriptorDataAvailable(ctx, descriptor);
}

DisplayModeCapability StripVisualization::capability() const {
    return DisplayModeCapability{true, false, false};
}

QStringList StripVisualization::legacyModeIds() const {
    return {};
}

QVector<StripComponent> StripVisualization::componentsFor(const SignalDescriptor& descriptor) const {
    QVector<StripComponent> components;
    if (descriptor.channels.isEmpty()) {
        components.push_back(StripComponent::Auto);
        return components;
    }

    for (const ChannelDescriptor& channel : descriptor.channels) {
        const std::optional<StripComponent> component = componentForChannel(channel);
        if (component && !components.contains(*component))
            components.push_back(*component);
    }
    if (components.isEmpty())
        components.push_back(StripComponent::Auto);
    return components;
}

std::optional<StripComponent> StripComponentForLegacyMode(const QString& modeId) {
    if (modeId == QStringLiteral("strip.vector.magnitude"))
        return StripComponent::VectorMagnitude;
    if (modeId == QStringLiteral("strip.tensor.T0"))
        return StripComponent::TensorT0;
    if (modeId == QStringLiteral("strip.tensor.T1"))
        return StripComponent::TensorT1;
    if (modeId == QStringLiteral("strip.tensor.T2"))
        return StripComponent::TensorT2;
    if (modeId == QStringLiteral("strip.tensor.component"))
        return StripComponent::TensorComponent;
    if (isStripMode(modeId))
        return StripComponent::Auto;
    return std::nullopt;
}

QString LegacyModeIdForStripComponent(StripComponent component) {
    switch (component) {
    case StripComponent::Auto:
        return QStringLiteral("strip.scalar");
    case StripComponent::VectorX:
    case StripComponent::VectorY:
    case StripComponent::VectorZ:
        return QStringLiteral("strip.vector.component");
    case StripComponent::VectorMagnitude:
        return QStringLiteral("strip.vector.magnitude");
    case StripComponent::TensorT0:
        return QStringLiteral("strip.tensor.T0");
    case StripComponent::TensorT1:
        return QStringLiteral("strip.tensor.T1");
    case StripComponent::TensorT2:
        return QStringLiteral("strip.tensor.T2");
    case StripComponent::TensorComponent:
        return QStringLiteral("strip.tensor.component");
    }
    return {};
}

bool StripModeWantsChannel(const SignalDescriptor& descriptor,
                           const QString& modeId,
                           const ChannelDescriptor& channel) {
    if (!isStripMode(modeId))
        return false;

    if (modeSelectsVectorComponents(modeId)) {
        const std::optional<StripComponent> component = componentForChannel(channel);
        return component == StripComponent::VectorX
            || component == StripComponent::VectorY
            || component == StripComponent::VectorZ
            || (!component && channel.id != QStringLiteral("magnitude"));
    }

    const std::optional<StripComponent> requested = StripComponentForLegacyMode(modeId);
    if (!requested || *requested == StripComponent::Auto)
        return true;

    const std::optional<StripComponent> channelComponent = componentForChannel(channel);
    if (channelComponent)
        return *channelComponent == *requested;

    Q_UNUSED(descriptor);
    return true;
}

}  // namespace h5reader::model

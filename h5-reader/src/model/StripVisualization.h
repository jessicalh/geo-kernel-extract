#pragma once

#include "VisualizationDefinition.h"

#include <optional>

namespace h5reader::model {

class StripVisualization final : public VisualizationDefinition {
public:
    VisualizationType type() const override { return VisualizationType::TemporalStrip; }
    QString label() const override;
    DisplaySurface surface() const override { return DisplaySurface::Strip; }
    bool supports(const SignalDescriptor& descriptor) const override;
    bool isAvailable(const VisualizationContext& ctx,
                     const SignalDescriptor& descriptor) const override;
    DisplayModeCapability capability() const override;
    QStringList legacyModeIds() const override;

    QVector<StripComponent> componentsFor(const SignalDescriptor& descriptor) const;
};

std::optional<StripComponent> StripComponentForLegacyMode(const QString& modeId);
QString LegacyModeIdForStripComponent(StripComponent component);
bool StripModeWantsChannel(const SignalDescriptor& descriptor,
                           const QString& modeId,
                           const ChannelDescriptor& channel);

}  // namespace h5reader::model

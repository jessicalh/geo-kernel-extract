#pragma once

#include "VisualizationDefinition.h"

namespace h5reader::model {

class LagCurveVisualization final : public VisualizationDefinition {
public:
    VisualizationType type() const override { return VisualizationType::LagCurve; }
    QString label() const override;
    DisplaySurface surface() const override { return DisplaySurface::Panel; }
    bool supports(const SignalDescriptor& descriptor) const override;
    bool isAvailable(const VisualizationContext& ctx,
                     const SignalDescriptor& descriptor) const override;
    DisplayModeCapability capability() const override;
    QStringList legacyModeIds() const override;
};

}  // namespace h5reader::model

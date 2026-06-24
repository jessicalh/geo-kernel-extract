#pragma once

#include "VisualizationDefinition.h"

namespace h5reader::model {

// The Newman display form for the residue-dihedral metric: a panel of phi/psi
// d'Arsonval dials (live needle + whole-trajectory occupancy haze) instead of a
// line strip. Panel surface, like SequenceBar.
class NewmanVisualization final : public VisualizationDefinition {
public:
    VisualizationType type() const override { return VisualizationType::Newman; }
    QString label() const override;
    DisplaySurface surface() const override { return DisplaySurface::Panel; }
    bool supports(const SignalDescriptor& descriptor) const override;
    bool isAvailable(const VisualizationContext& ctx,
                     const SignalDescriptor& descriptor) const override;
    DisplayModeCapability capability() const override;
    QStringList legacyModeIds() const override;
};

}  // namespace h5reader::model

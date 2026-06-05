#pragma once

#include "VisualizationDefinition.h"

#include <QStringList>
#include <QVector>

#include <memory>
#include <vector>

namespace h5reader::model {

class TrajectorySignalCatalog;

class VisualizationRegistry {
public:
    // Non-QObject, immutable-after-construction registry. Construct once on
    // first GUI-thread use; all later reads are const and re-entrant.
    static const VisualizationRegistry& instance();

    QVector<const VisualizationDefinition*> definitions() const;
    QVector<const VisualizationDefinition*> supporting(const SignalDescriptor& descriptor) const;
    QVector<const VisualizationDefinition*> visibleOfferable(const VisualizationContext& ctx,
                                                             const SignalDescriptor& descriptor) const;
    QVector<const VisualizationDefinition*> trackedButHidden(const SignalDescriptor& descriptor) const;
    DisplayModeCapability capabilityForMode(const QString& modeId) const;
    const VisualizationDefinition* definitionForMode(const QString& modeId) const;
    QStringList unresolvedStaticModes(const TrajectorySignalCatalog& catalog) const;

private:
    VisualizationRegistry();

    std::vector<std::unique_ptr<VisualizationDefinition>> defs_;
};

}  // namespace h5reader::model

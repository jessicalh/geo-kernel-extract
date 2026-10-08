#pragma once

#include "../model/MetricGlossary.h"

#include <QStringList>
#include <optional>

namespace h5reader::app {

std::optional<model::MetricGlossaryEntry> AtomInspectorGlossary(
    const QStringList& fieldPath, const QString& shieldingSource);

}  // namespace h5reader::app

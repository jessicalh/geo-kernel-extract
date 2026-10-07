#pragma once

#include "DashboardSignal.h"

#include <QString>

#include <optional>

namespace h5reader::model {

struct MetricGlossaryEntry {
    QString meaning;
    QString calculation;
    QString origin;
};

std::optional<MetricGlossaryEntry> MetricGlossaryFor(const SignalDescriptor& descriptor);

}  // namespace h5reader::model

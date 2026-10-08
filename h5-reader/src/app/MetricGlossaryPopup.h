#pragma once

#include <QPoint>

class QWidget;
class QString;

namespace h5reader::model {
struct SignalDescriptor;
struct MetricGlossaryEntry;
}

namespace h5reader::app {

void ShowMetricGlossaryPopup(const QString& title,
                             const model::MetricGlossaryEntry& entry,
                             const QPoint& globalPosition,
                             QWidget* parent);

void ShowMetricGlossaryPopup(const model::SignalDescriptor& descriptor,
                             const QPoint& globalPosition,
                             QWidget* parent);

}  // namespace h5reader::app

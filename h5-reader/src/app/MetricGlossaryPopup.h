#pragma once

#include <QPoint>

class QWidget;

namespace h5reader::model {
struct SignalDescriptor;
}

namespace h5reader::app {

void ShowMetricGlossaryPopup(const model::SignalDescriptor& descriptor,
                             const QPoint& globalPosition,
                             QWidget* parent);

}  // namespace h5reader::app

// DashboardLogging — shared category for selected-metric/dashboard transitions.
//
// The dashboard transition trace spans model, controller, dock, and window
// code. Keep the category declared once so every transition lands under the
// same structured logger category.

#pragma once

#include <QLoggingCategory>

namespace h5reader::diagnostics {

Q_DECLARE_LOGGING_CATEGORY(cDash)

}  // namespace h5reader::diagnostics

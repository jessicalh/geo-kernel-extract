#include "QtConformation.h"
#include "QtProtein.h"

#include "../diagnostics/ErrorBus.h"

#include <QString>
#include <QStringList>

namespace h5reader::model {

namespace {

// Parse ring_geometry/fields into an index map. Known field names are
// center_{x,y,z}, normal_{x,y,z}, radius — the writer's current output.
// Unknown names are ignored (caller's IsValid() catches missing ones).
RingGeometryLayout ParseRingGeometryLayout(const AnalysisFile& h5) {
    RingGeometryLayout layout;
    const auto& fields = h5.ring_geometry.fields;
    for (int i = 0; i < static_cast<int>(fields.size()); ++i) {
        const auto& name = fields[i];
        if      (name == "center_x") layout.centerX = i;
        else if (name == "center_y") layout.centerY = i;
        else if (name == "center_z") layout.centerZ = i;
        else if (name == "normal_x") layout.normalX = i;
        else if (name == "normal_y") layout.normalY = i;
        else if (name == "normal_z") layout.normalZ = i;
        else if (name == "radius")   layout.radius  = i;
    }
    return layout;
}

QString JoinFields(const std::vector<std::string>& fields) {
    QStringList parts;
    parts.reserve(static_cast<int>(fields.size()));
    for (const auto& s : fields) parts << QString::fromStdString(s);
    return parts.join(QStringLiteral(", "));
}

}  // namespace

QtConformation::QtConformation(const QtProtein*                    protein,
                               std::shared_ptr<const AnalysisFile> h5)
    : protein_(protein), h5_(std::move(h5)) {
    ringLayout_ = ParseRingGeometryLayout(*h5_);

    // Report if the writer emitted ring_geometry but its column labels
    // don't match what overlays need. Ring-scoped overlays (ring polygon,
    // ring-current glyphs, butterfly isosurfaces) will silently fall
    // back to zero geometry if this is invalid — that used to be an
    // unexplained "rings don't render" without a log line. Warning is
    // the right severity: rendering is degraded, not a user action failure.
    // No rings is fine; we only flag the label-mismatch case.
    if (h5_->n_rings > 0 && !ringLayout_.IsValid()) {
        h5reader::diagnostics::ErrorBus::Report(
            h5reader::diagnostics::Severity::Warning,
            QStringLiteral("QtConformation"),
            QStringLiteral("ring_geometry column layout is incomplete — "
                           "ring-scoped overlays will render zero geometry"),
            QStringLiteral("n_rings=%1; fields=[%2]")
                .arg(h5_->n_rings)
                .arg(JoinFields(h5_->ring_geometry.fields)));
    }

    const size_t T = h5_->n_frames;
    frames_.reserve(T);
    for (size_t t = 0; t < T; ++t) {
        frames_.emplace_back(this, h5_.get(), t);
    }
}

double QtConformation::startTimePicoseconds() const {
    return h5_->meta.frame_times.empty() ? 0.0 : h5_->meta.frame_times.front();
}

double QtConformation::endTimePicoseconds() const {
    return h5_->meta.frame_times.empty() ? 0.0 : h5_->meta.frame_times.back();
}

}  // namespace h5reader::model

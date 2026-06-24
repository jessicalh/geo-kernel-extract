#include "NewmanPanel.h"

#include "../model/Conformation.h"
#include "../model/QtProtein.h"

#include <QPainter>
#include <QRectF>

#include <utility>

namespace h5reader::app {

NewmanPanel::NewmanPanel(QString title,
                         const model::QtProtein* protein,
                         const model::Conformation* conf,
                         std::size_t residue,
                         NewmanOccupancy occPhi,
                         NewmanOccupancy occPsi,
                         model::SignalBinding revealBinding)
    : title_(std::move(title)),
      protein_(protein),
      conf_(conf),
      residue_(residue),
      occPhi_(std::move(occPhi)),
      occPsi_(std::move(occPsi)),
      reveal_(std::move(revealBinding)),
      hasReveal_(!reveal_.descriptorId.isEmpty()) {}

void NewmanPanel::paint(QPainter& p, const PanelGeometry& geometry,
                        const PaintContext& context) const {
    paintPanelBackground(p, geometry);
    // Title only; each dial draws its own torsion readout.
    paintHeader(p, geometry, HeaderText{title_, QString()}, 0, kStripText, hasReveal_, false);

    const std::size_t frame =
        context.currentFrame >= 0 ? static_cast<std::size_t>(context.currentFrame) : 0u;
    const NewmanProjection phi = (protein_ && conf_)
        ? ComputeNewmanProjection(*protein_, *conf_, residue_, frame, NewmanKind::Phi)
        : NewmanProjection{};
    const NewmanProjection psi = (protein_ && conf_)
        ? ComputeNewmanProjection(*protein_, *conf_, residue_, frame, NewmanKind::Psi)
        : NewmanProjection{};

    // phi on the left, psi on the right -- the residue's backbone conformation.
    const QRectF plot = geometry.plot;
    const double halfW = plot.width() / 2.0;
    PaintNewman(p, QRectF(plot.left(), plot.top(), halfW, plot.height()),
                phi, occPhi_, /*darkBg=*/true);
    PaintNewman(p, QRectF(plot.left() + halfW, plot.top(), halfW, plot.height()),
                psi, occPsi_, /*darkBg=*/true);
}

PanelDisplayData NewmanPanel::displayData() const {
    PanelDisplayData d;
    d.kind = QStringLiteral("newman");
    d.title = title_;
    d.descriptorId = reveal_.descriptorId;
    d.seriesCount = 2;  // phi + psi
    if (protein_ && conf_) {
        for (NewmanKind k : {NewmanKind::Phi, NewmanKind::Psi}) {
            const NewmanProjection pr = ComputeNewmanProjection(*protein_, *conf_, residue_, 0, k);
            if (pr.valid)
                d.note(pr.torsionDeg);
        }
    }
    for (float v : occPhi_.bins) d.note(static_cast<double>(v));
    for (float v : occPsi_.bins) d.note(static_cast<double>(v));
    return d;
}

}  // namespace h5reader::app

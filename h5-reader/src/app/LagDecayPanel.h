// LagDecayPanel — autocorrelation curve over lag, animated by playback.
//
// Renders one polyline per kernel channel for a fixed atom row (same
// data shape as PowerSpectrumPanel; different x-axis: lag time). The
// "animated" twist: a vertical cursor walks across the lag axis in
// step with the PaintContext.currentFrame so the viewer sees the
// curve relax as playback advances. The curve itself is static (no
// per-frame data change).

#pragma once

#include "AbstractStripPanel.h"

#include "../model/DashboardSignal.h"
#include "../model/QtPerAtomChannelBuffers.h"
#include "../model/SignalDictionary.h"

#include <QColor>
#include <QString>

#include <cstddef>
#include <memory>

namespace h5reader::app {

class LagDecayPanel final : public AbstractStripPanel {
public:
    // Borrowing constructor: data lifetime managed elsewhere (KD case).
    LagDecayPanel(QString label,
                  const model::QtPerAtomChannelCurve* data,
                  std::size_t atomRow,
                  model::SignalBinding revealBinding);

    // Owning constructor: panel takes ownership of a synthesized
    // single-row, single-channel view (ReorientationalDynamics
    // adapts its per-bond-vector curve into this shape).
    LagDecayPanel(QString label,
                  std::unique_ptr<model::QtPerAtomChannelCurve> ownedData,
                  std::size_t atomRow,
                  model::SignalBinding revealBinding);

    bool hasRevealBinding() const override { return haveBinding_; }
    model::SignalBinding revealBinding() const override { return reveal_; }

    void paint(QPainter& p,
               const PanelGeometry& geometry,
               const PaintContext& context) const override;

    bool plotContains(const PanelGeometry& geometry, const QPoint& pos) const override {
        return geometry.plot.contains(pos);
    }

private:
    QString label_;
    const model::QtPerAtomChannelCurve* data_ = nullptr;
    std::unique_ptr<model::QtPerAtomChannelCurve> owned_;
    std::size_t atomRow_ = 0;
    model::SignalBinding reveal_;
    bool haveBinding_ = false;
};

}  // namespace h5reader::app

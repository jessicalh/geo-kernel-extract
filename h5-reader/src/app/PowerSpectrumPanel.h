// PowerSpectrumPanel — continuous line plot of power vs frequency,
// multi-channel overlay. Renders one polyline per kernel channel
// (KernelDynamics emits 13) for a fixed atom row.
//
// Static panel: no per-frame redraw. Click anywhere in the plot →
// reveals the bound atom (single anchor; the panel is per-atom).

#pragma once

#include "AbstractStripPanel.h"

#include "../model/DashboardSignal.h"
#include "../model/QtPerAtomChannelBuffers.h"
#include "../model/SignalDictionary.h"

#include <QColor>
#include <QString>

#include <cstddef>
#include <vector>

namespace h5reader::app {

class PowerSpectrumPanel final : public AbstractStripPanel {
public:
    PowerSpectrumPanel(QString label,
                       const model::QtPerAtomChannelCurve* data,
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
    QColor channelColor(std::size_t channelIndex) const;

    QString label_;
    const model::QtPerAtomChannelCurve* data_ = nullptr;
    std::size_t atomRow_ = 0;
    model::SignalBinding reveal_;
    bool haveBinding_ = false;
};

}  // namespace h5reader::app

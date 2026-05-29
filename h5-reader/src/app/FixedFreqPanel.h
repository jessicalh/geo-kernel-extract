// FixedFreqPanel — discrete-frequency line plot for ReorientationalDynamics
// J(ω) (5 KTB Larmor combinations). Renders one trace per selected bond
// vector as 5 log-x markers connected by a polyline.
//
// Distinct from PowerSpectrumPanel: PowerSpectrumPanel walks a dense
// frequency grid (kernel PSD); FixedFreqPanel exposes a sparse,
// externally-fixed grid. Bonded to a single bond-vector row via the
// reveal binding; the panel owns its data view.
//
// L-3b (2026-05-29).

#pragma once

#include "AbstractStripPanel.h"

#include "../model/QtBondVectorBuffers.h"
#include "../model/SignalDictionary.h"

#include <QColor>
#include <QString>

#include <cstddef>
#include <memory>

namespace h5reader::app {

class FixedFreqPanel final : public AbstractStripPanel {
public:
    // The panel owns a private copy of the J(ω) row + frequency axis
    // by taking the QtPerBondVectorFixedFreqBlock view by unique_ptr.
    // The view carries only the resolved row's payload (5 values + 5
    // frequencies); the controller builder slices the buffer down to
    // this minimal shape so the panel never holds a pointer back into
    // a producer-owned vector that could outlive the rebuild.
    FixedFreqPanel(QString label,
                   std::unique_ptr<model::QtPerBondVectorFixedFreqBlock> data,
                   std::size_t row,
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
    std::unique_ptr<model::QtPerBondVectorFixedFreqBlock> data_;
    std::size_t row_ = 0;
    model::SignalBinding reveal_;
    bool haveBinding_ = false;
};

}  // namespace h5reader::app

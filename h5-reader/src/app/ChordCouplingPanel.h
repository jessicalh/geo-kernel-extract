// ChordCouplingPanel — chord/arc diagram of a per-atom NxN correlation
// matrix (KernelCoherence's 13×13 Pearson). Channel labels arranged
// around a circle; arcs between coupled pairs with thickness=|r| and
// hue encoding sign. Heatmap was explicitly rejected (per execution
// rules) — the chord idiom reads at a glance.

#pragma once

#include "AbstractStripPanel.h"

#include "../model/DashboardSignal.h"
#include "../model/QtPerAtomChannelBuffers.h"
#include "../model/SignalDictionary.h"

#include <QString>

#include <cstddef>

namespace h5reader::app {

class ChordCouplingPanel final : public AbstractStripPanel {
public:
    ChordCouplingPanel(QString label,
                       const model::QtPerAtomMatrix* matrix,
                       std::size_t atomRow,
                       double threshold,
                       model::SignalBinding revealBinding);

    bool hasRevealBinding() const override { return haveBinding_; }
    model::SignalBinding revealBinding() const override { return reveal_; }

    void paint(QPainter& p,
               const PanelGeometry& geometry,
               const PaintContext& context) const override;

    // Square panel: chord reads best at 1:1 aspect.
    std::optional<double> preferredAspect() const override { return 1.0; }

private:
    QString label_;
    const model::QtPerAtomMatrix* matrix_ = nullptr;
    std::size_t atomRow_ = 0;
    double threshold_ = 0.2;
    model::SignalBinding reveal_;
    bool haveBinding_ = false;
};

}  // namespace h5reader::app

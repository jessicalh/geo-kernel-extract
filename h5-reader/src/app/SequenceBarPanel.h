// SequenceBarPanel — per-row scalar bars laid out along a sequence axis.
//
// The NMR-relaxation idiom: S² (or τ_e, R1, R2, NOE, …) along the
// residue sequence as a bar chart. Used by 4 of the 5 new TRs
// (IRedOrderParameter, ReorientationalDynamics, DihedralAutocorrelation,
// and ReorientationalDynamics's rate panels). v1 is hosted by IRed for
// the BondVector axis; the same panel handles per-residue and
// per-bond-vector inputs because the rows carry their own residue
// indices.
//
// The panel is decoupled from any specific TR buffer type — callers
// provide a flat vector of (residue_index, value, optional kind) rows
// and a binding-builder lambda that turns a clicked bar's row index
// into a SignalBinding (so the scene can focus on the residue / atom
// pair the bar represents).

#pragma once

#include "AbstractStripPanel.h"

#include "../model/DashboardSignal.h"
#include "../model/SignalDictionary.h"

#include <QColor>
#include <QString>

#include <cstdint>
#include <functional>
#include <optional>
#include <vector>

namespace h5reader::app {

struct SequenceBarRow {
    std::int32_t residue_index = 0;
    double value = 0.0;
    // Vector-kind discriminator (0 = no discriminator, 1=NH, 2=CaHa, 3=CO);
    // used to colour-stripe bars when a single residue has multiple
    // bond-vector rows on the same panel.
    std::uint8_t kind = 0;
};

class SequenceBarPanel final : public AbstractStripPanel {
public:
    using BindingForRow = std::function<model::SignalBinding(std::size_t row)>;

    SequenceBarPanel(QString label,
                     QString unit,
                     std::vector<SequenceBarRow> rows,
                     BindingForRow bindingForRow,
                     QColor barColor = QColor(115, 229, 214),
                     std::optional<double> yMin = std::nullopt,
                     std::optional<double> yMax = std::nullopt);

    bool hasRevealBinding() const override { return false; }
    model::SignalBinding revealBinding() const override { return {}; }

    void paint(QPainter& p,
               const PanelGeometry& geometry,
               const PaintContext& context) const override;

    bool plotContains(const PanelGeometry& geometry, const QPoint& pos) const override {
        return geometry.plot.contains(pos);
    }

    std::optional<model::SignalBinding>
        mousePressInPlot(QMouseEvent* event, const PanelGeometry& geometry) override;

    QString tooltipLine(int frame) const override;

private:
    std::optional<std::size_t> rowAtX(double x, const QRectF& plot) const;

    QString label_;
    QString unit_;
    std::vector<SequenceBarRow> rows_;
    BindingForRow bindingForRow_;
    QColor barColor_;
    std::optional<double> yMin_;
    std::optional<double> yMax_;
    std::int32_t residueMin_ = 0;
    std::int32_t residueMax_ = 0;
};

}  // namespace h5reader::app

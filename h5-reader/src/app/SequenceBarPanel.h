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

    // L-4 (2026-05-29): multi-channel overlay (auto-compose). One
    // SequenceBarPanel can carry a primary series plus N overlays.
    // The auto-compose path in DashboardDisplayController bundles
    // multiple Reorient scalar signals (s2/tau_e/r1/r2/noe) with the
    // same `static.bar.sequence` mode into ONE panel — the primary
    // owns the left y-axis; each overlay paints with its own colour,
    // sub-slotted across the residue tick, and (if units differ from
    // primary) renders y-axis labels on the right margin scaled to
    // the overlay's own value range. Overlay rows respect the same
    // kindOffsetFraction sub-slot machinery as the primary.
    struct OverlaySeries {
        std::vector<SequenceBarRow> rows;
        QColor color;
        QString label;
        QString unit;
        // Computed once at addOverlay time so paint() doesn't re-walk
        // the rows. nullopt → empty rows, the overlay still occupies
        // a slot but draws nothing.
        std::optional<double> valueMin;
        std::optional<double> valueMax;
    };
    void addOverlay(std::vector<SequenceBarRow> rows,
                    QColor color,
                    QString label,
                    QString unit);
    std::size_t overlayCount() const { return overlays_.size(); }

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
    std::vector<OverlaySeries> overlays_;
};

}  // namespace h5reader::app

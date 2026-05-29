#include "SequenceBarPanel.h"

#include <QFontMetricsF>
#include <QMouseEvent>
#include <QPainter>
#include <QPoint>
#include <QRectF>
#include <QString>

#include <algorithm>
#include <cmath>
#include <utility>

namespace h5reader::app {

SequenceBarPanel::SequenceBarPanel(QString label,
                                   QString unit,
                                   std::vector<SequenceBarRow> rows,
                                   BindingForRow bindingForRow,
                                   QColor barColor,
                                   std::optional<double> yMin,
                                   std::optional<double> yMax)
    : label_(std::move(label)),
      unit_(std::move(unit)),
      rows_(std::move(rows)),
      bindingForRow_(std::move(bindingForRow)),
      barColor_(barColor),
      yMin_(yMin),
      yMax_(yMax)
{
    if (!rows_.empty()) {
        residueMin_ = rows_.front().residue_index;
        residueMax_ = residueMin_;
        for (const SequenceBarRow& row : rows_) {
            residueMin_ = std::min(residueMin_, row.residue_index);
            residueMax_ = std::max(residueMax_, row.residue_index);
        }
    }
}

// Per-kind sub-slot offset within a residue slot. Reorient packs up
// to 3 rows per residue (NH=1, CaHa=2, CO=3). Bars at distinct kinds
// fan out across the slot; kind=0 (no discriminator, used by iRED
// and Dihedral) centres on the residue tick.
static double kindOffsetFraction(std::uint8_t kind)
{
    switch (kind) {
    case 1: return -0.27;   // NH   — left
    case 2: return  0.00;   // CaHa — centre
    case 3: return  0.27;   // CO   — right
    default: return 0.00;
    }
}

std::optional<std::size_t>
SequenceBarPanel::rowAtX(double x, const QRectF& plot) const
{
    if (rows_.empty() || residueMax_ <= residueMin_)
        return std::nullopt;
    const double span = static_cast<double>(residueMax_ - residueMin_);
    const double slotWidth = plot.width() / std::max(1.0, span + 1.0);
    // Hit-test against the actual bar centres (residue + kind offset).
    // Use squared distance in (residue, sub-slot) space.
    std::optional<std::size_t> best;
    double bestDist = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < rows_.size(); ++i) {
        const double xCenter = plot.left()
            + (rows_[i].residue_index - residueMin_) / std::max(1.0, span) * plot.width()
            + kindOffsetFraction(rows_[i].kind) * slotWidth;
        const double d = std::abs(x - xCenter);
        if (d < bestDist) {
            bestDist = d;
            best = i;
        }
    }
    return best;
}

QString SequenceBarPanel::tooltipLine(int /*frame*/) const {
    // The hover-info we want is per-bar, not per-frame; surfaced via the
    // tooltip the widget shows over the panel. Returning empty so the
    // widget's frame-keyed temporal tooltip composer doesn't include
    // this panel.
    return {};
}

std::optional<model::SignalBinding>
SequenceBarPanel::mousePressInPlot(QMouseEvent* event, const PanelGeometry& geometry)
{
    if (event->button() != Qt::LeftButton)
        return std::nullopt;
    const auto row = rowAtX(event->pos().x(), geometry.plot);
    if (!row || !bindingForRow_)
        return std::nullopt;
    return bindingForRow_(*row);
}

void SequenceBarPanel::paint(QPainter& p,
                             const PanelGeometry& geometry,
                             const PaintContext& context) const
{
    paintPanelBackground(p, geometry);

    // Determine y-range: caller's explicit bounds, else min/max from rows
    // with mild padding (0 baseline if values stay non-negative).
    ValueRange range;
    range.valid = !rows_.empty();
    if (yMin_ && yMax_) {
        range.min = *yMin_;
        range.max = *yMax_;
    } else if (range.valid) {
        range.min = range.max = rows_.front().value;
        for (const SequenceBarRow& row : rows_) {
            if (!std::isfinite(row.value)) continue;
            range.min = std::min(range.min, row.value);
            range.max = std::max(range.max, row.value);
        }
        if (range.min > 0.0) range.min = 0.0;
        range = padRange(range);
        if (yMin_) range.min = *yMin_;
        if (yMax_) range.max = *yMax_;
    }

    const QString title = unit_.isEmpty()
                              ? label_
                              : QStringLiteral("%1 (%2)").arg(label_, unit_);
    paintHeader(p,
                geometry,
                HeaderText{title, QStringLiteral("residues %1-%2")
                                       .arg(residueMin_ + 1).arg(residueMax_ + 1)},
                160,
                kStripTextMuted,
                /*hasBinding=*/false,
                /*revealHover=*/false);
    paintGrid(p, geometry.plot);
    paintYAxisLabels(p, geometry, range);

    if (!rows_.empty() && residueMax_ > residueMin_) {
        const QRectF& plot = geometry.plot;
        const double xSpan = static_cast<double>(residueMax_ - residueMin_);
        const double slotWidth = plot.width() / std::max(1.0, xSpan + 1.0);
        // Detect whether ANY row carries a non-zero kind discriminator.
        // When the panel has multi-kind rows (Reorient), narrow bars to
        // ~1/3 of the slot so the three sub-slot bars don't overlap;
        // single-kind panels (iRED, Dihedral) get a wider bar.
        bool hasMultiKind = false;
        for (const SequenceBarRow& r : rows_) {
            if (r.kind != 0) { hasMultiKind = true; break; }
        }
        const double barWidth = std::max(2.0,
            slotWidth * (hasMultiKind ? 0.22 : 0.75));

        p.setRenderHint(QPainter::Antialiasing, true);
        for (const SequenceBarRow& row : rows_) {
            if (!std::isfinite(row.value)) continue;
            const double xCenter = plot.left()
                + (row.residue_index - residueMin_) / std::max(1.0, xSpan) * plot.width()
                + kindOffsetFraction(row.kind) * slotWidth;
            const double yTop = yForValue(row.value, range, plot);
            const double yBase = (range.min <= 0.0 && range.max >= 0.0)
                                     ? yForValue(0.0, range, plot)
                                     : plot.bottom();
            const QRectF bar(xCenter - barWidth * 0.5,
                             std::min(yTop, yBase),
                             barWidth,
                             std::abs(yBase - yTop));
            // Colour-stripe by kind: distinguishes N-H vs Cα-Hα vs C=O
            // when one residue has multiple rows on the same panel.
            QColor fill = barColor_;
            switch (row.kind) {
            case 2: fill = fill.darker(120); break;
            case 3: fill = fill.lighter(120); break;
            default: break;
            }
            p.fillRect(bar, fill);
        }
    }

    // Residue ticks along the bottom — fewer than dense per-bar labels
    // so the axis stays legible at ~50 residues.
    if (residueMax_ > residueMin_) {
        p.setFont(monoFont(10));
        p.setPen(kStripTextMuted);
        const QRectF& plot = geometry.plot;
        for (int g = 0; g <= 4; ++g) {
            const double frac = g / 4.0;
            const double x = plot.left() + plot.width() * frac;
            const int tickResidue = residueMin_
                + static_cast<int>(std::round(frac * (residueMax_ - residueMin_)));
            const QRectF labelRect(x - 24.0, plot.bottom() + 4.0, 48.0, 16.0);
            p.drawText(labelRect,
                       Qt::AlignHCenter | Qt::AlignVCenter,
                       QStringLiteral("r%1").arg(tickResidue + 1));
        }
    }

    paintPlotBorder(p, geometry.plot);

    (void)context;  // sequence-bar is static; no per-frame redraw
}

}  // namespace h5reader::app

// AbstractStripPanel — contract for one stackable panel inside
// StripStackWidget. Promoted out of StripStackWidget.cpp's anonymous
// namespace so new panel subclasses (SequenceBarPanel,
// PowerSpectrumPanel, LagDecayPanel, ChordCouplingPanel) can live in
// their own .cpp files and share the same painter/geometry/header
// helpers used by TemporalStripPanel + SpectrumStripPanel.
//
// Subclasses implement paint() + hasRevealBinding()/revealBinding() and
// may optionally override:
//   - preferredAspect()       — let the stack letterbox to a fixed
//                               aspect (chord wants square).
//   - mousePressInPlot(...)   — new panels handle their own click
//   - mouseMoveInPlot(...)      gestures. Temporal panels do NOT use
//                               these — their drag-select-range gesture
//                               lives in StripStackWidget directly.
//
// The free helpers (panelRectForIndex, plotRectForPanel,
// drawRevealButton, paintHeader, paintGrid, ...) are the same helpers
// the existing strip panels use; making them public-in-header lets new
// panel subclasses keep the same look without duplicating paint code.

#pragma once

#include "../model/SignalDictionary.h"

#include <QColor>
#include <QFont>
#include <QPoint>
#include <QRectF>
#include <QSize>
#include <QString>

#include <optional>

class QMouseEvent;
class QPainter;

namespace h5reader::app {

class TimeViewportController;

// --- theme + fonts -------------------------------------------------------

// Color constants previously local to StripStackWidget.cpp. inline const
// (not constexpr — QColor's ctor isn't constexpr) gives one definition
// across the program.
inline const QColor kStripCanvas(17, 18, 23);
inline const QColor kStripPanel(24, 27, 31);
inline const QColor kStripPanelBorder(43, 48, 56);
inline const QColor kStripGrid(48, 54, 66);
inline const QColor kStripText(216, 217, 218);
inline const QColor kStripTextMuted(154, 160, 166);
inline const QColor kStripCursor(217, 26, 26);
inline const QColor kStripSelection(87, 148, 242, 55);
inline const QColor kStripHover(199, 208, 217);
inline const QColor kStripReveal(115, 229, 214);

QFont uiFont(int px, QFont::Weight weight = QFont::Normal);
QFont monoFont(int px);
QString fmtValue(double value, int decimals = 2);

// --- geometry + paint context ------------------------------------------

constexpr int kStripMinPanelHeight = 72;
constexpr int kStripPanelGap = 6;

struct StackGeometry {
    QSize viewportSize;
    int panelCount = 0;
};

struct PanelGeometry {
    QRectF panel;
    QRectF plot;
    QRectF reveal;
};

struct HeaderText {
    QString title;
    QString readout;
};

struct ValueRange {
    double min = 0.0;
    double max = 1.0;
    bool valid = false;
};

struct FrameWindow {
    int first = 0;
    int last = 0;
    bool valid() const { return first <= last; }
};

struct TimeScale {
    int first = 0;
    int last = 0;
    double xForFrame(int frame, const QRectF& plot) const;
    int frameAt(const QPoint& pos, const QRectF& plot) const;
};

struct PaintContext {
    const TimeViewportController* viewport = nullptr;
    int currentFrame = 0;
    bool hasHover = false;
    QPoint hoverPos;
    TimeScale time;
};

// --- helpers ------------------------------------------------------------

QRectF panelRectForIndex(const StackGeometry& stack, int panelIndex);
QRectF plotRectForPanel(const QRectF& panel);
QRectF revealRectForPanel(const QRectF& panel);

// Letterbox `panel` to `aspect` (width/height) centred. nullopt → full rect.
QRectF letterboxedPanelRect(const QRectF& panel, std::optional<double> aspect);

PanelGeometry panelGeometryForIndex(const StackGeometry& stack,
                                    int panelIndex,
                                    std::optional<double> preferredAspect = std::nullopt);

ValueRange padRange(ValueRange range);
double yForValue(double value, const ValueRange& range, const QRectF& plot);

void drawRevealButton(QPainter& p, const QRectF& r, bool hover);
void paintPanelBackground(QPainter& p, const PanelGeometry& geometry);
void paintHeader(QPainter& p,
                 const PanelGeometry& geometry,
                 const HeaderText& text,
                 int readoutWidth,
                 const QColor& readoutColor,
                 bool hasBinding,
                 bool revealHover);
void paintGrid(QPainter& p, const QRectF& plot);
void paintYAxisLabels(QPainter& p, const PanelGeometry& geometry, const ValueRange& range);
void paintPlotBorder(QPainter& p, const QRectF& plot);
void paintTimeTicks(QPainter& p, const QRectF& plot, const TimeScale& time);
void paintSpectrumTicks(QPainter& p, const QRectF& plot, double xMin, double xMax);
void paintSelectedTimeRange(QPainter& p, const PanelGeometry& geometry, const PaintContext& context);
void paintTemporalCursor(QPainter& p, const PanelGeometry& geometry, const PaintContext& context);

// --- abstract base ------------------------------------------------------

class AbstractStripPanel {
public:
    virtual ~AbstractStripPanel() = default;

    virtual void paint(QPainter& p,
                       const PanelGeometry& geometry,
                       const PaintContext& context) const = 0;

    virtual bool hasRevealBinding() const = 0;
    virtual model::SignalBinding revealBinding() const = 0;

    virtual QString tooltipLine(int /*frame*/) const { return {}; }

    // Letterbox aspect (width/height) within the assigned rect when set.
    // Default nullopt = use the full rect.
    virtual std::optional<double> preferredAspect() const { return std::nullopt; }

    // Whether the cursor over this panel hits its plot area. The
    // existing temporal-panel drag-select-range gesture in
    // StripStackWidget keys off this for "is the press / drag inside a
    // temporal plot." Defaults false (chord, sequence-bar, etc. opt in).
    virtual bool plotContains(const PanelGeometry&, const QPoint&) const {
        return false;
    }

    // Optional in-plot gesture hooks for new panel types. Default no-op
    // — temporal panels don't use these (their drag-select-range is in
    // StripStackWidget). Chord / sequence-bar opt in here for
    // click-arc-to-highlight, click-bar-to-focus, etc. Returning a
    // SignalBinding from mousePressInPlot tells the widget to emit
    // revealRequested with that binding (scene focuses on the bar's
    // atom/residue); returning nullopt means the panel handled the
    // click itself or ignored it.
    virtual std::optional<model::SignalBinding>
        mousePressInPlot(QMouseEvent*, const PanelGeometry&) { return std::nullopt; }
    virtual void mouseMoveInPlot(QMouseEvent*, const PanelGeometry&) {}

    bool revealContains(const PanelGeometry& geometry, const QPoint& pos) const {
        return hasRevealBinding() && geometry.reveal.contains(pos);
    }
};

}  // namespace h5reader::app

// QtAtomTimeSeriesDock — scalar-vs-frame line chart for the picked atom.
//
// Uses Qt6 Charts (QChart + QLineSeries + QValueAxis). A QComboBox at
// the top selects which per-atom per-frame scalar to plot; the chart
// redraws for every atom pick or scalar change, and a vertical
// cursor line follows the current frame during playback / scrub.
//
// The dock is tabified with QtAtomInspectorDock in the right dock
// area — same picked-atom context, two views.

#pragma once

#include "../model/QtConformation.h"
#include "../model/QtProtein.h"

#include <QDockWidget>
#include <QPointer>

#include <cstddef>

QT_BEGIN_NAMESPACE
namespace QtCharts {}
QT_END_NAMESPACE

class QComboBox;
class QLabel;

// Qt 6.2+ merges QtCharts namespace into QT_CHARTS_USE_NAMESPACE.
// Using the global namespace here works under both Qt 6.4 (our
// build target) and later.
class QChart;
class QChartView;
class QLineSeries;
class QValueAxis;

namespace h5reader::app {

class QtAtomTimeSeriesDock final : public QDockWidget {
    Q_OBJECT
public:
    explicit QtAtomTimeSeriesDock(QWidget* parent = nullptr);
    ~QtAtomTimeSeriesDock() override = default;

    // Bind the typed model once at load.
    void setContext(const model::QtProtein*      protein,
                    const model::QtConformation* conformation);

public slots:
    void setPickedAtom(std::size_t atomIdx);
    void setFrame(int t);
    void clearSelection();

private slots:
    void onScalarChanged(int comboIndex);

private:
    // Walk the conformation's frames, extract the currently-selected
    // scalar for atomIdx_ at each frame, populate the line series,
    // rescale the Y axis, and reposition the cursor.
    void rebuildSeries();

    // Reposition the vertical cursor line at the current frame; cheap.
    void updateCursor();

    QPointer<QComboBox>  scalarCombo_;
    QPointer<QChartView> chartView_;
    QPointer<QLabel>     valueLabel_;

    // Chart + axes + data/cursor series. Owned by QChart; raw pointers.
    QChart*      chart_       = nullptr;
    QLineSeries* dataSeries_  = nullptr;
    QLineSeries* cursorSeries_ = nullptr;
    QValueAxis*  xAxis_       = nullptr;
    QValueAxis*  yAxis_       = nullptr;

    const model::QtProtein*      protein_      = nullptr;
    const model::QtConformation* conformation_ = nullptr;
    bool                         hasSelection_ = false;
    std::size_t                  atomIdx_      = 0;
    int                          frame_        = 0;
};

}  // namespace h5reader::app

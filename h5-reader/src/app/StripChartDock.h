// StripChartDock — the trajectory strip-chart instrument.
//
// The over-time companion to the 3-D measurement overlay. Where the overlay
// holds the geometry in space, this holds the OBSERVABLE over the trajectory
// ("first, Rosalind Franklin had to go get some graph paper" — memory
// project_h5reader_killer_app_multiatom_compare_20260526).
//
// Two stacked panels:
//   * TIME domain — the derived geometry (distance / angle / dihedral over
//     frames, via model::Measure) with a scrolling fixed-width window (the cure
//     for the old dock's "scrunch as you go"), a Fit-all toggle, major+minor
//     science grids, a red dashed playhead, and a digital value readout.
//   * FREQUENCY domain — the FFT power spectrum of that series
//     (model::ComputePowerSpectrum) so a periodic motion shows up as a peak;
//     the dominant period is read out. The dihedral is phase-unwrapped before
//     the transform so a ±180°-straddling oscillation reads as one clean peak.
//     Recomputed on selection change (the whole-trajectory spectrum is
//     frame-independent), so it carries no playhead.
//
// Seeing periodicity in a dihedral is the headline — directly observable and
// assumption-free. A future analysis tab can host the worthier 3-D methods
// (PCA / dihedral-PCA collective modes; recurrence / RQA) when wanted; see
// notes/PLANNED_ANALYSIS_METHODS.md.
//
// Engine: Qt Charts (QChartView), deliberately NOT a second VTK surface
// (decision 2026-05-27). The app's heavy GPU load is the molecule scene's
// QVTKOpenGLNativeWidget; a 2-D dockable chart neither needs VTK's muscle nor a
// second OpenGL context (which would demand Qt::AA_ShareOpenGLContexts and
// dockable-context-recreation handling). The FFT math is Eigen-only (no VTK).
//
// Trajectory-only, like the time-series dock: a single pose has no frame axis.

#pragma once

#include <QDockWidget>
#include <QPointer>

class QCheckBox;
class QLabel;

// Qt 6.2+ merges the QtCharts namespace into the global namespace, so these
// forward declarations work as-is (same as QtAtomTimeSeriesDock).
class QChart;
class QChartView;
class QLineSeries;
class QValueAxis;

namespace h5reader::model {
class QtProtein;
class Conformation;
class AtomSelection;
}

namespace h5reader::app {

class StripChartDock final : public QDockWidget {
    Q_OBJECT
public:
    explicit StripChartDock(QWidget* parent = nullptr);
    ~StripChartDock() override = default;

    // Bind the typed model once at load. conformation may be a single pose; the
    // dock then shows an explanatory placeholder rather than an empty axis.
    void setContext(const model::QtProtein* protein, model::Conformation* conformation);

    // Bind the selection whose derived geometry this charts. The dock ACONNECTs
    // its changed()/cleared() itself (same pattern as SelectionDock).
    void setSelection(model::AtomSelection* selection);

public slots:
    // Move the playhead + slide the scrolling window + refresh the readout.
    // Cheap: the series + spectrum are fixed; per frame only the cursor moves.
    void setFrame(int t);

    // Selection membership changed — rebuild the geometry series AND its spectrum.
    void onSelectionChanged();

    // Selection emptied — clear both panels.
    void clearSelection();

private slots:
    void onFitToggled(bool on);

private:
    void rebuildGeometrySeries();  // time series + FFT, on selection change
    void updateCursorAndWindow();  // playhead + scrolling window, per frame

    QPointer<QChartView> chartView_;
    QPointer<QCheckBox>  fitAllBox_;
    QPointer<QLabel>     readout_;

    // FFT (frequency-domain) panel, stacked below the time-domain chart.
    QPointer<QChartView> fftView_;
    QPointer<QLabel>     fftReadout_;

    // Owned by chart_ / chartView_ (Qt parent ownership); raw pointers, the
    // same convention as QtAtomTimeSeriesDock.
    QChart*      chart_        = nullptr;
    QLineSeries* geomSeries_   = nullptr;
    QLineSeries* cursorSeries_ = nullptr;
    QValueAxis*  xAxis_        = nullptr;
    QValueAxis*  yAxis_        = nullptr;

    QChart*      fftChart_  = nullptr;
    QLineSeries* fftSeries_ = nullptr;
    QValueAxis*  fftXAxis_  = nullptr;
    QValueAxis*  fftYAxis_  = nullptr;

    const model::QtProtein*        protein_ = nullptr;
    QPointer<model::Conformation>  conformation_;
    QPointer<model::AtomSelection> selection_;

    int    frame_      = 0;
    bool   fitAll_     = false;
    bool   haveSeries_ = false;
    double yMin_       = 0.0;
    double yMax_       = 1.0;
};

}  // namespace h5reader::app

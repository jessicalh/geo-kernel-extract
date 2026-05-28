// StripChartDock — the trajectory strip-chart instrument.
//
// The over-time companion to the 3-D measurement overlay. Where the overlay
// holds the geometry in space, this holds the OBSERVABLE over the trajectory
// ("first, Rosalind Franklin had to go get some graph paper" — memory
// project_h5reader_killer_app_multiatom_compare_20260526).
//
// A REAL strip chart (rearchitected 2026-05-27). The earlier version built the
// WHOLE per-frame series up front and slid a cursor over a finished curve — that
// cannot scale to the microsecond (~1e6-frame) trajectories ahead, and it is not
// a strip chart. Now each plotted series is an authoritative value buffer WE own
// (model::ChannelBuffer — a std::vector grown one append per playback frame);
// the PAST stays, the FUTURE is never drawn, and there is NO decimation of the
// stored data (a 1e6-frame buffer is ≈ 8 MB). See StripChartChannel.h for the
// buffer/source seam and why those are plain data, not QObjects.
//
// Panels (each a "track" = an owned buffer + its source + Qt Charts handles):
//   * STRUCTURAL — the selection's geometry (distance / angle / dihedral via
//     model::Measure), the always-present track.
//   * SHIELDING  — (added with the DFT channel) the focus atom's DFT σ; a gap
//     where no DFT frame exists, never a faked zero.
//   * FREQUENCY  — the FFT power spectrum of the structural track's collected
//     values (model::ComputePowerSpectrum); a periodic motion shows as a peak.
//
// Engine: Qt Charts fed only the VISIBLE WINDOW sliced from the authoritative
// buffer — deliberately NOT a second VTK surface (the molecule scene's
// QVTKOpenGLNativeWidget is the heavy GPU surface; a second context would demand
// Qt::AA_ShareOpenGLContexts), and deliberately NOT a software QLineSeries of
// all 1e6 points. The buffer + source are the durable seam; the Qt Charts
// presentation is the swap-out half when a richer widget lands (memory: keep the
// data, swap the view). Multi-frame-only: a single pose has no frame axis.

#pragma once

#include "../model/DftSigmaStrips.h"
#include "../model/SignalDictionary.h"
#include "../model/StripChartChannel.h"  // ChannelBuffer, ChannelSource

#include <QDockWidget>
#include <QPointF>
#include <QPointer>
#include <QString>
#include <QVector>

#include <memory>
#include <optional>

class QCheckBox;
class QColor;
class QLabel;
class QSpinBox;
class QVBoxLayout;

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
class DftShieldingStore;
}

namespace h5reader::app {

class TimeViewportController;
class StripStackWidget;

class StripChartDock final : public QDockWidget {
    Q_OBJECT
public:
    explicit StripChartDock(QWidget* parent = nullptr);
    ~StripChartDock() override = default;

    // Bind the typed model once at load. conformation may be a single pose; the
    // dock then shows an explanatory placeholder rather than an empty axis.
    void setContext(const model::QtProtein* protein, model::Conformation* conformation);

    // Bind the selection whose derived geometry this charts. Safe to call more
    // than once: any prior selection's signals are disconnected first (ACONNECT
    // cannot pass Qt::UniqueConnection, so the re-bind guard is explicit).
    void setSelection(model::AtomSelection* selection);

    // Bind the DFT shielding provider (present only when the run has a dft/
    // campaign). The shielding panel appears only once a store is set; its
    // channel follows the selection's FOCUS atom (σ of one nucleus over the
    // trajectory). nullptr leaves the panel hidden (no-DFT runs).
    void setDftStore(model::DftShieldingStore* store);

    // Bind shared time viewport state. Playback still appends values through
    // setFrame(); this controller owns the visible range/follow mode so future
    // panels can stay synchronized without talking to this dock directly.
    void setTimeViewport(TimeViewportController* viewport);

signals:
    void revealRequested(const model::SignalBinding& binding);

public slots:
    // Playback advanced / scrubbed. Forward motion APPENDS the new frames'
    // values to each track's buffer; backward motion only moves the cursor (the
    // past is kept). Then each track renders its visible window.
    void setFrame(int t);

    // Selection membership changed — rebuild the channel set: reset the buffers,
    // rebind the sources, backfill 0..currentFrame so the new series shows its
    // history up to now (never the future).
    void onSelectionChanged();

    // Selection emptied — clear all tracks.
    void clearSelection();

private slots:
    void onFitToggled(bool on);
    void onResidueModeToggled(bool on);
    void onVisibleRangeChanged(int first, int last);

private:
    // One plotted series: its authoritative buffer (the durable half), the
    // source that fills it, and the Qt Charts handles (the throwaway half).
    // active == this track currently has a channel bound.
    struct Track {
        model::ChannelBuffer buffer;
        model::ChannelSource source;
        bool                 active = false;
        std::optional<model::SignalBinding> binding;

        // Owned by chart (Qt parent ownership); raw pointers, the QtCharts
        // convention used across the reader's docks.
        QChart*      chart   = nullptr;
        QChartView*  view    = nullptr;
        QLineSeries* series  = nullptr;  // the visible-window slice only
        QLineSeries* cursor  = nullptr;  // red dashed playhead
        QValueAxis*  xAxis   = nullptr;
        QValueAxis*  yAxis   = nullptr;
        QLabel*      readout = nullptr;  // digital value at the current frame
    };

    // Build a track's panel (chart + data/cursor series + axes + readout) into
    // `into`. Called once per panel in the ctor; the channel is bound later.
    void buildTrack(Track& tr, QVBoxLayout* into, const QColor& traceColor, int stretch);

    void rebuildChannels();             // selection changed: rebind sources + backfill
    void extendTo(int t);               // append source(f) for f in (lastFrame, t]
    bool extendDenseTrackTo(Track& tr, int t);
    void renderTrack(Track& tr, int t); // visible-window slice + cursor + readout
    void refreshStackTracks();          // expose active buffers to the custom stack widget
    void rebuildFft();                  // FFT over the structural track's valid prefix

    // Bind/rebuild the DFT shielding channel to the current focus atom (called
    // from rebuildChannels when a store is present).
    void bindDftChannel();
    void clearDftSignalStrips();
    void bindResidueDihedralChannels();
    void clearTrack(Track& tr, const QString& readout = QStringLiteral("—"));
    void clearResidueDihedralChannels();
    QString atomDisplayLabel(std::size_t atomIdx) const;
    QString residueDisplayLabel(std::size_t residueIdx) const;
    QString selectionTupleLabel(const std::vector<std::size_t>& atoms) const;
    void updateViewportReadout(int first, int last);
    std::optional<model::SignalBinding> currentDftBinding() const;
    std::optional<model::SignalBinding> sourceBindingForFft() const;

    // ----- Tracks -----
    Track structural_;  // panel 1: the selection's geometry (always present)
    Track dft_;         // panel 2: DFT σ of the focus atom (shown only with a store)
    std::vector<Track> residueDihedrals_;
    std::unique_ptr<model::DftSigmaAtomTimeStrip> dftTimeStrip_;
    std::unique_ptr<model::DftSigmaAtomFftStrip> dftFftStrip_;

    // ----- Frequency-domain panel (Qt Charts; small N, bound to structural_) --
    QChart*      fftChart_   = nullptr;
    QChartView*  fftView_    = nullptr;
    QLineSeries* fftSeries_  = nullptr;
    QValueAxis*  fftXAxis_   = nullptr;
    QValueAxis*  fftYAxis_   = nullptr;
    QPointer<QLabel> fftReadout_;

    QPointer<QCheckBox> fitAllBox_;
    QPointer<QCheckBox> followBox_;
    QPointer<QCheckBox> residueModeBox_;
    QPointer<QSpinBox> windowFramesSpin_;
    QPointer<QLabel> viewportReadout_;
    QPointer<StripStackWidget> stackWidget_;

    const model::QtProtein*            protein_ = nullptr;
    QPointer<model::Conformation>      conformation_;
    QPointer<model::AtomSelection>     selection_;
    QPointer<model::DftShieldingStore> dftStore_;
    QPointer<TimeViewportController>   timeViewport_;

    int  frame_       = 0;
    bool fitAll_      = false;
    bool residueMode_ = false;
    int  fftRecomputeAtSize_ = 0;  // throttle: next valid-count at which to redo the FFT
    QVector<QPointF> fftPoints_;
    QString fftReadoutText_ = QStringLiteral("—");
    QVector<QPointF> dftFftPoints_;
    QString dftFftReadoutText_ = QStringLiteral("—");
};

}  // namespace h5reader::app

// StripStackWidget -- custom QPainter surface for stacked trajectory strips.
//
// Data ownership stays outside the widget. Each track points at a
// model::ChannelBuffer owned by DashboardDisplayController and handed through
// DashboardStripDock. This widget is only the renderer and gesture surface:
// visible range, selected range, and playback frame requests all go through
// TimeViewportController.

#pragma once

#include <QColor>
#include <QPointF>
#include <QPointer>
#include <QString>
#include <QVector>
#include <QWidget>

#include "../model/SignalDictionary.h"
#include "AbstractStripPanel.h"

#include <memory>
#include <vector>

namespace h5reader::model {
struct ChannelBuffer;
}

namespace h5reader::app {

class TimeViewportController;

class StripStackWidget final : public QWidget {
    Q_OBJECT

public:
    struct Track {
        const model::ChannelBuffer* buffer = nullptr;
        QColor color;
        bool hasBinding = false;
        model::SignalBinding binding;
    };

    struct SpectrumTrack {
        const QVector<QPointF>* points = nullptr;
        QString label;
        QString xUnit;
        QString yUnit;
        QString readout;
        QColor color;
        bool hasBinding = false;
        model::SignalBinding binding;
    };

    explicit StripStackWidget(QWidget* parent = nullptr);
    ~StripStackWidget() override = default;

    void setTracks(QVector<Track> tracks);
    void setSpectrumTracks(QVector<SpectrumTrack> tracks);

    // setOwnedPanels — canonical entry for heterogeneous panel types
    // (SequenceBarPanel, PowerSpectrumPanel, LagDecayPanel,
    // ChordCouplingPanel, ...) that don't ride the Track / SpectrumTrack
    // data shapes. These panels own their own data and render after the
    // temporal + spectrum tracks in the stack. New panel subclasses opt
    // in here rather than overloading Track.
    void setOwnedPanels(std::vector<std::unique_ptr<AbstractStripPanel>> panels);

    void setTimeViewport(TimeViewportController* viewport);
    void setCurrentFrame(int frame);
    // Backwards-compatible per-section counts after L-2b unified the
    // storage into a single panels_ vector. trackCount() returns the
    // number of temporal panels (front of panels_), spectrumTrackCount()
    // the middle, ownedPanelCount() the trailing.
    int trackCount() const { return static_cast<int>(n_temporal_); }
    int spectrumTrackCount() const { return static_cast<int>(n_spectrum_); }
    int ownedPanelCount() const {
        return static_cast<int>(panels_.size() - n_temporal_ - n_spectrum_);
    }

signals:
    void revealRequested(const model::SignalBinding& binding);

protected:
    void paintEvent(QPaintEvent* event) override;
    void mousePressEvent(QMouseEvent* event) override;
    void mouseMoveEvent(QMouseEvent* event) override;
    void mouseReleaseEvent(QMouseEvent* event) override;
    void mouseDoubleClickEvent(QMouseEvent* event) override;
    void leaveEvent(QEvent* event) override;
    void wheelEvent(QWheelEvent* event) override;

private:
    QRectF trackRect(int index) const;
    QRectF plotRect(const QRectF& trackRect) const;
    QRectF revealRect(const QRectF& trackRect) const;
    int panelCount() const;
    void updateMinimumHeight();
    bool timePlotContains(const QPoint& pos) const;
    bool revealAt(const QPoint& pos, model::SignalBinding* binding) const;
    int frameAt(const QPoint& pos) const;
    QString tooltipText(int frame) const;

    // L-2b (2026-05-29): one ordered panels_ vector replaces the prior
    // tracks_ / spectrumTracks_ / ownedPanels_ three-bucket storage.
    // Layout: indices [0..n_temporal_) hold TemporalStripPanel instances,
    // [n_temporal_..n_temporal_+n_spectrum_) hold SpectrumStripPanel,
    // the trailing section holds caller-provided owned panels. Each
    // setter rebuilds only its section, preserving the others' order.
    // TemporalStripPanel + SpectrumStripPanel now own their Track /
    // SpectrumTrack by value (was const-ref) so the panel survives a
    // setTracks() reallocation.
    std::vector<std::unique_ptr<AbstractStripPanel>> panels_;
    std::size_t n_temporal_ = 0;
    std::size_t n_spectrum_ = 0;
    QPointer<TimeViewportController> viewport_;
    int currentFrame_ = 0;
    bool hasHover_ = false;
    QPoint hoverPos_;
    bool selecting_ = false;
    bool dragSelecting_ = false;
    bool pressWasFollowing_ = false;
    int selectionAnchor_ = 0;
    QPoint pressPos_;
};

}  // namespace h5reader::app

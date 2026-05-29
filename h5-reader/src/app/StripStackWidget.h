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
    void setTimeViewport(TimeViewportController* viewport);
    void setCurrentFrame(int frame);
    int trackCount() const { return tracks_.size(); }
    int spectrumTrackCount() const { return spectrumTracks_.size(); }

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

    QVector<Track> tracks_;
    QVector<SpectrumTrack> spectrumTracks_;
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

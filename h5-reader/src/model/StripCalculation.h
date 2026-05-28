// StripCalculation -- common Qt contract for dashboard strip calculations.
//
// Concrete calculations own their binding, frame sampling, durable strip buffer,
// and render preparation. The dashboard owns StripCalculation instances, not
// transient chart widgets; source readers only provide current-frame data.

#pragma once

#include "SignalDictionary.h"

#include <QObject>
#include <QPointF>
#include <QString>
#include <QVector>

#include <cstddef>

namespace h5reader::model {

struct StripSpec {
    QString id;
    QString label;
    SignalKey signal;
    SignalAnchorKind anchorKind = SignalAnchorKind::None;
    SignalDomainKind domainKind = SignalDomainKind::FrameTime;
};

struct StripRenderData {
    QString label;
    QString xUnit;
    QString yUnit;
    QString readout;
    SignalDomainKind domainKind = SignalDomainKind::FrameTime;
    QVector<QPointF> points;
};

class StripCalculation : public QObject {
    Q_OBJECT

public:
    explicit StripCalculation(QObject* parent = nullptr) : QObject(parent) {}
    ~StripCalculation() override = default;

    virtual StripSpec spec() const = 0;
    virtual bool canBind(const SignalBinding& binding) const = 0;
    virtual void bind(const SignalBinding& binding) = 0;

    // Frame-time strips append only up to frame. Frequency-domain strips may
    // recompute their points from the same underlying signal over the collected
    // prefix/window, but they still do not read future trajectory data.
    virtual void extendToFrame(std::size_t frame) = 0;
    virtual StripRenderData renderData() const = 0;

signals:
    void bindingChanged();
    void bufferChanged(std::size_t frame);
};

}  // namespace h5reader::model

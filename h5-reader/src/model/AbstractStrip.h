// AbstractStrip -- common contract for dashboard strip instances.
//
// Concrete strips own their sampling/render preparation. They bind to a
// SignalKey/anchor and expose render data in a shared shape so the dashboard
// can stack many different strip types without each one inventing UI rules.

#pragma once

#include "SignalDictionary.h"

#include <QPointF>
#include <QString>
#include <QVector>

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

class AbstractStrip {
public:
    virtual ~AbstractStrip() = default;

    virtual StripSpec spec() const = 0;
    virtual bool canBind(const SignalBinding& binding) const = 0;
    virtual void bind(const SignalBinding& binding) = 0;

    // Frame-time strips append only up to frame. Frequency-domain strips may
    // recompute their points from the same underlying signal over the collected
    // prefix/window, but they still do not read future trajectory data.
    virtual void extendToFrame(std::size_t frame) = 0;
    virtual StripRenderData renderData() const = 0;
};

}  // namespace h5reader::model

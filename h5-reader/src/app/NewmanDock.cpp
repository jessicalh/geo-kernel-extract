#include "NewmanDock.h"

#include "../diagnostics/ConnectionAuditor.h"
#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"

#include "../model/Conformation.h"
#include "../model/QtProtein.h"

#include <QFontMetricsF>
#include <QHBoxLayout>
#include <QLabel>
#include <QPainter>
#include <QPaintEvent>
#include <QPen>
#include <QToolButton>
#include <QVBoxLayout>
#include <QWidget>

#include <algorithm>
#include <cmath>

namespace h5reader::app {

namespace {

QString kindName(NewmanKind k) {
    return k == NewmanKind::Phi ? QStringLiteral("phi") : QStringLiteral("psi");
}

void drawCenteredText(QPainter& p, const QPointF& at, const QString& text) {
    const QFontMetricsF fm(p.font());
    const double w = fm.horizontalAdvance(text);
    const double h = fm.height();
    p.drawText(QRectF(at.x() - w / 2.0, at.y() - h / 2.0, w, h), Qt::AlignCenter, text);
}

}  // namespace

// ---------------------------------------------------------------------------
// NewmanView
// ---------------------------------------------------------------------------

NewmanView::NewmanView(QWidget* parent) : QWidget(parent) {
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("NewmanView"));
    setMinimumSize(minimumSizeHint());
}

void NewmanView::setData(const NewmanProjection& drawn, const NewmanProjection& other) {
    drawn_ = drawn;
    other_ = other;
    has_   = true;
    update();
}

void NewmanView::setEmpty(const QString& hint) {
    has_  = false;
    hint_ = hint;
    update();
}

void NewmanView::paintEvent(QPaintEvent*) {
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing, true);
    const QRectF area = rect();
    p.fillRect(area, QColor(252, 252, 252));

    if (!has_ || !drawn_.valid) {
        p.setPen(QColor(120, 120, 120));
        const QString msg = !has_
            ? hint_
            : QStringLiteral("%1 not defined here:\n%2").arg(kindName(drawn_.kind), drawn_.invalidReason);
        p.drawText(area.adjusted(10, 10, -10, -10), Qt::AlignCenter | Qt::TextWordWrap, msg);
        return;
    }

    const QPointF c(area.center().x(), area.center().y() - 6.0);
    const double  R = std::min(area.width(), area.height()) * 0.30;
    auto pt = [&](double radius, double ang) {
        return QPointF(c.x() + radius * std::cos(ang), c.y() - radius * std::sin(ang));
    };

    const QColor frontCol(35, 90, 175), backCol(170, 55, 35);
    const QColor refFront(20, 140, 70), refBack(205, 120, 0);

    // Back atom = circle; front atom = centre hub.
    p.setPen(QPen(QColor(110, 110, 110), 1.4));
    p.setBrush(Qt::NoBrush);
    p.drawEllipse(c, R, R);
    p.setBrush(QColor(70, 70, 70));
    p.setPen(Qt::NoPen);
    p.drawEllipse(c, 3.2, 3.2);

    // Back spokes: from the rim outward.
    for (const NewmanSpoke& s : drawn_.spokes) {
        if (s.front) continue;
        const QColor col = s.reference ? refBack : backCol;
        p.setPen(QPen(col, s.reference ? 2.8 : 1.7));
        p.drawLine(pt(R, s.angleRad), pt(R * 1.34, s.angleRad));
        p.setPen(col);
        drawCenteredText(p, pt(R * 1.52, s.angleRad), s.label);
    }
    // Front spokes: from the hub to just inside the rim.
    for (const NewmanSpoke& s : drawn_.spokes) {
        if (!s.front) continue;
        const QColor col = s.reference ? refFront : frontCol;
        p.setPen(QPen(col, s.reference ? 2.8 : 1.7));
        p.drawLine(c, pt(R * 0.96, s.angleRad));
        p.setPen(col);
        drawCenteredText(p, pt(R * 0.60, s.angleRad), s.label);
    }

    // Front/back atom identity + the torsion readout.
    p.setPen(QColor(40, 40, 40));
    QFont fhdr = p.font();
    fhdr.setPointSizeF(fhdr.pointSizeF() * 0.95);
    p.setFont(fhdr);
    p.drawText(area.adjusted(6, 4, -6, 0), Qt::AlignLeft | Qt::AlignTop,
               QStringLiteral("front %1  -  back %2").arg(drawn_.frontLabel, drawn_.backLabel));

    QFont fbig = p.font();
    fbig.setPointSizeF(fbig.pointSizeF() * 1.6);
    fbig.setBold(true);
    p.setFont(fbig);
    p.setPen(Qt::black);
    p.drawText(QRectF(area.left(), area.bottom() - 30, area.width(), 26),
               Qt::AlignHCenter | Qt::AlignVCenter,
               QStringLiteral("%1 = %2 deg").arg(kindName(drawn_.kind))
                   .arg(drawn_.torsionDeg, 0, 'f', 1));
}

// ---------------------------------------------------------------------------
// NewmanDock
// ---------------------------------------------------------------------------

NewmanDock::NewmanDock(QWidget* parent)
    : QDockWidget(QStringLiteral("Newman"), parent) {
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("NewmanDock"));

    auto* container = new QWidget(this);
    auto* lay = new QVBoxLayout(container);
    lay->setContentsMargins(6, 6, 6, 6);
    lay->setSpacing(4);

    auto* top = new QHBoxLayout();
    caption_ = new QLabel(QStringLiteral("Select an atom"), container);
    caption_->setWordWrap(true);
    toggle_ = new QToolButton(container);
    toggle_->setText(QStringLiteral("phi / psi"));
    toggle_->setToolTip(QStringLiteral("Switch the drawn torsion between backbone phi and psi."));
    top->addWidget(caption_, 1);
    top->addWidget(toggle_, 0);
    lay->addLayout(top);

    view_ = new NewmanView(container);
    lay->addWidget(view_, 1);
    setWidget(container);

    ACONNECT(toggle_, &QToolButton::clicked, this, [this]() {
        kind_ = (kind_ == NewmanKind::Phi) ? NewmanKind::Psi : NewmanKind::Phi;
        recompute();
    });

    recompute();
}

void NewmanDock::setContext(const model::QtProtein* protein, const model::Conformation* conf) {
    ASSERT_THREAD(this);
    protein_ = protein;
    conf_    = conf;
    residue_ = -1;
    recompute();
}

void NewmanDock::setFocusAtom(std::size_t atomIdx) {
    ASSERT_THREAD(this);
    if (!protein_ || atomIdx >= protein_->atomCount()) {
        residue_ = -1;
    } else {
        const int r = protein_->atom(atomIdx).residueIndex;
        residue_ = (r >= 0 && static_cast<std::size_t>(r) < protein_->residueCount()) ? r : -1;
    }
    recompute();
}

void NewmanDock::setFrame(int frame) {
    ASSERT_THREAD(this);
    frame_ = frame;
    recompute();
}

void NewmanDock::clear() {
    ASSERT_THREAD(this);
    residue_ = -1;
    recompute();
}

void NewmanDock::recompute() {
    ASSERT_THREAD(this);
    if (!view_)
        return;
    if (!protein_ || !conf_ || residue_ < 0) {
        if (caption_) caption_->setText(QStringLiteral("Select an atom"));
        view_->setEmpty(QStringLiteral("Select an atom to see its backbone torsion."));
        return;
    }
    const std::size_t r = static_cast<std::size_t>(residue_);
    const std::size_t f = frame_ >= 0 ? static_cast<std::size_t>(frame_) : 0u;
    const NewmanProjection phi = ComputeNewmanProjection(*protein_, *conf_, r, f, NewmanKind::Phi);
    const NewmanProjection psi = ComputeNewmanProjection(*protein_, *conf_, r, f, NewmanKind::Psi);

    const NewmanProjection& drawn = (kind_ == NewmanKind::Phi) ? phi : psi;
    const NewmanProjection& other = (kind_ == NewmanKind::Phi) ? psi : phi;
    view_->setData(drawn, other);

    if (caption_) {
        auto fmt = [](const NewmanProjection& p) {
            return p.valid ? QStringLiteral("%1 %2").arg(kindName(p.kind))
                                 .arg(p.torsionDeg, 0, 'f', 1)
                           : QStringLiteral("%1 n/a").arg(kindName(p.kind));
        };
        caption_->setText(QStringLiteral("%1 %2   |   %3   %4")
                              .arg(phi.residueLabel)
                              .arg(residue_)
                              .arg(fmt(phi), fmt(psi)));
    }
}

}  // namespace h5reader::app

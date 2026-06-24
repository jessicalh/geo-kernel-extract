#include "NewmanPaint.h"

#include <QChar>
#include <QColor>
#include <QFont>
#include <QFontMetricsF>
#include <QPainter>
#include <QPainterPath>
#include <QPen>
#include <QPointF>
#include <QRectF>

#include <algorithm>
#include <cmath>

namespace h5reader::app {

namespace {

// Glyphs as code points so the source stays pure ASCII (no BOM / /utf-8 needed).
const QChar kPhi = QChar(char16_t(0x03C6));  // phi
const QChar kPsi = QChar(char16_t(0x03C8));  // psi
const QChar kDeg = QChar(char16_t(0x00B0));  // degree
const QChar kArrow = QChar(char16_t(0x2192));  // rightwards arrow (sight axis)

void drawCenteredText(QPainter& p, const QPointF& at, const QString& text) {
    const QFontMetricsF fm(p.font());
    const double w = fm.horizontalAdvance(text);
    const double h = fm.height();
    p.drawText(QRectF(at.x() - w / 2.0, at.y() - h / 2.0, w, h), Qt::AlignCenter, text);
}

}  // namespace

QString NewmanKindGlyph(NewmanKind kind) {
    return QString(kind == NewmanKind::Phi ? kPhi : kPsi);
}

void PaintNewman(QPainter& p, const QRectF& area, const NewmanProjection& drawn,
                 const NewmanOccupancy& occupancy, bool darkBg) {
    p.setRenderHint(QPainter::Antialiasing, true);

    const QColor ink        = darkBg ? QColor(225, 225, 225) : QColor(40, 40, 40);
    const QColor mutedInk    = darkBg ? QColor(155, 155, 155) : QColor(120, 120, 120);
    const QColor discPen     = darkBg ? QColor(150, 150, 150) : QColor(110, 110, 110);
    const QColor hubCol      = darkBg ? QColor(205, 205, 205) : QColor(70, 70, 70);
    const QColor arcPen      = darkBg ? QColor(170, 170, 170) : QColor(150, 150, 150);
    const QColor readoutInk  = darkBg ? QColor(245, 245, 245) : QColor(0, 0, 0);

    if (!drawn.valid) {
        p.setPen(mutedInk);
        p.drawText(area.adjusted(10, 10, -10, -10), Qt::AlignCenter | Qt::TextWordWrap,
                   QStringLiteral("%1 not defined here:\n%2")
                       .arg(NewmanKindGlyph(drawn.kind), drawn.invalidReason));
        return;
    }

    const QPointF c(area.center().x(), area.center().y() - 6.0);
    const double  R = std::min(area.width(), area.height()) * 0.32;
    auto pt = [&](double radius, double ang) {
        return QPointF(c.x() + radius * std::cos(ang), c.y() - radius * std::sin(ang));
    };

    // Occupancy haze (the d'Arsonval dial scale): a translucent pastel wedge per
    // angular bin, alpha by density, drawn FIRST so the disc + the live needle
    // sit on top. The torsion runs [-180, 180) == QPainter angle (0 at 3 o'clock,
    // CCW+), the same frame as the torsion arc + the back-reference spoke.
    if (!occupancy.bins.empty()) {
        float maxD = 0.0F;
        for (float d : occupancy.bins) maxD = std::max(maxD, d);
        if (maxD > 0.0F) {
            const int    n    = static_cast<int>(occupancy.bins.size());
            const double step = 360.0 / n;
            const QRectF disc(c.x() - R, c.y() - R, 2 * R, 2 * R);
            const QColor haze = darkBg ? QColor(150, 165, 245) : QColor(120, 140, 225);
            p.setPen(Qt::NoPen);
            for (int i = 0; i < n; ++i) {
                const double d = occupancy.bins[static_cast<std::size_t>(i)] / maxD;
                if (d <= 0.0) continue;
                QColor col = haze;
                col.setAlphaF(0.12 + 0.46 * d);
                p.setBrush(col);
                QPainterPath wedge;
                wedge.moveTo(c);
                wedge.arcTo(disc, -180.0 + i * step, step);
                wedge.closeSubpath();
                p.drawPath(wedge);
            }
        }
    }

    const QColor frontCol(35, 90, 175), backCol(170, 55, 35);
    const QColor refFront(20, 140, 70), refBack(205, 120, 0);

    // Back atom = circle; front atom = centre hub.
    p.setPen(QPen(discPen, 1.4));
    p.setBrush(Qt::NoBrush);
    p.drawEllipse(c, R, R);
    p.setBrush(hubCol);
    p.setPen(Qt::NoPen);
    p.drawEllipse(c, 3.2, 3.2);

    // Torsion arc: front reference (0) -> back reference (torsion), in the disc.
    {
        const double arcR = R * 0.44;
        const QRectF arcRect(c.x() - arcR, c.y() - arcR, 2 * arcR, 2 * arcR);
        p.setBrush(Qt::NoBrush);
        p.setPen(QPen(arcPen, 1.5));
        p.drawArc(arcRect, 0, static_cast<int>(std::lround(drawn.torsionDeg * 16.0)));
    }

    QFont lblFont = p.font();
    lblFont.setPointSizeF(lblFont.pointSizeF() * 1.06);
    p.setFont(lblFont);

    // Back spokes: from the rim outward. The reference spoke is the live needle.
    for (const NewmanSpoke& s : drawn.spokes) {
        if (s.front) continue;
        const QColor col = s.reference ? refBack : backCol;
        p.setPen(QPen(col, s.reference ? 2.8 : 1.7));
        p.drawLine(pt(R, s.angleRad), pt(R * 1.34, s.angleRad));
        p.setPen(col);
        drawCenteredText(p, pt(R * 1.52, s.angleRad), s.label);
    }
    // Front spokes: from the hub to just inside the rim.
    for (const NewmanSpoke& s : drawn.spokes) {
        if (!s.front) continue;
        const QColor col = s.reference ? refFront : frontCol;
        p.setPen(QPen(col, s.reference ? 2.8 : 1.7));
        p.drawLine(c, pt(R * 0.96, s.angleRad));
        p.setPen(col);
        drawCenteredText(p, pt(R * 0.66, s.angleRad), s.label);
    }

    // Sight-axis identity (front -> back), compact, top-left. The front atom is
    // the hub, the back atom the rim; spoke colour (front blue, back red) carries
    // the rest, so no separate legend block is needed (it only collided here).
    p.setPen(ink);
    QFont fhdr = p.font();
    fhdr.setPointSizeF(fhdr.pointSizeF() * 0.95);
    p.setFont(fhdr);
    p.drawText(area.adjusted(6, 4, -6, 0), Qt::AlignLeft | Qt::AlignTop,
               QStringLiteral("%1 %2 %3").arg(drawn.frontLabel, QString(kArrow), drawn.backLabel));

    QFont fbig = p.font();
    fbig.setPointSizeF(fbig.pointSizeF() * 1.6);
    fbig.setBold(true);
    p.setFont(fbig);
    p.setPen(readoutInk);
    p.drawText(QRectF(area.left(), area.bottom() - 30, area.width(), 26),
               Qt::AlignHCenter | Qt::AlignVCenter,
               QStringLiteral("%1 = %2%3").arg(NewmanKindGlyph(drawn.kind))
                   .arg(drawn.torsionDeg, 0, 'f', 1).arg(kDeg));
}

}  // namespace h5reader::app

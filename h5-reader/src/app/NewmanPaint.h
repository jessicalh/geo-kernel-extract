// NewmanPaint -- the reusable Newman-diagram drawing core, extracted from the
// (retired) NewmanDock so the dashboard's Newman display panel can share it.
// Pure painting of a NewmanProjection; no widget/model deps beyond the
// projection struct. Source stays ASCII (glyphs via QChar code points).

#pragma once

#include "NewmanProjection.h"

#include <QString>

#include <vector>

class QPainter;
class QRectF;

namespace h5reader::app {

// phi / psi glyph for a kind (captions / headers).
QString NewmanKindGlyph(NewmanKind kind);

// Whole-trajectory occupancy of the torsion angle: a histogram over the full
// circle [-180, 180), bin i covering [-180 + i*360/N, ...). Values are a
// relative density (the painter normalises to its own max). Empty => no haze.
struct NewmanOccupancy {
    std::vector<float> bins;  // angular density around the circle
};

// Paint `drawn` (the live, current-frame projection) into `area`, behind it a
// translucent pastel "haze" of `occupancy` (the d'Arsonval dial: the live spoke
// is the needle, the haze is where the angle lives over the trajectory). Does
// NOT fill the background (the caller owns it). darkBg selects the neutral ink
// palette. An invalid projection draws its reason.
void PaintNewman(QPainter& p, const QRectF& area, const NewmanProjection& drawn,
                 const NewmanOccupancy& occupancy, bool darkBg);

}  // namespace h5reader::app

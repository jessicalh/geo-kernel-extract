#include "RamaRegion.h"

#include <cmath>

namespace h5reader::rediscover {

namespace {

double deg(double rad) { return rad * 57.29577951308232; }

double wrapDeg(double d) {
    while (d <= -180.0) d += 360.0;
    while (d > 180.0) d -= 360.0;
    return d;
}

bool between(double v, double lo, double hi) { return v >= lo && v <= hi; }

}  // namespace

RowRamaRegion ClassifyRowRama(double phiRadians, double psiRadians) {
    if (!std::isfinite(phiRadians) || !std::isfinite(psiRadians))
        return RowRamaRegion::Undefined;
    const double phi = wrapDeg(deg(phiRadians));
    const double psi = wrapDeg(deg(psiRadians));
    if (between(phi, -170.0, -30.0) && (between(psi, 75.0, 180.0) || between(psi, -180.0, -120.0)))
        return RowRamaRegion::Beta;
    if (between(phi, -100.0, -35.0) && between(psi, 75.0, 180.0))
        return RowRamaRegion::PPII;
    if (between(phi, -170.0, -20.0) && between(psi, -90.0, 55.0))
        return RowRamaRegion::AlphaR;
    if (between(phi, 20.0, 120.0) && between(psi, -45.0, 100.0))
        return RowRamaRegion::AlphaL;
    return RowRamaRegion::Other;
}

const char* NameForRowRama(RowRamaRegion r) {
    switch (r) {
    case RowRamaRegion::AlphaR: return "alphaR";
    case RowRamaRegion::Beta: return "beta";
    case RowRamaRegion::PPII: return "PPII";
    case RowRamaRegion::AlphaL: return "alphaL";
    case RowRamaRegion::Other: return "other";
    case RowRamaRegion::Undefined: return "undefined";
    }
    return "undefined";
}

}  // namespace h5reader::rediscover

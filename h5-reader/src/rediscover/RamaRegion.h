#pragma once

namespace h5reader::rediscover {

enum class RowRamaRegion : int {
    AlphaR = 0,
    Beta = 1,
    PPII = 2,
    AlphaL = 3,
    Other = 4,
    Undefined = 5,
};

RowRamaRegion ClassifyRowRama(double phiRadians, double psiRadians);
const char* NameForRowRama(RowRamaRegion r);

}  // namespace h5reader::rediscover

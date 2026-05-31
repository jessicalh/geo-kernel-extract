// ResidentIndexes — day-one immutable indexes paired with RunData.

#pragma once

#include "RingGeometryCache.h"
#include "SpatialIndexSet.h"
#include "TemporalIndex.h"
#include "TypedAtomIndex.h"

namespace h5reader::rediscover {

class RunData;

struct ResidentIndexes {
    TypedAtomIndex typedAtoms;
    SpatialIndexSet spatial;
    RingGeometryCache ringGeometry;
    TemporalIndex temporal;
};

ResidentIndexes BuildResidentIndexes(const RunData& run);

}  // namespace h5reader::rediscover

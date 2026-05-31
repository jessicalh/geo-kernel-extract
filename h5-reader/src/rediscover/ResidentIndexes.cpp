#include "ResidentIndexes.h"

#include "RunData.h"

namespace h5reader::rediscover {

ResidentIndexes BuildResidentIndexes(const RunData& run) {
    ResidentIndexes idx;
    if (run.protein) idx.typedAtoms = TypedAtomIndex(*run.protein);
    idx.spatial = SpatialIndexSet(run);
    idx.ringGeometry = RingGeometryCache(run);
    idx.temporal = TemporalIndex(run.conformation ? run.conformation->frameCount() : 0);
    return idx;
}

}  // namespace h5reader::rediscover

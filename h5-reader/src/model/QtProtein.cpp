#include "QtProtein.h"

namespace h5reader::model {

int QtProtein::ringCountByType(RingTypeIndex t) const {
    int n = 0;
    for (const auto& ring : rings_) {
        if (ring && ring->TypeIndex() == t) ++n;
    }
    return n;
}

}  // namespace h5reader::model

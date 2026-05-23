// QtProtein implementation — ring-type counting helper.

#include "QtProtein.h"

namespace h5reader::model {

int QtProtein::ringCountByType(RingTypeIndex t) const {
    if (!topology_)
        return 0;
    int n = 0;
    for (std::size_t i = 0; i < topology_->ringCount(); ++i) {
        const QtRing& r = topology_->ringAt(i);
        if (r.TypeIndex() == t)
            ++n;
    }
    return n;
}

}  // namespace h5reader::model

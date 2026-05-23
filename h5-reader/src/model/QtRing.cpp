// QtRing factory — the ONE place RingTypeIndex → concrete subclass
// dispatch lives. All 9 ring chemistries (8 aromatic + 1 saturated).

#include "QtRing.h"

namespace h5reader::model {

std::unique_ptr<QtRing> CreateQtRing(RingTypeIndex idx) {
    switch (idx) {
    case RingTypeIndex::PheBenzene:
        return std::make_unique<QtPheBenzeneRing>();
    case RingTypeIndex::TyrPhenol:
        return std::make_unique<QtTyrPhenolRing>();
    case RingTypeIndex::TrpBenzene:
        return std::make_unique<QtTrpBenzeneRing>();
    case RingTypeIndex::TrpPyrrole:
        return std::make_unique<QtTrpPyrroleRing>();
    case RingTypeIndex::TrpPerimeter:
        return std::make_unique<QtIndolePerimeterRing>();
    case RingTypeIndex::HisImidazole:
        return std::make_unique<QtHisImidazoleRing>();
    case RingTypeIndex::HidImidazole:
        return std::make_unique<QtHidImidazoleRing>();
    case RingTypeIndex::HieImidazole:
        return std::make_unique<QtHieImidazoleRing>();
    case RingTypeIndex::ProPyrrolidine:
        return std::make_unique<QtProPyrrolidineRing>();
    }
    // Unknown ordinal — caller logs via ErrorBus.
    return nullptr;
}

}  // namespace h5reader::model

// QtBond — one covalent bond in the protein topology.
//
// Decoded from /topology at H5 load. Connectivity is static across the
// trajectory (the extractor determines bonds once at conformation 0);
// bond length, midpoint, and direction change per frame and are
// computed on demand in QtFrame, not stored here.

#pragma once

#include "Types.h"

#include <cstddef>

namespace h5reader::model {

struct QtBond {
    size_t        atomIndexA   = 0;                      // into QtProtein.atoms()
    size_t        atomIndexB   = 0;
    BondOrder     order        = BondOrder::Unknown;
    BondCategory  category     = BondCategory::Unknown;

    // Convenience query methods — typed, not string-based.
    bool IsPeptide()   const { return order == BondOrder::Peptide; }
    bool IsAromatic()  const { return order == BondOrder::Aromatic
                                    || category == BondCategory::Aromatic; }
    bool IsDisulfide() const { return category == BondCategory::Disulfide; }
    bool IsBackbone()  const {
        return category == BondCategory::PeptideCO
            || category == BondCategory::PeptideCN
            || category == BondCategory::BackboneOther;
    }
};

}  // namespace h5reader::model

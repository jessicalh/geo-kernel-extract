#pragma once
//
// GeometryResult: ring, bond, and global geometry of a conformation.
//

#include "ConformationResult.h"
#include "ProteinConformation.h"
#include "Protein.h"

namespace nmr {

class GeometryResult : public ConformationResult {
public:
    std::string Name() const override { return "GeometryResult"; }
    std::vector<std::type_index> Dependencies() const override { return {}; }

    static std::unique_ptr<GeometryResult> Compute(ProteinConformation& conf);

    const RingGeometry& RingGeometryAt(size_t ring_index) const;
    double BondLengthAt(size_t bond_index) const;
    Vec3 BondMidpointAt(size_t bond_index) const;
    Vec3 BondDirectionAt(size_t bond_index) const;

private:
    ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr

#pragma once
//
// GeometryResult: ring geometry, bond geometry, global geometry.
// Requires: nothing. First result attached to any conformation.
//

#include "ConformationResult.h"
#include "ProteinConformation.h"
#include "Protein.h"

namespace nmr {

class GeometryResult : public ConformationResult {
public:
    std::string Name() const override { return "GeometryResult"; }
    std::vector<std::type_index> Dependencies() const override { return {}; }

    // Factory: compute geometry and return ready-to-attach result.
    static std::unique_ptr<GeometryResult> Compute(ProteinConformation& conf);

    // Query methods
    const RingGeometry& RingGeometryAt(size_t ring_index) const;
    double BondLengthAt(size_t bond_index) const;
    Vec3 BondMidpointAt(size_t bond_index) const;
    Vec3 BondDirectionAt(size_t bond_index) const;

private:
    ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr

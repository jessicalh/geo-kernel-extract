#pragma once
//
// Shared geometry owner for per-(atom, aromatic-ring) rows.
//
// Ring calculators write different payload fields on RingNeighbourhood, but
// they all describe the same target/ring geometry.  This helper owns those
// shared fields so calculator execution order cannot change the azimuthal
// frame or leave a newly-created row with the default (1, 0) azimuth.
//


#include "ConformationAtom.h"
#include "Ring.h"

#include <cstddef>

namespace nmr {

RingNeighbourhood& EnsureRingNeighbourGeometry(
    ConformationAtom& atom,
    const RingGeometry& geom,
    const Ring& ring,
    size_t ring_index,
    const Vec3& atom_pos);

}  // namespace nmr

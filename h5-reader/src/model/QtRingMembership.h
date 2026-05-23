// QtRingMembership — one per-(ring, vertex-atom) record, mirroring
// ring_membership.npy (6 fields, 12 bytes per record).
//
// Each row binds an atom to a ring (the ring_id is the absolute row in
// rings.npy, atom_index references the atom axis), records the atom's
// position in the canonical ring walk (ring_atom_order), and tags
// whether the atom is a structural vertex or a substituent (currently
// always vertex; substituent is reserved for future extension per
// src/TopologySidecar.cpp:421).
//
// Ring-walk-order matters for the polygon overlay (vertices must be
// drawn in cyclic order, not arbitrary order) and for the BS/HM
// volumetric kernel evaluator (Johnson-Bovey two-loop integral
// orientation needs the canonical walk).

#pragma once

#include <cstdint>

namespace h5reader::model {

struct QtRingMembership {
    int32_t ringId = -1;       // absolute row index in rings.npy
    int32_t atomIndex = -1;    // into QtProtein.atoms()
    int8_t ringAtomOrder = 0;  // 0..ringSize-1, canonical walk order
    bool isVertex = false;
    bool isSubstituent = false;
};

}  // namespace h5reader::model

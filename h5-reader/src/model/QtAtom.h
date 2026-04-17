// QtAtom — per-atom identity + enrichment, decoded from the H5 /atoms
// and /topology groups at load time.
//
// Identity: what the atom IS. Never changes frame-to-frame.
// Position, per-frame tensors, per-frame vectors — NOT stored here.
// Those live as slab-view accessors on QtFrame.
//
// All enum fields are typed (Element, AtomRole, Hybridisation). The
// one string field — h5AtomName — is DISPLAY ONLY. No physics code
// compares atom names as strings. If you need to know whether this
// atom is an alpha-carbon, read `role == AtomRole::BackboneCA`, not
// the string. See notes/SCOPE.md, feedback_oo_modeling.

#pragma once

#include "Types.h"

#include <QString>
#include <cstddef>

namespace h5reader::model {

struct QtAtom {
    // ---- Identity (decoded from atoms/) ----
    Element         element        = Element::Unknown;
    AtomRole        role           = AtomRole::Unknown;
    Hybridisation  hybridisation   = Hybridisation::Unassigned;
    int             residueIndex   = -1;    // into QtProtein.residues()
    int             nBonded        = 0;     // number of covalent neighbours
    int             graphDistRing  = -1;    // bond-hops to nearest aromatic ring vertex
    int             graphDistN     = -1;    // bond-hops to nearest N
    int             graphDistO     = -1;    // bond-hops to nearest O

    // ---- Enrichment booleans (decoded from atoms/) ----
    bool            isBackbone          = false;
    bool            isAmideH            = false;
    bool            isAlphaH            = false;
    bool            isMethyl            = false;
    bool            isAromaticH         = false;
    bool            isOnAromaticResidue = false;
    bool            isHBondDonor        = false;
    bool            isHBondAcceptor     = false;
    bool            parentIsSp2         = false;
    bool            isConjugated        = false;

    // ---- Derived topology (decoded from topology/) ----
    int             parentAtomIndex = -1;  // for H atoms, heavy-atom parent
    double          enegSum1         = 0.0;   // Pauling units, 1-bond sum
    double          enegSum2         = 0.0;   // 2-bond sum
    int             nPiBonds3        = 0;     // pi bonds within 3 hops
    int             bfsToNearestRing = -1;    // BFS hops to nearest ring atom
    double          bfsDecay         = 0.0;   // exp decay weight

    // ---- Static charges (decoded from atoms/) ----
    double          partialCharge    = 0.0;   // ff14SB or CHARMM force field charge
    double          vdwRadius        = 0.0;   // Angstroms

    // ---- Display ----
    // The literal string the extractor wrote. CHARMM-flavoured for
    // GROMACS trajectories ("HN" for backbone amide, "HSD" variant
    // atom names for HID, etc.). DISPLAY ONLY. Never dispatch physics
    // off this value.
    QString         h5AtomName;
};

}  // namespace h5reader::model

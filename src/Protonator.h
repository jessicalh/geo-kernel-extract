#pragma once
//
// Protonator: the DECISION step for protonation.
//
// Given a protein, a conformation (for geometry), and a pH, returns a
// ProtonationState with per-residue decisions. Does NOT modify the
// protein or conformation. Does NOT add hydrogen atoms. Does NOT
// assign charges.
//
// The build step (tleap, pdb2gmx) is downstream: it takes a
// ProtonationState and produces a new Protein with hydrogen atoms
// and charges. The Protonator decides; the builder applies.
//
// Implementations:
//   PropkaProtonator  — PROPKA 3.5: pKa from structure, Henderson-Hasselbalch
//   KamlProtonator    — KaML-CBTrees: ML-predicted pKa (future)
//   TleapProtonator   — tleap default protonation per force field (future)
//

#include "ProtonationState.h"
#include <string>

namespace nmr {

class Protein;
class ProteinConformation;

struct ProtonationResult {
    ProtonationState state;
    bool ok = false;
    std::string error;

    bool Ok() const { return ok; }
};


class Protonator {
public:
    virtual ~Protonator() = default;

    // Decide protonation for this protein at this pH.
    // Reads geometry from the conformation to predict pKa values.
    // Returns a ProtonationState with per-residue decisions.
    // Does NOT modify the protein or conformation.
    virtual ProtonationResult Protonate(
        const Protein& protein,
        const ProteinConformation& conf,
        double pH) = 0;

    virtual std::string ToolName() const = 0;
    virtual std::string ToolVersion() const = 0;
};

}  // namespace nmr

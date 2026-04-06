#pragma once
//
// DsspResult: secondary structure, phi/psi, SASA from DSSP.
//
// Uses libdssp (Joosten et al. 2011, from Kabsch & Sander 1983).
// Per-residue: secondary structure char, phi (rad), psi (rad), SASA (A^2),
// H-bond acceptor and donor partners.
//
// Dependencies: none (DSSP needs only positions).
//

#include "ConformationResult.h"
#include "ProteinConformation.h"
#include <vector>

namespace nmr {

struct DsspResidue {
    char secondary_structure = 'C';  // H/G/I/E/B/T/S/C
    double phi = 0.0;                // radians
    double psi = 0.0;                // radians
    double sasa = 0.0;               // Angstroms^2

    // H-bond partner indices (into protein residue list, SIZE_MAX if none)
    struct HBondPartner {
        size_t residue_index = SIZE_MAX;
        double energy = 0.0;         // kcal/mol
    };
    HBondPartner acceptors[2];       // residues whose C=O accepts this N-H
    HBondPartner donors[2];          // residues whose N-H donates to this C=O
};


class DsspResult : public ConformationResult {
public:
    std::string Name() const override { return "DsspResult"; }
    std::vector<std::type_index> Dependencies() const override { return {}; }

    // Factory: run DSSP on the conformation and return the result.
    // Returns nullptr on failure (writes diagnostic to stderr).
    static std::unique_ptr<DsspResult> Compute(ProteinConformation& conf);

    // Physics query methods
    char SecondaryStructure(size_t residue_index) const;
    double Phi(size_t residue_index) const;
    double Psi(size_t residue_index) const;
    double SASA(size_t residue_index) const;

    // All per-residue data
    const std::vector<DsspResidue>& AllResidues() const { return residues_; }

private:
    std::vector<DsspResidue> residues_;
};

}  // namespace nmr

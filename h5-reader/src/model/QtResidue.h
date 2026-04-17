// QtResidue — one amino acid at a sequence position.
//
// Decoded from the H5 /residues group plus a pass over /atoms to
// build the backbone index cache. Protonation variant is inferred
// from ring subclass identity (HIS tautomers) and atom presence
// (ASP/GLU/LYS/TYR/CYS) — never from residue-name string comparison.
//
// Access patterns downstream: typed index into atoms via `.N`, `.CA`,
// `.C`, `.O`, `.H`, `.HA`, `.CB`. SIZE_MAX means absent (PRO has no
// amide H; GLY has no CB; terminal residues drop N or C).

#pragma once

#include "Types.h"

#include <QString>
#include <climits>
#include <cstddef>
#include <vector>

namespace h5reader::model {

struct QtResidue {
    // Identity
    AminoAcid          aminoAcid    = AminoAcid::Unknown;
    int                residueNumber = 0;         // PDB sequence number
    QString            chainId;                   // e.g. "A" (addressing, not identity)

    // Protonation state — inferred typed, not parsed from strings.
    ProtonationVariant variant = ProtonationVariant::Default;

    // Atom membership. Indices into QtProtein.atoms().
    std::vector<size_t> atomIndices;

    // Backbone index cache — SIZE_MAX means absent.
    // Populated at load time from AtomRole matching:
    //   BackboneN  → N,  BackboneCA → CA,  BackboneC → C,  BackboneO → O,
    //   AmideH     → H,  AlphaH     → HA (or HA2 for GLY).
    // CB is the first non-aromatic sidechain C encountered.
    static constexpr size_t NONE = SIZE_MAX;
    size_t N  = NONE;
    size_t CA = NONE;
    size_t C  = NONE;
    size_t O  = NONE;
    size_t H  = NONE;
    size_t HA = NONE;
    size_t CB = NONE;

    bool HasAmideH() const { return H != NONE; }
    bool IsTerminalN() const { return N == NONE; }   // missing N → N-terminal fragment
    bool IsTerminalC() const { return C == NONE; }
};

}  // namespace h5reader::model

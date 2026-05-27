#pragma once
//
// ReduceProtonation: thin wrapper around Richardson lab reduce (reducelib).
//
// Adds hydrogen atoms to a PDB structure with PDB-standard naming.
// This is the ONLY place reduce internals are touched. The rest of
// our codebase sees PDB strings in and protonated PDB strings out.
//
// reduce handles:
//   - HIS side-chain proton placement (HD1/HE2 pattern from optimization)
//   - NQH flip optimization (ASN, GLN, HIS ring flips)
//   - OH/SH hydrogen placement and rotation
//   - Standard PDB hydrogen naming (H, HA, HD1, HE2, etc.)
//
// Protein::ResolveProtonationStates recognises reduce's output without
// modification because it uses PDB-standard H atom names.
//

#include <string>

namespace nmr {

// Protonate a PDB structure using reduce (Richardson lab).
//
// Input:  PDB file content as string (heavy atoms, standard naming)
// Output: protonated PDB content as string (heavy + H atoms, standard naming)
//
// Returns empty string on failure.
// The het dictionary path is resolved from REDUCE_HET_DICT or the CMake default.
std::string ProtonateWithReduce(const std::string& pdb_content);

}  // namespace nmr

#pragma once
//
// PDB file reader: parses PDB/mmCIF into a fully constructed Protein.
//
// BuildFromPdb is the ONLY entry point for PDB files. It protonates
// with reduce (Richardson lab) and assigns ff14SB charges. Every
// protein in this system is protonated and charged — these are
// preconditions for the physics we compute.
//
// Uses cif++ for parsing. AminoAcidType provides ring definitions.
// reduce adds H atoms. Atom-name strings from cif++ are assigned raw
// to Atom::pdb_atom_name — no NamingRegistry translation on this
// path. AMBER-convention input is assumed (the documented contract
// for --pdb is reduce-protonated PDBs, which emit AMBER names). See
// KNOWN_BUGS.md Bug 1 for what fails when this assumption is wrong.
//

#include "Protein.h"
#include "BuildResult.h"
#include <string>

namespace nmr {

// Build a fully protonated, charge-assigned protein from a bare PDB.
// Protonates with reduce at the given pH, assigns ff14SB charges.
// The PDB is assumed to NOT be protonated. reduce strips any existing
// H atoms and rebuilds them cleanly.
BuildResult BuildFromPdb(const std::string& path, double pH = 7.0);

// Build from an already-protonated PDB (e.g., reduce output saved to disk,
// or a PDB known to contain all H atoms). Skips reduce, assigns charges.
// Use this for test data and pre-processed inputs.
BuildResult BuildFromProtonatedPdb(const std::string& path);

}  // namespace nmr

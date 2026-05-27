#pragma once
//
// PDB file reader: parses PDB input via cif++ into a fully constructed Protein.
//
// BuildFromPdb protonates with reduce (Richardson lab) and assigns
// ff14SB charges. BuildFromProtonatedPdb skips reduce for already-
// protonated input. Every protein from this loader is protonated and
// charged — these are preconditions for the physics we compute.
//
// Uses cif++ for parsing. reduce adds H atoms. Atom names from cif++
// pass through the NamingApplicator before they are stored on Atom.
//

#include "Protein.h"
#include "BuildResult.h"
#include <string>

namespace nmr {

// Build a fully protonated, charge-assigned protein from a bare PDB.
// Protonates with reduce, records the requested pH in the build context,
// and assigns ff14SB charges.
// The PDB is assumed to NOT be protonated. reduce strips any existing
// H atoms and rebuilds them cleanly.
BuildResult BuildFromPdb(const std::string& path, double pH = 7.0);

// Build from an already-protonated PDB (e.g., reduce output saved to disk,
// or a PDB known to contain all H atoms). Skips reduce, assigns charges.
// Use this for test data and pre-processed inputs.
BuildResult BuildFromProtonatedPdb(const std::string& path);

}  // namespace nmr

#pragma once
//
// PDB reader. BuildFromPdb protonates with reduce before parsing;
// BuildFromProtonatedPdb parses the input as-is. Both paths assign ff14SB
// charges through ResolveAmberChargeSource.
//

#include "Protein.h"
#include "BuildResult.h"
#include <string>

namespace nmr {

// reduce strips existing ATOM hydrogens before rebuilding them.
BuildResult BuildFromPdb(const std::string& path, double pH = 7.0);

// Skips reduce and parses the input directly.
BuildResult BuildFromProtonatedPdb(const std::string& path);

}  // namespace nmr

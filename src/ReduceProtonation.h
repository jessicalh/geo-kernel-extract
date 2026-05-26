#pragma once
//
// Wrapper around reducelib. Input and output are whole PDB files as strings;
// reduce globals are configured in the implementation before each call.
//

#include <string>

namespace nmr {

// Returns empty string on failure.
std::string ProtonateWithReduce(const std::string& pdb_content);

}  // namespace nmr

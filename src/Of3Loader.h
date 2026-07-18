#pragma once
//
// Of3Loader: load a protein from an OpenFold pose prepared with tleap/ff14SB.
//
// This is the --of3 mode's structural loader. It is a peer of the ORCA
// structural loader (OrcaRunLoader.h / BuildFromOrca) and shares the exact
// same prmtop heavy-lifting parse, but reads geometry from the AMBER restart
// (rst7/inpcrd) that tleap emits beside the prmtop, carries HONEST OpenFold
// provenance, and has NO DFT/NMR field.
//
// Required:
//   - inpcrd: AMBER rst7 coordinates from the tleap/ff14SB preparation
//     (tleap writes this next to the prmtop, so no xyz conversion is needed)
//   - prmtop: AMBER topology (atom names, residues, elements, charges/radii)
//
// Optional/provenance:
//   - PDB: upstream structure path, recorded on ProteinBuildContext
//   - tleap script: recorded on ProteinBuildContext
//   - prediction_method: REQUIRED conformation provenance, supplied at the
//     mode boundary ("OpenFold+tleap"). BuildFromOf3 fails clearly if it is
//     empty rather than silently recording false provenance.
//
// Deliberately there is NO nmr_out_path: DFT is not part of OF3's object
// shape. A co-located {root}_nmr.out is irrelevant and never opened.
//
// Loading path:
//   prmtop provides the protonated atom list (names, residues, element);
//   the inpcrd provides positions in the same (prmtop) atom order. Missing or
//   unreadable prmtop is a hard load error; there is no canonical-residue /
//   flat-ff14SB-table fallback. Charges resolve through the shared Branch-1
//   resolver over the input prmtop (PrmtopChargeSource); the returned
//   BuildResult::charges is Branch-1 prmtop authority and no
//   runtime-preparation fallback is permitted for a missing prmtop.
//
// The definition lives in OrcaRunLoader.cpp so that OF3 reuses the existing
// file-local heavy parse (LoadWithPrmtop / prmtop section readers) verbatim
// rather than duplicating it, alongside a small file-local rst7 reader. The
// ORCA entry point (and its xyz reader) is untouched.
//

#include "BuildResult.h"

#include <string>

namespace nmr {

struct Of3Input {
    std::string pdb_path;           // Upstream PDB path (provenance, optional)
    std::string inpcrd_path;        // AMBER rst7/inpcrd coordinates (required)
    std::string prmtop_path;        // AMBER topology + charges/radii (required)
    std::string tleap_script_path;  // tleap input script (provenance, optional)
    std::string prediction_method;  // conformation provenance (required, "OpenFold+tleap")
    // NO nmr_out_path: OF3 carries no DFT/NMR field by construction.
};

// Load a protonated Protein from an OpenFold/tleap-prepared pose.
//
// Requires prmtop_path (topology) and inpcrd_path (rst7 positions, prmtop
// atom order). Charges from PrmtopChargeSource (Branch 1). Net charge from
// prmtop charge sum. Requires a non-empty prediction_method (honest
// provenance guard).
//
// The resulting Protein has one PredictionConformation with the inpcrd
// positions, tagged with input.prediction_method. Charges are wrapped in the
// BuildResult's ChargeSource. No OrcaShieldingResult is ever attached.
BuildResult BuildFromOf3(const Of3Input& input);

}  // namespace nmr

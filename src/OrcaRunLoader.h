#pragma once
//
// OrcaRunLoader: load a protein from an ORCA/tleap-prepared pose.
//
// Required:
//   - XYZ: protonated coordinates from tleap/ORCA geometry
//   - prmtop: AMBER topology (atom names, residues, elements, charges)
//
// Optional/provenance:
//   - PDB: upstream structure path, recorded on ProteinBuildContext
//   - NMR output: per-atom shielding tensors (loaded separately)
//   - tleap artifacts: script, log, amber PDB, inpcrd
//
// Loading path:
//   prmtop provides the protonated atom list (names, residues, element);
//   XYZ provides positions. Missing or unreadable prmtop is a hard load
//   error; there is no canonical-residue fallback on ORCA paths.
//
// BuildFromOrca produces a Protein with the full protonated atom list and
// one conformation with XYZ positions. Charges and NMR tensors attach
// separately as ConformationResults from the supplied ORCA-path files.
//

#include "Protein.h"
#include "BuildResult.h"
#include "PdbFileReader.h"
#include <string>

namespace nmr {

struct OrcaRunFiles {
    std::string pdb_path;           // Upstream PDB path (provenance)
    std::string xyz_path;           // tleap-protonated coordinates (required)
    std::string prmtop_path;        // AMBER topology (required)
    std::string nmr_out_path;       // ORCA NMR shielding output (for OrcaShieldingResult)
    std::string tleap_script_path;  // tleap input script (provenance, optional)
};


// Load a protonated Protein from an ORCA/tleap-prepared pose.
//
// Requires prmtop_path: topology from prmtop, positions from XYZ.
// Charges from PrmtopChargeSource. Net charge from prmtop charge sum.
//
// The resulting Protein has one PredictionConformation with XYZ positions.
// Charges are wrapped in the BuildResult's ChargeSource.
BuildResult BuildFromOrca(const OrcaRunFiles& files);

}  // namespace nmr

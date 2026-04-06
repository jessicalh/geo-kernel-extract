#pragma once
//
// OrcaRunLoader: load a protein from an ORCA DFT run.
//
// An ORCA run always has:
//   - PDB: AlphaFold/crystal structure (heavy atoms, sequence identity)
//   - XYZ: protonated coordinates from tleap (all atoms, ORCA geometry)
//
// It may also have:
//   - prmtop: AMBER topology (authoritative atom names, residues, charges)
//   - NMR output: per-atom shielding tensors (loaded separately)
//   - tleap artifacts: script, log, amber PDB, inpcrd
//
// Two loading paths:
//   WITH prmtop (all 723 proteins): prmtop provides the protonated atom list
//     (names, residues, element). XYZ provides positions.
//     As of 2026-04-02, prmtop regenerated via tLeap ff14SB for the 253
//     proteins that were previously missing them. All proteins now use
//     this path.
//   WITHOUT prmtop (legacy fallback): PDB provides sequence. AminoAcidType
//     provides the canonical atom list per residue (including H). XYZ
//     provides positions. Atoms matched by residue structure and element.
//
// Both paths produce: a Protein with the full protonated atom list and
// one conformation with XYZ positions. Charges and NMR tensors attach
// separately as ConformationResults from whatever source is available.
//

#include "Protein.h"
#include "BuildResult.h"
#include "PdbFileReader.h"
#include <string>

namespace nmr {

struct OrcaRunFiles {
    std::string pdb_path;           // AlphaFold/crystal PDB (sequence identity)
    std::string xyz_path;           // tleap-protonated coordinates (required)
    std::string prmtop_path;        // AMBER topology (optional — empty if not available)
    std::string nmr_out_path;       // ORCA NMR shielding output (for OrcaShieldingResult)
    std::string tleap_script_path;  // tleap input script (provenance, optional)
};


// Load a protonated Protein from an ORCA run.
//
// If prmtop_path is provided: topology from prmtop, positions from XYZ.
// Charges from PrmtopChargeSource. Net charge from prmtop charge sum.
//
// The resulting Protein has one PredictionConformation with XYZ positions.
// Charges are wrapped in the BuildResult's ChargeSource.
BuildResult BuildFromOrca(const OrcaRunFiles& files);

}  // namespace nmr

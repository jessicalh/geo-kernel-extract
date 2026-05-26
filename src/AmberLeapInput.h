#pragma once
//
// Pure builders for the PDB body, LEaP script body, and extractor-to-PRMTOP
// residue mapping used by AmberPreparedChargeSource. They perform no I/O.
//

#include "AmberChargeResolver.h"
#include <cstddef>
#include <iosfwd>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace nmr {

class Protein;
class ProteinConformation;

namespace amber_leap {

// Per generated-PRMTOP residue: which extractor residue it corresponds
// to, or NONE_FOR_CAP if it is an inserted ACE/NME/NHE cap.
struct ResidueAmberMapping {
    static constexpr size_t NONE_FOR_CAP =
        std::numeric_limits<size_t>::max();
    std::vector<size_t> extractor_index_for_prmtop_residue;
};


struct LeapScriptInputs {
    std::string pdb_path;
    std::string prmtop_path;
    std::string inpcrd_path;
    // 1-based PRMTOP residue indices for SG-SG bonds. Built by the
    // caller from DetectDisulfides + the ResidueAmberMapping.
    std::vector<std::pair<size_t, size_t>> disulfide_residue_pairs_1based;
};


// Detect inter-CYS disulfide bonds by direct SG-SG distance check on the
// typed CYS residues. Returns extractor residue index pairs (0-based,
// each with a < b).
//
// Any pair of CYS residues whose SG atoms are within
// max_ss_distance_angstroms of each other (default 2.5 Angstrom) is renamed
// CYX and bonded explicitly in the LEaP script. It does not depend on
// CovalentTopology::Resolve or any external bond-perception library;
// AmberPreparedChargeSource is methodologically self-contained for
// AMBER charge preparation.
std::vector<std::pair<size_t, size_t>> DetectDisulfides(
    const Protein& protein,
    const ProteinConformation& conf,
    double max_ss_distance_angstroms = 2.5);


// Walks the Protein's residues in order; each residue's name is the
// AMBER unit name derived from typed state:
//   - HIS variant 0/1/2 → HID/HIE/HIP
//   - HIS no variant → HIE (ff14SB has no neutral generic HIS)
//   - CYS in a Disulfide bond → CYX (overrides any other variant)
//   - CYS variant 1 → CYM
//   - ASP variant 0 → ASH; GLU variant 0 → GLH; LYS variant 0 → LYN
//   - others → canonical three-letter
// Atom records are emitted in extractor atom-index order. ACE/NME caps under
// UseCappedFragmentsForUnsupportedTerminalVariants are marked NONE_FOR_CAP
// in ResidueAmberMapping.
void GenerateAmberPdb(const Protein& protein,
                      const ProteinConformation& conf,
                      AmberPreparationPolicy policy,
                      const AmberFlatTableCoverageVerdict& verdict,
                      std::ostream& pdb_out,
                      ResidueAmberMapping& map_out);


// Emits a LEaP input script. Body is exactly:
//
//   source leaprc.protein.ff14SB
//   set default PBRadii mbondi2
//   mol = loadPdb <pdb_path>
//   bond mol.<ri>.SG mol.<rj>.SG    (one line per disulfide pair)
//   saveamberparm mol <prmtop_path> <inpcrd_path>
//   quit
void GenerateLeapScript(const LeapScriptInputs& inputs,
                        std::ostream& script_out);

}  // namespace amber_leap
}  // namespace nmr

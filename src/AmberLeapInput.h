#pragma once
//
// AmberLeapInput: builders for deterministic LEaP input artifacts
// (PDB body and LEaP script body) from typed Protein state. No I/O of
// their own — callers stream output to a file.
//
// These functions are the testable boundary for AmberPreparedChargeSource:
//   GenerateAmberPdb       — typed PDB generator (residue naming derives
//                            from protonation_variant_index;
//                            disulfide-CYS detected via DetectDisulfides
//                            and emitted as CYX; ACE/NME caps inserted
//                            under UseCappedFragmentsForUnsupportedTerminalVariants
//                            for unsupported terminal variants).
//   GenerateLeapScript     — leaprc.protein.ff14SB + mbondi2 +
//                            loadPdb + bond mol.<ri>.SG mol.<rj>.SG +
//                            saveamberparm. Pure function on its inputs.
//   DetectDisulfides       — typed pairwise SG-SG distance check on CYS
//                            residues.
//   ResidueAmberMapping    — extractor↔PRMTOP residue index map; cap
//                            residues marked NONE_FOR_CAP.
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

// Per generated-PRMTOP residue: which extractor residue does it correspond
// to, or NONE_FOR_CAP if it is an inserted ACE/NME cap.
//
// Built during GenerateAmberPdb; consumed by the atom-mapping
// step of AmberPreparedChargeSource::LoadCharges.
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
// `max_ss_distance_angstroms` of each other (default 2.5 Å) is renamed
// CYX and bonded explicitly in the LEaP script. It does not depend on
// CovalentTopology::Resolve or any external bond-perception library;
// AmberPreparedChargeSource is methodologically self-contained for
// AMBER charge preparation.
std::vector<std::pair<size_t, size_t>> DetectDisulfides(
    const Protein& protein,
    const ProteinConformation& conf,
    double max_ss_distance_angstroms = 2.5);


// Emits a PDB body suitable for `loadPdb` in tleap.
//
// Walks the Protein's residues in order; each residue's name is the
// AMBER unit name derived from typed state:
//   - HIS variant 0/1/2 → HID/HIE/HIP
//   - HIS no variant → HIE
//   - CYS in a Disulfide bond → CYX (overrides any other variant)
//   - CYS variant 1 → CYM
//   - ASP variant 0 → ASH; GLU variant 0 → GLH; LYS variant 0 → LYN
//   - others → canonical three-letter
//
// Atom records emitted in extractor atom-index order. TER between chains.
//
// Without caps, every extractor residue maps 1:1 to a PRMTOP residue.
// Under UseCappedFragmentsForUnsupportedTerminalVariants, ACE/NME cap
// residues may be inserted and marked NONE_FOR_CAP.
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

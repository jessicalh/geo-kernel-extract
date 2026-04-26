#pragma once
//
// CcdValidator: read-only sanity check of a finalized Protein against the
// PDB Chemical Component Dictionary (CCD) via gemmi. Two roles:
//
//   DescribeAtom(residue_three_letter, atom_name) — single-atom CCD lookup.
//     Used by Protein::StampAtomTopology to enrich the per-atom diagnostic
//     when an atom name is not in our IUPAC table, with what CCD says about
//     it. PDB string boundary stays at the parameters; the typed flag
//     (CcdAtomDescription::found) is the signal calculator code reads.
//
//   Check(protein) — whole-protein validator. Walks every residue and atom,
//     records findings (atom not in CCD, atom missing from Protein, element
//     mismatch, residue type unknown to CCD). Result is a CcdValidationReport
//     that can be written to summary.npy + findings.csv.
//
// Both roles share a process-global CcdStore that lazy-loads
// data/ccd/components.cif.gz on first access. CCD path is resolved by
// RuntimeEnvironment::CcdPath() — TOML > env var > NMR_DATA_DIR default.
//
// No data flows from gemmi into Atom::topology or any other typed field.
// The validator enriches diagnostics. AminoAcidType remains the runtime
// authority for symbolic identity.
//

#include <string>
#include <vector>

namespace nmr {

class Protein;

// Description of one atom according to CCD. Strings only here — calculator
// code checks `found` and reads typed fields.
struct CcdAtomDescription {
    std::string atom_name;        // CCD canonical name (echoed for diagnostics)
    std::string element_symbol;   // "C", "N", "H", "O", ...
    bool        is_aromatic = false;
    bool        leaving_atom = false;
    bool        found = false;    // false = no CCD entry for this (residue, atom)
};

// One disagreement between Protein-side topology and CCD.
struct CcdValidationFinding {
    enum Kind {
        AtomNotInCcd,             // Protein has an atom name CCD does not list
        AtomMissingFromProtein,   // CCD lists an atom name the Protein lacks
        ElementMismatch,          // names match, elements disagree
        ResidueTypeUnknownToCcd   // CCD has no entry for this three-letter code
    };
    Kind        kind = AtomNotInCcd;
    std::string residue_three_letter;
    int         residue_position = 0;
    std::string chain_id;
    std::string atom_name;        // empty for residue-level findings
    std::string detail;           // human-readable detail
};

struct CcdValidationReport {
    int n_atoms_checked = 0;
    int n_residues_checked = 0;
    std::string ccd_path;          // resolved path, for diagnostics
    std::vector<CcdValidationFinding> findings;
    // Note: load is an environment requirement — failure to load aborts
    // before any report is constructed (see CcdStore::Load). A report
    // exists iff the CCD was successfully loaded.
};

class CcdValidator {
public:
    // Lookup a single atom in CCD. Returns CcdAtomDescription{found=false}
    // when the residue is not in CCD or the atom name is not listed for it.
    static CcdAtomDescription DescribeAtom(const std::string& residue_three_letter,
                                            const std::string& atom_name);

    // Walk a finalized Protein, comparing each atom against CCD records.
    // Read-only: never mutates the Protein. Logs each finding via
    // OperationLog::Warning at emit time; the returned report is the
    // structured form for NPY/CSV emission.
    static CcdValidationReport Check(const Protein& protein);
};

// Emit summary.npy + findings.csv to <output_dir>/ccd_validation/.
// Returns the number of files written (0 on failure).
int WriteCcdValidationDiagnostics(const CcdValidationReport& report,
                                   const std::string& output_dir);

}  // namespace nmr

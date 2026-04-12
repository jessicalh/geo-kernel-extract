#pragma once
//
// JobSpec: command-line to library handoff for all use cases.
//
// ParseJobSpec reads argc/argv into a JobSpec struct. ValidateJobSpec
// checks that every required file exists and reports what is missing
// in plain language. Both nmr_extract and nmr-viewer call these.
//
// Root-name convention (--orca, --mutant):
//   --root /path/to/A0A7C5FAR6_WT expands to:
//     {root}.xyz        — required
//     {root}.prmtop     — required
//     {root}_nmr.out    — optional (DFT tensors)
//
// The validator tells the caller exactly what was found and what
// was not, with full paths.  No regex. No directory scanning.
//

#include "OrcaRunLoader.h"       // OrcaRunFiles
#include "GromacsEnsembleLoader.h" // FleetPaths
#include <string>
#include <vector>

namespace nmr {

enum class JobMode {
    Pdb,             // bare PDB — protonate with reduce
    ProtonatedPdb,   // already-protonated PDB — detect state from H atoms
    Orca,            // single ORCA DFT (root → .xyz + .prmtop + optional _nmr.out)
    Mutant,          // WT + ALA ORCA pair (two roots)
    Fleet,           // GROMACS ensemble (TPR + poses dir)
    None             // parse failed or --help
};


struct JobSpec {
    JobMode mode = JobMode::None;

    // -- PDB modes --
    std::string pdb_path;
    double pH = 7.0;

    // -- ORCA mode --
    OrcaRunFiles orca_files;

    // -- Mutant mode --
    OrcaRunFiles wt_files;
    OrcaRunFiles ala_files;

    // -- Fleet mode --
    FleetPaths fleet_paths;

    // -- Common --
    std::string output_dir;    // empty = no feature output (viewer mode)
    std::string config_path;   // TOML calculator parameter overrides
    bool skip_mopac   = false;   // --no-mopac:   skip PM7+MOZYME (and derived calcs)
    bool skip_apbs    = false;   // --no-apbs:    skip Poisson-Boltzmann
    bool skip_coulomb = false;   // --no-coulomb:  skip vacuum Coulomb EFG (APBS preferred)
    std::string aimnet2_model_path;  // --aimnet2 MODEL: path to .jpt file (empty = skip)

    // -- Diagnostics --
    std::vector<std::string> warnings;   // non-fatal (e.g., "NMR .out not found")
    std::string error;                   // fatal parse error

    bool Ok() const { return error.empty() && mode != JobMode::None; }
};


// Parse command-line arguments into a JobSpec.
// On parse failure, error is set and mode is None.
// Does NOT check file existence — that is ValidateJobSpec's job.
JobSpec ParseJobSpec(int argc, char* argv[]);


// Check that every required file exists.  Returns false and populates
// spec.error with a specific diagnostic on the first fatal problem.
// Optional files that are absent get a warning in spec.warnings.
// Call AFTER ParseJobSpec, BEFORE building.
bool ValidateJobSpec(JobSpec& spec);


// Print usage to stderr.
void PrintJobSpecUsage(const char* program_name);


}  // namespace nmr

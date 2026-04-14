#include "JobSpec.h"
#include "OperationLog.h"
#include "CalculatorConfig.h"

#include <cstring>
#include <cstdio>
#include <filesystem>

namespace fs = std::filesystem;

namespace nmr {

// ============================================================================
// Argv helpers
// ============================================================================

static bool HasFlag(int argc, char* argv[], const char* name) {
    for (int i = 1; i < argc; ++i)
        if (std::strcmp(argv[i], name) == 0) return true;
    return false;
}

static std::string GetArg(int argc, char* argv[], const char* name) {
    for (int i = 1; i < argc - 1; ++i)
        if (std::strcmp(argv[i], name) == 0) return argv[i + 1];
    return "";
}


// ============================================================================
// Root expansion: deterministic, no regex, no scanning
//
//   root → root.xyz        (required)
//   root → root.prmtop     (required)
//   root → root_nmr.out    (optional — DFT tensors)
// ============================================================================

static OrcaRunFiles ExpandOrcaRoot(const std::string& root) {
    OrcaRunFiles f;
    f.xyz_path     = root + ".xyz";
    f.prmtop_path  = root + ".prmtop";
    f.nmr_out_path = root + "_nmr.out";
    return f;
}


// ============================================================================
// ParseJobSpec
// ============================================================================

JobSpec ParseJobSpec(int argc, char* argv[]) {
    JobSpec spec;

    if (argc < 2 || HasFlag(argc, argv, "--help") || HasFlag(argc, argv, "-h")) {
        spec.mode = JobMode::None;
        return spec;
    }

    // Common optional flags
    spec.output_dir  = GetArg(argc, argv, "--output");
    spec.config_path = GetArg(argc, argv, "--config");
    spec.skip_mopac   = HasFlag(argc, argv, "--no-mopac");
    spec.skip_apbs    = HasFlag(argc, argv, "--no-apbs");
    spec.skip_coulomb = HasFlag(argc, argv, "--no-coulomb");
    spec.aimnet2_model_path = GetArg(argc, argv, "--aimnet2");

    // ---- Mode dispatch ----

    if (HasFlag(argc, argv, "--pdb")) {
        spec.mode = JobMode::Pdb;
        spec.pdb_path = GetArg(argc, argv, "--pdb");
        std::string pH_str = GetArg(argc, argv, "--pH");
        if (!pH_str.empty()) spec.pH = std::atof(pH_str.c_str());

        if (spec.pdb_path.empty()) {
            spec.error = "--pdb requires a file path: --pdb FILE";
        }
        return spec;
    }

    if (HasFlag(argc, argv, "--protonated-pdb")) {
        spec.mode = JobMode::ProtonatedPdb;
        spec.pdb_path = GetArg(argc, argv, "--protonated-pdb");

        if (spec.pdb_path.empty()) {
            spec.error = "--protonated-pdb requires a file path: --protonated-pdb FILE";
        }
        return spec;
    }

    if (HasFlag(argc, argv, "--orca")) {
        spec.mode = JobMode::Orca;
        std::string root = GetArg(argc, argv, "--root");

        if (root.empty()) {
            spec.error = "--orca requires --root NAME\n"
                         "  root expands to: {root}.xyz, {root}.prmtop, {root}_nmr.out";
            return spec;
        }
        spec.orca_files = ExpandOrcaRoot(root);
        return spec;
    }

    if (HasFlag(argc, argv, "--mutant")) {
        spec.mode = JobMode::Mutant;
        std::string wt_root  = GetArg(argc, argv, "--wt");
        std::string ala_root = GetArg(argc, argv, "--ala");

        if (wt_root.empty() || ala_root.empty()) {
            spec.error = "--mutant requires --wt NAME --ala NAME\n"
                         "  each root expands to: {root}.xyz, {root}.prmtop, {root}_nmr.out";
            return spec;
        }
        spec.wt_files  = ExpandOrcaRoot(wt_root);
        spec.ala_files = ExpandOrcaRoot(ala_root);
        return spec;
    }

    if (HasFlag(argc, argv, "--fleet")) {
        spec.error = "--fleet mode removed (2026-04-12). Use --trajectory for "
                     "GROMACS ensemble extraction: --trajectory --tpr FILE --xtc FILE";
        return spec;
    }

    if (HasFlag(argc, argv, "--trajectory")) {
        spec.mode = JobMode::Trajectory;
        spec.traj_tpr = GetArg(argc, argv, "--tpr");
        spec.traj_xtc = GetArg(argc, argv, "--xtc");
        spec.traj_edr = GetArg(argc, argv, "--edr");
        spec.analysis = HasFlag(argc, argv, "--analysis");

        if (spec.traj_tpr.empty()) {
            spec.error = "--trajectory requires --tpr FILE (full-system topology)";
            return spec;
        }
        if (spec.traj_xtc.empty()) {
            spec.error = "--trajectory requires --xtc FILE (full-system trajectory)";
            return spec;
        }
        return spec;
    }

    spec.error = std::string("unknown mode: ") + argv[1] +
                 "\n  valid modes: --pdb, --protonated-pdb, --orca, --mutant, --trajectory";
    return spec;
}


// ============================================================================
// ValidateJobSpec — check every file, report what is missing
// ============================================================================

// Check one required file.  Returns false and sets spec.error on failure.
static bool RequireFile(JobSpec& spec, const std::string& path,
                        const char* description) {
    if (path.empty()) {
        spec.error = std::string(description) + ": no path specified";
        OperationLog::Error("JobSpec", spec.error);
        return false;
    }
    if (!fs::exists(path)) {
        spec.error = std::string(description) + ": file not found: " + path;
        OperationLog::Error("JobSpec", spec.error);
        return false;
    }
    OperationLog::Info(LogFileIO, "JobSpec",
        std::string(description) + ": " + path + " (ok)");
    return true;
}

// Check one required directory.
static bool RequireDir(JobSpec& spec, const std::string& path,
                       const char* description) {
    if (path.empty()) {
        spec.error = std::string(description) + ": no path specified";
        OperationLog::Error("JobSpec", spec.error);
        return false;
    }
    if (!fs::is_directory(path)) {
        spec.error = std::string(description) + ": directory not found: " + path;
        OperationLog::Error("JobSpec", spec.error);
        return false;
    }
    OperationLog::Info(LogFileIO, "JobSpec",
        std::string(description) + ": " + path + " (ok)");
    return true;
}

// Check one optional file.  Logs a warning if absent, clears the path.
static void OptionalFile(JobSpec& spec, std::string& path,
                         const char* description) {
    if (path.empty()) return;
    if (!fs::exists(path)) {
        std::string msg = std::string(description) + ": not found: " + path +
                          " — will proceed without it";
        spec.warnings.push_back(msg);
        OperationLog::Warn("JobSpec", msg);
        path.clear();
    } else {
        OperationLog::Info(LogFileIO, "JobSpec",
            std::string(description) + ": " + path + " (ok)");
    }
}

// Validate an ORCA root's expanded files.
static bool ValidateOrcaFiles(JobSpec& spec, OrcaRunFiles& files,
                              const char* label) {
    char desc[128];

    std::snprintf(desc, sizeof(desc), "%s XYZ coordinates", label);
    if (!RequireFile(spec, files.xyz_path, desc)) return false;

    std::snprintf(desc, sizeof(desc), "%s AMBER prmtop (charge source)", label);
    if (!RequireFile(spec, files.prmtop_path, desc)) return false;

    std::snprintf(desc, sizeof(desc), "%s ORCA NMR shielding tensors", label);
    OptionalFile(spec, files.nmr_out_path, desc);

    return true;
}


bool ValidateJobSpec(JobSpec& spec) {
    if (!spec.Ok()) return false;

    OperationLog::Info(LogFileIO, "JobSpec",
        "validating " + std::string(
            spec.mode == JobMode::Pdb           ? "--pdb" :
            spec.mode == JobMode::ProtonatedPdb  ? "--protonated-pdb" :
            spec.mode == JobMode::Orca           ? "--orca" :
            spec.mode == JobMode::Mutant         ? "--mutant" :
            spec.mode == JobMode::Trajectory     ? "--trajectory" : "unknown") +
        " job");

    // Config TOML — optional
    if (!spec.config_path.empty()) {
        if (!fs::exists(spec.config_path)) {
            spec.error = "calculator config TOML not found: " + spec.config_path;
            OperationLog::Error("JobSpec", spec.error);
            return false;
        }
        OperationLog::Info(LogFileIO, "JobSpec",
            "config TOML: " + spec.config_path + " (ok)");
    }

    // AIMNet2 model — optional, but if specified must exist
    if (!spec.aimnet2_model_path.empty()) {
        if (!RequireFile(spec, spec.aimnet2_model_path,
                         "AIMNet2 .jpt model")) return false;
        OperationLog::Info(LogFileIO, "JobSpec",
            "AIMNet2 model: " + spec.aimnet2_model_path + " (ok)");
    }

    // Output directory — create if needed (non-fatal if empty = viewer mode)
    if (!spec.output_dir.empty()) {
        OperationLog::Info(LogFileIO, "JobSpec",
            "output directory: " + spec.output_dir);
    }

    switch (spec.mode) {

    case JobMode::Pdb:
        return RequireFile(spec, spec.pdb_path, "PDB file (will protonate with reduce)");

    case JobMode::ProtonatedPdb:
        return RequireFile(spec, spec.pdb_path,
                           "pre-protonated PDB (protonation detected from H atoms)");

    case JobMode::Orca:
        return ValidateOrcaFiles(spec, spec.orca_files, "ORCA");

    case JobMode::Mutant:
        if (!ValidateOrcaFiles(spec, spec.wt_files, "WT")) return false;
        if (!ValidateOrcaFiles(spec, spec.ala_files, "ALA")) return false;
        return true;

    case JobMode::Trajectory:
        if (!RequireFile(spec, spec.traj_tpr,
                         "full-system TPR (topology + charges)")) return false;
        if (!RequireFile(spec, spec.traj_xtc,
                         "full-system XTC (protein + water + ions)")) return false;
        if (!spec.traj_edr.empty()) {
            if (!RequireFile(spec, spec.traj_edr,
                             "GROMACS .edr energy file")) return false;
        }
        return true;

    case JobMode::None:
        spec.error = "no mode specified";
        return false;
    }

    return false;  // unreachable
}


// ============================================================================
// Usage
// ============================================================================

void PrintJobSpecUsage(const char* prog) {
    fprintf(stderr,
        "Usage:\n"
        "  %s --pdb FILE [--pH N] [--config FILE] --output DIR\n"
        "      Load a bare PDB, protonate with reduce at the given pH\n"
        "      (default 7.0), assign ff14SB charges, run all calculators.\n"
        "\n"
        "  %s --protonated-pdb FILE [--config FILE] --output DIR\n"
        "      Load a PDB that already has hydrogen atoms (e.g. from reduce,\n"
        "      tleap, or GROMACS). Detect protonation state from H atoms.\n"
        "      Assign ff14SB charges.\n"
        "\n"
        "  %s --orca --root NAME [--config FILE] --output DIR\n"
        "      Load an ORCA DFT run. Root expands to:\n"
        "        {root}.xyz        coordinates from tleap (required)\n"
        "        {root}.prmtop     AMBER topology — charge source (required)\n"
        "        {root}_nmr.out    ORCA NMR shielding tensors (optional)\n"
        "\n"
        "  %s --mutant --wt NAME --ala NAME [--config FILE] --output DIR\n"
        "      Load a WT + ALA mutant pair. Each root expands as above.\n"
        "      Runs both conformations and computes WT-ALA delta tensors.\n"
        "\n"
        "  %s --trajectory --tpr FILE --xtc FILE [--edr FILE] --output DIR\n"
        "      Process a full-system GROMACS trajectory (protein + water + ions).\n"
        "      TPR: full-system topology and protein construction. XTC: full-system\n"
        "      trajectory. EDR: optional energy file.\n"
        "      Pass 1: scan all frames with lightweight calculators, accumulate stats.\n"
        "      Pass 2: extract selected frames with full calculators, write NPY.\n"
        "      Writes per-atom trajectory catalog (atom_catalog.csv).\n"
        "\n"
        "  %s --trajectory --analysis --tpr FILE --xtc FILE --output DIR\n"
        "      Analysis mode: single pass, all calculators on every sampled frame.\n"
        "      Writes exhaustive analysis H5 with full per-atom time series,\n"
        "      topology, and physics-group decomposition. ~1.2 GB per protein.\n"
        "\n"
        "Common options:\n"
        "  --output DIR     Output directory for NPY feature arrays (required for CLI)\n"
        "  --config FILE    TOML file with calculator parameter overrides\n"
        "  --no-mopac       Skip MOPAC semiempirical (and MopacCoulomb, MopacMcConnell)\n"
        "  --no-apbs        Skip APBS Poisson-Boltzmann solvated fields\n"
        "  --no-coulomb     Skip vacuum Coulomb EFG (APBS is preferred for electrostatics)\n"
        "  --aimnet2 FILE   AIMNet2 .jpt model for neural network charges + EFG\n"
        "  --help, -h       Show this message\n",
        prog, prog, prog, prog, prog, prog);
}


}  // namespace nmr

//
// nmr_extract: command-line entry point for all use cases.
//
// Parses command line via JobSpec, validates all inputs upfront,
// then dispatches to the appropriate builder + OperationRunner path.
//
// Usage: nmr_extract --help
//

#include "JobSpec.h"
#include "PdbFileReader.h"
#include "OrcaRunLoader.h"
#include "GromacsEnsembleLoader.h"
#include "OperationRunner.h"
#include "ConformationResult.h"
#include "OperationLog.h"
#include "RuntimeEnvironment.h"
#include "CalculatorConfig.h"

#include <cstdio>
#include <filesystem>

namespace fs = std::filesystem;
using namespace nmr;


// ============================================================================
// Use case A: bare PDB → protonate → run → write
// ============================================================================

static int RunPdb(const JobSpec& spec) {
    OperationLog::Info(LogFileIO, "nmr_extract",
        "PDB mode: " + spec.pdb_path + " pH=" + std::to_string(spec.pH));

    auto build = BuildFromPdb(spec.pdb_path, spec.pH);
    if (!build) {
        fprintf(stderr, "ERROR: %s\n", build.error.c_str());
        OperationLog::Error("nmr_extract", build.error);
        return 1;
    }

    auto& conf = build.protein->Conformation();

    RunOptions opts;
    opts.charge_source = build.charges.get();
    opts.net_charge = build.net_charge;

    auto result = OperationRunner::Run(conf, opts);
    if (!result.Ok()) {
        fprintf(stderr, "ERROR: %s\n", result.error.c_str());
        OperationLog::Error("nmr_extract", result.error);
        return 1;
    }

    fs::create_directories(spec.output_dir);
    int arrays = ConformationResult::WriteAllFeatures(conf, spec.output_dir);
    fprintf(stderr, "Wrote %d arrays to %s\n", arrays, spec.output_dir.c_str());
    OperationLog::Info(LogFileIO, "nmr_extract",
        "wrote " + std::to_string(arrays) + " arrays to " + spec.output_dir);
    return 0;
}


// ============================================================================
// Use case A': pre-protonated PDB → run → write
// ============================================================================

static int RunProtonatedPdb(const JobSpec& spec) {
    OperationLog::Info(LogFileIO, "nmr_extract",
        "protonated-pdb mode: " + spec.pdb_path);

    auto build = BuildFromProtonatedPdb(spec.pdb_path);
    if (!build) {
        fprintf(stderr, "ERROR: %s\n", build.error.c_str());
        OperationLog::Error("nmr_extract", build.error);
        return 1;
    }

    auto& conf = build.protein->Conformation();

    RunOptions opts;
    opts.charge_source = build.charges.get();
    opts.net_charge = build.net_charge;

    auto result = OperationRunner::Run(conf, opts);
    if (!result.Ok()) {
        fprintf(stderr, "ERROR: %s\n", result.error.c_str());
        OperationLog::Error("nmr_extract", result.error);
        return 1;
    }

    fs::create_directories(spec.output_dir);
    int arrays = ConformationResult::WriteAllFeatures(conf, spec.output_dir);
    fprintf(stderr, "Wrote %d arrays to %s\n", arrays, spec.output_dir.c_str());
    OperationLog::Info(LogFileIO, "nmr_extract",
        "wrote " + std::to_string(arrays) + " arrays to " + spec.output_dir);
    return 0;
}


// ============================================================================
// Use case B: ORCA DFT → run → write
// ============================================================================

static int RunOrca(const JobSpec& spec) {
    OperationLog::Info(LogFileIO, "nmr_extract",
        "orca mode: xyz=" + spec.orca_files.xyz_path);

    auto build = BuildFromOrca(spec.orca_files);
    if (!build) {
        fprintf(stderr, "ERROR: %s\n", build.error.c_str());
        OperationLog::Error("nmr_extract", build.error);
        return 1;
    }

    auto& conf = build.protein->Conformation();

    RunOptions opts;
    opts.charge_source = build.charges.get();
    opts.net_charge = build.net_charge;
    if (!spec.orca_files.nmr_out_path.empty())
        opts.orca_nmr_path = spec.orca_files.nmr_out_path;

    auto result = OperationRunner::Run(conf, opts);
    if (!result.Ok()) {
        fprintf(stderr, "ERROR: %s\n", result.error.c_str());
        OperationLog::Error("nmr_extract", result.error);
        return 1;
    }

    fs::create_directories(spec.output_dir);
    int arrays = ConformationResult::WriteAllFeatures(conf, spec.output_dir);
    fprintf(stderr, "Wrote %d arrays to %s\n", arrays, spec.output_dir.c_str());
    OperationLog::Info(LogFileIO, "nmr_extract",
        "wrote " + std::to_string(arrays) + " arrays to " + spec.output_dir);
    return 0;
}


// ============================================================================
// Use case C: WT + ALA mutant pair → run both → delta → write
// ============================================================================

static int RunMutant(const JobSpec& spec) {
    OperationLog::Info(LogFileIO, "nmr_extract",
        "mutant mode: wt=" + spec.wt_files.xyz_path +
        " ala=" + spec.ala_files.xyz_path);

    auto wt_build  = BuildFromOrca(spec.wt_files);
    auto ala_build = BuildFromOrca(spec.ala_files);
    if (!wt_build)  { fprintf(stderr, "WT ERROR: %s\n",  wt_build.error.c_str());  return 1; }
    if (!ala_build) { fprintf(stderr, "ALA ERROR: %s\n", ala_build.error.c_str()); return 1; }

    auto& wt_conf  = wt_build.protein->Conformation();
    auto& ala_conf = ala_build.protein->Conformation();

    RunOptions wt_opts;
    wt_opts.charge_source = wt_build.charges.get();
    wt_opts.net_charge = wt_build.net_charge;
    if (!spec.wt_files.nmr_out_path.empty())
        wt_opts.orca_nmr_path = spec.wt_files.nmr_out_path;

    RunOptions ala_opts;
    ala_opts.charge_source = ala_build.charges.get();
    ala_opts.net_charge = ala_build.net_charge;
    if (!spec.ala_files.nmr_out_path.empty())
        ala_opts.orca_nmr_path = spec.ala_files.nmr_out_path;

    auto result = OperationRunner::RunMutantComparison(
        wt_conf, wt_opts, ala_conf, ala_opts);
    if (!result.Ok()) {
        fprintf(stderr, "ERROR: %s\n", result.error.c_str());
        OperationLog::Error("nmr_extract", result.error);
        return 1;
    }

    fs::create_directories(spec.output_dir);
    int arrays = ConformationResult::WriteAllFeatures(wt_conf, spec.output_dir);
    fprintf(stderr, "Wrote %d arrays to %s\n", arrays, spec.output_dir.c_str());
    OperationLog::Info(LogFileIO, "nmr_extract",
        "wrote " + std::to_string(arrays) + " arrays to " + spec.output_dir);
    return 0;
}


// ============================================================================
// Use case D: GROMACS fleet → run all frames → write
// ============================================================================

static int RunFleet(const JobSpec& spec) {
    OperationLog::Info(LogFileIO, "nmr_extract",
        "fleet mode: poses=" + spec.fleet_paths.sampled_poses_dir +
        " tpr=" + spec.fleet_paths.tpr_path);

    auto build = BuildFromGromacs(spec.fleet_paths);
    if (!build) {
        fprintf(stderr, "ERROR: %s\n", build.error.c_str());
        OperationLog::Error("nmr_extract", build.error);
        return 1;
    }

    RunOptions opts;
    opts.charge_source = build.charges.get();
    opts.net_charge = build.net_charge;

    auto results = OperationRunner::RunEnsemble(*build.protein, opts);

    fs::create_directories(spec.output_dir);
    int total_arrays = 0;
    for (size_t i = 0; i < build.protein->ConformationCount(); ++i) {
        auto& conf = build.protein->ConformationAt(i);
        std::string subdir;
        if (i < build.pose_names.size() && !build.pose_names[i].empty()) {
            // Carry the source PDB name through, append our frame tag
            char tag[32];
            std::snprintf(tag, sizeof(tag), "_frame%03zu", i + 1);
            subdir = spec.output_dir + "/" + build.pose_names[i] + tag;
        } else {
            char frame_dir[512];
            std::snprintf(frame_dir, sizeof(frame_dir), "%s/frame_%03zu",
                          spec.output_dir.c_str(), i + 1);
            subdir = frame_dir;
        }
        fs::create_directories(subdir);
        total_arrays += ConformationResult::WriteAllFeatures(conf, subdir);
    }

    fprintf(stderr, "Wrote %d arrays across %zu frames to %s\n",
            total_arrays, build.protein->ConformationCount(),
            spec.output_dir.c_str());
    OperationLog::Info(LogFileIO, "nmr_extract",
        "wrote " + std::to_string(total_arrays) + " arrays across " +
        std::to_string(build.protein->ConformationCount()) + " frames to " +
        spec.output_dir);
    return 0;
}


// ============================================================================
// main
// ============================================================================

int main(int argc, char* argv[]) {
    RuntimeEnvironment::Load();

    auto spec = ParseJobSpec(argc, argv);

    if (spec.mode == JobMode::None) {
        if (!spec.error.empty())
            fprintf(stderr, "ERROR: %s\n\n", spec.error.c_str());
        PrintJobSpecUsage(argv[0]);
        return spec.error.empty() ? 0 : 1;
    }

    if (!ValidateJobSpec(spec)) {
        fprintf(stderr, "ERROR: %s\n", spec.error.c_str());
        return 1;
    }

    // Print warnings (e.g., missing optional NMR .out)
    for (const auto& w : spec.warnings)
        fprintf(stderr, "WARNING: %s\n", w.c_str());

    // Require --output for CLI
    if (spec.output_dir.empty()) {
        fprintf(stderr, "ERROR: --output DIR required\n");
        PrintJobSpecUsage(argv[0]);
        return 1;
    }

    // Load calculator config if specified
    if (!spec.config_path.empty())
        CalculatorConfig::Load(spec.config_path);

    switch (spec.mode) {
        case JobMode::Pdb:           return RunPdb(spec);
        case JobMode::ProtonatedPdb: return RunProtonatedPdb(spec);
        case JobMode::Orca:          return RunOrca(spec);
        case JobMode::Mutant:        return RunMutant(spec);
        case JobMode::Fleet:         return RunFleet(spec);
        case JobMode::None:          return 1;  // unreachable
    }
    return 1;
}

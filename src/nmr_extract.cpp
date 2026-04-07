//
// nmr_extract: command-line entry point for all use cases.
//
// Usage:
//   nmr_extract --pdb FILE [--pH N] --output DIR
//   nmr_extract --orca --xyz FILE --prmtop FILE [--nmr FILE] --output DIR
//   nmr_extract --mutant --wt-xyz FILE --wt-prmtop FILE [--wt-nmr FILE]
//                        --ala-xyz FILE --ala-prmtop FILE [--ala-nmr FILE]
//                        --output DIR
//   nmr_extract --fleet POSES_DIR TPR_PATH --output DIR
//
// All file paths are fully qualified. No directory scanning.
//

#include "PdbFileReader.h"
#include "OrcaRunLoader.h"
#include "GromacsEnsembleLoader.h"
#include "OperationRunner.h"
#include "ConformationResult.h"
#include "OperationLog.h"
#include "RuntimeEnvironment.h"

#include <cstdio>
#include <cstring>
#include <string>
#include <filesystem>

namespace fs = std::filesystem;
using namespace nmr;


static void PrintUsage(const char* prog) {
    fprintf(stderr,
        "Usage:\n"
        "  %s --pdb FILE [--pH N] --output DIR\n"
        "  %s --orca --xyz FILE --prmtop FILE [--nmr FILE] --output DIR\n"
        "  %s --mutant --wt-xyz FILE --wt-prmtop FILE [--wt-nmr FILE]\n"
        "              --ala-xyz FILE --ala-prmtop FILE [--ala-nmr FILE]\n"
        "              --output DIR\n"
        "  %s --fleet POSES_DIR TPR_PATH --output DIR\n",
        prog, prog, prog, prog);
}


// Helper: parse a named argument from argv. Returns "" if not found.
static std::string GetArg(int argc, char* argv[], const char* name) {
    for (int i = 1; i < argc - 1; ++i) {
        if (std::strcmp(argv[i], name) == 0)
            return argv[i + 1];
    }
    return "";
}


// ============================================================================
// Use case A: bare PDB → protonate → run → write
// ============================================================================

static int RunPdb(const std::string& pdb_path, double pH,
                  const std::string& output_dir) {
    fprintf(stderr, "Building from PDB: %s (pH=%.1f)\n", pdb_path.c_str(), pH);

    auto build = BuildFromPdb(pdb_path, pH);
    if (!build) {
        fprintf(stderr, "ERROR: %s\n", build.error.c_str());
        return 1;
    }

    auto& conf = build.protein->Conformation();

    RunOptions opts;
    opts.charge_source = build.charges.get();
    opts.net_charge = build.net_charge;

    auto result = OperationRunner::Run(conf, opts);
    if (!result.Ok()) {
        fprintf(stderr, "ERROR: %s\n", result.error.c_str());
        return 1;
    }

    fs::create_directories(output_dir);
    int arrays = ConformationResult::WriteAllFeatures(conf, output_dir);
    fprintf(stderr, "Wrote %d arrays to %s\n", arrays, output_dir.c_str());
    return 0;
}


// ============================================================================
// Use case B: ORCA DFT → run → write
// ============================================================================

static int RunOrca(const OrcaRunFiles& files,
                   const std::string& output_dir) {
    fprintf(stderr, "Building from ORCA: xyz=%s prmtop=%s\n",
            files.xyz_path.c_str(), files.prmtop_path.c_str());

    auto build = BuildFromOrca(files);
    if (!build) {
        fprintf(stderr, "ERROR: %s\n", build.error.c_str());
        return 1;
    }

    auto& conf = build.protein->Conformation();

    RunOptions opts;
    opts.charge_source = build.charges.get();
    opts.net_charge = build.net_charge;
    if (!files.nmr_out_path.empty())
        opts.orca_nmr_path = files.nmr_out_path;

    auto result = OperationRunner::Run(conf, opts);
    if (!result.Ok()) {
        fprintf(stderr, "ERROR: %s\n", result.error.c_str());
        return 1;
    }

    fs::create_directories(output_dir);
    int arrays = ConformationResult::WriteAllFeatures(conf, output_dir);
    fprintf(stderr, "Wrote %d arrays to %s\n", arrays, output_dir.c_str());
    return 0;
}


// ============================================================================
// Use case C: WT + ALA mutant pair → run both → delta → write
// ============================================================================

static int RunMutant(const OrcaRunFiles& wt_files, const OrcaRunFiles& ala_files,
                     const std::string& output_dir) {
    fprintf(stderr, "Building mutant pair:\n  WT:  xyz=%s prmtop=%s\n  ALA: xyz=%s prmtop=%s\n",
            wt_files.xyz_path.c_str(), wt_files.prmtop_path.c_str(),
            ala_files.xyz_path.c_str(), ala_files.prmtop_path.c_str());

    auto wt_build = BuildFromOrca(wt_files);
    auto ala_build = BuildFromOrca(ala_files);
    if (!wt_build)  { fprintf(stderr, "WT ERROR: %s\n", wt_build.error.c_str()); return 1; }
    if (!ala_build) { fprintf(stderr, "ALA ERROR: %s\n", ala_build.error.c_str()); return 1; }

    auto& wt_conf = wt_build.protein->Conformation();
    auto& ala_conf = ala_build.protein->Conformation();

    RunOptions wt_opts;
    wt_opts.charge_source = wt_build.charges.get();
    wt_opts.net_charge = wt_build.net_charge;
    if (!wt_files.nmr_out_path.empty())
        wt_opts.orca_nmr_path = wt_files.nmr_out_path;

    RunOptions ala_opts;
    ala_opts.charge_source = ala_build.charges.get();
    ala_opts.net_charge = ala_build.net_charge;
    if (!ala_files.nmr_out_path.empty())
        ala_opts.orca_nmr_path = ala_files.nmr_out_path;

    auto result = OperationRunner::RunMutantComparison(
        wt_conf, wt_opts, ala_conf, ala_opts);
    if (!result.Ok()) {
        fprintf(stderr, "ERROR: %s\n", result.error.c_str());
        return 1;
    }

    fs::create_directories(output_dir);
    int arrays = ConformationResult::WriteAllFeatures(wt_conf, output_dir);
    fprintf(stderr, "Wrote %d arrays to %s\n", arrays, output_dir.c_str());
    return 0;
}


// ============================================================================
// Use case D (trajectory): GROMACS fleet → run all frames → write
// ============================================================================

static int RunFleet(const std::string& poses_dir, const std::string& tpr_path,
                    const std::string& output_dir) {
    fprintf(stderr, "Building from GROMACS: %s\n", poses_dir.c_str());

    FleetPaths paths;
    paths.sampled_poses_dir = poses_dir;
    paths.tpr_path = tpr_path;

    auto build = BuildFromGromacs(paths);
    if (!build) {
        fprintf(stderr, "ERROR: %s\n", build.error.c_str());
        return 1;
    }

    RunOptions opts;
    opts.charge_source = build.charges.get();
    opts.net_charge = build.net_charge;

    auto results = OperationRunner::RunEnsemble(*build.protein, opts);

    fs::create_directories(output_dir);
    int total_arrays = 0;
    for (size_t i = 0; i < build.protein->ConformationCount(); ++i) {
        auto& conf = build.protein->ConformationAt(i);
        char frame_dir[512];
        std::snprintf(frame_dir, sizeof(frame_dir), "%s/frame_%03zu",
                      output_dir.c_str(), i + 1);
        fs::create_directories(frame_dir);
        total_arrays += ConformationResult::WriteAllFeatures(conf, frame_dir);
    }

    fprintf(stderr, "Wrote %d arrays across %zu frames to %s\n",
            total_arrays, build.protein->ConformationCount(),
            output_dir.c_str());
    return 0;
}


// ============================================================================
// main
// ============================================================================

int main(int argc, char* argv[]) {
    RuntimeEnvironment::Load();

    if (argc < 4) {
        PrintUsage(argv[0]);
        return 1;
    }

    std::string mode = argv[1];
    std::string output_dir = GetArg(argc, argv, "--output");

    if (output_dir.empty()) {
        fprintf(stderr, "ERROR: --output DIR required\n");
        PrintUsage(argv[0]);
        return 1;
    }

    if (mode == "--pdb") {
        std::string pH_str = GetArg(argc, argv, "--pH");
        double pH = pH_str.empty() ? 7.0 : std::atof(pH_str.c_str());
        return RunPdb(argv[2], pH, output_dir);
    }

    if (mode == "--orca") {
        OrcaRunFiles files;
        files.xyz_path = GetArg(argc, argv, "--xyz");
        files.prmtop_path = GetArg(argc, argv, "--prmtop");
        files.nmr_out_path = GetArg(argc, argv, "--nmr");
        if (files.xyz_path.empty()) {
            fprintf(stderr, "ERROR: --xyz FILE required for --orca\n");
            return 1;
        }
        return RunOrca(files, output_dir);
    }

    if (mode == "--mutant") {
        OrcaRunFiles wt_files;
        wt_files.xyz_path = GetArg(argc, argv, "--wt-xyz");
        wt_files.prmtop_path = GetArg(argc, argv, "--wt-prmtop");
        wt_files.nmr_out_path = GetArg(argc, argv, "--wt-nmr");

        OrcaRunFiles ala_files;
        ala_files.xyz_path = GetArg(argc, argv, "--ala-xyz");
        ala_files.prmtop_path = GetArg(argc, argv, "--ala-prmtop");
        ala_files.nmr_out_path = GetArg(argc, argv, "--ala-nmr");

        if (wt_files.xyz_path.empty() || ala_files.xyz_path.empty()) {
            fprintf(stderr, "ERROR: --wt-xyz and --ala-xyz required for --mutant\n");
            return 1;
        }
        return RunMutant(wt_files, ala_files, output_dir);
    }

    if (mode == "--fleet" && argc >= 5) {
        return RunFleet(argv[2], argv[3], output_dir);
    }

    fprintf(stderr, "ERROR: unknown mode '%s'\n", mode.c_str());
    PrintUsage(argv[0]);
    return 1;
}

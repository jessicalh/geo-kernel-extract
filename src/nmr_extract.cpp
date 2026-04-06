//
// nmr_extract: command-line entry point for all use cases.
//
// Usage:
//   nmr_extract --pdb FILE [--pH N] --output DIR
//   nmr_extract --orca DIR --output DIR
//   nmr_extract --mutant WT_DIR ALA_DIR --output DIR
//   nmr_extract --fleet POSES_DIR TPR_PATH --output DIR
//
// See spec/USE_CASES.md for the full description.
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
        "  %s --orca DIR --output DIR\n"
        "  %s --mutant WT_DIR ALA_DIR --output DIR\n"
        "  %s --fleet POSES_DIR TPR_PATH --output DIR\n",
        prog, prog, prog, prog);
}


// Find ORCA run files in a directory by convention:
//   *.pdb, *.xyz, *.prmtop, *_nmr.out
static OrcaRunFiles FindOrcaFiles(const std::string& dir) {
    OrcaRunFiles files;
    for (const auto& entry : fs::directory_iterator(dir)) {
        std::string name = entry.path().filename().string();
        std::string ext = entry.path().extension().string();
        if (ext == ".pdb")    files.pdb_path = entry.path().string();
        if (ext == ".xyz")    files.xyz_path = entry.path().string();
        if (ext == ".prmtop") files.prmtop_path = entry.path().string();
        if (name.find("_nmr.out") != std::string::npos)
            files.nmr_out_path = entry.path().string();
    }
    return files;
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

static int RunOrca(const std::string& orca_dir,
                   const std::string& output_dir) {
    auto files = FindOrcaFiles(orca_dir);
    if (files.xyz_path.empty()) {
        fprintf(stderr, "ERROR: no .xyz file in %s\n", orca_dir.c_str());
        return 1;
    }
    fprintf(stderr, "Building from ORCA: %s\n", orca_dir.c_str());

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

static int RunMutant(const std::string& wt_dir, const std::string& ala_dir,
                     const std::string& output_dir) {
    auto wt_files = FindOrcaFiles(wt_dir);
    auto ala_files = FindOrcaFiles(ala_dir);
    if (wt_files.xyz_path.empty() || ala_files.xyz_path.empty()) {
        fprintf(stderr, "ERROR: missing .xyz in WT or ALA directory\n");
        return 1;
    }
    fprintf(stderr, "Building mutant pair: WT=%s ALA=%s\n",
            wt_dir.c_str(), ala_dir.c_str());

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
    std::string output_dir;
    double pH = 7.0;

    // Parse --output from end of args
    for (int i = 1; i < argc - 1; ++i) {
        if (std::strcmp(argv[i], "--output") == 0) {
            output_dir = argv[i + 1];
        }
        if (std::strcmp(argv[i], "--pH") == 0) {
            pH = std::atof(argv[i + 1]);
        }
    }

    if (output_dir.empty()) {
        fprintf(stderr, "ERROR: --output DIR required\n");
        PrintUsage(argv[0]);
        return 1;
    }

    if (mode == "--pdb") {
        return RunPdb(argv[2], pH, output_dir);
    }
    if (mode == "--orca") {
        return RunOrca(argv[2], output_dir);
    }
    if (mode == "--mutant" && argc >= 5) {
        return RunMutant(argv[2], argv[3], output_dir);
    }
    if (mode == "--fleet" && argc >= 5) {
        return RunFleet(argv[2], argv[3], output_dir);
    }

    fprintf(stderr, "ERROR: unknown mode '%s'\n", mode.c_str());
    PrintUsage(argv[0]);
    return 1;
}

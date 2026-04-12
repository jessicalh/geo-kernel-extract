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
#include "GromacsProtein.h"
#include "OperationRunner.h"
#include "ConformationResult.h"
#include "AIMNet2Result.h"
#include "OperationLog.h"
#include "RuntimeEnvironment.h"
#include "CalculatorConfig.h"
#include "xtc_reader.h"

#include <cstdio>
#include <filesystem>

namespace fs = std::filesystem;
using namespace nmr;

// AIMNet2 model: loaded once in main, shared across all use cases.
static std::unique_ptr<AIMNet2Model> g_aimnet2_model;


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
    opts.skip_mopac   = spec.skip_mopac;
    opts.skip_apbs    = spec.skip_apbs;
    opts.skip_coulomb = spec.skip_coulomb;
    opts.aimnet2_model = g_aimnet2_model.get();

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
    opts.skip_mopac   = spec.skip_mopac;
    opts.skip_apbs    = spec.skip_apbs;
    opts.skip_coulomb = spec.skip_coulomb;
    opts.aimnet2_model = g_aimnet2_model.get();

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
    opts.skip_mopac   = spec.skip_mopac;
    opts.skip_apbs    = spec.skip_apbs;
    opts.skip_coulomb = spec.skip_coulomb;
    opts.aimnet2_model = g_aimnet2_model.get();
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
    wt_opts.skip_mopac   = spec.skip_mopac;
    wt_opts.skip_apbs    = spec.skip_apbs;
    wt_opts.skip_coulomb = spec.skip_coulomb;
    wt_opts.aimnet2_model = g_aimnet2_model.get();
    if (!spec.wt_files.nmr_out_path.empty())
        wt_opts.orca_nmr_path = spec.wt_files.nmr_out_path;

    RunOptions ala_opts;
    ala_opts.charge_source = ala_build.charges.get();
    ala_opts.net_charge = ala_build.net_charge;
    ala_opts.skip_mopac   = spec.skip_mopac;
    ala_opts.skip_apbs    = spec.skip_apbs;
    ala_opts.skip_coulomb = spec.skip_coulomb;
    ala_opts.aimnet2_model = g_aimnet2_model.get();
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
    opts.skip_mopac   = spec.skip_mopac;
    opts.skip_apbs    = spec.skip_apbs;
    opts.skip_coulomb = spec.skip_coulomb;
    opts.aimnet2_model = g_aimnet2_model.get();

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
// Use case E: full-system GROMACS trajectory
//
// Pass 1 — scan: lightweight calculators on all frames, accumulate
//          per-atom Welford stats.  No NPY.
// Finalize — write atom_catalog.csv, select frames.
// Pass 2 — extract: full calculators on selected frames, write NPY.
//
// Process control stays here. Data + computation on GromacsProtein.
// ============================================================================

static int RunTrajectory(const JobSpec& spec) {
    OperationLog::Info(LogFileIO, "nmr_extract",
        "trajectory mode: tpr=" + spec.traj_tpr +
        " xtc=" + spec.traj_xtc + " ref=" + spec.traj_ref);

    // Build trajectory manager
    GromacsProtein gp;
    if (!gp.BuildFromTrajectory(spec.traj_tpr, spec.traj_ref)) {
        fprintf(stderr, "ERROR: %s\n", gp.error().c_str());
        return 1;
    }

    // Load trajectory
    if (!gp.LoadTrajectory(spec.traj_xtc)) {
        fprintf(stderr, "ERROR: %s\n", gp.error().c_str());
        return 1;
    }
    fprintf(stderr, "Loaded %zu frames, %zu protein atoms\n",
            gp.FrameCount(), gp.AtomCount());

    // ── Pass 1: scan ─────────────────────────────────────────────

    RunOptions scan_opts;
    scan_opts.charge_source = gp.charges();
    scan_opts.skip_mopac   = true;
    scan_opts.skip_apbs    = true;
    scan_opts.skip_coulomb = true;
    scan_opts.skip_dssp    = true;

    fprintf(stderr, "\n=== Pass 1: scan %zu frames ===\n", gp.FrameCount());
    for (size_t fi = 0; fi < gp.FrameCount(); ++fi) {
        if (fi % 100 == 0)
            fprintf(stderr, "  scan %zu / %zu\n", fi, gp.FrameCount());
        gp.ScanFrame(fi, scan_opts);
    }

    // ── Finalize: catalog + frame selection ───────────────────────

    fs::create_directories(spec.output_dir);
    std::string catalog_path = spec.output_dir + "/atom_catalog.csv";
    gp.WriteCatalog(catalog_path);
    fprintf(stderr, "\nWrote %s\n", catalog_path.c_str());

    auto selected = gp.SelectFrames(200);
    fprintf(stderr, "Selected %zu frames for extraction\n", selected.size());

    if (selected.empty()) {
        fprintf(stderr, "Scan complete — no frames selected.\n");
        return 0;
    }

    // ── Pass 2: extract selected frames ──────────────────────────

    RunOptions extract_opts;
    extract_opts.charge_source = gp.charges();
    extract_opts.skip_mopac   = spec.skip_mopac;
    extract_opts.skip_apbs    = spec.skip_apbs;
    extract_opts.skip_coulomb = spec.skip_coulomb;
    extract_opts.aimnet2_model = g_aimnet2_model.get();
    if (!spec.traj_edr.empty())
        extract_opts.edr_path = spec.traj_edr;

    fprintf(stderr, "\n=== Pass 2: extract %zu frames ===\n", selected.size());
    int total_arrays = 0;
    for (size_t si = 0; si < selected.size(); ++si) {
        size_t fi = selected[si];
        fprintf(stderr, "  extract frame %zu (%zu/%zu)\n",
                fi, si + 1, selected.size());
        int n = gp.ExtractFrame(fi, extract_opts, spec.output_dir);
        if (n >= 0) total_arrays += n;
    }

    fprintf(stderr, "\nDone: %d arrays, %zu frames, %s\n",
            total_arrays, selected.size(), spec.output_dir.c_str());
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

    // Load AIMNet2 model if specified (once, shared across all conformations)
    if (!spec.aimnet2_model_path.empty()) {
        g_aimnet2_model = AIMNet2Model::Load(spec.aimnet2_model_path);
        if (!g_aimnet2_model) {
            fprintf(stderr, "ERROR: Failed to load AIMNet2 model: %s\n",
                    spec.aimnet2_model_path.c_str());
            return 1;
        }
    }

    switch (spec.mode) {
        case JobMode::Pdb:           return RunPdb(spec);
        case JobMode::ProtonatedPdb: return RunProtonatedPdb(spec);
        case JobMode::Orca:          return RunOrca(spec);
        case JobMode::Mutant:        return RunMutant(spec);
        case JobMode::Fleet:         return RunFleet(spec);
        case JobMode::Trajectory:    return RunTrajectory(spec);
        case JobMode::None:          return 1;  // unreachable
    }
    return 1;
}

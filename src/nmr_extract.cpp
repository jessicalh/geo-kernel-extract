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
#include "TrajectoryProtein.h"
#include "Trajectory.h"
#include "RunConfiguration.h"
#include "Session.h"
#include "errors.h"
#include "OperationRunner.h"
#include "ConformationResult.h"
#include "AIMNet2Result.h"
#include "OperationLog.h"
#include "RuntimeEnvironment.h"
#include "CalculatorConfig.h"
#include "GromacsFrameHandler.h"

#include <highfive/H5File.hpp>

#include <cstdio>
#include <filesystem>
#include <set>

namespace fs = std::filesystem;
using namespace nmr;

// Resources that every RunXxx needs live on `Session`, constructed
// in main and passed by const reference. No global singletons here.


// ============================================================================
// Use case A: bare PDB → protonate → run → write
// ============================================================================

static int RunPdb(const JobSpec& spec, const Session& session) {
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
    opts.aimnet2_model = session.Aimnet2Model();

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

static int RunProtonatedPdb(const JobSpec& spec, const Session& session) {
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
    opts.aimnet2_model = session.Aimnet2Model();

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

static int RunOrca(const JobSpec& spec, const Session& session) {
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
    opts.aimnet2_model = session.Aimnet2Model();
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

static int RunMutant(const JobSpec& spec, const Session& session) {
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
    wt_opts.aimnet2_model = session.Aimnet2Model();
    if (!spec.wt_files.nmr_out_path.empty())
        wt_opts.orca_nmr_path = spec.wt_files.nmr_out_path;

    RunOptions ala_opts;
    ala_opts.charge_source = ala_build.charges.get();
    ala_opts.net_charge = ala_build.net_charge;
    ala_opts.skip_mopac   = spec.skip_mopac;
    ala_opts.skip_apbs    = spec.skip_apbs;
    ala_opts.skip_coulomb = spec.skip_coulomb;
    ala_opts.aimnet2_model = session.Aimnet2Model();
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
// Use case D: full-system GROMACS trajectory
//
// TrajectoryProtein owns Protein + trajectory-scope model.
// Trajectory owns file paths + EDR and drives the run.
// RunConfiguration describes the run shape.
// ============================================================================

static int RunTrajectory(const JobSpec& spec, const Session& session) {
    OperationLog::Info(LogFileIO, "nmr_extract",
        "trajectory mode: dir=" + spec.traj_dir);

    // ── TrajectoryProtein: TPR → topology + protein + charges ────
    TrajectoryProtein tp;
    if (!tp.BuildFromTrajectory(spec.traj_dir)) {
        fprintf(stderr, "ERROR: %s\n", tp.Error().c_str());
        return 1;
    }

    // ── Trajectory: file paths + EDR preload ─────────────────────
    Trajectory traj(spec.traj_xtc, spec.traj_tpr,
                    fs::path(spec.traj_dir) / "md.edr");

    // ── Run shape: PerFrameExtractionSet; aimnet2 comes from session ─
    RunConfiguration config = RunConfiguration::PerFrameExtractionSet();

    // ── Drive ────────────────────────────────────────────────────
    const Status s = traj.Run(tp, config, session, /*extras=*/{},
                              /*output_dir=*/spec.output_dir);
    if (s != kOk) {
        fprintf(stderr, "ERROR: Trajectory::Run returned status 0x%x\n", s);
        return 1;
    }

    // ── Emit output H5 ───────────────────────────────────────────
    fs::create_directories(spec.output_dir);
    const std::string h5_path = spec.output_dir + "/trajectory.h5";
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        traj.WriteH5(file);
        tp.WriteH5(file);
    }
    fprintf(stderr, "Wrote %s (%zu frames, %zu atoms, %zu selections)\n",
            h5_path.c_str(), traj.FrameCount(),
            tp.AtomCount(), traj.Selections().Count());

    return 0;
}


// ============================================================================
// Analysis mode: --trajectory --analysis
//
// Stubbed pending dissolution of AnalysisWriter into a family of
// *TimeSeriesTrajectoryResult classes (each owning a DenseBuffer,
// each writing its own H5 group).
// ============================================================================

static int RunAnalysis(const JobSpec& spec) {
    (void)spec;
    fprintf(stderr,
        "ERROR: --trajectory --analysis is disabled pending the\n"
        "dissolution of AnalysisWriter into per-Result H5 emitters.\n"
        "Use --trajectory (no --analysis) for the Trajectory::Run path.\n");
    return 1;
}


// ============================================================================
// main
// ============================================================================

int main(int argc, char* argv[]) {
    // ── Session: one named object that holds process-wide resources ─
    // Loads RuntimeEnvironment, OperationLog channel config, emits
    // session-start log line. CalculatorConfig + AIMNet2 model are
    // loaded below after CLI parse (each depends on spec fields).
    Session session;
    if (session.LoadFromToml() != kOk) {
        fprintf(stderr, "ERROR: session load: %s\n",
                session.LastError().c_str());
        return 1;
    }

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

    // Load calculator config: explicit --config, or default from data dir.
    if (!spec.config_path.empty()) {
        CalculatorConfig::Load(spec.config_path);
    } else {
        std::string default_config =
            std::string(NMR_DATA_DIR) + "/calculator_params.toml";
        if (fs::exists(default_config))
            CalculatorConfig::Load(default_config);
    }

    // AIMNet2 model: CLI --aimnet2 takes priority; TOML fallback.
    // Session holds the loaded model for the rest of the process.
    if (spec.aimnet2_model_path.empty())
        spec.aimnet2_model_path =
            CalculatorConfig::GetString("aimnet2_model_path");
    if (!spec.aimnet2_model_path.empty()) {
        if (session.LoadAimnet2Model(spec.aimnet2_model_path) != kOk) {
            fprintf(stderr, "ERROR: %s\n", session.LastError().c_str());
            return 1;
        }
    }

    switch (spec.mode) {
        case JobMode::Pdb:           return RunPdb(spec, session);
        case JobMode::ProtonatedPdb: return RunProtonatedPdb(spec, session);
        case JobMode::Orca:          return RunOrca(spec, session);
        case JobMode::Mutant:        return RunMutant(spec, session);
        case JobMode::Trajectory:
            return spec.analysis ? RunAnalysis(spec)
                                 : RunTrajectory(spec, session);
        case JobMode::None:          return 1;  // unreachable
    }
    return 1;
}

/// @file
/// nmr_extract: command-line entry point for all supported modes.
///
/// Parses argv via @ref nmr::cli::Parse, validates inputs, then
/// dispatches to the typed mode runner via @c std::visit.

#include "AIMNet2Result.h"
#include "CalculatorConfig.h"
#include "CategoryInfoProjection.h"
#include "Cli/Parse.h"
#include "Cli/PrintUsage.h"
#include "Cli/Validate.h"
#include "ConformationResult.h"
#include "FrameNpyEmitter.h"
#include "FramePdbEmitter.h"
#include "GromacsFrameHandler.h"
#include "OperationLog.h"
#include "OperationRunner.h"
#include "OrcaRunLoader.h"
#include "PdbFileReader.h"
#include "RunConfiguration.h"
#include "RuntimeEnvironment.h"
#include "Session.h"
#include "TopologySidecar.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"
#include "errors.h"

#include <highfive/H5File.hpp>

#include <cstdio>
#include <exception>
#include <filesystem>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>

namespace fs = std::filesystem;
using namespace nmr;

namespace {

fs::path ResolveConfigRelativePath(const std::string& path,
                                   const fs::path& config_path) {
    fs::path p(path);
    if (p.empty() || p.is_absolute() || config_path.empty()) return p;
    fs::path base = config_path.parent_path();
    if (base.empty()) base = ".";
    return base / p;
}

/// Populate the common cross-mode @ref RunOptions slots. Mode-specific
/// fields (charges, orca_nmr_path) are filled by the caller.
RunOptions MakeBaseOpts(const Session& session, bool mopac) {
    RunOptions opts;
    opts.skip_mopac           = !mopac;
    opts.skip_apbs            = false;  // APBS is the canonical electrostatics — always on.
    opts.skip_coulomb         = true;   // Home-rolled vacuum Coulomb retired; APBS supersedes.
    opts.aimnet2_model        = session.Aimnet2Model();
    opts.tripeptide_dft_table = session.TripeptideDftTablePtr();
    opts.larsen_hbond_grid    = session.LarsenHBondGridPtr();
    return opts;
}

/// Emit the per-protein sidecars (CategoryInfoProjection + TopologySidecar)
/// to @p output_dir. Returns true on success, prints error and returns
/// false otherwise. Single-mode runners share this.
bool WriteSidecars(const Protein& protein, const std::string& output_dir) {
    const int cat  = CategoryInfoProjection::WriteFeatures(protein, output_dir);
    const int topo = TopologySidecar::WriteFeatures(protein, output_dir);
    if (cat != 1 || topo != 5) {
        std::fprintf(stderr,
            "ERROR: incomplete sidecar emission "
            "(atoms_category=%d/1, topology=%d/5) -- disk full or permission?\n",
            cat, topo);
        return false;
    }
    return true;
}

/// Emit WriteAllFeatures and report. Single-mode runners share this.
int WriteFeaturesAndReport(const ProteinConformation& conf, const std::string& output_dir) {
    const int arrays = ConformationResult::WriteAllFeatures(conf, output_dir);
    std::fprintf(stderr, "Wrote %d arrays to %s\n", arrays, output_dir.c_str());
    OperationLog::Info(LogFileIO, "nmr_extract",
        "wrote " + std::to_string(arrays) + " arrays to " + output_dir);
    return 0;
}

}  // namespace


// ============================================================================
// Mode runners
// ============================================================================

static int RunPdb(const cli::PdbMode& mode, const cli::CommonOptions& common,
                  const Session& session) {
    OperationLog::Info(LogFileIO, "nmr_extract",
        "PDB mode: " + mode.pdb.string() + " pH=" + std::to_string(mode.pH));

    auto build = BuildFromPdb(mode.pdb.string(), mode.pH);
    if (!build) {
        std::fprintf(stderr, "ERROR: %s\n", build.error.c_str());
        OperationLog::Error("nmr_extract", build.error);
        return 1;
    }

    auto& conf = build.protein->Conformation();
    RunOptions opts        = MakeBaseOpts(session, mode.mopac);
    opts.charge_source     = build.charges.get();
    opts.net_charge        = build.net_charge;

    auto result = OperationRunner::Run(conf, opts);
    if (!result.Ok()) {
        std::fprintf(stderr, "ERROR: %s\n", result.error.c_str());
        OperationLog::Error("nmr_extract", result.error);
        return 1;
    }

    const std::string output_dir = common.output_dir.string();
    fs::create_directories(output_dir);
    if (!WriteSidecars(*build.protein, output_dir)) return 1;
    return WriteFeaturesAndReport(conf, output_dir);
}


static int RunProtonatedPdb(const cli::ProtonatedPdbMode& mode,
                            const cli::CommonOptions& common,
                            const Session& session) {
    OperationLog::Info(LogFileIO, "nmr_extract",
        "protonated-pdb mode: " + mode.pdb.string());

    auto build = BuildFromProtonatedPdb(mode.pdb.string());
    if (!build) {
        std::fprintf(stderr, "ERROR: %s\n", build.error.c_str());
        OperationLog::Error("nmr_extract", build.error);
        return 1;
    }

    auto& conf = build.protein->Conformation();
    RunOptions opts    = MakeBaseOpts(session, mode.mopac);
    opts.charge_source = build.charges.get();
    opts.net_charge    = build.net_charge;

    auto result = OperationRunner::Run(conf, opts);
    if (!result.Ok()) {
        std::fprintf(stderr, "ERROR: %s\n", result.error.c_str());
        OperationLog::Error("nmr_extract", result.error);
        return 1;
    }

    const std::string output_dir = common.output_dir.string();
    fs::create_directories(output_dir);
    if (!WriteSidecars(*build.protein, output_dir)) return 1;
    return WriteFeaturesAndReport(conf, output_dir);
}


static int RunOrca(const cli::OrcaMode& mode, const cli::CommonOptions& common,
                   const Session& session) {
    OperationLog::Info(LogFileIO, "nmr_extract",
        "orca mode: xyz=" + mode.files.xyz_path);

    auto build = BuildFromOrca(mode.files);
    if (!build) {
        std::fprintf(stderr, "ERROR: %s\n", build.error.c_str());
        OperationLog::Error("nmr_extract", build.error);
        return 1;
    }

    auto& conf = build.protein->Conformation();
    RunOptions opts    = MakeBaseOpts(session, mode.mopac);
    opts.charge_source = build.charges.get();
    opts.net_charge    = build.net_charge;
    if (!mode.files.nmr_out_path.empty()) opts.orca_nmr_path = mode.files.nmr_out_path;

    auto result = OperationRunner::Run(conf, opts);
    if (!result.Ok()) {
        std::fprintf(stderr, "ERROR: %s\n", result.error.c_str());
        OperationLog::Error("nmr_extract", result.error);
        return 1;
    }

    const std::string output_dir = common.output_dir.string();
    fs::create_directories(output_dir);
    if (!WriteSidecars(*build.protein, output_dir)) return 1;
    return WriteFeaturesAndReport(conf, output_dir);
}


static int RunMutant(const cli::MutantMode& mode, const cli::CommonOptions& common,
                     const Session& session) {
    OperationLog::Info(LogFileIO, "nmr_extract",
        "mutant mode: wt=" + mode.wt.xyz_path + " ala=" + mode.ala.xyz_path);

    auto wt_build  = BuildFromOrca(mode.wt);
    auto ala_build = BuildFromOrca(mode.ala);
    if (!wt_build)  { std::fprintf(stderr, "WT ERROR: %s\n",  wt_build.error.c_str());  return 1; }
    if (!ala_build) { std::fprintf(stderr, "ALA ERROR: %s\n", ala_build.error.c_str()); return 1; }

    auto& wt_conf  = wt_build.protein->Conformation();
    auto& ala_conf = ala_build.protein->Conformation();

    RunOptions wt_opts        = MakeBaseOpts(session, mode.mopac);
    wt_opts.charge_source     = wt_build.charges.get();
    wt_opts.net_charge        = wt_build.net_charge;
    if (!mode.wt.nmr_out_path.empty())  wt_opts.orca_nmr_path = mode.wt.nmr_out_path;

    RunOptions ala_opts       = MakeBaseOpts(session, mode.mopac);
    ala_opts.charge_source    = ala_build.charges.get();
    ala_opts.net_charge       = ala_build.net_charge;
    if (!mode.ala.nmr_out_path.empty()) ala_opts.orca_nmr_path = mode.ala.nmr_out_path;

    auto result = OperationRunner::RunMutantComparison(
        wt_conf, wt_opts, ala_conf, ala_opts);
    if (!result.Ok()) {
        std::fprintf(stderr, "ERROR: %s\n", result.error.c_str());
        OperationLog::Error("nmr_extract", result.error);
        return 1;
    }

    const std::string output_dir = common.output_dir.string();
    fs::create_directories(output_dir);
    if (!WriteSidecars(*wt_build.protein, output_dir)) return 1;
    return WriteFeaturesAndReport(wt_conf, output_dir);
}


static int RunTrajectory(const cli::TrajectoryMode& mode,
                         const cli::CommonOptions& common,
                         const Session& session) {
    const std::string traj_dir = mode.dir.string();
    OperationLog::Info(LogFileIO, "nmr_extract", "trajectory mode: dir=" + traj_dir);

    const std::string output_dir = common.output_dir.string();
    if (output_dir.empty()) {
        std::fprintf(stderr, "ERROR: --trajectory requires --output DIR\n");
        return 1;
    }

    TrajectoryProtein tp;
    if (!tp.BuildFromTrajectory(traj_dir)) {
        std::fprintf(stderr, "ERROR: %s\n", tp.Error().c_str());
        return 1;
    }

    // Per-frame sidecars always emit in trajectory mode — one set per
    // dispatched frame, cadence governed solely by --stride. PDBs to
    // output/pdbs, the full per-conformation NPY set to output/npys.
    {
        std::string stem = fs::path(traj_dir).filename().string();
        if (stem.empty()) {
            stem = fs::path(traj_dir).parent_path().filename().string();
        }
        FramePdbEmitter::Config pcfg;
        pcfg.output_dir = fs::path(output_dir) / "pdbs";
        pcfg.stem       = stem;
        FramePdbEmitter::Configure(tp.ProteinRef(), std::move(pcfg));

        FrameNpyEmitter::Config ncfg;
        ncfg.output_dir = fs::path(output_dir) / "npys";
        FrameNpyEmitter::Configure(tp.ProteinRef(), std::move(ncfg));
    }

    const auto files = cli::TrajectoryInputFiles::FromProductionDir(mode.dir);
    Trajectory traj(files.trr, files.tpr, files.edr);

    RunConfiguration config = mode.mopac
        ? RunConfiguration::FullFatFrameExtraction()
        : RunConfiguration::PerFrameExtractionSet();
    // The single cadence knob: process every mode.stride-th TRR frame.
    config.SetStride(mode.stride);
    if (mode.mopac) {
        config.MutablePerFrameRunOptions().net_charge = tp.NetCharge();
    }

    const Status s = traj.Run(tp, config, session, /*extras=*/{}, /*output_dir=*/output_dir);
    if (s != kOk) {
        std::fprintf(stderr, "ERROR: Trajectory::Run returned status 0x%x\n", s);
        return 1;
    }

    fs::create_directories(output_dir);
    const std::string h5_path = output_dir + "/trajectory.h5";
    {
        HighFive::File file(h5_path, HighFive::File::Truncate);
        traj.WriteH5(file);
        tp.WriteH5(file);
    }
    std::fprintf(stderr, "Wrote %s (%zu frames, %zu atoms, %zu selections)\n",
            h5_path.c_str(), traj.FrameCount(),
            tp.AtomCount(), traj.Selections().Count());
    return 0;
}


// ============================================================================
// main
// ============================================================================

// All mode handlers use return codes (errors.h), not exceptions — but
// parsing, config/model loading, std::filesystem and allocation can still
// throw. RunExtract holds the real body; main is the outermost boundary that
// turns any escaped exception into a fail-loud diagnostic + nonzero exit
// rather than std::terminate with no message.
static int RunExtract(int argc, char* argv[]) {
    Session session;
    if (session.LoadFromToml() != kOk) {
        std::fprintf(stderr, "ERROR: session load: %s\n", session.LastError().c_str());
        return 1;
    }

    auto parse = cli::Parse(argc, argv);

    if (parse.help_requested) {
        cli::PrintUsage(argv[0]);
        return 0;
    }
    if (!parse.error.empty()) {
        std::fprintf(stderr, "ERROR: %s\n\n", parse.error.c_str());
        cli::PrintUsage(argv[0]);
        return 1;
    }

    // ParseResult's contract (Parse.h): on success spec and common are set
    // and error is empty. We've already returned on help_requested and on a
    // non-empty error, so this is the success path. Guard the invariant
    // anyway -- a contract violation in Parse would otherwise dereference an
    // empty optional (UB); fail loud instead.
    if (!parse.spec || !parse.common) {
        std::fprintf(stderr,
            "internal: parser reported success but left spec/common unset\n");
        return 1;
    }
    cli::ModeSpec      spec   = std::move(*parse.spec);
    cli::CommonOptions common = std::move(*parse.common);

    // Load CalculatorConfig before AIMNet2 fallback so the TOML key
    // is available for resolution.
    fs::path loaded_config_path;
    if (!common.config_path.empty()) {
        loaded_config_path = common.config_path;
        CalculatorConfig::Load(common.config_path.string());
    } else {
        loaded_config_path = fs::path(NMR_DATA_DIR) / "calculator_params.toml";
        if (fs::exists(loaded_config_path)) {
            CalculatorConfig::Load(loaded_config_path.string());
        } else {
            loaded_config_path.clear();
        }
    }

    // AIMNet2 fallback: CLI flag wins, then NMR_AIMNET2_MODEL, then
    // calculator_params.toml aimnet2_model_path. TOML-relative model paths
    // resolve from the TOML file directory, so committed configs stay
    // relocatable across machines and containers.
    if (common.aimnet2_model_path.empty()) {
        const char* env_model = std::getenv("NMR_AIMNET2_MODEL");
        if (env_model && *env_model) common.aimnet2_model_path = env_model;
    }
    if (common.aimnet2_model_path.empty()) {
        const std::string toml_default = CalculatorConfig::GetString("aimnet2_model_path");
        if (!toml_default.empty()) {
            common.aimnet2_model_path =
                ResolveConfigRelativePath(toml_default, loaded_config_path);
        }
    }
    if (common.aimnet2_model_path.empty()) {
        std::fprintf(stderr,
            "ERROR: AIMNet2 model required: pass --aimnet2 PATH, set "
            "NMR_AIMNET2_MODEL, or set aimnet2_model_path in "
            "calculator_params.toml.\n");
        return 1;
    }

    if (const std::string verror = cli::Validate(spec, common); !verror.empty()) {
        std::fprintf(stderr, "ERROR: %s\n", verror.c_str());
        return 1;
    }

    for (const auto& w : parse.warnings) {
        std::fprintf(stderr, "WARNING: %s\n", w.c_str());
    }

    if (common.output_dir.empty()) {
        std::fprintf(stderr, "ERROR: --output DIR required\n");
        cli::PrintUsage(argv[0]);
        return 1;
    }

    if (session.LoadAimnet2Model(common.aimnet2_model_path.string()) != kOk) {
        std::fprintf(stderr, "ERROR: %s\n", session.LastError().c_str());
        return 1;
    }
    if (session.LoadTripeptideDftTable() != kOk) {
        std::fprintf(stderr, "ERROR: %s\n", session.LastError().c_str());
        return 1;
    }
    if (session.LoadLarsenHBondGrid() != kOk) {
        std::fprintf(stderr, "ERROR: %s\n", session.LastError().c_str());
        return 1;
    }

    return std::visit(
        [&](const auto& m) -> int {
            using T = std::decay_t<decltype(m)>;
            if constexpr (std::is_same_v<T, cli::PdbMode>)            return RunPdb(m, common, session);
            else if constexpr (std::is_same_v<T, cli::ProtonatedPdbMode>) return RunProtonatedPdb(m, common, session);
            else if constexpr (std::is_same_v<T, cli::OrcaMode>)          return RunOrca(m, common, session);
            else if constexpr (std::is_same_v<T, cli::MutantMode>)        return RunMutant(m, common, session);
            else if constexpr (std::is_same_v<T, cli::TrajectoryMode>)    return RunTrajectory(m, common, session);
            else {
                std::fprintf(stderr, "internal: unhandled mode variant\n");
                return 1;
            }
        },
        spec);
}

int main(int argc, char* argv[]) {
    try {
        return RunExtract(argc, argv);
    } catch (const std::exception& e) {
        // stderr first: e.what() is a const char*, so this diagnostic needs no
        // allocation and lands even under memory pressure. The structured log
        // call (whose arg build could allocate) follows, itself guarded so a
        // bad_alloc there cannot turn a clean exit-1 into std::terminate.
        std::fprintf(stderr, "FATAL: unhandled exception: %s\n", e.what());
        try {
            OperationLog::Error("nmr_extract::main",
                                std::string("unhandled exception: ") + e.what());
        } catch (...) { /* stderr already carried it */ }
        return 1;
    } catch (...) {
        std::fprintf(stderr, "FATAL: unhandled non-std exception\n");
        try {
            OperationLog::Error("nmr_extract::main", "unhandled non-std exception");
        } catch (...) { /* stderr already carried it */ }
        return 1;
    }
}

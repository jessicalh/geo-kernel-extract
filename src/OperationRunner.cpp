#include "OperationRunner.h"
#include "Protein.h"

#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "EnrichmentResult.h"
#include "DsspResult.h"
#include "ChargeAssignmentResult.h"
#include "MopacResult.h"
#include "MopacCoulombResult.h"
#include "MopacMcConnellResult.h"
#include "ApbsFieldResult.h"
#include "BiotSavartResult.h"
#include "HaighMallionResult.h"
#include "McConnellResult.h"
#include "RingSusceptibilityResult.h"
#include "PiQuadrupoleResult.h"
#include "DispersionResult.h"
#include "CoulombResult.h"
#include "HBondResult.h"
#include "OrcaShieldingResult.h"
#include "MutationDeltaResult.h"
#include "AIMNet2Result.h"
#include "SasaResult.h"
#include "GromacsEnergyResult.h"
#include "BondedEnergyResult.h"
#include "WaterFieldResult.h"
#include "HydrationShellResult.h"
#include "HydrationGeometryResult.h"
#include "EeqResult.h"
#include "OperationLog.h"

namespace nmr {


// Helper: attach a result, log failure, accumulate name list.
static bool Attach(ProteinConformation& conf,
                   std::unique_ptr<ConformationResult> result,
                   const char* name,
                   RunResult& out) {
    if (!result) {
        out.error = std::string(name) + " computation returned null";
        return false;
    }
    if (!conf.AttachResult(std::move(result))) {
        out.error = std::string(name) + " failed to attach";
        return false;
    }
    out.attached.push_back(name);
    return true;
}


// =================================================================
// Run: the standard sequence.
//
// Order:
//   Tier 0 — foundation (no inter-dependencies)
//   Tier 0.5 — external tools (need charges)
//   Tier 1 — 8 classical calculators
//   Tier 2 — DFT comparison (optional)
// =================================================================

RunResult OperationRunner::Run(ProteinConformation& conf,
                               const RunOptions& opts) {
    RunResult out;

    OperationLog::Scope scope("OperationRunner::Run",
        "atoms=" + std::to_string(conf.AtomCount()));

    // --- Tier 0: foundation ---

    if (!Attach(conf, GeometryResult::Compute(conf),
                "GeometryResult", out)) return out;
    if (!Attach(conf, SpatialIndexResult::Compute(conf),
                "SpatialIndexResult", out)) return out;
    if (!Attach(conf, EnrichmentResult::Compute(conf),
                "EnrichmentResult", out)) return out;

    // DSSP
    if (!opts.skip_dssp) {
        auto dssp = DsspResult::Compute(conf);
        if (dssp) {
            Attach(conf, std::move(dssp), "DsspResult", out);
        } else {
            OperationLog::Error("OperationRunner",
                "DSSP failed — HBond calculator will be skipped");
        }
    }

    // Charges
    if (opts.charge_source) {
        if (!Attach(conf, ChargeAssignmentResult::Compute(conf, *opts.charge_source),
                    "ChargeAssignmentResult", out)) return out;
    }

    // --- Tier 0.5: external tools (when charges available) ---

    if (conf.HasResult<ChargeAssignmentResult>()) {
        // MOPAC: PM7+MOZYME semiempirical charges + bond orders
        if (!opts.skip_mopac) {
            auto mopac = MopacResult::Compute(conf, opts.net_charge);
            if (mopac) {
                Attach(conf, std::move(mopac), "MopacResult", out);
            } else {
                OperationLog::Error("OperationRunner",
                    "MOPAC failed (atoms=" +
                    std::to_string(conf.AtomCount()) + ")");
            }
        }

        // APBS: solvated Poisson-Boltzmann E-field
        if (!opts.skip_apbs) {
            auto apbs = ApbsFieldResult::Compute(conf);
            if (apbs) {
                Attach(conf, std::move(apbs), "ApbsFieldResult", out);
            } else {
                OperationLog::Error("OperationRunner",
                    "APBS failed — solvated fields unavailable");
            }
        }
    }

    // --- Tier 1: 8 classical calculators ---

    if (!Attach(conf, BiotSavartResult::Compute(conf),
                "BiotSavartResult", out)) return out;
    if (!Attach(conf, HaighMallionResult::Compute(conf),
                "HaighMallionResult", out)) return out;
    if (!Attach(conf, McConnellResult::Compute(conf),
                "McConnellResult", out)) return out;
    if (!Attach(conf, RingSusceptibilityResult::Compute(conf),
                "RingSusceptibilityResult", out)) return out;
    if (!Attach(conf, PiQuadrupoleResult::Compute(conf),
                "PiQuadrupoleResult", out)) return out;
    if (!Attach(conf, DispersionResult::Compute(conf),
                "DispersionResult", out)) return out;

    if (conf.HasResult<ChargeAssignmentResult>() && !opts.skip_coulomb) {
        if (!Attach(conf, CoulombResult::Compute(conf),
                    "CoulombResult", out)) return out;
    }

    if (conf.HasResult<MopacResult>()) {
        if (!Attach(conf, MopacCoulombResult::Compute(conf),
                    "MopacCoulombResult", out)) return out;
        if (!Attach(conf, MopacMcConnellResult::Compute(conf),
                    "MopacMcConnellResult", out)) return out;
    }

    if (conf.HasResult<DsspResult>()) {
        if (!Attach(conf, HBondResult::Compute(conf),
                    "HBondResult", out)) return out;
    }

    // SASA: per-atom Shrake-Rupley (needs SpatialIndexResult)
    if (!Attach(conf, SasaResult::Compute(conf),
                "SasaResult", out)) return out;

    // EEQ: geometry-dependent charges (Caldeweyher 2019, pure C++/Eigen).
    // No dependencies beyond protein geometry and element types.
    if (!Attach(conf, EeqResult::Compute(conf),
                "EeqResult", out)) return out;

    // AIMNet2: neural network charges + EFG (geometry-only, CUDA)
    // FAILURE POLICY: if model is loaded, AIMNet2 MUST succeed.
    // Silent degradation on a 4-week fleet run means no thesis.
    if (opts.aimnet2_model) {
        auto aimnet2 = AIMNet2Result::Compute(conf, *opts.aimnet2_model);
        if (!Attach(conf, std::move(aimnet2), "AIMNet2Result", out))
            return out;
    }

    // Explicit solvent calculators (trajectory path with full-system .xtc)
    // Explicit solvent calculators: if solvent data is provided, these MUST succeed.
    // A fleet run with silently missing water features is unrecoverable.
    if (opts.solvent && !opts.solvent->Empty()) {
        if (!Attach(conf, WaterFieldResult::Compute(conf, *opts.solvent),
                    "WaterFieldResult", out)) return out;
        if (!Attach(conf, HydrationShellResult::Compute(conf, *opts.solvent),
                    "HydrationShellResult", out)) return out;
        // HydrationGeometryResult: SASA-normal reference frame for water polarisation.
        // Depends on SasaResult (attached above) and SolventEnvironment.
        if (!Attach(conf, HydrationGeometryResult::Compute(conf, *opts.solvent),
                    "HydrationGeometryResult", out)) return out;
    }

    // GROMACS energy: from preloaded EDR via run context (O(1) per frame)
    if (opts.frame_energy) {
        if (!Attach(conf, GromacsEnergyResult::Compute(conf, *opts.frame_energy),
                    "GromacsEnergyResult", out)) return out;
    }

    // Bonded energy: per-atom decomposition from force field parameters
    if (opts.bonded_params) {
        if (!Attach(conf, BondedEnergyResult::Compute(conf, *opts.bonded_params),
                    "BondedEnergyResult", out)) return out;
    }

    // --- Tier 2: DFT comparison (optional) ---

    if (!opts.orca_nmr_path.empty()) {
        auto orca = OrcaShieldingResult::Compute(conf, opts.orca_nmr_path);
        if (orca) {
            Attach(conf, std::move(orca), "OrcaShieldingResult", out);
        } else {
            OperationLog::Error("OperationRunner",
                "ORCA shielding load failed for " + opts.orca_nmr_path);
        }
    }

    OperationLog::Info(LogCalcOther, "OperationRunner",
        "Run complete: " + std::to_string(out.attached.size()) + " results");

    return out;
}


// =================================================================
// RunMutantComparison: use case C.
// =================================================================

RunResult OperationRunner::RunMutantComparison(
        ProteinConformation& wt_conf,
        const RunOptions& wt_opts,
        ProteinConformation& ala_conf,
        const RunOptions& ala_opts) {

    RunResult out;

    OperationLog::Scope scope("OperationRunner::RunMutantComparison",
        "wt_atoms=" + std::to_string(wt_conf.AtomCount()) +
        " ala_atoms=" + std::to_string(ala_conf.AtomCount()));

    // Run standard sequence on WT
    auto wt_result = Run(wt_conf, wt_opts);
    if (!wt_result.Ok()) {
        out.error = "WT: " + wt_result.error;
        return out;
    }
    out.attached.insert(out.attached.end(),
        wt_result.attached.begin(), wt_result.attached.end());

    // Run standard sequence on ALA
    auto ala_result = Run(ala_conf, ala_opts);
    if (!ala_result.Ok()) {
        out.error = "ALA: " + ala_result.error;
        return out;
    }

    // Mutation delta (WT - ALA), attaches to WT
    if (wt_conf.HasResult<OrcaShieldingResult>() &&
        ala_conf.HasResult<OrcaShieldingResult>()) {

        auto delta = MutationDeltaResult::Compute(wt_conf, ala_conf);
        if (delta) {
            Attach(wt_conf, std::move(delta), "MutationDeltaResult", out);
        } else {
            OperationLog::Error("OperationRunner",
                "MutationDeltaResult computation failed");
        }
    }

    OperationLog::Info(LogCalcOther, "OperationRunner",
        "MutantComparison complete: " +
        std::to_string(out.attached.size()) + " results on WT");

    return out;
}


// =================================================================
// RunEnsemble: use case D with trajectories.
// =================================================================

std::vector<RunResult> OperationRunner::RunEnsemble(
        Protein& protein,
        const RunOptions& opts) {

    std::vector<RunResult> results;
    results.reserve(protein.ConformationCount());

    for (size_t i = 0; i < protein.ConformationCount(); ++i) {
        auto& conf = protein.ConformationAt(i);
        results.push_back(Run(conf, opts));
    }

    return results;
}


}  // namespace nmr

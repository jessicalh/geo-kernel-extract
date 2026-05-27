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
#include "TripeptideBackboneShieldingResult.h"
#include "LarsenHBondShieldingResult.h"
#include "TripeptideNeighborShieldingResult.h"
#include "AIMNet2ChargeResponseGradientResult.h"
#include "PlanarGeometryResult.h"
#include "SasaResult.h"
#include "GromacsEnergyResult.h"
#include "GromacsFramePullResult.h"
#include "BondedEnergyResult.h"
#include "WaterFieldResult.h"
#include "HydrationShellResult.h"
#include "HydrationGeometryResult.h"
#include "EeqResult.h"
#include "OperationLog.h"

namespace nmr {


// Helper: attach a result, log failure, accumulate name list.
// Attach itself just moves the pointer; compute timing belongs to the
// caller or to the result's own Compute implementation.
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

// Timed compute + attach. The Scope emits [BEGIN] and [END] with
// elapsed ms through OperationLog for the calculators routed here.
template<typename F>
static bool TimedAttach(ProteinConformation& conf, const char* name,
                        RunResult& out, F&& compute) {
    OperationLog::Scope scope(name);
    auto result = compute();
    return Attach(conf, std::move(result), name, out);
}


// =================================================================
// Run: the standard sequence.
//
// Order:
//   Tier 0 — foundation
//   Tier 0.5 — external tools (need charges)
//   Tier 1 — calculator stack
//   Tier 2 — DFT comparison (optional)
// =================================================================

RunResult OperationRunner::Run(ProteinConformation& conf,
                               const RunOptions& opts) {
    RunResult out;

    OperationLog::Scope scope("OperationRunner::Run",
        "atoms=" + std::to_string(conf.AtomCount()));

    // --- Tier 0: foundation ---

    if (!TimedAttach(conf, "GeometryResult", out, [&]{
            return GeometryResult::Compute(conf); })) return out;
    if (!TimedAttach(conf, "SpatialIndexResult", out, [&]{
            return SpatialIndexResult::Compute(conf); })) return out;
    if (!TimedAttach(conf, "EnrichmentResult", out, [&]{
            return EnrichmentResult::Compute(conf); })) return out;

    // Planar geometry — DSSP-style soft attach. Compute returns
    // nullptr on stub-fixture paths that don't populate substrate;
    // we log and continue rather than aborting the pipeline. Mirrors
    // the existing pattern for DSSP/MOPAC/APBS (calculators that can
    // legitimately fail under specific input shapes), per the
    // PATTERNS.md "Diagnostic error messages" rule.
    {
        OperationLog::Scope pgr_scope("PlanarGeometryResult");
        auto pgr = PlanarGeometryResult::Compute(conf);
        if (pgr) {
            Attach(conf, std::move(pgr), "PlanarGeometryResult", out);
        } else {
            OperationLog::Error("OperationRunner",
                "PlanarGeometryResult failed — pyramidalization, "
                "omega, aromatic_chi2, and pucker NPYs will be absent.");
        }
    }

    // DSSP
    if (!opts.skip_dssp) {
        OperationLog::Scope dssp_scope("DsspResult");
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
        if (!TimedAttach(conf, "ChargeAssignmentResult", out, [&]{
                return ChargeAssignmentResult::Compute(conf, *opts.charge_source); })) return out;
    }

    // --- Tier 0.5: external tools (when charges available) ---

    if (conf.HasResult<ChargeAssignmentResult>()) {
        if (!opts.skip_mopac) {
            OperationLog::Scope mopac_scope("MopacResult");
            auto mopac = MopacResult::Compute(conf, opts.net_charge);
            if (mopac) {
                Attach(conf, std::move(mopac), "MopacResult", out);
            } else {
                OperationLog::Error("OperationRunner",
                    "MOPAC failed (atoms=" +
                    std::to_string(conf.AtomCount()) + ")");
            }
        }

        if (!opts.skip_apbs) {
            OperationLog::Scope apbs_scope("ApbsFieldResult");
            auto apbs = ApbsFieldResult::Compute(conf);
            if (apbs) {
                Attach(conf, std::move(apbs), "ApbsFieldResult", out);
            } else {
                OperationLog::Error("OperationRunner",
                    "APBS failed — solvated fields unavailable");
            }
        }
    }

    // --- Tier 1: calculator stack ---

    if (!TimedAttach(conf, "BiotSavartResult", out, [&]{
            return BiotSavartResult::Compute(conf); })) return out;
    if (!TimedAttach(conf, "HaighMallionResult", out, [&]{
            return HaighMallionResult::Compute(conf); })) return out;
    if (!TimedAttach(conf, "McConnellResult", out, [&]{
            return McConnellResult::Compute(conf); })) return out;
    if (!TimedAttach(conf, "RingSusceptibilityResult", out, [&]{
            return RingSusceptibilityResult::Compute(conf); })) return out;
    if (!TimedAttach(conf, "PiQuadrupoleResult", out, [&]{
            return PiQuadrupoleResult::Compute(conf); })) return out;
    if (!TimedAttach(conf, "DispersionResult", out, [&]{
            return DispersionResult::Compute(conf); })) return out;

    if (conf.HasResult<ChargeAssignmentResult>() && !opts.skip_coulomb) {
        if (!TimedAttach(conf, "CoulombResult", out, [&]{
                return CoulombResult::Compute(conf); })) return out;
    }

    if (conf.HasResult<MopacResult>()) {
        if (!TimedAttach(conf, "MopacCoulombResult", out, [&]{
                return MopacCoulombResult::Compute(conf); })) return out;
        if (!TimedAttach(conf, "MopacMcConnellResult", out, [&]{
                return MopacMcConnellResult::Compute(conf); })) return out;
    }

    if (conf.HasResult<DsspResult>()) {
        if (!TimedAttach(conf, "HBondResult", out, [&]{
                return HBondResult::Compute(conf); })) return out;
    }

    if (!TimedAttach(conf, "SasaResult", out, [&]{
            return SasaResult::Compute(conf); })) return out;
    if (!TimedAttach(conf, "EeqResult", out, [&]{
            return EeqResult::Compute(conf); })) return out;

    // AIMNet2: FAILURE POLICY: if model is loaded, AIMNet2 MUST succeed.
    // ChargeResponseGradient runs unconditionally after AIMNet2Result;
    // chains via Dependencies() for attach ordering and runs its own
    // forward+backward pass. Trajectory mode supplies the model through
    // Trajectory::Run's per-frame RunOptions.
    if (opts.aimnet2_model) {
        if (!TimedAttach(conf, "AIMNet2Result", out, [&]{
                return AIMNet2Result::Compute(conf, *opts.aimnet2_model); })) return out;
        if (!TimedAttach(conf, "AIMNet2ChargeResponseGradientResult", out, [&]{
                return AIMNet2ChargeResponseGradientResult::Compute(
                    conf, *opts.aimnet2_model); })) return out;
    }

    // Explicit solvent calculators: MUST succeed if solvent provided.
    if (opts.solvent && !opts.solvent->Empty()) {
        if (!TimedAttach(conf, "WaterFieldResult", out, [&]{
                return WaterFieldResult::Compute(conf, *opts.solvent); })) return out;
        if (!TimedAttach(conf, "HydrationShellResult", out, [&]{
                return HydrationShellResult::Compute(conf, *opts.solvent); })) return out;
        if (!TimedAttach(conf, "HydrationGeometryResult", out, [&]{
                return HydrationGeometryResult::Compute(conf, *opts.solvent); })) return out;
    }

    // GROMACS frame pull: catch-all for what the trajectory frame
    // yielded. Source-available capture (velocities, box matrix today;
    // future fields land here, not as direct fields on
    // ProteinConformation). No dependencies; runs whenever either
    // pointer is non-null. See GromacsFramePullResult.h.
    if (opts.velocities || opts.box_matrix) {
        if (!TimedAttach(conf, "GromacsFramePullResult", out, [&]{
                return GromacsFramePullResult::Compute(conf, opts.velocities, opts.box_matrix); })) return out;
    }

    // GROMACS energy: from preloaded EDR via run context (O(1) per frame)
    if (opts.frame_energy) {
        if (!TimedAttach(conf, "GromacsEnergyResult", out, [&]{
                return GromacsEnergyResult::Compute(conf, *opts.frame_energy); })) return out;
    }

    // Bonded energy: per-atom decomposition from force field parameters
    if (opts.bonded_params) {
        if (!TimedAttach(conf, "BondedEnergyResult", out, [&]{
                return BondedEnergyResult::Compute(conf, *opts.bonded_params); })) return out;
    }

    // Tripeptide DFT shielding: ProCS15 tensorcs15 backbone σ_BB^i +
    // Larsen-2015-Eq.3 neighbour correction Δσ_BB^{i±1}. Always-on
    // when the table is loaded (Session populates from
    // ~/.nmr_tools.toml [databases].tensorcs15). Per-residue
    // perception failures inside the calculator are logged via
    // OperationLog::Warn but do not abort the run; misses degrade
    // gracefully (per-atom tripeptide_bb_has_match=false stays default).
    if (opts.tripeptide_dft_table) {
        if (!TimedAttach(conf, "TripeptideBackboneShieldingResult", out, [&]{
                return TripeptideBackboneShieldingResult::Compute(
                    conf, *opts.tripeptide_dft_table); })) return out;
        if (!TimedAttach(conf, "TripeptideNeighborShieldingResult", out, [&]{
                return TripeptideNeighborShieldingResult::Compute(
                    conf, *opts.tripeptide_dft_table); })) return out;
    }

    // LarsenHBondShieldingResult: Larsen 2015 ProCS15 H-bond terms via
    // direct DFT grid lookup. Spatial-enumeration based — H-bond
    // candidates come from SpatialIndexResult rather than DSSP, which
    // covers both amide H AND Hα donors uniformly and matches Larsen's
    // geometric H-bond criterion. Runs in parallel with the kernel-
    // form HBondResult (methods accumulate; feedback_methods_accumulate).
    if (opts.larsen_hbond_grid) {
        if (!TimedAttach(conf, "LarsenHBondShieldingResult", out, [&]{
                return LarsenHBondShieldingResult::Compute(
                    conf, *opts.larsen_hbond_grid); })) return out;
    }

    // --- Tier 2: DFT comparison (optional) ---

    if (!opts.orca_nmr_path.empty()) {
        auto orca = OrcaShieldingResult::Compute(conf, opts.orca_nmr_path);
        if (orca) {
            Attach(conf, std::move(orca), "OrcaShieldingResult", out);
        } else {
            // Fail loud: a non-empty --orca path that fails to load is an
            // error, not a silent WT-only success (OI-001). The caller asked
            // for ORCA shielding by supplying the path.
            out.error = "ORCA shielding load failed for " + opts.orca_nmr_path;
            OperationLog::Error("OperationRunner", out.error);
            return out;
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
// RunEnsemble: run a shared option set over existing conformations.
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

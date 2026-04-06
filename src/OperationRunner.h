#pragma once
//
// OperationRunner: the ONE home for all ordered sequences.
//
// THIS IS ONLY FOR ORDER. The conformation is the buffer where
// results accumulate. The runner does not hold state, does not
// cache intermediate results, does not decide what to compute.
// It sequences operations. If a step's prerequisites are on the
// conformation, the step runs. If not, it is skipped. Nothing
// runs backwards. Nothing runs twice. These are the ONLY
// variations — the presence or absence of upstream results.
//
// Every use case (UI, CLI, batch, training) calls one of
// these methods. The conformation arrives with its protein
// identity determined (protonated, topology resolved). The
// runner fills it with computed results in dependency order.
//
// All ordered sequences live here. Future agents add new
// sequences here, not in new files.
//
// Replaces both Pipeline.h and CalculationRunner.h.
//

#include "ProteinConformation.h"
#include "ChargeSource.h"
#include <vector>
#include <string>

namespace nmr {

class Protein;

struct RunOptions {
    // Charges — required for real physics.
    // Null = no Coulomb, no APBS, no MOPAC.
    const ChargeSource* charge_source = nullptr;
    int net_charge = 0;

    // Skip DSSP (and therefore HBond).
    bool skip_dssp = false;

    // DFT: load ORCA shielding tensors after calculators.
    std::string orca_nmr_path;
};

struct RunResult {
    std::vector<std::string> attached;
    std::string error;
    bool Ok() const { return error.empty(); }
};


class OperationRunner {
public:

    // =================================================================
    // Run: the standard sequence (use cases A, B, D-per-item).
    //
    // Tier 0: Geometry, SpatialIndex, Enrichment, DSSP
    //         Charges (if charge_source provided)
    //         MOPAC (if charges available)
    //         APBS (if charges)
    // Tier 1: 8 classical calculators
    //         Coulomb (if charges), HBond (if DSSP)
    // Tier 2: ORCA DFT (if orca_nmr_path provided)
    // =================================================================

    static RunResult Run(ProteinConformation& conf,
                         const RunOptions& opts);


    // =================================================================
    // RunMutantComparison: use case C.
    //
    // Runs the standard sequence on BOTH WT and ALA conformations,
    // then computes MutationDeltaResult (attaches to WT).
    // =================================================================

    static RunResult RunMutantComparison(
        ProteinConformation& wt_conf,
        const RunOptions& wt_opts,
        ProteinConformation& ala_conf,
        const RunOptions& ala_opts);


    // =================================================================
    // RunEnsemble: use case D with trajectories.
    //
    // Runs the standard sequence on every conformation of a protein.
    // Returns per-frame results.
    // =================================================================

    static std::vector<RunResult> RunEnsemble(
        Protein& protein,
        const RunOptions& opts);
};

}  // namespace nmr

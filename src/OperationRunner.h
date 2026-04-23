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
#include "SolventEnvironment.h"
#include <vector>
#include <string>

namespace nmr {

class Protein;

struct AIMNet2Model;  // forward declaration (AIMNet2Result.h)
struct BondedParameters;  // forward declaration (BondedEnergyResult.h)
struct GromacsEnergy;     // forward declaration (GromacsEnergyResult.h)

struct RunOptions {
    // Charges — required for real physics.
    // Null = no Coulomb, no APBS, no MOPAC.
    const ChargeSource* charge_source = nullptr;
    int net_charge = 0;

    // Skip DSSP (and therefore HBond).
    bool skip_dssp = false;

    // Skip MOPAC semiempirical (and therefore MopacCoulomb, MopacMcConnell).
    bool skip_mopac = false;

    // Skip APBS Poisson-Boltzmann.
    bool skip_apbs = false;

    // Skip vacuum Coulomb EFG (home-rolled, O(N*k), 25s at 4800 atoms).
    // APBS is faster and solvated — preferred for all production work.
    // Coulomb retained for special comparison tests only.
    bool skip_coulomb = false;

    // AIMNet2: loaded model for neural network charges + EFG.
    // Null = skip AIMNet2. Loaded once, shared across all conformations.
    AIMNet2Model* aimnet2_model = nullptr;

    // Per-frame energy from GROMACS .edr (preloaded by Trajectory).
    // Null = skip GromacsEnergyResult. O(1) per frame.
    // Set by GromacsFrameHandler via trajectory->EnergyAtTime() after
    // each frame.
    const GromacsEnergy* frame_energy = nullptr;

    // Explicit solvent: water + ion positions for this frame.
    // Null = no solvent calculators (protein-only trajectory).
    // Set by the full-system trajectory reader.
    const SolventEnvironment* solvent = nullptr;

    // Bonded force field parameters from TPR for per-atom energy
    // decomposition. Null = skip BondedEnergyResult.
    // Owned by TrajectoryProtein, borrowed here per frame.
    const BondedParameters* bonded_params = nullptr;

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
    // Tier 1: classical calculators
    //         Coulomb (if charges and not skip_coulomb), HBond (if DSSP)
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

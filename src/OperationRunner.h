#pragma once
//
// OperationRunner: the ONE home for all ordered sequences.
//
// THIS IS ONLY FOR ORDER. The conformation is the buffer where
// results accumulate. The runner does not hold state and does not
// cache intermediate results. It applies one ordered compute policy:
// RunOptions and already-attached prerequisites decide which optional
// steps run. Nothing runs backwards. Nothing runs twice. These are the
// variations — options plus the presence or absence of upstream results.
//
// In-tree dispatchers call these methods with a conformation whose
// protein identity is determined (protonated, topology resolved). The
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
class  TripeptideDftTable; // forward declaration (TripeptideDftTable.h)
class  LarsenHBondGrid;    // forward declaration (LarsenHBondGrid.h)

struct RunOptions {
    // Charges — required for charge-dependent calculators.
    // Null = no Coulomb, no APBS, no MOPAC.
    const ChargeSource* charge_source = nullptr;
    int net_charge = 0;

    // Skip DSSP (and therefore HBond).
    bool skip_dssp = false;

    // Skip MOPAC semiempirical (and therefore MopacCoulomb, MopacMcConnell).
    bool skip_mopac = false;

    // Skip APBS Poisson-Boltzmann.
    bool skip_apbs = false;

    // Skip vacuum Coulomb EFG (home-rolled, O(N*k)).
    // Retired from production: APBS is the canonical electrostatics, so
    // every single-frame mode and PerFrameExtractionSet set this true.
    // The only consumer that sets it false is FullFatFrameExtraction,
    // where CoulombResult feeds the MOPAC-vs-FF14SB reconciliation probe.
    bool skip_coulomb = false;

    // AIMNet2: loaded model for neural network charges + EFG.
    // Null = skip AIMNet2. Loaded once, shared across all conformations.
    // When set, OperationRunner::Run also attaches
    // AIMNet2ChargeResponseGradientResult automatically (single
    // forward+backward pass).
    AIMNet2Model* aimnet2_model = nullptr;

    // Per-frame energy from GROMACS .edr (preloaded by Trajectory).
    // Null = skip GromacsEnergyResult. O(1) per frame.
    // Set by Trajectory::Run via EnergyAtTime() after each frame.
    const GromacsEnergy* frame_energy = nullptr;

    // Explicit solvent: water + ion positions for this frame.
    // Null = no solvent calculators (protein-only trajectory).
    // Set by the full-system trajectory reader.
    const SolventEnvironment* solvent = nullptr;

    // Bonded force field parameters from TPR for per-atom energy
    // decomposition. Null = skip BondedEnergyResult.
    // Owned by TrajectoryProtein, borrowed here per frame.
    const BondedParameters* bonded_params = nullptr;

    // Per-frame TRR data from GromacsFrameHandler (set per frame from
    // TrajectoryEnv by Trajectory::Run). Null = no velocities / box for
    // this frame (XTC-only legacy or non-trajectory load).
    // Read by GromacsFramePullResult, which stashes them on its result
    // alongside derived geometry — the catch-all for "everything gromacs
    // gave us at frame pull time."
    const std::vector<Vec3>* velocities = nullptr;
    const Eigen::Matrix3d*   box_matrix = nullptr;

    // DFT: load ORCA shielding tensors after calculators.
    std::string orca_nmr_path;

    // Tripeptide DFT lookup table (ProCS15 tensorcs15 replica).
    // Null = skip tripeptide backbone + neighbor shielding results.
    // When set, OperationRunner::Run attaches both
    // TripeptideBackboneShieldingResult and
    // TripeptideNeighborShieldingResult. Owned by Session for the
    // process lifetime; borrowed here per run.
    const TripeptideDftTable* tripeptide_dft_table = nullptr;

    // LarsenHBondGrid pointer. Borrowed; Session owns.
    // Null = skip LarsenHBondShieldingResult.
    // When set, OperationRunner::Run attaches the calculator after the
    // tripeptide block. Methods accumulate side-by-side with the
    // kernel-form HBondResult (feedback_methods_accumulate).
    const LarsenHBondGrid* larsen_hbond_grid = nullptr;
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
    // Foundation: Geometry, SpatialIndex, Enrichment, PlanarGeometry,
    //             DSSP, Charges.
    // Calculator stack: classical shielding, Sasa, Eeq, optional
    //             MOPAC/APBS/Coulomb/HBond/AIMNet2, solvent, GROMACS,
    //             tripeptide, and Larsen resources.
    // DFT comparison: ORCA shielding tensors if orca_nmr_path is provided.
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
    // RunEnsemble: run a shared option set over existing conformations.
    //
    // Runs the standard sequence on every conformation of a protein.
    // Returns one RunResult per conformation.
    // =================================================================

    static std::vector<RunResult> RunEnsemble(
        Protein& protein,
        const RunOptions& opts);
};

}  // namespace nmr

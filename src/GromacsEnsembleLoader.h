#pragma once
//
// GromacsEnsembleLoader: load a protein ensemble from GROMACS fleet data.
//
// Input: one protein directory from sampled_poses/ plus the production
// TPR file from run_params_export/.
//
// The TPR (read via libgromacs) is the authority for:
//   - Atom names (CHARMM naming)
//   - Residue names (partially CHARMM: HSP preserved, HSE→"HIS")
//   - Per-atom charges (CHARMM36m)
//   - Per-atom element (from atomic number)
//   - Residue boundaries (from molecule topology)
//
// The pose PDBs provide positions only. Topology is from the TPR and
// the first pose's geometry (for bond detection via CovalentTopology).
//
// ToolContext is DECLARED as Charmm, not inferred.
// ForceField is DECLARED by the caller, not read from the data.
//
// Output: a Protein with N MDFrameConformations (typically 10).
// Charges are returned separately — the caller attaches them via
// ChargeAssignmentResult with whatever ChargeSource is appropriate.
//

#include "Protein.h"
#include "BuildResult.h"
#include "OperationRunner.h"
#include <string>
#include <vector>

namespace nmr {

struct FleetPaths {
    std::string sampled_poses_dir;   // e.g., .../sampled_poses/1A6J_5789/
    std::string tpr_path;            // e.g., .../run_params_export/1A6J_5789/prod.tpr
    ForceField force_field = ForceField::CHARMM36m;  // declared, not inferred
};


// Load a protein ensemble from fleet data.
//
// Reads the TPR via libgromacs for topology and charges, ensemble.json
// for frame metadata, and pose PDBs for positions. The first pose
// provides geometry for CovalentTopology::Resolve. All poses become
// MDFrameConformations.
//
// Charges from the TPR are wrapped in a PreloadedChargeSource inside
// the returned BuildResult. Net charge computed from the TPR charge sum.
BuildResult BuildFromGromacs(const FleetPaths& paths);

// Build a Protein from TPR topology only (no poses, no conformations).
// Returns a Protein with residues + atoms but no conformations.
// Caller provides positions separately for FinalizeConstruction.
// Charges are in the BuildResult.
// Used by BuildFromGromacs (fleet path). The trajectory path uses
// FullSystemReader::BuildProtein() which shares the same TPR parse.
BuildResult BuildProteinFromTpr(const std::string& tpr_path,
                                const std::string& protein_id,
                                ForceField force_field = ForceField::CHARMM36m);


// Run all frames through the standard sequence.
// Delegates to OperationRunner::RunEnsemble.
// Each MDFrameConformation is independent.
std::vector<RunResult> RunAllFrames(
    Protein& protein, const RunOptions& opts);

}  // namespace nmr

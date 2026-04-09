#pragma once
//
// BuildResult: the common output of all protein builders.
//
// Every loading path (PDB, ORCA, GROMACS) produces a BuildResult.
// The protein is fully constructed: all atoms present (including H
// if protonated), topology resolved, at least one conformation created.
// The ChargeSource matches the protein's atoms. The net_charge is
// the integer formal charge determined by the protonation state.
//
// Builders:
//   BuildFromPdb     — PDB file, optionally protonated via reduce
//   BuildFromOrca    — ORCA DFT run (prmtop + XYZ), already protonated
//   BuildFromGromacs — GROMACS ensemble (TPR + pose PDBs), already protonated
//

#include "ChargeSource.h"
#include <memory>
#include <string>
#include <vector>

namespace nmr {

class Protein;

struct BuildResult {
    std::unique_ptr<Protein> protein;       // fully constructed, topology resolved
    std::unique_ptr<ChargeSource> charges;  // matched to this protein's atoms
    int net_charge = 0;                     // formal charge (integer)
    bool ok = false;
    std::string error;

    // Per-conformation names from source PDB filenames (sans extension).
    // Empty for PDB/ORCA builds. Populated by BuildFromGromacs when
    // PDBs have descriptive names (e.g. fes-sampler output).
    // RunFleet uses these as output subdirectory bases when present.
    std::vector<std::string> pose_names;

    // Convenience: check success
    bool Ok() const { return ok && protein != nullptr; }
    explicit operator bool() const { return Ok(); }
};

}  // namespace nmr

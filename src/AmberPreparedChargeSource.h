#pragma once
//
// AmberPreparedChargeSource: produces an authoritative AMBER PRMTOP at
// run time via tleap when the flat ff14SB parameter file cannot
// represent the protein, then reads CHARGE/RADII back onto the
// extractor's atoms.
//
// Pipeline (LoadCharges):
//   1. Hard preconditions (std::abort on violation; re-tleap protection):
//        - protein.HasForceFieldCharges() must be false
//        - protein.BuildContext().prmtop_path must be empty
//   2. Resolve tleap binary (config_.tleap_path or RuntimeEnvironment::Tleap()).
//   3. Make a deterministic work directory under TmpDir.
//   4. Generate input PDB from typed Protein state via
//      amber_leap::GenerateAmberPdb (handles AMBER residue naming for
//      HID/HIE/HIP, CYX from typed disulfide detection, ASH/GLH/CYM/LYN
//      variant naming, and ACE/NME caps under
//      UseCappedFragmentsForUnsupportedTerminalVariants).
//   5. Generate LEaP script via amber_leap::GenerateLeapScript.
//   6. Run tleap as a subprocess; parse the resulting PRMTOP for
//      ATOM_NAME / RESIDUE_LABEL / RESIDUE_POINTER / CHARGE / RADII.
//   7. Map PRMTOP atoms back to extractor atom indices via the
//      ResidueAmberMapping built in step 4 (NONE_FOR_CAP-tagged
//      residues are dropped from the mapping).
//
// Configurable via AmberSourceConfig:
//   - preparation_policy (FailOnUnsupported | UseCappedFragments | UseStockTermini)
//   - tleap_path (overrides RuntimeEnvironment::Tleap() if non-empty)
//   - work_dir   (overrides RuntimeEnvironment::TempFilePath if non-empty)
//

#include "AmberChargeResolver.h"
#include "AmberLeapInput.h"
#include "ChargeSource.h"
#include <memory>
#include <string>
#include <vector>

namespace nmr {

class Protein;

class AmberPreparedChargeSource : public ChargeSource {
public:
    AmberPreparedChargeSource(const Protein& protein,
                              AmberPreparationPolicy policy,
                              AmberFlatTableCoverageVerdict reason,
                              AmberSourceConfig config);

    ForceField SourceForceField() const override {
        return ForceField::Amber_ff14SB;
    }
    ChargeModelKind Kind() const override {
        return ChargeModelKind::AmberPreparedPrmtop;
    }
    std::string Describe() const override;

    std::vector<AtomChargeRadius> LoadCharges(
        const Protein& protein,
        const ProteinConformation& conf,
        std::string& error_out) const override;

    // Inspection accessors (used by tests in step 4 and by the
    // LoadCharges body in step 5).
    AmberPreparationPolicy Policy() const { return policy_; }
    const AmberFlatTableCoverageVerdict& Reason() const { return reason_; }
    const AmberSourceConfig& Config() const { return config_; }

    // Generated PDB body for the given conformation. Side-effect:
    // populates the cached ResidueAmberMapping.
    std::string GeneratedPdb(const ProteinConformation& conf) const;

    // Generated LEaP script body for the given paths. Disulfide
    // bond pairs are derived from the cached ResidueAmberMapping;
    // call GeneratedPdb first or via LoadCharges so the mapping
    // exists.
    std::string GeneratedLeapScript(const std::string& pdb_path,
                                     const std::string& prmtop_path,
                                     const std::string& inpcrd_path) const;

    const amber_leap::ResidueAmberMapping& ResidueMapping() const {
        return residue_mapping_;
    }

private:
    const Protein& protein_;
    AmberPreparationPolicy policy_;
    AmberFlatTableCoverageVerdict reason_;
    AmberSourceConfig config_;
    mutable amber_leap::ResidueAmberMapping residue_mapping_;
    mutable std::vector<std::pair<size_t, size_t>>
        disulfide_extractor_pairs_;
    mutable bool mapping_built_ = false;

    // 1-based PRMTOP residue indices for each disulfide pair, derived
    // from disulfide_extractor_pairs_ projected through residue_mapping_.
    std::vector<std::pair<size_t, size_t>>
        DisulfidePairs1Based() const;
};

}  // namespace nmr

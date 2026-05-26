#pragma once
//
// Runtime AMBER PRMTOP preparation for proteins the flat ff14SB table cannot
// represent. LoadCharges requires a protein without existing force-field
// charges and without build_context.prmtop_path, writes a generated PDB and
// LEaP script, runs tleap, then maps CHARGE/RADII back to extractor atoms.
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

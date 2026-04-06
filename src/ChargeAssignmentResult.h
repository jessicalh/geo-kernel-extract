#pragma once
//
// ChargeAssignmentResult: per-atom partial charges and VdW radii.
//
// Dependencies: none.
//
// Charges come from a typed ChargeSource: ff14SB param file, CHARMM36m
// topology (.tpr), AMBER prmtop, or stub. The ChargeSource determines
// the force field; this result stores the numbers on ConformationAtom.
//
// For protonation variants (HID/HIE/HIP, ASH, GLH, CYX, LYN, ARN):
// the ff14SB param file path uses the residue's protonation_variant_index
// to select the correct charge set. Other sources (tpr, prmtop) carry
// the protonation in their topology already.
//

#include "ConformationResult.h"
#include "ProteinConformation.h"
#include <unordered_map>
#include <string>

namespace nmr {

class ChargeSource;

class ChargeAssignmentResult : public ConformationResult {
public:
    std::string Name() const override { return "ChargeAssignmentResult"; }
    std::vector<std::type_index> Dependencies() const override { return {}; }

    // Typed factory: assign charges from any ChargeSource.
    // This is the preferred path. The ChargeSource determines the
    // force field and data format.
    static std::unique_ptr<ChargeAssignmentResult> Compute(
        ProteinConformation& conf,
        const ChargeSource& source);

    // ff14SB param file factory (convenience, delegates to typed path).
    static std::unique_ptr<ChargeAssignmentResult> Compute(
        ProteinConformation& conf,
        const std::string& param_file_path);

    // No-arg Compute removed 2026-04-03.
    // Every protein has real charges. No stub fallback.

    // Query methods
    double ChargeAt(size_t atom_index) const;
    double RadiusAt(size_t atom_index) const;

    // Diagnostics
    double TotalCharge() const { return total_charge_; }
    int AssignedCount() const { return assigned_count_; }
    int UnassignedCount() const { return unassigned_count_; }
    const std::string& Source() const { return source_; }

    // PDB LOADING BOUNDARY: load ff14SB params file into a lookup table.
    // Key is "RESNAME ATOMNAME", value is (charge, radius).
    // Public because ParamFileChargeSource uses it.
    struct ParamEntry {
        double charge = 0.0;
        double radius = 0.0;
    };
    static std::unordered_map<std::string, ParamEntry> LoadParamFile(
        const std::string& path);

    // Map from standard amino acid to ff14SB variant residue name.
    // PDB LOADING BOUNDARY: this is where the typed variant_index becomes
    // the string residue name for the parameter file lookup.
    // Public because ParamFileChargeSource uses it.
    static std::string VariantResidueName(AminoAcid aa, int protonation_variant_index);

private:
    const ProteinConformation* conf_ = nullptr;
    double total_charge_ = 0.0;
    int assigned_count_ = 0;
    int unassigned_count_ = 0;
    std::string source_;
};

}  // namespace nmr

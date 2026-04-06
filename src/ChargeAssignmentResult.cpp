#include "ChargeAssignmentResult.h"
#include "ChargeSource.h"
#include "Protein.h"
#include "OperationLog.h"
#include <fstream>
#include <sstream>
#include <cmath>

namespace nmr {


// ============================================================================
// Typed factory: assign charges from any ChargeSource.
// ============================================================================

std::unique_ptr<ChargeAssignmentResult> ChargeAssignmentResult::Compute(
        ProteinConformation& conf,
        const ChargeSource& source) {

    OperationLog::Scope scope("ChargeAssignmentResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()) +
        " source=" + source.Describe());

    const Protein& protein = conf.ProteinRef();

    std::string error;
    auto charges = source.LoadCharges(protein, conf, error);
    if (charges.empty()) {
        OperationLog::Error("ChargeAssignmentResult::Compute",
            "charge source failed: " + error);
        return nullptr;
    }

    if (charges.size() != conf.AtomCount()) {
        OperationLog::Error("ChargeAssignmentResult::Compute",
            "charge source returned " + std::to_string(charges.size()) +
            " entries, expected " + std::to_string(conf.AtomCount()));
        return nullptr;
    }

    auto result = std::make_unique<ChargeAssignmentResult>();
    result->conf_ = &conf;
    result->source_ = source.Describe();

    double total_q = 0.0;
    int assigned = 0;
    int unassigned = 0;

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        auto& ca = conf.MutableAtomAt(ai);
        ca.partial_charge = charges[ai].charge;
        ca.vdw_radius = charges[ai].radius;
        total_q += charges[ai].charge;

        if (std::abs(charges[ai].charge) > 1e-10)
            assigned++;
        else
            unassigned++;
    }

    result->total_charge_ = total_q;
    result->assigned_count_ = assigned;
    result->unassigned_count_ = unassigned;

    OperationLog::Info(LogCharges, "ChargeAssignmentResult::Compute",
        "source=" + source.Describe() +
        " assigned=" + std::to_string(assigned) +
        " unassigned=" + std::to_string(unassigned) +
        " total_charge=" + std::to_string(total_q));

    return result;
}

// ============================================================================
// PDB LOADING BOUNDARY: load ff14SB parameter file.
//
// File format (one line per atom):
//   RESNAME ATOMNAME CHARGE LJ_EPSILON RADIUS
//
// Lines starting with # are comments. Blank lines are skipped.
// Key is "RESNAME ATOMNAME" -- string lookup at the loading boundary.
// After this function, charges are plain doubles.
// ============================================================================

std::unordered_map<std::string, ChargeAssignmentResult::ParamEntry>
ChargeAssignmentResult::LoadParamFile(const std::string& path) {
    std::unordered_map<std::string, ParamEntry> params;

    std::ifstream in(path);
    if (!in.is_open()) {
        OperationLog::Error("ChargeAssignmentResult::LoadParamFile",
            "cannot open " + path);
        return params;
    }

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string resname, atomname;
        double charge, epsilon, radius;
        if (!(iss >> resname >> atomname >> charge >> epsilon >> radius)) continue;

        // PDB LOADING BOUNDARY: string key for lookup
        std::string key = resname + " " + atomname;
        params[key] = {charge, radius};
    }

    OperationLog::Info(LogCharges, "ChargeAssignmentResult::LoadParamFile",
        "loaded " + std::to_string(params.size()) + " entries from " + path);

    return params;
}


// ============================================================================
// Map AminoAcid enum + protonation_variant_index to ff14SB residue name.
//
// PDB LOADING BOUNDARY: this is where the typed protonation variant becomes
// the string residue name for the parameter file lookup. After this,
// everything is numbers.
//
// The variant names in the ff14SB file: HIE, HID, HIP, ASH, GLH, CYX, LYN.
// Default (no variant): standard three-letter code (HIS, ASP, GLU, CYS, LYS).
// ============================================================================

std::string ChargeAssignmentResult::VariantResidueName(
        AminoAcid aa, int protonation_variant_index) {

    if (protonation_variant_index < 0) {
        // No variant assigned -- use standard name.
        // For HIS without variant, ff14SB file has HIS entries but
        // they may be incomplete. Default to HIE (most common tautomer).
        if (aa == AminoAcid::HIS) return "HIE";
        return ThreeLetterCodeForAminoAcid(aa);
    }

    // PDB LOADING BOUNDARY: variant index to ff14SB residue name
    switch (aa) {
        case AminoAcid::HIS:
            // 0 = HID, 1 = HIE, 2 = HIP (from AminoAcidType.cpp)
            if (protonation_variant_index == 0) return "HID";
            if (protonation_variant_index == 1) return "HIE";
            if (protonation_variant_index == 2) return "HIP";
            return "HIE";
        case AminoAcid::ASP:
            // 0 = ASH
            if (protonation_variant_index == 0) return "ASH";
            return "ASP";
        case AminoAcid::GLU:
            // 0 = GLH
            if (protonation_variant_index == 0) return "GLH";
            return "GLU";
        case AminoAcid::CYS:
            // 0 = CYX, 1 = CYM
            if (protonation_variant_index == 0) return "CYX";
            return "CYS";
        case AminoAcid::LYS:
            // 0 = LYN
            if (protonation_variant_index == 0) return "LYN";
            return "LYS";
        case AminoAcid::ARG:
            // 0 = ARN (deprotonated, very rare)
            if (protonation_variant_index == 0) return "ARN";
            return "ARG";
        default:
            return ThreeLetterCodeForAminoAcid(aa);
    }
}


// ============================================================================
// Compute: ff14SB from param file (convenience, delegates to typed path)
// ============================================================================

std::unique_ptr<ChargeAssignmentResult> ChargeAssignmentResult::Compute(
        ProteinConformation& conf,
        const std::string& param_file_path) {

    ParamFileChargeSource source(param_file_path);
    return Compute(conf, source);
}


// No-arg Compute (stub) removed 2026-04-03.


double ChargeAssignmentResult::ChargeAt(size_t atom_index) const {
    return conf_->AtomAt(atom_index).partial_charge;
}

double ChargeAssignmentResult::RadiusAt(size_t atom_index) const {
    return conf_->AtomAt(atom_index).vdw_radius;
}

}  // namespace nmr

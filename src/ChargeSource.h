#pragma once
//
// ChargeSource: typed abstraction for where per-atom charges come from.
//
// The force field determines charges, naming, and VdW parameters.
// Different force fields produce different numbers for the same atom.
// This is NOT a string dispatch — each source is a distinct type that
// knows how to provide charges for a Protein's atoms.
//
// Implementations:
//   ParamFileChargeSource  — ff14SB from flat parameter file (current)
//   GmxTprChargeSource     — CHARMM36m (or any FF) from GROMACS .tpr
//   PrmtopChargeSource     — AMBER prmtop (future, needs cpptraj)
//   StubChargeSource       — uniform test charges
//
// ChargeAssignmentResult::Compute takes a ChargeSource, not a path string.
//

#include "Types.h"
#include <string>
#include <vector>

namespace nmr {

class Protein;
class ProteinConformation;


// Force field identity. Recorded in ProteinBuildContext as provenance.
enum class ForceField {
    Amber_ff14SB,
    Amber_ff19SB,
    CHARMM36m,
    Unknown
};

inline const char* ForceFieldName(ForceField ff) {
    switch (ff) {
        case ForceField::Amber_ff14SB: return "ff14SB";
        case ForceField::Amber_ff19SB: return "ff19SB";
        case ForceField::CHARMM36m:    return "CHARMM36m";
        case ForceField::Unknown:      return "unknown";
    }
    return "unknown";
}


// Per-atom charge and radius from whatever source.
struct AtomChargeRadius {
    double charge = 0.0;    // elementary charges (e)
    double radius = 0.0;    // Angstroms (VdW radius for APBS)
};


// Abstract charge source. Each implementation knows how to read
// charges from its data format and match them to protein atoms.
class ChargeSource {
public:
    virtual ~ChargeSource() = default;

    // Which force field produced these charges.
    virtual ForceField SourceForceField() const = 0;

    // Human-readable description for logging.
    virtual std::string Describe() const = 0;

    // Load charges for all atoms in this protein/conformation.
    // Returns one entry per atom (parallel to protein's atom list).
    // Empty vector on failure (error reported via error_out).
    virtual std::vector<AtomChargeRadius> LoadCharges(
        const Protein& protein,
        const ProteinConformation& conf,
        std::string& error_out) const = 0;
};


// ============================================================================
// ff14SB from flat parameter file (the current implementation, now typed)
// ============================================================================

class ParamFileChargeSource : public ChargeSource {
public:
    explicit ParamFileChargeSource(const std::string& param_file_path)
        : path_(param_file_path) {}

    ForceField SourceForceField() const override { return ForceField::Amber_ff14SB; }
    std::string Describe() const override { return "ff14SB:" + path_; }

    std::vector<AtomChargeRadius> LoadCharges(
        const Protein& protein,
        const ProteinConformation& conf,
        std::string& error_out) const override;

private:
    std::string path_;
};


// ============================================================================
// CHARMM36m (or any FF) from GROMACS .tpr binary topology.
// Uses `gmx dump` to extract per-atom charges and masses.
// ============================================================================

class GmxTprChargeSource : public ChargeSource {
public:
    // gmx_binary: path to the gmx executable (from RuntimeEnvironment)
    // tpr_path: path to the .tpr file
    GmxTprChargeSource(const std::string& gmx_binary,
                       const std::string& tpr_path,
                       ForceField ff = ForceField::CHARMM36m)
        : gmx_(gmx_binary), tpr_(tpr_path), ff_(ff) {}

    ForceField SourceForceField() const override { return ff_; }
    std::string Describe() const override {
        return std::string(ForceFieldName(ff_)) + ":tpr:" + tpr_;
    }

    std::vector<AtomChargeRadius> LoadCharges(
        const Protein& protein,
        const ProteinConformation& conf,
        std::string& error_out) const override;

private:
    std::string gmx_;
    std::string tpr_;
    ForceField ff_;
};


// ============================================================================
// AMBER prmtop (ff14SB/ff19SB authoritative charges).
// Parses the Fortran-formatted text file directly — no cpptraj dependency.
// Reads %FLAG CHARGE (AMBER internal units, divide by 18.2223 for e)
// and %FLAG RADII (PB radii in Angstroms).
// ============================================================================

class PrmtopChargeSource : public ChargeSource {
public:
    explicit PrmtopChargeSource(const std::string& prmtop_path,
                                ForceField ff = ForceField::Amber_ff14SB)
        : path_(prmtop_path), ff_(ff) {}

    ForceField SourceForceField() const override { return ff_; }
    std::string Describe() const override {
        return std::string(ForceFieldName(ff_)) + ":prmtop:" + path_;
    }

    std::vector<AtomChargeRadius> LoadCharges(
        const Protein& protein,
        const ProteinConformation& conf,
        std::string& error_out) const override;

private:
    std::string path_;
    ForceField ff_;

    // AMBER prmtop charge units: e * 18.2223 (sqrt(kcal/mol * A))
    static constexpr double AMBER_CHARGE_FACTOR = 18.2223;
};


// No StubChargeSource: every protein enters the system protonated
// and charged; there is no fallback to fake charges.



// ============================================================================
// Pre-loaded charges (e.g., from TPR via libgromacs).
// The caller already has the per-atom values; this wraps them in the
// ChargeSource interface so ChargeAssignmentResult can consume them.
// ============================================================================

class PreloadedChargeSource : public ChargeSource {
public:
    PreloadedChargeSource(std::vector<AtomChargeRadius> charges,
                          ForceField ff)
        : charges_(std::move(charges)), ff_(ff) {}

    ForceField SourceForceField() const override { return ff_; }
    std::string Describe() const override {
        return std::string(ForceFieldName(ff_)) + ":preloaded";
    }

    std::vector<AtomChargeRadius> LoadCharges(
        const Protein& protein,
        const ProteinConformation& conf,
        std::string& error_out) const override;

private:
    std::vector<AtomChargeRadius> charges_;
    ForceField ff_;
};


}  // namespace nmr

#pragma once
//
// PropkaProtonator: pKa prediction via PROPKA 3.5.
//
// Calls the propka3 binary on a temp PDB (heavy atoms only — PROPKA
// works on heavy atoms), parses the .pka output for per-residue pKa
// values, and applies Henderson-Hasselbalch at the target pH to produce
// a ProtonationState.
//
// The naming boundary is in Protonate(): PROPKA output uses PDB residue
// names and numbers. We match to our typed Protein residues by sequence
// number and amino acid type. After matching, everything is typed.
//
// PROPKA citation: Olsson et al. (2011) JCTC 7:525.
//

#include "Protonator.h"

namespace nmr {

struct PkaResult {
    std::string residue_type;   // "ASP", "GLU", "HIS", "CYS", "TYR", "LYS"
    int residue_number;
    std::string chain_id;
    double pKa;
};


class PropkaProtonator : public Protonator {
public:
    ProtonationResult Protonate(
        const Protein& protein,
        const ProteinConformation& conf,
        double pH) override;

    std::string ToolName() const override { return "PROPKA"; }
    std::string ToolVersion() const override { return "3.5.1"; }

    // Low-level: predict pKa values without making protonation decisions.
    // Useful for analysis ("what does PROPKA think?") without committing
    // to a ProtonationState.
    static std::vector<PkaResult> PredictPka(
        const Protein& protein,
        const ProteinConformation& conf,
        std::string& error_out);

private:
    // Write a temp PDB for PROPKA (heavy atoms only).
    static std::string WriteTempPdb(
        const Protein& protein,
        const ProteinConformation& conf,
        const std::string& dir);

    // Parse a PROPKA .pka output file.
    static std::vector<PkaResult> ParsePkaFile(const std::string& path);

    // Apply Henderson-Hasselbalch to pKa predictions at a given pH.
    static ProtonationState ApplyHendersonHasselbalch(
        const Protein& protein,
        const std::vector<PkaResult>& pkas,
        double pH);
};

}  // namespace nmr

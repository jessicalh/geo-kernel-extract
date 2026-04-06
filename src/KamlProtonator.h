#pragma once
//
// KamlProtonator: ML-predicted pKa via KaML-CBTrees.
//
// Calls the KaML-CBtree.py script on a temp PDB (heavy atoms only),
// parses the _pka.csv output, and applies Henderson-Hasselbalch.
//
// KaML uses CatBoost models trained on experimental pKa data.
// Structure-dependent: different conformations can yield different pKa.
// Good for re-evaluating protonation on trajectory frames.
//
// ~80% success rate on ARM64 — callers should handle failure gracefully.
//
// Citation: Reis & Machuqueiro (2024) JCIM.
//

#include "Protonator.h"
#include "PropkaProtonator.h"  // PkaResult shared type

namespace nmr {

class KamlProtonator : public Protonator {
public:
    ProtonationResult Protonate(
        const Protein& protein,
        const ProteinConformation& conf,
        double pH) override;

    std::string ToolName() const override { return "KaML"; }
    std::string ToolVersion() const override { return "KaML-CBTrees"; }

    // Low-level: predict pKa values without making protonation decisions.
    static std::vector<PkaResult> PredictPka(
        const Protein& protein,
        const ProteinConformation& conf,
        std::string& error_out);

private:
    // Write temp PDB (same as PropkaProtonator — heavy atoms only).
    static std::string WriteTempPdb(
        const Protein& protein,
        const ProteinConformation& conf,
        const std::string& dir);

    // Parse KaML _pka.csv output.
    // Format: "Residue, pKa\nASP-43,3.245\nHIS-24,6.789\n..."
    static std::vector<PkaResult> ParsePkaCsv(const std::string& path);
};

}  // namespace nmr

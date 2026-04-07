#pragma once
//
// OrcaShieldingResult: DFT shielding tensors from ORCA NMR calculation.
//
// Dependencies: none (loaded from external ORCA output file).
//
// Stores per-atom shielding tensors (diamagnetic, paramagnetic, total)
// as Mat3 + SphericalTensor on ConformationAtom. This is a feature input
// for calibration against DFT and for T2 residual analysis.
//
// Each protein gets its own OrcaShieldingResult on its own conformation.
// WT and mutant are separate proteins with separate results. Comparison
// between them is done by MutantProteinConformationComparison, not here.
//

#include "ConformationResult.h"
#include "ProteinConformation.h"
#include <string>

namespace nmr {

class OrcaShieldingResult : public ConformationResult {
public:
    std::string Name() const override { return "OrcaShieldingResult"; }
    std::vector<std::type_index> Dependencies() const override { return {}; }

    // Factory: parse ORCA NMR output and store tensors on conformation atoms.
    // The NMR output atom ordering must match the conformation's atom ordering
    // (both come from the same tleap-protonated geometry).
    static std::unique_ptr<OrcaShieldingResult> Compute(
        ProteinConformation& conf,
        const std::string& nmr_out_path);

    // Query methods
    const std::string& Source() const { return source_; }
    int ParsedAtomCount() const { return parsed_count_; }

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
    std::string source_;
    int parsed_count_ = 0;
};

}  // namespace nmr

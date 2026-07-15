#pragma once
//
// MopacMcConnellResult: compatibility handle for the McConnell BO channel.
//
// The new McConnell model is unified in McConnellResult:
//   fixed channel: sum D(r) Qhat
//   BO channel:    sum Wiberg_bond_order * D(r) Qhat
//
// This result is retained so static and trajectory consumers keep the
// MOPAC-specific BO projection and its sparse/conditional source semantics.
// McConnellResult owns the computation; this calculator reads back the
// ConformationAtom::mopac_mc_* values under their own MOPAC names.
//

#include "ConformationResult.h"
#include "Types.h"

#include <vector>
#include <typeindex>

namespace nmr {

class ProteinConformation;

// NOTE: vestigial. The unified McConnell BO channel uses
// "mcconnell_bond_anisotropy_cutoff"; the old
// "mopac_mcconnell_bond_anisotropy_cutoff" config key is intentionally
// absent. This constant is referenced only in comments/docs, never in
// the math path.
constexpr double MOPAC_MCCONNELL_CUTOFF_A = 10.0;


class MopacMcConnellResult : public ConformationResult {
public:
    std::string Name() const override { return "MopacMcConnellResult"; }

    std::vector<std::type_index> Dependencies() const override;

    // Factory: attach compatibility result after McConnellResult has populated
    // ConformationAtom::mopac_mc_* legacy projections from the BO channel.
    static std::unique_ptr<MopacMcConnellResult> Compute(
        ProteinConformation& conf);

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    const ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr

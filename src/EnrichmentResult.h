#pragma once
//
// EnrichmentResult: assign AtomRole and categorical booleans to every atom.
//
// Dependencies: none (uses bond graph and backbone index cache from Protein).
//
// Assignment logic uses TYPED properties only:
//   - Backbone index cache (res.N, res.CA, etc.) for backbone roles
//   - Ring atom_indices membership for aromatic roles
//   - Bond graph traversal for hydrogen parent classification
//   - Element enum for sidechain roles
//
// No string comparisons after the PDB loading boundary.
//

#include "ConformationResult.h"
#include "ProteinConformation.h"
#include <map>
#include <vector>
#include <typeindex>

namespace nmr {

class EnrichmentResult : public ConformationResult {
public:
    std::string Name() const override { return "EnrichmentResult"; }
    std::vector<std::type_index> Dependencies() const override { return {}; }

    // Factory: assign roles and categorical booleans to all atoms
    static std::unique_ptr<EnrichmentResult> Compute(ProteinConformation& conf);

    // Pre-built collection: atoms grouped by role
    const std::map<AtomRole, std::vector<size_t>>& AtomsByRole() const {
        return atoms_by_role_;
    }

    // Diagnostics
    int UnknownCount() const { return unknown_count_; }

private:
    ProteinConformation* conf_ = nullptr;
    std::map<AtomRole, std::vector<size_t>> atoms_by_role_;
    int unknown_count_ = 0;
};

}  // namespace nmr

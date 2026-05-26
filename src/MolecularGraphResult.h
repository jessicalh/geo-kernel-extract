#pragma once
// BFS-based through-bond features from Protein's bond graph.

#include "ConformationResult.h"
#include "ProteinConformation.h"
#include "SpatialIndexResult.h"
#include <typeindex>

namespace nmr {

// BFS decay: exp(-d / DECAY_LENGTH) where d is graph distance in bonds
constexpr double BFS_DECAY_LENGTH = 4.0;

class MolecularGraphResult : public ConformationResult {
public:
    std::string Name() const override { return "MolecularGraphResult"; }

    std::vector<std::type_index> Dependencies() const override {
        return { typeid(SpatialIndexResult) };
    }

    static std::unique_ptr<MolecularGraphResult> Compute(ProteinConformation& conf);

private:
    ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr

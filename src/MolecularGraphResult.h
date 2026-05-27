#pragma once
//
// MolecularGraphResult: BFS-based through-bond features from the bond graph.
//
// Dependencies: SpatialIndexResult.
//
// Uses the bond graph (Protein's bond list) for BFS traversal.
// Per atom: graph_dist_ring, graph_dist_N, graph_dist_O,
// eneg_sum_1, eneg_sum_2, n_pi_bonds_3, is_conjugated,
// bfs_to_nearest_ring_atom, bfs_decay.
//
// The BFS uses the TYPED bond graph from Protein -- atom.bond_indices
// gives bonds, each bond gives the other atom. No string work.
//

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

    // Factory: compute BFS distances and through-bond features
    static std::unique_ptr<MolecularGraphResult> Compute(ProteinConformation& conf);

private:
    ProteinConformation* conf_ = nullptr;
};

}  // namespace nmr

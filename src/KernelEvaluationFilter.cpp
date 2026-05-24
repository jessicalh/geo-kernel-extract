#include "KernelEvaluationFilter.h"
#include "Protein.h"

namespace nmr {


// ============================================================================
// RingBondedExclusionFilter constructor.
//
// Walks: Protein → RingAt(ri) → atom_indices → AtomAt(vi) → bond_indices
//        → BondAt(bi) → atom_index_a / atom_index_b
//
// For each ring, the exclusion set contains:
//   - Every vertex atom of the ring
//   - Every atom covalently bonded to any vertex
//
// This is the same topology walk that DispersionResult::BondedToVertices
// performs. Centralised here so all ring calculators share one path
// through the protein's bond graph.
// ============================================================================

RingBondedExclusionFilter::RingBondedExclusionFilter(const Protein& protein) {
    size_t n_rings = protein.RingCount();
    ring_bonded_.resize(n_rings);

    for (size_t ri = 0; ri < n_rings; ++ri) {
        const Ring& ring = protein.RingAt(ri);
        auto& bonded = ring_bonded_[ri];

        for (size_t vi : ring.atom_indices) {
            bonded.insert(vi);  // the vertex itself

            const auto& atom = protein.AtomAt(vi);
            for (size_t bi : atom.bond_indices) {
                const auto& bond = protein.BondAt(bi);
                bonded.insert(bond.atom_index_a);
                bonded.insert(bond.atom_index_b);
            }
        }
    }
}


}  // namespace nmr

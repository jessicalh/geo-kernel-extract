#include "MolecularGraphResult.h"
#include "Protein.h"
#include "OperationLog.h"
#include <queue>
#include <set>
#include <cmath>
#include <limits>

namespace nmr {

// ============================================================================
// Multi-source BFS: distance from each atom to the nearest atom in a source set.
// Returns -1 for atoms not reachable from any source.
// ============================================================================

static std::vector<int> BfsFromSet(
        const std::set<size_t>& sources,
        const Protein& protein) {
    std::vector<int> dist(protein.AtomCount(), -1);
    std::queue<size_t> queue;

    for (size_t s : sources) {
        dist[s] = 0;
        queue.push(s);
    }

    while (!queue.empty()) {
        size_t current = queue.front();
        queue.pop();
        int d = dist[current];

        for (size_t bi : protein.AtomAt(current).bond_indices) {
            const Bond& bond = protein.BondAt(bi);
            size_t other = (bond.atom_index_a == current)
                ? bond.atom_index_b : bond.atom_index_a;
            if (dist[other] < 0) {
                dist[other] = d + 1;
                queue.push(other);
            }
        }
    }
    return dist;
}


std::unique_ptr<MolecularGraphResult> MolecularGraphResult::Compute(
        ProteinConformation& conf) {

    OperationLog::Scope scope("MolecularGraphResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()));

    auto result = std::make_unique<MolecularGraphResult>();
    result->conf_ = &conf;

    const Protein& protein = conf.ProteinRef();

    // ---------------------------------------------------------------
    // Build source sets for multi-source BFS
    // ---------------------------------------------------------------

    // Ring atom set (from typed ring atom_indices)
    std::set<size_t> ring_atoms;
    for (size_t ri = 0; ri < protein.RingCount(); ++ri) {
        for (size_t ai : protein.RingAt(ri).atom_indices) {
            ring_atoms.insert(ai);
        }
    }

    // Nitrogen atom set
    std::set<size_t> nitrogen_atoms;
    for (size_t ai = 0; ai < protein.AtomCount(); ++ai) {
        if (protein.AtomAt(ai).element == Element::N) {
            nitrogen_atoms.insert(ai);
        }
    }

    // Oxygen atom set
    std::set<size_t> oxygen_atoms;
    for (size_t ai = 0; ai < protein.AtomCount(); ++ai) {
        if (protein.AtomAt(ai).element == Element::O) {
            oxygen_atoms.insert(ai);
        }
    }

    // ---------------------------------------------------------------
    // Multi-source BFS: distance to nearest ring atom, N, O
    // ---------------------------------------------------------------
    std::vector<int> dist_ring = BfsFromSet(ring_atoms, protein);
    std::vector<int> dist_N = BfsFromSet(nitrogen_atoms, protein);
    std::vector<int> dist_O = BfsFromSet(oxygen_atoms, protein);

    // ---------------------------------------------------------------
    // Build pi-bond set: bonds with aromatic or double bond character
    // ---------------------------------------------------------------
    std::set<size_t> pi_bond_indices;
    for (size_t bi = 0; bi < protein.BondCount(); ++bi) {
        const Bond& bond = protein.BondAt(bi);
        if (bond.order == BondOrder::Double ||
            bond.order == BondOrder::Aromatic ||
            bond.order == BondOrder::Peptide) {
            pi_bond_indices.insert(bi);
        }
    }

    // ---------------------------------------------------------------
    // Per-atom features
    // ---------------------------------------------------------------
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        auto& ca = conf.MutableAtomAt(ai);
        const Atom& identity = protein.AtomAt(ai);

        // Graph distances from BFS
        ca.graph_dist_ring = dist_ring[ai];
        ca.graph_dist_N = dist_N[ai];
        ca.graph_dist_O = dist_O[ai];

        // bfs_to_nearest_ring_atom: same as graph_dist_ring
        ca.bfs_to_nearest_ring_atom = dist_ring[ai];

        // bfs_decay: exp(-d / DECAY_LENGTH), 0 if unreachable
        if (dist_ring[ai] >= 0) {
            ca.bfs_decay = std::exp(-static_cast<double>(dist_ring[ai]) / BFS_DECAY_LENGTH);
        } else {
            ca.bfs_decay = 0.0;
        }

        // eneg_sum_1: sum of electronegativity of atoms bonded to this atom
        // (1-bond neighbourhood)
        double eneg1 = 0.0;
        for (size_t bi : identity.bond_indices) {
            const Bond& bond = protein.BondAt(bi);
            size_t other = (bond.atom_index_a == ai)
                ? bond.atom_index_b : bond.atom_index_a;
            eneg1 += protein.AtomAt(other).Electronegativity();
        }
        ca.eneg_sum_1 = eneg1;

        // eneg_sum_2: sum of electronegativity of atoms 2 bonds away
        double eneg2 = 0.0;
        for (size_t bi : identity.bond_indices) {
            const Bond& bond = protein.BondAt(bi);
            size_t nb1 = (bond.atom_index_a == ai)
                ? bond.atom_index_b : bond.atom_index_a;
            for (size_t bi2 : protein.AtomAt(nb1).bond_indices) {
                const Bond& bond2 = protein.BondAt(bi2);
                size_t nb2 = (bond2.atom_index_a == nb1)
                    ? bond2.atom_index_b : bond2.atom_index_a;
                if (nb2 != ai) {  // skip going back to self
                    eneg2 += protein.AtomAt(nb2).Electronegativity();
                }
            }
        }
        ca.eneg_sum_2 = eneg2;

        // n_pi_bonds_3: number of pi bonds within 3 bonds
        // Uses a small BFS limited to depth 3 from this atom, counting pi bonds
        int pi_count = 0;
        {
            std::vector<int> local_dist(protein.AtomCount(), -1);
            std::queue<size_t> q;
            local_dist[ai] = 0;
            q.push(ai);
            while (!q.empty()) {
                size_t cur = q.front();
                q.pop();
                int d = local_dist[cur];
                if (d >= 3) continue;

                for (size_t bi : protein.AtomAt(cur).bond_indices) {
                    const Bond& b = protein.BondAt(bi);
                    size_t other = (b.atom_index_a == cur)
                        ? b.atom_index_b : b.atom_index_a;
                    if (local_dist[other] < 0) {
                        local_dist[other] = d + 1;
                        q.push(other);
                    }
                    // Count pi bond if it hasn't been counted yet
                    // (check both directions to avoid double counting)
                    if (pi_bond_indices.count(bi) > 0 && local_dist[cur] == d) {
                        // Only count if this is the first time we traverse this bond
                        if (local_dist[other] == d + 1) {
                            pi_count++;
                        }
                    }
                }
            }
        }
        ca.n_pi_bonds_3 = pi_count;

        // is_conjugated: part of a chain of alternating single/double bonds.
        // Simple heuristic: atom has at least one pi bond and at least 2 bonds total
        bool has_pi = false;
        for (size_t bi : identity.bond_indices) {
            if (pi_bond_indices.count(bi) > 0) {
                has_pi = true;
                break;
            }
        }
        ca.is_conjugated = has_pi && identity.bond_indices.size() >= 2;
    }

    // Diagnostics
    int ring_atoms_at_zero = 0;
    int unreachable = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        if (conf.AtomAt(ai).graph_dist_ring == 0) ring_atoms_at_zero++;
        if (conf.AtomAt(ai).graph_dist_ring < 0) unreachable++;
    }

    OperationLog::Info(LogResultAttach, "MolecularGraphResult::Compute",
        "ring_atoms=" + std::to_string(ring_atoms_at_zero) +
        " unreachable_from_ring=" + std::to_string(unreachable));

    return result;
}

}  // namespace nmr

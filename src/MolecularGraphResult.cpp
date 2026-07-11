#include "MolecularGraphResult.h"
#include "Protein.h"
#include "OperationLog.h"
#include "NpyWriter.h"
#include <queue>
#include <set>
#include <cmath>
#include <limits>

namespace nmr {

// ============================================================================
// Multi-source BFS: distance from each atom to the nearest atom in a source set.
// Returns -1 for atoms not reachable from any source.
// ============================================================================

struct BfsResult {
    std::vector<int> distance;
    std::vector<int> nearest_source;
};

static BfsResult BfsFromSet(
        const std::set<size_t>& sources,
        const Protein& protein) {
    BfsResult result;
    result.distance.assign(protein.AtomCount(), -1);
    result.nearest_source.assign(protein.AtomCount(), -1);
    std::queue<size_t> queue;

    for (size_t s : sources) {
        result.distance[s] = 0;
        result.nearest_source[s] = static_cast<int>(s);
        queue.push(s);
    }

    while (!queue.empty()) {
        size_t current = queue.front();
        queue.pop();
        int d = result.distance[current];

        for (size_t bi : protein.AtomAt(current).bond_indices) {
            const Bond& bond = protein.BondAt(bi);
            size_t other = (bond.atom_index_a == current)
                ? bond.atom_index_b : bond.atom_index_a;
            if (result.distance[other] < 0) {
                result.distance[other] = d + 1;
                result.nearest_source[other] = result.nearest_source[current];
                queue.push(other);
            }
        }
    }
    return result;
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

    // Aromatic ring atom set (from typed ring atom_indices)
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
    // Multi-source BFS: distance to nearest aromatic ring atom, N, O
    // ---------------------------------------------------------------
    const BfsResult ring_bfs = BfsFromSet(ring_atoms, protein);
    const BfsResult nitrogen_bfs = BfsFromSet(nitrogen_atoms, protein);
    const BfsResult oxygen_bfs = BfsFromSet(oxygen_atoms, protein);

    // ---------------------------------------------------------------
    // Build pi-like bond set: aromatic, double, and peptide bonds.
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
        ca.graph_dist_ring = ring_bfs.distance[ai];
        ca.graph_dist_N = nitrogen_bfs.distance[ai];
        ca.graph_dist_O = oxygen_bfs.distance[ai];

        // bfs_to_nearest_ring_atom: same as graph_dist_ring.
        ca.bfs_to_nearest_ring_atom = ring_bfs.distance[ai];
        ca.nearest_ring_atom_index = ring_bfs.nearest_source[ai];

        // bfs_decay: exp(-d / DECAY_LENGTH), 0 if unreachable
        if (ring_bfs.distance[ai] >= 0) {
            ca.bfs_decay = std::exp(
                -static_cast<double>(ring_bfs.distance[ai]) /
                BFS_DECAY_LENGTH);
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

        // n_pi_bonds_3: pi-like BFS tree edges within depth 3.
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
                    // Count only forward frontier edges, so the reverse
                    // traversal does not count the same bond again.
                    if (pi_bond_indices.count(bi) > 0 && local_dist[cur] == d) {
                        if (local_dist[other] == d + 1) {
                            pi_count++;
                        }
                    }
                }
            }
        }
        ca.n_pi_bonds_3 = pi_count;

        // is_conjugated: local pi-like bond heuristic, not a full
        // alternating-bond path search.
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

int MolecularGraphResult::WriteFeatures(const ProteinConformation& conf,
                                        const std::string& output_dir) const {
    const size_t N = conf.AtomCount();
    std::vector<int32_t> integer_data(N * 6, 0);
    std::vector<double> float_data(N * 3, 0.0);
    std::vector<double> compatibility(N * 9, 0.0);
    for (size_t i = 0; i < N; ++i) {
        const ConformationAtom& atom = conf.AtomAt(i);
        integer_data[i*6 + 0] = atom.graph_dist_ring;
        integer_data[i*6 + 1] = atom.graph_dist_N;
        integer_data[i*6 + 2] = atom.graph_dist_O;
        integer_data[i*6 + 3] = atom.n_pi_bonds_3;
        integer_data[i*6 + 4] = atom.is_conjugated ? 1 : 0;
        integer_data[i*6 + 5] = atom.nearest_ring_atom_index;
        float_data[i*3 + 0] = atom.eneg_sum_1;
        float_data[i*3 + 1] = atom.eneg_sum_2;
        float_data[i*3 + 2] = atom.bfs_decay;

        compatibility[i*9 + 0] = atom.graph_dist_ring;
        compatibility[i*9 + 1] = atom.graph_dist_N;
        compatibility[i*9 + 2] = atom.graph_dist_O;
        compatibility[i*9 + 3] = atom.eneg_sum_1;
        compatibility[i*9 + 4] = atom.eneg_sum_2;
        compatibility[i*9 + 5] = atom.n_pi_bonds_3;
        compatibility[i*9 + 6] = atom.is_conjugated ? 1.0 : 0.0;
        compatibility[i*9 + 7] = atom.nearest_ring_atom_index;
        compatibility[i*9 + 8] = atom.bfs_decay;
    }

    int written = 0;
    if (NpyWriter::WriteInt32(output_dir + "/molecular_graph_int.npy",
                              integer_data.data(), N, 6)) ++written;
    if (NpyWriter::WriteFloat64(output_dir + "/molecular_graph_float.npy",
                                float_data.data(), N, 3)) ++written;
    if (NpyWriter::WriteFloat64(output_dir + "/molecular_graph.npy",
                                compatibility.data(), N, 9)) ++written;
    return written;
}

}  // namespace nmr

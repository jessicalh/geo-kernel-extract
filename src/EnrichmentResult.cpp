#include "EnrichmentResult.h"
#include "Protein.h"
#include "OperationLog.h"
#include <set>

namespace nmr {

std::unique_ptr<EnrichmentResult> EnrichmentResult::Compute(
        ProteinConformation& conf) {

    OperationLog::Scope scope("EnrichmentResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()));

    auto result = std::make_unique<EnrichmentResult>();
    result->conf_ = &conf;

    const Protein& protein = conf.ProteinRef();

    // ---------------------------------------------------------------
    // Pre-build sets used for role assignment (all typed, no strings)
    // ---------------------------------------------------------------

    // Set of atom indices that are members of aromatic rings
    std::set<size_t> aromatic_atom_set;
    for (size_t ri = 0; ri < protein.RingCount(); ++ri) {
        for (size_t ai : protein.RingAt(ri).atom_indices) {
            aromatic_atom_set.insert(ai);
        }
    }

    // Per-atom backbone membership using backbone index cache
    // An atom is backbone if it matches a backbone slot of its own residue
    std::vector<bool> is_bb(protein.AtomCount(), false);
    for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
        const Residue& res = protein.ResidueAt(ri);
        if (res.N  != Residue::NONE) is_bb[res.N]  = true;
        if (res.CA != Residue::NONE) is_bb[res.CA] = true;
        if (res.C  != Residue::NONE) is_bb[res.C]  = true;
        if (res.O  != Residue::NONE) is_bb[res.O]  = true;
        if (res.H  != Residue::NONE) is_bb[res.H]  = true;
        if (res.HA != Residue::NONE) is_bb[res.HA] = true;
    }

    // ---------------------------------------------------------------
    // Assign roles to every atom
    // ---------------------------------------------------------------

    int unknown_count = 0;

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const Atom& identity = protein.AtomAt(ai);
        const Residue& res = protein.ResidueAt(identity.residue_index);
        auto& ca = conf.MutableAtomAt(ai);

        AtomRole role = AtomRole::Unknown;

        // --- Heavy atom roles ---
        if (identity.element != Element::H) {
            // Check backbone identity via backbone index cache
            if (ai == res.N) {
                role = AtomRole::BackboneN;
            } else if (ai == res.CA) {
                role = AtomRole::BackboneCA;
            } else if (ai == res.C) {
                role = AtomRole::BackboneC;
            } else if (ai == res.O) {
                role = AtomRole::BackboneO;
            }
            // Check aromatic ring membership
            else if (aromatic_atom_set.count(ai) > 0) {
                if (identity.element == Element::N) {
                    role = AtomRole::AromaticN;
                } else {
                    // C in aromatic ring
                    role = AtomRole::AromaticC;
                }
            }
            // Sidechain by element
            else {
                switch (identity.element) {
                    case Element::C: role = AtomRole::SidechainC; break;
                    case Element::N: role = AtomRole::SidechainN; break;
                    case Element::O: role = AtomRole::SidechainO; break;
                    case Element::S: role = AtomRole::SidechainS; break;
                    default: role = AtomRole::Unknown; break;
                }
            }
        }
        // --- Hydrogen roles ---
        else {
            // AmideH: atom index matches residue backbone H cache
            if (ai == res.H) {
                role = AtomRole::AmideH;
            }
            // AlphaH: atom index matches residue HA cache
            else if (ai == res.HA) {
                role = AtomRole::AlphaH;
            }
            // For other H, check parent atom
            else if (identity.parent_atom_index != SIZE_MAX) {
                size_t parent = identity.parent_atom_index;
                const Atom& parent_identity = protein.AtomAt(parent);

                // AromaticH: parent is in the aromatic atom set
                if (aromatic_atom_set.count(parent) > 0) {
                    role = AtomRole::AromaticH;
                }
                // MethylH: parent is sp3 C with exactly 3 H bonds
                else if (parent_identity.element == Element::C) {
                    // Count H bonded to parent
                    int h_count = 0;
                    for (size_t bi : parent_identity.bond_indices) {
                        const Bond& bond = protein.BondAt(bi);
                        size_t other = (bond.atom_index_a == parent)
                            ? bond.atom_index_b : bond.atom_index_a;
                        if (protein.AtomAt(other).element == Element::H) {
                            h_count++;
                        }
                    }
                    if (h_count == 3) {
                        role = AtomRole::MethylH;
                    } else {
                        role = AtomRole::OtherH;
                    }
                }
                // HydroxylH: parent is O bonded to C (SER OG, THR OG1, TYR OH)
                else if (parent_identity.element == Element::O) {
                    bool parent_bonded_to_C = false;
                    for (size_t bi : parent_identity.bond_indices) {
                        const Bond& bond = protein.BondAt(bi);
                        size_t other = (bond.atom_index_a == parent)
                            ? bond.atom_index_b : bond.atom_index_a;
                        if (protein.AtomAt(other).element == Element::C) {
                            parent_bonded_to_C = true;
                            break;
                        }
                    }
                    if (parent_bonded_to_C) {
                        role = AtomRole::HydroxylH;
                    } else {
                        role = AtomRole::OtherH;
                    }
                }
                else {
                    role = AtomRole::OtherH;
                }
            }
            else {
                // No parent assigned -- fallback
                role = AtomRole::OtherH;
            }
        }

        ca.role = role;

        // --- Categorical booleans ---
        ca.is_backbone = is_bb[ai];
        ca.is_amide_H = (role == AtomRole::AmideH);
        ca.is_alpha_H = (role == AtomRole::AlphaH);
        ca.is_methyl = (role == AtomRole::MethylH);
        ca.is_aromatic_H = (role == AtomRole::AromaticH);
        ca.is_on_aromatic_residue = res.IsAromatic();

        // H-bond donor: hydrogen bonded to N or O
        if (identity.element == Element::H &&
            identity.parent_atom_index != SIZE_MAX) {
            Element parent_elem = protein.AtomAt(identity.parent_atom_index).element;
            ca.is_hbond_donor = (parent_elem == Element::N || parent_elem == Element::O);
        } else {
            ca.is_hbond_donor = false;
        }

        // H-bond acceptor: N or O with lone pair (heavy atoms)
        ca.is_hbond_acceptor = (identity.element == Element::N ||
                                identity.element == Element::O);

        // --- Hybridisation (heavy atoms only) ---
        // Approximated from bond environment and role:
        //   sp2: aromatic atoms, backbone C (carbonyl), backbone N (peptide)
        //   sp3: all other heavy C/N/O atoms (CA, CB, sidechain CH2, etc.)
        //   Unassigned: H, S, and anything we can't classify
        if (identity.element != Element::H && identity.element != Element::S) {
            if (aromatic_atom_set.count(ai) > 0) {
                ca.hybridisation = Hybridisation::sp2;
            } else if (ai == res.C) {
                ca.hybridisation = Hybridisation::sp2;  // carbonyl C
            } else if (ai == res.N) {
                ca.hybridisation = Hybridisation::sp2;  // peptide N
            } else {
                ca.hybridisation = Hybridisation::sp3;
            }
        }

        // parent_is_sp2: for H atoms, check parent hybridisation.
        // Heavy atoms are processed before their hydrogens in the
        // residue atom list, so parent hybridisation is already set.
        if (identity.element == Element::H &&
            identity.parent_atom_index != SIZE_MAX) {
            size_t parent = identity.parent_atom_index;
            bool parent_sp2 = false;
            if (aromatic_atom_set.count(parent) > 0)
                parent_sp2 = true;
            else if (parent == res.C)
                parent_sp2 = true;
            else if (parent == res.N)
                parent_sp2 = true;
            ca.parent_is_sp2 = parent_sp2;
        } else {
            ca.parent_is_sp2 = false;
        }

        // Track unknown
        if (role == AtomRole::Unknown) unknown_count++;

        // Build atoms_by_role
        result->atoms_by_role_[role].push_back(ai);
    }

    result->unknown_count_ = unknown_count;

    // Log summary
    std::string summary;
    for (const auto& kv : result->atoms_by_role_) {
        if (!summary.empty()) summary += " ";
        summary += std::to_string(static_cast<int>(kv.first)) + ":" +
                   std::to_string(kv.second.size());
    }

    OperationLog::Info(LogResultAttach, "EnrichmentResult::Compute",
        "roles_assigned total=" + std::to_string(conf.AtomCount()) +
        " unknown=" + std::to_string(unknown_count) +
        " [" + summary + "]");

    return result;
}

}  // namespace nmr

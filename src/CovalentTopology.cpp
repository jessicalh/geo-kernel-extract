#include "CovalentTopology.h"
#include "Atom.h"
#include "Residue.h"
#include "Ring.h"
#include "OperationLog.h"

#include <set>
#include <sstream>
#include <cstdio>

#ifdef HAS_OPENBABEL
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/builder.h>
#include <openbabel/obiter.h>
#endif

namespace nmr {


#ifdef HAS_OPENBABEL
static std::vector<std::pair<size_t, size_t>> DetectBondsViaOpenBabel(
        const std::vector<std::unique_ptr<Atom>>& atoms,
        const std::vector<Vec3>& positions) {

    OpenBabel::OBMol mol;
    mol.BeginModify();

    for (size_t i = 0; i < atoms.size(); ++i) {
        OpenBabel::OBAtom* ob_atom = mol.NewAtom();
        ob_atom->SetAtomicNum(AtomicNumberForElement(atoms[i]->element));
        ob_atom->SetVector(positions[i].x(), positions[i].y(), positions[i].z());
    }
    mol.EndModify();

    mol.ConnectTheDots();
    mol.PerceiveBondOrders();

    std::vector<std::pair<size_t, size_t>> pairs;
    for (auto bond_iter = mol.BeginBonds(); bond_iter != mol.EndBonds(); ++bond_iter) {
        OpenBabel::OBBond* ob_bond = *bond_iter;
        size_t a = ob_bond->GetBeginAtomIdx() - 1;
        size_t b = ob_bond->GetEndAtomIdx() - 1;
        if (a < atoms.size() && b < atoms.size())
            pairs.push_back({std::min(a, b), std::max(a, b)});
    }
    return pairs;
}
#endif


std::unique_ptr<CovalentTopology> CovalentTopology::Resolve(
        const std::vector<std::unique_ptr<Atom>>& atoms,
        const std::vector<Residue>& residues,
        const std::vector<Vec3>& positions,
        double bond_tolerance) {

    auto topo = std::make_unique<CovalentTopology>();
    const size_t n_atoms = atoms.size();

    topo->bond_indices_.resize(n_atoms);
    topo->h_parent_.assign(n_atoms, SIZE_MAX);

    std::vector<bool> is_backbone(n_atoms, false);
    for (const auto& res : residues) {
        if (res.N  != Residue::NONE) is_backbone[res.N]  = true;
        if (res.CA != Residue::NONE) is_backbone[res.CA] = true;
        if (res.C  != Residue::NONE) is_backbone[res.C]  = true;
        if (res.O  != Residue::NONE) is_backbone[res.O]  = true;
        if (res.H  != Residue::NONE) is_backbone[res.H]  = true;
        if (res.HA != Residue::NONE) is_backbone[res.HA] = true;
        if (res.CB != Residue::NONE) is_backbone[res.CB] = true;
    }

    std::vector<std::pair<size_t, size_t>> bond_pairs;

#ifdef HAS_OPENBABEL
    bond_pairs = DetectBondsViaOpenBabel(atoms, positions);
    fprintf(stderr, "CovalentTopology::Resolve: OpenBabel detected %zu bonds.\n",
            bond_pairs.size());
#else
    for (size_t i = 0; i < n_atoms; ++i) {
        for (size_t j = i + 1; j < n_atoms; ++j) {
            double dist = (positions[i] - positions[j]).norm();
            double threshold = atoms[i]->CovalentRadius() +
                               atoms[j]->CovalentRadius() + bond_tolerance;
            if (dist > threshold) continue;
            bond_pairs.push_back({i, j});
        }
    }
#endif

    for (const auto& [i, j] : bond_pairs) {
        double dist = (positions[i] - positions[j]).norm();

        Bond bond;
        bond.atom_index_a = i;
        bond.atom_index_b = j;

        const Atom& a = *atoms[i];
        const Atom& b = *atoms[j];

        if (a.element == Element::S && b.element == Element::S) {
            bond.order = BondOrder::Single;
            bond.category = BondCategory::Disulfide;
        }
        else {
            bool a_bb = is_backbone[i];
            bool b_bb = is_backbone[j];

            if (a_bb && b_bb) {
                bool is_peptide_co = false;
                bool is_peptide_cn = false;

                const Residue& res_a = residues[a.residue_index];
                const Residue& res_b = residues[b.residue_index];

                if (a.residue_index == b.residue_index) {
                    if ((i == res_a.C && j == res_a.O) ||
                        (i == res_a.O && j == res_a.C)) {
                        is_peptide_co = true;
                    }
                }

                if (!is_peptide_co && a.residue_index != b.residue_index) {
                    if ((i == res_a.C && j == res_b.N) ||
                        (i == res_a.N && j == res_b.C)) {
                        is_peptide_cn = true;
                    }
                }

                if (is_peptide_co) {
                    bond.order = BondOrder::Double;
                    bond.category = BondCategory::PeptideCO;
                } else if (is_peptide_cn) {
                    bond.order = BondOrder::Peptide;
                    bond.category = BondCategory::PeptideCN;
                } else {
                    bond.order = BondOrder::Single;
                    bond.category = BondCategory::BackboneOther;
                }
            }
            else if (!a_bb && !b_bb) {
                if (((a.element == Element::C && b.element == Element::O) ||
                     (a.element == Element::O && b.element == Element::C)) &&
                    dist < 1.35) {
                    bond.order = BondOrder::Double;
                    bond.category = BondCategory::SidechainCO;
                } else {
                    bond.order = BondOrder::Single;
                    bond.category = BondCategory::SidechainOther;
                }
            }
            else {
                bond.order = BondOrder::Single;
                bond.category = BondCategory::SidechainOther;
            }
        }

        size_t bond_idx = topo->bonds_.size();
        topo->bonds_.push_back(bond);

        topo->bond_indices_[i].push_back(bond_idx);
        topo->bond_indices_[j].push_back(bond_idx);
    }

    for (size_t i = 0; i < n_atoms; ++i) {
        if (atoms[i]->element != Element::H) continue;

        double min_dist = 1e30;
        for (size_t bi : topo->bond_indices_[i]) {
            const Bond& bnd = topo->bonds_[bi];
            size_t other = (bnd.atom_index_a == i) ? bnd.atom_index_b : bnd.atom_index_a;
            if (atoms[other]->element == Element::H) continue;
            double d = (positions[i] - positions[other]).norm();
            if (d < min_dist) {
                min_dist = d;
                topo->h_parent_[i] = other;
            }
        }
    }

    return topo;
}


void CovalentTopology::TagAromaticBonds(
        const std::vector<std::unique_ptr<Ring>>& rings) {

    std::set<size_t> aromatic_atoms;
    for (const auto& ring : rings) {
        for (size_t ai : ring->atom_indices) {
            aromatic_atoms.insert(ai);
        }
    }
    if (aromatic_atoms.empty()) return;

    for (Bond& bond : bonds_) {
        if (bond.category == BondCategory::Disulfide) continue;
        if (aromatic_atoms.count(bond.atom_index_a) > 0 &&
            aromatic_atoms.count(bond.atom_index_b) > 0) {
            bond.order = BondOrder::Aromatic;
            bond.category = BondCategory::Aromatic;
        }
    }
}


std::string CovalentTopology::OverrideDisulfides(
        const std::vector<DisulfidePair>& pairs) {

    // Track which authority pairs got applied so we can detect leftover
    // geometric Disulfide bonds that aren't in the authority list.
    std::set<std::pair<size_t, size_t>> authority_pair_set;
    auto canonical_pair = [](size_t a, size_t b) {
        return std::pair<size_t, size_t>(std::min(a, b), std::max(a, b));
    };

    int n_force_added = 0;
    int n_overrides = 0;

    for (const auto& dp : pairs) {
        const size_t a = dp.atom_index_sg_a;
        const size_t b = dp.atom_index_sg_b;

        if (a >= bond_indices_.size() || b >= bond_indices_.size()) {
            std::ostringstream err;
            err << "CovalentTopology::OverrideDisulfides: atom index "
                << "out of range — got pair (" << a << ", " << b
                << ") with bond_indices_.size()=" << bond_indices_.size();
            return err.str();
        }
        if (a == b) {
            std::ostringstream err;
            err << "CovalentTopology::OverrideDisulfides: degenerate "
                << "DisulfidePair both indices = " << a;
            return err.str();
        }

        authority_pair_set.insert(canonical_pair(a, b));

        size_t found_index = SIZE_MAX;
        for (size_t bi : bond_indices_[a]) {
            const Bond& bnd = bonds_[bi];
            if ((bnd.atom_index_a == a && bnd.atom_index_b == b) ||
                (bnd.atom_index_a == b && bnd.atom_index_b == a)) {
                found_index = bi;
                break;
            }
        }

        if (found_index != SIZE_MAX) {
            Bond& bnd = bonds_[found_index];
            if (bnd.category != BondCategory::Disulfide) {
                ++n_overrides;
            }
            bnd.category = BondCategory::Disulfide;
            bnd.order    = BondOrder::Single;
        } else {
            // Force-added bonds use the same canonical endpoint order
            // as the geometric path.
            Bond bnd;
            bnd.atom_index_a = std::min(a, b);
            bnd.atom_index_b = std::max(a, b);
            bnd.order    = BondOrder::Single;
            bnd.category = BondCategory::Disulfide;
            const size_t new_idx = bonds_.size();
            bonds_.push_back(bnd);
            bond_indices_[a].push_back(new_idx);
            bond_indices_[b].push_back(new_idx);
            ++n_force_added;
        }
    }

    int n_demoted = 0;
    for (Bond& bnd : bonds_) {
        if (bnd.category != BondCategory::Disulfide) continue;
        if (authority_pair_set.count(canonical_pair(
                bnd.atom_index_a, bnd.atom_index_b)) == 0) {
            bnd.category = BondCategory::SidechainOther;
            ++n_demoted;
        }
    }

    if (n_overrides > 0 || n_force_added > 0 || n_demoted > 0) {
        std::ostringstream msg;
        msg << "applied " << pairs.size() << " authority pairs ("
            << n_overrides << " overrides, "
            << n_force_added << " force-added, "
            << n_demoted << " demoted from geometric)";
        if (n_demoted > 0) {
            OperationLog::Warn(
                "CovalentTopology::OverrideDisulfides", msg.str());
        } else {
            OperationLog::Info(LogCalcOther,
                "CovalentTopology::OverrideDisulfides", msg.str());
        }
    }

    return "";
}


}  // namespace nmr

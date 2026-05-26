#include "Protein.h"
#include "AminoAcidType.h"
#include "LegacyAmberTopology.h"
#include "ChargeSource.h"
#include "NamingRegistry.h"
#include "OperationLog.h"
#include "generated/LegacyAmberSemanticTables.h"
#include <algorithm>
#include <map>
#include <set>
#include <vector>
#include <cstdio>
#include <cstdlib>

namespace nmr {

size_t Protein::AddAtom(std::unique_ptr<Atom> atom) {
    size_t idx = atoms_.size();
    atoms_.push_back(std::move(atom));
    return idx;
}

size_t Protein::AddResidue(Residue residue) {
    size_t idx = residues_.size();
    residues_.push_back(std::move(residue));
    return idx;
}

ProteinConformation& Protein::AddConformation(
    std::vector<Vec3> positions,
    std::string description)
{
    auto conf = std::make_unique<ProteinConformation>(
        this, std::move(positions), std::move(description));
    conformations_.push_back(std::move(conf));
    return *conformations_.back();
}

CrystalConformation& Protein::AddCrystalConformation(
    std::vector<Vec3> positions,
    double resolution, double r_factor,
    double temperature, std::string pdb_id)
{
    auto conf = std::make_unique<CrystalConformation>(
        this, std::move(positions), resolution, r_factor,
        temperature, std::move(pdb_id));
    crystal_index_ = conformations_.size();
    conformations_.push_back(std::move(conf));
    return static_cast<CrystalConformation&>(*conformations_.back());
}

PredictionConformation& Protein::AddPrediction(
    std::vector<Vec3> positions,
    std::string method,
    double confidence)
{
    auto conf = std::make_unique<PredictionConformation>(
        this, std::move(positions), std::move(method), confidence);
    prediction_indices_.push_back(conformations_.size());
    conformations_.push_back(std::move(conf));
    return static_cast<PredictionConformation&>(*conformations_.back());
}

MDFrameConformation& Protein::AddMDFrame(
    std::vector<Vec3> positions,
    int walker, double time_ps, double weight,
    double rmsd_nm, double rg_nm)
{
    auto conf = std::make_unique<MDFrameConformation>(
        this, std::move(positions), walker, time_ps, weight, rmsd_nm, rg_nm);
    md_frame_indices_.push_back(conformations_.size());
    conformations_.push_back(std::move(conf));
    return static_cast<MDFrameConformation&>(*conformations_.back());
}

PredictionConformation& Protein::PredictionAt(size_t i) {
    return static_cast<PredictionConformation&>(*conformations_[prediction_indices_[i]]);
}

const PredictionConformation& Protein::PredictionAt(size_t i) const {
    return static_cast<const PredictionConformation&>(*conformations_[prediction_indices_[i]]);
}

MDFrameConformation& Protein::MDFrameAt(size_t i) {
    return static_cast<MDFrameConformation&>(*conformations_[md_frame_indices_[i]]);
}

const MDFrameConformation& Protein::MDFrameAt(size_t i) const {
    return static_cast<const MDFrameConformation&>(*conformations_[md_frame_indices_[i]]);
}

DerivedConformation& Protein::AddDerived(
    const ProteinConformation& /*parent*/,
    std::string description,
    std::vector<Vec3> positions)
{
    auto conf = std::make_unique<DerivedConformation>(
        this, std::move(positions), std::move(description));
    conformations_.push_back(std::move(conf));
    return static_cast<DerivedConformation&>(*conformations_.back());
}

ProteinConformation& Protein::Conformation() {
    if (conformations_.empty()) {
        fprintf(stderr, "FATAL: Protein::Conformation() -- no conformations.\n");
        std::abort();
    }
    return *conformations_[0];
}

const ProteinConformation& Protein::Conformation() const {
    if (conformations_.empty()) {
        fprintf(stderr, "FATAL: Protein::Conformation() -- no conformations.\n");
        std::abort();
    }
    return *conformations_[0];
}

CrystalConformation& Protein::CrystalConf() {
    if (crystal_index_ == SIZE_MAX) {
        fprintf(stderr, "FATAL: Protein::CrystalConf() -- no crystal conformation.\n");
        std::abort();
    }
    return static_cast<CrystalConformation&>(*conformations_[crystal_index_]);
}

const CrystalConformation& Protein::CrystalConf() const {
    if (crystal_index_ == SIZE_MAX) {
        fprintf(stderr, "FATAL: Protein::CrystalConf() -- no crystal conformation.\n");
        std::abort();
    }
    return static_cast<const CrystalConformation&>(*conformations_[crystal_index_]);
}

const ProteinTopology& Protein::TopologyBase() const {
    if (!protein_topology_) {
        fprintf(stderr, "FATAL: Protein::TopologyBase() -- no topology.\n");
        std::abort();
    }
    return *protein_topology_;
}

const LegacyAmberTopology& Protein::LegacyAmber() const {
    return TopologyAs<LegacyAmberTopology>();
}

size_t Protein::BondCount() const {
    return protein_topology_ ? LegacyAmber().BondCount() : 0;
}

const Bond& Protein::BondAt(size_t i) const {
    return LegacyAmber().BondAt(i);
}

const std::vector<Bond>& Protein::Bonds() const {
    return LegacyAmber().BondList();
}

const CovalentTopology& Protein::BondTopology() const {
    return LegacyAmber().Bonds();
}


bool Protein::BackboneConnected(size_t residue_a_idx,
                                 size_t residue_b_idx) const {
    auto next = BackboneSuccessor(residue_a_idx);
    return next.has_value() && *next == residue_b_idx;
}

std::optional<size_t> Protein::BackbonePredecessor(size_t residue_idx) const {
    if (residue_idx >= residues_.size()) return std::nullopt;
    const Residue& r = residues_[residue_idx];
    if (r.N == Residue::NONE) return std::nullopt;

    const LegacyAmberTopology& topo = LegacyAmber();
    std::optional<size_t> first_match;
    for (size_t bond_idx : topo.BondIndicesFor(r.N)) {
        const Bond& bond = topo.BondAt(bond_idx);
        if (!bond.IsPeptideBond()) continue;
        const size_t other = (bond.atom_index_a == r.N)
            ? bond.atom_index_b : bond.atom_index_a;
        if (other >= atoms_.size()) continue;
        const size_t other_res_idx = atoms_[other]->residue_index;
        if (other_res_idx == residue_idx) continue;
        if (other_res_idx >= residues_.size()) continue;
        if (residues_[other_res_idx].C != other) continue;
        if (first_match.has_value()) {
            OperationLog::Warn(
                "Protein::BackbonePredecessor",
                "Multiple PeptideCN bonds off residue " +
                std::to_string(residue_idx) + " backbone N atom; "
                "first match (residue " + std::to_string(*first_match) +
                ") returned. Loader may have mistagged a cross-link.");
            return first_match;
        }
        first_match = other_res_idx;
    }
    return first_match;
}

std::optional<size_t> Protein::BackboneSuccessor(size_t residue_idx) const {
    if (residue_idx >= residues_.size()) return std::nullopt;
    const Residue& r = residues_[residue_idx];
    if (r.C == Residue::NONE) return std::nullopt;

    const LegacyAmberTopology& topo = LegacyAmber();
    std::optional<size_t> first_match;
    for (size_t bond_idx : topo.BondIndicesFor(r.C)) {
        const Bond& bond = topo.BondAt(bond_idx);
        if (!bond.IsPeptideBond()) continue;
        const size_t other = (bond.atom_index_a == r.C)
            ? bond.atom_index_b : bond.atom_index_a;
        if (other >= atoms_.size()) continue;
        const size_t other_res_idx = atoms_[other]->residue_index;
        if (other_res_idx == residue_idx) continue;
        if (other_res_idx >= residues_.size()) continue;
        if (residues_[other_res_idx].N != other) continue;
        if (first_match.has_value()) {
            OperationLog::Warn(
                "Protein::BackboneSuccessor",
                "Multiple PeptideCN bonds off residue " +
                std::to_string(residue_idx) + " backbone C atom; "
                "first match (residue " + std::to_string(*first_match) +
                ") returned. Loader may have mistagged a cross-link.");
            return first_match;
        }
        first_match = other_res_idx;
    }
    return first_match;
}

size_t Protein::RingCount() const {
    return protein_topology_ ? LegacyAmber().AromaticRingCount() : 0;
}

const Ring& Protein::RingAt(size_t i) const {
    return LegacyAmber().AromaticRingAt(i);
}

const std::vector<std::unique_ptr<Ring>>& Protein::Rings() const {
    return LegacyAmber().AromaticRingList();
}

size_t Protein::SaturatedRingCount() const {
    return protein_topology_ ? LegacyAmber().SaturatedRingCount() : 0;
}

const Ring& Protein::SaturatedRingAt(size_t i) const {
    return LegacyAmber().SaturatedRingAt(i);
}

const std::vector<std::unique_ptr<Ring>>& Protein::SaturatedRings() const {
    return LegacyAmber().SaturatedRingList();
}

const ForceFieldChargeTable& Protein::ForceFieldCharges() const {
    if (!force_field_charges_) {
        fprintf(stderr, "FATAL: Protein::ForceFieldCharges() -- no loaded force-field charges.\n");
        std::abort();
    }
    return *force_field_charges_;
}

void Protein::SetForceFieldCharges(
        std::unique_ptr<ForceFieldChargeTable> charges) {
    if (!charges) {
        fprintf(stderr, "FATAL: Protein::SetForceFieldCharges(nullptr).\n");
        std::abort();
    }
    if (charges->AtomCount() != AtomCount()) {
        fprintf(stderr,
            "FATAL: ForceFieldChargeTable atom count %zu != protein atom count %zu.\n",
            charges->AtomCount(), AtomCount());
        std::abort();
    }
    force_field_charges_ = std::move(charges);
}

bool Protein::PrepareForceFieldCharges(
        const ChargeSource& source,
        const ProteinConformation& conf,
        std::string& error_out) {
    auto table = ForceFieldChargeTable::Build(source, *this, conf, error_out);
    if (!table) return false;
    // ForceFieldChargeTable is the sole authority for charge-source identity.
    SetForceFieldCharges(std::move(table));
    return true;
}

void Protein::FinalizeConstruction(const std::vector<Vec3>& positions,
                                    LegacyAmberInvariants invariants,
                                    double bond_tolerance) {
    // The first backbone cache is name-based; the typed substrate pass
    // at the end overwrites it when semantic tables are available.
    ResolveResidueTerminalStates();
    CacheResidueBackboneIndices();
    ResolveProtonationStates(/*bonds=*/nullptr);

    // Loaders without source force-field data pass an empty invariant pack.
    auto bonds = CovalentTopology::Resolve(atoms_, residues_,
                                           positions, bond_tolerance);

    // An explicit readback authority with zero disulfides still demotes
    // geometry-inferred disulfides; absence of authority leaves geometry intact.
    if (invariants.has_disulfide_authority) {
        const std::string err =
            bonds->OverrideDisulfides(invariants.disulfide_pairs);
        if (!err.empty()) {
            std::fprintf(stderr,
                "FATAL: Protein::FinalizeConstruction "
                "OverrideDisulfides: %s\n", err.c_str());
            std::abort();
        }
    }

    // CYS/CYX depends on the final disulfide bond list, so this pass
    // must precede substrate composition.
    ResolveProtonationStates(bonds.get());

    // ComposeAtomSemantic needs Hydrogen parent names already cached on Atom.
    for (size_t i = 0; i < atoms_.size(); ++i) {
        atoms_[i]->bond_indices = bonds->BondIndicesFor(i);
        atoms_[i]->parent_atom_index = bonds->HydrogenParentOf(i);
    }

    // Re-apply naming after variant resolution. Atom does not retain
    // the original naming source, so pass 2 treats current atom names
    // as AmberFf14SBCanonical inputs.
    {
        const auto& applicator = GlobalNamingApplicator();
        for (Residue& res : residues_) {
            if (!res.protonation_state_resolved) continue;
            if (res.protonation_variant_index < 0) continue;
            const AminoAcidType& aatype = res.AminoAcidInfo();
            if (static_cast<size_t>(res.protonation_variant_index)
                    >= aatype.variants.size()) continue;

            std::vector<std::string> input_names;
            std::vector<std::string> parent_names;
            input_names.reserve(res.atom_indices.size());
            parent_names.reserve(res.atom_indices.size());
            for (size_t ai : res.atom_indices) {
                if (ai >= atoms_.size() || !atoms_[ai]) {
                    input_names.push_back("");
                    parent_names.push_back("");
                    continue;
                }
                input_names.push_back(atoms_[ai]->pdb_atom_name);
                const size_t pai = atoms_[ai]->parent_atom_index;
                parent_names.push_back(
                    (pai != SIZE_MAX && pai < atoms_.size() && atoms_[pai])
                        ? atoms_[pai]->pdb_atom_name : std::string{});
            }

            // Structural N-termini use the ff14SB charged default; this
            // code does not infer neutral termini.
            TerminalState terminal_state = TerminalState::Internal;
            switch (res.terminal_state) {
                case ResidueTerminalState::NTerminus:
                case ResidueTerminalState::NAndCTerminus:
                    terminal_state = TerminalState::NtermCharged;
                    break;
                case ResidueTerminalState::CTerminus:
                    terminal_state = TerminalState::CtermDeprotonated;
                    break;
                case ResidueTerminalState::Internal:
                case ResidueTerminalState::Unknown:
                    terminal_state = TerminalState::Internal;
                    break;
            }

            const auto canonical_names = applicator.ApplyResidue(
                input_names,
                parent_names,
                res.type,
                res.protonation_variant_index,
                terminal_state,
                NamingSource::AmberFf14SBCanonical,
                res.sequence_number,
                res.chain_id);

            for (size_t k = 0; k < res.atom_indices.size(); ++k) {
                const size_t ai = res.atom_indices[k];
                if (ai >= atoms_.size() || !atoms_[ai]) continue;
                atoms_[ai]->pdb_atom_name = canonical_names[k];
            }
        }
    }

    // Stub fixtures with empty PDB names produce an empty semantic table.
    std::vector<AtomSemanticTable> atom_semantic =
        ComposeAtomSemantic(atoms_, residues_, *bonds);

    // Empty atom_semantic yields an empty RingTopology.
    auto rings = RingTopology::ConstructFromSubstrate(
        residues_, atom_semantic, *bonds);

    // Aromatic categories are assigned only after ring construction.
    bonds->TagAromaticBonds(rings->Aromatic());

    protein_topology_ = std::make_unique<LegacyAmberTopology>(
        atoms_.size(), residues_.size(), std::move(bonds),
        std::move(invariants), std::move(atom_semantic),
        std::move(rings));

    CacheResidueBackboneIndices_Typed();
}

// Uses residue insertion order within each chain, not the bond graph.
void Protein::ResolveResidueTerminalStates() {
    std::map<std::string, std::vector<size_t>> by_chain;
    for (size_t ri = 0; ri < residues_.size(); ++ri) {
        residues_[ri].terminal_state = ResidueTerminalState::Internal;
        by_chain[residues_[ri].chain_id].push_back(ri);
    }

    for (const auto& kv : by_chain) {
        const auto& indices = kv.second;
        if (indices.empty()) continue;

        if (indices.size() == 1) {
            residues_[indices.front()].terminal_state =
                ResidueTerminalState::NAndCTerminus;
            continue;
        }

        residues_[indices.front()].terminal_state =
            ResidueTerminalState::NTerminus;
        residues_[indices.back()].terminal_state =
            ResidueTerminalState::CTerminus;
    }
}

void Protein::ResolveProtonationStates(const CovalentTopology* bonds) {
    // First pass resolves H-present variants; second pass can also see
    // CYS disulfides through the bond graph.
    std::set<size_t> disulfide_sg;
    if (bonds != nullptr) {
        for (const Bond& bond : bonds->Bonds()) {
            if (bond.category == BondCategory::Disulfide) {
                disulfide_sg.insert(bond.atom_index_a);
                disulfide_sg.insert(bond.atom_index_b);
            }
        }
    }

    for (auto& res : residues_) {
        const AminoAcidType& aatype = res.AminoAcidInfo();
        if (!aatype.is_titratable || aatype.variants.empty()) continue;

        if (res.protonation_variant_index >= 0) {
            res.protonation_state_resolved = true;
        }

        std::map<std::string, size_t> name_to_idx;
        bool has_any_H = false;
        for (size_t ai : res.atom_indices) {
            const Atom& atom = *atoms_[ai];
            name_to_idx[atom.pdb_atom_name] = ai;
            if (atom.element == Element::H) has_any_H = true;
        }

        int variant_idx = res.protonation_variant_index;
        bool resolved = res.protonation_state_resolved;

        if (res.type == AminoAcid::HIS) {
            bool has_HD1 = name_to_idx.find("HD1") != name_to_idx.end();
            bool has_HE2 = name_to_idx.find("HE2") != name_to_idx.end();

            if (has_HD1 && has_HE2) {
                variant_idx = 2;  // HIP
                resolved = true;
            } else if (has_HD1) {
                variant_idx = 0;  // HID
                resolved = true;
            } else if (has_HE2) {
                variant_idx = 1;  // HIE
                resolved = true;
            } else if (has_any_H) {
                resolved = true;
            }
        }
        else if (res.type == AminoAcid::ASP) {
            bool has_HD2 = name_to_idx.find("HD2") != name_to_idx.end();
            if (has_HD2) {
                variant_idx = 0;  // ASH
                resolved = true;
            } else if (has_any_H) {
                resolved = true;  // charged ASP default
            }
        }
        else if (res.type == AminoAcid::GLU) {
            bool has_HE2 = name_to_idx.find("HE2") != name_to_idx.end();
            if (has_HE2) {
                variant_idx = 0;  // GLH
                resolved = true;
            } else if (has_any_H) {
                resolved = true;  // charged GLU default
            }
        }
        else if (res.type == AminoAcid::CYS) {
            auto sg_it = name_to_idx.find("SG");
            if (sg_it != name_to_idx.end() &&
                disulfide_sg.count(sg_it->second) > 0) {
                variant_idx = 0;  // CYX
                resolved = true;
            } else {
                bool has_HG = name_to_idx.find("HG") != name_to_idx.end();
                if (has_HG || has_any_H) {
                    resolved = true;  // free CYS default unless CYX detected
                }
            }
        }
        else if (res.type == AminoAcid::LYS) {
            bool has_HZ1 = name_to_idx.find("HZ1") != name_to_idx.end();
            bool has_HZ2 = name_to_idx.find("HZ2") != name_to_idx.end();
            bool has_HZ3 = name_to_idx.find("HZ3") != name_to_idx.end();

            if (has_HZ1 && has_HZ2 && has_HZ3) {
                resolved = true;  // charged LYS default
            } else if (has_HZ1 || has_HZ2) {
                variant_idx = 0;  // LYN
                resolved = true;
            } else if (has_any_H) {
                resolved = true;
            }
        }
        else if (res.type == AminoAcid::ARG) {
            resolved = true;  // charged ARG default; ARN is not inferred from names
        }
        else if (res.type == AminoAcid::TYR) {
            bool has_HH = name_to_idx.find("HH") != name_to_idx.end();
            if (!has_HH && has_any_H) {
                variant_idx = 0;  // TYM
                resolved = true;
            } else if (has_HH) {
                resolved = true;  // neutral TYR default
            }
        }

        res.protonation_variant_index = variant_idx;
        res.protonation_state_resolved = resolved;
    }
}

// Name-based construction pass used before AtomSemanticTable exists.
void Protein::CacheResidueBackboneIndices() {
    for (auto& res : residues_) {
        for (size_t ai : res.atom_indices) {
            const std::string& name = atoms_[ai]->pdb_atom_name;
            if      (name == "N")   res.N  = ai;
            else if (name == "CA")  res.CA = ai;
            else if (name == "C")   res.C  = ai;
            else if (name == "O")   res.O  = ai;
            else if (name == "H" || name == "HN")  res.H  = ai;
            else if (name == "HA" || name == "HA2") res.HA = ai;
            else if (name == "CB")  res.CB = ai;
        }

        const AminoAcidType& aatype = res.AminoAcidInfo();
        for (int ci = 0; ci < aatype.chi_angle_count && ci < 4; ++ci) {
            const ChiAngleDef& def = aatype.chi_angles[ci];
            for (int j = 0; j < 4; ++j) {
                for (size_t ai : res.atom_indices) {
                    if (atoms_[ai]->pdb_atom_name == def.atoms[j]) {
                        res.chi[ci].a[j] = ai;
                        break;
                    }
                }
            }
        }
    }
}

// Substrate-driven cache. Without AtomSemanticTable, the earlier
// name-based cache remains in place.
void Protein::CacheResidueBackboneIndices_Typed() {
    if (!protein_topology_) return;
    const LegacyAmberTopology& topo = LegacyAmber();
    if (!topo.HasAtomSemantic()) return;

    for (size_t res_idx = 0; res_idx < residues_.size(); ++res_idx) {
        Residue& res = residues_[res_idx];

        // The typed pass may legitimately leave a slot at NONE
        // (Pro res.H, Gly res.CB).
        res.N  = Residue::NONE;
        res.CA = Residue::NONE;
        res.C  = Residue::NONE;
        res.O  = Residue::NONE;
        res.H  = Residue::NONE;
        res.HA = Residue::NONE;
        res.CB = Residue::NONE;

        for (size_t ai : res.atom_indices) {
            if (ai >= topo.AtomSemantic().size()) continue;
            const AtomSemanticTable& sem = topo.SemanticAt(ai);
            switch (sem.backbone_role) {
                case BackboneRole::Nitrogen:        res.N  = ai; break;
                case BackboneRole::AlphaCarbon:     res.CA = ai; break;
                case BackboneRole::CarbonylCarbon:  res.C  = ai; break;
                case BackboneRole::CarbonylOxygen:  res.O  = ai; break;
                case BackboneRole::AmideHydrogen:   res.H  = ai; break;
                case BackboneRole::AlphaHydrogen:   res.HA = ai; break;
                case BackboneRole::None:            break;
            }
        }

        // Gly HA2/HA3 are distinguished by DiastereotopicIndex, not
        // BackboneRole; the cache picks HA2.
        if (res.type == AminoAcid::GLY) {
            AtomMechanicalIdentity gly_ha2_id{
                Element::H, Locant::Alpha, BranchAddress{},
                DiastereotopicIndex::Position2, BackboneRole::None
            };
            std::vector<size_t> matches = topo.ResidueAtomsWithIdentity(
                res_idx, gly_ha2_id, residues_);
            if (!matches.empty()) res.HA = matches[0];
        }

        // Gly has no C-beta identity in the substrate table.
        AtomMechanicalIdentity cb_id{
            Element::C, Locant::Beta, BranchAddress{},
            DiastereotopicIndex::None, BackboneRole::None
        };
        std::vector<size_t> cb_matches = topo.ResidueAtomsWithIdentity(
            res_idx, cb_id, residues_);
        if (!cb_matches.empty()) res.CB = cb_matches[0];

        // Chi definitions use canonical heavy-atom names, so parent_name
        // is unused when converting them to AtomMechanicalIdentity.
        const AminoAcidType& aatype = res.AminoAcidInfo();
        for (int ci = 0; ci < aatype.chi_angle_count && ci < 4; ++ci) {
            const ChiAngleDef& def = aatype.chi_angles[ci];
            for (int j = 0; j < 4; ++j) {
                res.chi[ci].a[j] = Residue::NONE;
                const std::string atom_name = def.atoms[j];
                if (atom_name.empty()) continue;

                const Element elem =
                    topology_generated::ElementFromAtomName(atom_name);
                if (elem == Element::Unknown) continue;

                const AtomMechanicalIdentity id =
                    topology_generated::ComputeAtomMechanicalIdentity(
                        elem, atom_name, /*parent_name=*/"");

                std::vector<size_t> matches =
                    topo.ResidueAtomsWithIdentity(res_idx, id, residues_);
                if (!matches.empty()) res.chi[ci].a[j] = matches[0];
            }
        }
    }
}


}  // namespace nmr

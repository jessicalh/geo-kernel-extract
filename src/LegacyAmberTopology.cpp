#include "LegacyAmberTopology.h"

#include "Atom.h"
#include "AminoAcidType.h"
#include "Bond.h"
#include "generated/LegacyAmberSemanticTables.h"

#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>

namespace nmr {

LegacyAmberTopology::LegacyAmberTopology(
        size_t atom_count,
        size_t residue_count,
        std::unique_ptr<CovalentTopology> bonds,
        LegacyAmberInvariants invariants,
        std::vector<AtomSemanticTable> atom_semantic,
        std::unique_ptr<RingTopology> rings)
    : atom_count_(atom_count)
    , residue_count_(residue_count)
    , bonds_(std::move(bonds))
    , mass_(std::move(invariants.mass))
    , ff_atom_type_index_(std::move(invariants.ff_atom_type_index))
    , ptype_(std::move(invariants.ptype))
    , atomtype_string_(std::move(invariants.atomtype_string))
    , exclusions_(std::move(invariants.exclusions))
    , fudge_qq_(invariants.fudge_qq)
    , rep_pow_(invariants.rep_pow)
    , atnr_(invariants.atnr)
    , num_non_perturbed_(invariants.num_non_perturbed)
    , atom_semantic_(std::move(atom_semantic))
    , rings_(std::move(rings)) {
    if (!bonds_) {
        std::fprintf(stderr,
            "FATAL: LegacyAmberTopology requires a CovalentTopology.\n");
        std::abort();
    }
    if (!atom_semantic_.empty() && atom_semantic_.size() != atom_count_) {
        std::fprintf(stderr,
            "FATAL: LegacyAmberTopology atom_semantic size %zu != atom_count %zu.\n",
            atom_semantic_.size(), atom_count_);
        std::abort();
    }
    // Accessors assume rings_ is non-null even for empty-substrate paths.
    if (!rings_) {
        rings_ = std::make_unique<RingTopology>();
    }
}


const AtomSemanticTable&
LegacyAmberTopology::SemanticAt(size_t atom_index) const {
    if (atom_semantic_.empty()) {
        std::fprintf(stderr,
            "FATAL: LegacyAmberTopology::SemanticAt: atom_semantic not populated. "
            "Caller must gate on HasAtomSemantic() — stub fixtures (atoms with "
            "empty pdb_atom_name) leave the substrate empty. See "
            "spec/plan/bones/topology-encoding-dependencies-2026-05-05.md §H.5.\n");
        std::abort();
    }
    if (atom_index >= atom_semantic_.size()) {
        std::fprintf(stderr,
            "FATAL: LegacyAmberTopology::SemanticAt: atom_index %zu out of "
            "range (size %zu).\n", atom_index, atom_semantic_.size());
        std::abort();
    }
    return atom_semantic_[atom_index];
}


std::vector<size_t>
LegacyAmberTopology::ResidueAtomsWithIdentity(
        size_t residue_index,
        const AtomMechanicalIdentity& identity,
        const std::vector<Residue>& residues) const {
    if (residue_index >= residues.size()) {
        std::fprintf(stderr,
            "FATAL: LegacyAmberTopology::ResidueAtomsWithIdentity: "
            "residue_index %zu out of range (size %zu).\n",
            residue_index, residues.size());
        std::abort();
    }
    std::vector<size_t> matches;
    if (atom_semantic_.empty()) return matches;
    const Residue& res = residues[residue_index];
    for (size_t ai : res.atom_indices) {
        if (ai >= atom_semantic_.size()) continue;
        const AtomSemanticTable& sem = atom_semantic_[ai];
        AtomMechanicalIdentity sem_id{
            sem.element, sem.locant, sem.branch, sem.di_index,
            sem.backbone_role
        };
        if (sem_id == identity) matches.push_back(ai);
    }
    return matches;
}


size_t
LegacyAmberTopology::AtomWithRole(size_t residue_index,
                                  BackboneRole role,
                                  const std::vector<Residue>& residues) const {
    if (residue_index >= residues.size()) {
        std::fprintf(stderr,
            "FATAL: LegacyAmberTopology::AtomWithRole: residue_index %zu "
            "out of range (size %zu).\n", residue_index, residues.size());
        std::abort();
    }
    if (atom_semantic_.empty()) return Residue::NONE;
    const Residue& res = residues[residue_index];
    for (size_t ai : res.atom_indices) {
        if (ai >= atom_semantic_.size()) continue;
        if (atom_semantic_[ai].backbone_role == role) return ai;
    }
    return Residue::NONE;
}


namespace {

nmr::TerminalState NTerminalStateForResidue(const nmr::Residue& res) {
    // AMBER ff14SB default at neutral pH: NtermCharged (NH3+).
    if (res.terminal_state == nmr::ResidueTerminalState::NTerminus ||
        res.terminal_state == nmr::ResidueTerminalState::NAndCTerminus) {
        return nmr::TerminalState::NtermCharged;
    }
    return nmr::TerminalState::Internal;
}

nmr::TerminalState CTerminalStateForResidue(const nmr::Residue& res) {
    // AMBER ff14SB default at neutral pH: CtermDeprotonated (COO-).
    if (res.terminal_state == nmr::ResidueTerminalState::CTerminus ||
        res.terminal_state == nmr::ResidueTerminalState::NAndCTerminus) {
        return nmr::TerminalState::CtermDeprotonated;
    }
    return nmr::TerminalState::Internal;
}

[[noreturn]] void FatalSubstrateMiss(const char* kind,
                                     size_t atom_index,
                                     const Residue& res,
                                     const Atom& atom,
                                     const AtomMechanicalIdentity& ident,
                                     std::uint8_t variant_idx,
                                     int cap_state_int) {
    std::fprintf(stderr,
        "FATAL: ComposeAtomSemantic: %s lookup miss\n"
        "  atom_index = %zu\n"
        "  residue    = %s seq %d chain '%s'\n"
        "  atom_name  = '%s'\n"
        "  identity   = element=%u/locant=%u/branch={%u,%u}/di=%u/role=%u\n"
        "  variant_idx= %u\n"
        "  cap_state  = %d\n"
        "Bundle A canonicalisation + Bundle B post-protonation re-pass should "
        "have caught any naming variance upstream. This is either a substrate "
        "gap or a NamingRegistry rule gap. See spec/plan/"
        "topology-encoding-dependencies-2026-05-05.md §H.5.\n",
        kind, atom_index,
        res.AminoAcidInfo().three_letter_code,
        res.sequence_number,
        res.chain_id.c_str(),
        atom.pdb_atom_name.c_str(),
        static_cast<unsigned>(ident.element),
        static_cast<unsigned>(ident.locant),
        static_cast<unsigned>(ident.branch.outer),
        static_cast<unsigned>(ident.branch.inner),
        static_cast<unsigned>(ident.di_index),
        static_cast<unsigned>(ident.backbone_role),
        static_cast<unsigned>(variant_idx),
        cap_state_int);
    std::abort();
}

}  // namespace


std::vector<AtomSemanticTable>
ComposeAtomSemantic(const std::vector<std::unique_ptr<Atom>>& atoms,
                    const std::vector<Residue>& residues,
                    const CovalentTopology& bonds) {
    namespace gen = nmr::topology_generated;

    // Named atoms on AminoAcid::Unknown residues must fall through to
    // the explicit unsupported-residue check below.
    bool has_real_atom_names = false;
    for (const Residue& res : residues) {
        for (size_t ai : res.atom_indices) {
            if (ai < atoms.size() && !atoms[ai]->pdb_atom_name.empty()) {
                has_real_atom_names = true;
                break;
            }
        }
        if (has_real_atom_names) break;
    }
    if (!has_real_atom_names) return {};

    for (size_t ri = 0; ri < residues.size(); ++ri) {
        const Residue& res = residues[ri];
        if (res.type != AminoAcid::Unknown) continue;
        size_t named_atom_count = 0;
        for (size_t ai : res.atom_indices) {
            if (ai < atoms.size() && !atoms[ai]->pdb_atom_name.empty()) {
                ++named_atom_count;
            }
        }
        if (named_atom_count == 0) continue;

        std::string first_name;
        for (size_t ai : res.atom_indices) {
            if (ai < atoms.size() && !atoms[ai]->pdb_atom_name.empty()) {
                first_name = atoms[ai]->pdb_atom_name;
                break;
            }
        }
        std::fprintf(stderr,
            "FATAL: ComposeAtomSemantic: AminoAcid::Unknown residue at "
            "index %zu (sequence %d, chain '%s') carries %zu named atom(s) "
            "(first: '%s'). The standard-20 substrate has no row for "
            "non-standard residues; composition cannot proceed without a "
            "default row leaking to downstream calculators. The load path "
            "is unsupported for non-standard residues; refuse before "
            "substrate composition. See spec/plan/bones/topology-encoding-"
            "dependencies-2026-05-05.md §H.5 (fail-loud discipline) and "
            "codex-review Finding 4.\n",
            ri,
            res.sequence_number,
            res.chain_id.c_str(),
            named_atom_count,
            first_name.c_str());
        std::abort();
    }

    std::vector<AtomSemanticTable> result;
    result.resize(atoms.size());

    for (const Residue& res : residues) {
        if (res.type == AminoAcid::Unknown) continue;

        const nmr::TerminalState n_state = NTerminalStateForResidue(res);
        const nmr::TerminalState c_state = CTerminalStateForResidue(res);

        const std::uint8_t variant_idx =
            (res.protonation_state_resolved &&
             res.protonation_variant_index >= 0)
            ? static_cast<std::uint8_t>(res.protonation_variant_index)
            : gen::kBaseVariantIdx;

        for (size_t ai : res.atom_indices) {
            if (ai >= atoms.size()) continue;
            const Atom& atom = *atoms[ai];
            const std::string& name = atom.pdb_atom_name;

            // Parent name disambiguates hydrogens such as HG21 vs HG2.
            std::string parent_name;
            if (atom.parent_atom_index != SIZE_MAX &&
                atom.parent_atom_index < atoms.size()) {
                parent_name = atoms[atom.parent_atom_index]->pdb_atom_name;
            }

            const gen::ParsedAtomName parsed =
                gen::ParseAtomName(name, parent_name);

            nmr::AtomMechanicalIdentity ident;
            ident.element       = atom.element;  // typed-Element authority
            ident.locant        = parsed.locant;
            ident.branch        = parsed.branch;
            ident.di_index      = parsed.di_index;
            ident.backbone_role = parsed.backbone_role;

            // Methyl Hs share one generated substrate row, so runtime
            // identity also clears di_index when the parent has >=3 Hs.
            if (atom.element == nmr::Element::H &&
                ident.di_index != nmr::DiastereotopicIndex::None &&
                atom.parent_atom_index != SIZE_MAX &&
                atom.parent_atom_index < atoms.size()) {
                int parent_h_count = 0;
                const Atom& parent = *atoms[atom.parent_atom_index];
                for (size_t bi : parent.bond_indices) {
                    const Bond& bond = bonds.BondAt(bi);
                    const size_t other =
                        (bond.atom_index_a == atom.parent_atom_index)
                            ? bond.atom_index_b : bond.atom_index_a;
                    if (other < atoms.size() &&
                        atoms[other]->element == nmr::Element::H) {
                        ++parent_h_count;
                    }
                }
                if (parent_h_count >= 3) {
                    ident.di_index = nmr::DiastereotopicIndex::None;
                }
            }

            if (parsed.is_cap_only_n || parsed.is_cap_only_c) {
                const nmr::TerminalState cap_state =
                    parsed.is_cap_only_n ? n_state
                    : parsed.is_cap_only_c ? c_state
                    : nmr::TerminalState::Internal;

                const nmr::AtomSemanticTable* cap =
                    gen::LookupCap(cap_state, ident);
                if (cap == nullptr) {
                    FatalSubstrateMiss("LookupCap (cap-only)",
                                       ai, res, atom, ident, variant_idx,
                                       static_cast<int>(cap_state));
                }
                result[ai] = *cap;
                continue;
            }

            const nmr::AtomSemanticTable* base =
                gen::LookupBy(res.type, variant_idx, ident);
            if (base == nullptr) {
                FatalSubstrateMiss("LookupBy (chain)", ai, res, atom,
                                   ident, variant_idx, /*cap_state=*/-1);
            }
            result[ai] = *base;

            // Apply cap deltas field-by-field so terminal backbone rows
            // retain the chain identity fields not overridden by the cap.
            const nmr::BackboneRole bb_role = ident.backbone_role;
            if (gen::IsBackboneCapOverlayAtom(bb_role, n_state, c_state)) {
                const nmr::TerminalState cap_state =
                    (bb_role == nmr::BackboneRole::Nitrogen) ? n_state
                                                              : c_state;
                const nmr::AtomSemanticTable* cap =
                    gen::LookupCap(cap_state, ident);
                if (cap == nullptr) {
                    FatalSubstrateMiss("LookupCap (backbone-cap overlay)",
                                       ai, res, atom, ident, variant_idx,
                                       static_cast<int>(cap_state));
                }
                gen::ApplyCapDelta(result[ai], *cap);
            }
        }
    }

    return result;
}

}  // namespace nmr

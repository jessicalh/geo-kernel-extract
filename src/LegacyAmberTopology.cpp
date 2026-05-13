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
    // Default: empty RingTopology if caller didn't supply one. This
    // keeps stub-fixture call sites simple — they pass `{}` (or omit)
    // and get an empty topology (AromaticCount() == 0,
    // SaturatedCount() == 0). The accessors never need to null-check.
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


// ============================================================================
// ComposeAtomSemantic
//
// Per spec/plan/bones/topology-encoding-dependencies-2026-05-05.md §H.5. The
// composition rule:
//
//   1. Parse the (canonical) atom name + its heavy-atom parent's name
//      ONCE per atom via ParseAtomName (lifted into
//      src/generated/LegacyAmberSemanticTables.h post-Bundle-A). The
//      typed flags on the parse result drive cap-only / backbone /
//      chain dispatch below; the AtomMechanicalIdentity tuple used
//      for LookupBy / LookupCap is built inline from the same parse
//      output. No second ParseAtomName call, no string-taking
//      IsCapOnlyAtomName(name) re-read after the parser boundary
//      (codex-review Finding 3).
//   2. Methyl-H pseudoatom collapse: substrate clears
//      DiastereotopicIndex on methyl Hs (PseudoatomKind::M); the runtime
//      composition replicates this via bond-graph methyl detection
//      (parent has 3+ H neighbours).
//   3. Cap-only atoms (H1/H2/H3 / OXT / HXT) resolve via LookupCap with
//      the residue's terminal state. nullptr is FATAL — the post-
//      protonation re-canonicalisation pass should have caught any
//      naming variance upstream.
//   4. Standard chain atoms resolve via LookupBy(residue, variant_idx,
//      identity). nullptr is FATAL on the same grounds.
//   5. Backbone-cap overlay: terminal-residue backbone N (NTERM) or
//      backbone C/O (CTERM) get cap-delta overlay on top of the chain
//      row via ApplyCapDelta (field-level rule per §H.5; whole-row
//      assignment is forbidden). cap_delta nullptr is acceptable for
//      this overlay path (no cap-side override for this backbone role
//      at this state); chain row stays as-is.
//
// Stub-fixture guard: tests that build proteins from raw element-and-
// position fixtures (e.g. tests/test_coulomb_result.cpp) have empty
// pdb_atom_name on every atom. Detect this and return empty (the
// legitimate "no substrate populated" signal). Calculators that need
// substrate gate on HasAtomSemantic().
//
// The "fail-loudly" policy: post-Bundle-A canonicalisation + post-
// Bundle-B post-protonation re-canonicalisation are designed to make
// every atom's (canonical) name match the substrate. Any remaining
// LookupBy/LookupCap miss is a substrate gap or a naming-rule gap;
// crashing here surfaces it instead of letting a calculator silently
// see a default-constructed AtomSemanticTable. The stash@{0}
// degrade-and-skip predecessor is rejected.
// ============================================================================

namespace {

nmr::TerminalState NTerminalStateForResidue(const nmr::Residue& res) {
    // AMBER ff14SB default at neutral pH: NtermCharged (NH3+).
    // Future: branch on a NTERM_NEUTRAL variant if/when one lands.
    if (res.terminal_state == nmr::ResidueTerminalState::NTerminus ||
        res.terminal_state == nmr::ResidueTerminalState::NAndCTerminus) {
        return nmr::TerminalState::NtermCharged;
    }
    return nmr::TerminalState::Internal;
}

nmr::TerminalState CTerminalStateForResidue(const nmr::Residue& res) {
    // AMBER ff14SB default at neutral pH: CtermDeprotonated (COO-).
    // Future: branch on a CTERM_PROTONATED variant when added.
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

    // Stub-fixture guard. Tests that construct proteins with raw
    // element-and-position fixtures (no PDB names) leave pdb_atom_name
    // empty on every atom. Substrate doesn't apply; return empty vector.
    // Calculators gate on HasAtomSemantic().
    //
    // Codex Finding F4 (2026-05-06): the guard's predicate is "ANY atom
    // anywhere in the protein has a non-empty pdb_atom_name." It does
    // NOT pre-skip Unknown residues. An all-Unknown-but-named-atoms
    // protein hits this guard with has_real_atom_names = true, falls
    // through to the Unknown-residue fail-loud loop below, and aborts
    // there — which is the correct behaviour. The previous formulation
    // ("ignore Unknown residues, look for named atoms only on standard
    // residues") let an all-Unknown-named-atoms protein silently return
    // empty, masking a real chemistry error.
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

    // Fail-loud on AminoAcid::Unknown residues that carry named atoms
    // (codex-review Finding 4). The standard-20 substrate has no row
    // for Unknown residues; previously the per-residue loop simply
    // skipped them and the result vector kept default-constructed
    // entries on those atom slots — which `HasAtomSemantic()` reports
    // as populated, producing default `Element::Unknown` rows that
    // calculators silently consumed. The fail-loud discipline:
    // unsupported residues abort before composition rather than leak
    // default rows downstream. Toy-fixture stub-proteins (atoms with
    // empty pdb_atom_name) are caught above by the stub guard.
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

        // First named atom for diagnostic context.
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

        // Variant index: kBaseVariantIdx unless protonation has been
        // resolved AND the residue carries a typed variant.
        const std::uint8_t variant_idx =
            (res.protonation_state_resolved &&
             res.protonation_variant_index >= 0)
            ? static_cast<std::uint8_t>(res.protonation_variant_index)
            : gen::kBaseVariantIdx;

        for (size_t ai : res.atom_indices) {
            if (ai >= atoms.size()) continue;
            const Atom& atom = *atoms[ai];
            const std::string& name = atom.pdb_atom_name;

            // Parent-atom name for H disambiguation (e.g. HG21 on CG2
            // vs HG2 on CG). parent_atom_index is SIZE_MAX for non-H
            // atoms; ParseAtomName treats empty as "no parent."
            std::string parent_name;
            if (atom.parent_atom_index != SIZE_MAX &&
                atom.parent_atom_index < atoms.size()) {
                parent_name = atoms[atom.parent_atom_index]->pdb_atom_name;
            }

            // PARSE ONCE per atom (codex-review Finding 3 — "parse once,
            // then no string work in composition"). The typed flags on
            // `parsed` drive cap-only / backbone / chain dispatch below;
            // strings DIE at this call's return.
            const gen::ParsedAtomName parsed =
                gen::ParseAtomName(name, parent_name);

            nmr::AtomMechanicalIdentity ident;
            ident.element       = atom.element;  // typed-Element authority
            ident.locant        = parsed.locant;
            ident.branch        = parsed.branch;
            ident.di_index      = parsed.di_index;
            ident.backbone_role = parsed.backbone_role;

            // Methyl-H pseudoatom collapse: substrate clears
            // DiastereotopicIndex on methyl Hs (PseudoatomKind::M) so
            // HE1/HE2/HE3 of MET, HD11/HD12/HD13 of LEU, etc. share a
            // single identity row. Detect via the bond graph (parent
            // has 3+ H neighbours) and clear di_index to match.
            // Mirrors tools/topology/build_semantic_tables.cpp where
            // e.di_index = None when e.pseudoatom.kind == M.
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

            // Dispatch from typed flags (no second ParseAtomName call,
            // no IsCapOnlyAtomName(string) re-read).
            if (parsed.is_cap_only_n || parsed.is_cap_only_c) {
                // Cap-only atoms live only in cap tables. The cap state
                // selecting the table comes from the atom's typed
                // family flag.
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

            // Standard chain atom: must resolve via LookupBy.
            const nmr::AtomSemanticTable* base =
                gen::LookupBy(res.type, variant_idx, ident);
            if (base == nullptr) {
                FatalSubstrateMiss("LookupBy (chain)", ai, res, atom,
                                   ident, variant_idx, /*cap_state=*/-1);
            }
            result[ai] = *base;

            // Backbone-cap overlay (FIELD-LEVEL via ApplyCapDelta).
            //
            // Codex Finding CC1 (2026-05-06): when
            // IsBackboneCapOverlayAtom() selects this branch, the
            // backbone role + terminal state combination is one the
            // substrate is required to populate (Nitrogen at NTERM*,
            // CarbonylCarbon / CarbonylOxygen at CTERM*). A nullptr
            // cap row here is a substrate gap — the chain row would
            // silently keep its internal-residue chemistry while the
            // load-time terminal state says it should have cap-delta
            // overlays applied. Treat nullptr as FATAL with full
            // mechanical-identity context, mirroring the chain-side
            // FatalSubstrateMiss path. Test-gap note: constructing a
            // synthetic fixture that triggers this assertion requires
            // modifying the (auto-generated) substrate table to omit
            // the row; the architectural fix is the value, and the
            // failure becomes obvious any time the substrate is
            // extended for a new variant or terminal state without
            // the matching cap-table population.
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

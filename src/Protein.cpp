#include "Protein.h"
#include "AminoAcidType.h"
#include "LegacyAmberTopology.h"
#include "ChargeSource.h"
#include "NamingRegistry.h"
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


// ============================================================================
// Conformation factory methods
// ============================================================================

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
    // parent is accepted for future use (tracking derivation lineage)
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


// ============================================================================
// Explicit topology / charge contract access
// ============================================================================

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

// ============================================================================
// Ring access — delegated through RingTopology on LegacyAmberTopology.
// Bundle C / Slice B (2026-05-07): mirrors the bond delegation above.
// `RingCount()` and `SaturatedRingCount()` return 0 when no topology
// is attached (pre-FinalizeConstruction state); the `*At(i)` and
// `*s()` accessors require a topology and abort if called without one
// (matches BondAt(i)).
// ============================================================================

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
    // ForceFieldChargeTable already records source_force_field_, kind_,
    // and source_description_ at construction. ProteinBuildContext is
    // not the second authority on charge identity.
    SetForceFieldCharges(std::move(table));
    return true;
}


// ============================================================================
// FinalizeConstruction: every loader must call this after adding all atoms
// and residues. Ensures backbone indices, bonds, and rings are all detected.
// ============================================================================

void Protein::FinalizeConstruction(const std::vector<Vec3>& positions,
                                    LegacyAmberInvariants invariants,
                                    double bond_tolerance) {
    // Layer 2: symbolic topology (no geometry needed). The first
    // CacheResidueBackboneIndices() call here is string-matched and
    // feeds the first ResolveProtonationStates pass (HIS/LYS/etc.
    // detection from explicit H presence). The typed
    // CacheResidueBackboneIndices_Typed() pass at the end of this
    // function overwrites the cache from the substrate.
    ResolveResidueTerminalStates();
    CacheResidueBackboneIndices();
    ResolveProtonationStates(/*bonds=*/nullptr);

    // Layer 3: geometric topology + FF-numerical invariants. The
    // value-pack is moved into the topology's plain fields and goes out
    // of scope after this call. Loaders without source FF data pass {}.
    //
    // Bundle C / Slice B (2026-05-07): rings are no longer an input
    // here. The aromatic-bond categorisation moved to a post-Resolve
    // overlay (CovalentTopology::TagAromaticBonds), called below
    // after substrate-driven ring construction.
    auto bonds = CovalentTopology::Resolve(atoms_, residues_,
                                           positions, bond_tolerance);

    // Apply pdb2gmx's authoritative disulfide pairing (TPR bonded list
    // + rtp comments). Geometric SG-SG inference at this point in
    // CovalentTopology agrees with chemistry on standard fixtures, but
    // the override establishes authority direction: GROMACS decided,
    // we read.
    //
    // Gated on has_disulfide_authority, NOT on !disulfide_pairs.empty():
    // an authority that says "zero disulfides" is meaningful — any
    // geometric Disulfide tag in CovalentTopology must be demoted, not
    // preserved. Falsy authority (PDB load) leaves geometric inference
    // as source of truth.
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

    // Second protonation pass against the now-finalized covalent
    // topology (CYS -> CYX from the disulfide bond list). Run BEFORE
    // LegacyAmberTopology construction so the substrate-composition
    // step sees the final variant_index for every residue.
    ResolveProtonationStates(bonds.get());

    // Copy connectivity onto Atom for convenient calculator access.
    // Done before substrate composition because ComposeAtomSemantic
    // uses atom.parent_atom_index for the H-atom parent-name lookup
    // that ParseAtomName needs.
    for (size_t i = 0; i < atoms_.size(); ++i) {
        atoms_[i]->bond_indices = bonds->BondIndicesFor(i);
        atoms_[i]->parent_atom_index = bonds->HydrogenParentOf(i);
    }

    // Post-protonation second applicator pass: now that variant_index
    // is resolved per residue (LYS-labelled-LYN, HIS variants, etc.),
    // walk every resolved-variant residue and rewrite atom names with
    // the now-known variant_index in NamingContext. The applicator's
    // sibling-aware predicates make this pass IDEMPOTENT on residues
    // whose Pass-1 names are already canonical for the resolved variant
    // (canonical LYN siblings {HZ2,HZ3,no HZ1} do NOT match the LYN
    // pre-Markley pattern and the shift rules don't fire; canonical
    // LYS-NH3+ siblings {HZ1,HZ2,HZ3} match LysAmmoniumHzPassThrough
    // which preserves the input). The 1Z9B fleet variance — LYS-labelled
    // residues with LYN chemistry — is captured during Pass 1 already
    // (loader source = CifppPdbInput; siblings {HZ1,HZ2,no HZ3} fire
    // the LYN shift rules); this Pass-2 call is idempotent on those.
    //
    // No source tag is preserved across passes: Pass-1 source was
    // recorded onto Atom::pdb_atom_name as the canonical output; Pass 2
    // operates on canonical strings tagged AmberFf14SBCanonical (the
    // applicator's pass-through branch returns canonical inputs
    // unchanged). Rules that fire on Pass 2 are sibling-aware shift
    // rules that examine the (now-canonical) sibling set.
    {
        const auto& applicator = GlobalNamingApplicator();
        for (Residue& res : residues_) {
            if (!res.protonation_state_resolved) continue;
            if (res.protonation_variant_index < 0) continue;
            const AminoAcidType& aatype = res.AminoAcidInfo();
            if (static_cast<size_t>(res.protonation_variant_index)
                    >= aatype.variants.size()) continue;

            // Snapshot input names + parent names from the current
            // pdb_atom_name across the residue.
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

            // Map ResidueTerminalState -> TerminalState. The applicator
            // accepts the typed SemanticEnums.h TerminalState; the
            // mapping is conservative (NTerminus -> NtermCharged is
            // the AMBER ff14SB default; NTermNeutral is reserved for
            // PROPKA-driven distinctions not yet exposed).
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

    // Compose the per-atom AtomSemanticTable substrate. Empty for stub
    // calculator-physics fixtures (no PDB names); populated otherwise.
    std::vector<AtomSemanticTable> atom_semantic =
        ComposeAtomSemantic(atoms_, residues_, *bonds);

    // Substrate-driven ring construction (Bundle C / Slice B). Reads
    // typed RingPosition slots from the substrate; produces aromatic
    // rings (PHE/TYR/HIS-variants/TRP-{benzene,pyrrole,9}) and
    // saturated rings (Pro pyrrolidine) in canonical cyclic walk
    // order. Stub fixtures (empty atom_semantic) get an empty
    // RingTopology. See spec/plan/ring-investigation-2026-05-06/.
    auto rings = RingTopology::ConstructFromSubstrate(
        residues_, atom_semantic, *bonds);

    // Aromatic-bond tagging overlay: now that rings exist, tag bonds
    // whose both endpoints sit in any aromatic ring. Mirrors the
    // OverrideDisulfides post-Resolve pattern.
    bonds->TagAromaticBonds(rings->Aromatic());

    // Construct the final LegacyAmberTopology with substrate + rings.
    protein_topology_ = std::make_unique<LegacyAmberTopology>(
        atoms_.size(), residues_.size(), std::move(bonds),
        std::move(invariants), std::move(atom_semantic),
        std::move(rings));

    // Typed CacheResidueBackboneIndices: overwrite res.{N, CA, C, O,
    // H, HA, CB} with substrate-driven indices from
    // LegacyAmber().AtomSemantic(). On stub fixtures (no substrate),
    // this is a no-op and the string-matched cache from the first
    // pass stays.
    CacheResidueBackboneIndices_Typed();
}


// ============================================================================
// ResolveResidueTerminalStates
//
// Construction-boundary interpretation of polymer end state. This records
// which residues are structurally first and last within each chain so
// force-field adapters can choose their own terminal templates explicitly.
// ============================================================================

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


// ============================================================================
// ResolveProtonationStates
//
// Construction-boundary string interpretation. This prepares residue
// protonation/variant state before calculators run. ProtonationDetectionResult
// reports this state; it does not perform identity resolution.
// ============================================================================

void Protein::ResolveProtonationStates(const CovalentTopology* bonds) {
    // bonds == nullptr  : first pass. HIS/LYS/TYR/ASP/GLU variants
    //                     from explicit H presence; CYS variant
    //                     deferred (needs the bond graph).
    // bonds != nullptr  : second pass. Adds CYS -> CYX from
    //                     BondCategory::Disulfide entries.
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


// ============================================================================
// DetectAromaticRings: REMOVED Bundle C / Slice B (2026-05-07).
//
// Ring construction is now substrate-driven via
// RingTopology::ConstructFromSubstrate, called from
// FinalizeConstruction after ComposeAtomSemantic. Rings live on
// LegacyAmberTopology (parallel to bonds on CovalentTopology) and
// are reached through the delegating accessors above.
// ============================================================================


// ============================================================================
// CacheResidueBackboneIndices
//
// PDB LOADING BOUNDARY: string comparisons are used here to match PDB atom
// names ("N", "CA", "C", "O", "H", "HA", "CB") against each residue's atoms.
// This populates typed backbone index fields (res.N, res.CA, etc.) which are
// atom INDICES, not strings. After this function, backbone identity is
// determined entirely by these typed indices -- no further string work.
// ============================================================================

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

        // PDB LOADING BOUNDARY (continued): match chi angle atom names
        // against the residue's atoms. After this, chi angles are typed
        // index tuples — no further string work.
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


// ============================================================================
// CacheResidueBackboneIndices_Typed
//
// Substrate-driven backbone-index cache. Reads BackboneRole + Locant +
// DiastereotopicIndex from the per-atom AtomSemanticTable populated by
// ComposeAtomSemantic, and overwrites res.{N, CA, C, O, H, HA, CB}
// with typed indices.
//
// Special cases:
//   - Glycine HA: Gly's HA2/HA3 carry Locant::Alpha + DiastereotopicIndex
//     (BackboneRole stays None per the parser convention; the typed
//     identity already disambiguates them). Pick HA2 (Position2) to
//     match the string-matched cache's prior assignment of res.HA = HA2.
//   - Proline H: PRO's chain table drops the backbone amide H per the
//     substrate dependencies §H.10 (Pro is a secondary amine); no atom
//     has BackboneRole::AmideHydrogen, so res.H stays Residue::NONE.
//   - CB: the substrate carries Locant::Beta + Element::C + branch{0,0}
//     for the canonical CB. Glycine has no CB atom in its table, so
//     res.CB stays Residue::NONE for Gly.
//
// Chi-angle resolver: STAYS string-matched. Audit Hotspot 2; separate
// slice. Calculators that consume chi angles read res.chi[i] which
// still works against the AminoAcidType chi_angles atom-name list.
//
// Stub-fixture path: when LegacyAmber().HasAtomSemantic() is false
// (atoms with empty pdb_atom_name; calculator-physics tests), this
// function is a no-op and the string-matched cache from the first
// pass stays (which is fine — those tests don't carry residues with
// real atom names anyway).
// ============================================================================

void Protein::CacheResidueBackboneIndices_Typed() {
    if (!protein_topology_) return;
    const LegacyAmberTopology& topo = LegacyAmber();
    if (!topo.HasAtomSemantic()) return;

    for (size_t res_idx = 0; res_idx < residues_.size(); ++res_idx) {
        Residue& res = residues_[res_idx];

        // Reset the backbone slots before substrate-driven repopulation.
        // The first-pass cache wrote them from strings; the typed pass
        // is the authoritative version and may legitimately leave a
        // slot at NONE (Pro res.H, Gly res.CB).
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

        // Glycine special case. HA2/HA3 carry Locant::Alpha + di_index;
        // BackboneRole::None. Pick HA2 (Position2) for res.HA.
        if (res.type == AminoAcid::GLY) {
            AtomMechanicalIdentity gly_ha2_id{
                Element::H, Locant::Alpha, BranchAddress{},
                DiastereotopicIndex::Position2, BackboneRole::None
            };
            std::vector<size_t> matches = topo.ResidueAtomsWithIdentity(
                res_idx, gly_ha2_id, residues_);
            if (!matches.empty()) res.HA = matches[0];
        }

        // CB cache: typed identity for the canonical sidechain Cβ
        // (Locant::Beta, Element::C, branch{0,0}, di_index=None).
        // Gly has no Cβ atom in its substrate table, so the lookup
        // returns empty and res.CB stays NONE.
        AtomMechanicalIdentity cb_id{
            Element::C, Locant::Beta, BranchAddress{},
            DiastereotopicIndex::None, BackboneRole::None
        };
        std::vector<size_t> cb_matches = topo.ResidueAtomsWithIdentity(
            res_idx, cb_id, residues_);
        if (!cb_matches.empty()) res.CB = cb_matches[0];

        // Chi-angle resolver: STAYS string-matched (Audit Hotspot 2;
        // separate slice). Re-resolve here so the typed pass produces
        // the same chi indices the string-matched first pass produced.
        const AminoAcidType& aatype = res.AminoAcidInfo();
        for (int ci = 0; ci < aatype.chi_angle_count && ci < 4; ++ci) {
            const ChiAngleDef& def = aatype.chi_angles[ci];
            for (int j = 0; j < 4; ++j) {
                res.chi[ci].a[j] = Residue::NONE;
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


}  // namespace nmr

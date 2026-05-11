// LarsenResidue.cpp -- bond-graph perception of Larsen ProCS15
// tripeptide DFT records.
//
// The bond-graph + K=3 WL algorithm was prototyped in
//   scripts/perceive_larsen_tripeptide.py
// and validated against all 20 DB (tripeptide, frame_type) combinations
// and the AAA Gaussian log (AAA_4_54_nmr.log byte-for-byte match with
// DB row 21755). The Python script remains a useful diagnostic but is
// NOT a normative spec — it still returns canonical atom names from
// `match_piece` (predates the typed-substrate work). The C++ object
// model is the runtime authority: MatchPiece returns canonical NODE
// INDICES into the CanonicalPiece, canonical identities are stamped
// by looking up the generated topology table
// (`StampChainIdentitiesViaTable`), and slot dispatch happens by
// typed BackboneRole/Locant enum — never by name comparison.
//
// Design doc: spec/plan/larsen-residue-design-2026-05-11.md.

#include "LarsenResidue.h"

#include "AminoAcidType.h"
#include "OperationLog.h"
#include "TripeptideDftTable.h"
#include "generated/LegacyAmberSemanticTables.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <stack>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace nmr {

namespace {

using AdjacencySet = std::vector<std::set<int>>;


// ============================================================================
// Bond perception
// ============================================================================

double BondCutoffAngstrom(Element a, Element b) {
    // Order-insensitive lookup. Cutoffs match the Python POC.
    Element lo = (static_cast<int>(a) <= static_cast<int>(b)) ? a : b;
    Element hi = (lo == a) ? b : a;
    if (lo == Element::H || hi == Element::H) {
        if (lo == Element::H && hi == Element::H) return 0.0;
        if (hi == Element::S) return 1.50;
        return 1.30;
    }
    // Heavy-heavy.
    if (lo == Element::S || hi == Element::S) {
        if (lo == Element::S && hi == Element::S) return 2.30;
        return 2.10;
    }
    if (lo == Element::O && hi == Element::O) return 1.60;
    return 1.85;
}


double Distance(const Vec3& a, const Vec3& b) {
    return (a - b).norm();
}


AdjacencySet BuildBondGraph(const std::vector<TripeptideDftAtom>& atoms) {
    const int n = static_cast<int>(atoms.size());
    AdjacencySet adj(n);
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const double cutoff =
                BondCutoffAngstrom(atoms[i].element, atoms[j].element);
            if (cutoff <= 0.0) continue;
            if (Distance(atoms[i].position, atoms[j].position) < cutoff) {
                adj[i].insert(j);
                adj[j].insert(i);
            }
        }
    }
    return adj;
}


// ============================================================================
// Amide identification
// ============================================================================

// Returns (carbonyl_C_idx, amide_N_idx) pairs for all peptide-amide bonds.
// Filters out sidechain primary amides (ASN ND2, GLN NE2) by requiring the
// amide N to have at least one heavy neighbour besides the carbonyl C.
std::vector<std::pair<int, int>> FindPeptideAmides(
        const std::vector<TripeptideDftAtom>& atoms,
        const AdjacencySet& adj) {
    std::vector<std::pair<int, int>> out;
    const int n = static_cast<int>(atoms.size());

    // First, find all amide-shape candidates (sp2 C with C=O and C-N).
    std::vector<std::pair<int, int>> candidates;
    for (int i = 0; i < n; ++i) {
        if (atoms[i].element != Element::C) continue;
        if (adj[i].size() != 3) continue;  // sp2 carbonyl C has 3 neighbours

        bool has_double_o = false;
        for (int j : adj[i]) {
            if (atoms[j].element == Element::O &&
                Distance(atoms[i].position, atoms[j].position) < 1.32) {
                has_double_o = true; break;
            }
        }
        if (!has_double_o) continue;

        for (int j : adj[i]) {
            if (atoms[j].element != Element::N) continue;
            const double d = Distance(atoms[i].position, atoms[j].position);
            if (d > 1.20 && d < 1.50) {
                candidates.emplace_back(i, j);
            }
        }
    }

    // Filter: peptide amide N has at least one heavy non-carbonyl-C neighbour.
    for (const auto& [c, nj] : candidates) {
        int n_heavy_other = 0;
        for (int k : adj[nj]) {
            if (k == c) continue;
            if (atoms[k].element != Element::H) ++n_heavy_other;
        }
        if (n_heavy_other >= 1) out.emplace_back(c, nj);
    }

    return out;
}


// ============================================================================
// Connected components + piece ordering
// ============================================================================

std::vector<std::vector<int>> ConnectedComponents(const AdjacencySet& adj,
                                                    int n) {
    std::vector<std::vector<int>> comps;
    std::vector<bool> visited(n, false);
    for (int start = 0; start < n; ++start) {
        if (visited[start]) continue;
        std::vector<int> comp;
        std::stack<int> stk; stk.push(start);
        while (!stk.empty()) {
            int cur = stk.top(); stk.pop();
            if (visited[cur]) continue;
            visited[cur] = true;
            comp.push_back(cur);
            for (int nbr : adj[cur]) if (!visited[nbr]) stk.push(nbr);
        }
        std::sort(comp.begin(), comp.end());
        comps.push_back(std::move(comp));
    }
    return comps;
}


// Order the 5 pieces along the directed amide chain ACE → NCapAla → Central
// → CCapAla → NME. Returns a list of 5 pieces in chain order.
std::vector<std::vector<int>> OrderPiecesAlongAmideChain(
        const std::vector<std::vector<int>>& pieces,
        const std::vector<std::pair<int, int>>& amides) {
    int total_atoms = 0;
    for (const auto& p : pieces) total_atoms += static_cast<int>(p.size());
    std::vector<int> piece_of(total_atoms, -1);
    for (int pi = 0; pi < static_cast<int>(pieces.size()); ++pi) {
        for (int ai : pieces[pi]) piece_of[ai] = pi;
    }
    // Directed edges: piece(C) -> piece(N) per amide.
    std::map<int, int> out_edges;
    std::map<int, int> in_edges;
    for (const auto& [c, nj] : amides) {
        out_edges[piece_of[c]] = piece_of[nj];
        in_edges[piece_of[nj]] = piece_of[c];
    }
    // Chain start: piece with outgoing but no incoming amide.
    int start = -1;
    for (int pi = 0; pi < static_cast<int>(pieces.size()); ++pi) {
        if (out_edges.count(pi) && !in_edges.count(pi)) { start = pi; break; }
    }
    if (start < 0) return {};  // cyclic or disconnected — fail

    std::vector<std::vector<int>> chain;
    int cur = start;
    chain.push_back(pieces[cur]);
    while (out_edges.count(cur)) {
        cur = out_edges.at(cur);
        chain.push_back(pieces[cur]);
    }
    return chain;
}


// ============================================================================
// Canonical chemistry per piece
// ============================================================================

struct CanonicalAtom {
    std::string             name;
    Element                 element;
    AtomMechanicalIdentity  identity;   // pre-computed at canonical build time
};

struct CanonicalPiece {
    std::string                          variant_name;
    std::vector<CanonicalAtom>           atoms;

    // Build-time adjacency keyed by atom name. Used by AddBond /
    // FirstHeavyNeighbour during CanonicalPiece construction; NOT
    // consulted at runtime (perception's hot path runs over
    // adj_by_idx). Names die here — runtime identity comes from the
    // typed AtomMechanicalIdentity carried on each CanonicalAtom.
    std::map<std::string, std::set<std::string>> adj_by_name;

    // Index-keyed adjacency, populated by FinalizeAdjacency() at the
    // end of each builder (CanonicalAce / CanonicalNme /
    // CanonicalResidue / CanonicalHisVariant). The WL signature
    // computation and MatchPiece operate on this — no name lookups
    // inside perception's match path.
    std::vector<std::set<int>>           adj_by_idx;
};


// Translate adj_by_name → adj_by_idx after canonical construction.
// Called once per CanonicalPiece (cached via function-local statics),
// not in the per-record hot path.
void FinalizeAdjacency(CanonicalPiece& p) {
    std::map<std::string, int> name_to_idx;
    for (int i = 0; i < static_cast<int>(p.atoms.size()); ++i) {
        name_to_idx[p.atoms[i].name] = i;
    }
    p.adj_by_idx.assign(p.atoms.size(), {});
    for (const auto& [a_name, nbrs] : p.adj_by_name) {
        auto a_it = name_to_idx.find(a_name);
        if (a_it == name_to_idx.end()) continue;
        for (const std::string& b_name : nbrs) {
            auto b_it = name_to_idx.find(b_name);
            if (b_it == name_to_idx.end()) continue;
            p.adj_by_idx[a_it->second].insert(b_it->second);
        }
    }
}


void AddBond(CanonicalPiece& p, const std::string& a, const std::string& b) {
    auto& A = p.adj_by_name[a];
    auto& B = p.adj_by_name[b];
    // Confirm both atoms are in the piece.
    bool has_a = false, has_b = false;
    for (const auto& at : p.atoms) {
        if (at.name == a) has_a = true;
        if (at.name == b) has_b = true;
    }
    if (!has_a || !has_b) return;
    A.insert(b); B.insert(a);
}


void AddAtom(CanonicalPiece& p, const std::string& name, Element e) {
    p.atoms.push_back({name, e, AtomMechanicalIdentity{}});
    p.adj_by_name[name];  // ensure key exists
}


// Forward declarations — ACE/NME identity emitters + the two stampers.
// Cap stamping is hand-coded (ACE/NME are not in the standard-20 substrate
// table); chain stamping (NCapAla / CCapAla / Central) routes through the
// generated topology table so the LarsenResidue side speaks the same typed
// vocabulary as ComposeAtomSemantic on the protein side.
AtomMechanicalIdentity AceCapIdentity(const std::string& canonical_name);
AtomMechanicalIdentity NmeCapIdentity(const std::string& canonical_name);
void StampCapIdentities(CanonicalPiece& p, LarsenResidue::Kind kind);
void StampChainIdentitiesViaTable(CanonicalPiece& p,
                                    AminoAcid aa,
                                    std::uint8_t variant_idx);


// Build canonical chemistry for ACE cap.
CanonicalPiece CanonicalAce() {
    CanonicalPiece p;
    AddAtom(p, "C",    Element::C);
    AddAtom(p, "O",    Element::O);
    AddAtom(p, "CH3",  Element::C);
    AddAtom(p, "HH31", Element::H);
    AddAtom(p, "HH32", Element::H);
    AddAtom(p, "HH33", Element::H);
    AddBond(p, "C",   "O");
    AddBond(p, "C",   "CH3");
    AddBond(p, "CH3", "HH31");
    AddBond(p, "CH3", "HH32");
    AddBond(p, "CH3", "HH33");
    StampCapIdentities(p, LarsenResidue::Kind::AceCap);
    FinalizeAdjacency(p);
    return p;
}


// Build canonical chemistry for NME cap.
CanonicalPiece CanonicalNme() {
    CanonicalPiece p;
    AddAtom(p, "N",    Element::N);
    AddAtom(p, "H",    Element::H);
    AddAtom(p, "CH3",  Element::C);
    AddAtom(p, "HH31", Element::H);
    AddAtom(p, "HH32", Element::H);
    AddAtom(p, "HH33", Element::H);
    AddBond(p, "N",   "H");
    AddBond(p, "N",   "CH3");
    AddBond(p, "CH3", "HH31");
    AddBond(p, "CH3", "HH32");
    AddBond(p, "CH3", "HH33");
    StampCapIdentities(p, LarsenResidue::Kind::NmeCap);
    FinalizeAdjacency(p);
    return p;
}


// Bonds beyond chi-walk that the AminoAcidType representation doesn't
// encode. Mirror of the Python POC's EXTRA_BONDS table.
const std::vector<std::pair<std::string, std::string>>&
ExtraBondsFor(const std::string& code) {
    static const std::map<std::string, std::vector<std::pair<std::string, std::string>>> kTable = {
        {"ARG", {{"CZ","NH1"},{"CZ","NH2"}}},
        {"ASN", {{"CG","ND2"}}},
        {"ASP", {{"CG","OD2"}}},
        {"GLN", {{"CD","NE2"}}},
        {"GLU", {{"CD","OE2"}}},
        {"HIS", {{"CG","CD2"},{"CD2","NE2"},{"NE2","CE1"},{"CE1","ND1"}}},
        {"ILE", {{"CB","CG2"}}},
        {"LEU", {{"CG","CD2"}}},
        {"PHE", {{"CG","CD2"},{"CD1","CE1"},{"CD2","CE2"},
                 {"CE1","CZ"},{"CE2","CZ"}}},
        {"PRO", {{"CD","N"}}},
        {"THR", {{"CB","CG2"}}},
        {"TRP", {{"CG","CD2"},{"CD1","NE1"},{"NE1","CE2"},{"CD2","CE2"},
                 {"CD2","CE3"},{"CE2","CZ2"},{"CE3","CZ3"},{"CZ2","CH2"},
                 {"CZ3","CH2"}}},
        {"TYR", {{"CG","CD2"},{"CD1","CE1"},{"CD2","CE2"},
                 {"CE1","CZ"},{"CE2","CZ"},{"CZ","OH"}}},
        {"VAL", {{"CB","CG2"}}},
    };
    static const std::vector<std::pair<std::string, std::string>> kEmpty;
    auto it = kTable.find(code);
    return (it == kTable.end()) ? kEmpty : it->second;
}


// Hydrogen attachment rules. Returns the heavy-atom parent for a given
// H atom name, or empty string if not handled. Mirrors the Python POC's
// canonical_bonds_of_residue H-attach section.
std::string ParentHeavyForH(const std::string& h_name,
                              const std::set<std::string>& atoms_in_piece) {
    auto has = [&](const std::string& n) { return atoms_in_piece.count(n); };
    // Polar Hs handled at end of function (HE→NE, HZ→CZ, HH→OH, etc.).
    if (h_name == "H" || h_name == "HN")    return "N";
    if (h_name == "HA")                      return "CA";
    if (h_name == "HA2" || h_name == "HA3")  return "CA";  // Gly
    if (h_name.rfind("HB", 0) == 0)         return "CB";
    if (h_name.rfind("HG1", 0) == 0 && has("CG1")) return "CG1";
    if (h_name.rfind("HG1", 0) == 0 && has("OG1")) return "OG1";  // Thr HG1
    if (h_name.rfind("HG2", 0) == 0 && has("CG2")) return "CG2";
    if (h_name.rfind("HG2", 0) == 0 && has("CG"))  return "CG";   // most cases
    if (h_name.rfind("HG3", 0) == 0 && has("CG3")) return "CG3";
    if (h_name.rfind("HG3", 0) == 0 && has("CG"))  return "CG";
    if (h_name == "HG" && has("OG")) return "OG";  // Ser
    if (h_name == "HG" && has("SG")) return "SG";  // Cys
    if (h_name == "HG" && has("CG")) return "CG";  // Leu
    if (h_name.rfind("HD1", 0) == 0 && has("CD1")) return "CD1";
    if (h_name.rfind("HD2", 0) == 0 && has("CD2")) return "CD2";
    if (h_name.rfind("HD2", 0) == 0 && has("ND2")) return "ND2";  // Asn
    if (h_name.rfind("HD", 0) == 0 && has("CD"))   return "CD";
    if (h_name.rfind("HE1", 0) == 0 && has("CE1")) return "CE1";
    if (h_name.rfind("HE1", 0) == 0 && has("NE1")) return "NE1";  // Trp
    if (h_name.rfind("HE2", 0) == 0 && has("CE2")) return "CE2";
    if (h_name.rfind("HE2", 0) == 0 && has("NE2")) return "NE2";  // Gln/HIE
    if (h_name.rfind("HE3", 0) == 0 && has("CE3")) return "CE3";
    if (h_name.rfind("HE", 0) == 0 && has("CE"))   return "CE";
    if (h_name == "HE" && has("NE"))               return "NE";    // Arg HE
    if (h_name.rfind("HZ1", 0) == 0 && has("NZ"))  return "NZ";
    if (h_name.rfind("HZ2", 0) == 0 && has("NZ"))  return "NZ";
    if (h_name.rfind("HZ3", 0) == 0 && has("NZ"))  return "NZ";
    if (h_name.rfind("HZ2", 0) == 0 && has("CZ2")) return "CZ2";
    if (h_name.rfind("HZ3", 0) == 0 && has("CZ3")) return "CZ3";
    if (h_name == "HZ" && has("CZ"))               return "CZ";   // Phe
    if (h_name.rfind("HH1", 0) == 0)               return "NH1";
    if (h_name.rfind("HH2", 0) == 0 && has("NH2")) return "NH2";
    if (h_name.rfind("HH2", 0) == 0 && has("CH2")) return "CH2";  // Trp HH2
    if (h_name == "HH" && has("OH"))               return "OH";   // Tyr
    if (h_name == "HD1" && has("ND1"))             return "ND1";  // HID/HIP
    return "";
}


// HIS canonical variants (HID/HIE/HIP differ in which ring N carries H).
CanonicalPiece CanonicalHisVariant(const std::string& variant_name) {
    CanonicalPiece p;
    p.variant_name = variant_name;
    AddAtom(p, "N",   Element::N);
    AddAtom(p, "CA",  Element::C);
    AddAtom(p, "C",   Element::C);
    AddAtom(p, "O",   Element::O);
    AddAtom(p, "H",   Element::H);
    AddAtom(p, "HA",  Element::H);
    AddAtom(p, "CB",  Element::C);
    AddAtom(p, "HB2", Element::H);
    AddAtom(p, "HB3", Element::H);
    AddAtom(p, "CG",  Element::C);
    AddAtom(p, "ND1", Element::N);
    if (variant_name == "HID" || variant_name == "HIP") AddAtom(p, "HD1", Element::H);
    AddAtom(p, "CD2", Element::C);
    AddAtom(p, "HD2", Element::H);
    AddAtom(p, "CE1", Element::C);
    AddAtom(p, "HE1", Element::H);
    AddAtom(p, "NE2", Element::N);
    if (variant_name == "HIE" || variant_name == "HIP") AddAtom(p, "HE2", Element::H);

    AddBond(p, "N","CA"); AddBond(p, "CA","C"); AddBond(p, "C","O");
    AddBond(p, "N","H"); AddBond(p, "CA","HA"); AddBond(p, "CA","CB");
    AddBond(p, "CB","HB2"); AddBond(p, "CB","HB3"); AddBond(p, "CB","CG");
    AddBond(p, "CG","ND1"); AddBond(p, "CG","CD2");
    AddBond(p, "CD2","NE2"); AddBond(p, "NE2","CE1"); AddBond(p, "CE1","ND1");
    AddBond(p, "CD2","HD2"); AddBond(p, "CE1","HE1");
    if (variant_name == "HID" || variant_name == "HIP")
        AddBond(p, "HD1", "ND1");
    if (variant_name == "HIE" || variant_name == "HIP")
        AddBond(p, "HE2", "NE2");

    // adj_by_idx must exist before chain stamping so the methyl-collapse
    // step can count parent-H neighbours via the canonical bond graph.
    FinalizeAdjacency(p);

    // Resolve variant_idx from the variant name. Convention is shared
    // with Residue::protonation_variant_index and PerceiveLarsenTripeptide
    // (0=HID, 1=HIE, 2=HIP) — checked against AminoAcidType variant
    // ordering by ValidateVariantIndices() at test startup.
    std::uint8_t variant_idx = topology_generated::kBaseVariantIdx;
    if      (variant_name == "HID") variant_idx = 0;
    else if (variant_name == "HIE") variant_idx = 1;
    else if (variant_name == "HIP") variant_idx = 2;
    StampChainIdentitiesViaTable(p, AminoAcid::HIS, variant_idx);
    return p;
}


// Build canonical chemistry for a standard amino acid (uses
// AminoAcidType's atoms + chi_angles + ExtraBondsFor + H attachment rules).
CanonicalPiece CanonicalResidue(AminoAcid aa) {
    const AminoAcidType& aat = GetAminoAcidType(aa);
    CanonicalPiece p;

    // Atoms.
    std::set<std::string> atoms_set;
    for (const auto& a : aat.atoms) {
        AddAtom(p, a.name, a.element);
        atoms_set.insert(a.name);
    }

    // Backbone bonds.
    AddBond(p, "N","CA"); AddBond(p, "CA","C"); AddBond(p, "C","O");
    if (atoms_set.count("H"))  AddBond(p, "N","H");
    if (atoms_set.count("HA")) AddBond(p, "CA","HA");
    if (atoms_set.count("HA2")) AddBond(p, "CA","HA2");
    if (atoms_set.count("HA3")) AddBond(p, "CA","HA3");
    if (atoms_set.count("CB")) AddBond(p, "CA","CB");

    // Chi-chain bonds.
    for (const auto& chi : aat.chi_angles) {
        for (int i = 0; i < 3; ++i) {
            AddBond(p, chi.atoms[i], chi.atoms[i + 1]);
        }
    }

    // Extras.
    const std::string three_letter(aat.three_letter_code);
    for (const auto& [a, b] : ExtraBondsFor(three_letter)) {
        AddBond(p, a, b);
    }

    // Hydrogen attachments.
    for (const auto& a : aat.atoms) {
        if (a.element != Element::H) continue;
        const std::string parent = ParentHeavyForH(a.name, atoms_set);
        if (!parent.empty()) AddBond(p, a.name, parent);
    }

    // adj_by_idx must exist before chain stamping so the methyl-collapse
    // step can count parent-H neighbours via the canonical bond graph.
    FinalizeAdjacency(p);

    // Non-titratable standard residues use the base (non-variant) row;
    // titratable residues that perception drives via a non-HIS path
    // (none currently, but the protein side stays consistent) also land
    // here at the base index. HIS variants are stamped from CanonicalHisVariant.
    StampChainIdentitiesViaTable(p, aa, topology_generated::kBaseVariantIdx);
    return p;
}


// ============================================================================
// Graph isomorphism matching (K=3 Weisfeiler-Lehman)
// ============================================================================
//
// Why multi-round WL: a single-round (element, degree, neighbour-multiset)
// signature collapses chemically-distinct atoms whose 1-hop neighbourhood
// happens to look identical. The canonical case is PHE/TYR aromatic
// rings: CD1, CD2, CE1, CE2, CZ all carry signature
// (C, 3, [(C,3), (C,3), (H,1)]) — 5 chemically-different atoms in a
// single class. zip-assigning canonical names within that class scrambles
// tensors across locants (perceived CE2's tensor lands on protein CD1).
//
// K=3 WL refinement: each round's signature is the prior round's
// signature combined with the sorted multiset of neighbours' prior-round
// signatures. After K=3, CG/CD/CE/CZ are all distinguished because
// their 3-hop neighbourhoods differ (CG's neighbour CB has degree 4;
// CZ's neighbours are aromatic Cs with degree-3 only; etc.). The
// remaining residual classes are the truly-graph-automorphic pairs:
// {CD1, CD2}, {CE1, CE2}, {HD1, HD2}, {HE1, HE2}, and analogous Trp /
// Arg / Asn pairs. Those need spatial tiebreak at match time
// (AssembleCentralTyped's nearest-spatial within the relaxed candidate
// set), not graph refinement — no K resolves them.
//
// Cost: K=3 over ~50 atoms × ~20 residues = ~3000 string-tuple
// comparisons per perception call. Negligible.

using Signature = std::tuple<Element, int, std::vector<std::pair<Element, int>>>;

// One round's signature = (element, degree, sorted-multiset of
// (neighbour-element, neighbour-degree)). The "prior signatures"
// argument carries round k-1's per-atom signature; for round 0 it's
// empty and we fall back to (element, degree) only.
using SigKey = std::vector<std::uint64_t>;  // K signatures per atom, hashed

// Stable hash of a Signature for use as a numeric per-round id.
std::uint64_t HashSignature(const Signature& s) {
    // FNV-1a 64-bit over the tuple serialised as bytes.
    auto fnv = [](std::uint64_t h, std::uint64_t v) {
        h ^= v;
        h *= 0x100000001b3ull;
        return h;
    };
    std::uint64_t h = 0xcbf29ce484222325ull;
    h = fnv(h, static_cast<std::uint64_t>(std::get<0>(s)));
    h = fnv(h, static_cast<std::uint64_t>(std::get<1>(s)));
    for (const auto& [e, d] : std::get<2>(s)) {
        h = fnv(h, static_cast<std::uint64_t>(e));
        h = fnv(h, static_cast<std::uint64_t>(d));
    }
    return h;
}


// K=3 WL per-canonical-atom signature ids, keyed by canonical node
// index. Each intermediate round's hash is (element, degree, sorted-
// multiset of neighbours' prior-round hashes); round 0 uses
// (element, 0) as the seed. Walks adj_by_idx — no name lookups.
std::vector<std::uint64_t> CanonicalWLSignatures(
        const CanonicalPiece& p, int rounds = 3) {
    const int n = static_cast<int>(p.atoms.size());
    std::vector<std::uint64_t> cur(n);
    for (int i = 0; i < n; ++i) {
        Signature seed{p.atoms[i].element, 0, {}};
        cur[i] = HashSignature(seed);
    }
    for (int r = 0; r < rounds; ++r) {
        std::vector<std::uint64_t> next(n);
        for (int i = 0; i < n; ++i) {
            std::vector<std::pair<Element, int>> nbr_sig;
            for (int j : p.adj_by_idx[i]) {
                // Neighbour's prior-round signature truncated to 31 bits
                // and paired with its element gives the WL refinement
                // discriminator. Order-insensitive (sorted multiset).
                nbr_sig.emplace_back(
                    p.atoms[j].element,
                    static_cast<int>(cur[j] & 0x7fffffff));
            }
            std::sort(nbr_sig.begin(), nbr_sig.end());
            Signature s{p.atoms[i].element,
                         static_cast<int>(p.adj_by_idx[i].size()), nbr_sig};
            next[i] = HashSignature(s);
        }
        cur = std::move(next);
    }
    return cur;
}


// Same as CanonicalWLSignatures but for the perceived atoms.
std::map<int, std::uint64_t> PerceivedWLSignatures(
        const std::vector<int>& piece_atoms,
        const AdjacencySet& sub_adj,
        const std::vector<TripeptideDftAtom>& atoms,
        int rounds = 3) {
    std::map<int, std::uint64_t> cur;
    for (int i : piece_atoms) {
        Signature seed{atoms[i].element, 0, {}};
        cur[i] = HashSignature(seed);
    }
    for (int r = 0; r < rounds; ++r) {
        std::map<int, std::uint64_t> next;
        for (int i : piece_atoms) {
            std::vector<std::pair<Element, int>> nbr_sig;
            for (int j : sub_adj[i]) {
                const std::uint64_t prior = cur.at(j);
                nbr_sig.emplace_back(atoms[j].element,
                                      static_cast<int>(prior & 0x7fffffff));
            }
            std::sort(nbr_sig.begin(), nbr_sig.end());
            Signature s{atoms[i].element,
                         static_cast<int>(sub_adj[i].size()), nbr_sig};
            next[i] = HashSignature(s);
        }
        cur = std::move(next);
    }
    return cur;
}


// Match a perceived piece's atoms to canonical-piece NODE INDICES by
// K=3 WL signature class. Returns a PieceMatch with:
//
//   - canon_idx_by_dft: perceived atom dft_idx → canonical node index
//     (position in `canon.atoms`). The canonical node carries the
//     typed AtomMechanicalIdentity; consumers copy that identity at
//     EmitPiece time. Names are not part of the matching contract;
//     they live only inside CanonicalPiece construction.
//   - canonical_ambiguous: perceived atom dft_idx → bool, true iff
//     the assignment was made within a WL signature class of size ≥ 2
//     (a graph-automorphic pair such as PHE CD1/CD2 that K rounds of
//     WL cannot split). Singleton classes (the common case after K=3,
//     including chemistry-distinct branches like ILE CG1/CG2) flag
//     false — BranchAddress + DiastereotopicIndex bind strictly
//     downstream. Multi-atom classes flag true — the assembler drops
//     those two fields and resolves by nearest-spatial.
//
// Empty canon_idx_by_dft on failure (signature-set / cardinality
// mismatch).
struct PieceMatch {
    std::map<int, int>  canon_idx_by_dft;
    std::map<int, bool> canonical_ambiguous;
};

PieceMatch MatchPiece(
        const std::vector<int>& piece_atoms,
        const std::vector<TripeptideDftAtom>& full_atoms,
        const AdjacencySet& full_adj,
        const CanonicalPiece& canon) {

    // Restrict adjacency to atoms inside the piece (the cut adjacency).
    std::set<int> piece_set(piece_atoms.begin(), piece_atoms.end());
    AdjacencySet sub_adj(full_adj.size());
    for (int i : piece_atoms) {
        for (int j : full_adj[i]) {
            if (piece_set.count(j)) sub_adj[i].insert(j);
        }
    }

    // K=3 WL signatures, both sides.
    const std::vector<std::uint64_t> canon_sigs =
        CanonicalWLSignatures(canon);
    const std::map<int, std::uint64_t> perc_sigs =
        PerceivedWLSignatures(piece_atoms, sub_adj, full_atoms);

    // Group canonical nodes by WL signature.
    std::map<std::uint64_t, std::vector<int>> canon_by_sig;
    for (int ci = 0; ci < static_cast<int>(canon.atoms.size()); ++ci) {
        canon_by_sig[canon_sigs[ci]].push_back(ci);
    }

    // Group perceived atoms.
    std::map<std::uint64_t, std::vector<int>> perceived_by_sig;
    for (int i : piece_atoms) {
        perceived_by_sig[perc_sigs.at(i)].push_back(i);
    }

    // Sig sets must match; per-class cardinality must match.
    if (canon_by_sig.size() != perceived_by_sig.size()) return {};
    for (const auto& [sig, canon_ids] : canon_by_sig) {
        auto it = perceived_by_sig.find(sig);
        if (it == perceived_by_sig.end()) return {};
        if (it->second.size() != canon_ids.size()) return {};
    }

    // Map perceived → canonical node within each signature class. For
    // singleton classes (the common case after K=3 WL), this is a
    // unique assignment and canonical_ambiguous=false — the canonical
    // node's typed identity (including BranchAddress and
    // DiastereotopicIndex) binds strictly downstream. For multi-atom
    // classes (residual symmetric pairs like CD1/CD2), the within-
    // class assignment is by container order; canonical_ambiguous=true
    // tells AssembleCentralTyped to drop BranchAddress +
    // DiastereotopicIndex from the equality check and resolve by
    // nearest-spatial within the equivalence class.
    PieceMatch out;
    for (const auto& [sig, canon_ids] : canon_by_sig) {
        const auto& perc_ids = perceived_by_sig.at(sig);
        const bool ambiguous = canon_ids.size() > 1;
        for (std::size_t k = 0; k < canon_ids.size(); ++k) {
            out.canon_idx_by_dft[perc_ids[k]]     = canon_ids[k];
            out.canonical_ambiguous[perc_ids[k]]  = ambiguous;
        }
    }
    return out;
}


// ============================================================================
// ACE / NME identity emission (hand-coded since they're not standard residues).
// ============================================================================

AtomMechanicalIdentity AceCapIdentity(const std::string& canonical_name) {
    AtomMechanicalIdentity id;
    if (canonical_name == "C") {
        id.element = Element::C;
        id.backbone_role = BackboneRole::CarbonylCarbon;
    } else if (canonical_name == "O") {
        id.element = Element::O;
        id.backbone_role = BackboneRole::CarbonylOxygen;
    } else if (canonical_name == "CH3") {
        id.element = Element::C;
    } else if (canonical_name.rfind("HH3", 0) == 0) {
        id.element = Element::H;
    }
    return id;
}


AtomMechanicalIdentity NmeCapIdentity(const std::string& canonical_name) {
    AtomMechanicalIdentity id;
    if (canonical_name == "N") {
        id.element = Element::N;
        id.backbone_role = BackboneRole::Nitrogen;
    } else if (canonical_name == "H") {
        id.element = Element::H;
        id.backbone_role = BackboneRole::AmideHydrogen;
    } else if (canonical_name == "CH3") {
        id.element = Element::C;
    } else if (canonical_name.rfind("HH3", 0) == 0) {
        id.element = Element::H;
    }
    return id;
}


// Find the first heavy neighbour of `atom_name` in the canonical piece's
// bond graph. Used as the `parent_name` argument to ParseAtomName, which
// disambiguates H atoms with branch indices (e.g., HG21 on Thr CG2).
// Returns empty string if no heavy neighbour.
std::string FirstHeavyNeighbour(const std::string& atom_name,
                                  const CanonicalPiece& p) {
    auto it = p.adj_by_name.find(atom_name);
    if (it == p.adj_by_name.end()) return "";
    for (const std::string& n : it->second) {
        for (const auto& a : p.atoms) {
            if (a.name == n && a.element != Element::H) return n;
        }
    }
    return "";
}


// Stamp typed identity onto ACE / NME cap atoms. Hand-coded because
// the cap atom names ("CH3", "HH31") aren't in the standard-20
// substrate table; the generated topology tables only cover the 20
// canonical residues. The cap identities below are synthesised to
// give the cap's carbonyl-C / amide-N / etc. the same typed
// BackboneRole/Locant vocabulary the protein side uses for chain
// atoms — so cap atoms can still match protein-side N/CA/C/O slots
// via typed dispatch in EmitPiece and AssembleAlaCap.
void StampCapIdentities(CanonicalPiece& p, LarsenResidue::Kind kind) {
    for (auto& a : p.atoms) {
        if (kind == LarsenResidue::Kind::AceCap) {
            a.identity = AceCapIdentity(a.name);
        } else if (kind == LarsenResidue::Kind::NmeCap) {
            a.identity = NmeCapIdentity(a.name);
        }
    }
}


// Stamp typed identity onto each chain-piece atom by looking up the
// authoritative AtomSemanticTable row in the generated topology table.
// This matches the protein side's `ComposeAtomSemantic` discipline
// (`src/LegacyAmberTopology.cpp:233`) so both substrates speak the
// same typed vocabulary:
//
//   1. Parse name + heavy-parent-name → ParsedAtomName.
//   2. Build a tentative AtomMechanicalIdentity from parsed flags.
//   3. Apply methyl-H pseudoatom collapse: if the H's heavy parent has
//      3+ H neighbours in the canonical bond graph (a methyl), clear
//      DiastereotopicIndex (mirrors the substrate generator's
//      `e.di_index = None when e.pseudoatom.kind == M` rule).
//   4. LookupBy(amino_acid, variant_idx, collapsed_identity) on the
//      generated table; the returned row's identity is the
//      authoritative typed tuple.
//
// A lookup miss is FATAL — same discipline as ComposeAtomSemantic:
// a missing row means perception's canonical-chemistry definition
// doesn't match the generated table, which is a project-side bug
// (either the chemistry author missed an atom, or the substrate
// generator produced a non-matching row). Either way we fail loud
// rather than emit a default-constructed identity into perception's
// hot path.
void StampChainIdentitiesViaTable(CanonicalPiece& p,
                                    AminoAcid aa,
                                    std::uint8_t variant_idx) {
    namespace gen = topology_generated;
    for (int i = 0; i < static_cast<int>(p.atoms.size()); ++i) {
        auto& a = p.atoms[i];

        // Parse name + parent-name (parent is empty for non-H atoms).
        const std::string parent = (a.element == Element::H)
            ? FirstHeavyNeighbour(a.name, p) : "";
        const gen::ParsedAtomName parsed = gen::ParseAtomName(a.name, parent);

        AtomMechanicalIdentity ident;
        ident.element       = a.element;
        ident.locant        = parsed.locant;
        ident.branch        = parsed.branch;
        ident.di_index      = parsed.di_index;
        ident.backbone_role = parsed.backbone_role;

        // Methyl-H collapse: parent has 3+ H neighbours in the
        // canonical bond graph → methyl → clear DI. Mirrors
        // ComposeAtomSemantic lines 356-382.
        if (a.element == Element::H &&
            ident.di_index != DiastereotopicIndex::None &&
            !parent.empty()) {
            int parent_idx = -1;
            for (int k = 0; k < static_cast<int>(p.atoms.size()); ++k) {
                if (p.atoms[k].name == parent) { parent_idx = k; break; }
            }
            if (parent_idx >= 0) {
                int h_children = 0;
                for (int nbr : p.adj_by_idx[parent_idx]) {
                    if (p.atoms[nbr].element == Element::H) ++h_children;
                }
                if (h_children >= 3) {
                    ident.di_index = DiastereotopicIndex::None;
                }
            }
        }

        // Generated topology table is the authority. The row's
        // identity equals the (post-collapse) lookup key by
        // construction; we still copy from the row so any future
        // table-side adjustment to the canonical identity propagates
        // here automatically.
        const AtomSemanticTable* row =
            gen::LookupBy(aa, variant_idx, ident);
        if (row == nullptr) {
            std::fprintf(stderr,
                "FATAL: LarsenResidue StampChainIdentitiesViaTable "
                "lookup miss: aa=%s variant_idx=%u atom_name='%s' "
                "identity=(element=%u locant=%u branch={%u,%u} "
                "di=%u role=%u). The canonical chemistry definition "
                "in LarsenResidue.cpp does not match the generated "
                "topology table at src/generated/LegacyAmberSemanticTables.cpp "
                "— update one to match the other. Positive coverage "
                "for this path lives in LarsenResiduePerceptionTest"
                ".AllCombinationsPerceiveCleanly: that test walks all "
                "20 (residue, frame_type) combinations and the AAA "
                "Gaussian log; if any LookupBy returned null the test "
                "would abort here instead of pass. A passing test "
                "suite means this FATAL is unreachable from the "
                "standard-20 substrate; reaching it indicates a "
                "development-time edit to one side without the other.\n",
                GetAminoAcidType(aa).three_letter_code,
                static_cast<unsigned>(variant_idx),
                a.name.c_str(),
                static_cast<unsigned>(ident.element),
                static_cast<unsigned>(ident.locant),
                static_cast<unsigned>(ident.branch.outer),
                static_cast<unsigned>(ident.branch.inner),
                static_cast<unsigned>(ident.di_index),
                static_cast<unsigned>(ident.backbone_role));
            std::abort();
        }
        ident.element       = row->element;
        ident.locant        = row->locant;
        ident.branch        = row->branch;
        ident.di_index      = row->di_index;
        ident.backbone_role = row->backbone_role;
        a.identity = ident;
    }
}


// ============================================================================
// Per-piece emission
// ============================================================================

bool EmitPiece(LarsenResidue::Kind kind,
                const std::vector<int>& piece_atoms,
                const std::vector<TripeptideDftAtom>& full_atoms,
                const AdjacencySet& full_adj,
                const CanonicalPiece& canon,
                LarsenResidue& out) {
    if (static_cast<int>(piece_atoms.size()) !=
        static_cast<int>(canon.atoms.size())) {
        return false;
    }

    PieceMatch match = MatchPiece(piece_atoms, full_atoms, full_adj, canon);
    if (match.canon_idx_by_dft.empty()) return false;

    out.kind = kind;
    out.atoms.clear();
    out.atoms.reserve(piece_atoms.size());

    // Reset slot cache.
    out.N_idx = out.H_idx = out.CA_idx = out.HA_idx
              = out.CB_idx = out.C_idx = out.O_idx = -1;

    // Build local-idx ordering preserving the DFT-row order within the piece.
    std::vector<int> ordered = piece_atoms;
    std::sort(ordered.begin(), ordered.end(),
               [&](int a, int b) {
                   return full_atoms[a].atom_idx < full_atoms[b].atom_idx;
               });

    for (int local_idx = 0; local_idx < static_cast<int>(ordered.size());
          ++local_idx) {
        const int dft_idx = ordered[local_idx];
        const TripeptideDftAtom& a = full_atoms[dft_idx];

        // Read typed identity directly from the assigned canonical
        // node. Names die at canonical construction — runtime
        // identity is the AtomMechanicalIdentity stamped via
        // ParseAtomName at CanonicalPiece build time.
        const int canon_idx = match.canon_idx_by_dft.at(dft_idx);
        const CanonicalAtom& canon_atom = canon.atoms[canon_idx];

        LarsenResidue::PerAtom pa;
        pa.dft_atom_idx     = a.atom_idx;
        pa.element          = a.element;
        pa.position         = a.position;
        pa.shielding_tensor = a.shielding_tensor;
        pa.isotropic        = a.isotropic;
        pa.anisotropy       = a.anisotropy;
        pa.t2_components    = a.t2_components;
        pa.canonical_assignment_ambiguous =
            match.canonical_ambiguous.at(dft_idx);
        pa.identity         = canon_atom.identity;
        out.atoms.push_back(std::move(pa));

        // Slot cache — typed dispatch on the canonical identity. No
        // name comparisons. BackboneRole carries the peptide-amide
        // slot semantics for N / H / CA / HA / C / O; CB is the
        // (Element::C, Locant::Beta) atom on every standard residue
        // except GLY (which has no CB and falls through). For ACE
        // (no N/CA/H/HA/CB) and NME (no CA/HA/CB/C/O), the typed
        // identity simply doesn't carry the unrelated roles, so the
        // unused slots stay -1.
        const auto& id = canon_atom.identity;
        if (id.backbone_role == BackboneRole::Nitrogen)       out.N_idx  = local_idx;
        if (id.backbone_role == BackboneRole::AmideHydrogen)  out.H_idx  = local_idx;
        if (id.backbone_role == BackboneRole::AlphaCarbon)    out.CA_idx = local_idx;
        if (id.backbone_role == BackboneRole::AlphaHydrogen)  out.HA_idx = local_idx;
        if (id.backbone_role == BackboneRole::CarbonylCarbon) out.C_idx  = local_idx;
        if (id.backbone_role == BackboneRole::CarbonylOxygen) out.O_idx  = local_idx;
        if (id.element == Element::C && id.locant == Locant::Beta)
            out.CB_idx = local_idx;
    }

    // Emit bond list (heavy-heavy + heavy-H edges within the piece, in
    // local-index space).
    std::map<int, int> dft_to_local;
    for (int li = 0; li < static_cast<int>(ordered.size()); ++li) {
        dft_to_local[ordered[li]] = li;
    }
    out.bonds.clear();
    std::set<int> piece_set(piece_atoms.begin(), piece_atoms.end());
    for (int i : piece_atoms) {
        for (int j : full_adj[i]) {
            if (!piece_set.count(j)) continue;
            if (i < j) {
                LarsenResidue::Bond b;
                b.a = dft_to_local[i];
                b.b = dft_to_local[j];
                b.order = 1;  // heuristic; refinement deferred
                out.bonds.push_back(b);
            }
        }
    }

    return true;
}


}  // anonymous namespace


// ============================================================================
// LarsenResidue methods
// ============================================================================

int LarsenResidue::LookupByIdentity(const AtomMechanicalIdentity& id) const {
    for (int i = 0; i < static_cast<int>(atoms.size()); ++i) {
        if (atoms[i].identity == id) return i;
    }
    return -1;
}


bool LarsenResidue::HasAllRequiredSlots() const {
    switch (kind) {
        case Kind::AceCap:
            // ACE has no BB N/CA/C, but it does have a carbonyl C+O.
            // Required slots: C, O (the carbonyl pair); 6 atoms total.
            return atoms.size() == 6;
        case Kind::NmeCap:
            // NME has BB N + H but no carbonyl C+O. 6 atoms total.
            return N_idx >= 0 && H_idx >= 0 && atoms.size() == 6;
        case Kind::NCapAla:
        case Kind::CCapAla:
            // ALA caps: 10 atoms with N/H/CA/HA/CB/C/O all populated.
            return N_idx >= 0 && H_idx >= 0 && CA_idx >= 0 && HA_idx >= 0
                && CB_idx >= 0 && C_idx >= 0 && O_idx >= 0
                && atoms.size() == 10;
        case Kind::Central:
            // Central residue: N/CA/C/O always; H present except PRO; HA
            // present except GLY (which has HA2/HA3); CB present except GLY.
            if (N_idx < 0 || CA_idx < 0 || C_idx < 0 || O_idx < 0) return false;
            if (residue != AminoAcid::PRO && H_idx < 0) return false;
            if (residue != AminoAcid::GLY && HA_idx < 0) return false;
            if (residue != AminoAcid::GLY && CB_idx < 0) return false;
            return true;
    }
    return false;
}


// ============================================================================
// LarsenTripeptide methods
// ============================================================================

int LarsenTripeptide::TotalAtoms() const {
    return static_cast<int>(ace.atoms.size() + n_cap.atoms.size()
                           + central.atoms.size() + c_cap.atoms.size()
                           + nme.atoms.size());
}


std::pair<const LarsenResidue*, int>
LarsenTripeptide::FindByDftIdx(int dft_atom_idx) const {
    for (const LarsenResidue* p : {&ace, &n_cap, &central, &c_cap, &nme}) {
        for (int i = 0; i < static_cast<int>(p->atoms.size()); ++i) {
            if (p->atoms[i].dft_atom_idx == dft_atom_idx) return {p, i};
        }
    }
    return {nullptr, -1};
}


// ============================================================================
// PerceiveLarsenTripeptide
// ============================================================================

// Emit a structured perception-failure warning. Carries calc_id +
// tripeptide + frame_type so failures can be grouped in the audit log.
static void LogPerceptionFailure(const TripeptideDftRecord& rec,
                                   const std::string& reason) {
    OperationLog::Warn(
        "PerceiveLarsenTripeptide",
        "calc_id=" + std::to_string(rec.calc_id) +
        " tripeptide=" + rec.tripeptide +
        " frame_type=" + rec.frame_type +
        " reason: " + reason);
}


std::optional<LarsenTripeptide> PerceiveLarsenTripeptide(
        const TripeptideDftRecord& rec,
        AminoAcid expected_central,
        int his_variant_hint) {
    if (rec.atoms.empty()) {
        LogPerceptionFailure(rec, "rec.atoms is empty");
        return std::nullopt;
    }

    AdjacencySet adj = BuildBondGraph(rec.atoms);
    auto amides      = FindPeptideAmides(rec.atoms, adj);
    if (amides.size() != 4) {
        LogPerceptionFailure(rec,
            "expected 4 peptide amides, found " +
            std::to_string(amides.size()));
        return std::nullopt;
    }

    // Cut amides → 5 components.
    AdjacencySet cut = adj;
    for (const auto& [c, nj] : amides) {
        cut[c].erase(nj); cut[nj].erase(c);
    }
    auto comps = ConnectedComponents(cut, static_cast<int>(rec.atoms.size()));
    if (comps.size() != 5) {
        std::string sizes;
        for (const auto& c : comps) {
            if (!sizes.empty()) sizes += ",";
            sizes += std::to_string(c.size());
        }
        LogPerceptionFailure(rec,
            "expected 5 connected components after amide cut, found " +
            std::to_string(comps.size()) + " (sizes=[" + sizes + "])");
        return std::nullopt;
    }

    auto ordered = OrderPiecesAlongAmideChain(comps, amides);
    if (ordered.size() != 5) {
        LogPerceptionFailure(rec,
            "OrderPiecesAlongAmideChain returned " +
            std::to_string(ordered.size()) +
            " pieces — chain ordering failed (cyclic or disconnected)");
        return std::nullopt;
    }

    LarsenTripeptide out;
    const std::array<LarsenResidue::Kind, 5> kinds = {
        LarsenResidue::Kind::AceCap,
        LarsenResidue::Kind::NCapAla,
        LarsenResidue::Kind::Central,
        LarsenResidue::Kind::CCapAla,
        LarsenResidue::Kind::NmeCap,
    };
    const std::array<LarsenResidue*, 5> targets = {
        &out.ace, &out.n_cap, &out.central, &out.c_cap, &out.nme
    };

    // Cached canonical chemistry. String parsing happens once at first
    // call per residue-kind, never per record. Thread-safe per C++11
    // function-local static init rules.
    auto cached_ace = []() -> const CanonicalPiece& {
        static const CanonicalPiece p = CanonicalAce(); return p;
    };
    auto cached_nme = []() -> const CanonicalPiece& {
        static const CanonicalPiece p = CanonicalNme(); return p;
    };
    auto cached_residue = [](AminoAcid aa) -> const CanonicalPiece& {
        static const auto table = []() {
            std::array<CanonicalPiece, 21> t;
            for (int i = 0; i < 21; ++i) {
                t[i] = CanonicalResidue(static_cast<AminoAcid>(i));
            }
            return t;
        }();
        int idx = static_cast<int>(aa);
        if (idx < 0 || idx >= 21) idx = 20;
        return table[idx];
    };
    auto cached_his = [](const char* variant) -> const CanonicalPiece& {
        static const CanonicalPiece hid = CanonicalHisVariant("HID");
        static const CanonicalPiece hie = CanonicalHisVariant("HIE");
        static const CanonicalPiece hip = CanonicalHisVariant("HIP");
        if (std::string(variant) == "HID") return hid;
        if (std::string(variant) == "HIE") return hie;
        return hip;
    };

    static const std::array<const char*, 5> kKindNames = {
        "AceCap", "NCapAla", "Central", "CCapAla", "NmeCap"
    };

    for (int k = 0; k < 5; ++k) {
        const auto kind = kinds[k];
        const auto& piece = ordered[k];
        LarsenResidue& tgt = *targets[k];
        const std::string kind_label =
            std::string("piece[") + std::to_string(k) + "/" +
            kKindNames[k] + "]";

        if (kind == LarsenResidue::Kind::AceCap) {
            if (!EmitPiece(kind, piece, rec.atoms, adj, cached_ace(), tgt)) {
                LogPerceptionFailure(rec, kind_label + " EmitPiece failed "
                    "(piece_n=" + std::to_string(piece.size()) + ")");
                return std::nullopt;
            }
            tgt.residue = AminoAcid::Unknown;
        } else if (kind == LarsenResidue::Kind::NmeCap) {
            if (!EmitPiece(kind, piece, rec.atoms, adj, cached_nme(), tgt)) {
                LogPerceptionFailure(rec, kind_label + " EmitPiece failed "
                    "(piece_n=" + std::to_string(piece.size()) + ")");
                return std::nullopt;
            }
            tgt.residue = AminoAcid::Unknown;
        } else if (kind == LarsenResidue::Kind::NCapAla ||
                   kind == LarsenResidue::Kind::CCapAla) {
            if (!EmitPiece(kind, piece, rec.atoms, adj,
                            cached_residue(AminoAcid::ALA), tgt)) {
                LogPerceptionFailure(rec, kind_label + " EmitPiece failed "
                    "(piece_n=" + std::to_string(piece.size()) + ")");
                return std::nullopt;
            }
            tgt.residue = AminoAcid::ALA;
        } else {  // Central
            tgt.residue = expected_central;
            if (expected_central == AminoAcid::HIS) {
                // Strict variant enforcement:
                //   - his_variant_hint ≥ 0: ONLY the hinted variant is
                //     accepted. Mismatch → fail loud, perception
                //     returns nullopt. This is the contract the header
                //     documents and the only way to guarantee a HID
                //     protein doesn't silently consume HIE/HIP DB rows.
                //   - his_variant_hint == -1 (no hint): try HID/HIE/HIP
                //     in canonical order, accept the first match.
                //     Calculator sites that don't have a protein-side
                //     variant available rely on this; producing SOME
                //     LarsenTripeptide is better than declining.
                static const std::array<const char*, 3> kCanonOrder = {
                    "HID", "HIE", "HIP"
                };
                const char* hinted = nullptr;
                if      (his_variant_hint == 0) hinted = "HID";
                else if (his_variant_hint == 1) hinted = "HIE";
                else if (his_variant_hint == 2) hinted = "HIP";

                std::vector<const char*> try_order;
                if (hinted) {
                    try_order.push_back(hinted);
                } else {
                    for (const char* v : kCanonOrder) try_order.push_back(v);
                }

                bool any_match = false;
                std::string tried;
                int matched_variant_idx = -1;
                for (const char* v : try_order) {
                    const auto& canon = cached_his(v);
                    if (!tried.empty()) tried += ",";
                    tried += std::string(v) + "(n=" +
                              std::to_string(canon.atoms.size()) + ")";
                    if (static_cast<int>(piece.size()) !=
                        static_cast<int>(canon.atoms.size())) continue;
                    if (EmitPiece(kind, piece, rec.atoms, adj, canon, tgt)) {
                        any_match = true;
                        // Canonical variant index: 0=HID, 1=HIE, 2=HIP.
                        if      (std::string(v) == "HID") matched_variant_idx = 0;
                        else if (std::string(v) == "HIE") matched_variant_idx = 1;
                        else                              matched_variant_idx = 2;
                        break;
                    }
                }
                if (!any_match) {
                    if (hinted) {
                        LogPerceptionFailure(rec, kind_label +
                            " HIS hinted variant " + std::string(hinted) +
                            " does not match perceived piece_n=" +
                            std::to_string(piece.size()) +
                            "; declining (hint=" +
                            std::to_string(his_variant_hint) +
                            " tried " + tried +
                            ") — DB row's protonation form differs "
                            "from protein residue.protonation_variant_index");
                    } else {
                        LogPerceptionFailure(rec, kind_label +
                            " no HIS variant matched perceived piece_n=" +
                            std::to_string(piece.size()) +
                            " (no hint provided; tried " + tried + ")");
                    }
                    return std::nullopt;
                }
                out.central_variant_index = matched_variant_idx;
            } else {
                const auto& canon = cached_residue(expected_central);
                if (!EmitPiece(kind, piece, rec.atoms, adj, canon, tgt)) {
                    LogPerceptionFailure(rec, kind_label +
                        " EmitPiece failed (expected_residue=" +
                        std::string(GetAminoAcidType(expected_central).three_letter_code) +
                        " canonical_n=" + std::to_string(canon.atoms.size()) +
                        " piece_n=" + std::to_string(piece.size()) + ")");
                    return std::nullopt;
                }
            }
        }

        if (!tgt.HasAllRequiredSlots()) {
            LogPerceptionFailure(rec, kind_label +
                " missing required slots after emit (residue=" +
                std::string(GetAminoAcidType(tgt.residue).three_letter_code) +
                ")");
            return std::nullopt;
        }
    }

    return out;
}


}  // namespace nmr

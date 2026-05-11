// LarsenResidue.cpp -- bond-graph perception of Larsen ProCS15
// tripeptide DFT records.
//
// Implementation mirrors the Python POC at
//   scripts/perceive_larsen_tripeptide.py
// which was validated against all 20 DB (tripeptide, frame_type)
// combinations AND against the original Larsen Gaussian log
// (AAA_4_54_nmr.log byte-for-byte match with DB row 21755).
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
    std::map<std::string, std::set<std::string>> adj_by_name;
};


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


// Forward declarations — ACE/NME identity emitters + StampCanonicalIdentities
// live below; the Canonical{Ace,Nme,Residue,HisVariant} functions call them.
AtomMechanicalIdentity AceCapIdentity(const std::string& canonical_name);
AtomMechanicalIdentity NmeCapIdentity(const std::string& canonical_name);
void StampCanonicalIdentities(CanonicalPiece& p, LarsenResidue::Kind kind);


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
    StampCanonicalIdentities(p, LarsenResidue::Kind::AceCap);
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
    StampCanonicalIdentities(p, LarsenResidue::Kind::NmeCap);
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
    StampCanonicalIdentities(p, LarsenResidue::Kind::Central);
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

    StampCanonicalIdentities(p, LarsenResidue::Kind::Central);
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


// K=3 WL per-canonical-atom signature ids. Returns a map name → final
// (K=3) signature hash. Each intermediate round's hash is a function of
// (element, degree, sorted-multiset of neighbours' prior-round hashes).
// Round 0 uses (element, 0) as the seed.
std::map<std::string, std::uint64_t> CanonicalWLSignatures(
        const CanonicalPiece& p, int rounds = 3) {
    // Round 0: bare (element, 0) hashed.
    std::map<std::string, std::uint64_t> cur;
    for (const auto& a : p.atoms) {
        Signature seed{a.element, 0, {}};
        cur[a.name] = HashSignature(seed);
    }
    // Subsequent rounds.
    for (int r = 0; r < rounds; ++r) {
        std::map<std::string, std::uint64_t> next;
        for (const auto& a : p.atoms) {
            std::vector<std::pair<Element, int>> nbr_sig;
            const auto& nbrs = p.adj_by_name.at(a.name);
            for (const std::string& n : nbrs) {
                // Treat each neighbour's prior-round signature as a
                // pseudo-(Element, degree) pair so we can re-use the
                // Signature tuple shape. The "element" slot carries the
                // neighbour's actual element; the "degree" slot carries
                // the truncated prior-round hash. Together they
                // discriminate by both immediate chemistry and recursive
                // neighbourhood structure.
                Element ne = Element::Unknown;
                for (const auto& b : p.atoms) {
                    if (b.name == n) { ne = b.element; break; }
                }
                const std::uint64_t prior = cur.at(n);
                nbr_sig.emplace_back(ne, static_cast<int>(prior & 0x7fffffff));
            }
            std::sort(nbr_sig.begin(), nbr_sig.end());
            Signature s{a.element, static_cast<int>(nbrs.size()), nbr_sig};
            next[a.name] = HashSignature(s);
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


// Match a perceived piece's atoms to canonical atom names by K=3 WL
// signature class. Returns {perceived_idx -> canonical name} or empty
// map on failure. Within remaining equivalence classes (graph-automorphic
// atom pairs like aromatic CD1/CD2), the assignment is by container
// order — the downstream nearest-spatial in AssembleCentralTyped's
// relaxed candidate path resolves the within-pair swap.
std::map<int, std::string> MatchPiece(
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
    const auto canon_sigs = CanonicalWLSignatures(canon);
    const auto perc_sigs = PerceivedWLSignatures(piece_atoms, sub_adj,
                                                   full_atoms);

    // Group canonical atoms by WL signature.
    std::map<std::uint64_t, std::vector<std::string>> canon_by_sig;
    for (const auto& a : canon.atoms) {
        canon_by_sig[canon_sigs.at(a.name)].push_back(a.name);
    }

    // Group perceived atoms.
    std::map<std::uint64_t, std::vector<int>> perceived_by_sig;
    for (int i : piece_atoms) {
        perceived_by_sig[perc_sigs.at(i)].push_back(i);
    }

    // Sig sets must match; per-class cardinality must match.
    if (canon_by_sig.size() != perceived_by_sig.size()) return {};
    for (const auto& [sig, names] : canon_by_sig) {
        auto it = perceived_by_sig.find(sig);
        if (it == perceived_by_sig.end()) return {};
        if (it->second.size() != names.size()) return {};
    }

    // Map perceived → canonical within each signature class. For
    // singleton classes (the common case after K=3 WL), this is a
    // unique assignment. For the residual symmetric pairs (CD1/CD2 etc.),
    // the assignment is by container order; the downstream relaxed-
    // match in AssembleCentralTyped resolves any swap by nearest-spatial.
    std::map<int, std::string> out;
    for (const auto& [sig, names] : canon_by_sig) {
        const auto& idxs = perceived_by_sig.at(sig);
        for (size_t k = 0; k < names.size(); ++k) {
            out[idxs[k]] = names[k];
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


// Translate a canonical (residue, atom_name) pair into AtomMechanicalIdentity
// via the existing ParseAtomName authority.
//
// Called ONCE per canonical atom at CanonicalPiece build time (not per
// perceived atom per record). The pre-computed identity is stored on
// `CanonicalAtom::identity` so the per-record perception loop only does
// typed-enum copies — no string work inside the hot path.
AtomMechanicalIdentity CanonicalIdentity(const std::string& atom_name,
                                          const std::string& parent_name) {
    const auto parsed = topology_generated::ParseAtomName(atom_name, parent_name);
    AtomMechanicalIdentity id;
    id.element       = parsed.element;
    id.locant        = parsed.locant;
    id.branch        = parsed.branch;
    id.di_index      = parsed.di_index;
    id.backbone_role = parsed.backbone_role;
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


// Stamp typed identity onto every CanonicalAtom in the piece. Called
// after the bond graph is fully built (so FirstHeavyNeighbour works).
// For ACE/NME caps, identity is hand-coded (the atom names like "CH3",
// "HH31" don't parse cleanly through ParseAtomName which is built for
// standard-20 residue atoms).
void StampCanonicalIdentities(CanonicalPiece& p,
                                LarsenResidue::Kind kind) {
    for (auto& a : p.atoms) {
        if (kind == LarsenResidue::Kind::AceCap) {
            a.identity = AceCapIdentity(a.name);
        } else if (kind == LarsenResidue::Kind::NmeCap) {
            a.identity = NmeCapIdentity(a.name);
        } else {
            const std::string parent = (a.element == Element::H)
                                        ? FirstHeavyNeighbour(a.name, p) : "";
            a.identity = CanonicalIdentity(a.name, parent);
        }
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

    auto mapping = MatchPiece(piece_atoms, full_atoms, full_adj, canon);
    if (mapping.empty()) return false;

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
        const std::string& cname = mapping.at(dft_idx);

        LarsenResidue::PerAtom pa;
        pa.dft_atom_idx     = a.atom_idx;
        pa.element          = a.element;
        pa.position         = a.position;
        pa.shielding_tensor = a.shielding_tensor;
        pa.isotropic        = a.isotropic;
        pa.anisotropy       = a.anisotropy;
        pa.t2_components    = a.t2_components;

        // Typed identity — copied from the pre-stamped canonical atom.
        // No ParseAtomName call inside this loop; strings stay at the
        // CanonicalPiece-construction boundary (once per residue, not
        // per record per atom).
        for (const auto& ca : canon.atoms) {
            if (ca.name == cname) { pa.identity = ca.identity; break; }
        }
        out.atoms.push_back(std::move(pa));

        // Slot cache (BB roles).
        if (cname == "N")  out.N_idx  = local_idx;
        if (cname == "H")  out.H_idx  = local_idx;
        if (cname == "CA") out.CA_idx = local_idx;
        if (cname == "HA") out.HA_idx = local_idx;
        if (cname == "CB") out.CB_idx = local_idx;
        if (cname == "C")  out.C_idx  = local_idx;
        if (cname == "O")  out.O_idx  = local_idx;
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

#include "RingTopology.h"

#include "Atom.h"
#include "Bond.h"
#include "CovalentTopology.h"
#include "Residue.h"
#include "Ring.h"

#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
#include <utility>

namespace nmr {

namespace {

// Read the RingPositionLabel for a given ring kind from any of the
// three RingPosition slots. TRP bridge atoms carry both 5-ring and
// 6-ring labels, so walkers must search all slots.
RingPositionLabel LabelForKind(const AtomSemanticTable& sem,
                               RingSystemKind kind) {
    if (sem.ring_position.primary.ring == kind) {
        return sem.ring_position.primary.position;
    }
    if (sem.ring_position.secondary.ring == kind) {
        return sem.ring_position.secondary.position;
    }
    if (sem.ring_position.tertiary.ring == kind) {
        return sem.ring_position.tertiary.position;
    }
    return RingPositionLabel::NotInRing;
}

bool AreBonded(size_t a, size_t b, const CovalentTopology& bonds) {
    for (size_t bi : bonds.BondIndicesFor(a)) {
        const Bond& bond = bonds.BondAt(bi);
        if ((bond.atom_index_a == a && bond.atom_index_b == b) ||
            (bond.atom_index_a == b && bond.atom_index_b == a)) {
            return true;
        }
    }
    return false;
}

size_t FindWithLabel(const std::vector<size_t>& candidates,
                     const std::vector<AtomSemanticTable>& sem,
                     RingSystemKind kind,
                     RingPositionLabel label) {
    for (size_t ai : candidates) {
        if (LabelForKind(sem[ai], kind) == label) return ai;
    }
    return SIZE_MAX;
}

[[noreturn]] void FatalSubstrateGap(const char* what,
                                    RingSystemKind kind,
                                    int residue_seq) {
    std::fprintf(stderr,
        "FATAL: RingTopology::ConstructFromSubstrate: %s for ring "
        "kind=%d at residue seq %d. Substrate is required to provide "
        "complete labelling for ring construction; no string fallback.\n",
        what, static_cast<int>(kind), residue_seq);
    std::abort();
}

// PHE benzene / TYR phenol walk:
//   Ipso (CG) -> Ortho1 (CD1) -> Meta1 (CE1) -> Para (CZ) ->
//   Meta2 (CE2) -> Ortho2 (CD2)
std::vector<size_t> WalkSixRingByLabel(
        const std::vector<size_t>& candidates,
        const std::vector<AtomSemanticTable>& sem,
        RingSystemKind kind,
        int residue_seq) {

    const std::array<RingPositionLabel, 6> walk_labels = {
        RingPositionLabel::Ipso,
        RingPositionLabel::Ortho1,
        RingPositionLabel::Meta1,
        RingPositionLabel::Para,
        RingPositionLabel::Meta2,
        RingPositionLabel::Ortho2,
    };

    std::vector<size_t> walk;
    walk.reserve(6);
    for (RingPositionLabel label : walk_labels) {
        const size_t ai = FindWithLabel(candidates, sem, kind, label);
        if (ai == SIZE_MAX) {
            FatalSubstrateGap("six-ring missing canonical label",
                              kind, residue_seq);
        }
        walk.push_back(ai);
    }
    return walk;
}

// HIS imidazole walk:
//   Ipso (CG) -> Heteroatom-Locant-Delta (ND1) -> PyrroleAlpha (CE1)
//   -> Heteroatom-Locant-Epsilon (NE2) -> PyrroleBeta (CD2)
// The two Heteroatom_* atoms (ND1, NE2) are disambiguated by Locant.
// This is variant-stable: the substrate generator emits Locant::Delta
// for ND1 and Locant::Epsilon for NE2 in HID/HIE/HIP alike; only the
// Heteroatom_NH/_NoH label flips between them across variants.
std::vector<size_t> WalkHisImidazole(
        const std::vector<size_t>& candidates,
        const std::vector<AtomSemanticTable>& sem,
        int residue_seq) {

    const RingSystemKind kind = RingSystemKind::Imidazole_His;

    auto find_heteroatom_with_locant = [&](Locant loc) -> size_t {
        for (size_t ai : candidates) {
            const RingPositionLabel pos = LabelForKind(sem[ai], kind);
            if ((pos == RingPositionLabel::Heteroatom_NH ||
                 pos == RingPositionLabel::Heteroatom_NoH) &&
                sem[ai].locant == loc) {
                return ai;
            }
        }
        return SIZE_MAX;
    };

    const size_t ipso      = FindWithLabel(candidates, sem, kind,
                                           RingPositionLabel::Ipso);
    const size_t delta_n   = find_heteroatom_with_locant(Locant::Delta);
    const size_t alpha     = FindWithLabel(candidates, sem, kind,
                                           RingPositionLabel::PyrroleAlpha);
    const size_t epsilon_n = find_heteroatom_with_locant(Locant::Epsilon);
    const size_t beta      = FindWithLabel(candidates, sem, kind,
                                           RingPositionLabel::PyrroleBeta);

    if (ipso == SIZE_MAX || delta_n == SIZE_MAX || alpha == SIZE_MAX ||
        epsilon_n == SIZE_MAX || beta == SIZE_MAX) {
        FatalSubstrateGap("HIS imidazole walk incomplete",
                          kind, residue_seq);
    }

    return {ipso, delta_n, alpha, epsilon_n, beta};
}

// TRP-5 (pyrrole) walk:
//   Ipso (CG) -> PyrroleBeta (CD1) -> Heteroatom_NH (NE1) ->
//   BridgeFusion-bonded-to-NE1 (CE2) ->
//   BridgeFusion-bonded-to-CG  (CD2)
// The two BridgeFusion atoms (CE2 and CD2) share the same label.
// Disambiguation by bond-graph adjacency: CE2 is the BridgeFusion
// bonded to NE1 (the Heteroatom_NH); CD2 is the BridgeFusion bonded
// to CG (the Ipso). The 5-ring's outer edge runs CG-CD1-NE1-CE2-CD2-CG;
// the bridge bond CE2-CD2 is the shared edge with the 6-ring.
std::vector<size_t> WalkIndolePyrrole(
        const std::vector<size_t>& candidates,
        const std::vector<AtomSemanticTable>& sem,
        const CovalentTopology& bonds,
        int residue_seq) {

    const RingSystemKind kind = RingSystemKind::Indole_Trp_5;

    const size_t ipso = FindWithLabel(candidates, sem, kind,
                                      RingPositionLabel::Ipso);
    const size_t beta = FindWithLabel(candidates, sem, kind,
                                      RingPositionLabel::PyrroleBeta);
    const size_t het  = FindWithLabel(candidates, sem, kind,
                                      RingPositionLabel::Heteroatom_NH);

    if (ipso == SIZE_MAX || beta == SIZE_MAX || het == SIZE_MAX) {
        FatalSubstrateGap("TRP-5 missing Ipso / PyrroleBeta / "
                          "Heteroatom_NH", kind, residue_seq);
    }

    size_t bridge_adj_het = SIZE_MAX;
    size_t bridge_adj_ipso = SIZE_MAX;
    for (size_t ai : candidates) {
        if (LabelForKind(sem[ai], kind) !=
            RingPositionLabel::BridgeFusion) continue;
        const bool to_het  = AreBonded(ai, het,  bonds);
        const bool to_ipso = AreBonded(ai, ipso, bonds);
        if (to_het && !to_ipso) {
            bridge_adj_het = ai;
        } else if (to_ipso && !to_het) {
            bridge_adj_ipso = ai;
        }
    }

    if (bridge_adj_het == SIZE_MAX || bridge_adj_ipso == SIZE_MAX) {
        FatalSubstrateGap("TRP-5 BridgeFusion disambiguation failed",
                          kind, residue_seq);
    }

    return {ipso, beta, het, bridge_adj_het, bridge_adj_ipso};
}

// TRP-6 (benzene) walk:
//   BridgeFusion-bonded-to-Ortho1 (CD2) ->
//   BridgeFusion-bonded-to-Ortho2 (CE2) ->
//   Ortho2 (CZ2) -> Meta2 (CH2) -> Meta1 (CZ3) -> Ortho1 (CE3)
// Two BridgeFusion atoms; disambiguation by bond-graph adjacency.
// CD2 sits between CE2 (across the bridge) and CE3 (the Ortho1
// neighbour); CE2 sits between CD2 and CZ2 (the Ortho2 neighbour).
std::vector<size_t> WalkIndoleBenzene(
        const std::vector<size_t>& candidates,
        const std::vector<AtomSemanticTable>& sem,
        const CovalentTopology& bonds,
        int residue_seq) {

    const RingSystemKind kind = RingSystemKind::Indole_Trp_6;

    const size_t ortho1 = FindWithLabel(candidates, sem, kind,
                                        RingPositionLabel::Ortho1);
    const size_t ortho2 = FindWithLabel(candidates, sem, kind,
                                        RingPositionLabel::Ortho2);
    const size_t meta1  = FindWithLabel(candidates, sem, kind,
                                        RingPositionLabel::Meta1);
    const size_t meta2  = FindWithLabel(candidates, sem, kind,
                                        RingPositionLabel::Meta2);

    if (ortho1 == SIZE_MAX || ortho2 == SIZE_MAX ||
        meta1 == SIZE_MAX || meta2 == SIZE_MAX) {
        FatalSubstrateGap("TRP-6 missing Ortho/Meta labels",
                          kind, residue_seq);
    }

    size_t bridge_adj_ortho1 = SIZE_MAX;
    size_t bridge_adj_ortho2 = SIZE_MAX;
    for (size_t ai : candidates) {
        if (LabelForKind(sem[ai], kind) !=
            RingPositionLabel::BridgeFusion) continue;
        const bool to_o1 = AreBonded(ai, ortho1, bonds);
        const bool to_o2 = AreBonded(ai, ortho2, bonds);
        if (to_o1 && !to_o2) {
            bridge_adj_ortho1 = ai;
        } else if (to_o2 && !to_o1) {
            bridge_adj_ortho2 = ai;
        }
    }

    if (bridge_adj_ortho1 == SIZE_MAX ||
        bridge_adj_ortho2 == SIZE_MAX) {
        FatalSubstrateGap("TRP-6 BridgeFusion disambiguation failed",
                          kind, residue_seq);
    }

    return {bridge_adj_ortho1, bridge_adj_ortho2,
            ortho2, meta2, meta1, ortho1};
}

// TRP-9 (indole 9-atom perimeter) walk:
//   5/Ipso (CG) -> 5/PyrroleBeta (CD1) -> 5/Heteroatom_NH (NE1) ->
//   5/BridgeFusion-bonded-to-NE1 (CE2) -> 6/Ortho2 (CZ2) ->
//   6/Meta2 (CH2) -> 6/Meta1 (CZ3) -> 6/Ortho1 (CE3) ->
//   5/BridgeFusion-bonded-to-CG (CD2)
//
// The tertiary slot uniformly carries (Indole_Trp_9, PerimeterMember)
// for all nine atoms — that label set is information-free for
// walking. The walk uses each atom's PRIMARY-slot label (which the
// 5-ring atoms own as 5/* and the non-bridge 6-ring atoms own as 6/*),
// with bond-graph adjacency disambiguating the two BridgeFusion atoms.
// Bridge atoms are selected through their 5-ring BridgeFusion labels.
// Case 1995, J. Biomol. NMR 6, 341-346 is the citation for the
// 9-atom indole pi current circuit.
std::vector<size_t> WalkIndolePerimeter(
        const std::vector<size_t>& candidates,
        const std::vector<AtomSemanticTable>& sem,
        const CovalentTopology& bonds,
        int residue_seq) {

    const RingSystemKind k5 = RingSystemKind::Indole_Trp_5;
    const RingSystemKind k6 = RingSystemKind::Indole_Trp_6;

    const size_t ipso   = FindWithLabel(candidates, sem, k5,
                                        RingPositionLabel::Ipso);
    const size_t pbeta  = FindWithLabel(candidates, sem, k5,
                                        RingPositionLabel::PyrroleBeta);
    const size_t het    = FindWithLabel(candidates, sem, k5,
                                        RingPositionLabel::Heteroatom_NH);
    const size_t ortho1 = FindWithLabel(candidates, sem, k6,
                                        RingPositionLabel::Ortho1);
    const size_t ortho2 = FindWithLabel(candidates, sem, k6,
                                        RingPositionLabel::Ortho2);
    const size_t meta1  = FindWithLabel(candidates, sem, k6,
                                        RingPositionLabel::Meta1);
    const size_t meta2  = FindWithLabel(candidates, sem, k6,
                                        RingPositionLabel::Meta2);

    if (ipso == SIZE_MAX || pbeta == SIZE_MAX || het == SIZE_MAX ||
        ortho1 == SIZE_MAX || ortho2 == SIZE_MAX ||
        meta1 == SIZE_MAX || meta2 == SIZE_MAX) {
        FatalSubstrateGap("TRP-9 perimeter substrate incomplete",
                          RingSystemKind::Indole_Trp_9, residue_seq);
    }

    size_t bridge_adj_het = SIZE_MAX;
    size_t bridge_adj_ipso = SIZE_MAX;
    for (size_t ai : candidates) {
        if (LabelForKind(sem[ai], k5) !=
            RingPositionLabel::BridgeFusion) continue;
        const bool to_het  = AreBonded(ai, het,  bonds);
        const bool to_ipso = AreBonded(ai, ipso, bonds);
        if (to_het && !to_ipso) {
            bridge_adj_het = ai;
        } else if (to_ipso && !to_het) {
            bridge_adj_ipso = ai;
        }
    }

    if (bridge_adj_het == SIZE_MAX || bridge_adj_ipso == SIZE_MAX) {
        FatalSubstrateGap("TRP-9 BridgeFusion disambiguation failed",
                          RingSystemKind::Indole_Trp_9, residue_seq);
    }

    return {ipso, pbeta, het, bridge_adj_het,
            ortho2, meta2, meta1, ortho1,
            bridge_adj_ipso};
}

// Pro pyrrolidine walk:
//   ProRingNitrogen (N) -> ProRingAlphaCarbon (CA) -> ProRingBeta (CB)
//   -> ProRingPuckerPivot (CG) -> ProRingDelta (CD)
// ProPyrrolidineRing::Intensity() returns 0.0 because this is a
// saturated heterocycle.
std::vector<size_t> WalkProPyrrolidine(
        const std::vector<size_t>& candidates,
        const std::vector<AtomSemanticTable>& sem,
        int residue_seq) {

    const RingSystemKind kind = RingSystemKind::Pyrrolidine_Pro;

    const std::array<RingPositionLabel, 5> walk_labels = {
        RingPositionLabel::ProRingNitrogen,
        RingPositionLabel::ProRingAlphaCarbon,
        RingPositionLabel::ProRingBeta,
        RingPositionLabel::ProRingPuckerPivot,
        RingPositionLabel::ProRingDelta,
    };

    std::vector<size_t> walk;
    walk.reserve(5);
    for (RingPositionLabel label : walk_labels) {
        const size_t ai = FindWithLabel(candidates, sem, kind, label);
        if (ai == SIZE_MAX) {
            FatalSubstrateGap("Pro pyrrolidine walk incomplete",
                              kind, residue_seq);
        }
        walk.push_back(ai);
    }
    return walk;
}

// HIS and Pro are intentionally outside this mapping.
RingTypeIndex AromaticTypeFromKind(RingSystemKind kind) {
    switch (kind) {
        case RingSystemKind::Benzene_Phe:    return RingTypeIndex::PheBenzene;
        case RingSystemKind::Benzene_Tyr:    return RingTypeIndex::TyrPhenol;
        case RingSystemKind::Indole_Trp_5:   return RingTypeIndex::TrpPyrrole;
        case RingSystemKind::Indole_Trp_6:   return RingTypeIndex::TrpBenzene;
        case RingSystemKind::Indole_Trp_9:   return RingTypeIndex::TrpPerimeter;
        case RingSystemKind::Imidazole_His:
        case RingSystemKind::Pyrrolidine_Pro:
        case RingSystemKind::NotInRing:
            break;
    }
    std::fprintf(stderr,
        "FATAL: AromaticTypeFromKind: kind %d is not an aromatic "
        "ring system; HIS dispatches via variant_index, Pro via "
        "the saturated path.\n", static_cast<int>(kind));
    std::abort();
}

// variant_index: 0=HID, 1=HIE, 2=HIP. Unresolved and HIP both use the
// ambiguous HisImidazole type.
RingTypeIndex HisRingTypeFromVariant(int variant_index) {
    switch (variant_index) {
        case 0: return RingTypeIndex::HidImidazole;
        case 1: return RingTypeIndex::HieImidazole;
        case 2: return RingTypeIndex::HisImidazole;
        default: return RingTypeIndex::HisImidazole;
    }
}

}  // anonymous namespace

std::vector<size_t> RingTopology::CanonicalCyclicWalk(
        RingSystemKind kind,
        const std::vector<size_t>& candidates,
        const std::vector<AtomSemanticTable>& atom_semantic,
        const CovalentTopology& bonds) {

    constexpr int residue_seq = -1;

    switch (kind) {
        case RingSystemKind::Benzene_Phe:
        case RingSystemKind::Benzene_Tyr:
            return WalkSixRingByLabel(candidates, atom_semantic,
                                      kind, residue_seq);
        case RingSystemKind::Imidazole_His:
            return WalkHisImidazole(candidates, atom_semantic,
                                    residue_seq);
        case RingSystemKind::Indole_Trp_5:
            return WalkIndolePyrrole(candidates, atom_semantic,
                                     bonds, residue_seq);
        case RingSystemKind::Indole_Trp_6:
            return WalkIndoleBenzene(candidates, atom_semantic,
                                     bonds, residue_seq);
        case RingSystemKind::Indole_Trp_9:
            return WalkIndolePerimeter(candidates, atom_semantic,
                                       bonds, residue_seq);
        case RingSystemKind::Pyrrolidine_Pro:
            return WalkProPyrrolidine(candidates, atom_semantic,
                                      residue_seq);
        case RingSystemKind::NotInRing:
            break;
    }
    std::fprintf(stderr,
        "FATAL: RingTopology::CanonicalCyclicWalk: unhandled "
        "RingSystemKind=%d.\n", static_cast<int>(kind));
    std::abort();
}


std::unique_ptr<RingTopology> RingTopology::ConstructFromSubstrate(
        const std::vector<Residue>& residues,
        const std::vector<AtomSemanticTable>& atom_semantic,
        const CovalentTopology& bonds) {

    auto topo = std::unique_ptr<RingTopology>(new RingTopology());

    if (atom_semantic.empty()) return topo;

    // Track each TRP residue's aromatic_-vector indices for the
    // benzene <-> pyrrole fused-partner post-pass. The 9-atom
    // perimeter is not linked because fused_partner_index is binary.
    struct TrpRingIndices {
        size_t benzene_idx = SIZE_MAX;
        size_t pyrrole_idx = SIZE_MAX;
    };
    std::map<size_t, TrpRingIndices> trp_rings;

    auto is_heavy = [&](size_t ai) {
        return atom_semantic[ai].element != Element::H;
    };

    for (size_t ri = 0; ri < residues.size(); ++ri) {
        const Residue& res = residues[ri];

        // RingTypeIndex order controls output ring order. For TRP it
        // gives benzene, pyrrole, then 9-atom perimeter; enum order
        // would put the 5-ring before the 6-ring.
        std::vector<std::pair<RingSystemKind, RingTypeIndex>> kinds_with_type;
        std::set<RingSystemKind> seen;
        auto note_kind = [&](RingSystemKind kind) {
            if (kind == RingSystemKind::NotInRing) return;
            if (seen.count(kind) > 0) return;
            seen.insert(kind);

            RingTypeIndex type;
            if (kind == RingSystemKind::Imidazole_His) {
                type = HisRingTypeFromVariant(
                    res.protonation_variant_index);
            } else if (kind == RingSystemKind::Pyrrolidine_Pro) {
                type = RingTypeIndex::ProPyrrolidine;
            } else {
                type = AromaticTypeFromKind(kind);
            }
            kinds_with_type.push_back({kind, type});
        };

        for (size_t ai : res.atom_indices) {
            if (!is_heavy(ai)) continue;
            const RingPosition& rp = atom_semantic[ai].ring_position;
            note_kind(rp.primary.ring);
            note_kind(rp.secondary.ring);
            note_kind(rp.tertiary.ring);
        }
        if (kinds_with_type.empty()) continue;

        std::sort(kinds_with_type.begin(), kinds_with_type.end(),
            [](const auto& a, const auto& b) {
                return static_cast<int>(a.second) <
                       static_cast<int>(b.second);
            });

        for (const auto& [kind, type] : kinds_with_type) {
            std::vector<size_t> heavy;
            for (size_t ai : res.atom_indices) {
                if (!is_heavy(ai)) continue;
                if (LabelForKind(atom_semantic[ai], kind) !=
                    RingPositionLabel::NotInRing) {
                    heavy.push_back(ai);
                }
            }
            if (heavy.empty()) continue;

            std::vector<size_t> walk = CanonicalCyclicWalk(
                kind, heavy, atom_semantic, bonds);

            auto ring = CreateRing(type);
            ring->atom_indices = std::move(walk);
            ring->parent_residue_index = ri;
            ring->parent_residue_number = res.sequence_number;

            if (type == RingTypeIndex::ProPyrrolidine) {
                topo->saturated_.push_back(std::move(ring));
            } else {
                const size_t aromatic_idx = topo->aromatic_.size();
                if (res.type == AminoAcid::TRP) {
                    auto& trp = trp_rings[ri];
                    if (type == RingTypeIndex::TrpBenzene) {
                        trp.benzene_idx = aromatic_idx;
                    } else if (type == RingTypeIndex::TrpPyrrole) {
                        trp.pyrrole_idx = aromatic_idx;
                    }
                }
                topo->aromatic_.push_back(std::move(ring));
            }
        }
    }

    for (const auto& kv : trp_rings) {
        const TrpRingIndices& trp = kv.second;
        if (trp.benzene_idx != SIZE_MAX &&
            trp.pyrrole_idx != SIZE_MAX) {
            topo->aromatic_[trp.benzene_idx]->fused_partner_index =
                trp.pyrrole_idx;
            topo->aromatic_[trp.pyrrole_idx]->fused_partner_index =
                trp.benzene_idx;
        }
    }

    return topo;
}

}  // namespace nmr

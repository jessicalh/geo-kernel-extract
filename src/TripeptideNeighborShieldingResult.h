#pragma once
//
// TripeptideNeighborShieldingResult: Δσ_BB^{i±1} per Larsen 2015 Eq 3.
//
// Larsen 2015 Eq 3 (read directly from
// references/larsen-2015-procs15-dft-chemical-shift-predictor.pdf p.3):
//
//   Δσ_BB^{i-1}(residue i)
//     = σ_BB^{i-1}(φ^{i-1}, ψ^{i-1}, χ_1^{i-1}, χ_2^{i-1}, …)
//       − σ^A(φ_std, ψ_std)
//
//   with φ_std = -120°, ψ_std = 140° for all φ and ψ in the AAA
//   reference tripeptide.
//
// Critical reading per the Cβ/Val example in the paper:
//
//   "if residue i is a Ser and residue i − 1 is a Val then the effect
//    of the Val side-chain on the Cβ chemical shielding of the Ser
//    residue is computed as the difference in the chemical shielding
//    of the Cβ atom in the **C-terminal Ala residue** computed for an
//    AVA and AAA tripeptide."
//
// So σ_BB^{i-1} is the σ at the C-terminal ALA atoms of the AXA
// tripeptide where X = residue (i-1). The C-terminal ALA represents
// residue (i-1)+1 = i in the protein. The AAA reference at standard
// angles is queried once. The DIFFERENCE in σ at the C-terminal ALA
// atoms between AXA and AAA is the (i-1) sidechain's effect on
// residue i.
//
// Symmetric for Δσ_BB^{i+1}: query AYA where Y = (i+1), read σ at the
// N-TERMINAL ALA atoms (which represents residue (i+1)-1 = i).
//
// Algorithm:
//
//   1. AAA reference: query AAA at (φ=-120, ψ=140), one-shot cached.
//   2. For each residue i in the protein:
//      a. For δ ∈ {-1, +1} where (i+δ) exists in the same chain:
//         i.   Get (i+δ)'s identity, actual φ/ψ/χ.
//         ii.  Query DB: A-(i+δ)-A at (i+δ)'s actual angles.
//         iii. Identify the flanking ALA cap of (i+δ)'s tripeptide:
//                  δ = -1: C-terminal ALA (atoms after central, before NME)
//                  δ = +1: N-terminal ALA (atoms after ACE, before central)
//         iv.  Identify the SAME flanking ALA in the AAA reference.
//         v.   Kabsch align: flanking ALA's N/CA/C → protein i's N/CA/C
//              in BOTH the (i+δ) tripeptide and AAA reference → R_(i+δ),
//              R_AAA.
//         vi.  For each flanking ALA atom k:
//                 σ_AXA_k → R_(i+δ) σ_AXA_k R_(i+δ)^T (aligned to protein)
//                 σ_AAA_k → R_AAA  σ_AAA_k R_AAA^T  (aligned to protein)
//                 Δσ_k = σ_AXA_k - σ_AAA_k
//              Match the flanking ALA atom k to a protein atom in
//              residue i by element + nearest distance (post-rotation),
//              accumulate Δσ_k on that protein atom.
//
// Larsen's assumption: the effect of i±1's side-chain on residue i
// is dominated by the chemistry change at the flanking ALA position,
// and the φ/ψ/χ of residue i itself plus the nature of residue i are
// secondary. This is the "AXA-scan reuse" trick — no new DFT data
// needed beyond what's already in tensorcs15.
//
// Per-atom storage on ConformationAtom (defined in ConformationAtom.h):
//   tripeptide_neighbor_shielding_tensor       Mat3 ppm
//   tripeptide_neighbor_shielding_spherical    SphericalTensor ppm
//   tripeptide_neighbor_has_match              bool
//
// Sum of (i-1) + (i+1) contributions is what's stored — Larsen's
// six-term decomposition treats them as independent additive terms.
// Per-side breakdown can be derived in a separate calculator if a
// downstream consumer needs it.
//
// frame_type discriminator: both AXA and AAA rows are typically
// gaussian_standard_orientation (OPBE). If the flanking residue is
// SER (i-1 = S), the AXA query returns an ASA row with
// frame_type=orca_input_orientation. Δσ then mixes OPBE (AAA) and
// PBE (ASA) — methods caveat applies (project_serine_pbe_discontinuity).
//

#include "ConformationResult.h"
#include "Types.h"
#include "TripeptideDftTable.h"

#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class ProteinConformation;


class TripeptideNeighborShieldingResult : public ConformationResult {
public:
    std::string Name() const override {
        return "TripeptideNeighborShieldingResult";
    }

    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<TripeptideNeighborShieldingResult> Compute(
        ProteinConformation& conf,
        const TripeptideDftTable& table);

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

    // Per-residue diagnostics (length = ResidueCount()).
    struct ResidueMatch {
        // i-1 contribution
        int    prev_calc_id        = 0;       // 0 = miss
        double prev_backbone_rmsd  = 0.0;     // Å
        std::string prev_frame_type;
        int    prev_n_atoms_matched = 0;
        // i+1 contribution
        int    next_calc_id        = 0;
        double next_backbone_rmsd  = 0.0;
        std::string next_frame_type;
        int    next_n_atoms_matched = 0;
    };
    const std::vector<ResidueMatch>& ResidueMatches() const {
        return residue_matches_;
    }

    int ResiduesWithAnyNeighbor() const { return residues_any_; }
    int AtomsAccumulated()        const { return atoms_accumulated_; }

private:
    const ProteinConformation* conf_ = nullptr;
    std::vector<ResidueMatch>  residue_matches_;
    int residues_any_       = 0;
    int atoms_accumulated_  = 0;
};


}  // namespace nmr

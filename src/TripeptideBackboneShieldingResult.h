#pragma once
//
// TripeptideBackboneShieldingResult: per-atom DFT shielding from
// ProCS15 tripeptide DB lookup, aligned to the protein conformation.
//
// Algorithm (port of gotham TripeptideAssembler with substrate gates +
// chi3/chi4 extension):
//
//   For each protein residue r with valid backbone N/CA/C cache and
//   well-defined phi/psi:
//     1. Compute φ, ψ, χ₁..χₙ from cached atom indices on the
//        ProteinConformation.
//     2. Query TripeptideDftTable.QueryNearest(letter, φ, ψ, χ₁..χ₄)
//        for the matching DFT pose. n_chi is set per residue type.
//     3. On a chi-specific miss, fall back to fewer chi axes
//        (drop χ₄, then χ₃, etc.) until the WHERE clause matches.
//     4. Kabsch-align tripeptide central N/CA/C onto the protein's
//        N/CA/C → R, t.
//     5. Apply the same R to the σ tensors (R σ Rᵀ) at every central-
//        residue atom.
//     6. Apply a sidechain re-rotation around the CA-CB axis
//        (Rodrigues with the angle that takes the aligned trip CB to
//        the protein CB direction) — fixes the 20° grid coarseness on
//        χ₁ while preserving the backbone alignment.
//     7. Match each rotated tripeptide central atom to the nearest
//        protein atom of the same element within 5 Å, take the
//        rotated σ tensor, store on ConformationAtom.
//
// Substrate gates: BackboneRole and matched-atom-element checks via
// LegacyAmberTopology / Protein::AtomAt(i).element. NO string
// comparisons on atom names. Per PATTERNS.md anti-string-traversal.
//
// Per-atom storage on ConformationAtom (T2 sacred):
//   tripeptide_bb_shielding_tensor       Mat3 ppm
//   tripeptide_bb_shielding_spherical    SphericalTensor ppm
//   tripeptide_bb_match_distance         Å
//   tripeptide_bb_has_match              bool
//   tripeptide_bb_method_tag             0|1|2 (frame_type provenance)
//
// Per-residue storage on the result (diagnostics + ML features):
//   residue_match_calc_id                int
//   residue_match_backbone_rmsd          Å
//   residue_match_frame_type             "" or "gaussian_…" / "orca_…"
//   residue_phi_actual / psi_actual      degrees
//   residue_phi_db / psi_db              degrees (rounded grid)
//   residue_chi_db                       per-chi grid points
//

#include "ConformationResult.h"
#include "Types.h"
#include "TripeptideDftTable.h"

#include <array>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class ProteinConformation;


class TripeptideBackboneShieldingResult : public ConformationResult {
public:
    std::string Name() const override {
        return "TripeptideBackboneShieldingResult";
    }

    std::vector<std::type_index> Dependencies() const override;

    // Factory. Returns nullptr only on hard structural errors (zero
    // atoms, table not connected); per-residue misses are handled
    // internally with NaN / has_match=false on the per-atom record.
    static std::unique_ptr<TripeptideBackboneShieldingResult> Compute(
        ProteinConformation& conf,
        const TripeptideDftTable& table);

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

    // ── Per-residue diagnostics (length = ResidueCount()) ───────────

    struct ResidueMatch {
        int    calc_id        = 0;        // 0 = miss
        double backbone_rmsd  = 0.0;      // Å, post-Kabsch
        double ca_match_dist  = 0.0;      // Å, CA matched-atom distance
        std::string frame_type;            // discriminator from DB row
        // Actual (protein) angles
        double phi_actual = 0.0, psi_actual = 0.0;
        double chi_actual[4] = {0.0, 0.0, 0.0, 0.0};
        // DB-grid-rounded angles
        int phi_db = 0, psi_db = 0;
        int chi_db[4] = {0, 0, 0, 0};
        // Atom-mapping stats
        int n_atoms_matched   = 0;
        int n_chi_axes_used   = 0;        // 0..4 — how deep the lookup went
    };

    const std::vector<ResidueMatch>& ResidueMatches() const {
        return residue_matches_;
    }

    // Aggregate stats over the matched residues.
    int ResiduesAttempted() const { return residues_attempted_; }
    int ResiduesMatched()   const { return residues_matched_; }
    int ResiduesFailed()    const { return residues_failed_; }
    int AtomsAssigned()     const { return atoms_assigned_; }
    double MeanBackboneRmsd() const { return mean_backbone_rmsd_; }
    double MaxBackboneRmsd()  const { return max_backbone_rmsd_; }

private:
    const ProteinConformation* conf_ = nullptr;

    std::vector<ResidueMatch> residue_matches_;

    int    residues_attempted_   = 0;
    int    residues_matched_     = 0;
    int    residues_failed_      = 0;
    int    atoms_assigned_       = 0;
    double mean_backbone_rmsd_   = 0.0;
    double max_backbone_rmsd_    = 0.0;
};


}  // namespace nmr

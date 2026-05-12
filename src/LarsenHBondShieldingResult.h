#pragma once
//
// LarsenHBondShieldingResult: per-atom Larsen 2015 ProCS15 H-bond term
// shielding contributions, computed via direct DFT grid lookup against
// LarsenHBondGrid.
//
// Implements Larsen 2015 Eq 5 four-term decomposition:
//
//   Δσ_HB^i  = Δσ_1°HB(rHO, θ, ρ)   + Δσ_2°HB(rOH, θO, ρO)     (amide H donor)
//   Δσ_HαB^i = Δσ_1°HαB(rHαO, θ, ρ) + Δσ_2°HαB(rOHα, θO, ρO)   (Hα donor)
//
// Plus the water term Δσ_w = 2.07 ppm isotropic on amide H atoms that
// have no geometric H-bond candidate at all (solvent-exposed amides
// per Larsen's NMA-water complex DFT).
//
// Methods accumulate (feedback_methods_accumulate): this calculator
// runs side-by-side with the kernel-form HBondResult. Both emit their
// own NPYs. The kernel-vs-grid per-atom-type residual is itself a
// methodological coordinate.
//
// Enumeration: spatial-search-driven over SpatialIndexResult. For each
// donor candidate (substrate-typed amide H or any α-hydrogen via
// IsAnyAlphaHydrogen — GLY HA2 and HA3 enumerate as separate donors),
// find candidate acceptor O atoms within kSpatialCutoff_A. Each
// candidate O is classified (Backbone / SidechainCarbonyl / Hydroxyl /
// Carboxylate) via ClassifyAcceptor; the (donor_class, acceptor_class)
// pair routes to one of 6 Larsen grids. Larsen's framework is
// geometric (not DSSP-energy-based), and the spatial sweep IS the
// H-bond finder.
//
// Per-atom-type contribution dispatch follows Larsen 2015 Table 2,
// encoded as `LarsenContribDispatch::Applies` below. Cβ row is
// intentionally all-false per Larsen — the calculator still emits a
// diagnostic Cβ tensor (should be near-zero) to verify the parser →
// loader → rotation pipeline produces what the physics predicts.
//
// Output per-atom fields land on ConformationAtom: larsen_hbond_*.
// NPYs emitted by WriteFeatures.
//
// Structural deviation from PATTERNS.md §6: grid-lookup calculators
// like this one have no separable kernel × parameter stage — the DFT
// grid output IS the calibrated ppm shielding directly. Pattern 11
// (both Mat3 AND SphericalTensor stored) is satisfied at every
// per-class field on ConformationAtom.
//
// Sign verification (PATTERNS.md §22): transitive via the
// LarsenHBondGrid loader test `TightHBondGeometryGivesNonzeroHα`,
// which locks the loader's sign at a known geometry. End-to-end
// calibrated-sign test against Larsen's published 1UBQ .procs / DFT
// log lands in the reality-check pass alongside SDK comparison.

#include "ConformationResult.h"
#include "LarsenHBondGrid.h"
#include "Types.h"

#include <cstdint>
#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class ProteinConformation;


// ============================================================================
// LarsenContribDispatch — Larsen 2015 Table 2 encoded as a constexpr
// ============================================================================
//
// Per-(target_atom_type, contribution_term) "applies/doesn't apply".
// One cell per row; one cell per column. A wrong cell silently zeroes
// or doubles an entire contribution class.
//
// Cβ row is INTENTIONALLY all-false per Larsen — we emit Cβ tensors
// as a diagnostic. A non-zero Cβ contribution in production is a
// methodological signal worth reporting (feedback_methods_accumulate).
//
namespace LarsenContribDispatch {

enum class TargetAtom : std::uint8_t {
    N = 0, CA, CB, C, HA, HN, Count
};

enum class Term : std::uint8_t {
    Primary_HB = 0,   // 1° amide-H donor effect, applies to donor residue i
    Secondary_HB,     // 2° amide-H donor effect, applies to acceptor residue i+1
    Primary_HaB,      // 1° Hα donor effect, applies to donor residue i
    Secondary_HaB,    // 2° Hα donor effect, applies to acceptor residue i+1
    RingCurrent,      // Δσ_RC, applied separately (existing calculators)
    Water,            // Δσ_w = 2.07 ppm, applied separately below
    Count
};

// Table 2 from Larsen 2015 PeerJ (DOI 10.7717/peerj.1344):
//                  1HB     2HB    1HaB   2HaB    RC     w
// N                ✓       ✓      ·      ✓       ·      ·
// Cα               ·       ✓      ·      ·       ·      ·
// Cβ               ·       ·      ·      ·       ·      ·   (all-zero diagnostic)
// C'               ·       ✓      ·      ·       ·      ·
// Hα               ✓       ✓      ✓      ✓       ✓      ·
// HN               ✓       ✓      ✓      ✓       ✓      ✓
constexpr bool Applies(TargetAtom t, Term term) {
    constexpr bool table[6][6] = {
        //  1HB    2HB    1HaB   2HaB   RC     w
        {  true,  true, false,  true, false, false }, // N
        { false,  true, false, false, false, false }, // CA
        { false, false, false, false, false, false }, // CB
        { false,  true, false, false, false, false }, // C
        {  true,  true,  true,  true,  true, false }, // HA
        {  true,  true,  true,  true,  true,  true }, // HN
    };
    return table[static_cast<int>(t)][static_cast<int>(term)];
}

}  // namespace LarsenContribDispatch


// ============================================================================
// LarsenHBondShieldingResult
// ============================================================================

class LarsenHBondShieldingResult : public ConformationResult {
public:
    std::string Name() const override {
        return "LarsenHBondShieldingResult";
    }

    std::vector<std::type_index> Dependencies() const override;

    // Factory. Returns nullptr only on hard structural errors (zero
    // atoms, grid not loaded). Per-pair classification failures are
    // logged via GeometryChoiceBuilder and skipped silently.
    static std::unique_ptr<LarsenHBondShieldingResult> Compute(
        ProteinConformation& conf,
        const LarsenHBondGrid& grid);

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

    // ── Per-pair diagnostic record (kept on the result for future
    // ── per-pair NPY emission and adversarial review).
    struct PairRecord {
        std::size_t        donor_atom_idx = 0;   // donor H atom
        std::size_t        acceptor_atom_idx = 0;  // acceptor O atom
        std::size_t        donor_residue_idx = 0;
        std::size_t        acceptor_residue_idx = 0;
        HBondDonorClass    donor_class    = HBondDonorClass::AmideHydrogen;
        HBondAcceptorClass acceptor_class = HBondAcceptorClass::BackboneCarbonyl;
        double             r_angstrom = 0.0;
        double             theta_deg  = 0.0;
        double             rho_deg    = 0.0;
        double             isotropic_total = 0.0;
        bool               any_corner_imputed = false;
    };
    const std::vector<PairRecord>& Pairs() const { return pairs_; }

    // Aggregate stats.
    //
    // PairsFound counts pairs the grid path successfully processed
    // (geometry computed, grid hit, tensors accumulated).
    // PairsGridSkipped counts geometric H-bond candidates the grid
    // path SKIPPED (out-of-range θ, grid miss, no i+1 mapping for the
    // 2° term on C-terminus acceptor). The two together sum to the
    // total geometric H-bond candidate count.
    // AmideHsUnboundWithWater counts amide Hs that received the
    // Δσ_w = 2.07 ppm term — gated on "ZERO geometric H-bond
    // candidates found." A grid-skipped pair does NOT trigger spurious
    // water-term assignment.
    int    PairsFound()             const { return static_cast<int>(pairs_.size()); }
    int    PairsGridSkipped()       const { return pairs_grid_skipped_; }
    int    AtomsWithContribution()  const { return atoms_with_contribution_; }
    int    AmideHsUnboundWithWater() const { return amide_hs_unbound_; }

private:
    const ProteinConformation* conf_ = nullptr;

    std::vector<PairRecord> pairs_;
    int atoms_with_contribution_ = 0;
    int amide_hs_unbound_ = 0;
    int pairs_grid_skipped_ = 0;
};


}  // namespace nmr

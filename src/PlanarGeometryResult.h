#pragma once
//
// PlanarGeometryResult: per-frame conformation companion to the
// LegacyAmber substrate's typed planar-group / ring-position fields.
//
// Optional downstream readers (RingPuckerTimeSeriesTrajectoryResult,
// DihedralTimeSeriesTrajectoryResult) read these getters per-frame;
// they attach conditionally, so absence is captured by the
// source-attached gate.
//
// The substrate carries the typed *classification* (PlanarGroupKind,
// PlanarStereo, RingPosition). This calculator carries the actual
// *deviation from canonical* in each frame.
//
// Per Amendment 2026-05-08(a) (spec/PLANNED_CALCULATORS_2026-04-22.md).
// Four quantities, all derived from positions only — no electronic
// structure work, negligible cost relative to any kernel calculator.
//
// Convention pins:
//
//   - **Peptide-bond planarity** — ω dihedral
//     Cα(i)–C(i)–N(successor)–Cα(successor) per residue, with
//     deviation from π.
//     ω and ω_deviation are emitted for every well-defined peptide
//     bond (covalent C(i)–N(successor) per `protein.BackboneSuccessor`),
//     INCLUDING X→Pro bonds. cis/trans isomerism at X-Pro is real
//     conformational signal — use the `omega_is_xpro` mask to flag
//     those rows for consumer-side interpretation, not NaN-fill.
//     NaN at: chain-break boundary (no covalent C(i)-N(successor)),
//     C-terminus, missing backbone-cache atoms. Convention: ω in
//     radians, deviation wrapped to [−π, π].
//
//   - **sp2 pyramidalization** — signed out-of-plane displacement of
//     each atom whose AtomSemanticTable::planar_group != None from
//     the plane of its three bonded neighbours. CHARMM/AMBER
//     improper-dihedral convention: sign by the right-hand rule on
//     the cross product of the first two sorted neighbour vectors
//     relative to the neighbour centroid (neighbours are sorted
//     ascending by atom index, so the sign is build-stable; see
//     ThreeBondedNeighbours). Per-atom field on ConformationAtom
//     (zero for non-planar atoms).
//
//   - **Aromatic χ₂** — per aromatic ring, the parent residue's χ₂
//     dihedral (Cα–Cβ–Cγ–Cδ for PHE/TYR; Cα–Cβ–Cγ–Nδ1 for HIS;
//     Cα–Cβ–Cγ–Cδ1 for TRP). Per Akke & Weininger 2023 J. Phys.
//     Chem. B 127, 591 (M17): Phe/Tyr ring flips are 180° rotamer
//     transitions of χ₂. **Important honesty caveat**: a single MD
//     frame gives the instantaneous χ₂ value; flip *kinetics* require
//     exchange-sensitive NMR rate measurements and are NOT measurable
//     from one frame. Field name `aromatic_chi2` reflects that — not
//     `ring_flip_state`.
//
//   - **Cremer-Pople pucker** — per saturated ring (Pro pyrrolidine
//     today; future saturated rings handled the same way), the
//     puckering amplitude Q (Å) and phase angle θ (degrees).
//     Cremer & Pople 1975 J. Am. Chem. Soc. 97, 1354 is the
//     canonical reference for the (Q, θ) parameters; for 5-rings,
//     θ mod 72° gives the envelope/twist endo/exo classification.
//     Atoms in canonical cyclic order via
//     RingTopology::CanonicalCyclicWalk (already at construction
//     time per the Bundle C / Slice B substrate).
//
// Lifecycle: ConformationResult subclass.
// Dependencies: GeometryResult (positions + ring geometry),
//               EnrichmentResult (substrate semantics finalised).
//

#include "ConformationResult.h"
#include "Types.h"

#include <memory>
#include <string>
#include <typeindex>
#include <vector>

namespace nmr {

class ProteinConformation;


class PlanarGeometryResult : public ConformationResult {
public:
    std::string Name() const override { return "PlanarGeometryResult"; }

    std::vector<std::type_index> Dependencies() const override;

    // Factory: compute all four quantities. Returns nullptr only on
    // structural errors (zero atoms, missing substrate semantics).
    // NaN-fill is the legitimate "not applicable" signal in the
    // per-residue and per-ring vectors.
    static std::unique_ptr<PlanarGeometryResult> Compute(
        ProteinConformation& conf);

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

    // ── Per-residue queries ────────────────────────────────────────
    const std::vector<double>& OmegaActual() const { return omega_actual_; }
    const std::vector<double>& OmegaDeviation() const { return omega_deviation_; }
    const std::vector<uint8_t>& OmegaIsXpro() const { return omega_is_xpro_; }

    // ── Per-ring queries ───────────────────────────────────────────
    const std::vector<double>& AromaticChi2() const { return aromatic_chi2_; }
    const std::vector<double>& PuckerQ() const { return pucker_Q_; }
    const std::vector<double>& PuckerTheta() const { return pucker_theta_; }

private:
    const ProteinConformation* conf_ = nullptr;

    // Per-residue: indexed by Protein residue index. NaN when there is
    // no backbone successor or either CA cache is missing. X→Pro bonds
    // are emitted with their actual ω value
    // (cis/trans isomerism is a real signal, not a deviation), and
    // the `omega_is_xpro_` mask flags those rows for the consumer to
    // interpret. Length = ResidueCount() for all three vectors.
    std::vector<double>  omega_actual_;
    std::vector<double>  omega_deviation_;
    std::vector<uint8_t> omega_is_xpro_;   // 1 if successor is Pro, else 0

    // Per aromatic ring: indexed by LegacyAmberTopology::AromaticRingAt
    // index. NaN if parent residue has no χ₂ defined. Length =
    // AromaticRingCount().
    std::vector<double> aromatic_chi2_;

    // Per saturated ring: indexed by LegacyAmberTopology::SaturatedRingAt
    // index. Q in Å, θ in degrees [0, 360). Length =
    // SaturatedRingCount().
    std::vector<double> pucker_Q_;
    std::vector<double> pucker_theta_;
};

}  // namespace nmr

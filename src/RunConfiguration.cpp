#include "RunConfiguration.h"
#include "TrajectoryProtein.h"
#include "TrajectoryResult.h"

// Per-frame ConformationResult types used by the named trajectory
// configurations. Include full types so type_index values are valid at
// factory call time.
#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "EnrichmentResult.h"
#include "DsspResult.h"
#include "BiotSavartResult.h"
#include "HaighMallionResult.h"
#include "McConnellResult.h"
#include "SasaResult.h"
#include "ChargeAssignmentResult.h"
#include "ApbsFieldResult.h"
#include "RingSusceptibilityResult.h"
#include "PiQuadrupoleResult.h"
#include "DispersionResult.h"
#include "HBondResult.h"
#include "EeqResult.h"
#include "AIMNet2Result.h"
#include "WaterFieldResult.h"
#include "HydrationShellResult.h"
#include "HydrationGeometryResult.h"
#include "GromacsEnergyResult.h"
#include "GromacsEnergyTimeSeriesTrajectoryResult.h"
#include "BondedEnergyTimeSeriesTrajectoryResult.h"
#include "WaterFieldTimeSeriesTrajectoryResult.h"
#include "WaterFieldWelfordTrajectoryResult.h"
#include "HydrationGeometryTimeSeriesTrajectoryResult.h"
#include "HydrationGeometryWelfordTrajectoryResult.h"
#include "HydrationShellTimeSeriesTrajectoryResult.h"
#include "HydrationShellWelfordTrajectoryResult.h"
#include "DihedralTimeSeriesTrajectoryResult.h"
#include "RingNeighbourhoodTrajectoryStats.h"
#include "RmsdTrackingTrajectoryResult.h"
#include "RmsdSpikeSelectionTrajectoryResult.h"
#include "DftPoseCoordinatorTrajectoryResult.h"
#include "ChiRotamerSelectionTrajectoryResult.h"
#include "DihedralBinTransitionTrajectoryResult.h"
#include "Dssp8TimeSeriesTrajectoryResult.h"
#include "Dssp8TransitionTrajectoryResult.h"
#include "RingPuckerTimeSeriesTrajectoryResult.h"
#include "JCouplingTimeSeriesTrajectoryResult.h"
#include "GromacsFramePullResult.h"
#include "BondedEnergyResult.h"
#include "TripeptideNeighborShieldingResult.h"

// Concrete TrajectoryResults populating the canonical configurations.
#include "BsWelfordTrajectoryResult.h"
#include "HmWelfordTrajectoryResult.h"
#include "McConnellWelfordTrajectoryResult.h"
#include "EeqWelfordTrajectoryResult.h"
#include "SasaWelfordTrajectoryResult.h"
#include "HBondCountWelfordTrajectoryResult.h"
#include "BsShieldingTimeSeriesTrajectoryResult.h"
#include "HmShieldingTimeSeriesTrajectoryResult.h"
#include "McConnellShieldingTimeSeriesTrajectoryResult.h"
#include "PiQuadrupoleShieldingTimeSeriesTrajectoryResult.h"
#include "RingSusceptibilityShieldingTimeSeriesTrajectoryResult.h"
#include "DispersionShieldingTimeSeriesTrajectoryResult.h"
#include "HBondShieldingTimeSeriesTrajectoryResult.h"
#include "SasaTimeSeriesTrajectoryResult.h"
#include "AIMNet2ChargeTimeSeriesTrajectoryResult.h"
#include "AIMNet2EmbeddingTimeSeriesTrajectoryResult.h"
#include "AIMNet2ChargeResponseGradientResult.h"
#include "AIMNet2ChargeResponseGradientTimeSeriesTrajectoryResult.h"
#include "AIMNet2ChargeResponseGradientWelfordTrajectoryResult.h"
#include "ApbsEfieldTimeSeriesTrajectoryResult.h"
#include "ApbsEfgTimeSeriesTrajectoryResult.h"
#include "MopacChargeWelfordTrajectoryResult.h"
#include "MopacBondOrderWelfordTrajectoryResult.h"
#include "MopacCoulombShieldingTimeSeriesTrajectoryResult.h"
#include "MopacMcConnellShieldingTimeSeriesTrajectoryResult.h"
#include "MopacVsFf14SbReconciliationTrajectoryResult.h"
#include "TripeptideBackboneShieldingTimeSeriesTrajectoryResult.h"
#include "TripeptideBackboneResidualVecTimeSeriesTrajectoryResult.h"
#include "TripeptideNeighborShieldingTimeSeriesTrajectoryResult.h"
#include "TripeptideNeighborResidualVecPrevTimeSeriesTrajectoryResult.h"
#include "TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult.h"
#include "TripeptideBackboneMethodTagTimeSeriesTrajectoryResult.h"
#include "LarsenHBondWaterTermTimeSeriesTrajectoryResult.h"
#include "LarsenHBondCountTimeSeriesTrajectoryResult.h"
#include "LarsenHBond1pHBShieldingTimeSeriesTrajectoryResult.h"
#include "LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult.h"
#include "LarsenHBond1pHaBShieldingTimeSeriesTrajectoryResult.h"
#include "LarsenHBond2pHaBShieldingTimeSeriesTrajectoryResult.h"
#include "BsAnomalousAtomMarkerTrajectoryResult.h"
#include "BsT0AutocorrelationTrajectoryResult.h"
#include "KernelDynamicsTrajectoryResult.h"
#include "KernelCoherenceTrajectoryResult.h"
#include "ReorientationalDynamicsTrajectoryResult.h"
#include "IRedOrderParameterTrajectoryResult.h"
#include "DihedralAutocorrelationTrajectoryResult.h"
#include "BondLengthStatsTrajectoryResult.h"
#include "PositionsTimeSeriesTrajectoryResult.h"

#include <type_traits>
#include <typeindex>

namespace nmr {

namespace {

// Register one or more TrajectoryResults onto a shape. The run
// "produces" these results; argument order is attach order is dispatch
// order, so it is load-bearing. Centralizes the unique_ptr<Derived> ->
// unique_ptr<TrajectoryResult> upcast and the deferral lambda the call
// sites would otherwise repeat per result.
template <class... TRs>
void Produces(RunConfiguration& c) {
    static_assert((std::is_base_of_v<TrajectoryResult, TRs> && ...),
                  "Produces<> takes TrajectoryResult subclasses");
    (c.AddTrajectoryResultFactory(
         [](const TrajectoryProtein& tp) -> std::unique_ptr<TrajectoryResult> {
             return TRs::Create(tp);
         }),
     ...);
}

}  // namespace


// ── PerFrameExtractionSet ────────────────────────────────────────
//
// Production canonical trajectory shape. Full classical stack + APBS +
// AIMNet2 every dispatched frame. MOPAC skipped (FullFat only);
// vacuum Coulomb skipped because APBS is the active electrostatic field.
//
// Stride is intentionally NOT set here: it defaults to 1 (every frame)
// and is the caller's single knob (CLI --stride → SetStride). A buried
// SetStride(2) here used to silently halve every production run while the
// CLI advertised "stride 1" at the emit layer — removed 2026-05-31.

RunConfiguration RunConfiguration::PerFrameExtractionSet() {
    RunConfiguration c;
    c.SetName("PerFrameExtractionSet");

    c.per_frame_opts_.skip_mopac   = true;   // sparse-frame only
    c.per_frame_opts_.skip_apbs    = false;
    c.per_frame_opts_.skip_coulomb = true;   // APBS supersedes
    c.per_frame_opts_.skip_dssp    = false;

    // Mandatory per frame; Phase 4 returns kConfigRequiresAimnet2
    // if the Session has no model loaded.
    c.SetRequiresAimnet2(true);

    // Required ConformationResult set for Phase 4 validation.
    c.RequireConformationResult(typeid(GeometryResult));
    c.RequireConformationResult(typeid(SpatialIndexResult));
    c.RequireConformationResult(typeid(EnrichmentResult));
    c.RequireConformationResult(typeid(DsspResult));
    c.RequireConformationResult(typeid(ChargeAssignmentResult));
    c.RequireConformationResult(typeid(ApbsFieldResult));
    c.RequireConformationResult(typeid(BiotSavartResult));
    c.RequireConformationResult(typeid(HaighMallionResult));
    c.RequireConformationResult(typeid(McConnellResult));
    c.RequireConformationResult(typeid(RingSusceptibilityResult));
    c.RequireConformationResult(typeid(PiQuadrupoleResult));
    c.RequireConformationResult(typeid(DispersionResult));
    c.RequireConformationResult(typeid(HBondResult));
    c.RequireConformationResult(typeid(SasaResult));
    c.RequireConformationResult(typeid(EeqResult));
    c.RequireConformationResult(typeid(AIMNet2Result));
    c.RequireConformationResult(typeid(AIMNet2ChargeResponseGradientResult));
    c.RequireConformationResult(typeid(WaterFieldResult));
    c.RequireConformationResult(typeid(HydrationShellResult));
    c.RequireConformationResult(typeid(HydrationGeometryResult));
    c.RequireConformationResult(typeid(GromacsFramePullResult));
    c.RequireConformationResult(typeid(GromacsEnergyResult));
    c.RequireConformationResult(typeid(BondedEnergyResult));
    // Note: TripeptideBackbone/NeighborShieldingResult are NOT required
    // here. They are conditionally attached by OperationRunner when the
    // tensorcs15 DSN is configured (project precondition). The
    // corresponding TimeSeries TRs capture whatever is in the source
    // ConformationAtom field — zero SphericalTensor if the calc didn't
    // run. Dep checks for always-present-when-DSN-configured deps would
    // be cruft that breaks fleet runs that legitimately lack DSN.

    // Attach order is dispatch order — the argument order below is
    // load-bearing. BsWelford comes first so downstream TRs that
    // cross-read its fields (BsAnomalousAtomMarker) see fresh values.

    // ── Welford accumulators (per-atom mean/variance over frames) ──
    Produces<BsWelfordTrajectoryResult,
             HmWelfordTrajectoryResult,
             McConnellWelfordTrajectoryResult,
             EeqWelfordTrajectoryResult,
             SasaWelfordTrajectoryResult,
             HBondCountWelfordTrajectoryResult>(c);

    // ── Per-frame classical shielding time-series ──
    Produces<BsShieldingTimeSeriesTrajectoryResult,
             HmShieldingTimeSeriesTrajectoryResult,
             McConnellShieldingTimeSeriesTrajectoryResult,
             PiQuadrupoleShieldingTimeSeriesTrajectoryResult,
             RingSusceptibilityShieldingTimeSeriesTrajectoryResult,
             DispersionShieldingTimeSeriesTrajectoryResult,
             HBondShieldingTimeSeriesTrajectoryResult,
             SasaTimeSeriesTrajectoryResult>(c);

    // ── AIMNet2 (charge, embedding, charge-response gradient) ──
    Produces<AIMNet2ChargeTimeSeriesTrajectoryResult,
             AIMNet2EmbeddingTimeSeriesTrajectoryResult,
             AIMNet2ChargeResponseGradientTimeSeriesTrajectoryResult,
             AIMNet2ChargeResponseGradientWelfordTrajectoryResult>(c);

    // ── APBS Poisson–Boltzmann (E-field + EFG) ──
    Produces<ApbsEfieldTimeSeriesTrajectoryResult,
             ApbsEfgTimeSeriesTrajectoryResult>(c);

    // The MOPAC family (MopacCoulomb / MopacMcConnell shielding, charge
    // and bond-order Welford, the MOPAC-vs-FF14SB reconciliation) is NOT
    // here: PerFrameExtractionSet skips MOPAC and vacuum Coulomb, so those
    // source calcs never attach per frame and the TRs would be dead-code.
    // They live in FullFatFrameExtraction, which enables both. Decision
    // 2026-05-21 per science + math adversarial review H3.

    // ── Tripeptide backbone + neighbor (ProCS15) shielding + residuals ──
    Produces<TripeptideBackboneShieldingTimeSeriesTrajectoryResult,
             TripeptideBackboneResidualVecTimeSeriesTrajectoryResult,
             TripeptideNeighborShieldingTimeSeriesTrajectoryResult,
             TripeptideNeighborResidualVecPrevTimeSeriesTrajectoryResult,
             TripeptideNeighborResidualVecNextTimeSeriesTrajectoryResult,
             TripeptideBackboneMethodTagTimeSeriesTrajectoryResult>(c);

    // ── Larsen H-bond family (water term, count, 1°/2° HB + HαB) ──
    Produces<LarsenHBondWaterTermTimeSeriesTrajectoryResult,
             LarsenHBondCountTimeSeriesTrajectoryResult,
             LarsenHBond1pHBShieldingTimeSeriesTrajectoryResult,
             LarsenHBond2pHBShieldingTimeSeriesTrajectoryResult,
             LarsenHBond1pHaBShieldingTimeSeriesTrajectoryResult,
             LarsenHBond2pHaBShieldingTimeSeriesTrajectoryResult>(c);

    // ── Biot–Savart-derived per-atom markers (read BsWelford, above) ──
    Produces<BsAnomalousAtomMarkerTrajectoryResult,
             BsT0AutocorrelationTrajectoryResult>(c);

    // ── Dynamics observables ──
    // The instrument: per-atom autocorrelation + power spectrum of every
    // shielding kernel (KernelDynamics) and the zero-lag kernel-kernel
    // coherence (KernelCoherence). The model-free layer: backbone
    // reorientational S^2 / tau_e / TCF (Reorientational), reference-free
    // iRED order parameters, and torsional circular ACF
    // (DihedralAutocorrelation). None cross-reads another, so order here is
    // free; KernelDynamics/Coherence depend on the classical kernel
    // ConformationResults required above.
    Produces<KernelDynamicsTrajectoryResult,
             KernelCoherenceTrajectoryResult,
             ReorientationalDynamicsTrajectoryResult,
             IRedOrderParameterTrajectoryResult,
             DihedralAutocorrelationTrajectoryResult>(c);

    // ── Geometry + energy bookkeeping ──
    Produces<BondLengthStatsTrajectoryResult,
             PositionsTimeSeriesTrajectoryResult,
             GromacsEnergyTimeSeriesTrajectoryResult,
             BondedEnergyTimeSeriesTrajectoryResult>(c);

    // ── Explicit solvent: water field + hydration geometry/shell ──
    Produces<WaterFieldTimeSeriesTrajectoryResult,
             WaterFieldWelfordTrajectoryResult,
             HydrationGeometryTimeSeriesTrajectoryResult,
             HydrationGeometryWelfordTrajectoryResult,
             HydrationShellTimeSeriesTrajectoryResult,
             HydrationShellWelfordTrajectoryResult>(c);

    // ── Backbone/sidechain conformation (dihedrals, DSSP8, ring pucker) ──
    Produces<DihedralTimeSeriesTrajectoryResult,
             DihedralBinTransitionTrajectoryResult,
             Dssp8TimeSeriesTrajectoryResult,
             Dssp8TransitionTrajectoryResult,
             RingPuckerTimeSeriesTrajectoryResult>(c);

    // ── Scalar couplings + ring spatial neighbourhood ──
    Produces<JCouplingTimeSeriesTrajectoryResult,
             RingNeighbourhoodTrajectoryStats>(c);

    // ── Frame selection (RMSD tracking/spikes, chi rotamer, DFT pose) ──
    // ChiRotamerSelection feeds DftPoseCoordinator (which reads its
    // SelectionBag at Finalize) and emits sidechain-rotamer DFT pose
    // candidates alongside RmsdSpikeSelection's whole-protein spikes.
    Produces<RmsdTrackingTrajectoryResult,
             RmsdSpikeSelectionTrajectoryResult,
             ChiRotamerSelectionTrajectoryResult,
             DftPoseCoordinatorTrajectoryResult>(c);

    return c;
}


// ── FullFatFrameExtraction ───────────────────────────────────────
//
// PerFrameExtractionSet with MOPAC and vacuum Coulomb enabled. MOPAC runs
// on every dispatched frame (the single --stride governs the cadence);
// there is no separate MOPAC stride.
//
// MOPAC-family ConformationResult sources attach conditionally inside
// OperationRunner when MOPAC runs and succeeds; the MOPAC-family
// TrajectoryResults use HasResult gates rather than hard Phase-4
// RequireConformationResult dependencies.

RunConfiguration RunConfiguration::FullFatFrameExtraction() {
    RunConfiguration c = PerFrameExtractionSet();
    c.SetName("FullFatFrameExtraction");

    // Both MOPAC AND vacuum Coulomb attach here — required for the
    // MOPAC-family TRs (TR5-TR9 of the 13-TR plan) and for TR9's
    // cross-source MopacCoulomb-vs-Coulomb reconciliation. APBS
    // stays on (inherited from PerFrameExtractionSet) for the
    // hybrid APBS-Mopac calibration probe.
    c.per_frame_opts_.skip_mopac   = false;
    c.per_frame_opts_.skip_coulomb = false;

    // ── MOPAC family ── only meaningful once skip_mopac and skip_coulomb
    // are false (above), so these live here rather than in
    // PerFrameExtractionSet. MopacVsFf14SbReconciliation is the cross-
    // source MopacCoulomb-vs-FF14SB-Coulomb probe (TR9).
    Produces<MopacChargeWelfordTrajectoryResult,
             MopacBondOrderWelfordTrajectoryResult,
             MopacCoulombShieldingTimeSeriesTrajectoryResult,
             MopacMcConnellShieldingTimeSeriesTrajectoryResult,
             MopacVsFf14SbReconciliationTrajectoryResult>(c);

    return c;
}

}  // namespace nmr

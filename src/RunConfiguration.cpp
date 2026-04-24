#include "RunConfiguration.h"
#include "TrajectoryProtein.h"
#include "TrajectoryResult.h"

// Per-frame ConformationResult types that TrajectoryResults in this
// session (BsWelford) depend on. Include full types so type_index
// values are valid at factory call time.
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
#include "BondedEnergyResult.h"

// Concrete TrajectoryResults populating the canonical configurations.
#include "BsWelfordTrajectoryResult.h"
#include "BsShieldingTimeSeriesTrajectoryResult.h"
#include "BsAnomalousAtomMarkerTrajectoryResult.h"
#include "BsT0AutocorrelationTrajectoryResult.h"
#include "BondLengthStatsTrajectoryResult.h"
#include "PositionsTimeSeriesTrajectoryResult.h"

#include <typeindex>

namespace nmr {


// ── ScanForDftPointSet ───────────────────────────────────────────
//
// Cheap per-frame (no MOPAC, APBS, Coulomb, AIMNet2). Intended for
// rotamer / RMSD / bin-crossing detection to choose frames for
// FullFatFrameExtraction. The ConformationResults listed below are
// what OperationRunner attaches under these opts; TRs declaring any
// of them as Dependencies() must find them here (Phase 4 validation).

RunConfiguration RunConfiguration::ScanForDftPointSet() {
    RunConfiguration c;
    c.SetName("ScanForDftPointSet");

    // Per-frame options: cheap set.
    c.per_frame_opts_.skip_mopac   = true;
    c.per_frame_opts_.skip_apbs    = true;
    c.per_frame_opts_.skip_coulomb = true;
    c.per_frame_opts_.skip_dssp    = false;   // dihedral bins need phi/psi/chi

    // ConformationResults that must run each frame for the attached
    // TrajectoryResults to have valid inputs.
    c.RequireConformationResult(typeid(GeometryResult));
    c.RequireConformationResult(typeid(SpatialIndexResult));
    c.RequireConformationResult(typeid(EnrichmentResult));
    c.RequireConformationResult(typeid(DsspResult));
    c.RequireConformationResult(typeid(BiotSavartResult));
    c.RequireConformationResult(typeid(SasaResult));

    // Attached TrajectoryResults: BsWelford. The scan-selection TRs
    // (ChiRotamerSelection et al.) are a pending-decision item —
    // see spec/pending_decisions_20260423.md item 2.
    c.AddTrajectoryResultFactory(
        [](const TrajectoryProtein& tp) -> std::unique_ptr<TrajectoryResult> {
            return BsWelfordTrajectoryResult::Create(tp);
        });

    return c;
}


// ── PerFrameExtractionSet ────────────────────────────────────────
//
// Production canonical for the 685-protein fleet. Full classical
// stack + APBS + AIMNet2 every frame. MOPAC skipped (FullFat only);
// vacuum Coulomb skipped (APBS supersedes at N > 1000 atoms).
// Stride 2: 25 ns × 1250 frames → 625 sampled.

RunConfiguration RunConfiguration::PerFrameExtractionSet() {
    RunConfiguration c;
    c.SetName("PerFrameExtractionSet");

    c.per_frame_opts_.skip_mopac   = true;   // sparse-frame only
    c.per_frame_opts_.skip_apbs    = false;
    c.per_frame_opts_.skip_coulomb = true;   // APBS supersedes
    c.per_frame_opts_.skip_dssp    = false;

    // Production stride: 25 ns × 1250 frames × stride 2 → 625 sampled.
    c.SetStride(2);

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
    c.RequireConformationResult(typeid(WaterFieldResult));
    c.RequireConformationResult(typeid(HydrationShellResult));
    c.RequireConformationResult(typeid(HydrationGeometryResult));
    c.RequireConformationResult(typeid(GromacsEnergyResult));
    c.RequireConformationResult(typeid(BondedEnergyResult));

    // Attach order is dispatch order. BsWelford runs first so
    // downstream TRs that cross-read its fields
    // (BsAnomalousAtomMarker) see fresh values.
    c.AddTrajectoryResultFactory(
        [](const TrajectoryProtein& tp) -> std::unique_ptr<TrajectoryResult> {
            return BsWelfordTrajectoryResult::Create(tp);
        });
    c.AddTrajectoryResultFactory(
        [](const TrajectoryProtein& tp) -> std::unique_ptr<TrajectoryResult> {
            return BsShieldingTimeSeriesTrajectoryResult::Create(tp);
        });
    c.AddTrajectoryResultFactory(
        [](const TrajectoryProtein& tp) -> std::unique_ptr<TrajectoryResult> {
            return BsAnomalousAtomMarkerTrajectoryResult::Create(tp);
        });
    c.AddTrajectoryResultFactory(
        [](const TrajectoryProtein& tp) -> std::unique_ptr<TrajectoryResult> {
            return BsT0AutocorrelationTrajectoryResult::Create(tp);
        });
    c.AddTrajectoryResultFactory(
        [](const TrajectoryProtein& tp) -> std::unique_ptr<TrajectoryResult> {
            return BondLengthStatsTrajectoryResult::Create(tp);
        });
    c.AddTrajectoryResultFactory(
        [](const TrajectoryProtein& tp) -> std::unique_ptr<TrajectoryResult> {
            return PositionsTimeSeriesTrajectoryResult::Create(tp);
        });

    return c;
}


// ── FullFatFrameExtraction ───────────────────────────────────────
//
// PerFrameExtractionSet with MOPAC enabled. Meant for a selected-frame
// subset (DFT pose set, harvester checkpoints) — every-frame MOPAC
// on a 25 ns trajectory is ~15 h/protein. Selected-frame mechanism
// is a pending-decision item (spec/pending_decisions_20260423.md).
//
// MOPAC-family ConformationResult deps (MopacResult,
// MopacCoulombResult, MopacMcConnellResult) are a pending-decision
// item — see spec/pending_decisions_20260423.md item 3.

RunConfiguration RunConfiguration::FullFatFrameExtraction() {
    RunConfiguration c = PerFrameExtractionSet();
    c.SetName("FullFatFrameExtraction");

    c.per_frame_opts_.skip_mopac = false;

    return c;
}

}  // namespace nmr

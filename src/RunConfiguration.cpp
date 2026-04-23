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
// Cheap per-frame: no MOPAC, no vacuum Coulomb, no APBS, no AIMNet2.
// The intended use is rotamer / RMSD / bin-crossing detection to
// choose frames for FullFatFrameExtraction. BiotSavart runs because
// ring proximity matters for selection; BsWelford aggregates it.
//
// Per-frame ConformationResults that actually run are determined by
// OperationRunner based on RunOptions — the required set below is a
// declaration for dependency validation, not a drive list.

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

    // Attached TrajectoryResults. This session: only BsWelford.
    // Follow-up sessions populate the rest per Appendix F:
    //   DihedralBinTransitionTrajectoryResult
    //   RmsdTrackingTrajectoryResult
    //   SasaWelfordTrajectoryResult
    //   ChiRotamerSelectionTrajectoryResult (mixin)
    //   DftPoseCoordinatorTrajectoryResult
    c.AddTrajectoryResultFactory(
        [](const TrajectoryProtein& tp) -> std::unique_ptr<TrajectoryResult> {
            return BsWelfordTrajectoryResult::Create(tp);
        });

    return c;
}


// ── PerFrameExtractionSet ────────────────────────────────────────
//
// The production canonical configuration for the 685-protein fleet
// extraction. Full classical stack per frame, MOPAC skipped (sparse-
// frame only under FullFatFrameExtraction), Coulomb skipped (APBS
// supersedes — faster and solvated for N > 1000 atoms).
//
// Stride 2 is the production sample: 25 ns trajectories at 1250
// frames → 625 sampled. The fleet runs ~685 of these; standardising
// here means the shape downstream consumers (calibration, analysis)
// see is uniform across proteins.
//
// Per-frame ConformationResults: every type OperationRunner actually
// attaches under these opts, given that trajectory runs always
// provide opts.solvent (from XTC), opts.frame_energy (from preloaded
// EDR), opts.bonded_params (from TPR), and the session's AIMNet2
// model (mandatory per RequiresAimnet2 below). Declaring the full
// set here lets future TrajectoryResults legitimately depend on any
// of them via Dependencies(); Phase 4's validation then catches a
// configuration drift rather than letting a TR fail silently on
// stale zeros.

RunConfiguration RunConfiguration::PerFrameExtractionSet() {
    RunConfiguration c;
    c.SetName("PerFrameExtractionSet");

    c.per_frame_opts_.skip_mopac   = true;   // sparse-frame only
    c.per_frame_opts_.skip_apbs    = false;
    c.per_frame_opts_.skip_coulomb = true;   // APBS supersedes
    c.per_frame_opts_.skip_dssp    = false;

    // Production stride: 25 ns × 1250 frames × stride 2 → 625 sampled.
    c.SetStride(2);

    // AIMNet2 is MANDATORY per frame — neural-network Hirshfeld
    // charges + 256-dim embedding. Trajectory::Run Phase 4 returns
    // kConfigRequiresAimnet2 if the caller's Session has no model
    // loaded (prevents silent-switch-off when a model path is
    // missing).
    c.SetRequiresAimnet2(true);

    // ── Canonical per-frame ConformationResult set ──
    // Every type OperationRunner attaches under the opts above. The
    // order here is not load-bearing (dispatch uses attach order
    // within OperationRunner itself); declared as a set for
    // validation.
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

    // ── Canonical per-trajectory TrajectoryResults ──
    // The two worked examples currently exercise the two canonical
    // shapes: BsWelfordTrajectoryResult (always-valid-mid-stream
    // rollup on a ConformationAtom field) and
    // PositionsTimeSeriesTrajectoryResult (Finalize-only
    // DenseBuffer<Vec3> time series with no ConformationResult dep).
    // The rest of Appendix F lands in follow-up sessions — one
    // AddTrajectoryResultFactory line per class.
    // Attach order is dispatch order. BsWelford runs first on each
    // frame so its accumulator state is fresh when downstream
    // Results (BsAnomalousAtomMarker) read from tp.AtomAt(i) during
    // the same frame.
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
// PerFrameExtractionSet plus MOPAC family on the frames that get run.
// Framing note (WIP §6): this configuration is meant to run only on a
// selected-frame subset (DFT pose set or μs harvester checkpoints) —
// every-frame MOPAC on a 25 ns trajectory is ~15 hours per protein.
// Selection happens at the caller layer (frame-filter XTC from a
// ScanForDftPointSet run, or RunContext::SetSelectedFrameIndices in
// a future iteration).

RunConfiguration RunConfiguration::FullFatFrameExtraction() {
    RunConfiguration c = PerFrameExtractionSet();
    c.SetName("FullFatFrameExtraction");

    c.per_frame_opts_.skip_mopac = false;

    // MOPAC-family ConformationResult dependencies would be added
    // here once their types are referenced (MopacResult,
    // MopacCoulombResult, MopacMcConnellResult). Deferred until the
    // MOPAC-family TrajectoryResults land in a follow-up session;
    // including the typeids without attached TrajectoryResults that
    // need them is harmless (the set just lists more than required).

    // No additional TrajectoryResults in this session.

    return c;
}

}  // namespace nmr

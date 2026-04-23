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

// First concrete TrajectoryResult.
#include "BsWelfordTrajectoryResult.h"

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
// Full classical stack per frame (analysis mode today). Produces the
// exhaustive per-frame time-series H5 and the per-atom Welford
// rollups consumed by calibration.
//
// Per-frame ConformationResults (required): the full classical set
// minus MOPAC. MOPAC is optional here — selected frames only, and
// only under FullFatFrameExtraction.

RunConfiguration RunConfiguration::PerFrameExtractionSet() {
    RunConfiguration c;
    c.SetName("PerFrameExtractionSet");

    c.per_frame_opts_.skip_mopac   = true;   // sparse-frame only
    c.per_frame_opts_.skip_apbs    = false;
    c.per_frame_opts_.skip_coulomb = true;   // APBS supersedes
    c.per_frame_opts_.skip_dssp    = false;
    // AIMNet2 is MANDATORY per frame — the neural-network charge +
    // 256-dim embedding calculator that drives AIMNet2-derived
    // features. Caller must supply the model via
    // RunContext::SetAimnet2Model; Trajectory::Run Phase 2 throws
    // otherwise.
    c.SetRequiresAimnet2(true);

    c.RequireConformationResult(typeid(GeometryResult));
    c.RequireConformationResult(typeid(SpatialIndexResult));
    c.RequireConformationResult(typeid(EnrichmentResult));
    c.RequireConformationResult(typeid(DsspResult));
    c.RequireConformationResult(typeid(ChargeAssignmentResult));
    c.RequireConformationResult(typeid(ApbsFieldResult));
    c.RequireConformationResult(typeid(BiotSavartResult));
    c.RequireConformationResult(typeid(HaighMallionResult));
    c.RequireConformationResult(typeid(McConnellResult));
    c.RequireConformationResult(typeid(SasaResult));

    // This session: only BsWelford.
    // Follow-up sessions populate per Appendix F: positions time
    // series, all *TimeSeriesTrajectoryResult family, Welford
    // rollups for MC/HM/RS/PQ/Disp/AIMNet2/SASA/Water/EEQ/DSSP/
    // bonded energy, etc.
    c.AddTrajectoryResultFactory(
        [](const TrajectoryProtein& tp) -> std::unique_ptr<TrajectoryResult> {
            return BsWelfordTrajectoryResult::Create(tp);
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

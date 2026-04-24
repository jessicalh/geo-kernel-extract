// RunConfiguration.cpp
//
// Per WIP_OBJECT_MODEL.md §6. The three named configurations. Each
// factory reads as intent; "what does this run do?" is answered by
// reading one function (spec §6).
//
// The traj_factories_ list here carries ONLY the one worked-example
// TrajectoryResult this sandbox implements (BsWelfordTrajectoryResult).
// Appendix F's full catalogue of ~30 TrajectoryResult classes names
// what the real factories carry; that is the implementation checklist,
// not the framework — see README.

#include "RunConfiguration.h"

#include "BsWelfordTrajectoryResult.h"

RunConfiguration RunConfiguration::ScanForDftPointSet() {
    RunConfiguration c;
    c.name_ = "ScanForDftPointSet";

    // Per-frame: very cheap. No APBS, no MOPAC, no Coulomb, no AIMNet2.
    c.per_frame_opts_.skip_dssp = false;   // need DSSP for Ramachandran bins
    c.per_frame_opts_.skip_mopac = true;
    c.per_frame_opts_.skip_apbs = true;
    c.per_frame_opts_.skip_coulomb = true;
    c.per_frame_opts_.aimnet2_model = nullptr;

    // The ConformationResult types that the scan-mode per-frame pipeline
    // runs, per §6 body. Real required_conf_result_types_ contains
    // typeid(GeometryResult), typeid(SpatialIndexResult), typeid(DsspResult),
    // typeid(SasaResult). The sandbox forward-declares none of these
    // so the set stays empty — the shape is what matters for the framework.

    // Real factories per §6: DihedralBinTransitionTrajectoryResult,
    // RmsdTrackingTrajectoryResult, DftPoseSelectionTrajectoryResult,
    // plus the selection-emitting ones per Appendix F.
    //
    // The sandbox keeps the list empty for ScanForDftPointSet; it does
    // not implement those Results.

    return c;
}

RunConfiguration RunConfiguration::PerFrameExtractionSet() {
    RunConfiguration c;
    c.name_ = "PerFrameExtractionSet";

    // Per-frame: the full classical stack, no MOPAC, per §6.
    c.per_frame_opts_.skip_dssp = false;
    c.per_frame_opts_.skip_mopac = true;
    c.per_frame_opts_.skip_apbs = false;
    c.per_frame_opts_.skip_coulomb = false;

    // Real required_conf_result_types_ includes typeid(BiotSavartResult)
    // among ~20 others per §6. BsWelfordTrajectoryResult's Dependencies()
    // returns typeid(BiotSavartResult); Trajectory::Run Phase 2 consults
    // this set. The sandbox cannot register typeid(BiotSavartResult)
    // here because Stubs.h forward-declares BiotSavartResult without
    // a body, and typeid on a forward-declared class requires the
    // complete type. This is noted in README.

    c.traj_factories_ = {
        [](const TrajectoryProtein& tp) {
            return BsWelfordTrajectoryResult::Create(tp);
        },
        // Appendix F's remaining entries (PositionsTimeSeries,
        // HmWelford, McConnellWelford, Coulomb*, APBS*, Water*, ...)
        // go here in the real factory. Not implemented in sandbox.
    };

    return c;
}

RunConfiguration RunConfiguration::FullFatFrameExtraction() {
    RunConfiguration c = PerFrameExtractionSet();  // inherit, per §6
    c.name_ = "FullFatFrameExtraction";

    // Add MOPAC.
    c.per_frame_opts_.skip_mopac = false;

    // Real code also inserts typeid(MopacResult), typeid(MopacCoulombResult),
    // typeid(MopacMcConnellResult) into required_conf_result_types_,
    // and pushes MopacCoulombTimeSeries, MopacMcConnellTimeSeries,
    // MopacVsFf14SbReconciliation onto traj_factories_. Not implemented
    // in sandbox.

    return c;
}

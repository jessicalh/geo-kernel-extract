# Next Session Prompt (SUPERSEDED — see next_session_2026-04-11.md)

**This prompt is from before the 2026-04-10 session which settled the
calibration.  The MLP questions below are answered: per-element ridge
beats the MLP (0.818 vs 0.61).  Read next_session_2026-04-11.md instead.**

Read spec/INDEX.md first (it tells you what to read next). Then read
learn/EXPERIMENTS.md — the full experiment log including the 2026-04-09
session. Then read learn/docs/secondary_analysis_2026-04-09.md and
learn/docs/workbench_physics_analysis_2026-04-09.md.

## Where we are

The calibration pipeline works: ridge R²=0.35, MLP R²=0.61 on T2
across 723 proteins. We built secondary analysis tools (5 Python +
R graphics), an R analytical workbench, and a hidden layer probe.

## The key finding from the last session

The MLP hidden space collapses to 2-3 dimensions out of 64. The
model learns a near-constant relative weight vector across all atoms.
The kernel self-gating (magnitude / (magnitude + threshold)) provides
all the spatial differentiation — 9 effective dimensions in the gated
weights, all from the kernels' own magnitude variation.

This means the 0.35→0.61 improvement is mostly gating, not learned
scalar-dependent kernel weighting. The MLP is not using the 261
scalars to modulate weights per-atom. It found one good set of
relative importances and lets the gating handle the rest.

## The open question

Is the MLP FAILING to learn scalar-dependent weights, or is one
weight vector plus gating genuinely optimal?

Test this by:
1. **Tiny hidden layer (hidden_mix=4 or 8).** If R² barely changes
   from 0.61, the MLP's learned content is truly low-dimensional
   and the hidden=64 was always wasted capacity.
2. **Disable gating entirely.** Set gate = 1.0 for all kernels,
   force the MLP to do all the spatial differentiation through
   scalar-dependent weights. If R² drops to ~0.35 (ridge level),
   gating IS the mechanism. If it stays near 0.61, the MLP can
   learn the same thing without gating.
3. **Per-protein weight analysis.** Extract the gated weight vector
   per protein (average across atoms in each protein). If the
   weight vectors cluster into a few types (PHE-dominated proteins
   vs TRP-dominated vs HIE-dominated), the MLP learned a protein-
   type classifier, not per-atom context.

## Also pending

- The AzimuthalExtraction is running (~110 of 723 proteins done as
  of 2026-04-09 evening). When complete, run the full 723-protein
  workbench with the 91+new-scalars config to get definitive numbers.
- Stratified forward selection per ring type (PHE-only vs HIE-only
  atoms) — tells you whether the kernel ordering changes by ring type.
- The pose selection idea is written up in
  /shared/2026Thesis/nmr-training/POSE_SELECTION_FROM_CALIBRATION.md
- Viewer annotation exporter (learn/src/secondary/viewer_export.py)
  ready but not yet tested against the actual viewer.

## Files to know

- learn/src/secondary/ — 7 Python tools (divergence, ring_strata,
  kernel_structure, ablation, scalar_ablation, viewer_export,
  export_for_r, hidden_probe)
- learn/R/ — 8 R scripts (common, divergence, ring_strata,
  kernel_structure, ablation, scalar_ablation, workbench,
  hypothesis_e)
- learn/runs/ — 5 experiment runs from 2026-04-09
- spec/AZIMUTHAL_ANGLE.md — C++ spec for cos(phi), sin(phi)
- learn/docs/viewer_overlays.md — VTK overlay ideas

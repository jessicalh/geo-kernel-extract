# Codex — prior-art research: trajectory display stabilisation / structural superposition

## (The standard qt6-cpp + VTK-docs preamble is prepended above. For THIS task you ALSO have explicit BROAD web-search license — see below. You will NOT change any code: this is research-and-write only.)

## You have explicit broad web-search license for this task
Search the open web freely — official docs, source repos, papers, mailing lists, forum
threads — for how established molecular-structure / MD tools stabilise a trajectory for
display and compute structural superposition. **Cite every claim with a URL**; quote the
tool's own docs or source where you can. "I don't know / the docs don't say" beats a
confident guess (project ethos: no handwaving about external systems).

## Why we're asking
The h5-reader displays one MD trajectory and removes global tumbling + translation so the
viewer sees internal motion, not box drift. Current approach: fit each frame (Kabsch) to
**frame 0**, then temporally smooth the (R,T) sequence over a window. We just diagnosed a
bug — smoothing R and T *independently*, with the rotation applied **about the origin**
while the molecule sits ~93 Å away, walks the whole molecule across the screen ("duck
walk"; measured ~40 Å excursion by the trajectory end). We are about to rebuild it as:
**(A)** pin the centroid — rotate about the molecule's own centroid and place it on the
reference centroid; **(B)** fit to an iterative converged-**MEAN** structure over the whole
trajectory instead of frame 0. Before building, the lead wants to know what the mature
tools do — "PyMOL must do something." Use them as reference behaviour to check our plan.

## Research questions — answer each, with citations
1. **Reference structure.** Do tools fit each frame to (a) the first frame, (b) a chosen
   reference, or (c) an iteratively-computed AVERAGE/mean structure? Document the
   iterative average-structure superposition specifically:
   - GROMACS: `gmx rmsf` average structure; `gmx trjconv -fit rot+trans`; how `-fit`
     picks/uses its reference; any iterate-to-average option.
   - MDAnalysis: `align.AverageStructure` + `align.AlignTraj` (the documented "align to the
     mean, iterate" recipe) — iteration count, convergence criterion.
   - Theseus (Theobald & Wuttke): maximum-likelihood superposition to an average — what it
     adds over plain least-squares.
2. **Translation / centroid handling.** Is translation removed by centre-of-geometry or
   centre-of-mass? Is the standard to CENTRE both structures (subtract centroid) BEFORE
   computing AND applying the rotation — i.e. rotation is about the centroid, never the
   origin? Confirm across Kabsch / Kabsch–Umeyama, GROMACS, MDAnalysis `_fit_to`. This is
   the crux of our duck-walk fix — verify the established form.
3. **Fit selection & weighting (the "gyroscope").** What atom set is fitted (all atoms,
   backbone, Cα, heavy atoms, a rigid core)? Critically: does any tool DOWN-WEIGHT or
   iteratively reject mobile atoms so flexible loops/sidechains don't drag the rotation
   (our all-atom precession)? Look at Theseus ML weighting, GROMACS fit groups, MDAnalysis
   weighted/`select` fits, iterative outlier rejection. Is "fit on a stable core" the
   standard answer to apparent precession?
4. **Temporal smoothing of the alignment.** Do any trajectory VIEWERS temporally smooth the
   alignment transform over a time window (like our windowed quaternion mean)? Or is the
   universal approach a purely per-frame fit to a stable (mean) reference with NO temporal
   coupling? Does PyMOL's `smooth` or VMD's trajectory smoothing average COORDINATES — and
   if so, does that distort internal motion? We want to know whether windowed
   transform-smoothing is a thing mature tools do at all, or whether a good reference
   removes the need for it.
5. **PyMOL.** `intra_fit` / `intra_rms_cur` (reference across states), `fit` / `align` /
   `super`, and the `smooth` command (what it actually does to coordinates); any movie /
   trajectory stabilisation.
6. **VMD.** RMSD Trajectory Tool / RMSD Visualizer / `measure fit` — fit to a reference
   frame vs an average; how it centres. (ChimeraX `mseries`/`align` too, if relevant.)

## Deliverable
Write `notes/STABILISATION_PRIORART_2026-06-05.md`: a concise per-tool summary (PyMOL, VMD,
GROMACS, MDAnalysis, Theseus, ChimeraX if relevant), each claim cited with a URL. End with a
RECOMMENDATION section answering: (1) frame-0 vs iterative-mean reference — what's standard;
(2) centroid-pinning (rotate about the centroid) — confirm it's the universal form; (3) fit
weighting / stable-core — is that the standard cure for all-atom precession; (4) temporal
smoothing — does anyone do it, or does a stable reference make it unnecessary. Flag anything
that contradicts our planned A+B fix. Do NOT change code. Do NOT run git.

# Stabilisation Prior Art - 2026-06-05

Research-only note. No code changes.

## Local h5-reader context verified

- The display wrapper is documented as caching one rigid transform per frame, computed by raw Kabsch and then temporally smoothed before `atomPosition()` applies `R * raw + T`: `src/model/TransformedConformation.h:20`, `src/model/TransformedConformation.h:22`, `src/model/TransformedConformation.h:80`.
- Startup/default behaviour is documented as `FitSubset` over the typed backbone subset, reference frame 0: `src/model/TransformedConformation.h:36`, `src/model/TransformedConformation.h:37`.
- The implementation caches reference positions from `referenceFrame_`, computes a raw transform sequence, then calls `smoothTransformSequence(raw)`: `src/model/TransformedConformation.cpp:140`, `src/model/TransformedConformation.cpp:157`, `src/model/TransformedConformation.cpp:265`.
- The current smoothing averages rotation quaternions and translations independently (`qSum += ...`, `tSum += raw[k].T`, then `out.T = tSum / count`): `src/model/TransformedConformation.cpp:298`, `src/model/TransformedConformation.cpp:308`, `src/model/TransformedConformation.cpp:321`.
- The underlying Kabsch math itself is centroid-form, not origin-form: it computes `cc` and `cr`, subtracts them before SVD, and sets `T = cr - R * cc`: `src/app/FitTargetMath.h:189`, `src/app/FitTargetMath.h:201`, `src/app/FitTargetMath.h:202`, `src/app/FitTargetMath.h:227`, `src/app/FitTargetMath.h:228`.

## PyMOL

- `intra_fit` is a state-to-state stabilisation command: PyMOL says it "fits all states of an object to an atom selection in the specified state" and takes `selection` plus target `state`, so the reference is a chosen state, not an automatically computed average. Source: https://pymol.org/pymol-command-ref.html lines 3203-3223.
- PyMOL's own trajectory movie tutorial recommends `intra_fit SampleTrajectory`; it says this "superposes all states on the first state and stops it from floating around." Source: https://www.pymol.org/tutorials/moviemaking/ lines 153-158.
- `intra_rms_cur` is not stabilisation: PyMOL says it calculates RMS values across states relative to the indicated state "without performing any fitting." Source: https://pymol.org/pymol-command-ref.html lines 3264-3286.
- `align` and `super` are not mean-trajectory stabilisers; they are pairwise structural superposition commands with iterative outlier rejection. PyMOL says `align` performs sequence alignment, structural superposition, and refinement cycles to reject outliers; it exposes `cutoff` and `cycles`, and says setting `cycles=0` disables outlier rejection. Source: https://pymol.org/pymol-command-ref.html lines 368-416. PyMOL says `super` performs residue-based pairwise alignment, structural superposition, and refinement cycles to reject outliers. Source: https://pymol.org/pymol-command-ref.html lines 6891-6917.
- PyMOL's `smooth` is coordinate smoothing, not transform smoothing. The command reference says it "performs a window average of coordinate states" and is used to suppress high-frequency vibrations in MD trajectories. Source: https://pymol.org/pymol-command-ref.html lines 6602-6626. The movie tutorial applies `smooth` after `intra_fit` as a separate coordinate filter. Source: https://www.pymol.org/tutorials/moviemaking/ lines 153-162.

## VMD

- VMD's `measure fit` returns a 4x4 matrix that, when applied to `selection1`, minimizes weighted RMSD to `selection2`; the `weight` argument can be `none`, an atom-selection keyword such as `mass`, or a per-atom list. Source: https://www.ks.uiuc.edu/Research/vmd/current/ug/node138.html lines 3-5 and 18.
- VMD's `measure center` returns the geometric center by default, with a `weight` option; using `mass` as the weight gives center-of-mass style behaviour. Source: https://www.ks.uiuc.edu/Research/vmd/current/ug/node138.html lines 3-8.
- VMD's RMSD Trajectory Tool documents explicit fit selections: `Trace` adds `name CA`, `Backbone` adds `name C CA N` or `name C CA N O`, and `noh` excludes hydrogens. Source: https://www.ks.uiuc.edu/Research/vmd/plugins/rmsdtt/ lines 35-40.
- RMSDTT supports weighted RMSD/alignment with fields including `mass`, `charge`, `beta`, and `occupancy`. Source: https://www.ks.uiuc.edu/Research/vmd/plugins/rmsdtt/ lines 49-52.
- RMSDTT aligns frames to a selected reference molecule/frame, but its average reference is for RMSD calculation rather than alignment: it says `ALIGN` uses the reference molecule, and "ALIGN is not available if ... Average is used as reference." Source: https://www.ks.uiuc.edu/Research/vmd/plugins/rmsdtt/ lines 53-60.
- I found VMD documentation for per-frame fitting, average positions (`measure avpos`), and frame playback, but not for temporal smoothing of the rigid alignment transform. Sources checked: `measure avpos` and `measure fit` at https://www.ks.uiuc.edu/Research/vmd/current/ug/node138.html lines 6 and 18; RMSDTT trajectory options at https://www.ks.uiuc.edu/Research/vmd/plugins/rmsdtt/ lines 61-70.

## GROMACS

- `gmx trjconv -fit rot+trans` fits each frame to a reference structure from the structure file; the current docs list `-fit` choices as `none`, `rot+trans`, `rotxy+transxy`, `translation`, `transxy`, and `progressive`. Source: https://manual.gromacs.org/current/onlinehelp/gmx-trjconv.html lines 485-486 and 615-617.
- GROMACS has a progressive fit option, but the docs describe it as fitting the first timeframe to the structure-file reference and each subsequent timeframe to the previously fitted structure, not as iterating to a converged mean structure. Source: https://manual.gromacs.org/current/onlinehelp/gmx-trjconv.html lines 485-486.
- `gmx rmsf` computes RMSF after optional fitting to a reference frame supplied with `-s`; it can write average coordinates with `-ox`, but its documented `-[no]fit` option is least-squares superposition before RMSF, not an iterate-to-average alignment pipeline. Source: https://manual.gromacs.org/current/onlinehelp/gmx-rmsf.html lines 497-499 and 580-582.
- `trjconv -center` uses the geometric center of a user-selected group for box centering, separate from least-squares fitting. Source: https://manual.gromacs.org/current/onlinehelp/gmx-trjconv.html lines 547-548.
- For movie smoothing, GROMACS points to `gmx filter` low-pass filtering to reduce aliasing of high-frequency motions; this is documented as filtering the trajectory/movie frames, not as smoothing the alignment transform. Source: https://manual.gromacs.org/current/onlinehelp/gmx-trjconv.html line 550.

## MDAnalysis

- `AlignTraj` fits a whole trajectory to a chosen reference `Universe`/frame using a selection; `select` defaults to `all`, and `weights` can be `mass`, `None`, or a per-atom array. Source: https://docs.mdanalysis.org/stable/documentation_pages/analysis/align.html lines 240-258.
- `AverageStructure` RMS-aligns a trajectory to a reference and calculates average coordinates; it exposes `ref_frame`, `select`, and `weights`. Source: https://docs.mdanalysis.org/stable/documentation_pages/analysis/align.html lines 331-369.
- `iterative_average` is the clearest mature-tool support for the planned mean-reference path: it iteratively computes an optimal reference that is itself the average structure after RMSD alignment, with defaults `niter=100` and `eps=1e-6`; convergence is when the RMSD distance between reference and average is below `eps`. Source: https://docs.mdanalysis.org/stable/documentation_pages/analysis/align.html lines 513-529.
- MDAnalysis explicitly documents the centroid-pinned transform form in `_fit_to`: input coordinates are already shifted to a center at the origin, and all moved atoms are transformed as `X' = R (X - Xbar) + Xbar_ref`. Source: https://docs.mdanalysis.org/stable/documentation_pages/analysis/align.html lines 564-594.
- The source follows the same form: `AlignTraj` stores `ref_atoms.positions - ref_com`, each frame uses `mobile_atoms.positions - mobile_com`, then passes both centers into `_fit_to`; `AverageStructure` accumulates positions after `_fit_to`. Source: https://docs.mdanalysis.org/stable/_modules/MDAnalysis/analysis/align.html lines 837-857 and 1063-1087.

## Theseus

- Theseus is explicitly an average/simultaneous-superposition tool, not a frame-0 trajectory viewer. The Theobald Lab page says it "simultaneously superimposes multiple macromolecular structures" with maximum likelihood rather than conventional least-squares. Source: https://www.theobaldlab.org/theseus lines 60-76.
- What Theseus adds over plain least-squares is directly relevant to flexible-loop precession: the lab page says ML "downweights variable regions of the superposition and corrects for correlations among atoms." Source: https://www.theobaldlab.org/theseus lines 73-77.
- Theseus also handles homologous proteins with gaps/missing data by using an ML algorithm that includes all data instead of discarding gap-aligned residues. Source: https://www.theobaldlab.org/theseus lines 77-78.

## ChimeraX

- ChimeraX `align` performs least-squares fitting of `matchatoms` onto `refatoms`, and for trajectories it has `each coordset` semantics: "match each coordinate set in matchatoms ... separately." Source: https://www.rbvi.ucsf.edu/chimerax/docs/user/commands/align.html.
- ChimeraX `matchmaker` fits one atom per aligned residue, CA for amino acids and C4' or P for nucleic acids; it points users to `align` for a different atom set. Source: https://www.rbvi.ucsf.edu/chimerax/docs/user/commands/matchmaker.html.
- ChimeraX Matchmaker has iterative fit behaviour that excludes bad regions: the tool docs say incorrect alignment portions "tend to be excluded from the fit during iteration." Source: https://www.cgl.ucsf.edu/chimerax/docs/user/tools/matchmaker.html.
- ChimeraX trajectory/model playback commands (`coordset`, `mseries`) play coordinate sets or model series; the docs found describe playback, per-coordinate-set alignment/RMSD, and morph interpolation, not windowed smoothing of the rigid alignment transform. Sources: https://rbvi.ucsf.edu/chimerax/docs/user/commands/coordset.html, https://www.cgl.ucsf.edu/chimerax/docs/user/commands/mseries.html, https://www.rbvi.ucsf.edu/chimerax/docs/user/commands/morph.html.

## Kabsch / Umeyama centroid form

- The standard rigid/similarity least-squares transform is centroid-pinned. Eigen's `umeyama()` source computes source and destination means, demeans both point sets, computes SVD on the demeaned covariance, and sets translation to `dst_mean - R * src_mean` when scaling is off. Source: https://eigen.googlesource.com/mirror/+/refs/heads/master/Eigen/src/Geometry/Umeyama.h lines 103-141.
- The Kabsch-Umeyama literature frames the orientation-preserving rigid-motion problem as reducing to constrained orthogonal Procrustes; the accessible arXiv abstract for the NIST Journal article states this reduction explicitly. Source: https://arxiv.org/abs/1902.03138 lines 37-40.
- This matches MDAnalysis `_fit_to`, which applies `R (X - Xbar) + Xbar_ref`, and the local h5-reader `ComputeSubsetTransform`, which sets `T = cr - R * cc`. Sources: https://docs.mdanalysis.org/stable/documentation_pages/analysis/align.html lines 586-594; `src/app/FitTargetMath.h:227`.

## Recommendation

1. **Reference: frame 0 vs iterative mean.** Mature viewers and command-line tools commonly fit to a chosen reference frame/structure: PyMOL `intra_fit` targets a state, VMD/RMSDTT targets a reference molecule/frame, GROMACS `trjconv -fit` targets the structure file, and ChimeraX `align` targets explicit reference atoms. Sources: https://pymol.org/pymol-command-ref.html lines 3203-3223; https://www.ks.uiuc.edu/Research/vmd/plugins/rmsdtt/ lines 53-64; https://manual.gromacs.org/current/onlinehelp/gmx-trjconv.html lines 485-486; https://www.rbvi.ucsf.edu/chimerax/docs/user/commands/align.html. However, iterative/converged mean reference is also established in analysis-grade tools: MDAnalysis `iterative_average` and Theseus ML simultaneous superposition. Sources: https://docs.mdanalysis.org/stable/documentation_pages/analysis/align.html lines 513-529; https://www.theobaldlab.org/theseus lines 60-77. Planned fix B is supported, but should not be represented as universal viewer default.

2. **Centroid pinning.** The centroid-pinned form is universal in the sources checked: MDAnalysis documents `R (X - Xbar) + Xbar_ref`; Eigen/Umeyama computes translation from destination mean minus rotated source mean; h5-reader's raw Kabsch already does `T = cr - R * cc`. Sources: https://docs.mdanalysis.org/stable/documentation_pages/analysis/align.html lines 586-594; https://eigen.googlesource.com/mirror/+/refs/heads/master/Eigen/src/Geometry/Umeyama.h lines 103-141; `src/app/FitTargetMath.h:227`. Planned fix A is not contradicted; it is the standard form.

3. **Fit selection and weighting.** The standard cure for all-atom precession is not temporal smoothing; it is choosing a stable fit set and/or robust weighting/rejection. Examples: VMD exposes trace/backbone/no-H selections and weights; MDAnalysis exposes `select` and `weights`; PyMOL `align`/`super` reject outliers; ChimeraX Matchmaker excludes bad fit regions during iteration; Theseus downweights variable regions. Sources: https://www.ks.uiuc.edu/Research/vmd/plugins/rmsdtt/ lines 35-52; https://docs.mdanalysis.org/stable/documentation_pages/analysis/align.html lines 252-258; https://pymol.org/pymol-command-ref.html lines 368-416 and 6891-6917; https://www.cgl.ucsf.edu/chimerax/docs/user/tools/matchmaker.html; https://www.theobaldlab.org/theseus lines 73-77.

4. **Temporal smoothing.** I found no mature-tool evidence for windowed smoothing of the rigid-body alignment transform as the primary stabilisation method. The documented smoothing operations are coordinate filters or coordinate interpolation: PyMOL `smooth` averages coordinate states, and GROMACS points to low-pass filtering for smooth movies. Sources: https://pymol.org/pymol-command-ref.html lines 6602-6626; https://manual.gromacs.org/current/onlinehelp/gmx-trjconv.html line 550. Recommendation: after centroid-pinned fitting to a stable mean/core reference, make transform smoothing off by default or strictly cosmetic; if retained, smooth a full centroid-consistent transform and recompute/apply translation as `T = cref - R * ccur`, not as an independently averaged world-space translation.

5. **Contradictions to A+B.** No source contradicts centroid pinning. The only tension is that frame/chosen-reference fitting remains the common viewer default, while iterative-mean fitting is better documented in analysis tools (MDAnalysis, Theseus) than in basic trajectory viewers. Sources: https://www.pymol.org/tutorials/moviemaking/ lines 153-158; https://docs.mdanalysis.org/stable/documentation_pages/analysis/align.html lines 513-529; https://www.theobaldlab.org/theseus lines 60-77.

# Applied-Maths Audit: Rediscover Analysis

> **Historical audit record — not current truth (trued 2026-06-04).** Superseded
> by the June 4 maths walk/checklist/fix postmortems; preserve as provenance.

Date: 2026-06-01. Scope: read-only audit of the rediscover emitters, fitters,
and emitted CSV/NPY substrates. I did not inspect any parallel-agent audit output.
`learn/` / Stage-1 normalisation is taken as archived and unverified, per brief.

## Findings And Correct Treatments

1. Buckingham E-field is a valid local-axis scalar projection, not an unframed
   equivariant scalar.

   Evidence: the C++ emitter builds a per-atom backbone frame and stores
   `efield_local = frame.ToLocal(eLab)` and `E_proj = efield_local.z()`
   (`src/rediscover/BuckinghamEfield.cpp:161-198`). The local-frame z axes are
   chemistry-defined and signful by stratum, e.g. N: N->CA, C: C->O, O: O->C,
   HA: CA->HA (`src/rediscover/LocalFrameBasis.h:110-154`,
   `src/rediscover/LocalFrameBasis.cpp:80-135`).

   Correct treatment: report the scalar Buckingham model as
   `sigma_iso_s = a_s E_z(local,s) + b_s |E|^2`, per target stratum and per
   explicit frame convention. Do not describe `E_z` as a standalone invariant of
   an l=1 field; the reference axis is the atom local frame. Coefficient signs
   are frame-convention signs and should not be compared across strata without
   restating the z-axis convention.

2. Buckingham T1 is emitted but not mathematically audited; it is silently
   outside the fitted result.

   Evidence: the sink writes raw `dft_total_T1_*` columns and a T1 sidecar, but
   marks the status as `unverified_emit_only` (`src/rediscover/BuckinghamEfieldSink.cpp:87-89`,
   `src/rediscover/BuckinghamEfieldSink.cpp:149-155`). The fitter states that
   T1/T2 payloads are ignored and deletes the loaded field sidecar after a shape
   check (`src/rediscover/analysis/buckingham_efield_t0.py:13-14`,
   `src/rediscover/analysis/buckingham_efield_t0.py:212-222`).

   Correct treatment: first pin the T1 convention. The rediscover decomposition
   defines T1 as the antisymmetric pseudovector `(A_yz, A_zx, A_xy)`
   (`src/rediscover/SphericalBasis.cpp:13-17`); the model type notes this is
   axial/even parity, not polar (`src/model/Types.h:523-526`). A polar electric
   field to axial T1 map is not an O(3)-equivariant scalar multiple without a
   pseudoscalar or local chiral axes. With the atom local frame present, fit and
   report a local-frame vector response, e.g. `T1_local = A_s E_local`, or a
   symmetry-restricted submodel, after emitting verified local T1. Until then:
   no T1 null/positive claim.

3. APBS field units and FF14SB charge-field units are mixed.

   Evidence: APBS E-field and EFG are catalogued as `V/Angstrom` and
   `V/Angstrom^2` (`src/rediscover/Catalog.cpp:76-80`). The broad FF14SB field
   reducer computes `field += (-q/r^3) * disp_local` without the Coulomb
   prefactor (`src/rediscover/BroadBackboneSink.h:96-128`), while the charge T2
   EFG kernel uses `kCoulombKe = 14.3996` (`src/rediscover/BroadBackbone.cpp:31`,
   `src/rediscover/BroadBackbone.cpp:118-122`).

   Correct treatment: either multiply the FF14SB field by the same electrostatic
   prefactor and emit it in `V/Angstrom`, or label it explicitly as `e/Angstrom^2`.
   Correlations are scale-invariant, but fitted coefficients and APBS-vs-FF14SB
   comparisons are not.

4. EFG 2e->2e Schur scalar gate is implemented, but the emitted target/features
   are lab-frame only; the within-atom de-meaning is therefore not the intended
   local-frame dynamic decomposition.

   Evidence: the EFG emitter intentionally passes an invalid local frame to
   `BuildTarget` (`src/rediscover/EfgFeature.cpp:83-87`), emits APBS EFG in the
   H5/library T2 basis (`src/rediscover/EfgFeature.cpp:111-115`), and writes
   `target.total_decomp.T2`, i.e. lab-frame DFT T2
   (`src/rediscover/EfgFeatureSink.cpp:145-146`). The fitter maps feature and
   target with the frozen C matrix (`src/rediscover/analysis/equiv_t2_efg_e3nn.py:194-195`)
   and implements Schur's scalar gate as `pred = g(|EFG|) * feature`
   (`src/rediscover/analysis/equiv_t2_efg_e3nn.py:148-160`). It then de-means
   the lab-frame target per atom (`src/rediscover/analysis/equiv_t2_efg_e3nn.py:199`)
   and de-means predictions inside the forward pass (`src/rediscover/analysis/equiv_t2_efg_e3nn.py:156-160`).

   Correct treatment: emit and fit `efg_feature_local_T2` and
   `efg_target_local_T2` using the same backbone frames used by broad/Buckingham,
   plus a `local_frame_valid` flag. Then the scalar gate can be tested on
   co-rotating components. If the desired model is truly lab-frame
   `T2_lab = g | EFG_lab`, then the static atom baseline must also rotate in lab
   and cannot be removed by subtracting a lab-frame component mean.

5. The frozen library-to-e3nn T2 basis is mathematically sound in this substrate.

   Evidence: rediscover's library T2 basis is Frobenius-isometric in order
   `[xy, yz, zz, xz, xx-yy]` (`src/rediscover/SphericalBasis.cpp:27-35`).
   `get_C()` returns the frozen orthogonal map and `lib_to_e3nn` is a constant
   matmul (`src/rediscover/analysis/change_of_basis.py:162-187`). My numerical
   check found `max |C.T C - I| = 1.11e-16`, determinant `1.0`, and norm
   preservation errors <= `1.8e-13` on the emitted broad/EFG T2 arrays.

   Correct treatment: keep using the frozen C. Do not rederive/project tensors
   in model paths.

6. Broad-backbone axis handling is now present and mostly correct, but the fitter
   still has stale schema-gap prose and train/test leakage in normalisation.

   Evidence: the broad source schema includes ring normals and bond axes
   (`src/rediscover/BroadBackboneSink.cpp:87-96`), and writes them for ring/bond
   rows (`src/rediscover/BroadBackboneSink.cpp:175-195`). The fitter detects the
   axis columns (`src/rediscover/analysis/equiv_t2_backbone_e3nn.py:207-215`) and
   builds `Y2(disp)`, `Y2(axis)`, and the `1o x 1o -> 2e` cross term with e3nn
   (`src/rediscover/analysis/equiv_t2_backbone_e3nn.py:381-391`). My scan of all
   emitted axis rows found 384,672 ring normals and 5,024,736 bond axes, with no
   norms outside [0.999, 1.001]. However, the file header still says those axes
   are not emitted (`src/rediscover/analysis/equiv_t2_backbone_e3nn.py:45-65`).
   The same fitter computes per-kind feature means/stds over all sources before
   the frame split (`src/rediscover/analysis/equiv_t2_backbone_e3nn.py:375-377`).

   Correct treatment: treat the current axes path as active and remove the stale
   caveat from reports. Compute normalisation statistics from train groups only,
   per stratum and source kind, then apply them to held-out groups. This is not a
   physics recompute; it is leakage-free preprocessing of emitted features.

7. The equivariant T2 fitters' per-atom de-meaning leaks test information and
   turns leave-atoms-out into a within-component diagnostic, not new-atom
   prediction.

   Evidence: broad target de-meaning uses all frames of each atom before the
   split (`src/rediscover/analysis/equiv_t2_backbone_e3nn.py:327-331`), and the
   model subtracts prediction means over all groups in `forward`
   (`src/rediscover/analysis/equiv_t2_backbone_e3nn.py:279-284`). EFG has the
   same pattern (`src/rediscover/analysis/equiv_t2_efg_e3nn.py:156-168`,
   `src/rediscover/analysis/equiv_t2_efg_e3nn.py:199-200`). Broad LOAO further
   subtracts the held-out atom's own target/predicted means
   (`src/rediscover/analysis/equiv_t2_backbone_e3nn.py:478-486`).

   Correct treatment: for frame-split within-atom scoring, estimate atom means
   and feature normalisers on train frames only. For leave-atoms-out, either
   report it explicitly as "shape after oracle centering of the held-out atom" or
   add a between-atom model that predicts the held-out atom baseline from emitted
   atom/source summaries. Do not call oracle-centered LOAO an absolute
   held-out-atom predictor.

8. Charge T2 distillation uses the field radial power (`q/r^2`) where the emitted
   fixed T2 charge kernel uses the EFG radial power (`q/r^3`).

   Evidence: the fixed broad charge T2 kernel is an EFG-like traceless tensor
   `q * (3 d d / r^5 - I/r^3)` times the Coulomb prefactor
   (`src/rediscover/BroadBackbone.cpp:107-123`). The distillation comparator for
   charge instead tests `gate_vs_q_r^-2`
   (`src/rediscover/analysis/backbone_distill_evidence.py:310-339`), and PySR
   exposes only `r` and `source_q_e`
   (`src/rediscover/analysis/backbone_pysr_distill.py:127-128`).

   Correct treatment: for a T2 point-charge/EFG law, the radial gate multiplying
   unit `Y2(disp_hat)` should be proportional to `q/r^3` (up to the chosen
   normalisation and Coulomb prefactor). Reserve `q/r^2` for vector electric
   field / scalar Buckingham T0/T1 analyses.

9. Charge "dipole" is an origin-dependent descriptor, not a translation-invariant
   physical shielding law.

   Evidence: the charge-dipole cell emits only `mu = sum q_i disp_local`
   (`src/rediscover/ChargeDipoleNeighborhood.cpp:188-191`), and the downstream
   script fits only `mu_norm` or `mu_xyz`
   (`src/rediscover/analysis/look_charge_dipole.py:122-126`). The manifest states
   the charge-multipole origin is the target atom
   (`src/rediscover/OutputManifest.cpp:44-48`).

   Correct treatment: if reporting multipoles, specify origin, source exclusion,
   net charge, dipole, and traceless quadrupole convention. For electrostatic
   shielding mechanisms, prefer the directly emitted field and EFG sums:
   `E = -k_e sum q_i disp_i/r_i^3` and
   `EFG = k_e sum q_i (3 disp_i disp_i/r_i^5 - I/r_i^3)`, converted to local T1/T2
   or scalar projections as appropriate. `mu` can be a descriptive target-origin
   moment, not a field law or a null result for charge response.

10. The de-circularising evidence document is stale relative to the current
    broad substrate.

    Evidence: `BACKBONE_LAW_EVIDENCE.md` says fixed literature T2 kernels were
    not emitted (`src/rediscover/analysis/BACKBONE_LAW_EVIDENCE.md:122-128`).
    The current broad sink emits aggregate total/ring/bond/charge T2 kernel
    sidecars (`src/rediscover/BroadBackboneSink.cpp:139-147`) and writes their
    CSV/NPY payloads (`src/rediscover/BroadBackboneSink.cpp:219-249`). The emitted
    `/tmp/rdc-broad-backbone-axes/broad_backbone_literature_kernel_t2_correlation.csv`
    contains the current de-circularised correlations.

    Correct treatment: use the emitted fixed-form kernel T2 arrays for the
    de-circularised correlation check. Report it as "fixed functional form,
    no fitted source geometry", not as a fixed ppm coefficient unless a literature
    shielding coefficient is also fixed.

11. Current R2/correlation reporting is incomplete for T2.

    Evidence: broad T2 R2 pools all five components into one denominator
    (`src/rediscover/analysis/equiv_t2_backbone_e3nn.py:287-288`) and reports a
    separate `|T2|` correlation (`src/rediscover/analysis/equiv_t2_backbone_e3nn.py:449-451`).
    The literature-kernel table already has `component_r`, `magnitude_r`, and
    `mean_row_cosine`, showing that magnitude alone can disagree with component
    alignment.

    Correct treatment: for every T2 result, report pooled Frobenius R2, per-T2
    component R2/correlation, vector row cosine or component correlation, and
    `|T2|` correlation. A good magnitude correlation with poor component alignment
    is not a tensor law.

## Complete Variance-Decomposition Design

Do this per mechanism, target (`sigma_iso` and T2), stratum
(`N/CA/C/O/HN/HA`, with `HA2/HA3` flagged separately when retained), and split.
Do not use the archived Stage-1 normalisation wholesale; normalise each component
and each decomposition term on the relevant train data.

For target component `c` of atom `i` frame `t`, let `y_itc` be scalar sigma or a
T2 component in the orthonormal library/e3nn basis. Let `n_i` be frames for atom
`i`, `mu_c` the grand mean, and `bar_y_ic` the atom mean.

Descriptive variance identity:

```text
SS_total,c   = sum_i sum_t (y_itc - mu_c)^2
SS_between,c = sum_i n_i (bar_y_ic - mu_c)^2
SS_within,c  = sum_i sum_t (y_itc - bar_y_ic)^2
SS_total,c = SS_between,c + SS_within,c
```

Fit/report three models:

```text
between:  bar_y_ic - mu_c       ~ B_m,c(atom/source summaries)
within:   y_itc - bar_y_ic      ~ W_m,c(frame-varying emitted features)
combined: y_itc                 ~ mu_c + Bhat_i,c + What_it,c
```

For leakage-free frame-split reporting, estimate `bar_y_ic`, feature means/stds,
and all ridge/MLP normalisers from train frames only; apply them to held-out
frames. For atom-split reporting, either predict `Bhat_i,c` for held-out atoms
or label the metric as oracle-centered within-shape only.

Normalisation and metrics:

```text
R2_between,c = 1 - sum_i n_i (bar_y_ic - yhat_B_ic)^2 / SS_between,c
R2_within,c  = 1 - sum_i sum_t (yW_itc - yhat_W_itc)^2 / SS_within,c
R2_total,c   = 1 - sum_i sum_t (y_itc - yhat_itc)^2 / SS_total,c
```

For T2, report per-component values plus pooled Frobenius versions where the
denominator is summed over components. Also report component correlation,
row-cosine, and `|T2|` correlation. Report variance shares
`SS_between/SS_total` and `SS_within/SS_total` so a mechanism can be
static-strong, dynamic-strong, both, or neither. Do not infer the axis from one
case: the current fixed charge kernel is dynamic-strong for N T2, while
Buckingham sigma_iso is near-null on the same within axis.

Effective N:

```text
N_eff_between = number of atoms with valid data
N_eff_within,c = sum_i n_i * (1 - rho_i,c) / (1 + rho_i,c)
```

Clip each AR(1) term to `[1, n_i]`; report both atom count and row-equivalent
`N_eff_within,c`. For T2, compute this per component and optionally for the
within-vector norm; keep `HA2/HA3` thin flags.

## Solvation Caveat

The DFT target is r2SCAN/def2-SVP CPCM(water), protein-only. The electrostatic
kernels mix FF14SB vacuum Coulomb sums, APBS Poisson-Boltzmann fields/EFGs, and
some explicit-MD hydration-derived signals. This mismatch is a plausible reason
that electrostatic kernels underperform or disagree in coefficient scale. It is
a disclosable interpretation caveat, not an applied-math fix or a work item for
this audit.

## Numerical Checks Run

Read-only checks on emitted CSV/NPY data:

- Frozen C: `max |C.T C - I| = 1.11e-16`, determinant `1.0`; norm preservation
  max error was `8.9e-16` on random vectors, `1.7e-13` on broad target-local T2,
  and `1.7e-13` on EFG target T2.
- Sidecar alignment: broad field sidecar vs CSV max error `5.0e-10`;
  Buckingham field sidecar vs CSV max error `5.0e-12`; EFG feature magnitude
  sidecar vs CSV max error `5.0e-12`.
- Axis vectors: all emitted ring normal and bond axis norms were within
  [0.999, 1.001].
- Variance shares from emitted targets show de-meaning discards material signal.
  Broad local T2 between/within shares include N `0.465/0.535`, CA `0.547/0.453`,
  C `0.333/0.667`, O `0.343/0.657`, HN `0.582/0.418`, HA `0.682/0.318`.
  Buckingham sigma_iso between/within shares include N `0.771/0.229`, CA
  `0.543/0.457`, C `0.081/0.919`, O `0.235/0.765`, HN `0.438/0.562`, HA
  `0.526/0.474`.
- Lag-1 check: broad local-frame T2 component median lag-1 was low/moderate
  (`N 0.074`, `CA 0.087`, `C 0.052`, `O 0.052`, `HN 0.162`, `HA 0.271`), so the
  random frame split is defensible there as a correlational gate. EFG lab-frame
  T2 component median lag-1 was high (`N 0.861`, `CA 0.761`, `C 0.866`,
  `O 0.808`, `HN 0.782`, `HA 0.783`), so the current EFG lab-frame split is not
  decorrelated; emit local EFG/target or use blocked splits.
- Buckingham scalar OLS on `/tmp/rdc-buckingham` reproduced near-null
  within-atom sigma results: N `R2=0.005`, CA `0.035`, C `0.016`, O `0.003`,
  HN `0.104`, HA `0.021`. The emitted T1 norm was not zero
  (median `4.87`, p95 `33.6` ppm in the unverified convention), confirming T1 is
  incomplete rather than absent.

No production code changed. Nothing committed. This report is
`src/rediscover/APPLIED_MATHS_AUDIT_codex.md`.

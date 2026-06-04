# Maths Audit Checklist (2026-06-04) — the open agenda for "are the numbers real?"

**Disposition: VERIFY TOGETHER (lead + Claude, NOT codex); do NOT fix piecemeal.** The unified combine +
every LOAO/between number are PROVISIONAL until this closes. A SECOND audit pass after, if the first
surfaces more. Sources: `POSTMORTEM_MATHS_WALK_2026-06-04` (adversarial walk — mostly sound + 3 issues),
the COMBINE NOTE (NOW.md / STATE.md), the `PHYSICS_ARCHITECTURE_UNIFICATION` fold-ins.

## The 3 concrete issues (maths-walk)
1. **DFT dia+para validated on T0 only, not T1/T2** (`src/io/DftShieldingLoader.cpp:95`). Trace equality ≠
   component-level dia+para=total. Affects the dia/para SPLIT targets, NOT the total-T2 combine (which we
   fit). → add a full-tensor dia+para=total check; trust total, distrust the split until then.
   **[RESOLVED-BENIGN 2026-06-04 — `POSTMORTEM_DIAPARA_CHECK`]:** full-tensor check on build4 targets shows
   dia+para=total holds at ALL components (T0/T1/T2) to ORCA print-rounding (max 1.0e-3 ppm; T2 1.63e-3 = the
   isometric-basis projection of 1e-3, not drift; inverse 3×3 max 1.0e-3). The split targets are SOUND; the
   T0-only loader validation was sufficient in practice; total never implicated. (A full-tensor check in the
   loader is belt-and-suspenders, not required.)
2. **[CODE-CONFIRMED — the consequential one] the "LOAO" path centers the held-out atom by its OWN mean**
   (`analysis/stage2_law_fits.py:663-665, 1681-1682`; fresh-look verified) → it measures WITHIN-atom
   modulation, NOT between-atom recovery. **CONSEQUENCE: 1P9J has NO clean between axis; the 720-WT is the
   ONLY between instrument.** Every "LOAO/between" number — unified 0.26, charge 0.38, ring 0.17 — is a
   mislabeled WITHIN number; do NOT quote as between. Decide together: is a true LOAO (predict the held-out
   atom's MEAN, not its modulation) even the right between-test on ONE protein, or is between deferred
   entirely to the 720-WT?
   **[RESOLVED 2026-06-04 — `POSTMORTEM_TRUE_LOAO` / `b6e4d2e`]:** true-LOAO run → 1P9J TRUE between-atom
   recovery is ~null (charge 0.036, 0.03-class; ring −1.0; unified −105, overfit, null p0.70/z0.06). Between
   numbers RETRACTED; between / transferability / combine-DEPTH defers ENTIRELY to the 720-WT. **CODE FIX**
   landing (`bh5f0e7ve`): the LOAO/between path → true between-atom recovery (reuses the validated `true_loao_*`
   machinery); full re-fit HELD per lead.
3. **ring/broad e3nn: de-mean over ALL groups + unpurged random frame splits** (`analysis/equiv_t2_e3nn.py`,
   `analysis/equiv_t2_backbone_e3nn.py`) → possible leakage. The EFG e3nn path (`analysis/equiv_t2_efg_e3nn.py`)
   is clean (blocked/purged) = the template to align to. Affects the equivariant-path numbers in "3 paths agree."
   **[FIXED 2026-06-04 — `POSTMORTEM_E3NN_PROTOCOL_FIX` / `b804wd9rr`, code-only]:** ring/broad e3nn aligned to
   the clean EFG protocol (shared `analysis/e3nn_protocol.py`: blocked/purged split, train-only
   centering/normalization, `center_mask=g_tr`); all-group de-mean + random unpurged splits REMOVED; structural
   checks green, NO e3nn re-run. The clean-vs-leaky recovery VERDICT pends the HELD e3nn re-run.

## The combine: VALIDITY vs ATTRIBUTION (do not conflate)
- **VALIDITY** (could the combine be a mirage): regularization + effective DOF, the right NULL for a combine
  (structured vs shuffle-target), basis-invariance (change_of_basis / per-type sums), held-out recovery.
  **Collinearity is NOT a validity gate** — it is EXPECTED + confirmatory for shadows-of-one-object (recovery
  = projection onto the spanned column space, invariant to collinearity); DESCRIBE it.
- **ATTRIBUTION** (collinearity-affected, soft): the per-shadow drop-one marginals ("field +0.198, mc +0.169").
  Report RECOVERY as the robust claim; describe the collinearity; treat per-shadow attribution as soft.

## Fold-ins (PHYSICS_ARCHITECTURE_UNIFICATION)
- Grade statistical position against the **fixed-eigenstructure null** (the `(3cos²−1)` magic-angle shape,
  node at cos²θ=1/3), not a generic fit family.
- Make path-agreement a **coefficient agreement in PHYSICAL UNITS** on a clean stratum, NOT "the three R²s close."

## The settled within-axis results the audit does NOT threaten
charge q/r³ (within 0.28, z451), ring current (within 0.28, z155), the unified combine's WITHIN recovery
(0.43, z263, field+McConnell, not charge-in-a-coat). The defensible cookies. The audit is about the BETWEEN
story + the combine's DEPTH/attribution claim — not these.

## Reporting standard (binds the audit)
Statistical position + determinability, never survives/overfit. Probability where N earns it; case-study
where tiny (ring-5). 1P9J = WITHIN instrument. ~0.03 ≈ null; ~0.2 = something-or-trash decided by
determinability. ([[feedback_law_as_statistical_position]], [[feedback_transparent_cutoffs]].)

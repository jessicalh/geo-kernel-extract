# Deferred ledger — the couponed items

**Purpose.** The durable home for everything we **deliberately did NOT fix in place** — deferred,
"not ours," parked, passed to another author, or flagged for a later pass. The rule: **every time a coupon
is spent (a thing flagged-and-set-aside instead of fixed), it lands here**, with what / why / owner /
where. Nothing important should live only in chat. Jessica is told about each as it's added.

Status legend: **[parked]** needs Jessica's vet before firing · **[reader-pass]** owned by the
comprehensive h5-reader downstream-fix pass ([[project_reader_downstream_fix_no_cruft]]) · **[design]** a
human design call · **[author]** belongs to another author's concurrent work · **[follow-up]** a near-term
to-do on this work · **[robustness]** low-priority hardening.

---

## Parity / labels

- **[parked] The `1o`→`1e` parity sweep.** Every *shielding-kernel* H5 time-series writer still advertising
  `0e+1o+2e` is mislabelled — a shielding tensor's antisymmetric `T1` is the **axial** pseudovector →
  `1e`, not polar `1o`. Ring (BS/HM/ring-chi) and McConnell are now fixed; the **stragglers** are
  **π-quadrupole, dispersion, H-bond, Larsen, Tripeptide** (`PiQuadrupoleShieldingTimeSeriesTrajectoryResult.cpp:104`,
  `DispersionShieldingTimeSeriesTrajectoryResult.cpp:104`, `HBondShieldingTimeSeriesTrajectoryResult.cpp:104`,
  + Larsen/Tripeptide writers). The genuine **polar E-fields stay `1o`** (Coulomb/MOPAC/APBS/water,
  AIMNet2 charge-response-gradient). *Confirmed still-needed by the ring-review divergence (codex caught
  it; opus rationalised it away; codex was right). Scope question: several stragglers are CAGED kernels —
  a parity-label fix is a correctness fix, not a feature update, but it touches caged files; Jessica
  decides whether the sweep includes them or leaves them caged-mislabelled.* 2026-06-07.
- **[author] `mopac_dipole_components` is internally contradictory** — `irreps="1e"` with `parity="odd"`
  (`python/nmr_extract/_catalog.py:208`). A dipole is polar → should be `1o`/odd. Same class of bug the
  ring patch fixed on `bs_total_B`. Belongs to the **concurrent MOPAC-full work**, not the ring change.
  2026-06-07 (opus ring review).

## Reader-pass inbox (the h5-reader comprehensive downstream-fix pass owns)

- **[reader-pass] h5-reader stale generated catalog.** `QtFieldCatalog.gen.h` is generated from
  `python/nmr_extract/_catalog.py` (2026-05-26) and still maps old stems (`coulomb_shielding`,
  `mopac_coulomb_shielding`) — the new `coulomb_efg` etc. are invisible until regenerated
  (`python3 h5-reader/scripts/gen_field_catalog.py`). Generated + 1:1 provenance + graceful-absent → a
  clean wording/regen, not a functional break.
- **[reader-pass] Internal C++ `*_shielding_contribution` field names** (`coulomb_shielding_contribution`,
  `mopac_coulomb_shielding_contribution`) — the emitted NPYs are renamed `*_efg`, but the internal field
  names (and the **MOPAC-vs-FF14SB reconciliation diagnostic** that reads them) still say "shielding."
  Renaming ripples into the reconciliation result → deferred whole.
- **[reader-pass] Legacy `ArraySpec` wrappers** kept for reading pre-rename dirs (`coulomb_shielding`,
  `mopac_coulomb_shielding`, `mc_shielding`, `mc_category_T2`, `mopac_mc_*` → `*_legacy`, `is_feature=False`).
  Cruft by the no-cruft rule; harmless, swept in the reader pass.
- **[reader-pass] Docs / golden baselines** still carry old names; the **golden/smoke re-bless** is a
  separate wholesale act (below).

## Design calls (Jessica)

- **[design] B-field diagnostics: producer vs reader.** `bs_ring_B_field` / `bs_ring_B_cylindrical` /
  `hm_ring_B_field` are sparse `(P,3)`, ~`24P` bytes each, **~15.5% of the `ring_contributions` payload**,
  scaling on the atom-ring-pair axis. The reader re-evaluates BS/HM closed-form, so it can **recompute**
  them for vetting. Emitted now (per "all five"), but a real revisit if trajectory/fleet storage pressure
  matters: keep emitting, or move to reader-recompute. 2026-06-07.
- **[design] Golden/smoke wholesale re-bless.** Deferred until binary-compat matters again (release /
  fleet). The smoke baseline is entangled across EFG + McConnell + MOPAC-full + broad_backbone, so the
  re-bless is one wholesale act, not per-patch ([[feedback_dont_overhold_builds_bless]]).

## Ring to-builds (follow-up on this work)

- **[RESOLVED 2026-06-07] HM name → "HM-style" (option c); bond-sum benchmark DROPPED.** Earning the
  literal "Haigh–Mallion" name by implementing the bond-sum and self-correlating it against our
  surface-integral is **self-test ≠ validation** ([[feedback_consistency_not_validation]]) — Jessica
  killed it ("using code we didn't produce to test the thing we wrote"; the two paths landed at −2.5× + a
  sign flip). We keep the surface-integral, labelled **"HM-style / HM-inspired tensor."** Re-openable later
  via (a) acquire HM 1979/1980 + test our T0 vs *their* published factors, or (b) an algebraic-equivalence
  proof — neither is held/done; not now.
- **[design — open] #4 BS↔HM two-path: regression guard, framing TBD.** The agreement is *internal
  consistency* (both paths ours), valuable as a CTest-on-fixture **regression guard** (flip the report-only
  test to assert; derived 0.36/1.0 band) — but **not** a validation of the physics. Open call: also *emit*
  it as a labelled consistency-diagnostic, or keep test-only? Stale hazard: the test reconstructs HM on a
  different sign path than production (`test:311-315` vs `HaighMallion:293-296`) — any emit must use
  production tensors.

## Part-1 decisions

- **[follow-up] X–H ablation.** The `mcconnell_include_xh_sources` toggle (default 1) is a built scaffold;
  Part 1 decides keep/drop the C–H/N–H/O–H source bonds vs DFT.

## McConnell rhombic (value-gated; spec landed 2026-06-07)

Full scope: `kernel_design/mcconnell_rhombic_spec.md`. Jessica wants it; additive, drops nothing.

- **[follow-up] The rhombic shape build (#5).** Extend `Q̂` for C=O + C=C with the sp²-plane-normal
  rhombic block (`McConnellResult.cpp:230`, axial today); emit the unit shape, scale learned. Buildable
  without the cited value (geometry, additive). Pairs with #3/#4 as the McConnell+ring tensor-richness
  pass; schedule after #2 (same emit surface). Build routes to **codex**.
- **[RESOLVED 2026-06-07] Rhombic C=O sourced + pinned; C=C dropped.** C=O rhombic LANDED + sign-resolved:
  Hooper & Kaiser 1965 (10.1139/v65-318) + Abraham & Ainger 1999 (10.1039/A808908F) held + provenance-
  verified. **Pinned: Hooper EF-corrected acetamide A `(−5.4, +4.0, −14)` ×10⁻⁶ cm³ mol⁻¹**, Abraham-
  anchored sign (high-conf), magnitude medium → **sensitivity-reported** (lead sign-off). **C=C DROPPED**
  — moot (no non-aromatic C=C in the 20 residues; aromatic C=C is ring's, zeroed). Rhombic = C=O-only.
- **[follow-up — contingent] £40 Williamson–Asakura acquisition.** The pinned magnitude is medium-confidence
  (Hooper's EF-corrected fit poor). **If** the rhombic shows promise *and* the result proves magnitude-
  sensitive, spend ~£40 on the WA peptide primary to firm it. Don't pay until it matters; the sensitivity
  report defends the pin until then. (Jessica, 2026-06-07.)
- **[author — scholarship] `MCCONNELL_DCHI_LITERATURE.md` hygiene.** Its "Source Anchors" list unheld
  primaries (Hooper-Kaiser, Abraham, WA…) as if they back the table — relabel POINTER-ONLY and add the
  two genuinely-held axial anchors (Pauling 1979, Worcester 1978) it omits.

## Robustness (low-priority)

- **[robustness] B-field conditional-zero in partial pipelines.** `ConformationResult` writes the B-field
  rows for every `ring_neighbours` row regardless of whether BS/HM ran. Production always runs them first
  (real values); a *partial/custom* run (e.g. ring-chi without BS) would write default zeros
  (`ConformationResult.cpp:95`, `ConformationAtom.h:37`). Edge case, production-fine. 2026-06-07 (codex ring review).

## Fleet-size + .LGS (2026-06-08)

- **[design — real] Complete-emit heavy arrays blow the fleet budget.** COO arrays (mc_bond_neighbors 63MB,
  spatial_neighbors 35MB on Ubq) + MOPAC high-dim matrices are huge per-protein; at fleet scale (676)
  mc_bond_neighbors alone ~42GB vs the <15GB budget. Revisit (cap/top-K/opt-in) BEFORE any fleet run.
  Already excluded from golden + gitignored (BLESS_NOTES).
- **[follow-up] .LGS emitter (drafting).** KISS: one root (=--output), `<protein_id>_<timestamp>.lgs`,
  minimal/derived rest. Codex draft (CalcsetManifestEmitter; trajectory + single-pose; mutant may stub).
  Schema spec/CALCSET_MANIFEST.md v1; port h5-reader/tools/lgs_write.py. Review the draft next session.
- **[follow-up] Golden rename-orphan git rm.** 8 renamed-away baselines (coulomb_shielding, mc_category_T2/
  scalars/shielding + mopac_ versions) still in blessed -> smoke reports them "Missing from run". `git rm`
  them in a fresh context (deferred per lead, not late-session).

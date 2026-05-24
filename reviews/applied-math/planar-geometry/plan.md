# Fix plan — PlanarGeometryResult.{h,cpp}

## 1. Summary

The file tells a coherent story. The header is a four-quantity contract
(ω, sp2 pyramidalization, aromatic χ₂, Cremer-Pople pucker), each with a
convention-pinned doc block; `Compute` runs four numbered stages in physical
order; the four helpers are each self-contained with a units/sign signpost.
Both reviews agree there is **no algorithmic correctness bug**. The signs and
numbers are earned by iteration and traced cleanly to consumers (see §4) — the
pass moves **no numbers**.

What the pass **will** touch (readability only):
- Trim design-history / process-archaeology prose out of header and `.cpp`
  comments down to signposts (keeping the literature anchors and the
  hard-won normal-orientation rationale).
- Improve the worst internal-symbol opacity in the Cremer-Pople block
  (`cs`/`sn`/`R1`/`R2`/`Qcos`/`Qsin`) and a few local names elsewhere.
- Fix one inaccurate comment (`A == G` for planar sp2).
- Add the two direct includes (`<array>`, `<algorithm>`) the file uses
  transitively.

What the pass **will not** touch:
- No output names: `pyramidalization`, `omega_actual`, `omega_deviation`,
  `omega_is_xpro`, `aromatic_chi2`, `pucker_Q`, `pucker_theta` (NPY/H5/SDK
  contract) stay exactly as they are.
- No public method names (`OmegaActual`, `PuckerQ`, `PuckerTheta`,
  `AromaticChi2`, etc.) — these are read cross-result by the TR layer.
- No algorithm, no sign, no value, no refactor of stage structure.
- The `pucker_valid` diagnostic count semantics — see §4 / §5; this is a
  log-line label decision, not a number.

---

## 2. Review-finding ledger

Every finding from `codex.md` and `claude.md`, one row each. (No
`codex-correctness.md` present in this directory.)

### codex.md

| # | Finding (loc) | Disposition |
|---|---|---|
| C1 | h:6 cross-result-read block leads the header before the science | **adopted** → §3 E1 (trim to short note, keep below contract is optional; minimal: shorten) |
| C2 | h:16 "landed in the 2026-05-05 → 2026-05-08 topology slice" interrupts | **adopted** → §3 E2 |
| C3 | h:27/36-39 omega historical-correction prose forces current-vs-historical read | **adopted** → §3 E3 |
| C4 | cpp:107 Cremer-Pople block mixes derivation + warning + review history; move "adversarial review 2026-05-09" out | **adopted** → §3 E8 (drop the dated parenthetical only; keep the anti-parallel rationale) |
| C5 | cpp:153 degeneracy guard explanation too long | **declined** → claude.md (CL10) judges the literature anchor (Ho & Cornilescu Q≈0.35–0.45) earns its keep; I agree — it is grounded science, not archaeology. Trim only the dated "math-review MED-1, 2026-05-19" tag (§3 E10). |
| C6 | cpp:246 stage comments longer than code, esp. omega 278-309 | **adopted (partial)** → §3 E12 (trim omega block's process prose; keep the bond-graph-connectivity rationale, which is load-bearing) |
| C7 | cpp:70 `A,B,C,D,G,n` hide central/neighbours/centroid/normal | **adopted (partial)** → §3 E6. Rename `G`→`centroid`, `n`/`n_hat`→`normal`/`normal_hat`. **Decline** renaming `A,B,C,D` (CL4 calls them conventional + documented for an improper-dihedral helper; agree). |
| C8 | cpp:121 `R1`/`R2` → `sin_weighted_basis`/`cos_weighted_basis` | **adopted** → §3 E7 (use `sin_basis`/`cos_basis` to echo the comment's R'₁/R'₂) |
| C9 | cpp:142 `cs`/`sn` worst compression; neg-sine hidden | **adopted** → §3 E9 |
| C10 | cpp:150 `Qcos`/`Qsin` → `q2_cos_theta`/`q2_sin_theta` | **adopted** → §3 E9 (rename to `q2_cos`/`q2_sin`) |
| C11 | cpp:348 `ri` = aromatic ring idx (was residue idx in omega loop) | **adopted** → §3 E13 (`arom_i`) |
| C12 | cpp:377 `ri` = saturated ring idx | **adopted** → §3 E14 (`sat_i`) |
| C13 | cpp:264 `nb` too compressed for sign-sensitive ordering | **adopted** → §3 E5 (`neighbours`) |
| C14 | cpp:272 `pyr` hides signed-Å | **declined** → local, one-line scope, immediately assigned to `ca.pyramidalization`; the helper name already carries it. Low value, declined as churn. |
| C15 | cpp:35 `Dihedral` fuse → named intermediates `normal_123` etc. | **declined** → CL3.2 calls this the textbook b1/b2/b3→n1/n2/m1→atan2 sequence, signposted; agree it is clean as-is. The DihedralTS sibling uses the identical compact form; renaming here only diverges them. |
| C16 | cpp:134 pucker projection collapses z/phase/sign | **adopted** → covered by §3 E9 (named projection sums) |
| C17 | cpp:330 `WrapPi(omega - M_PI)` → named `deviation_from_trans` | **declined** → the result is stored directly into `omega_deviation_[ri]`; the lhs name already carries it. Adding a temp is fragmentation. |
| C18 | cpp:35 `Dihedral` → `SignedDihedralRadians` | **declined** → block comment states signed/radians/range; CL4 agrees clean. Renaming diverges from the bit-identical DihedralTS sibling helper of the same name (load-bearing parallel, see §4). |
| C19 | cpp:70 `Pyramidalization` → `SignedOutOfPlaneDisplacement` | **declined** → block comment states signed/out-of-plane/Å; rename is verbose churn, no consumer reads the file-local helper name. |
| C20 | cpp:100 `CremerPople5Ring` → `ComputeCremerPople5RingPucker` | **declined** → CL4.1 notes the `PuckerCP{Q,theta}` return struct carries the return; acceptable as-is. |
| C21 | cpp:177 `ThreeBondedNeighbours` → `SortedThreeBondedNeighbourAtoms` | **declined** → "sorted" is an implementation detail of the sign-stability mechanism, not the function's purpose (return the 3 neighbours). Will instead signpost the sort (§3 E11). |
| C22 | h:103 getter names clear | **no-op** (affirming) |
| C23 | h:2 file header more history than signpost | **adopted** → §3 E1–E4 (consolidated header trim) |
| C24 | cpp:49 WrapPi comment mostly history → `// wrap to [-pi,pi]` | **adopted (partial)** → §3 E4. Trim the dated "MED-2, 2026-05-19" tag but **keep** the bit-identity-with-DihedralTS rationale (CL5.1 — it is load-bearing, see §4). |
| C25 | cpp:66 "A == G for planar sp2" not reliable | **adopted** → §3 E6 (comment correctness fix; duplicate of C29) |
| C26 | cpp:187 neighbour-sort comment too process-heavy ("720-protein fleet") | **adopted** → §3 E11 |
| C27 | cpp:287 `feedback_huxley_data_discipline` / PATTERNS refs break through-line | **adopted** → §3 E12 (drop the policy-doc citations from the comment) |
| C28 | cpp:465 X-Pro write comment longer than needed | **declined** → it documents the serialized `omega_is_xpro` field semantics at the write site (consumer-facing contract: "consumer must apply the mask"). This earns its length at an output boundary. Low value to trim, risk of dropping contract prose. |
| C29 | cpp:67 comment correctness: A lies in neighbour plane, not at centroid | **adopted** → §3 E6 |
| C30 | cpp:177 missing `<array>` include | **adopted** → §3 E15 |
| C31 | cpp:193 missing `<algorithm>` include | **adopted** → §3 E15 |
| C32 | cpp:390 `pucker_valid` counts finite Q even when theta NaN | **adopted (label fix)** → §3 E16. This is a diagnostic-log label, not an output value (see §4). Rename the log token to say what it counts; do not change the count. |

### claude.md

| # | Finding (loc) | Disposition |
|---|---|---|
| CL1 | cpp:142 `cs`/`sn` read as generic temporaries until renamed 8 lines later | **duplicate** of C9 → §3 E9 |
| CL2 | cpp:145 `phi` reused as different angle (4π vs 2π) with no signpost on doubling | **adopted** → §3 E9 (add `// m=2 (2nd-harmonic) phase` signpost) |
| CL3 | h:6-14 cross-result block dense provenance at top of header | **duplicate** of C1 → §3 E1 |
| CL4 | cpp:142 `cs`/`sn` too terse | **duplicate** of C9 |
| CL5 | cpp:73,121-129 `n`,`R1`,`R2`,`n_hat`; `n`→`mean_plane_normal` | **adopted** → §3 E7 (R1/R2 via C8; n→`mean_plane_normal`, n_hat→`mean_plane_normal_hat`) |
| CL6 | cpp:70-71 `A,B,C,D` conventional, documented; no change | **adopted (as decline of rename)** → consistent with C7 partial |
| CL7 | cpp:319-320 `res_i`/`res_ip` ("ip"=i-plus-one) cryptic → `res_next` | **adopted** → §3 E12b (`res_next`) |
| CL8 | cpp:100-165 CremerPople5Ring three stages run flat; add captions | **adopted** → §3 E7/E9 (`// mean-plane normal`, `// (Q,θ) projection`) |
| CL9 | cpp:35-46 Dihedral clean | **no-op** (affirming; supports decline of C15/C18) |
| CL10 | cpp:152-164 well sequenced; guard at hazard | **no-op** (affirming) |
| CL11 | cpp:49-53 WrapPi comment half changelog → trim history, keep bit-identity | **duplicate** of C24 → §3 E4 |
| CL12 | cpp:153-160 MED-1 degeneracy comment long but lit anchor earns keep; leave | **adopted** → supports decline of C5; trim only the dated tag (§3 E10) |
| CL13 | cpp:107-120 normal block excellent signpost; model for the file | **no-op** (affirming; constrains C4 to dropping only the dated parenthetical) |
| CL14 | h:36-39 X-Pro "documentation error" history prose → drop to one-liner | **duplicate** of C3 → §3 E3 |
| CL15 | cpp:30-34,60-69,81-94 helper block captions right shape; clean | **no-op** (affirming) |
| CL16 | cpp:184 `b.atom_index_a == ai` possible sign/width mismatch | **declined (verified non-issue)** → `Bond::atom_index_a`/`_b` are `size_t` (src/Bond.h:15-16); `ai` is `size_t`. No sign/width comparison. See §4. |
| CL17 | cpp:193 sort makes sign build-stable but sign depends on index order not chemistry; confirm .h "first two neighbours" means post-sort | **adopted (comment clarify)** → §4 reason; §3 E6b: align the .h:43-45 wording so "first two neighbour vectors" reads as post-sort (lowest two atom indices). No code/sign change. |
| CL18 | cpp:354 confirm `chi[1]` is genuinely χ₂ not χ₁ | **declined (verified correct)** → `Residue::chi[4]` is 0-indexed (chi[0]=χ₁); every consumer (DihedralTS, DsspResult, ChiRotamer, DihedralBinTransition) reads chi[k] 0-indexed. `chi[1]`=χ₂. See §4. No change. |
| CL19 | cpp:379 6-ring saturated ring would vanish without a log line; acceptable | **declined (scope)** → reviewer self-classifies as no change needed given current scope (Pro pyrrolidine only). Agree; out of readability scope. Noted in §6. |
| CL20 | cpp:101 vs 379 double size guard ==5; harmless redundancy | **declined** → reviewer agrees fine; the inner guard makes the helper safe standalone. No change. |
| CL21 | no accumulation/bounds hazard in the four Compute loops | **no-op** (affirming) |

---

## 3. Edits that don't move numbers

All in `src/PlanarGeometryResult.{h,cpp}`. Comment trims keep every literature
anchor and every load-bearing rationale; they drop only dates, internal
review tags, ticket IDs, policy-doc citations, and changelog prose.

**Header (`.h`)**
- **E1** — h:6-14 — shrink the `CROSS-RESULT READ` block to a 2-3 line note:
  "Optional downstream readers (RingPuckerTimeSeriesTrajectoryResult,
  DihedralTimeSeriesTrajectoryResult) read these getters per-frame; they
  attach conditionally, so absence is captured by the source-attached gate."
  Drop the dated/PATTERNS-§17 framing and the attr-name list.
- **E2** — h:16-19 — keep "substrate carries the typed *classification*; this
  calculator carries the per-frame *deviation from canonical*"; drop the
  "landed in the 2026-05-05 → 2026-05-08 topology slice" clause.
- **E3** — h:36-39 — delete the "X-Pro NaN-fill described in earlier doc
  versions was a documentation error … corrected 2026-05-19 sweep" sentence.
  The live behavior ("ω emitted for X-Pro, use omega_is_xpro mask") is already
  stated above it.
- **E6b** — h:43-45 — adjust the pyramidalization sign sentence so "first two
  neighbour vectors" is unambiguously the two lowest-atom-index neighbours
  (post-sort), matching the code's `std::sort` + `nb[0],nb[1]` use. Wording
  only; the convention and sign are unchanged (see §4, CL17).

**`.cpp` — Pyramidalization helper**
- **E6** — cpp:65-68 — replace "For a perfectly planar sp2 site, A == G and
  the displacement is exactly zero" with: "For a perfectly planar sp2 site A
  lies in the neighbour plane, so `(A - centroid)·n̂` is zero." (C25/C29.)
- **E6** — cpp:70-77 — rename `G`→`centroid`, `n`→`normal`, `n_hat`→
  `normal_hat`. Keep `A,B,C,D` (documented improper-dihedral convention).

**`.cpp` — WrapPi**
- **E4** — cpp:49-53 — trim to: "Wrap to [-π, π] via `std::remainder`
  (single-call, IEEE round-half-to-even at ±π). Matches
  DihedralTimeSeriesTrajectoryResult::WrapPi bit-identically so omega_deviation
  agrees across the two producers." Drop the "historic while-loop" and
  "MED-2, 2026-05-19" history; **keep** the bit-identity rationale (§4).

**`.cpp` — Cremer-Pople helper** (the main readability win)
- **E7** — cpp:107-132 — keep the existing excellent normal-orientation
  comment (Eqs 11-12 + the anti-parallel trap). Rename `R1`→`sin_basis`,
  `R2`→`cos_basis`, `n`→`mean_plane_normal`, `n_hat`→`mean_plane_normal_hat`.
  Add a blank-line + `// mean-plane normal` caption to group the block.
- **E8** — cpp:120 — drop only the "(caught in adversarial review 2026-05-09)"
  parenthetical; keep the rest of the anti-parallel explanation verbatim.
- **E9** — cpp:134-151 — add `// (Q,θ) projection` caption; add a one-line
  `// m=2 (2nd-harmonic) phase, hence 4π j/5` signpost where `phi` doubles
  (CL2). Rename `cs`→`q2_cos_sum`, `sn`→`q2_sin_sum` (and keep the explicit
  `-z_j` so the negative-sine convention stays visible at the accumulation,
  matching the comment's `Q₂ sin θ₂ = -√(2/5) Σ …`). Rename `Qcos`→`q2_cos`,
  `Qsin`→`q2_sin`.
- **E10** — cpp:153-160 — keep the degeneracy-guard science and the Ho &
  Cornilescu 2000 Q≈0.35–0.45 Å anchor; drop only the "math-review MED-1,
  2026-05-19" tag.

**`.cpp` — ThreeBondedNeighbours**
- **E11** — cpp:187-192 — trim "Without this, a join of pyramidalization NPYs
  across the 720-protein fleet would silently mix sign conventions" to
  "Sort by atom index so the pyramidalization sign is build-stable and
  portable across structures (bond order in `bond_indices` is not stable)."
  Keep the mechanism; drop the fleet-anecdote framing.

**`.cpp` — Compute, stage 2 (omega)**
- **E12** — cpp:278-309 — keep the bond-graph-connectivity rationale
  (`BackboneSuccessor`, wrap-correct / insertion-codes / numbering-gaps /
  caps) — it is load-bearing science about why label heuristics are banned.
  Drop the `feedback_huxley_data_discipline` + PATTERNS.md citations (C27) and
  the "Without this gate the original loop walked raw ri+1 …" changelog
  sentence.
- **E12b** — cpp:319-320 — rename `res_ip`→`res_next` (and `res_i` may stay or
  become `res_curr`; minimal: `res_next`). Matches the `next_idx` local.

**`.cpp` — Compute, stages 3 & 4 (loop counters)**
- **E13** — cpp:348 — rename loop var `ri`→`arom_i` (aromatic ring index).
- **E14** — cpp:377 — rename loop var `ri`→`sat_i` (saturated ring index).

**`.cpp` — includes**
- **E15** — cpp:19-20 — add `#include <array>` (used at cpp:178) and
  `#include <algorithm>` (used at cpp:193) alongside `<cmath>`, `<limits>`.
  Self-containment; no behavior change.

**`.cpp` — diagnostic log label**
- **E16** — cpp:390/401 — the `pucker_valid` counter increments on finite `Q`
  even when `theta` is NaN (sub-amplitude degeneracy). This is a **log line**,
  not an output. Rename the printed token from `"pucker valid="` to
  `"pucker finite-Q="` (and optionally the local `pucker_valid`→`pucker_finiteQ`)
  so the diagnostic states what it counts. **Do not** change the increment
  condition — see §4 / §5.

---

## 4. Usage notes (the signs/values the reviews raised)

### Pyramidalization sign (CL17, C7, C26)
- **Reason discovered.** Sign is `(A - centroid)·n̂` where
  `n̂ = ((B-G)×(C-G)).normalized()` and B, C, D are the three bonded
  neighbours **sorted ascending by atom index** (`std::sort`, cpp:193). The
  sort exists so the cross-product orientation — and therefore the sign — is
  reproducible across builds and across proteins, because `Atom::bond_indices`
  ordering is not guaranteed stable. So the sign is a *stable, arbitrary-but-
  fixed* convention keyed to the two lowest-index neighbours, not a
  chemistry-derived absolute (e.g. it is not pinned to a specific substituent).
- **Consumers.** `ConformationAtom::pyramidalization` → `pyramidalization.npy`
  (catalog, float64, units Å, "0 for non-planar"). Read by `ui/MainWindow.cpp`
  (display only, signed Å) and `ui/RestServer.cpp` (`pyramidalization_A`).
  No consumer flips or reinterprets the sign; `learn/` uses it as a feature
  column. The TR layer does not carry pyramidalization (it is per-atom static
  emission, not in RingPuckerTS).
- **Producer/consumer agreement.** Agree. Every consumer treats it as a signed
  magnitude; none depends on an absolute chemical orientation. The convention
  is internally consistent and stable, which is exactly what the consumers
  need. **Coherent (expected).** The only action is the comment fix E6/E6b
  (the "A == G" line is wrong) and signposting the sort's purpose (E11).

### Cremer-Pople θ sign / negative-sine convention (C4, C9, C16, CL2)
- **Reason discovered.** `q2_sin_sum += -z_j·sin(4π j/5)` carries an explicit
  negative sine, matching the comment's `Q₂ sin θ₂ = -√(2/5) Σ z_j sin(…)`.
  Combined with the mean-plane normal built from the sin/cos-weighted basis
  (`R'₁×R'₂`, not the naive edge-cross accumulator), this fixes θ's phase to
  the canonical Cremer-Pople 1975 direction. The header CLAUDE.md record and
  the in-file comment both note the **θ inversion that was caught in
  adversarial review 2026-05-09** — i.e. this exact sign/normal-orientation
  pair is the *result* of iteration, the governing-prior case.
- **Consumers.** `pucker_theta.npy` (degrees, [0,360)); H5
  `/trajectory/ring_pucker_time_series/pucker_theta` via
  RingPuckerTimeSeriesTrajectoryResult, which **copies the per-frame value
  verbatim** (cpp:92) with no sign reinterpretation, and stamps
  `pucker_5ring_endvtwist = "theta mod 72 deg gives envelope/twist endo/exo"`.
  `h5-reader/QtSpecialBuffers.h` reads it for display. `learn/` consumes as a
  feature.
- **Producer/consumer agreement.** Agree. Consumers carry θ through unmodified
  and rely on `θ mod 72°` for the E/T classification, which depends on the
  exact phase the producer fixes. **Coherent (expected).** Action: keep the
  sign and the explicit `-z_j`; only rename `cs`/`sn`/`Qcos`/`Qsin` and drop
  the dated parenthetical (E8/E9).
- **Note (consumer-side comment bug, out of this file's scope):**
  `h5-reader/src/model/QtSpecialBuffers.h:375` comments `pucker_theta // rad`
  but the value is **degrees** (producer + catalog + H5 attr all say degrees).
  This is a stale comment on the consumer side; flag for the h5-reader owner,
  not fixed here (brief restricts edits to PlanarGeometryResult.{h,cpp}).

### omega / omega_deviation and the WrapPi bit-identity (C24, CL11)
- **Reason discovered.** `omega_deviation_ = WrapPi(omega - π)` with
  `WrapPi` = `std::remainder(a, 2π)` → range `[-π, π]`. The same `omega` and
  the same `WrapPi` are independently re-derived in
  `DihedralTimeSeriesTrajectoryResult` using a byte-for-byte identical
  `Dihedral` (atan2 formulation) and `WrapPi` (`std::remainder`). The catalog
  header documents a cross-result test asserting the two producers agree
  **bit-identically** at the wrap boundary. That is why both helpers must keep
  the exact same form — this is the load-bearing reason to **decline** renaming
  `Dihedral`/`WrapPi` (C15/C18) and to **keep** the bit-identity sentence in
  the comment (E4).
- **Consumers.** `omega_actual.npy`, `omega_deviation.npy`, `omega_is_xpro.npy`
  (catalog, per-residue). H5 dihedral_time_series group (DihedralTS, separate
  producer). `h5-reader/QtPerResidueBuffers.h` (omega/omega_deviation/
  omega_is_xpro). `omega_is_xpro` mask convention (set on the i-row carrying
  the peptide bond into i+1=Pro) is explicitly mirrored by DihedralTS
  (cpp:198-205 "we follow that convention").
- **Producer/consumer agreement.** Agree, and additionally cross-checked
  producer-vs-producer by the bit-identity test. **Coherent (expected).**
  Minor doc nit: catalog says `omega_deviation` is `(-π, π]` while
  `std::remainder` yields `[-π, π]`; both endpoints are measure-zero and
  consistent with "wrapped to ±π" — not worth a number-touching change, noted
  in §6.

### chi[1] = χ₂ (CL18)
- **Reason discovered.** `Residue::chi[4]` is 0-indexed; `chi[0]`=χ₁,
  `chi[1]`=χ₂. Verified against every other consumer that reads `chi[k]`
  (DihedralTimeSeriesTrajectoryResult, DsspResult, ChiRotamerSelectionTR,
  DihedralBinTransitionTR) — all 0-indexed. So `res.chi[1]` at cpp:354 is
  genuinely χ₂, matching the `aromatic_chi2_` field name and the Akke &
  Weininger 2023 anchor. **Coherent (expected).** No change.

### Bond field width (CL16)
- **Reason discovered.** `Bond::atom_index_a`/`atom_index_b` are `size_t`
  (src/Bond.h:15-16); the loop variable `ai` is `size_t` (Compute). The
  `b.atom_index_a == ai` comparison at cpp:184 is `size_t == size_t` — no
  signed/width mismatch. **Non-issue.** No change.

### pucker_valid counts finite-Q-but-NaN-θ (C32)
- **Reason discovered.** Under the sub-amplitude guard, `CremerPople5Ring`
  returns `{Q, NaN}` — a *finite, meaningful* amplitude with an undefined
  phase. `pucker_valid` increments on `!isnan(cp.Q)`, so it counts rings with
  a finite Q regardless of θ. This is correct as a "how many rings produced a
  pucker amplitude" diagnostic, but the printed label `"pucker valid="`
  overstates it (a reader could think θ is also valid). It is a **log line
  only** — it is not serialized and feeds no output. **Action:** relabel the
  diagnostic to `"pucker finite-Q="` (E16); do not change the count. This
  moves no number and keeps the (intentional) finite-Q semantics.

---

## 5. Bug-by-exhaustion candidates

**None.** No sign or value finding survived tracing:

- Pyramidalization sign — traced to every consumer (ui display, RestServer,
  learn feature, no TR carrier); all treat it as a stable signed magnitude;
  none requires an absolute chemical orientation. Coherent.
- Cremer-Pople θ — traced producer → RingPuckerTS (verbatim copy) → H5 attr
  (θ mod 72 classification) → h5-reader/learn; the negative-sine + sin/cos-
  basis normal is the documented, iterated fix. Coherent.
- omega_deviation — producer-vs-producer bit-identity test plus consumer trace.
  Coherent.
- `pucker_valid` — a diagnostic label, not a value; relabel only (E16), no
  number moves.

`<array>`/`<algorithm>` are missing direct includes (currently transitive).
This is a latent self-containment defect, not an algorithmic bug; E15 fixes it
without touching behavior.

---

## 6. Questions & Ambiguities

1. **omega_deviation range wording.** Code (`std::remainder`) yields `[-π, π]`;
   the SDK catalog string says `(-π, π]`. Both are fine in practice (boundary
   is measure-zero). Should the catalog string be normalized to `[-π, π]` for
   accuracy? That edit is in `_catalog.py` (a contract description, not a key)
   and outside this file — flag for the SDK owner, not planned here.
2. **h5-reader `pucker_theta` "rad" comment.** `QtSpecialBuffers.h:375` labels
   pucker_theta as radians; it is degrees everywhere on the producer side. A
   consumer-side comment bug. Out of scope for this file; should be raised with
   the h5-reader owner. Confirm whether to file it.
3. **Saturated-ring non-5-ring silent skip (CL19).** A future 6-membered
   saturated ring would be left NaN with no log line (cpp:379). Both reviewers
   call this acceptable at current scope (Pro pyrrolidine only). Confirm it
   stays out of scope, or whether a one-line `OperationLog` skip note is wanted
   — that would be a behavior/log addition, deferred pending the user's call.
4. **`res_i` rename.** E12b renames `res_ip`→`res_next`; leaving `res_i` as-is
   vs renaming to `res_curr` is a taste call. Proposed minimal: rename only
   `res_ip`. Confirm if `res_curr` is also wanted.

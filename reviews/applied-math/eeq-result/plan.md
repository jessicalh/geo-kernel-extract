# Fix plan — EeqResult.{h,cpp}

## 1. Summary

`EeqResult.h` and `EeqResult.cpp` tell a coherent story. The header is a
self-contained EEQ spec (energy functional → constrained minimisation →
parameter provenance → kernel choice → output contract). The `.cpp` is
explicitly staged in five numbered banners (CN → χ_eff → Coulomb matrix →
constrained solve → store) that trace the physics in textbook order, with
dense formulas pinned to citations at the point of use. Both reviews agree
the file is "close to the good case."

The friction is almost entirely **local naming inside the solve block**
(`u`, `v`, `q`, `gam_prod_inv`) and **two comments that overstate the math**
they label ("Enforce exact charge neutrality" before a later clamp; the
unused `n_fallback` with a contradictory comment). Those are exactly the
improvable internal-name / comment-conformance kind.

This pass will: improve solve-block and matrix local names, add 2–4 word
field-aliases / signposts where a textbook symbol does not tie to a struct
field, and correct the two overstating comments to state what the code
actually does. It will **not** touch: the algorithm, any number, the output
NPY/H5/SDK names (`eeq_charges`, `eeq_cn`, `eeq_welford` and its datasets),
or the recorded `GeometryChoice` key strings (those are serialized
provenance contract).

**Key traced finding (controls several review items):** charge neutrality
is enforced *before* the per-atom clamp, and the clamp can move the stored
sum off `Q_total`. Every consumer reads the **stored, post-clamp**
`eeq_charge` and none asserts or re-derives neutrality — so this is
**coherent (neutrality is a pre-clamp property; clamp is a rarely-fired
guard)**, and the only fix is to make the comment say so. See §4.

---

## 2. Review-finding ledger

Both reviews are readability passes. `codex-correctness.md` is referenced in
`codex.md` but is **not present** in this directory (only `codex.md` and
`claude.md` exist); its findings are therefore not separately ledgered. The
correctness items folded into `codex.md` §6 are ledgered below.

### codex.md

| # | Finding (loc) | Disposition |
|---|---|---|
| C1 | h:23 — output story lists only `eeq_charges.npy`; impl also writes `eeq_cn.npy` | **adopted** → §3 E1 |
| C2 | cpp:14 — large "D4 EEQ parameters" banner then 26–27 say params live elsewhere; collapse | **adopted (partial)** → §3 E2 (collapse the banner; do not move prose to PhysicalConstants.h — out of scope file) |
| C3 | cpp:51 — `n_fallback` computed, never used | **adopted** → §3 E3 (surface in the parameter `GeometryChoice` record; see also CL10) |
| C4 | cpp:153 — solve stage needs named intermediates (`A_inv_chi`, `A_inv_ones`, `constraint_scale`, `charges`) | **adopted** → §3 E5 (names chosen to match the banner's algebra; see §4) |
| C5 | cpp:187 — "Enforce exact charge neutrality" then later clamp can change sum; reword | **adopted** → §3 E6 + §4 (the central traced item) |
| C6 | cpp:32 — `N` → `atom_count` | **declined** → `N` is the matrix dimension used in `Eigen::MatrixXd(N,N)`, all banners, and the log line; `N` is the textbook symbol here and renaming it fragments the math-to-code tie. Conservative good practice keeps `N` for the linear-algebra dimension. |
| C7 | cpp:40 — `Q_total` → `total_charge_e` | **declined** → `Q_total` is the exact symbol in the header functional (`Σqᵢ = Q_total`) and the solve banner RHS; the tie is the point. |
| C8 | cpp:41 — `cn_k` → `cn_steepness` | **adopted** → §3 E4 (matches the TOML key `eeq_cn_steepness` and the Step-1 banner's `k`; low cost, improves the tie) |
| C9 | cpp:86 — `cn` → `coordination_number` | **declined** → `cn(i)` is the `Eigen::VectorXd` named exactly as the banner's `CN_i`; the field is even called `eeq_cn`. Verbose rename here loses the symbol tie. |
| C10 | cpp:95 — `R2` → `distance_bohr_sq` | **adopted** → §3 E7 (units/frame genuinely ambiguous; also unifies with the Coulomb block per C13) |
| C11 | cpp:135 — `A` → `coulomb_hardness_matrix` | **declined** → `A` is the augmented-system matrix named `A` in the Step-4 banner and the `LLT<>(A)` call; renaming breaks the banner tie. A signpost already exists (Step-3 banner). |
| C12 | cpp:174 — `u`,`v` are the worst names → `A_inv_chi_eff`, `A_inv_ones` | **adopted** → §3 E5 |
| C13 | cpp:184 `lambda` → `charge_constraint_lambda`; cpp:185 `q` → `charges` | **adopted (partial)** → §3 E5 renames `q`→`charges` only at storage-facing use; **declines** the `lambda` rename — `λ` is the Lagrange multiplier named `λ` in the banner and is used only on the adjacent two lines. |
| C14 | cpp:100/146 — duplicated distance calc; name the three quantities consistently | **adopted** → §3 E7 (naming only; no abstraction, as the reviewer themselves scoped) |
| C15 | cpp:146 — `gam_prod_inv` fuses semantics → `hardness_product_inverse` | **adopted** → §3 E8 |
| C16 | cpp:184 — split constraint into numerator/denominator/lambda/charges | **declined** → over-fragmentation. `denom` already names the denominator (line 178); the numerator is one term. Splitting two clean lines into four named temporaries is the "condescending fragmentation" the brief warns against. Renaming `q`→`charges` (E5) is the useful part and is adopted. |
| C17 | cpp:225 — stats seed from `q(0)` but report clamped `conf` charges; seed from stored or rename | **adopted** → §3 E9 (seed from stored charge; see §4 and CL9 — both reviews agree) |
| C18 | h:48 — `Compute` → `ComputeEeqCharges` | **declined** → `Compute` is the `ConformationResult` result-pattern convention used by every sibling calculator (`*Result::Compute`); renaming one breaks the pattern and the OperationRunner dispatch expectation. Cross-file convention, weighed and declined. |
| C19 | h:50 — `WriteFeatures` framework terminology, acceptable | **declined (no-op)** → reviewer concurs it's acceptable; it's an overridden virtual. |
| C20 | h:42 — `Name()` clear enough | **declined (no-op)** → no change requested. |
| C21 | cpp:38 — "Named constants from TOML (no buried literals)" is process prose → `// configuration` | **adopted** → §3 E10 |
| C22 | cpp:48 — signpost slightly long → `// atom parameters` | **declined** → current "Pre-cache per-atom parameters and positions in Bohr" carries the load-bearing "in Bohr" unit fact; trimming to "atom parameters" drops it. Keep. |
| C23 | cpp:69 — "GeometryChoice: parameter summary" clear | **declined (no-op)** |
| C24 | cpp:80 — block label could be terser | **declined** → "(error function counting)" names the CN method (Caldeweyher); it earns its length. |
| C25 | cpp:109 — good label, no change | **declined (no-op)** |
| C26 | cpp:128 — Ohno-Klopman justification reads as process prose → terse `// Ohno-Klopman approximation` | **adopted (partial / compress only)** → §3 E11. **Do not delete** the rationale — it documents a deliberate divergence from the cited reference that *moves the numbers* (claude CL8 agrees it is justified). Compress the last two sentences; keep the divergence statement. |
| C27 | cpp:137 — replace self-Coulomb diagonal block with `// self-Coulomb diagonal` | **adopted (partial)** → §3 E12 (keep the formula + citation, trim the "Makes the matrix..." explanatory sentence into the SPD note that already exists at the solve banner) |
| C28 | cpp:163 — SPD explanation "too confident" → `// Cholesky solve` | **declined** → the SPD justification is load-bearing (it explains *why* LLT is valid here and why the diagonal carries the `√(2/π)/r` term). The brief says comments conform to code, not that confident-but-correct comments be stripped. It is correct (diagonal dominance from the self-energy term); keep. |
| C29 | cpp:187 — too much prose, misleading after clamp → `// solve residual correction` | **adopted** → §3 E6 (same item as C5; reword to state it is the *pre-clamp* solve-residual removal) |
| C30 | cpp:210 — "traceable decision" wordy → `// charge clamp record` | **adopted** → §3 E13 |
| C31 | cpp:191/205 — post-clamp neutrality may break; check path or fix comment | **adopted (comment only)** → §4 + §3 E6. Traced: coherent, comment fix only. **Duplicate of C5/C29.** |
| C32 | cpp:225 — `q(0)` no empty-conformation guard; crash if `AtomCount()==0` | **→ §6 Question Q1** (need the upstream invariant; both reviews raise it — duplicate of CL11) |

### claude.md

| # | Finding (loc) | Disposition |
|---|---|---|
| CL1 | cpp:140 — diagonal comment uses `η_i`,`r_i` but code reads `params[i].gam`,`params[i].rad`; no in-line tie | **adopted** → §3 E14 (one-line field-alias / tie comment at first use) |
| CL2 | cpp:146-147 — confirm `gam` on off-diagonal is same η as diagonal (it is); good | **declined (no-op / confirmation)** → traced true; the E14 alias makes it visible. |
| CL3 | cpp:184-185 — sign chain matches banner exactly; reads cleanly | **declined (no-op)** → confirms coherence; recorded in §4. |
| CL4 | cpp:140,146 — `gam` invites confusion with the `gamma`/γ kernel 3 lines down; add local `eta_i` alias | **adopted** → §3 E14 (the alias also fixes CL1; resolves the γ-vs-gam clash) |
| CL5 | cpp:140 — `params[i].rad` ambiguous (radius/radians/rcov?), `rcov` is a sibling field; alias `r_gauss` | **adopted** → §3 E14 (Gaussian-charge-radius alias at first use) |
| CL6 | cpp:174 — `u`,`v` carry no meaning at use sites; acceptable given banner; minor | **adopted** → §3 E5 (duplicate of C12; the rename removes even the minor friction) |
| CL7 | cpp:40/41 — `Q_total`,`cn_k` fine given banner | **declined (no-op)** → agrees with our C7 decline; we adopt `cn_k`→`cn_steepness` (E4) only because it ties to the TOML key. |
| CL8 | cpp:128-133 — Ohno-Klopman-vs-Gaussian rationale justified (documents a number-moving divergence); keep, compress last two sentences | **adopted** → §3 E11 (duplicate of C26; claude's "keep but compress" is the disposition we take) |
| CL9 | cpp:225-231 — `q_min=q(0)` seeds from unclamped `q` while loop reads clamped stored charge; seed from stored | **adopted** → §3 E9 (duplicate of C17; both agree) |
| CL10 | cpp:55-56 — comment "S uses real params; Unknown uses fallback" but the `if` counts BOTH S and Unknown into `n_fallback`; comment contradicts code | **adopted** → §3 E3 (the comment is wrong about S; see §4 — relabel `n_fallback`→`n_special` and fix the comment) |
| CL11 | cpp:191-192 — neutrality shift divides by `N`; div-by-zero if `N==0`; no empty guard | **→ §6 Question Q1** (duplicate of C32) |
| CL12 | cpp:137 — verify "Caldeweyher 2019 Eq. 8" citation is the self-energy equation | **→ §6 Question Q2** (citation-number verification needs the paper; cannot confirm from the source tree) |
| CL13 | cpp:128-147 — hybrid Ohno-Klopman off-diagonal + Gaussian self-energy diagonal is a documented deliberate hybrid, not a bug; flag so nobody "fixes" the diagonal | **adopted** → §3 E11/E12 keep the divergence note explicit; §4 records it as coherent. |
| CL14 | cpp:117 — `1e-14` guard signpost is exactly right | **declined (no-op)** → praised; no change. |
| CL15 | cpp:187-190 — neutrality-shift comment clear and grounded; good | **partially superseded** → the comment is well-written *for the solve residual* but C5/C31 show it overstates "exact neutrality" relative to the later clamp. E6 keeps its grounded phrasing and adds the pre-clamp scope. |
| CL16 | cpp:184 — `q=-(u+λv)` sign self-consistent; can't independently confirm EEQ sign of χ from this file; check `D4EeqParams.chi` | **→ §6 Question Q3** (sign of `chi` in the parameter table; see §4 for what was traced) |
| CL17 | cpp:55 — `n_fallback` accumulated but never read/reported anywhere | **adopted** → §3 E3 (duplicate of C3; surface in the parameter record) |

---

## 3. Edits that don't move numbers

All line numbers are pre-edit (current `HEAD`).

- **E1** — `EeqResult.h:24` — add `eeq_cn.npy (N,) float64 — coordination
  number (intermediate, for traceability)` under the `Output:` block so the
  header matches `WriteFeatures` (which writes two NPYs). *(C1)*

- **E2** — `EeqResult.cpp:14-24` — collapse the 11-line "D4 EEQ parameters"
  banner to a 2–3 line note (source + DOI + "values in PhysicalConstants.h
  (`D4EeqParams`, `D4EeqParamsFor`)"). The line-26-27 note already says where
  the params live; remove the redundancy, keep the citation. Do not relocate
  prose into `PhysicalConstants.h` (different file, out of scope). *(C2)*

- **E3** — `EeqResult.cpp:51,55-57` — rename `n_fallback` →
  `n_special_params` (or `n_fallback_or_special`) and fix the comment to
  match: the counter increments for **both** `S` (real params, flagged as a
  watched element) and `Unknown` (true fallback). Surface the value in the
  existing `eeq_parameters` `GeometryChoice` record as
  `AddNumber(gc, "n_special_params", ...)` so it is no longer dead. *(C3,
  CL10, CL17)* — see §4 for why S is counted.

- **E4** — `EeqResult.cpp:41` (+ the lambda capture at :71-73 and the record
  key) — rename local `cn_k` → `cn_steepness` to match the TOML key
  `eeq_cn_steepness` and the recorded `GeometryChoice` field name
  `cn_steepness`. The recorded **key string** stays `"cn_steepness"` (already
  matches). *(C8, CL7)*

- **E5** — `EeqResult.cpp:174-185` — rename the solve-block locals to match
  the Step-4 banner's algebra:
  - `u` → `A_inv_chi` (solves `A·u = χ_eff`)
  - `v` → `A_inv_ones` (solves `A·v = 1`)
  - `q` → `charges` (the storage-facing vector)
  - update the banner at :159-161 to use the same names so banner and code
    read identically, OR keep the banner's `u`/`v`/`q` math symbols and add a
    one-line "`u ≡ A_inv_chi, v ≡ A_inv_ones`" tie. **Prefer renaming the code
    to the descriptive names and keeping the banner's math symbols** — the
    banner is the derivation, the code is the implementation; the tie line
    bridges them. Keep `lambda`/`λ` and `denom` as-is (adjacent, banner-named).
    *(C4, C12, C13-partial, CL6)*

- **E6** — `EeqResult.cpp:187-192` — reword the "Enforce exact charge
  neutrality" comment to state precisely what it does: it removes the
  **Cholesky solve residual** (proportional to cond(A)) so the *pre-clamp*
  sum equals `Q_total`; the per-atom clamp in Step 5 may move the *stored*
  sum off `Q_total` when it fires. Keep the existing grounded "residual
  proportional to cond(A) / uniform shift preserves differences" content.
  Suggested: `// Remove the Cholesky solve residual so the (pre-clamp) sum`
  `// equals Q_total. Uniform shift preserves per-atom differences;`
  `// residual is ~cond(A). NOTE: the Step-5 clamp can move the stored sum`
  `// off Q_total when it fires (see n_clamped).` *(C5, C29, C31, CL15)* —
  see §4.

- **E7** — `EeqResult.cpp:95,100,101` and `:145` — name the distance
  quantities consistently across the CN block and the Coulomb block:
  `R2` → `dist_bohr_sq`, `R` → `dist_bohr`, `rcov_sum` → `rcov_sum_bohr`.
  Apply the same `dist_bohr_sq` name at :145. Naming only, no shared helper.
  *(C10, C14)*

- **E8** — `EeqResult.cpp:146` — rename `gam_prod_inv` →
  `hardness_product_inv` (it is `1/(η_i·η_j)`). *(C15)*

- **E9** — `EeqResult.cpp:225` — seed `q_min`/`q_max` from the **stored**
  charge `conf.AtomAt(0).eeq_charge` instead of the unclamped `q(0)`, so the
  reported min/max are consistent with the values the loop body compares and
  with what every consumer reads. **This can change the reported
  `charge_min`/`charge_max` statistic only in the case where atom 0 is itself
  clamped and is the extreme** — i.e. it corrects a stats-reporting
  inconsistency, it does not change any stored charge or the algorithm. Flag
  to the human in §5 as the one "may move a *recorded number*" edit; the
  number it moves is a provenance statistic, not a physics output. *(C17,
  CL9)*

- **E10** — `EeqResult.cpp:38` — replace `// ── Named constants from TOML
  (no buried literals) ──` with `// ── Configuration (TOML) ──`. The
  "no buried literals" is process commentary. *(C21)*

- **E11** — `EeqResult.cpp:128-133` — compress the Ohno-Klopman rationale:
  keep the first sentences that state the deliberate divergence from the
  dftd4 Gaussian form and that **charges differ from exact dftd4**; compress
  the closing "For our purpose ... acceptable" to a single clause
  (e.g. `// geometry-responsive charges, not exact dftd4 reproduction`). Do
  **not** delete the divergence statement — it documents a number-moving
  choice. *(C26-partial, CL8, CL13)*

- **E12** — `EeqResult.cpp:137-139` — keep the diagonal formula + Eq.
  citation; trim the standalone "Makes the matrix diagonally dominant → SPD"
  sentence (it is restated at the solve banner :163-165). *(C27-partial,
  CL13)*

- **E13** — `EeqResult.cpp:210` — `// GeometryChoice: traceable decision —
  charge clamped` → `// GeometryChoice: charge-clamp record`. The recorded
  **key string** `"eeq charge clamp"` is provenance contract — leave it.
  *(C30)*

- **E14** — `EeqResult.cpp:140` — add a one-line field-alias tie at first
  use of the diagonal so the textbook symbols resolve to struct fields, e.g.
  `// gam = chemical hardness η_i; rad = Gaussian charge radius r_i`
  (`rcov`, used in Step 1, is the separate covalent radius). Optionally
  introduce local aliases `const double eta_i = params[i].gam;` /
  `const double r_gauss_i = params[i].rad;` at the diagonal to kill the
  `gam`-vs-`gamma` clash three lines down — **suggestion, weigh against the
  extra locals**; the comment tie alone may suffice. The `D4EeqParams`
  field names (`gam`, `rad`, `rcov`, `chi`, `kappa`) live in
  `PhysicalConstants.h` and are **not renamed** (cross-file struct contract;
  the comment/alias is the local remedy). *(CL1, CL4, CL5)*

**Not touched:** output names `eeq_charges`/`eeq_cn`/`eeq_welford` and its
datasets; `GeometryChoice` record key strings; `D4EeqParams` field names;
`Compute`/`WriteFeatures`/`Name`/`Dependencies` signatures; every numeric
literal and formula.

---

## 4. Usage notes (the traced reasons)

### 4a. Post-solve clamp vs enforced neutrality `Σq = Q_total` — **coherent (expected)**

**What the producer does, in order:** (1) solve the augmented system; (2)
`q.array() -= q_residual/N` at :191-192, which makes `Σq == Q_total`
*exactly* at that point; (3) in the Step-5 store loop (:200-221) each `qi` is
clamped to `±eeq_charge_clamp` (default **2.0 e**) *before* being written to
`atom.eeq_charge`. The clamp therefore happens **after** the neutrality shift,
so when any atom is clamped the **stored** sum no longer equals `Q_total`.

**Why this is coherent and not a bug — traced across every consumer:**

- `EeqWelfordTrajectoryResult::Compute` (`src/EeqWelfordTrajectoryResult.cpp:66`)
  reads `conf.AtomAt(i).eeq_charge` — the stored, post-clamp value — and runs
  per-atom Welford. It never sums charges, never checks neutrality.
- `learn/bones/GromacsProtein.cpp:282-285` reads stored `eeq_charge` /
  `eeq_cn` per atom (with a `==0.0` "did EEQ run" guard). No neutrality check.
- `learn/bones/AnalysisWriter.cpp:229,998` writes the stored per-atom
  `eeq_charge` straight to H5 `(T,N)`. No sum.
- `ui/src/RestServer.cpp:1023` and `ui/src/MainWindow.cpp:1385` display the
  stored per-atom value. Per-atom only.
- SDK `python/nmr_extract/_catalog.py:250` (`eeq_charges`, units `e`) and
  `_protein.py:1479` expose the stored `(N,)` array as-is. No neutrality
  invariant in the SDK `load()` path for this array.

**Conclusion:** producer and all consumers agree on one thing —
`eeq_charge` is the **stored per-atom charge**, consumed per-atom. *No
consumer requires the stored sum to equal `Q_total`.* The clamp is a guard
against pathologically large charges from an ill-conditioned solve; it fires
rarely (clamp at 2 e), is counted in `n_clamped`, recorded as a
`GeometryChoice` per clamped atom, and logged. Neutrality is therefore a
**pre-clamp solve property**, which is what the algorithm guarantees and what
the comment *should* say. The fix is **comment-only** (E6): state that the
shift removes the solve residual to make the pre-clamp sum exact, and that the
clamp can move the stored sum. The bug, if one insisted on exact stored
neutrality, would be a *spec* question, not a producer or consumer defect —
no side currently asserts that spec. **Coherent; comment fix only.**

### 4b. `q_min`/`q_max` seeded from unclamped `q(0)` — **bug-by-exhaustion, producer-side, statistic-only** (see §5)

Both reviews caught this and they are right; it is small and bounded. Treated
in §5 with the exhaustion.

### 4c. Sign chain `λ = -(Q + 1'u)/denom`, `q = -(u + λv)` — **coherent**

Traced against the Step-4 banner (:155-161): the augmented system is
`[A 1; 1' 0][q; λ] = [-χ_eff; Q]`. Block elimination of that system gives
exactly `u = A⁻¹χ_eff`, `v = A⁻¹1`, `λ = -(Q + 1'u)/(1'v)`,
`q = -(u + λv)`. Code and banner agree line-for-line (claude CL3, CL16
both confirm self-consistency). The **only** thing not confirmable from this
file is the sign convention of `χ` (`D4EeqParams.chi`) in the parameter table
— deferred to Q3.

### 4d. `n_fallback` counter — **coherent intent, contradictory comment**

The counter increments for `Element::Unknown` (genuine fallback to default
params) **and** `Element::S` (sulfur). Sulfur uses *real* D4 params, so the
comment "S uses real params; Unknown uses fallback" describes the elements
but mislabels the counter as "fallback" while deliberately including S. The
plausible reason S is watched: sulfur is a known soft spot for charge-equil
parameterisations (large, polarisable; common clamp culprit), so the run
wants a count of "atoms whose EEQ params merit watching." The fix (E3) is to
rename the counter to reflect "special/watched," fix the comment, and surface
the value (currently dead) in the parameter `GeometryChoice` record.

---

## 5. Bug-by-exhaustion candidates

**B1 — `q_min`/`q_max` statistic seeded from the pre-clamp `q(0)` (producer-side, statistic-only).** `src/EeqResult.cpp:225`.

- **Exhaustion.** The stats loop (:225-231) seeds `q_min = q_max = q(0)`
  (the *unclamped* Eigen vector) then iterates comparing against
  `conf.AtomAt(i).eeq_charge` (the *clamped, stored* value). For every atom
  *except index 0* the seed and the compared values are the same domain
  (stored). Only if atom 0 is itself clamped **and** is the running extreme
  can the reported `charge_min`/`charge_max` reflect the pre-clamp value of
  atom 0 — inconsistent with the stored array every consumer reads (traced in
  §4a: Welford, GromacsProtein, AnalysisWriter, RestServer, MainWindow, SDK
  all read stored). The recorded `charge_sum` (:228) already uses the stored
  value, so within the same record `sum` and `min/max` can disagree about
  what domain they describe.
- **Verdict.** Producer-side inconsistency in a **provenance statistic only**
  (the `eeq_charge_statistics` `GeometryChoice` record and the Info log). It
  does **not** touch any stored charge, NPY, or H5 per-atom value, and no
  consumer reads these min/max stats for physics.
- **Minimal fix (E9).** Seed both from the stored charge:
  `double q_min = conf.AtomAt(0).eeq_charge, q_max = q_min;`. This *may change
  the recorded `charge_min`/`charge_max` number* in the atom-0-clamped-extreme
  case — flagged because the brief reserves number-moving for outright bugs
  after exhaustion. This qualifies: it is a reporting bug, the fix is minimal,
  and it moves only a recorded statistic, never a physics output. Both reviews
  independently flagged it (C17, CL9).
  **APPROVED 2026-05-24 (Jessica): fix it.**

No other number-moving change is proposed. The clamp-vs-neutrality item (§4a)
is comment-only.

---

## 6. Questions & Ambiguities

- **Q1 — Empty-conformation guard (`N == 0`).** `src/EeqResult.cpp:191-192`
  divides by `N`, and :225 reads `q(0)` unguarded. Both reviews flag a
  potential div-by-zero / out-of-bounds if `AtomCount() == 0`. Is a non-empty
  conformation an upstream invariant (OperationRunner only attaches results
  to conformations with atoms)? If guaranteed, document the assumption near
  :32 (one line). If not, add an early `if (N == 0) return ...;`. **Need the
  upstream invariant before choosing** — this is a behaviour question, not a
  readability edit. *(C32, CL11)*

- **Q2 — Citation "Caldeweyher 2019 Eq. 8" for the Gaussian self-energy
  diagonal** (`:137`). Cannot verify the equation number from the source tree;
  the paper PDF would confirm. The comment is load-bearing (it justifies the
  hybrid diagonal). Leave as-is unless the human can confirm/correct the Eq.
  number. *(CL12)*

- **Q3 — Sign of `D4EeqParams.chi`.** The solve sign chain is internally
  self-consistent with the banner (§4c), but the EEQ convention for the sign
  of electronegativity entering `χ_eff` cannot be confirmed from `EeqResult`
  alone; it depends on how `chi` is stored in `PhysicalConstants.h`
  (`D4EeqParams`) relative to the dftd4 reference. No edit proposed; flagging
  per the governing prior that sign questions are resolved by tracing the
  parameter source, which is a separate (frozen) file. *(CL16)*

- **Q4 — E5 banner-vs-code naming preference.** §3 E5 recommends descriptive
  code locals (`A_inv_chi`, `A_inv_ones`, `charges`) with the math banner
  keeping `u`/`v`/`q` plus a one-line tie. The alternative is to rename the
  banner symbols too. Confirm the preferred convention (derivation-symbols in
  the banner, descriptive-names in code) before applying.

- **Q5 — `codex-correctness.md` absent.** `codex.md` references a sibling
  `codex-correctness.md` (earlier correctness-first pass) that is not in this
  directory. If it exists elsewhere, its findings should be ledgered too; this
  plan ledgers only the two reviews present.

# Fix plan — AIMNet2ChargeResponseGradientResult

Targets: `src/AIMNet2ChargeResponseGradientResult.h`, `src/AIMNet2ChargeResponseGradientResult.cpp`
Scope: readability only. No algorithm change, no number change, no refactor.

---

## 1. Summary

The file tells a coherent story end to end. The header states the objective
`L = sum_j q_j^2`, explains *why* the naive `sum(q)` objective is degenerate
(AIMNet2 enforces total-charge conservation → `d(sum q)/d(coord) ≈ 0`), and
labels the L2-of-charges scalar as a deliberate **single-pass proxy** — the
exact `d(q_i)/d(r_i)` diagonal would need N backward passes. **This is the
earned design reason, documented in three independent places** (the header
block, the inline objective comment at `.cpp:149-158`, and the SDK at
`_catalog.py:372-377` / `_tensors.py:807-818`). It is correct as a proxy and
must not be flagged as wrong. The `.cpp` walks: guard → record objective →
build inputs (duplicating AIMNet2Result's tensor build) → forward with grad →
slice sentinel → L2 objective → backward → finite pre-scan → store per-atom
vector+norm + log summary.

The friction is entirely cosmetic and falls in two buckets:
1. **Provenance prose embedded in source** — dated amendment tags, "PROMOTED
   FROM TEST FLAG", "Codex 2026-05-20 F1", "first-test choice". This is
   changelog/history that belongs in git/memory, not in an interface header or
   a local guard comment. Both reviews flag this consistently.
2. **Terse tensor names** — `N`, `N1`, `numbers_cpu`, `charge_cpu`,
   `charges_n`, `loss`, `v`, `s`. All local; improving them is a goal, not
   churn.

What the fix pass **will** touch: local variable names, comment wording,
2–4-word signposts before existing stages. What it **will not** touch: the
algorithm, the signs, the numbers, the output NPY names, the
`ConformationAtom` field names (`aimnet2_charge_response_gradient_vector` /
`_scalar`), the GeometryChoice record name, or the class name (see ledger).

**One naming inconsistency worth surfacing** (not a code bug): the C++ source
comments still call the quantity "polarisability" (`.h:3`, `.h:15`,
`.cpp:44 "polarisability backward"`), but the SDK was deliberately **renamed
away from "polarisability" on 2026-05-20 (commit 58594f5)** because the
emission is `∂(Σq²)/∂r`, NOT a Buckingham `α = ∂μ/∂E` (`_catalog.py:376-377`,
`_tensors.py:813-814`). The source comments are now stale relative to the
contract the SDK advertises. Aligning the comment language is a pure-comment
fix (§3). The `ui/` viewer still labels it "AIMNet2 polarisability"
(`MainWindow.cpp:1487, 1490`) — that is a `ui/`-owned file, out of scope for
this pass, but noted in §4 as downstream carry-through the human may want.

---

## 2. Review-finding ledger

### codex.md

| # | Finding (loc) | Disposition |
|---|---|---|
| C1 | `.h:19` — result name sounds like the exact diagonal `d(q_i)/d(r_i)`; add a "returned quantity" sentence | **adopted** → §3 E1 (the header already has the rationale; tighten to lead with the returned quantity). Class rename **declined** — see C13. |
| C2 | `.h:24` — "Per Amendment 2026-05-08(b)…" is audit-trail prose | **adopted** → §3 E2 (collapse to a terse signpost) |
| C3 | `.h:28` — "PROMOTED FROM TEST FLAG…" reads like release history | **adopted** → §3 E3 |
| C4 | `.cpp:66` — `objective_kind = 1.0` not self-explanatory; signpost it | **adopted** → §3 E4 (the `"L2_of_charges"` label string already disambiguates; add a 3-word comment) |
| C5 | `.cpp:149` — core math comment repeats header + "first-test choice" | **adopted** → §3 E5 (trim history, keep equation) |
| C6 | `.cpp:175` — non-finite pre-scan comment over-explains pipeline history ("Codex 2026-05-20 F1") | **adopted** → §3 E6 (keep the Welford-poisoning rationale, drop the dated tag) |
| C7 | `.cpp:40` — `N` weak for outside reader → `atom_count` | **declined** → §3 note N1 (see reasoning; `N` is used 9× as a loop bound and pairs with `N1`; renaming both is the better move but cost-weighed — proposed as optional E11) |
| C8 | `.cpp:80` — `N1` hides sentinel convention → `padded_atom_count` | **adopted** → §3 E7 |
| C9 | `.cpp:82` — `coord_cpu` no units/frame → `coords_angstrom_cpu` | **adopted (comment form)** → §3 E8. Rename **declined**; unit added as comment instead (positions come from `conf.PositionAt`, Å is the library convention; a comment is lower-risk than baking a unit assertion into the name). |
| C10 | `.cpp:92` — `numbers_cpu` → `atomic_numbers_cpu` | **adopted** → §3 E9 (also `num_acc` → `atomic_num_acc`) |
| C11 | `.cpp:105` — `charge_cpu` is total charge, not per-atom → `total_charge_cpu` | **adopted** → §3 E10 |
| C12 | `.cpp:159` — `charges_n` cryptic → `real_atom_charges` | **adopted** → §3 E12 |
| C13 | `.cpp:160` — `loss` generic → `charge_l2_objective` / `charge_norm_squared` | **adopted** → §3 E13 |
| C14 | `.cpp:196` — `v`/`s` → `charge_response_grad` / `grad_norm` | **adopted (partial)** → §3 E14. Claude (CL4) says `v`/`s` are fine as-is (local, obvious). Compromise: rename to `grad_vec` / `grad_norm` — short, unambiguous, doesn't bloat the storage lines. |
| C15 | `.cpp:149` (visible math) — make "drop sentinel charges" / "form L2 objective" visible stages | **adopted** → §3 E5/E12/E13 (the renames + trimmed comment achieve this without splitting code) |
| C16 | `.cpp:172` (stage labels) — add `// copy gradient`, `// finite-gradient guard`, `// store gradients`, `// norm summary` | **adopted (partial)** → §3 E6, E15. The finite-guard and the store/summary loop get one-line signposts; "copy gradient" on `.cpp:172` is self-evident from the cast line — declined to avoid over-labelling. |
| C17 | `.cpp:193` — storage + summary in one loop; add a signpost | **adopted** → §3 E15 (duplicate of CL5) |
| C18 | `.h:58` — class name underspecified; `AIMNet2ChargeL2CoordinateGradientResult` more literal, else state returned quantity in comment | **declined (rename) / adopted (comment)** → see §3 note R1. Rename declined: the class name is the H5/TR group lineage string (`AIMNet2ChargeResponseGradientResult.{…}` appears verbatim in two TimeSeries/Welford TR provenance strings, and matches the `aimnet2_charge_response_gradient*` NPY/field/H5-group family). It is effectively a contract surface. The "state the returned quantity" half is adopted via E1. |
| C19 | `.h:68` — `Compute` acceptable | **declined** (no action; concurs) |
| C20 | `.h:72` — `WriteFeatures` acceptable | **declined** (no action; concurs) |
| C21 | `.h:3` — align "charge-polarisation gradient"/"polarisability" with actual `dL/dr` | **adopted** → §3 E16 (the polarisability-language fix; see §1 and §4) |
| C22 | `.cpp:48` — element-guard comment is good | **declined** (no action; concurs) |
| C23 | `.cpp:74` — separator block heavier than needed → `// AIMNet2 inputs` | **adopted (partial)** → §3 E17 (trim the banner; keep the one sentence about gradient-on-coord + shared neighbour helper — that content is load-bearing) |
| C24 | `.cpp:119` — tighten `.detach()` lead label → `// differentiable coordinates` | **declined** → both reviews call this libtorch lore "worth keeping" (CL6); the current 3-line comment explains *why* `.detach()` is needed (leaf-tensor requirement). Tightening the lead loses the reason. Keep as-is. |
| C25 | `.cpp:137` — "Forward pass…" restates next line → `// AIMNet2 forward` | **adopted** → §3 E18 (trim the redundant clause; keep the "gradient tape records the graph" sentence — it is the non-obvious part) |
| C26 | `.cpp:149` — replace long rationale with `// L2 charge objective` + derivative | **duplicate** of C5 |
| C27 | `.cpp:175` — replace with `// finite-gradient guard`; remove dated tag | **duplicate** of C6 |
| C28 | `.cpp:221`/`:236` — WriteFeatures comments good; add units if known | **adopted** → §3 E19 (units are known: `e^2/Å` per `_catalog.py:379,381`; add to both comments) |
| C29 | `.cpp:162` — **correctness check**: do model parameters still `requires_grad`, accumulating unused grads on `loss.backward()`? | **see §6 Q1** — cannot settle from these two files; functionally harmless (only `coord_gpu.grad()` is read), at most a wasted-compute question. Logged as a question, not a fix. |

### claude.md

| # | Finding (loc) | Disposition |
|---|---|---|
| CL1 | `.cpp:105-106` — `charge_cpu`/`mol_idx_cpu` built with no explanation; add `// total system charge = 0` / `// single molecule: all atoms in batch 0` | **adopted** → §3 E10, E20 |
| CL2 | `.cpp:80,90` — `N1 = N+1` and the sentinel comment are 10 lines apart; move a 2-word note to line 80 | **adopted** → §3 E7 (the `padded_atom_count` rename + a `// +1 sentinel padding row` note at the declaration carries this) |
| CL3 | `.cpp:92` — `numbers_cpu`/`num_acc` → `atomic_numbers_cpu` | **duplicate** of C10 → §3 E9 |
| CL4 | `.cpp:159-160` — `charges_n` → `charges_real`/`charges_non_sentinel` | **duplicate** of C12 → §3 E12 |
| CL5 | `.cpp:80` — `N1` → `N_padded` | **duplicate** of C8 → §3 E7 (`padded_atom_count` adopted over `N_padded` for the spelled-out form) |
| CL6 | `.cpp:196-197` — `v`/`s` terse but fine as-is | **noted** — disagrees with C14; resolved in favour of a light rename (`grad_vec`/`grad_norm`) per §3 E14, splitting the difference |
| CL7 | `.cpp:159-162` — slice→square-sum→backward well grouped; good | **declined** (no action; concurs grouping is good) |
| CL8 | `.cpp:193-204` — reduction loop does two things; add `// store per-atom, accumulate norm stats` | **adopted** → §3 E15 (duplicate of C17) |
| CL9 | `Compute`/`WriteFeatures`/`Dependencies` conventional; clear | **declined** (no action; concurs) |
| CL10 | output names `_vector`/`_scalar` self-documenting; good | **declined** (no action; concurs — and these are protected output names regardless) |
| CL11 | `.cpp:149-158` — inline objective comment partly restates header; trim to equation + one line + pointer | **duplicate** of C5 → §3 E5 |
| CL12 | `.cpp:48-49` — element-guard comment exactly right | **duplicate** of C22 (no action) |
| CL13 | `.cpp:175-179` — drop the `Codex 2026-05-20 F1` tag, keep the Welford rationale | **duplicate** of C6 → §3 E6 |
| CL14 | `.cpp:119-121` — `.detach()` lore worth keeping; good | **adopted (as: keep)** — reinforces declining C24 |
| CL15 | `.h:24-34` — Amendment/PROMOTED/RequireConformationResult dates are changelog prose; one line of lifecycle suffices | **duplicate** of C2+C3 → §3 E2, E3 |
| CL16 | `.cpp:159` — **correctness**: `loss` sums over first N, excludes sentinel row; looks right | **confirmed coherent** → §4 (verified: `charges_gpu.slice(0,0,N)` drops row N) |
| CL17 | `.cpp:108,112,116` — **correctness**: squared vs raw cutoff dual use; check `BuildNeighbourMatrix` expects squared, model expects raw | **confirmed coherent** → §4 (traced to AIMNet2Result) |
| CL18 | `.cpp:172` — gradient cast to float64 before finite-scan/norm; good | **confirmed coherent** → §4 (matches `feedback_embedding_float32`: aim stored float32, all derived C++ math float64) |
| CL19 | `.cpp:204` — `mean_norm = sum_norm/N`; N≠0 guarded at line 42; OK | **confirmed coherent** → §4 |
| CL20 | `.cpp:122` — `requires_grad_(true)` on `coord_gpu` only; structural inputs left non-grad; consistent | **confirmed coherent** → §4 |
| CL21 | No sign/transpose/units bug; gradient is raw `dL/dr` magnitude, no claimed physical unit, consistent with "kernel not shielding" | **confirmed coherent** → §4 (note: `_catalog.py` *does* assign units `e^2/Å`; reconciled in §4) |

---

## 3. Edits that don't move numbers

`file:line` references are to the current source.

**Header (`AIMNet2ChargeResponseGradientResult.h`)**

- **E1** `.h:3-22` — Keep the objective/rationale block (it is the strongest
  part of the file per both reviews). Add one lead sentence stating the
  *returned* quantity plainly before the derivation, e.g.
  `// Returned quantity: the per-atom 3-vector dL/d(r_i) and its L2 norm,`
  `// NOT the exact charge-response diagonal d(q_i)/d(r_i).` (Addresses C1, C18-comment-half.)
- **E2** `.h:24-27` — Replace the dated "Per Amendment 2026-05-08(b) … PASSED on
  2026-05-09 …" block with a terse signpost, e.g.
  `// Prerequisite: the .jpt model must support requires_grad on coords.` Move
  the amendment/date/check-path provenance to the commit message or a memory
  entry. (C2, CL15.)
- **E3** `.h:28-40` — Replace "PROMOTED FROM TEST FLAG to always-on … (RunConfiguration.cpp:167)"
  with one lifecycle line, e.g.
  `// Always-on for the non-trajectory and trajectory pipelines; runs its own`
  `// grad-tracking forward (no NoGradGuard) then a single backward.` Keep the
  "depends on AIMNet2Result for attach ordering, does NOT share state" sentence
  — that is current behaviour, not history. The cost note (`~5-6 s`,
  state-sharing deferred) may stay or move to memory; recommend trimming to one
  clause. (C3, CL15.)
- **E16** `.h:3, .h:15` — Replace "charge-polarisation gradient" / "a
  charge-weighted per-atom polarisability" language with the contract-aligned
  phrasing already used by the SDK: charge-response gradient `∂(Σ q²)/∂r`, with
  an explicit `// NOT a Buckingham polarisability α = ∂μ/∂E` aside mirroring
  `_catalog.py:376-377`. (C21, and resolves the §1 staleness.)
- **R1 (note, not an edit)** — Class **not** renamed (C18). The name is a
  provenance/contract surface (see ledger C18). Returned-quantity clarity is
  delivered by E1 instead.

**Implementation (`AIMNet2ChargeResponseGradientResult.cpp`)**

- **E4** `.cpp:66-68` — Add a 3-word signpost on `objective_kind`, e.g.
  `// objective code: 1.0 = L2_of_charges (label string carries the name)`.
  Leave the value `1.0` and the `"L2_of_charges"` label unchanged. (C4.)
- **E5** `.cpp:149-158` — Trim the inline objective comment to the equation
  + the one-line "sum(q) is conservation-degenerate" reason + "see header for
  full rationale." Cut "first-test choice" and the repeated N-backward-passes
  sentence (already in the header). (C5, C15, C26, CL11.)
- **E6** `.cpp:175-179` — Reword the pre-scan comment to keep the
  *why* (a defined-but-NaN/Inf gradient would set TS mask=1 and poison the
  Welford running statistics) and drop the `Codex 2026-05-20 F1` tag. Optional
  short label `// finite-gradient guard` on the loop. (C6, C16, C27, CL13.)
- **E7** `.cpp:80` — Rename `N1` → `padded_atom_count`; add
  `// +1 sentinel padding row (stays zero)` at the declaration so the `+1` is
  explained where it appears, not 10 lines later. Carry the rename through all
  uses (`.cpp:82, 92, 106` tensor shapes). The standalone `// Row N stays zero
  (sentinel padding row).` at `.cpp:90` can then shorten or stay. (C8, CL2, CL5.)
- **E8** `.cpp:82` — Keep name `coord_cpu`; add `// positions in Å` on the
  declaration (rename to `coords_angstrom_cpu` declined — comment is
  lower-risk). (C9.)
- **E9** `.cpp:92-93` — Rename `numbers_cpu` → `atomic_numbers_cpu` and
  `num_acc` → `atomic_num_acc`. (C10, CL3.)
- **E10** `.cpp:105` — Rename `charge_cpu` → `total_charge_cpu`; add
  `// total system charge = 0 (neutral; required model input)`. (C11, CL1.)
- **E11 (optional)** `.cpp:40` — Rename `N` → `atom_count`. Cost: 9 use-sites
  in this function plus the `N1`→`padded_atom_count` pairing. Weighed call;
  recommend **adopting together with E7** so the pair reads `atom_count` /
  `padded_atom_count`. If the human prefers minimal diff, leave `N` (it is a
  near-universal loop-bound convention in this codebase). (C7.)
- **E12** `.cpp:159` — Rename `charges_n` → `real_atom_charges` (sentinel
  dropped). (C12, CL4.)
- **E13** `.cpp:160` — Rename `loss` → `charge_l2_objective`. Carry through to
  `.cpp:162` (`charge_l2_objective.backward()`). (C13, CL11.)
- **E14** `.cpp:196-197` — Rename `v` → `grad_vec`, `s` → `grad_norm`. Light
  touch (splits the C14/CL6 disagreement). (C14.)
- **E15** `.cpp:195` — Add one signpost before the loop:
  `// store per-atom gradient + accumulate norm stats for the log`. (C17, CL8.)
- **E17** `.cpp:74-79` — Replace the `// ---` banner with `// AIMNet2 inputs`
  and keep the one sentence: build duplicates AIMNet2Result's, gradient enabled
  on `coord`, neighbour matrices from the shared static helper. (C23.)
- **E18** `.cpp:137-138` — Trim "Forward pass." redundancy; keep "Gradient tape
  records the graph from coord_gpu through the forward to the charges output."
  (C25.)
- **E19** `.cpp:221, 236` — Add the unit to both WriteFeatures comments:
  `(N, 3) float64 — per-atom charge-response gradient dL/dr [e^2/Å]` and
  `(N,) float64 — L2 norm of the gradient vector [e^2/Å]`. Use the
  contract-aligned name (charge-response gradient, not "polarisability"). (C28.)
- **E20** `.cpp:106` — Add `// single molecule: all atoms in batch 0` on
  `mol_idx_cpu`. (CL1.)

**Not changed (protected / declined):**
- Output NPY filenames `aimnet2_charge_response_gradient.npy`,
  `aimnet2_charge_response_gradient_scalar.npy`.
- `ConformationAtom` fields `aimnet2_charge_response_gradient_vector` /
  `_scalar` (consumed in `ui/`, `h5-reader/`, both TR result classes).
- GeometryChoice record string `"charge_response_gradient_backward"` and label
  `"L2_of_charges"`.
- Class name `AIMNet2ChargeResponseGradientResult` (provenance/contract surface).
- `.detach()` comment (C24, kept per CL14).

---

## 4. Usage notes (the sign/value reasons, traced)

**The L2-of-charges proxy (the central design choice).** *Coherent, earned.*
The objective `L = Σ q_j²` is a deliberate one-backward-pass proxy, not the
exact `dq/dr` diagonal. Reason (documented and verified): `sum(q)` is fixed by
AIMNet2's charge-conservation projection, so `d(sum q)/d(coord) ≈ 0`; the L2
scalar is the cheapest single-pass objective with a non-trivial gradient. The
exact diagonal would cost N backward passes. **Do not flag this as wrong.** It
is stated in the header, at `.cpp:149-158`, and in the SDK (`_catalog.py:372`,
`_tensors.py:807`). Per the brief, the only action is making the comment state
the reason plainly (E5) and aligning the "polarisability" language (E16).

**Sentinel-row exclusion** (`.cpp:159`, CL16). *Coherent.* Tensors are built at
size `N+1`; row N is the zero sentinel padding row (AIMNet2's neighbour-matrix
convention, mirrored from AIMNet2Result). The objective slices
`charges_gpu.slice(0, 0, N)` so the sentinel's charge is excluded from `L`.
Producer-internal; no external consumer sees the sentinel row.

**Squared vs raw cutoff** (`.cpp:108-117`, CL17). *Coherent — traced to producer.*
`BuildNeighbourMatrix` takes a **squared** cutoff (`cutoff_sq`,
`cutoff_lr_sq`) and internally `std::sqrt`s it back (AIMNet2Result.cpp:110);
the model receives the **raw** `cutoff_lr` tensor (AIMNet2Result.cpp:248). This
Result reproduces AIMNet2Result's build line-for-line (`.cpp:108-117` vs
AIMNet2Result.cpp:240-248). Producer and the upstream sibling agree exactly;
the dual use is correct. (Could add a one-word note, but both reviews call it
"correct-looking" — left to the comment trim in E17.)

**float64 cast before finite-scan/norm** (`.cpp:172`, CL18). *Coherent.*
AIMNet2 charges/embeddings are float32 native (`feedback_embedding_float32`);
all derived C++ math is float64. Casting `coord_gpu.grad()` to `kFloat64`
before the `isfinite` scan and norm avoids float32 sentinel ambiguity. Matches
the project rule.

**`requires_grad_(true)` on coord only** (`.cpp:122`, CL20). *Coherent.*
`coord_gpu` is the single leaf differentiated against; `numbers`, `charge`,
`mol_idx`, `nbmat*` are integer/structural inputs correctly left non-grad. The
`.detach().requires_grad_(true)` chain guarantees a leaf tensor on-device
(setting `requires_grad` on a non-leaf raises in libtorch) — the comment
already explains this and is kept.

**`mean_norm = sum_norm / N`** (`.cpp:204`, CL19). *Coherent.* `N == 0` is
rejected at `.cpp:42`; no div-by-zero. `mean_norm`/`max_norm` are log-only
summary stats, not stored or serialized.

**Units `e^2/Å` and the "no claimed unit" reading** (CL21). *Reconciled.*
Claude read the gradient as "no claimed physical unit." In fact the SDK
**does** assign `units="e^2/Å"` (`_catalog.py:379, 381`) — dimensionally
`∂(charge²)/∂(length)`. Both are defensible framings (the value is a raw
geometric kernel, but the SDK records its dimension). No conflict: the producer
emits the raw `dL/dr`; the contract layer annotates its dimension. E19 adds the
unit to the WriteFeatures comments so the producer and contract agree visibly.

**Downstream consumers, traced (all read the same two fields, no sign/transform):**
- `python/nmr_extract/_catalog.py:378-381` — NPY/SDK contract. Already
  contract-correct (`d(sum q_j^2)/d(r_i)`, `e^2/Å`, `1o`/odd parity, renamed
  off "polarisability" 2026-05-20). **Aligned with producer.**
- `python/nmr_extract/_tensors.py:807-847` — `AIMNet2ChargeResponseGradient`
  wrapper; `norms` recomputes `||v||` from the vectors (consistent with the
  stored scalar "up to fp precision"). **Aligned.**
- `python/nmr_extract/_protein.py:1418-1421`,
  `_trajectory.py` (time-series + Welford groups) — pass-through. **Aligned.**
- `src/AIMNet2ChargeResponseGradient{TimeSeries,Welford}TrajectoryResult.cpp` —
  read `_vector`/`_scalar` verbatim into TR H5 groups; their provenance strings
  embed the producer class name and field names. **Aligned** (and the reason
  the class name / field names are treated as protected here).
- `h5-reader/src/model/QtFrameAtomView.h:66-67` — `std::optional<Vec3>` +
  `std::optional<double>`, read-only display. **Aligned.**
- `ui/src/MainWindow.cpp:1485-1493`, `ui/src/RestServer.cpp:1081-1082` —
  inspector display of `dL/dr` and `|dL/dr|`. **Aligned in value.** *Comment
  staleness:* `ui/` still labels it "AIMNet2 polarisability"
  (`MainWindow.cpp:1487, 1490`). This mirrors the stale C++ source comments; it
  is a `ui/`-owned display string, **out of scope for this pass**, but if E16
  lands the human may want a follow-up `ui/` label edit for consistency with
  the 2026-05-20 rename. Not a value/sign issue.

**Verdict on every sign/value the reviews raised: coherent.** No producer or
consumer disagreement on sign, transform, or value. The only producer/consumer
*mismatch* anywhere is the **naming/comment language** (source says
"polarisability"; contract says "charge-response gradient / NOT polarisability")
— a comment fix on the producer side, not a number fix.

---

## 5. Bug-by-exhaustion candidates

**None.** No sign, value, transform, or units bug survives tracing. Every
sign/value finding lands on *coherent* in §4. The "polarisability" naming
mismatch is a stale-comment issue, fixed by E16/E19 without touching any number.

---

## 6. Questions & Ambiguities

- **Q1 (codex C29) — do AIMNet2 model parameters still carry `requires_grad`,
  causing `backward()` to accumulate unused parameter gradients?** Cannot be
  settled from the two files under review. `AIMNet2Result.cpp:50` calls
  `module.eval()` (disables dropout/batchnorm training behaviour) but `.eval()`
  does **not** clear `requires_grad` on parameters. If the exported `.jpt`
  parameters default to `requires_grad=true`, `charge_l2_objective.backward()`
  would populate `.grad` on every model parameter. **Functionally harmless**
  for correctness — only `coord_gpu.grad()` is read, and the parameter grads
  are discarded when the graph is freed — so this is at most a small wasted-
  compute concern, not a value bug. Confirming requires inspecting the `.jpt`
  export or AIMNet2Model load path; out of scope for a readability pass.
  Flagging for the human, not proposing a fix.

- **Q2 — class rename `AIMNet2ChargeResponseGradientResult` →
  `…ChargeL2CoordinateGradientResult` (codex C18).** Declined here because the
  name appears verbatim in two TR provenance strings and matches the NPY/field/
  H5-group naming family, so it reads as a contract surface. If the human
  considers the class name *internal* (it is a C++ type, not a serialized
  string per se), the rename would be more literal but carries the documented
  carry-through. Surfacing the call rather than deciding it unilaterally.

- **Q3 — `ui/` "polarisability" labels** (`MainWindow.cpp:1487, 1490`,
  and the inspector header "AIMNet2 polarisability"). If E16 aligns the
  source/SDK language, the `ui/` display strings become the last place still
  saying "polarisability." `ui/` is a separate ownership boundary (read-only
  consumer of the library). Recommend a follow-up `ui/`-scoped edit, but not
  part of this plan's source changes.

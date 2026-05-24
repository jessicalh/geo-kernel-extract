# Fix plan — Coulomb EFG (`CoulombResult` + `MopacCoulombResult`)

Scope: readability-only fix pass for one algorithm — the vacuum/QM Coulomb
E-field + EFG sum — across both variants:

- `src/CoulombResult.{h,cpp}`   (ff14SB partial charges, spatial cutoff, APBS solvent delta)
- `src/MopacCoulombResult.{h,cpp}` (PM7 Mulliken charges, all-pairs, no solvent)

No source files are modified by this plan. This is the plan only.

---

## 1. Summary

Both files tell a coherent story: classify atoms (backbone / aromatic /
sidechain) → build primary bond directions → per-target sum of point-charge
field `E` and dipolar EFG kernel `V` → apply Coulomb constant → enforce
tracelessness → sanitise → clamp → store features and (Coulomb only) the
APBS solvent delta. The math kernel matches the header documentation
term-for-term, and the units block is accurate.

There is **one real defect, and it is a stale comment, not the code**: the
`CoulombResult.h` banner claims "sum over ALL atoms (no spatial cutoff)" and
"The physics is long-range," but the implementation iterates
`AtomsWithinRadius(pos_i, coulomb_efield_cutoff)` with a **configured 20.0 Å
cutoff** (`data/calculator_params.toml:31`). The cutoff is the intended,
iterated behaviour (it is a named, tuned TOML parameter); the comment is what
drifted. Per the brief ("comments conform to code unless exhaustion proves the
code wrong"), **the fix is to the comments**.

The remaining review findings are genuine readability improvements (terse
names carried through a long loop, a "frac" name that is really a projection,
cryptic `PackST_*` helpers, over-long in-loop essays) and two
audit-field-naming clarifications (`sources_beyond_cutoff`, the
not-stored sidechain EFG). None of them move a number.

**What this pass will touch:** comments, signposts, internal local-variable
names, the two file-static helper names, and one SDK docstring/comment.
**What it will not touch:** the algorithm, the numbers, any NPY field name or
column order, any H5 dataset name, any `_catalog.py` key, the `CoulombScalars`
SDK property names, the clamp/sanitise/traceless logic, the cutoff value, or
the all-pairs-vs-cutoff divergence between the two variants (it is intended —
see §4).

The MOPAC variant is the same kernel with a different charge source; every
shared edit below is applied to both files so they stay in lockstep.

---

## 2. Review-finding ledger

Every finding from `codex.md` and `claude.md`, with disposition. (`codex-correctness.md` is not present.)

### codex.md

| # | Finding | Disposition |
|---|---------|-------------|
| C1 | `CoulombResult.h:34` "no spatial cutoff" contradicts cutoff in `.cpp` | **adopted** → §3 E1 |
| C2 | `.cpp:30-33` formula presents all-source sum; loop is cutoff-limited | **adopted** → §3 E2 |
| C3 | `.cpp:107` "Coulomb sum with spatial cutoff" contradicts header/formula | **adopted** (resolved by E1+E2 making all three agree); banner kept as the *correct* one → §3 E2 |
| C4 | `.cpp:126-326` one long per-target loop; add stage block labels | **adopted (restrained)** → §3 E3 |
| C5 | `.cpp:139` `n_sidechain_aromatic_sources` appears before reason is clear | **adopted** → §3 E4 (diagnostic-only note) |
| C6 | `.cpp:299-325` solvent + shielding are late side-branches; add signposts | **adopted** → §3 E3 |
| C7 | `MopacCoulombResult.cpp:99-100` "Main N^2 Coulomb sum" matches loop — clean | noted, no edit (confirms §4 divergence is real) |
| C8 | `Mopac.cpp:109-264` same long-loop labelling opportunity | **adopted** → §3 E3 (mirrored) |
| C9 | `.cpp:94`/`:101` `d` hides frame/direction → `parent_to_hydrogen` / `atom_to_bond_neighbor` | **adopted** → §3 E5 |
| C10 | `.cpp:156` `r` underspecified → `source_to_target`/`source_to_field` | **adopted (as `r_vec`)** → §3 E6 (see also Cl1) |
| C11 | `.cpp:157` `r_mag` unitless → `source_distance_A` | **declined** — `r_mag` is locally unambiguous next to `r_vec`; the unit is stated in the kernel comment and header units block. Adding `_A` to one local while peers stay bare is inconsistent churn. |
| C12 | `.cpp:163` `E_j` → `source_E` | **declined** — `E_j` mirrors the header math `E_a(i)=…sum_j…`; the `_j` index is the meaning. Renaming breaks the comment↔code symbol match. |
| C13 | `.cpp:166` `V_j` → `source_EFG` | **declined** — same reason as C12; `V_ab` is the header symbol. |
| C14 | `.cpp:288` `E_hat` frame/unit implicit → `total_E_unit` | **declined** — `E_hat` is standard unit-vector notation and the line above computes it from `E_total`; clear in context. |
| C15 | `.cpp:289` `coulomb_E_backbone_frac` is a projection (V/Å), not a fraction | **adopted (comment/docstring only; field name protected)** → §3 E7, §4. Duplicate of Cl2. |
| C16 | `Mopac.cpp:87` `d` direction ambiguity | **adopted** → §3 E5 (mirrored) |
| C17 | `Mopac.cpp:135` `r` ambiguity | **adopted** → §3 E6 (mirrored) |
| C18 | `Mopac.cpp:142` `E_j` → `source_E` | **declined** — as C12 |
| C19 | `Mopac.cpp:145` `V_j` → `source_EFG` | **declined** — as C13 |
| C20 | `Mopac.cpp:246` `mopac_coulomb_E_backbone_frac` same "fraction" issue | **adopted (comment only; field protected)** → §3 E7 |
| C21 | `.cpp:159-167` add named intermediates `inv_r3`,`inv_r5`,`dyadic_rr` | **partially adopted** → §3 E8 keeps `r3`/`r5` (they match the header `r^3`/`r^5`) and adds the comment alignment; declines `dyadic_rr` as an extra local for a one-use `r*r.transpose()`. |
| C22 | `.cpp:188-207` unit conversion + detrace are good distinct stages — clean | noted, no edit |
| C23 | `.cpp:209-245` sanitise + clamp are guards mixed into math; add labels | **adopted** → §3 E3 |
| C24 | `.cpp:269-276` sidechain EFG accumulated but never stored; comment if intentional | **adopted** → §3 E9. Duplicate of Cl4. |
| C25 | `Mopac.cpp:138-146` same compressed kernel; name powers/dyad | **partially adopted** → §3 E8 (mirrored, same restraint as C21) |
| C26 | `Mopac.cpp:163-181` unit-convert then detrace clear — clean | noted, no edit |
| C27 | `Mopac.cpp:231-238` sidechain accumulated not stored; add note | **adopted** → §3 E9 (mirrored) |
| C28 | `.h:48-50` `EFieldAt`/`EFGAt`/`EFGSphericalAt` → `Total…` prefix | **declined** — these return the total by documented design; the header comment already says "total Coulomb E-field." Renaming three public query methods is a cross-file API change (h5-reader/ui mirror these names) for marginal gain. Note carry-through cost; declined as a weighed call. |
| C29 | `.cpp:381` `PackST_C` cryptic → `PackSphericalTensor9` | **adopted** → §3 E10. Duplicate of Cl3. |
| C30 | `Mopac.h:46-48` same total-vs-source ambiguity | **declined** — as C28 |
| C31 | `Mopac.cpp:296` `PackST_MCC` cryptic | **adopted** → §3 E10 (mirrored) |
| C32 | `.cpp:55-58` verbose → `// atom classes` | **adopted (trim)** → §3 E11 |
| C33 | `.cpp:83-87` long → `// bond projection axes` | **adopted (trim)** → §3 E11 |
| C34 | `.cpp:111-112` clear and grounded — acceptable | noted, no edit |
| C35 | `.cpp:198-200` longer than needed → `// traceless projection` | **declined** — the comment's value is *why* (each term traceless by Gauss's law; FP accumulation breaks it). That rationale is the non-obvious part and should stay. Trimming to "traceless projection" loses the reason. |
| C36 | `.cpp:209` spelling `sanitise` vs US norm | **declined** — repo uses British `sanitise` consistently (both files, `OperationLog`); changing it here would be the inconsistency. |
| C37 | `.cpp:307-324` too much prose inside compute loop; move/reduce | **adopted** → §3 E12 (relocate bulk to header, leave terse signpost). Duplicate of Cl5. |
| C38 | `Mopac.cpp:50-54` verbose → `// atom classes` | **adopted (trim)** → §3 E11 (mirrored) |
| C39 | `Mopac.cpp:79-81` good signpost — clean | noted, no edit |
| C40 | `Mopac.cpp:173-174` shorter → `// traceless projection` | **declined** — as C35 |
| C41 | `Mopac.cpp:251-252` concise and grounded — clean | noted, no edit |
| C42 | `.h:34`/`.cpp:141` doc-vs-code cutoff (correctness §6) | **adopted** → §3 E1/E2. Duplicate of C1/C2. |
| C43 | `.cpp:142-143` `sources_beyond_cutoff` means "outside radius," not "excluded from sum" | **adopted (comment + name; internal)** → §3 E13, §4. Note: the *recorded GeometryChoice key* `sources_beyond_cutoff` (`.cpp:254`) is an audit-output name — see §6 Q1. |
| C44 | `.cpp:241-244`/`Mopac.cpp:215-218` clamp rescales E but not EFG; check consumers | **adopted as Usage note (coherent, comment only)** → §4. Duplicate of Cl-clamp note. |
| C45 | `Mopac.cpp:124` all-pairs loop matches comment; no contradiction | noted — confirms §4 divergence is intended on the MOPAC side |
| C46 | headers: no concrete correctness issue beyond cutoff | noted, no edit |

### claude.md

| # | Finding | Disposition |
|---|---------|-------------|
| Cl1 | `.cpp:159-160` `r` is the vector but reads like scalar distance → `r_vec` | **adopted** → §3 E6. Duplicate of C10 (adopting claude's `r_vec` form). |
| Cl2 | `.cpp:289/:282` `coulomb_E_backbone_frac` is a V/Å projection, not a ratio; document at field/store | **adopted (comment/docstring; field name protected)** → §3 E7, §4. Duplicate of C15. |
| Cl3 | `.cpp:381`/`Mopac:296` `PackST_C`/`PackST_MCC` cryptic → `PackSphericalTensor` | **adopted** → §3 E10. Duplicate of C29/C31. |
| Cl4 | `.cpp:202` sidechain EFG detraced but never stored — dead, harmless | **adopted (note)** → §3 E9. Duplicate of C24. |
| Cl5 | `.cpp:307-325` 18-line essay stops the read; move bulk, leave 2-4 word signpost | **adopted** → §3 E12. Duplicate of C37. |
| Cl6 | `.cpp:34-35/.h vs :120,141` header "cutoff" vs comments "no cutoff" contradiction | **adopted** → §3 E1/E2. Duplicate of C1/C2. |
| Cl7 | `.cpp:108` banner contradicts line-34 comment | **adopted** → §3 E2. Duplicate of C3. |
| Cl8 | `.cpp:198-326` ~130-line post-loop tail buries the two physical outputs; one banner separating "physics outputs" from "diagnostics/derived" | **adopted** → §3 E3 |
| Cl9 | `Mopac.cpp:38-264` near-verbatim copy of Coulomb::Compute; add a 3-line "differs only in…" banner | **adopted** → §3 E14 |
| Cl10 | `.cpp:122-124,139,178,295-297` aromatic-source counts feed only a log/one scalar; "diagnostic only" note | **adopted** → §3 E4 |
| Cl11 | `.h:53` `SampleEFieldAt` clear — no change | noted, no edit |
| Cl12 | `.h:48-50` `EFieldAt`/`EFGAt`/`EFGSphericalAt` clear (per-atom quantity) — no change | noted — directly contradicts C28; **C28 declined**, claude's "clean" verdict adopted |
| Cl13 | Both — `E_j`,`V_j`,`q_j` carry units via comment; acceptable — no change | noted — supports declining C12/C13/C18/C19 |
| Cl14 | `.cpp:162-167`/`Mopac:141-146` E/V build well grouped — no change | noted (tension with C21's "compressed"; resolved toward keep — §3 E8 only aligns comments) |
| Cl15 | `.cpp:188-207` unit-convert→detrace correctly sequenced — no change | noted, no edit |
| Cl16 | `.cpp:278-297` ten derived scalars run together; blank-line + 2-word labels per sub-group | **adopted** → §3 E15 |
| Cl17 | `.cpp:209-226`/`Mopac:183-200` sanitise lambdas redefined each loop iter; harmless, reader expects hoisting | **declined as a code change** (no refactor/hoist per brief); **noted** in §4 so a reader knows it is intentional-and-harmless. |
| Cl18 | `Compute` static factory is convention — no change | noted, no edit |
| Cl19 | `.cpp:146` "Self-exclusion via filter framework (not inline check)" editorialises → "// self-exclusion" | **adopted (trim)** → §3 E16 |
| Cl20 | `.cpp:176-177` verbose aromatic-source comment → "// aromatic sidechain source" | **adopted (trim)** → §3 E16 |
| Cl21 | `.cpp:307-324` keep caveat ("T2 only; Buckingham T0 not included"), relocate rest | **adopted** → §3 E12. Duplicate of Cl5/C37. |
| Cl22 | `.cpp:393-394`/`Mopac:308-309` keep physics, drop rev date "EFG schema rev 2026-05-18" | **adopted** → §3 E17 |
| Cl23 | `Mopac.cpp:251-252` good terse signpost — no change | noted, no edit |
| Cl24 | `.cpp:284-286` backbone-projection comment accurate/helpful — no change | noted, no edit (this comment is *why* the name is "frac" but the value isn't — keep it; §4) |
| Cl-corr1 | `.cpp:139,178,297` `n_sidechain_aromatic_sources` has no proximity gate despite "near this atom"; stored per-target but is global | **adopted (comment fix; coherent)** → §3 E18, §4 |
| Cl-corr2 | `.cpp:114`/`Mopac:104` no explicit `r_mag` singularity guard in the per-charge math; relies on `MinDistanceFilter` | **resolved coherent** → §4. `MinDistanceFilter` rejects `distance < 0.1 Å` (`KernelEvaluationFilter.h:153`, `singularity_guard_distance=0.1`), so `r3`/`r5` cannot underflow. No edit. |
| Cl-corr3 | `.cpp:202` sidechain detrace dead — harmless | **adopted (note)** → §3 E9. Duplicate of C24/Cl4. |
| Cl-corr4 | `Mopac.cpp:124` raw N² vs Coulomb's `AtomsWithinRadius`; variants sum over different source sets, not term-for-term comparable | **resolved coherent (intended); comment fix** → §3 E19, §4 + §6 Q2 |
| Cl-corr5 | `.cpp:235`/`Mopac:209` clamp lambda captures `E_mag` pre-scale by value; correct — note for reader | noted; covered by §4 clamp note. No edit. |
| Cl-corr6 | Both — detrace happens after `COULOMB_KE` scaling and after accumulation; consistent with header — no issue | noted, no edit |

---

## 3. Edits that don't move numbers

Shared edits (E*) are applied to **both** `CoulombResult.cpp` and
`MopacCoulombResult.cpp` where the construct exists; line numbers cite the
Coulomb file with the Mopac counterpart in parentheses.

- **E1** `CoulombResult.h:34-35` — replace the false banner. From
  `// Coulomb cutoff: sum over ALL atoms (no spatial cutoff).` /
  `// N^2 is cheap for N < 1000. The physics is long-range.`
  to a statement of the real behaviour, e.g.
  `// Source set: all atoms within coulomb_efield_cutoff (20 Å default, TOML).`
  `// The 1/r^2 field is long-range but truncated at the configured radius.`

- **E2** `CoulombResult.cpp:30-35` (header formula block) — change the sum
  qualifier from `sum_{j!=i}` over an implied all-atom set to note the radius,
  e.g. append `(j within coulomb_efield_cutoff of i)` to the `E_a`/`V_ab`
  lines. Keep the `:107-109` banner "Coulomb sum with spatial cutoff" — it is
  the *correct* one; after E1/E2 the header, formula, and banner all agree.

- **E3** Stage signposts in the per-target loop (Cl8, C4, C6, C23). Add short
  block labels (single-line, no fragmentation) before the existing groups:
  `// --- per-source field + EFG sum ---` (before `for (size_t j ...`),
  `// --- units, tracelessness, sanitise, clamp ---` (before `*= COULOMB_KE`),
  and a single banner before line ~262 separating physics outputs from derived
  diagnostics: `// ===== store: physics outputs, then derived scalars =====`.
  Mirror in Mopac (`:109`, `:163`, `:222`).

- **E4** `CoulombResult.cpp:122` and `:139` — one-line "diagnostic only" note on
  `aromatic_source_count` and `n_sidechain_aromatic_sources`, e.g.
  `// diagnostic counts: feed a log line / one stored scalar, not the sum`
  (C5, Cl10).

- **E5** Rename the local `d` to its direction.
  `CoulombResult.cpp:94` `d` → `parent_to_hydrogen`;
  `:101` `d` → `atom_to_bond_neighbor` (Mopac `:87`, `:93`). Local-only;
  no consumer. (C9, C16)

- **E6** Rename the displacement vector `r` → `r_vec` and keep `r_mag` as its
  norm. `CoulombResult.cpp:156-160` (Mopac `:135-139`). Local-only; mirrors the
  existing `SampleEFieldAt` use of `d`/`r`. (C10, Cl1)

- **E7** Comment-only clarification at the backbone-projection store site that
  the **field name is a projection in V/Å, not a fraction**. `CoulombResult.cpp:289`
  (Mopac `:246`): keep the existing `:284-286` comment (it already explains the
  stability choice) and add one clause naming the unit, e.g.
  `// signed projection (V/Å), not a ratio — field name is historical`.
  **The `coulomb_E_backbone_frac` field name and the `coulomb_scalars` column
  order are NOT changed** (protected; see §4, §6 Q3). (C15, C20, Cl2)

- **E8** `CoulombResult.cpp:162-167` (Mopac `:141-146`) — leave `r3`/`r5` as-is
  (they match the header `r^3`/`r^5`); align the inline comments so the `E_a`
  and `V_ab` comments read with the same indices as the header math. No new
  intermediates, no `dyadic_rr`. (C21/C25 partial, consistent with Cl14.)

- **E9** `CoulombResult.cpp:202-207` / `:269-276` (Mopac `:181`/`:231-238`) —
  one-line note that `EFG_sidechain` is accumulated and detraced for symmetry
  but **not stored or exported** (only total/backbone/aromatic are written),
  e.g. `// EFG_sidechain: kept for symmetry; not stored or written to NPY`.
  (C24, C27, Cl4, Cl-corr3)

- **E10** Rename file-static `PackST_C` → `PackSphericalTensor9` (CoulombResult.cpp:381)
  and `PackST_MCC` → `PackSphericalTensor9` (MopacCoulombResult.cpp:296). Both
  are `static` (file-local, internal linkage) so no collision and no
  cross-file carry-through. Add a one-line layout comment:
  `// pack T0(1) + T1(3) + T2(5) → 9 doubles`. (C29, C31, Cl3)

- **E11** Trim the verbose section comments to terse signposts where the prose
  restates the code: `:55-58` → `// atom classes (backbone / aromatic)`;
  `:83-87` → `// primary bond direction (for E_bond_proj)` (Mopac `:50-54`,
  `:79-81` already concise). (C32, C33, C38)

- **E12** `CoulombResult.cpp:307-325` — relocate the design-narrative bulk of
  the 18-line WARNING comment into the `CoulombResult.h` header rationale
  (where the two-separate-kernels explanation belongs), and leave a 2–4 word
  in-loop signpost: `// T2 only; Buckingham T0 not included here`. (C37, Cl5, Cl21)

- **E13** `CoulombResult.cpp:142-143` — comment on `sources_beyond_cutoff` that
  it counts atoms **outside the search radius**, not atoms excluded from the
  final sum (filters + charge-floor skips reduce the summed set further):
  `// atoms outside coulomb_efield_cutoff; NOT the count dropped by filters/charge floor`.
  Consider renaming the **local** variable to `sources_outside_radius` for
  precision — see §6 Q1 for the recorded-key caveat. (C43)

- **E14** `MopacCoulombResult.cpp:38` (top of `Compute`) — add a 3-line banner
  noting it differs from `CoulombResult::Compute` only in: charge source
  (`mopac_charge` vs `partial_charge`), all-pairs vs cutoff source set, and
  no APBS solvent / no aromatic-source diagnostic count. (Cl9)

- **E15** `CoulombResult.cpp:278-297` — blank-line + 2-word labels between the
  three derived-scalar sub-groups: `// bond projection`, `// backbone alignment`,
  `// aromatic scalars`. Grouping only; no statements moved across the loop. (Cl16)

- **E16** Trim editorialising inline comments:
  `:146` "Self-exclusion via filter framework (not inline check)" → `// self-exclusion`;
  `:176-177` → `// aromatic sidechain source` (the proximity wording is removed
  here because there is no proximity test — see E18). (Cl19, Cl20)

- **E17** `CoulombResult.cpp:393-394` (Mopac `:308-309`) — keep the physics
  ("T2 only; symmetric per charge → T0/T1 structural zeros"), drop the
  `rev 2026-05-18` date stamp. (Cl22)

- **E18** `CoulombResult.cpp:176-178` — fix the "near this atom" wording on
  `n_sidechain_aromatic_sources`: there is no proximity gate; the count is over
  the in-radius source loop for target `i` but applies no extra distance test,
  so the comment must not claim proximity. Reword to
  `// aromatic-and-sidechain sources contributing to this target`. (Cl-corr1; see §4.)

- **E19** `MopacCoulombResult.cpp:99-101` — note in the existing "Main N^2"
  banner that this variant sums over **all pairs** (no `AtomsWithinRadius`
  cutoff), unlike `CoulombResult`, so the two source sets differ by design.
  Comment only. (Cl-corr4; see §4 and §6 Q2.)

---

## 4. Usage notes — the sign/value reasons (the real product)

**(a) The cutoff is real; the header comment is the bug.**
`coulomb_efield_cutoff = 20.0 Å` is a named, tuned TOML parameter
(`data/calculator_params.toml:31`, registered with description "Coulomb
E-field spatial cutoff" in `CalculatorConfig.cpp:54`). `CoulombResult::Compute`
iterates `spatial.AtomsWithinRadius(pos_i, coulomb_cutoff)` (`.cpp:141`). This
is iterated, intended behaviour — the field is long-range but truncated at
20 Å for cost. The header's "no spatial cutoff / physics is long-range" line is
stale prose. **Fix the comment (E1/E2), not the code.** Consumers
(h5-reader `QtAtomInspectorDock`, `QtAtomTimeSeriesDock`; SDK `coulomb_E`,
`coulomb_shielding`) read the stored result and never re-derive the source
set, so the comment fix has zero numerical reach.

**(b) `sources_beyond_cutoff` counts "outside radius," not "excluded from
the sum."** `.cpp:142-143` computes `n_atoms - 1 - neighbours.size()` — the
number of atoms outside the 20 Å search, *before* the `MinDistanceFilter`,
`SelfSourceFilter`, and charge-noise-floor skips further thin the summed set.
The audit field therefore under-reports exclusions relative to its name. This
is a **logging/audit clarity issue, not a physics error**: the field never
enters `E`/`V`, only a `GeometryChoice` record (`.cpp:248-257`). Producer is
internally consistent; the *name and comment* should say "outside radius"
(E13). See §6 Q1 for whether the recorded key may be renamed.

**(c) The clamp rescales E-field but not EFG — and that is correct.** Both
variants clamp `E_total` magnitude to `efield_magnitude_sanity_clamp`
(100 V/Å) and apply the same `scale` to the four E vectors only, never to the
EFG matrices (`.cpp:228-245`, Mopac `:202-219`). I traced both sides:
- The clamp is explicitly an **E-field magnitude** guard — it triggers on
  `E_total.norm()`, records `actual_E_magnitude` in V/Å, and is named
  `efield_magnitude_sanity_clamp`. It is not a tensor guard.
- The EFG has its **own, independent** protection: the per-term tracelessness
  is enforced by `project_traceless` and the `MinDistanceFilter` (0.1 Å)
  prevents `r5` underflow that would blow up `V`. EFG sanity is handled by
  `sanitise_mat` (NaN/Inf → 0), not by the E-field clamp.
- Consumers do not assume joint clamping: h5-reader/ui read `coulombETotal`
  and `coulombShielding` separately; the SDK exposes `coulomb_E` and
  `coulomb_efg_*`/`coulomb_shielding` as separate arrays. No consumer divides
  EFG by E or otherwise couples them.
So E and EFG are guarded **on their own terms**; rescaling E without EFG is the
intended, decoupled design. **No code change; the E3 "E-field clamp" signpost
makes the scope explicit.** (Resolves C44 / Cl-corr5.)

**(d) `coulomb_E_backbone_frac` is a signed projection (V/Å), not a fraction.**
`.cpp:287-289`: `E_backbone.dot(E_hat)` where `E_hat = E_total/|E_total|`. The
existing `:284-286` comment already explains this is bounded by `|E_backbone|`
and chosen for stability near cancellation (vs `|bb|/|total|`). The name is
historical. It flows to `coulomb_scalars.npy` column 2 → SDK
`CoulombScalars.E_backbone_frac` (`_tensors.py:457-459`). **The NPY column and
the SDK property name are protected contract** (the column index is what the
SDK reads). The fix is a comment at the store site (E7) and, optionally, a
docstring note on the SDK property; **not a rename**. See §6 Q3.

**(e) No per-charge singularity guard is needed in the loop.** The per-charge
math (`.cpp:156-167`, Mopac `:135-146`) has no inline `r_mag` guard, unlike
`SampleEFieldAt` (`.cpp:366`). It does not need one: `MinDistanceFilter`
(`KernelEvaluationFilter.h:150-154`) rejects any source with
`distance < singularity_guard_distance` (0.1 Å) before the math runs, so
`r3 = r_mag^3 >= 1e-3` and `r5 >= 1e-5` are bounded away from underflow.
`SampleEFieldAt` guards inline because it does **not** run through the filter
set (it samples at arbitrary grid points). Coherent; no edit. (Resolves Cl-corr2.)

**(f) MOPAC all-pairs vs Coulomb cutoff is an intended divergence.**
`CoulombResult` sums over `AtomsWithinRadius(…, 20 Å)`; `MopacCoulombResult`
sums over the raw `for j in 0..n_atoms` (true N², `Mopac.cpp:124`). The two
variants therefore use **different source sets** and are not strictly
term-for-term comparable, despite the MOPAC header saying "same kernel … same
decomposition" (which refers to the *kernel form*, not the *source set*). The
model learns separate `gamma`/`gamma_mopac`, so the divergence is absorbed by
calibration. I could not find a written rationale for why the MOPAC variant
omits the cutoff. **Treated as intended (comment fix E19), but flagged** —
§6 Q2.

**(g) `n_sidechain_aromatic_sources` has no proximity gate.** The comment says
"near this atom" but the count is incremented for every in-radius aromatic
non-backbone source with no extra distance test (`.cpp:176-178`). Since the
outer source set is already radius-limited to 20 Å, every counted source is
"within 20 Å," so the count is a per-target tally of aromatic-sidechain
sources inside the cutoff — coherent, but the "near" wording overstates a
proximity test that does not exist. It feeds only the stored scalar
`aromatic_n_sidechain_atoms` (not written to NPY in the current 4-scalar
schema) and is otherwise diagnostic. **Comment fix only (E18).** (Resolves
Cl-corr1.)

---

## 5. Bug-by-exhaustion candidates

**None.** Every sign/value the reviews raised landed as **coherent** after
tracing:
- cutoff → real parameter, stale comment (§4a);
- `sources_beyond_cutoff` → correct number, misleading name (§4b);
- clamp E-not-EFG → intended decoupled guard, no consumer couples them (§4c);
- `*_backbone_frac` → real projection, historical name, protected contract (§4d);
- missing inline singularity guard → covered by `MinDistanceFilter` (§4e);
- sidechain EFG detrace → dead-but-harmless, not exported (E9);
- MOPAC all-pairs → intended, calibration-absorbed (§4f, with open question Q2).

No edit in this plan moves a number.

---

## 6. Questions & Ambiguities

- **Q1 — `sources_beyond_cutoff` recorded GeometryChoice key.** §4b/E13 rename
  the **local int** freely. But the same string `"sources_beyond_cutoff"` is
  also the *recorded audit key* (`.cpp:254`, `AddNumber(gc,
  "sources_beyond_cutoff", …)`). Is the GeometryChoice key part of any
  downstream contract (forensics tables, h5-reader GeometryChoice viewer)? If
  it is read by name anywhere, the recorded key should stay and only the
  comment changes; if it is free, rename it to `sources_outside_radius` for
  honesty. I did not find a by-name consumer of this key, but GeometryChoice
  records are surfaced in the UI, so confirm before renaming the key.

- **Q2 — Why does the MOPAC variant omit the 20 Å cutoff?** (§4f) The Coulomb
  variant truncates at 20 Å; the MOPAC variant is full N². Is this a deliberate
  choice (MOPAC charges are conformation-dependent and a different source set
  was wanted), or an oversight from cloning before the cutoff was added to the
  ff14SB path? It does not change correctness (calibration absorbs it), but if
  the two should match, that is a one-line code change outside this readability
  pass and needs your call. Flagged, not changed.

- **Q3 — May the SDK property `CoulombScalars.E_backbone_frac` get a clearer
  name?** The NPY column order is hard contract and stays. The SDK *property*
  name (`_tensors.py:458`) is a softer surface but still a public API for the
  `learn/` calibration and any notebook/script reading the SDK. My default is
  to leave the property name and only add a docstring clause ("signed
  projection in V/Å, not a ratio"). Confirm you don't want a deprecation-style
  alias; I declined to propose one to avoid API churn.

- **Q4 — Review disagreement on the query-method names** (`EFieldAt` etc.):
  codex wants `Total…` prefixes (C28/C30); claude says they are already clear
  (Cl12). I sided with claude and declined (they are documented as returning
  the total, and h5-reader/ui mirror the names — cross-file cost for marginal
  gain). Flagging the disagreement in case you prefer the explicit prefix.

- **Q5 — Header essay relocation (E12).** I propose moving the
  two-separate-kernels narrative from the loop body into `CoulombResult.h`.
  Confirm the header is the right home (vs a `spec/` design note); the brief
  prefers not to bloat headers, but this rationale is genuinely physics and has
  no better in-tree home that a code reader would find.

# SASA — readability fix plan

Targets: `src/SasaResult.h`, `src/SasaResult.cpp`
Reviews: `codex.md` (gpt-5.5 xhigh), `claude.md` (opus). No `codex-correctness.md`.

## 1. Summary

`SasaResult` tells a coherent Shrake-Rupley story top to bottom: build a
Fibonacci-lattice sphere of directions → for each atom expand the radius by the
probe → place test points → test occlusion by neighbours' expanded spheres →
SASA = exposed_fraction × sphere_area, plus an averaged outward normal. The
header documents method, the two TOML parameters, and both output NPYs with
shape/units. Nothing reads as symbol soup.

Both reviews converge on the same low-stakes friction: a few `r_*` names that
silently mix "bare vdW" vs "vdW + probe" conventions, two verbose GeometryChoice
comments, the implicit ownership model (results live on atoms, not in this
object), and one magic-number dependency (`Element::S` = largest Bondi radius).

The fix pass will: improve a handful of internal local names, add a few terse
signposts, and tighten two comments. It will **not** touch the algorithm, the
numbers, the output names (`atom_sasa`, `sasa_normal`), the per-atom storage
model, or any consumer. No sign/value bug found — all consumers read both
outputs straight through (area as scalar Å², normal as a unit vector) with no
compensating transform.

## 2. Review-finding ledger

| # | Source | Finding | Disposition |
|---|--------|---------|-------------|
| C1 | codex §1 (h:38) | `AtomSASA`/`AllSASA` imply this object owns values; only member is `conf_` — add `// values stored on atoms` near `conf_` | **adopted** → E1 |
| C2 | codex §1/§2/§5 (cpp:62-65, :63, :71) | `max_vdw` lacks role/unit; rename `max_occluder_vdw_radius` + `// occluder search radius` | **adopted (partial)** → E4 (comment + targeted rename); see note in §4 on scope |
| C3 | codex §1 (cpp:69-85) | add block label before line 69, e.g. `// exposed surface samples` | **adopted** → E6 |
| C4 | codex §1 (cpp:88-97) | split end block with `// accessible area` (l88) and `// exposed normal` (l91) | **adopted (partial)** → E7 (terse `// accessible area`; normal block already commented at l91-92, tighten not duplicate) |
| C5 | codex §1/§5 (cpp:100-107, :54, :100) | `radii_source` recorded as numeric `0.0` with unit `"Bondi_1964"` — add `// metadata-only choice` | **declined** → this is the `AddNumber(name, value, unit)` provenance convention used identically across calculators; the "unit" slot is a free-text tag and `0.0` is a sentinel. Tightening the block comment (E8) makes the intent clear without flagging a non-issue. |
| C6 | codex §2 (cpp:19) | `FibonacciSphere(int n)` — rename `n` → `point_count` | **adopted** → E2 |
| C7 | codex §2/§4 (cpp:49, :19, :61) | `unit_sphere` holds directions → `sphere_directions`; `FibonacciSphere` returns directions → `FibonacciSphereDirections` | **adopted (partial)** → E2/E3: rename local `unit_sphere` → `sphere_directions`. **Decline** the function rename — `FibonacciSphere` is conventional and the header/comment already says "lattice on unit sphere"; the return type `vector<Vec3>` of unit directions is idiomatic. Weighed: low benefit, file-local but reads fine. |
| C8 | codex §2 (cpp:57) | prefer `atom_pos`/`expanded_atom_radius` over `pos_i`/`r_i` | **adopted (partial)** → E5: rename `r_i` → `r_expanded_i` (also covers Cl1 below). Keep `pos_i` — the `i`/`j` atom-index pairing is the clearest expression of the double loop. |
| C9 | codex §2 (cpp:67) | `exposed` → `exposed_count`/`exposed_samples` | **adopted** → E6b |
| C10 | codex §2 (cpp:72) | `occluded` is clear | **declined (no-op)** → reviewer confirms fine |
| C11 | codex §2/§4 (cpp:118, h:38) | `AtomSASA` doesn't state units (Å²) | **declined** → output-name style; project does not use unit-bearing accessor names (cf. consumers in ui/h5-reader). Reviewer marks "acceptable in this code style." |
| C12 | codex §3 (cpp:23-27) | name Fibonacci intermediates `z`, `radius_xy`, `azimuth` | **adopted (partial)** → E2b: introduce named `theta`/`phi` are already named; add `z` is implicit. Add named intermediates for the polar angle's vertical coordinate readability via a one-line comment rather than restructuring trig. See §3. |
| C13 | codex §3 (cpp:70-77) | add `// occlusion test` signpost | **adopted** → E6c |
| C14 | codex §3 (cpp:89) | introduce `exposed_fraction` named intermediate | **declined** → `sphere_area * exposed / n_points` is one readable line; an intermediate fragments a textbook formula (clarity bar: don't over-explain). The `// accessible area` label (E7) carries the meaning. |
| C15 | codex §4 (h:33) | `Compute` doesn't say it mutates atom fields → `ComputeAtomSasa` or comment | **adopted (comment form)** → E1 covers ownership at the class; add header note that Compute writes to atoms. Decline the rename: `Compute` is the framework-wide `ConformationResult` contract name. |
| C16 | codex §4 / claude §4 (cpp:15) | `BondiRadius` is a one-line passthrough to `BondiVdwRadius`; wrapper adds little | **declined** → harmless one-liner that localizes the Bondi-1964 citation comment; removing it spreads the citation across call sites. Both reviews stop short of demanding removal. No edit. |
| C17 | codex §5 (cpp:54) | `// GeometryChoice: record the parameters used` is process prose | **adopted** → E8 (tighten to `// record SASA parameters as provenance`) |
| C18 | codex §5 (cpp:100) | `// Record a single GeometryChoice summarising...` verbose | **adopted** → E8 |
| C19 | codex §5 (cpp:91-92) | good signpost; could be terser | **declined** → the two-line comment explains the buried-atom zero case, which is the non-obvious part; keep. |
| C20 | codex §5 (cpp:123, :135, :141) | l123 useful keep; l135 restates shape (harmless); l141 useful | **declined (no-op)** → reviewer says keep/acceptable for all three. |
| C21 | codex §6 (cpp:47-49,89) | divide-by-zero / invalid vector size if `sasa_n_points <= 0` | **declined (verified safe)** → see §4: TOML ships `sasa_n_points = 92`, an exact integer; no negative/zero path. Not a runtime hazard in the shipped config. Noted as a question only if config validation is desired. |
| C22 | codex §6 / claude §5 (cpp:63) | comment claims S is largest Bondi radius — verify against table | **verified true** → see §4. Adopt the guard-note suggestion (E4b: `// assumes S = max Bondi`). |
| Cl1 | claude §2 (cpp:59) | `r_i` is probe-expanded, name reads as bare radius → `r_expanded_i` or `// vdW + probe` | **adopted** → E5 (rename `r_i` → `r_expanded_i`); duplicate of C8 |
| Cl2 | claude §2 (cpp:63) | `max_vdw` is bare Bondi while `r_i` includes probe — inconsistent read; add `// bare vdW, probe added separately` | **adopted** → E4 (rename + note); duplicate of C2 |
| Cl3 | claude §2 (cpp:73-77) | `r_j` includes probe, unlike `max_vdw`; same `r_*` prefix two conventions | **adopted** → E5 also renames `r_j` → `r_expanded_j`; the `r_expanded_*` prefix now consistently means "vdW + probe" and `max_vdw` (bare) is visibly distinct |
| Cl4 | claude §3 (cpp:60) | `sphere_area` computed up top, consumed at l89 — slightly out of sequence, harmless | **declined (no-op)** → computing it once per atom near the radius is fine; moving it would not improve clarity. |
| Cl5 | claude §3 (cpp:62-65) | search-radius derivation well isolated; good | **declined (no-op)** → praise |
| Cl6 | claude §4 (cpp:15) | `BondiRadius` passthrough indirection | **duplicate** of C16 |
| Cl7 | claude §5 (cpp:54, :100) | both GeometryChoice comments restate code; shorten | **duplicate** of C17/C18 |
| Cl8 | claude §5 (cpp:63) | `// assumes S = max Bondi` guard note | **adopted** → E4b; duplicate of C22 |
| Cl9 | claude §5 (cpp:91-92, :141, :123) | concise, grounded — keep | **duplicate** of C19/C20 |
| Cl10 | claude §6 (cpp:63) | `max_vdw` recomputed every `i` iteration; hoist outside loop (no result change) | **declined (deferred)** → see §3 note. This is a micro-optimization (constant `BondiRadius(Element::S)`); hoisting is a legitimate no-number-change edit but lies outside readability and the brief forbids refactors beyond clarity. Flag as optional E9. |
| Cl11 | claude §6 (cpp:73-76) | per-neighbour `r_j` and `r_j*r_j` recomputed per test point; could precompute per atom | **declined** → accumulation/cost, not readability or a bug; precomputing per-atom neighbour radii is a perf refactor the brief excludes. |
| Cl12 | claude §6 (cpp:47) | `static_cast<double→int>` truncates not rounds | **declined (verified safe)** → see §4; TOML value is exactly 92. |
| Cl13 | claude §6 (cpp:77) | strict `<` occlusion; confirm boundary (dist == r_j) counts as exposed | **answered** → see §4. Strict `<` means a point exactly on a neighbour's expanded surface is counted exposed. Standard Shrake-Rupley convention; measure-zero set, no effect on results. No edit. |
| Cl14 | claude §6 (cpp:94) | buried-atom guard zeroes normal correctly; consistent with header l19 | **declined (no-op)** → confirms correct |

## 3. Edits that don't move numbers

All edits are internal (locals, comments). No output name, no signature on the
`ConformationResult` contract, no number changes.

- **E1** `SasaResult.h:42` — add `// SASA results are stored on the atoms (conf_), not in this object` next to `const ProteinConformation* conf_`. (C1, C15)
- **E2** `SasaResult.cpp:19` — rename param `n` → `point_count` in `FibonacciSphere` (and its uses at l20, 22, 23, 24). (C6)
- **E2b** `SasaResult.cpp:23` — keep the trig; `theta` (polar) and `phi` (azimuth) are already named. Optionally add `// polar angle / azimuth (Fibonacci spiral)` one-liner above the loop body if a signpost helps. Decline adding `z`/`radius_xy` intermediates — would split a recognizable spherical-coordinate formula. (C12 partial)
- **E3** `SasaResult.cpp:49,70,84` — rename local `unit_sphere` → `sphere_directions`. (C7)
- **E4** `SasaResult.cpp:63` — rename `max_vdw` → `max_occluder_vdw_radius`; change comment `// largest Bondi radius in table` per E4b. (C2, Cl2)
- **E4b** `SasaResult.cpp:63` — comment → `// bare Bondi vdW; assumes S = max in table (PhysicalConstants.h)`. (C22, Cl8, Cl2)
- **E5** `SasaResult.cpp:59,70,75,77` — rename `r_i` → `r_expanded_i` and `r_j` → `r_expanded_j`, so the `r_expanded_*` prefix consistently means "vdW + probe" and is visibly distinct from bare `max_occluder_vdw_radius`. (C8, Cl1, Cl3)
- **E6** `SasaResult.cpp:69` — add block signpost `// exposed-surface sampling`. (C3)
- **E6b** `SasaResult.cpp:67` — rename `exposed` → `exposed_count` (used at l83, 89). (C9)
- **E6c** `SasaResult.cpp:72` — add `// occlusion test` above the neighbour loop. (C13)
- **E7** `SasaResult.cpp:88` — add `// accessible area` above the `ca.atom_sasa = ...` assignment. Leave the existing normal comment (l91-92) as-is; it already labels the block and explains the buried-atom zero. (C4)
- **E8** `SasaResult.cpp:54` → `// record SASA parameters as provenance`; `SasaResult.cpp:100` → `// SASA parameter provenance`. (C17, C18, C5)
- **E9 (optional, perf, not readability)** `SasaResult.cpp:63` — hoist `max_occluder_vdw_radius = BondiRadius(Element::S)` above the `i` loop. No number change. Listed for the human's call; not required by the brief. (Cl10)

Naming note for the human: E2/E3/E4/E5/E6b are all file-local locals/params with
zero cross-file carry-through (verified: no consumer references these symbols —
see §4). They are pure-win internal renames per the brief's "improve internal
names where it helps."

## 4. Usage notes (sign/value tracing)

The reviews raised no sign findings; the value/correctness items traced as
follows. Outputs are `atom_sasa` (scalar Å²) and `sasa_normal` (unit Vec3).

**`atom_sasa` (Å², ≥ 0 by construction).** Producer:
`sphere_area * exposed_count / n_points`, a non-negative fraction times a
positive area. Consumers:
- `python/nmr_extract/_catalog.py:224` — `ArraySpec("atom_sasa", ... units="Å^2", mechanism="solvation")`, read as a plain `np.ndarray`.
- `python/nmr_extract/_protein.py:1571` — `sasa=get("atom_sasa")`, no transform.
- `ui/src/MainWindow.cpp:1453` — displayed as `"%.2f A^2"`, no transform.
- `ui/src/TrajectoryH5.cpp` — Welford mean/std and frame-0 read straight.
- `h5-reader` — `QtFrame::sasa()`, inspector/time-series dock display as `Å²`.
- `learn/bones/*` — Welford-accumulated as-is (`ga.sasa.Update(ca.atom_sasa, ...)`).

Producer and all consumers agree: scalar non-negative area in Å², no sign,
no scaling. **Coherent (expected).**

**`sasa_normal` (unit Vec3, outward).** Producer: `normal_sum / |normal_sum|`,
where `normal_sum` accumulates the *directions* (`sphere_directions[p]`, i.e.
unit lattice vectors pointing from the atom centre outward) of non-occluded test
points. Outward sign is intrinsic — every test point lies on the atom's own
expanded sphere, so its direction from the centre is outward by construction;
exposed (solvent-facing) points dominate the average, giving an outward
resultant. Buried atoms (`|normal_sum| ≤ near_zero_vector_norm_threshold`,
1e-10) get exactly `Vec3::Zero()`. Consumers:
- `_catalog.py:226` — `irreps="1e", tensor_rank=1, parity="odd"` (a polar vector — correct for a spatial direction), read as `VectorField` (N,3).
- `_protein.py:1572` — `sasa_normal=get("sasa_normal")`, no transform.
- `ui/MainWindow.cpp:1454`, `h5-reader QtAtomInspectorDock.cpp:265` — displayed via `vec3Str`/`AddVec3`, no sign flip.
- `learn/bones` — three Welford channels x/y/z, accumulated as-is.

Producer and consumers agree; no consumer negates or re-orients. The
`parity="odd"` SDK tag is consistent with a polar direction vector.
**Coherent (expected).**

**`Element::S` = largest Bondi radius (cpp:63).** Verified against
`src/PhysicalConstants.h:140-149`: `BondiVdwRadius` returns H 1.20, C 1.70,
N 1.55, O 1.52, S 1.80, default 1.70. S (1.80) is strictly the maximum.
The comment is **true today**; E4b adds a guard note pointing at the table so a
future element addition is flagged. **Coherent (expected).**

**`sasa_n_points` truncation (cpp:47).** `data/calculator_params.toml:68` ships
`sasa_n_points = 92` — an exact integer; `static_cast<int>(92.0)` is 92 with no
truncation loss. No zero/negative path in the shipped config. The divide-by-zero
concern (C21) is not reachable as configured. **Coherent (expected)**; a config
validator requiring `> 0` would be a defensive nicety, not a fix (see §6).

## 5. Bug-by-exhaustion candidates

None. No sign or value finding survived tracing; all consumers read both
outputs straight through with no compensating transform, and the two
numeric concerns (S-max, n_points truncation) verified safe against the
shipped table and TOML.

## 6. Questions & Ambiguities

1. **Config validation for `sasa_n_points` (C21, Cl12).** The shipped TOML
   value (92) is safe, but the code does not itself guard `n_points > 0` before
   `FibonacciSphere(n_points)` and the `/ n_points` divide. Is there a central
   `CalculatorConfig` validation layer that rejects non-positive
   `sasa_n_points`, or should a one-line guard be added? This is a defensive
   question, not a bug in the current run — flagging per the brief's "ask
   rather than guess."
2. **Optional perf hoist (E9 / Cl10).** Hoisting the constant
   `BondiRadius(Element::S)` out of the per-atom loop changes no number and is a
   legitimate clarity-neutral cleanup, but it edges toward "refactor." Include it
   in the fix pass, or leave it as the human's call? Listed as optional in §3.

# bonded-energy — claude review (readability focus)

- **Targets:** src/BondedEnergyResult.{h,cpp}
- **Model:** Claude (general-purpose agent, opus)
- **Date:** 2026-05-24
- **Brief:** `../_brief-claude.md`

---

## Verdict

**BondedEnergyResult.h** — Clean. The header tells a coherent story: a tagged-union `BondedInteraction`, the per-type parameter meanings spelled out in a units-annotated comment block, and a result class whose accessors name quantity and units. A force-field chemist reads it top-to-bottom without friction.

**BondedEnergyResult.cpp** — Mostly clean and follows a sensible story: small named `Eval*` functions per energy term (each prefaced by its functional form), a `DihedralAngle` helper, then a dispatch loop that gathers positions, evaluates, and splits energy per atom. The only genuinely dense spot is `EvalCMAP`'s index/weight setup, and there are a couple of single-letter / magic-number snags. Correctness is sound modulo two convention checks (CMAP grid axis order, dihedral sign convention) the file can't self-confirm.

---

## 1. Coherent story / readability

- `BondedEnergyResult.cpp:95-105` — `EvalCMAP` fuses angle→fractional-index→cell-corner→weight in one dense block; as written, a reader must mentally separate "where in the grid" from "how far into the cell." — Add a `// fractional grid coords` signpost at line 95 and `// cell corners + weights` at line 100.
- `BondedEnergyResult.cpp:59-61` — the `m = n1 × b̂2`, `sin_phi = m·n2` trick for signed dihedral is correct but opaque to anyone not recalling the standard construction. — One-line `// signed angle via reference axis` signpost.
- Clean elsewhere: the `Eval*` functions are textbook-legible, each headed by `E = ...`; the dispatch `switch` is uniform and easy to scan.

## 2. Naming carries meaning

- `BondedEnergyResult.cpp:95` — `dx` names a grid *angular step* (radians/cell), not a spatial x; misleads. — `cell_step_rad` or `dphi`.
- `BondedEnergyResult.cpp:96-97` — `fi`, `fj` (fractional grid indices) and `wi`, `wj` (interpolation weights) are terse but locally consistent and explained; borderline-acceptable, leave.
- `BondedEnergyResult.cpp:139,218` — loop index `k` shadows the conventional force-constant `k` used throughout the `Eval*` signatures; momentarily confusing in an energy file. — rename loop var to `a` (atom index).
- `BondedEnergyResult.cpp:136` — `ix` for an interaction is fine and consistent.

## 3. Visible math structure (grouping)

- `BondedEnergyResult.cpp:112-113` — the four-term bilinear sum is grouped well across two lines; clear.
- `BondedEnergyResult.cpp:100-105` — corner indices and weights are interleaved; grouping the two `int i/j` lines then the two `double w` lines (already mostly done) reads cleanly — minor.
- `BondedEnergyResult.cpp:148-211` — each `case` repeats the `PositionAt` gather inline; verbose but the repetition makes each term self-contained and scannable — acceptable, not a finding.

## 4. Function / method naming

- `BondedEnergyResult.cpp:19-86` — `EvalBond`/`EvalAngle`/`EvalUB`/`EvalProperDih`/`EvalImproperDih`/`EvalCMAP` all say what they compute and return energy; good. `DihedralAngle` says quantity + that it returns radians (comment); good.
- Header accessors (`BondEnergy`, `UBEnergy`, etc.) are clear and units-documented at line 86. Clean axis.

## 5. Comments as signposts

- `BondedEnergyResult.cpp:99` — `// Bilinear interpolation (sufficient — CMAP grids are smooth)` editorializes; the "sufficient — smooth" justification is process prose. — trim to `// bilinear interpolation`.
- `BondedEnergyResult.cpp:94` — `// Map angle [-π, π] to grid index [0, spacing-1]` is a good signpost; keep.
- `BondedEnergyResult.cpp:77` — `// Wrap difference to [-π, π]` good; keep.
- `BondedEnergyResult.cpp:95-105` — block lacks the two internal signposts noted in §1; otherwise comment discipline is good and terse throughout.

## 6. Correctness (secondary)

- `BondedEnergyResult.cpp:100-101` — `static_cast<int>(fi) % spacing`: if `phi == +π` then `fi == spacing` exactly and `i0` wraps to 0 with `wi == 0` — correct by wraparound, but check `fi` can't slightly exceed `spacing` (float round-up of `phi`) producing `i0 == 1`-off; a `% spacing` on a value barely above `spacing` is still fine here, so likely OK. Worth a `check`.
- `BondedEnergyResult.cpp:107-110` — grid indexed `grid[i_phi * spacing + j_psi]`; header comment (`.h:65`) says `[phi_idx * spacing + psi_idx]`, consistent. Check this matches GROMACS CMAP storage convention (row-major φ-major) at extraction time — convention the file can't confirm.
- `BondedEnergyResult.cpp:46-62` — `DihedralAngle` sign convention (IUPAC/GROMACS positive direction) can't be confirmed from this file; matters for `EvalProperDih`'s `cos(n·φ − φ0)`. — check against the extraction-side φ0 convention.
- `BondedEnergyResult.cpp:217` — `energy / ix.n_atoms` even split is a documented modeling choice (header line 12), not a bug; `n_atoms` is guaranteed ≥2 by construction so no div-by-zero.
- `BondedEnergyResult.cpp:31-34` — degeneracy guards present: angle clamps `cos_theta`, dihedral guards near-zero normals (`1e-10`), CMAP guards empty grid / `spacing<2`. Bond/UB have no zero-length guard but harmonic forms don't need one. Clean.
- `BondedEnergyResult.cpp:89` vs `196-199` — CMAP guards grid emptiness inside `EvalCMAP` and also bounds-checks `cmap_idx` at the call site; no out-of-range access. Clean.

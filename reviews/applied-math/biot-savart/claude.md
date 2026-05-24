# biot-savart — claude review (readability focus)

- **Targets:** src/BiotSavartResult.{h,cpp}
- **Model:** Claude (general-purpose agent, opus)
- **Date:** 2026-05-24
- **Brief:** `../_brief-claude.md`

---

## Verdict

**BiotSavartResult.h** — Reads cleanly. The header comment is an unusually good through-line: physics model → kernel construction → spherical decomposition → sign convention with a worked numeric check → units. A physicist can follow it once, top to bottom, and know exactly what the class produces and in what frame. The one self-contradiction (header says `G_ab = n_b * B_a` in one place and `-n_b * B_a` in another) is the only thing that forces a re-read.

**BiotSavartResult.cpp** — Mostly coherent and well-signposted. `WireSegmentField` and `JohnsonBoveyField` are exemplary: SI-at-the-boundary discipline is stated and visible, the formula is in the comment, guards are explained. `Compute` is where the story frays — it interleaves the physics (B → G → decompose) with three GeometryChoice recording side-quests, a manual linear-scan find-or-create, and two near-identical cylindrical-coordinate blocks computed twice. A reader must repeatedly step over bookkeeping to keep the physics thread. The inlined `CalculatorConfig::Get("...")` calls in `if` conditions are the single biggest local readability tax.

## 1. Coherent story / readability

- BiotSavartResult.cpp:53,57,165,291,321-324,370,376,400-402 — `CalculatorConfig::Get("literal string")` inlined inside conditionals fuses a config lookup, a magic string, and the comparison into one dense expression; as written, a reader must mentally resolve a runtime string-keyed lookup to know what threshold is being tested — hoist to a named `const double` at the top of the scope (e.g. `const double cutoff = CalculatorConfig::Get("ring_current_spatial_cutoff");`).
- BiotSavartResult.cpp:251-274 vs 286-295 — the cylindrical-frame decomposition (`z`, `d_plane`, `rho`) is computed twice with different local names (`z`/`rho` then `z_coord`/`rho_mag`); a reader cannot tell if the second block is intentionally distinct or accidental duplication — add a one-line signpost on the second block (`// re-derive ring-frame basis for B-field projection`) or compute once.
- BiotSavartResult.cpp:176-196 — the GeometryChoice "ring current" recording block (20 lines, a nested lambda containing a second sampler lambda that re-derives B and G) sits between "get the ring" and "use the ring", interrupting the physics thread; signpost is present but the block dwarfs the computation it precedes — consider a `// (provenance recording; physics resumes below)` marker at its end.
- BiotSavartResult.cpp:258 — `theta = atan2(d_plane.norm(), abs(z))` recomputes `d_plane.norm()` (already in `rho` one line up) and uses `abs(z)`, folding the polar angle into a single call whose branch choice (always positive z half) is not obvious — name it and comment `// polar angle from ring axis, sign-folded`.
- BiotSavartResult.cpp:238-244 — linear scan to find-or-create `RingNeighbourhood` is plumbing inlined into the physics loop; correct but distracts — a one-word signpost `// find-or-create` would let the eye skip it.

## 2. Naming carries meaning

- BiotSavartResult.cpp:44-45 — `a_m`, `b_m`, `r_m`, `dA_m`, `dB_m` are clear (`_m` = metres) but `dA`/`dB` collide visually with the B-field symbol and with each other; the header comment defines them, but in-code they read as cryptic — acceptable given the comment, minor.
- BiotSavartResult.cpp:90 — `halfI_A` is good (half current, amperes); `I_A` at line 45 likewise. No change.
- BiotSavartResult.cpp:272-273 — `cos_phi` / `sin_phi` are clear; the surrounding `d_hat`, `ref_hat`, `ref_plane`, `ref_norm` carry units/frame well. Clean.
- BiotSavartResult.h:64 / cpp:433 — `PackST` / `out` are terse but the comment names the 9-double layout; acceptable for a private packer.
- BiotSavartResult.cpp:302 `ti` (type index) and :170 `ri` (ring index), :161 `ai` (atom index) — single letters, but consistent and obvious in context. No change.

## 3. Visible math structure (grouping)

- BiotSavartResult.cpp:232-235 and :409-412 and :189-192 — the kernel build `G(a,b) = -n_b B_a PPM` appears three times as a raw double loop; each is correct and commented at one site, but the other two are bare — add the one-line `// G_ab = -n_b B_a (see header sign convention)` to the two uncommented copies so the shape is self-evident wherever it appears.
- BiotSavartResult.cpp:312-317 — the finalize sequence (accumulate B/G, decompose total, set shielding contribution) is a clean named-step block; note that `total_G_spherical` (315) and `bs_shielding_contribution` (317) both decompose the same `G_total` family — verify these are intended to be distinct fields, the duplication reads as possibly-redundant.
- BiotSavartResult.cpp:43-111 — `WireSegmentField` → `JohnsonBoveyField` staging is exemplary: build segments → sum loops → return Tesla, with unit conversion isolated at the boundary. Clean.
- BiotSavartResult.cpp:386-418 — `SampleShieldingAt` duplicates `SampleBFieldAt`'s ring-loop guards verbatim then adds the G build; the parallel structure is readable but the shared guard cascade (370-376 ≡ 400-402) is copy-paste — fine for code-complete, no action.

## 4. Function / method naming

- Clean axis. `Compute`, `SampleShieldingAt`, `SampleBFieldAt`, `WriteFeatures`, `WireSegmentField`, `JohnsonBoveyField`, `PackST` all say what they compute and (via comment/signature) what they return and in what units. No vague or misleading names.

## 5. Comments as signposts

- BiotSavartResult.cpp:117 — comment header says cutoff is `CalculatorConfig::Get("ring_current_spatial_cutoff") (15A)` — the parenthetical "15A" duplicates a runtime value that could drift from config; flag as a potential stale-comment hazard — drop the hardcoded "15A" or mark it "(default 15A)".
- BiotSavartResult.cpp:225-231 — the in-loop sign-convention comment with the worked example (`I=-12, G_T0=-0.116 → +1.40 ppm`) repeats the header verbatim; useful here but consider trimming to `// G_ab = -n_b B_a; sign per header` to avoid two copies drifting apart.
- BiotSavartResult.cpp:263-265 — azimuthal-angle comment is a genuinely good signpost (explains *why*: distinguishes N-side from C-side on HIE/TRP). Keep.
- BiotSavartResult.cpp:96,100 — `// Upper loop (z = +d)` / `// Lower loop (z = -d)` use `z` and `d` symbols not defined locally (offset is `lobe_offset_ang`); minor — `// upper loop (+lobe_offset)` is grounded.
- BiotSavartResult.cpp:294 — `B_phi (zero by axial symmetry choice)` is a precise, honest signpost for a non-obvious zero. Keep.
- BiotSavartResult.h:1-32 — exemplary header narrative; no change.

## 6. Correctness (secondary)

- BiotSavartResult.h:10 vs h:23 — header states `G_ab = n_b * B_a * PPM_FACTOR` (line 10, no minus) but then `G_ab = -n_b * B_a * PPM_FACTOR` (line 23); the code (cpp:192,235,412) uses the minus form — the line-10 statement is wrong/stale; fix to match the signed convention.
- BiotSavartResult.cpp:291 vs :269 — `rho_hat` guard uses `CalculatorConfig::Get("near_zero_vector_norm_threshold")` while the earlier in-plane block (269) uses a hardcoded `1e-10`; two different small-rho thresholds for the same geometric degeneracy — check whether the divergence is intentional.
- BiotSavartResult.cpp:370-376 vs Compute — `SampleBFieldAt`/`SampleShieldingAt` apply `distance < geom.radius` and a singularity guard but do NOT apply the `MinDistanceFilter`/`RingBondedExclusionFilter` that `Compute` does; the comment (348-354) states this is intentional (grid points aren't atoms) — consistent, not a bug, but confirm grid and per-atom shielding are meant to differ near rings.
- BiotSavartResult.cpp:303,457,467 — the `8` ring-type bound and `ti < 8` guard are magic numbers tied to the RingTypeIndex enum; if the enum grows, silent truncation occurs at accumulation (303) — check that 8 is enforced as `RingTypeIndex::Count` somewhere, else a future enum addition drops contributions silently.
- BiotSavartResult.cpp:315-317 — `total_G_spherical` and `bs_shielding_contribution` decompose `total_G_tensor` and `G_total` respectively; if `ca.total_G_tensor` carries prior-frame accumulation (`+=` at 314) while `G_total` is per-call, the two fields diverge across repeated Compute calls — check whether Compute is ever re-run on the same conformation.

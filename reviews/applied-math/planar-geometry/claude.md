# planar-geometry — claude review (readability focus)

- **Targets:** src/PlanarGeometryResult.{h,cpp}
- **Model:** Claude (general-purpose agent, opus)
- **Date:** 2026-05-24
- **Brief:** `../_brief-claude.md`

---

## Verdict

**PlanarGeometryResult.h** — Reads cleanly. The header is a four-quantity contract: each of ω, pyramidalization, χ₂, and Cremer-Pople pucker gets a convention-pinned doc block with units, frame, NaN semantics, and a literature anchor. A chemist can learn what the file produces without reading the .cpp. The cross-result-read preamble is the only dense spot, and it's metadata, not math.

**PlanarGeometryResult.cpp** — Tells a coherent story. Four helpers (Dihedral, WrapPi, Pyramidalization, CremerPople5Ring, ThreeBondedNeighbours) are each self-contained with a units/sign signpost, and `Compute` is four numbered, captioned stages in physical order. The only real friction is the Cremer-Pople block, where the comment math and the code's index/scale conventions need careful side-by-side reading to confirm they agree. Naming leans on single letters (A/B/C/D, R1/R2, cs/sn) but each is locally defined and documented, so the soup never spills past its helper.

---

## 1. Coherent story / readability

- .cpp:142 — `cs`/`sn` are the Q·cosθ and Q·sinθ accumulators but read like generic cos/sin temporaries until lines 150-151 rename them `Qcos`/`Qsin` — as written, a reader must scan ahead 8 lines to learn what `cs`/`sn` accumulate — name them `q_cos_acc`/`q_sin_acc` (or just compute into `Qcos`/`Qsin` directly).
- .cpp:145 — `phi` is reused as a different angle (`4π j/5` here vs `2π j/5` at line 125) in two adjacent loops with no signpost on the doubling — add `// 2nd harmonic (m=2 mode)` so the reader sees why the frequency doubled.
- .h:6-14 — the cross-result-read block is dense provenance prose (TR names, policy attrs, dates) at the top of the public header; correct but it's the first thing a reader hits — fine to keep, but consider it reads as machine-audit metadata, not science.
- .cpp clean otherwise: the numbered stage captions (1 pyramidalization, 2 ω, 3 χ₂, 4 pucker) give the through-line at a glance.

## 2. Naming carries meaning

- .cpp:142 — `cs`, `sn` — too terse for accumulators that carry physical meaning (the unscaled puckering projections) — `q_cos_acc`, `q_sin_acc`.
- .cpp:73,121-129 — `n`, `R1`, `R2`, `n_hat` are the Cremer-Pople basis vectors; `R1`/`R2` echo the comment's `R'₁`/`R'₂` so they're defensible, but `n`→`mean_plane_normal` would carry the frame.
- .cpp:70-71 — `A,B,C,D` for the central atom + 3 neighbours is conventional for an improper-dihedral helper and is documented in the block comment; acceptable, no change.
- .cpp:319-320 — `res_i` / `res_ip` (the "ip" = i-plus-one) is slightly cryptic — `res_next` reads at a glance and matches `next_idx`.
- .h naming clean: `omega_actual_`, `omega_deviation_`, `omega_is_xpro_`, `pucker_Q_`, `pucker_theta_` all carry quantity and the doc blocks carry units.

## 3. Visible math structure (grouping)

- .cpp:100-165 `CremerPople5Ring` — three stages (centroid → mean-plane normal → (Q,θ) projection) are present but run as one flat function body; the normal-construction block (121-132) and the projection block (142-164) would benefit from blank-line + 2-word captions (`// mean-plane normal`, `// (Q,θ) projection`) the way `Compute` captions its four stages.
- .cpp:35-46 `Dihedral` — clean: b1/b2/b3 → n1/n2/m1 → atan2 is the textbook sequence, signposted at 31-34.
- .cpp:72-77 `Pyramidalization` — clean: centroid → normal → degeneracy guard → signed projection, four obvious steps.
- .cpp:152-164 — Q computed, then degeneracy guard, then θ; well sequenced and the guard sits exactly at the hazard.

## 4. Function / method naming

- .cpp:100 `CremerPople5Ring` — says ring size and method but not the return; the `PuckerCP{Q,theta}` struct name carries it, so acceptable.
- .cpp:177 `ThreeBondedNeighbours` — predicate-style name returning bool + out-param reads correctly; clear.
- .cpp:35/54/70 `Dihedral`/`WrapPi`/`Pyramidalization` — each says what it computes; units live in the block comment. Clean.
- .h:96 `Compute` / .h:99 `WriteFeatures` — framework-conventional names, fine in context.

## 5. Comments as signposts

- .cpp:49-53 `WrapPi` — the comment is half signpost, half changelog ("replaces a historic while-loop", "MED-2, 2026-05-19", bit-identical-with-DihedralTS) — keep the bit-identity rationale, trim the history to `// wrap to [-π,π]; matches DihedralTS WrapPi (bit-identical at boundary)`.
- .cpp:153-160 — the MED-1 degeneracy comment is good grounded science (cites typical Pro Q range) but long; the literature anchor earns its keep, leave it.
- .cpp:107-120 `CremerPople5Ring` normal block — excellent signpost: states the equation AND the anti-parallel trap the naive edge-cross accumulator falls into; this is the model for the rest of the file.
- .h:36-39 — "The X-Pro NaN-fill described in earlier doc versions was a documentation error…" is history prose in a header convention block — drop to a one-liner or remove; the live behavior is already stated above it.
- .cpp:30-34, 60-69, 81-94 — block-comment captions on each helper are the right shape (units + sign convention up front). Clean.

## 6. Correctness (secondary)

- .cpp:184 — `b.atom_index_a == ai` compares a `decltype(atom_index_a)` to `size_t ai`; if `atom_index_a` is a narrower/signed type this is a sign/width comparison — check the Bond field type matches `size_t`.
- .cpp:193 — sorting neighbours by atom index makes the pyramidalization SIGN build-stable, but the sign then depends on index ordering rather than chemistry; the .h says "right-hand rule on first two neighbour vectors" — confirm the doc's "first two neighbours" means post-sort, else .h:43-45 description and code disagree.
- .cpp:354 — `res.chi[1]` is read as χ₂ (1-indexed χ₂ = `chi[1]`); confirm the `chi` array is 0-indexed-as-χ₁ so `chi[1]` is genuinely χ₂ and not χ₁ — the field name `aromatic_chi2_` makes this load-bearing.
- .cpp:379 — `SaturatedRingAt(ri)` rings with `atom_indices.size() != 5` are silently skipped (left NaN); correct for Pro today, but a 6-ring saturated ring would vanish without a log line — acceptable given current scope, no change needed.
- .cpp:101 vs .cpp:379 — double size guard (caller checks ==5, then `CremerPople5Ring` checks ==5 returning NaN); harmless redundancy, the inner guard makes the helper safe standalone. Fine.
- No accumulation/bounds hazard found in the four `Compute` loops; all index against `assign`-sized vectors and guard `parent >= N_res` / `CA == NONE` / `next_idx` before dereferencing.

# eeq-result — claude review (readability focus)

- **Targets:** src/EeqResult.{h,cpp}
- **Model:** Claude (general-purpose agent, opus)
- **Date:** 2026-05-24
- **Brief:** `../_brief-claude.md`

---

# Adversarial readability review — EeqResult.{h,cpp}

## Verdict

**EeqResult.h** — Clean. The header reads as a complete, self-contained spec: the energy functional, the constrained minimisation, the parameter provenance, the kernel choice, and the output contract are all stated before any code. A chemist who knows EEQ but not this codebase can read it once and know exactly what the class produces.

**EeqResult.cpp** — Genuinely coherent. The five numbered step-banners (CN → χ_eff → Coulomb matrix → constrained solve → store) trace the physics in textbook order, and the dense formulas are pinned to equations with units and citations at the point of use. This is close to the good case for a readability review: the story holds top-to-bottom. The few flags below are local naming mismatches between the comments' textbook symbols and the struct-field names, plus two comments that overstate or mismatch the math they label.

---

### 1. Coherent story / readability

- EeqResult.cpp:140 — diagonal comment says `η_i + √(2/π)/r_i` but code reads `params[i].gam + SQRT_2_OVER_PI / params[i].rad`; reader must trust that `gam`=η and `rad`=r_i with no in-line tie — the only friction in the file is this symbol↔field translation. — add `// gam=η, rad=Gaussian charge radius` once near first use.
- EeqResult.cpp:146-147 — `gam_prod_inv = 1/(η_i·η_j)`, but the Step-3 banner (line 122) writes the off-diagonal as `1/√(R² + 1/(η_i·η_j))` while the header (h:13) writes the same kernel as `1/√(R² + 1/(ηᵢ·ηⱼ))` — consistent with each other, good; just confirm `gam` is the same η used on the diagonal (it is) so the reader isn't tracking two hardness symbols.
- EeqResult.cpp:184-185 — sign chain `λ = -(Q+1'u)/denom` then `q = -(u+λv)` matches the banner (159-161) exactly; reads cleanly, no decode needed.

### 2. Naming carries meaning

- EeqResult.cpp:140,146 — `params[i].gam` is the chemical hardness η used everywhere in the comments as η; the name `gam` invites confusion with `gamma`/γ (the Ohno-Klopman kernel) used three lines down. — the field is external (PhysicalConstants.h) so can't rename, but flag: a local `const double eta_i = params[i].gam;` alias at the diagonal would kill the γ-vs-gam clash. (suggestion only; no math change.)
- EeqResult.cpp:140 — `params[i].rad` is the Gaussian charge radius r_i; `rad` is ambiguous (radius? radians? rcov?), and `rcov` is a sibling field — name does not by itself say which radius. — same alias remedy: `r_gauss`.
- EeqResult.cpp:174 — `u` and `v` are the two solve vectors (A⁻¹χ_eff, A⁻¹1); named in the banner but the bare `u`/`v` carry no meaning at the use sites (184-185). — acceptable given the adjacent banner; minor.
- EeqResult.cpp:40 — `Q_total` clear; `cn_k` (line 41) is the erfc steepness — fine given the Step-1 banner names it `k`.

### 3. Visible math structure

- Clean. The four-block augmented-system banner (153-161) names the algebra before the code executes it, and Step 4 sequences Cholesky → two solves → λ → q → neutrality-shift as named substeps. The shape is obvious.
- EeqResult.cpp:225-231 — the min/max/sum loop is a small unnamed stat block under a `// Charge statistics` banner; fine, but `q_min = q(0)` seeds from the *unclamped* `q` while the loop body reads the *clamped* `conf.AtomAt(i).eeq_charge` (line 227) — mixed sources for the same stat. — seed from `conf.AtomAt(0).eeq_charge` for consistency (see also Correctness).

### 4. Function / method naming

- Clean. `Compute` returns the populated result, `WriteFeatures` returns the file count (2) and the two emitted arrays are commented at each block (265, 273). `Name()`/`Dependencies()` are framework-conventional.

### 5. Comments as signposts

- EeqResult.cpp:55-56 — `// S uses real params; Unknown uses fallback` while the `if` counts BOTH S and Unknown into `n_fallback`; the comment says S is NOT a fallback yet S increments the fallback counter — comment contradicts the code it labels. — either exclude S from the count or rename to `n_special`/relabel the comment.
- EeqResult.cpp:137 — comment `// Diagonal: η_i + √(2/π)/r_i (Caldeweyher 2019 Eq. 8)` is a good signpost; verify the Eq. number — the self-energy term is commonly Eq. 8 in that paper, but "check Eq. 8" since the file cites it as load-bearing.
- EeqResult.cpp:128-133 — the Ohno-Klopman-vs-Gaussian rationale is 6 lines of prose; justified here (it documents a deliberate divergence from the cited reference that changes the numbers) — keep, but it's the one verbose block; could compress the last two sentences to `// geometry-responsive charges, not exact dftd4 reproduction`.
- EeqResult.cpp:117 — `// 1e-14: dftd4 guard against sqrt(0)` is exactly the right terse signpost. Good.
- EeqResult.cpp:187-190 — neutrality-shift comment is clear and grounded ("residual proportional to cond(A)"); good.

### 6. Correctness (secondary)

- EeqResult.cpp:128-147 — comment block (128-133) describes the off-diagonal as the Ohno-Klopman `1/√(R²+1/(η_iη_j))` AND the diagonal comment (137-139) describes a Gaussian self-energy `√(2/π)/r`; these are two different charge models mixed in one matrix. The code matches the comments, and the header (h:6-13) commits to Ohno-Klopman — so this is a documented deliberate hybrid, not a bug. Flag only so a future reader doesn't "fix" the diagonal to match the off-diagonal model.
- EeqResult.cpp:225 — `q_min`/`q_max` seeded from unclamped `q(0)` but compared against clamped `conf.AtomAt(i).eeq_charge`; if atom 0 is clamped and is the extreme, the reported min/max can be the pre-clamp value, inconsistent with the stored charges. — seed both from `conf.AtomAt(0).eeq_charge`.
- EeqResult.cpp:191-192 — uniform neutrality shift divides by `N`; if `N==0` this is div-by-zero, and the CN/matrix loops would also be vacuous. No guard for empty conformation. — check whether an upstream invariant guarantees `N>0`; if not, early-return on `N==0`.
- EeqResult.cpp:184 — sign convention: `q = -(u + λv)` yields charges from the augmented system `[A 1; 1' 0][q;λ]=[-χ_eff; Q]`. The banner (156-157) writes the RHS top block as `-χ_eff`, consistent. Convention looks self-consistent; can't independently confirm the EEQ sign of χ from this file alone — check against D4EeqParams sign of `chi`.
- EeqResult.cpp:55 — `n_fallback` is computed and logged-into via `++` but never read/reported anywhere (not in either GeometryChoice record nor the Info log). — dead accumulation; either surface it in the statistics record or drop it.

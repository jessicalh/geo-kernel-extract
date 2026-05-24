# eeq-result — codex review (readability focus)

- **Targets:** src/EeqResult.{h,cpp}
- **Model:** codex `gpt-5.5`, reasoning `xhigh`
- **Date:** 2026-05-24
- **Tokens:** 9,127
- **Brief:** `../_brief.md` (readability-first)
- **See also:** `codex-correctness.md` (earlier correctness-first pass)

---

**Verdict**
`src/EeqResult.h` is mostly coherent: the top comment tells the physical story clearly, but it is heavier than the tiny interface it documents and it is now slightly stale. `src/EeqResult.cpp` has a strong through-line because the computation is explicitly staged, but the solve block and parameter/matrix naming still require decoding algebra rather than reading intent. The file is not symbol soup overall, but several dense names and comments make a human reader carry too much context.

**1. Coherent Story / Readability**
`src/EeqResult.h:23` — output story says only `eeq_charges.npy`, but implementation also writes `eeq_cn.npy` — add `eeq_cn.npy` here or remove output details from the header.

`src/EeqResult.cpp:14` — large banner says “D4 EEQ parameters”, then lines 26–27 say the parameters live elsewhere — replace the banner with a short `// D4 EEQ parameter source` note or move this prose to `PhysicalConstants.h`.

`src/EeqResult.cpp:51` — `n_fallback` is computed but never used; as written, a reader must wonder whether a warning or summary path is missing — remove it or record/log it.

`src/EeqResult.cpp:153` — solve stage is mathematically compact but not narrative enough; as written, a reader must mentally derive the sign convention for `u`, `v`, `lambda`, and `q` — add named intermediates like `A_inv_chi`, `A_inv_ones`, `constraint_scale`, `charges`.

`src/EeqResult.cpp:187` — “Enforce exact charge neutrality” is followed by clamping later, which can change the final sum — reword to `// remove solve residual before clamping` or record final post-clamp neutrality explicitly.

**2. Naming Carries Meaning**
`src/EeqResult.cpp:32` — `N` is conventional but low-information for a top-level algorithm — `atom_count`.

`src/EeqResult.cpp:40` — `Q_total` carries meaning, but style differs from nearby C++ locals — `total_charge_e`.

`src/EeqResult.cpp:41` — `cn_k` hides what it controls — `cn_steepness`.

`src/EeqResult.cpp:86` — `cn` is understandable only after the comment — `coordination_number`.

`src/EeqResult.cpp:95` — `R2` omits frame/units — `distance_bohr_sq`.

`src/EeqResult.cpp:135` — `A` is mathematically standard but weak for a chemist reading once — `coulomb_hardness_matrix`.

`src/EeqResult.cpp:174` — `u` and `v` are the worst names in the file; they hide the block-elimination story — `A_inv_chi_eff` and `A_inv_ones`.

`src/EeqResult.cpp:184` — `lambda` is fine mathematically, but clearer as `charge_constraint_lambda`.

`src/EeqResult.cpp:185` — `q` is fine in formulas, but for storage code `charges` would reduce context switching.

**3. Visible Math Structure**
`src/EeqResult.cpp:100` — distance calculation is duplicated in both CN and Coulomb blocks; not asking for abstraction, but a local naming pattern would help — name the three quantities consistently as `distance_bohr_sq`, `distance_bohr`, `cov_radius_sum_bohr`.

`src/EeqResult.cpp:146` — `gam_prod_inv` fuses parameter semantics and kernel algebra — use a named intermediate such as `hardness_product_inverse`.

`src/EeqResult.cpp:184` — constraint projection is compressed into two lines — split into `constraint_numerator`, `constraint_denominator`, `charge_constraint_lambda`, then `charges`.

`src/EeqResult.cpp:225` — statistics block starts from `q(0)` but then reports clamped atom charges from `conf` — initialize from the stored first atom charge, or name this clearly as `stored_charge_min/max`.

**4. Function / Method Naming**
`src/EeqResult.h:48` — `Compute` is acceptable in the local result pattern, but not self-describing outside it — clearer standalone name would be `ComputeEeqCharges`.

`src/EeqResult.h:50` — `WriteFeatures` is probably framework terminology and acceptable if inherited; locally, it writes charges and CN.

`src/EeqResult.h:42` — `Name()` returning `"EeqResult"` is clear enough.

**5. Comments As Signposts**
`src/EeqResult.cpp:38` — “Named constants from TOML (no buried literals)” is process commentary — replace with `// configuration`.

`src/EeqResult.cpp:48` — good signpost, slightly long — `// atom parameters`.

`src/EeqResult.cpp:69` — “GeometryChoice: parameter summary” is clear enough.

`src/EeqResult.cpp:80` — good block label; could be terser — `// coordination numbers`.

`src/EeqResult.cpp:109` — good block label — no change needed.

`src/EeqResult.cpp:128` — long justification about deliberate simplification and acceptability reads like process prose — replace with a terse note such as `// Ohno-Klopman approximation`; keep validation rationale in docs, not the kernel.

`src/EeqResult.cpp:137` — useful physical signpost, but the “Makes the matrix...” sentence is explanatory prose — replace block with `// self-Coulomb diagonal`.

`src/EeqResult.cpp:163` — SPD explanation is too confident for an inline comment unless validated elsewhere — replace with `// Cholesky solve`.

`src/EeqResult.cpp:187` — too much prose and slightly misleading after clamp — replace with `// solve residual correction`.

`src/EeqResult.cpp:210` — good signpost, but “traceable decision” is wordy — `// charge clamp record`.

**6. Correctness**
`src/EeqResult.cpp:191` / `src/EeqResult.cpp:205` — possible final charge-neutrality issue: charges are neutralized before clamping, but clamping individual atoms can change `q_sum` away from `Q_total` — if exact post-clamp neutrality is required, check this path; otherwise make the comment say neutrality is only pre-clamp.

`src/EeqResult.cpp:225` — `q(0)` is read without an empty-conformation guard — if `AtomCount()` can ever be zero, this is an actual crash; if zero-atom conformations are impossible by contract, document that assumption near line 32.

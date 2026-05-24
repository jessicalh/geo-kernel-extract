# eeq-result — codex review

- **Targets:** src/EeqResult.{h,cpp}
- **Model:** codex `gpt-5.5`, reasoning `xhigh`
- **Date:** 2026-05-24
- **Tokens:** 32666
- **Brief:** `../_brief.md`

---

**1. Algorithm Correctness**
src/EeqResult.cpp:205 — post-solve clamping can break `Σq = Q_total` after neutrality was enforced — re-neutralize after clamp, or explicitly mark clamped outputs as non-neutral.

src/EeqResult.cpp:225 — `q_min`/`q_max` initialize from unclamped `q(0)` but loop reads clamped stored charges — initialize from `conf.AtomAt(0).eeq_charge`.

src/EeqResult.cpp:43 — negative `eeq_charge_clamp` would flip clamp signs at lines 205–207 — guard `charge_clamp > 0`.

src/EeqResult.cpp:146 — I don’t know whether the intended Ohno convention is geometric-mean hardness; code gives `γ(0)=sqrt(gam_i gam_j)` — check against the chosen EEQ validation fixture.

**2. Variable Naming Clarity**
src/EeqResult.cpp:40 — `Q_total` hides units — `total_charge_e`.

src/EeqResult.cpp:51 — `n_fallback` is misleading and unused; it also counts sulfur despite real parameters — `n_special_params` or remove.

src/EeqResult.cpp:59 — `pos` hides Bohr units — `pos_bohr`.

src/EeqResult.cpp:95 — `R2` hides units — `r2_bohr`.

src/EeqResult.cpp:135 — `A` hides the physical role — `charge_hessian`.

src/EeqResult.cpp:146 — `gam_prod_inv` is cryptic — `ohno_softening_sq`.

src/EeqResult.cpp:174 — `u`/`v` obscure solve meaning — `chi_response`, `constraint_response`.

**3. Grouping Of Math Operations**
src/EeqResult.cpp:115 — CN guard and EN shift are folded into one expression — name `cn_guarded` before `sqrt`.

src/EeqResult.cpp:146 — Ohno softening and kernel evaluation are compressed — split as `softening_sq` then `gamma_ij`.

src/EeqResult.cpp:184 — Lagrange multiplier solve is legible; no grouping issue.

**4. Method / Function Naming**
src/EeqResult.h:42 — `Name()` is framework-standard and clear — no issue.

src/EeqResult.h:48 — `Compute()` is generic but consistent with result classes — no issue.

src/EeqResult.cpp:259 — `WriteFeatures()` accurately describes output writing — no issue.

**5. Comments**
src/EeqResult.h:6 — energy comment omits the diagonal self-Coulomb term used at cpp:140 — replace with `// EEQ charge model`.

src/EeqResult.h:23 — output comment omits `eeq_cn.npy` — replace with `// EEQ outputs`.

src/EeqResult.cpp:38 — decorative banner is verbose — `// config constants`.

src/EeqResult.cpp:48 — decorative banner is verbose — `// positions in Bohr`.

src/EeqResult.cpp:80 — verbose step label — `// coordination numbers`.

src/EeqResult.cpp:119 — verbose step label — `// charge Hessian`.

src/EeqResult.cpp:137 — diagonal-dominance/SPD claim is stronger than shown — `// diagonal self term`.

src/EeqResult.cpp:187 — long explanatory neutrality comment — `// neutrality residual`.

src/EeqResult.cpp:210 — comment restates the record call — `// charge clamp`.

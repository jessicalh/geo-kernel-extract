# coulomb — codex review (readability focus)

- **Targets:** src/CoulombResult.{h,cpp} + src/MopacCoulombResult.{h,cpp}
- **Model:** codex `gpt-5.5`, reasoning `xhigh`
- **Date:** 2026-05-24
- **Tokens:** 17,036
- **Brief:** `../_brief.md` (readability-first)

---

**Verdict**

`src/CoulombResult.h` mostly tells a coherent top-level story: vacuum Coulomb field/EFG, source decomposition, solvent delta, units. The main break is that it explicitly says there is no cutoff, while the implementation uses a configured cutoff.

`src/CoulombResult.cpp` has the right computational order, but it reads partly as a production math kernel and partly as accumulated instrumentation. A reader can follow it, but they must reconcile contradictory cutoff claims and carry several terse names through a long per-atom loop.

`src/MopacCoulombResult.h` is coherent and easier than the Coulomb header: it explains the only real distinction, the charge source. No major narrative break.

`src/MopacCoulombResult.cpp` is the clearest of the two implementations because it is a straight all-pairs loop. Its main readability cost is inherited compression: terse `E_j`, `V_j`, `r`, `d`, `E_hat`, plus storage/output sections whose schemas are only partially named.

**1. Coherent Story / Readability**

`src/CoulombResult.h:34` — says “sum over ALL atoms (no spatial cutoff)” but `CoulombResult.cpp:107-120` uses `coulomb_efield_cutoff` and `AtomsWithinRadius` — make the header match reality, e.g. `// Coulomb source radius: configured spatial cutoff`.

`src/CoulombResult.cpp:30-33` — formula presents an all-source sum, but the actual loop is cutoff-limited — add a short qualification: `// evaluated over configured source radius`.

`src/CoulombResult.cpp:107` — “Coulomb sum with spatial cutoff” is clear, but it contradicts the header and earlier formula; this is the worst story break. As written, a reader must decide whether the physics is intentionally long-range or intentionally truncated.

`src/CoulombResult.cpp:126-326` — one long per-target loop combines kernel evaluation, classification, unit conversion, projection, sanitizing, clamping, audit logging, storage, derived features, APBS subtraction, and shielding notes — add terse block labels before each stage and consider grouping existing statements more visibly; no computation change needed.

`src/CoulombResult.cpp:139` — `n_sidechain_aromatic_sources` appears before the reader knows why it exists — rename or add `// aromatic source count` before it.

`src/CoulombResult.cpp:299-325` — solvent subtraction and shielding contribution are late side branches after storage has already begun — add signpost labels such as `// solvent delta` and `// EFG shielding tensor`.

`src/MopacCoulombResult.cpp:99-100` — “Main N^2 Coulomb sum” is coherent and matches the loop — clean on this point.

`src/MopacCoulombResult.cpp:109-264` — same long-loop readability issue as Coulomb, though less severe — clearer block labels would let a reader see: target setup → source sum → units → traceless projection → guards → store features.

**2. Naming Carries Meaning**

`src/CoulombResult.cpp:94` — `d` hides frame and direction — suggest `parent_to_hydrogen`.

`src/CoulombResult.cpp:101` — `d` hides frame and direction — suggest `atom_to_bond_neighbor`.

`src/CoulombResult.cpp:156` — `r` is physically central but underspecified — suggest `source_to_target` or `source_to_field`.

`src/CoulombResult.cpp:157` — `r_mag` is readable but unitless — suggest `source_distance_A`.

`src/CoulombResult.cpp:163` — `E_j` is compact but makes the reader remember it is one source’s raw field contribution — suggest `source_E`.

`src/CoulombResult.cpp:166` — `V_j` is domain-recognizable only if the reader already maps `V` to EFG — suggest `source_EFG`.

`src/CoulombResult.cpp:288` — `E_hat` is recognizable but frame/unit implicit — suggest `total_E_unit`.

`src/CoulombResult.cpp:289` — `coulomb_E_backbone_frac` is misleading: the value is a signed projection in V/A, not a fraction — suggest `coulomb_E_backbone_projection`.

`src/MopacCoulombResult.cpp:87` — `d` has the same direction ambiguity — suggest `parent_to_hydrogen`.

`src/MopacCoulombResult.cpp:135` — `r` has the same central ambiguity — suggest `source_to_target`.

`src/MopacCoulombResult.cpp:142` — `E_j` — suggest `source_E`.

`src/MopacCoulombResult.cpp:145` — `V_j` — suggest `source_EFG`.

`src/MopacCoulombResult.cpp:246` — `mopac_coulomb_E_backbone_frac` has the same “fraction but not fraction” issue — suggest `mopac_coulomb_E_backbone_projection`.

**3. Visible Math Structure**

`src/CoulombResult.cpp:159-167` — inverse powers, field term, and EFG kernel are visible but compressed — add named intermediates such as `inv_r3`, `inv_r5`, `dyadic_rr`; as written, a reader must parse algebra and units at once.

`src/CoulombResult.cpp:188-207` — unit conversion and traceless projection are good distinct stages — clean.

`src/CoulombResult.cpp:209-245` — sanitizing and clamping are algorithmic guards mixed into the math narrative — add labels `// finite-value guard` and `// E-field clamp`.

`src/CoulombResult.cpp:269-276` — total/backbone/aromatic EFG storage and spherical decomposition are structurally clear, but sidechain EFG is accumulated and never stored here — add a comment if intentionally omitted from features/storage.

`src/MopacCoulombResult.cpp:138-146` — same compressed kernel construction — name inverse-distance powers or dyadic product.

`src/MopacCoulombResult.cpp:163-181` — unit conversion then traceless projection is clear — clean.

`src/MopacCoulombResult.cpp:231-238` — total/backbone/aromatic decomposition is clear; sidechain is accumulated but not stored — add `// sidechain EFG not exported` if intentional.

**4. Function / Method Naming**

`src/CoulombResult.h:48` — `EFieldAt` returns total Coulomb E-field, not a selected component/source — suggest `TotalEFieldAt`.

`src/CoulombResult.h:49` — `EFGAt` returns total EFG — suggest `TotalEFGAt`.

`src/CoulombResult.h:50` — `EFGSphericalAt` returns total spherical EFG — suggest `TotalEFGSphericalAt`.

`src/CoulombResult.cpp:381` — `PackST_C` is cryptic and file-local but still opaque — suggest `PackSphericalTensor9`.

`src/MopacCoulombResult.h:46-48` — same total-vs-source ambiguity as Coulomb — suggest `TotalEFieldAt`, `TotalEFGAt`, `TotalEFGSphericalAt`.

`src/MopacCoulombResult.cpp:296` — `PackST_MCC` is cryptic — suggest `PackSphericalTensor9`.

**5. Comments As Signposts**

`src/CoulombResult.cpp:55-58` — verbose process wording — terse replacement: `// atom classes`.

`src/CoulombResult.cpp:83-87` — useful but long — terse replacement: `// bond projection axes`.

`src/CoulombResult.cpp:111-112` — clear and grounded — acceptable.

`src/CoulombResult.cpp:198-200` — useful but longer than needed in-loop — terse replacement: `// traceless projection`.

`src/CoulombResult.cpp:209` — good signpost, but spelling differs from likely codebase norms if US English is used — `// finite-value guard`.

`src/CoulombResult.cpp:307-324` — too much explanatory/process prose inside the compute loop — move the long warning to a header/design note or reduce in place to `// T2 EFG kernel only`.

`src/MopacCoulombResult.cpp:50-54` — useful but long — terse replacement: `// atom classes`.

`src/MopacCoulombResult.cpp:79-81` — good signpost — clean.

`src/MopacCoulombResult.cpp:173-174` — good content but can be shorter — `// traceless projection`.

`src/MopacCoulombResult.cpp:251-252` — concise and grounded — clean.

**6. Correctness**

`src/CoulombResult.h:34` / `src/CoulombResult.cpp:141` — possible documentation/spec bug: header says no cutoff, implementation excludes sources outside `coulomb_efield_cutoff`. I would not call this a physics bug from the snippet alone; check intended production semantics.

`src/CoulombResult.cpp:142-143` — `sources_beyond_cutoff` subtracts `neighbours.size()` before filters and charge-floor skips, so the audit count means “outside radius result set,” not “excluded from the final Coulomb sum.” If that is the intended log meaning, rename the recorded field or comment it.

`src/CoulombResult.cpp:241-244` / `src/MopacCoulombResult.cpp:215-218` — clamp rescales E-field components but not EFG. That may be intentional because it is an E-field clamp; check downstream consumers do not assume both were guarded together.

`src/MopacCoulombResult.cpp:124` — all-pairs loop matches the comment; no cutoff contradiction seen.

`src/CoulombResult.h`, `src/MopacCoulombResult.h` — no concrete correctness issue seen beyond the Coulomb cutoff contradiction.

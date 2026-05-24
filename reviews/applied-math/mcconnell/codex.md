# mcconnell — codex review (readability focus)

- **Targets:** src/McConnellResult.{h,cpp} + src/MopacMcConnellResult.{h,cpp}
- **Model:** codex `gpt-5.5`, reasoning `xhigh`
- **Date:** 2026-05-24
- **Tokens:** 20,994
- **Brief:** `../_brief.md` (readability-first)

---

**Verdict**
`src/McConnellResult.h` mostly tells a coherent story: it states the kernel, units, and stored quantities up front. The public query names are less clear than the physics comment, especially scalar-vs-tensor meaning.

`src/McConnellResult.cpp` is followable at the top level: setup → nearby bonds → filters → kernel → accumulation → spherical features. The story breaks hardest around the category `T2` totals, where names/comments say “K” or “dipolar kernel” but the code projects the symmetric traceless part of the full `M` tensor.

`src/MopacMcConnellResult.h` is coherent as “McConnell plus Wiberg bond-order weight.” It becomes slightly over-explanatory and overconfident about zero bond order.

`src/MopacMcConnellResult.cpp` mirrors the McConnell file clearly enough, and the bond-order weighting is visible. The same `T2`/`K` ambiguity remains, and the bond-order gate changes the population of “nearest” bonds without making that scope obvious.

**1. Coherent Story / Readability**
`src/McConnellResult.cpp:86-92` — dense tensor expression fuses all three physical terms; as written, a reader must decode indices, signs, and tensor roles at once — split into named local terms or add term labels beside each line.

`src/McConnellResult.cpp:198-240` — one switch mixes scalar sums, tensor category totals, and nearest-bond tracking — add signposts such as `// scalar sums`, `// tensor totals`, `// nearest CO/CN`.

`src/McConnellResult.cpp:267-285` — comment says “from symmetric dipolar kernel sums,” but code symmetrizes and de-traces accumulated full `M` tensors — say explicitly whether this is `T2(full M)` or accumulated dipolar `K`.

`src/MopacMcConnellResult.cpp:151-153` — bond-order filtering happens before geometry filtering; as written, the reader has to infer that “nearest” means nearest above the bond-order floor — add `// bond-order gate`.

`src/MopacMcConnellResult.cpp:251-262` — same `T2(full M)` vs `K` story break as McConnell — make the projected tensor’s source explicit.

**2. Naming Carries Meaning**
`src/McConnellResult.cpp:45-50` — `M_over_r3`, `K`, and especially `f` are compact symbols rescued only by comments — prefer `full_tensor_over_r3`, `dipolar_kernel`, `mcconnell_scalar`.

`src/McConnellResult.cpp:61` — `d` hides the frame — `midpoint_to_atom` would carry direction.

`src/McConnellResult.cpp:142` — `sidechain_sum` only accumulates `SidechainCO`, not all sidechain categories — use `sidechain_co_sum`.

`src/McConnellResult.cpp:143-146` — `M_total`/category totals are already divided by `r^3` — use `M_over_r3_total` or `shielding_tensor_total`.

`src/McConnellResult.cpp:153,206` — `best_co_direction` stores `kernel.direction`, i.e. midpoint-to-atom direction, not bond direction — use `best_co_midpoint_to_atom`.

`src/McConnellResult.cpp:275,279,283` — `K_backbone` etc. conflict with earlier `K` meaning “dipolar kernel” — use `backbone_T2_tensor` or `backbone_traceless_M`.

`src/MopacMcConnellResult.cpp:152-153` — `bo` and `zero_bo_skipped` are terse and threshold-based, not strictly zero — use `bond_order` and `below_floor_pairs_skipped`.

**3. Visible Math Structure**
`src/McConnellResult.cpp:74-92` — scalar, dipolar kernel, and full tensor are present, but the full tensor block is the only one not visibly decomposed into named terms — name the three terms from the formula.

`src/McConnellResult.cpp:269-272` — `project_traceless` is defined but unused, then the same operation is repeated inline — either use it consistently or remove it.

`src/McConnellResult.cpp:372-398` — feature widths `9`, `25`, and `6` are magic numbers tied to schema — add named dimensions or a terse schema label.

`src/MopacMcConnellResult.cpp:179-181` — bond-order weighting is nicely isolated and readable — clean on this sub-block.

`src/MopacMcConnellResult.cpp:303-333` — same schema magic numbers as McConnell — name dimensions or add a compact category-order comment.

**4. Function / Method Naming**
`src/McConnellResult.h:51` — `CategorySum` does not say it returns the McConnell scalar sum — use `CategoryScalarSum`.

`src/McConnellResult.h:52` — `NearestCOContribution` hides that this is scalar `f`, not tensor contribution — use `NearestCOScalarContribution`.

`src/McConnellResult.h:55` — `SampleShieldingAt` is acceptable in class context, but clearer as `SampleMcConnellShieldingAt`.

`src/McConnellResult.cpp:54` / `src/MopacMcConnellResult.cpp:47` — `ComputeBondKernel` is generic — `ComputeMcConnellBondKernel` would make copied/local use clearer.

`src/McConnellResult.cpp:362` / `src/MopacMcConnellResult.cpp:293` — `PackST_MC` and `PackST_MMC` are cryptic — use `PackSphericalTensor9`.

**5. Comments As Signposts**
`src/McConnellResult.cpp:119-124` — comment omits `MinDistanceFilter` — replace with `// distance/self/near-field filters`.

`src/McConnellResult.cpp:172` — `// ---- GeometryChoice: filter exclusion ----` is process-heavy — replace with `// filter exclusion`.

`src/McConnellResult.cpp:267-268` — comment is potentially misleading — replace with `// T2 projection of M` if that is intended.

`src/McConnellResult.cpp:330-331` — says `MCCONNELL_CUTOFF_A`, but code uses config value — replace with `// configured bond cutoff`.

`src/McConnellResult.cpp:346` — “skip if inside the bond” is less exact than the code — replace with `// near-field cutoff`.

`src/MopacMcConnellResult.h:13-16` — explanatory/model prose is longer than needed — replace with `// Wiberg order scales each bond contribution.`

`src/MopacMcConnellResult.cpp:94-96` — “electronically insignificant” is an interpretation, and the code uses a noise floor — replace with `// bond-order floor`.

**6. Correctness**
No confirmed correctness bug from the inlined text alone.

`src/McConnellResult.cpp:189` — possible contract risk: `direction_to_midpoint` receives `kernel.direction`, documented as midpoint-to-atom — check the field contract; if the field means atom-to-midpoint, this is a sign error.

`src/McConnellResult.cpp:267-285` and `src/MopacMcConnellResult.cpp:251-262` — possible feature-contract risk: comments imply dipolar `K` sums, code computes symmetric traceless projection of full `M` totals — check the expected feature definition.

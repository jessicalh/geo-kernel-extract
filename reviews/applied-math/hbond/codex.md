# hbond — codex review (readability focus)

- **Targets:** src/HBondResult.{h,cpp}
- **Model:** codex `gpt-5.5` `xhigh`
- **Date:** 2026-05-24
- **Tokens:** 11,305
- **Brief:** `../_brief.md`

---

**Verdict**

`src/HBondResult.h` mostly tells a coherent story: it states the physical object, source geometry, tensor formula, and exposed API clearly. The main breaks are small: `SampleAt`/`SampleShieldingAt` naming drift and slightly compressed notation in the formula.

`src/HBondResult.cpp` has a recognizable through-line: resolve DSSP H-bonds → filter atom/source evaluations → accumulate tensors → write features. The hard part is not the high-level order; it is local readability. Several passages make the reader decode duplicated donor/acceptor resolution logic, compact tensor algebra, and unlabeled bookkeeping before they can see the physics.

**1. Coherent Story / Readability**

`src/HBondResult.cpp:72` — `ComputeHBondKernel` jumps straight from geometry to a compact tensor loop; as written, a reader must map `d`, `r`, `cos_theta`, `f`, and `M_over_r3` back to the header formula while reading the loop — add short block labels such as `// field direction`, `// angular factor`, `// tensor kernel`.

`src/HBondResult.cpp:89` — `result.f` is computed but never used in this file; as written, a reader must wonder whether this scalar is physically important or leftover diagnostic state — either name it as a diagnostic angular factor or remove/avoid exposing it if unused.

`src/HBondResult.cpp:91` — the tensor expression fuses three meaningful mathematical terms into one assignment — split into named intermediates like `alignment_term`, `hbond_axis_term`, and `radial_term`.

`src/HBondResult.cpp:148` — the two DSSP resolution loops are logically symmetric but long enough that the reader must compare them line-by-line to confirm donor/acceptor roles — add signposts before each loop and use locally explicit names for the resolved endpoints.

`src/HBondResult.cpp:160` — sequence separation filtering is embedded inside resolution bookkeeping and `GeometryChoice` logging — add a terse `// sequence exclusion` label.

`src/HBondResult.cpp:182` — distance validity is a compound condition with two physically different meanings: singular geometry and maximum H-bond distance — add `// source distance gate`, or split into named booleans.

`src/HBondResult.cpp:316` — the kernel is computed before filter acceptance, so rejected pairs still carry a tensor result that is then discarded; this is probably necessary for distance context, but as written the order makes the reader ask why expensive math precedes filtering — add `// geometry for filters`.

`src/HBondResult.cpp:389` — donor/acceptor flagging is separated from H-bond resolution and atom accumulation, requiring another full scan after the main loop — add `// endpoint flags`, or move visually closer to per-atom bookkeeping if possible without changing behavior.

`src/HBondResult.cpp:428` — `PackST_HB` appears abruptly after the sampling function — add `// feature packing`.

**2. Naming Carries Meaning**

`src/HBondResult.h:17` — `h_b` and `h_a` in the formula are less explicit than the surrounding `h_hat` language — use `h_hat_a` / `h_hat_b`.

`src/HBondResult.cpp:44` — `donor_N` and `acceptor_O` are clear, but they hold atom indices, not atoms — consider `donor_N_atom` and `acceptor_O_atom`.

`src/HBondResult.cpp:49` — `h_hat` is compact but domain-readable once introduced — acceptable.

`src/HBondResult.cpp:66` — `M_over_r3` preserves the formula symbol but not the physical result — consider `kernel_tensor_Ainv3` or `shielding_kernel_Ainv3`.

`src/HBondResult.cpp:67` — `f` carries no meaning — if kept, use `angular_factor`.

`src/HBondResult.cpp:79` — `d` is conventional but its frame matters — consider `midpoint_to_atom`.

`src/HBondResult.cpp:80` — `r` is acceptable locally, but `atom_distance` would reduce decoding.

`src/HBondResult.cpp:148` — `ri`, `bi`, `acc_ri`, `don_ri` are compact; a one-pass reader must keep residue roles in memory — consider `residue_idx`, `bond_slot`, `acceptor_residue_idx`, `donor_residue_idx`.

`src/HBondResult.cpp:311` — `count_3_5` hides the configured radius behind a literal-looking name — consider `nearby_hbond_count`.

`src/HBondResult.cpp:428` — `PackST_HB` is terse and abbreviation-heavy — consider `PackHBondSphericalTensor`.

**3. Visible Math Structure**

`src/HBondResult.cpp:84` — distance, unit vector, angular factor, scalar factor, and tensor construction are present but not visually grouped — insert short labels for these stages.

`src/HBondResult.cpp:91` — tensor construction is the worst readability hotspot — expose the three formula terms as named pieces before summing.

`src/HBondResult.cpp:327` — sequence separation is mathematically simple but visually buried in context setup — label it as `// endpoint sequence gap`.

`src/HBondResult.cpp:380` — nearest tensor decomposition repeats the kernel call after nearest tracking — readable enough, but add `// nearest contribution`.

`src/HBondResult.cpp:395` — accumulated Cartesian tensor becomes spherical tensor without a signpost — add `// decompose total tensor`.

**4. Function / Method Naming**

`src/HBondResult.h:41` — `Compute` is standard for result classes and acceptable if this is the local convention.

`src/HBondResult.h:46` — `SampleShieldingAt` is clear about return purpose and spatial sampling.

`src/HBondResult.cpp:72` — `ComputeHBondKernel` is somewhat vague about return units and tensor form — consider `ComputeHBondTensorKernel`.

`src/HBondResult.cpp:428` — `PackST_HB` is not self-explanatory — use `PackHBondSphericalTensor`.

**5. Comments As Signposts**

`src/HBondResult.cpp:26` — large banner comment is accurate but verbose for production code — replace with `// resolved H-bond`.

`src/HBondResult.cpp:55` — large derivation banner is useful, but the framing is heavy — replace the banner with `// dipolar tensor kernel` plus the formula.

`src/HBondResult.cpp:117` — this comment explains an implementation dependency clearly; keep it.

`src/HBondResult.cpp:128` — Step 1 comment is helpful but long — keep the intent, shorten to `// resolve DSSP partners`.

`src/HBondResult.cpp:162` and `src/HBondResult.cpp:183` — `// ---- GeometryChoice: hbond resolution ----` is process-flavored and visually noisy — replace with `// record rejected H-bond`.

`src/HBondResult.cpp:276` — Step 2 comment is coherent but too explanatory — replace with `// evaluation filters`.

`src/HBondResult.cpp:293` — Step 3 comment is useful — could be shortened to `// accumulate atom tensors`.

`src/HBondResult.cpp:351` — “3.5A” comment conflicts with configurable radius naming — use `// nearby H-bond count`.

`src/HBondResult.cpp:418` — clear signpost; keep or shorten to `// near-field cutoff`.

**6. Correctness**

`src/HBondResult.h` — no concrete correctness issue found from the inlined text.

`src/HBondResult.cpp:148` — possible bounds assumption: `dssp.AllResidues()[ri]` assumes DSSP residue count matches `protein.ResidueCount()` — I do not know whether that invariant is guaranteed; check `DsspResult` construction.

`src/HBondResult.cpp:416` — redundant singularity check after `ComputeHBondKernel` already returns zero with `distance == 0` for singular points — not a bug, but it may mislead readers into thinking the second guard changes behavior.

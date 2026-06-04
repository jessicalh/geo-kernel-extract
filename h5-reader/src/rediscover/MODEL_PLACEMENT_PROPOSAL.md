# Model-placement proposal: one model in C++, an e3nn fitter, no Python recreation

Status trued 2026-06-04: the proposal boundary remains useful, but later work
implemented/evolved parts of it. The e3nn protocol fix and the bounded
interpolation result supersede pre-fix protocol details; see
`POSTMORTEM_E3NN_PROTOCOL_FIX_2026-06-04.md` and
`BUILD_INTERP_RESULT_2026-06-04.md`. Positive 1P9J leave-atoms-out/between
claims below are historical only under the true-LOAO retraction.

Status: **proposal for lead + user review. NO code changed in this pass; nothing
merged.** Branch `h5-reader-pysr-spike`. Written by an independent Opus
architect/reviewer, deliberately NOT the agent that authored the offending
`equiv_t2.py` / `lib_T2` end-run, so the judgment is not the offender rationalizing
its own shortcut.

This proposal answers ONE question: **what do we do** so that the model lives in
C++ (the spine), C++ exports only what Python needs, and the equivariant fitter
uses e3nn properly on that export — with the Python physics-recreation retired or
demoted to explicitly-labeled integrity tests, and **no half-fix**. It folds
codex's `ENGINE_TOTALITY_DESIGN.md` (fold/sink mechanics) into one coherent design
and supersedes it on the model-placement / e3nn question.

It enforces `feedback_no_python_physics_except_labeled_integrity_test`,
`feedback_model_is_spine`, and `feedback_functional_api_minimal_clarifying_abstraction`.

---

## 0. What I verified first (so the boundary claim is grounded, not asserted)

The whole proposal rests on one factual claim — **the C++ substrate already emits
everything the equivariant fitter needs, so every Python recompute is a pure
end-run, not a missing-export workaround.** Confirmed against the code:

- **Per-source geometry is already exported** (`RediscoverTypes.h:57-113`,
  `RecordSink.cpp:137-169`): each source row carries `disp_local_{x,y,z}` (Å,
  target→source in the recorded local frame), `r`, `cos_theta`, `dipolar`
  (`(3cos²θ−1)/r³`), and the orientation vectors `source_normal_local_{x,y,z}`
  (ring dipole axis) / `bond_axis_local_{x,y,z}` (McConnell), plus typed identity
  and `ring_intensity` (`QtRing::LiteratureIntensity`). For charge sources:
  `source_q_e`, `disp_local`, `r`.
- **The local frame is exported per record** (`RecordSink.cpp:100-109`):
  `frame_z/x/y`, `frame_valid`, `frame_anchor_atom_index`, `frame_variant`.
- **The target tensor is exported in BOTH frames and BOTH the library T0/T1/T2 and
  raw 3×3 form** (`RecordSink.cpp:114-129`, `RediscoverTypes.h:33-50`): `total_raw`
  (lab 3×3), `total_local` (rotated into the atom's local frame), and the
  library-basis decomposition `total_decomp.{T0,T1,T2}` from `DecomposeLibrary`.
- **The T2 target is ALSO an NPY sidecar, already shaped `(rows, 5)`**
  (`RecordSink.cpp:177-188, 280-303`): `rediscover_<case>_sources_target_T2.npy`,
  `..._sources_target_local_T2.npy`, `..._sources_bare_kernel_T2.npy`, and the
  aggregated equivalents.
- **The library T2 basis is fixed and pinned** (`SphericalBasis.cpp:29-37`):
  component order `[xy, yz, zz, xz, xx−yy]` (m = −2..+2) with isometric
  normalization `[√2, √2, √(3/2), √2, 1/√2]`, a faithful port of
  `nmr::SphericalTensor::Decompose`, fixture-pinned (`√6`, basis check 4.88e-8).
- **`equiv_t2.py`'s `lib_T2` (lines 32-39) is a byte-for-byte reimplementation of
  that exact C++ function** — its own comment says "matches DecomposeLibrary", and
  its own basis check (line 46) confirms `lib_T2(raw) == emitted dft_total_T2` to
  4.9e-8. It recomputes, in numpy, a tensor the producer already wrote to both a
  CSV column and an NPY. That is the end-run in one line.
- **e3nn 0.6.0 and torch 2.11.0+cu130 are installed in system Python**
  (`/usr/bin/python3`), NOT in `analysis/venv` (which is PySR/julia-only:
  pysr 1.5.10, juliacall, no torch, no e3nn). The torch-based offenders run under
  system python today; the e3nn fitter runs there too (with the
  `feedback_cuda_ld_path` `LD_LIBRARY_PATH` note). **Environment-split finding,
  section 7 open question Q1.**
- **`o3.spherical_harmonics("2e", x, normalize, normalization='component')`**
  returns the l=2 (`1x2e`, dim 5) irrep — the equivariant Y2 the fitter needs.
  The library↔e3nn reconciliation is a fixed 5×5 change-of-basis (section 2.3).

**Conclusion of the audit:** there is no export the fitter is missing. The Python
recreation exists for the author's convenience (plain numpy/torch over learning
e3nn's convention), which is exactly the offense the law names. This is not a
close call and there is no defensible half-fix.

---

## 1. The C++/Python boundary and the minimal export

### 1.1 The law, stated as a boundary

**C++ (the spine) owns and computes:** the typed protein model; per-source
geometry in a recorded local frame; the physics kernels and sums (`dipolar`,
`bare_T0`/`bare_kernel`, McConnell sums, charge multipoles); the spherical
decomposition (`DecomposeLibrary`, T0/T1/T2); the DFT target in raw, lab, and
local frames; the local-frame bases. It exports a method-agnostic substrate.

**Python (the fitter) owns and computes:** model fitting only — ridge / OLS on the
exported scalars and aggregates; the **e3nn equivariant network** on the exported
geometry; PySR symbolic distillation of a fitted readout. Python **reads** the
substrate (CSV + NPY sidecars + manifest), assembles design matrices and tensors,
trains, reports. Python builds its OWN equivariant featurization **inside e3nn**
(that Y2 is the network's architecture, computed by the library from the exported
vectors — not a numpy projection of the target).

**Python may NEVER recompute** (each is a rollback offense unless labeled +
separately-sourced per section 3):

- the spherical projection / `DecomposeLibrary` (`lib_T2`) — it is exported.
- the Pople ring-current kernel `(3cos²θ−1)/r³` or its intensity-weighted sum
  `Σ intensity·dipolar` — `dipolar` per source and `bare_T0` per atom are exported.
- the charge electric field `Σ q·d/r³` (`look_charge_dipole.py:129-134`) — the
  charge multipole is a C++ relationship (`charge_dipole`, `charge_quadrupole`);
  the field/EFG come from APBS (`buckingham_efield`, `efg`), exported.
- any kernel, frame, normal, or sum the substrate already carries.

### 1.2 The minimal export set for the two fitters

The export is **already almost exactly this**; the change is to make the
equivariant tensor inputs first-class NPY sidecars (not only CSV columns), so the
e3nn consumer never parses physics out of CSV text and never needs a numpy
projection. Per relationship the substrate emits:

**Scalar / ridge fitter (the aggregated row — already correct):**
`<case>_aggregated.csv` — identity, frame, `sum_dipolar_all`,
`sum_dipolar_producer_valid`, per-type sums, `cutoff_A`, the bare-kernel
cross-check, and the full DFT target (raw 3×3 + T0/T1/T2 + `total_local`).

**Equivariant fitter (the per-source set — already emitted; promote to NPY):**
the inputs an e3nn net consumes, keyed by `(row_id, source_row_id)`:

| Export | Current location | Proposed minimal form |
|---|---|---|
| `disp_local` (Å, vec3) | `_sources.csv` cols | per-source NPY `(n_src, 3)` |
| `source_normal_local` / `bond_axis_local` (vec3) | `_sources.csv` cols | per-source NPY `(n_src, 3)` |
| `r`, `cos_theta`, `dipolar`, `ring_intensity` (invariant scalars) | `_sources.csv` cols | per-source NPY `(n_src, k)` |
| `is_self_or_bonded`, type key, `ring_index`/`bond_index` | `_sources.csv` cols | per-source NPY (int) |
| group key `(atom_index, h5_row)` → `row_id` | derivable | per-source NPY `source_row_id → row_id` |
| **target T2 (library basis)** `(rows, 5)` | **already NPY** | keep `rediscover_<case>_*_target_T2.npy` |
| target `total_local` 3×3, T0, T1 | CSV + (T2 NPY) | add `(rows, 9)` local-tensor NPY |

This is the same content as today; the only structural change is **per-source
arrays move from CSV columns to columnar NPY** (which the fold/sink totality
design already wants for scale, section 4) so the fitter does
`np.load(...)` rather than `pd.read_csv(...).to_numpy().reshape(...)`. No new
physics is computed in C++ for this; it is a carrier change, not a model change.

**The minimality test (per `feedback_functional_api_minimal_clarifying_abstraction`):**
the equivariant net needs, per source, (i) the displacement vector, (ii) the
source orientation vector, (iii) invariant scalars for the radial MLP, (iv) a
group index to pool by, and per group the target tensor in a known basis/frame.
That is the whole list. Nothing in `equiv_t2.py` is computed from anything outside
it — which is the proof the recompute is gratuitous.

---

## 2. The equivariant model on e3nn

### 2.1 What replaces `equiv_t2.py`

`equiv_t2.py` is a hand-rolled 3-channel l=2 ansatz: per source it forms
`Y2(r̂)`, `Y2(n̂)`, `Y2sym(r̂,n̂)` **via the numpy `lib_T2`**, weights each by a
learned invariant radial MLP `R(r, intensity)`, and sum-pools (`index_add_`) to a
per-group T2, de-meaned per atom, against the de-meaned library-basis DFT T2. It
reached T2 R²=0.467, |T2| r=0.756 (FINDINGS.md, STATE.md).

The replacement keeps that physics ansatz exactly but builds it **inside e3nn**:

1. **Inputs** come from the C++ NPY export (section 1.2): per source `disp_local`
   (vec3), `source_normal_local` (vec3), invariant scalars `[r, ring_intensity]`,
   and `source_row_id → row_id`.
2. **Equivariant features** are e3nn spherical harmonics, computed by the library:
   `Y2_r = o3.spherical_harmonics("2e", disp_local, normalize=True, normalization="component")`
   and `Y2_n = o3.spherical_harmonics("2e", source_normal_local, ...)`. The mixed
   `Y2sym(r̂,n̂)` channel is an `o3.FullyConnectedTensorProduct(1o, 1o → 2e)` on the
   two unit vectors (the principled e3nn way to build the rank-2 cross term), or —
   if we want to match the ansatz exactly — `o3.TensorProduct` of `1o ⊗ 1o → 2e`.
3. **Radial weights** are an invariant MLP on `[r, intensity]` (and optionally an
   e3nn `soft_one_hot_linspace` radial embedding), producing the per-channel
   scalars that gate each `2e` feature — i.e. the same `R_A, R_B, R_C` of the
   ansatz, now as e3nn `weight` inputs to the tensor products.
4. **Permutation-invariant pooling** is `scatter`/`index_add` of the per-source
   `2e` contributions into per-group `2e`, which is the Deep-Sets sum the substrate
   was designed for and is equivariant because it sums equal-irrep features.
5. **Target** is the C++-exported `target_local_T2` NPY (already library basis),
   **after the change-of-basis to e3nn's `2e` convention** (section 2.3),
   de-meaned per atom (the same baseline strip).
6. **Loss / honesty** unchanged: MSE on the 5-vector; frame-split AND
   leave-atoms-out; report `|T2|` (rotation invariant) and per-component R².

This is strictly more capable than the hand-roll (e3nn gives correct, tested
Wigner-D equivariance, gated nonlinearities if we want depth, and the tensor
product for the cross term) and removes every line of numpy spherical-harmonic
math. The 3-channel ansatz is recoverable as the minimal e3nn instance, so we can
reproduce R²=0.467 as a parity check before extending.

### 2.2 Why e3nn and not the hand-roll (the line held)

The hand-roll's only defensible reason to have existed is a feasibility spike
("does the T2 carry signal" → R²=0.44). That justifies a throwaway prototype, not
a persisted unlabeled model that end-runs the declared dependency
(`sumpool_t0.py` literally calls itself "the precursor to the e3nn T2 model"). The
basis-convention excuse is thin: the library basis is real-SH / e3nn-style, the
reconciliation is one fixed 5×5 matrix, and e3nn exists precisely to get
equivariance right. **There is no "doing it our way is the only way" here.**

### 2.3 The library-basis ↔ e3nn irreps change-of-basis, and how to test it

This is the one genuinely technical seam and must be handled as a labeled,
fixture-pinned transform — NOT by re-deriving the projection in Python.

- The C++ target is in the **library T2 basis**: order `[xy, yz, zz, xz, xx−yy]`,
  normalization `[√2, √2, √(3/2), √2, 1/√2]`, Frobenius-norm preserving
  (`SphericalBasis.cpp`).
- e3nn's `2e` irrep is real SH in e3nn's component order and `normalization`
  convention, **with e3nn's internal axis convention (y, z, x)** — a known gotcha
  that will silently rotate everything if ignored.
- The map between them is a single **constant orthogonal 5×5 matrix** `C` such
  that `target_e3nn = C @ target_library`. It does not depend on data; it is
  determined once by the two conventions.

**How to derive `C` without recomputing physics:** feed a small set of *known*
symmetric traceless 3×3 tensors (the 5 basis tensors, plus a random one) through
(a) the C++ `DecomposeLibrary` path — already available as the exported
`target_T2` for a controlled fixture — and (b) e3nn's `o3` Cartesian-to-irrep map
(`o3.Irreps("2e")` with e3nn's documented rank-2 reduction, or
`spherical_harmonics` on chosen directions). `C` is read off by least squares on
those known pairs. This `C` lives in a small, **explicitly labeled** module
(`change_of_basis.py`, "library T2 basis ↔ e3nn 2e; derived from convention, not
a physics recompute") with a test, NOT inside the model forward pass as a
re-projection.

**The test (two independent checks):**
1. **Orthogonality / round-trip:** `C` is orthogonal (`Cᵀ C = I` up to the
   normalization), and `C⁻¹ @ (C @ t) == t` for random `t` to ~1e-10.
2. **Equivariance round-trip against C++:** take a random rotation `R`; rotate a
   raw 3×3 tensor in C++ (or rotate `disp_local`/`source_normal_local` and let the
   fixture re-extract), confirm the e3nn `2e` of the rotated tensor equals
   `D²(R) @ (e3nn 2e of the unrotated tensor)` where `D²` is e3nn's Wigner-D. This
   proves the C++ library tensor, transported through `C`, transforms under e3nn's
   `D²` — i.e. the basis match is correct AND equivariant. This is the rigorous
   replacement for `equiv_t2.py`'s one-line `np.nanmax(|lib_T2(raw) − emitted|)`
   check, and it is justified by an independent source (C++ vs e3nn, two
   implementations) so it is a *legitimate* labeled integrity test (section 3).

Note this also subsumes the existing `equiv_t2.py:46` basis sanity check —
demoted from "load-bearing line inside the model" to "one assertion in the
change-of-basis test."

### 2.4 Alternative equivariant libraries, in parallel (with trade-offs)

e3nn is the declared dependency and the recommendation. Alternatives, so the
choice is on the record:

| Library | Pros | Cons / trade-offs | Verdict |
|---|---|---|---|
| **e3nn (0.6.0, installed)** | Declared dep; mature `o3` irreps + Wigner-D + tensor products; exact match to the l=2 ansatz; PyTorch; CUDA present | Convention learning curve (the y-z-x axis gotcha, component order) — handled once by section 2.3 | **Recommend.** The line the brief draws. |
| `e3nn-jax` | Same algebra in JAX; faster autodiff/JIT for big nets | New runtime (JAX/CUDA stack); not installed; project is torch elsewhere | Reject for now — runtime churn, no gain at this scale (~7 coupled atoms). |
| MACE / NequIP / Allegro | Batteries-included equivariant *interatomic* models | Built for energies/forces, heavy, opinionated featurization; we want a thin per-source l=2 sum-pool, not an MPNN; would re-introduce a framework | Reject — over-framework for a 3-channel ansatz; `feedback_no_abstractions`. |
| Hand-roll (status quo) | Zero deps | The offense; unverified equivariance; numpy `lib_T2` end-run | **Reject — this is what we are removing.** |

If e3nn's convention handling ever proves painful, `e3nn-jax` is the fallback with
the *same* irrep semantics, so the change-of-basis test (2.3) ports unchanged.

---

## 3. Disposition of every offending script

Principle: a Python physics recompute survives ONLY if it is (a) explicitly
labeled as an integrity test AND (b) justified by an *independent* data source it
cross-checks against. Otherwise it is retired or rebuilt on the boundary.

| Script | What it recomputes | Disposition | Justifying independent source (if kept) |
|---|---|---|---|
| `equiv_t2.py` | `lib_T2` (== `DecomposeLibrary`) + hand-rolled Y2; IS the model | **Rebuilt on e3nn** (section 2). The model moves to `equiv_t2_e3nn.py` consuming the NPY export; `lib_T2` deleted from the model path. The basis check becomes the labeled `change_of_basis` test. | n/a (model, not a test) |
| `sumpool_t0.py` | Pople kernel readout `(3cos²−1)/r³` for "did the equation fall out" | **Rebuilt as the scalar precursor on the boundary** *or retire*: the sum-pool over exported `dipolar`/invariants is fine (that is fitting, not recompute); the *readout-vs-Pople* comparison is the recompute. Keep the readout comparison ONLY as a labeled "kernel-form integrity check," justified because the producer's `bare_T0` (Giessner-Prettre) is an *independent* analytic source the learned `g` is compared against. Label it; do not headline (it is the circular half per FINDINGS.md). | producer `bare_T0` analytic kernel (independent of the NN's input geometry path) |
| `sumpool_kernel.py` | Pople kernel + 1/r³ + angular readout comparisons | **Same as above** — the pooling fit stays; the explicit `(3cos²−1)/r³` / `1/r³` comparison arrays are demoted to a labeled "form-recovery integrity check" against `bare_T0`. The npz it writes (`sumpool_kernel_readout.npz`) is the input to `pysr_distill.py`, so it must keep emitting the learned `g` (allowed: that is the fitted readout) — but stop recomputing `pople_geo`/`pople_int` as model inputs. | producer `bare_T0` |
| `refine_kernel.py` | far-field fit of `g ~ A·intensity·(3cos²−1)/r³` + figure | **Keep as labeled coefficient-extraction/figure**, reading the learned `g` readout; the `(3cos²−1)/r³` it fits against is the analytic comparator (the literature-coefficient check FINDINGS.md still wants), so it is a legitimate labeled cross-check, not a parallel model. Rename intent in the header to "integrity/coefficient check." | literature Pople form (analytic, independent) |
| `look03_coefficient.py` | `K_int = Σ intensity·dipolar` in Python | **Rebuilt / retired.** The intensity-weighted sum is a C++ reducer output (`sum_dipolar_producer_valid` already; add a per-type intensity-weighted aggregate if needed). Python should read the exported aggregate, not re-sum `ts.ring_intensity * ts.dipolar`. The within-atom corr-vs-DFT analysis stays (that is analysis on exported columns). | n/a — recompute removed; analysis kept on exported aggregate |
| `look_charge_dipole.py` | `E = Σ q·d/r³` field in Python (lines 129-134) | **Rebuilt / retired.** The charge field/dipole is the C++ `charge_dipole` relationship (μ aggregated row) and APBS `buckingham_efield`. Python must read the exported `mu`/field, not recompute `q*d/r³`. The `mu_*` columns it also reads ARE exported (`WriteChargeDipoleAggregatedRow`) — keep that path; delete the hand-rolled field block. | n/a — recompute removed; the within-atom fit on exported `mu`/field stays |

**Non-offenders (read on the boundary, keep as-is):** `look01_ring_triangulate.py`,
`look02_self_vs_throughspace.py`, `credibility_check.py`,
`credibility2_instantaneous.py`, `diag_differencing.py`, `oracle_parity.py`,
`sumpool_mcconnell.py`, `pysr_distill.py` — these fit / analyze / distill from
exported columns and the learned readout; verify each only consumes exported data
and the fitted `g` (not a recomputed kernel) during migration. `pysr_distill.py`
distills the *learned* `g`, which is legitimate; its only dependency to clean is
that its input npz no longer carries recomputed Pople arrays as if they were data.

**The labeling standard (so "labeled integrity test" is not a loophole):** a kept
recompute must (1) name, in its module docstring, that it is an integrity/form
check and NOT a model input; (2) name the independent data source it validates
against (here always the producer's analytic `bare_T0` or the literature Pople
form, which are independent of the NN's geometry→shielding path); (3) feed nothing
into the model's prediction path. If any of the three fails, it is retired.

---

## 4. How this folds into and supersedes the fold/sink totality design

codex's `ENGINE_TOTALITY_DESIGN.md` correctly diagnoses the engine duplication
(`RunRelationship` hardwires `WriteSourceRows`/`WriteAggregatedRow`, forcing a
sibling `RunBroadBackbone`) and proposes the right minimal fix: parameterize the
loop tail by a typed `Reducer` + `Sink` over a closed set of record shapes
(`ScalarSourceSumRecord`, `Vec3SourceSumRecord`, `T2SourceSumRecord`,
`Vec3FeatureRecord`, `T2FeatureRecord`, etc.), with a per-source emission policy
(`none` / `compact_npy` / `debug_csv`) for scale. **That part stands and I adopt
it.** codex is good at the fold/sink mechanics and the enumeration is sound
(`RelationshipEngine.cpp:77-79` is exactly the hardwired tail it describes).

What that design did NOT carry — and where this proposal supersedes it — is the
**consumer side and the model-placement law**. The two become ONE design by adding:

1. **The `compact_npy` source payload IS the e3nn fitter's input contract.**
   codex's `compact_npy` was justified only by *scale* (charge rows hit
   1.7–5.5 GB). Reframe it: `compact_npy` columnar per-source arrays
   (`disp_local`, orientation vec3, invariant scalars, group index) are the
   *canonical equivariant-fitter inputs* (section 1.2). So the same carrier change
   serves both scale and the boundary — one decision, not two. The `T2FeatureRecord`
   / `T2SourceSumRecord` sidecars are the target/feature tensors e3nn consumes.
   This kills the residual reason Python reached for `lib_T2`: the tensor is on
   disk in a load-ready `(n, 5)` / `(n, 3)` NPY, never parsed from CSV.

2. **A consumer-side boundary rule added to the design's "deliberately cut" list.**
   codex listed what the C++ engine must not become (no registry, no superset
   carrier, no DSL). Add the symmetric Python rule: **the consumer must not
   recompute any kernel/projection/field the engine emits; the equivariant fitter
   is e3nn on the exported geometry.** This makes the totality design two-sided —
   the same `feedback_functional_api_minimal_clarifying_abstraction` razor applied
   to both the C++ fold/sink AND the Python fitter, designed across the whole known
   shape set at once (the 9 relationships × {ridge, e3nn} fitters), not patched per
   script.

3. **Resolve codex's open Q3/Q4 (L1 target, all-atom frame) consistently with the
   equivariant path.** Because the equivariant fit is tumbling-safe in the lab
   frame (SURFACE_DESIGN "T2 frames RESOLVED"), the e3nn consumer can take the
   lab-frame `target_T2` directly for field/charge items; the per-atom local frame
   stays only where it already exists (ring/HN). This is a consumer requirement
   that informs which target NPY the engine emits per relationship — folding the
   two designs' frame decisions together.

Net: **one design.** C++ side = codex's parameterized fold/sink with the typed
record/sink set and source-emission policy. Python side = read the
`aggregated`/`compact_npy`/T2-NPY export; ridge on aggregates; e3nn on per-source
geometry; no recompute. The `compact_npy` carrier is the hinge that makes both
true at once.

---

## 5. Migration and re-validation plan

Phased, each phase gated, byte-parity and broad output preserved. No C++ physics
change; the C++ work is carrier/fold-sink only (codex's plan), the Python work is
the e3nn fitter + script disposition.

**Phase 0 — review.** Lead + user sign off on this boundary + the e3nn plan + the
script dispositions. No code until then. (This file.)

**Phase 1 — e3nn fitter stands up against the CURRENT export (no C++ change yet).**
Build `change_of_basis.py` + its test (section 2.3) and `equiv_t2_e3nn.py`
consuming today's `_sources.csv` + `rediscover_<case>_*_target_local_T2.npy`.
Gate: **reproduce or beat** the hand-roll baseline — T2 R² ≈ 0.467 (≥, target
improve), |T2| r ≈ 0.756, frame-split AND leave-atoms-out, on the canonical
substrate (`/tmp/rediscover-out-v2/`). This proves the e3nn model on the existing
boundary before touching C++. **`equiv_t2.py` is retired the moment
`equiv_t2_e3nn.py` matches it.**

**Phase 2 — script dispositions (Python only).** Apply section 3: rebuild
`look03` / `look_charge_dipole` to read exported aggregates (delete the
`Σ intensity·dipolar` and `Σ q·d/r³` recomputes); relabel the kept integrity
checks (`sumpool_*`, `refine_kernel`) per the labeling standard; verify the
non-offenders only touch exported data. Gate: no script imports a recomputed
kernel into a model path; `pysr_distill` still distills the learned `g`.

**Phase 3 — C++ fold/sink totality (codex's plan), byte-parity preserved.**
Implement the parameterized loop tail + typed records + sink-policy. Gate
(unchanged from codex + the existing oracle): **ring/mc composed-vs-procedural
byte parity exact** (`oracle_parity.py`, the 4 CSVs 141000/20500/812205/26000 +
12 sidecar NPYs identical), broad-backbone counts/schema stable, 0% invalid
frames, per-atom-type σ_iso within correlate-not-match tolerance. Ring/mc keep
their legacy CSV during this phase (do not change byte-parity output in the same
step as the engine).

**Phase 4 — promote per-source arrays to `compact_npy` and repoint the e3nn
consumer.** Emit the columnar per-source NPY (section 1.2 / section 4) under the
sink policy; switch `equiv_t2_e3nn.py` to `np.load` the NPY instead of CSV. Gate:
the e3nn fit result is unchanged (same R²) reading NPY vs CSV — a pure carrier
swap. This is also where charge-source CSV stops being multi-GB.

**Phase 5 — extend e3nn to the new T2 relationships** (`efg`, `charge_quadrupole`,
field items) once their data flows, lab-frame-equivariant per SURFACE_DESIGN, using
the *same* e3nn model and change-of-basis. Gate: fail-loud on absent arrays;
per-atom rows aligned to DFT order; flagged-unverified T1 carried not dropped.

**Re-validation invariants held throughout:** ring/mc byte-parity (Phase 3 gate);
broad output stable (Phase 3 gate); equivariant result reproduced-or-improved on
e3nn (Phase 1 gate, re-checked Phase 4); reader owns H5, GUI untouched, library
not linked, never merged (branch discipline).

---

## 6. Half-fixes I am explicitly NOT taking, and the full-correct path instead

1. **NOT** "keep `lib_T2` but add a comment that it matches `DecomposeLibrary`."
   That re-labels the end-run without removing it. **Full fix:** the tensor is read
   from the exported NPY; `lib_T2` is deleted from the model path; the one
   surviving projection comparison is the labeled, two-implementation
   change-of-basis test.

2. **NOT** "wrap the hand-rolled Y2 in a function called `e3nn_compat` and call it
   done." Hand-rolled equivariant math is off the table regardless of naming.
   **Full fix:** `o3.spherical_harmonics` / `o3.TensorProduct` compute the
   features; equivariance is the library's, tested by Wigner-D round-trip.

3. **NOT** "the existing Python approach reaches R²=0.44, so it's fine / good
   enough." The bar is the project's law (one model in C++, fitter on a real
   library), not the metric. **Full fix:** e3nn must reproduce-or-beat 0.467 as the
   *gate*, but the reason to switch is correctness of equivariance + no end-run, not
   the number. (This is the exact rationalization the brief told me to stop on; I
   am stopping on it.)

4. **NOT** "only fix `equiv_t2.py`; leave `look_charge_dipole` / `look03`
   recomputing fields and kernel sums because they're 'just analysis.'" A field
   recompute is a kernel recompute. **Full fix:** every model-path recompute is
   removed; only labeled, independently-sourced integrity checks survive.

5. **NOT** "do the e3nn fitter but skip the change-of-basis test and trust the
   conventions line up." The y-z-x axis + component-order gotcha silently rotates
   everything. **Full fix:** the change-of-basis is derived once, fixture-pinned,
   and Wigner-D round-tripped against the C++ tensor.

6. **NOT** "build a new equivariant runner / framework in C++ to emit tensors."
   The tensors are already emitted; e3nn is the fitter. **Full fix:** the only C++
   change is codex's fold/sink carrier (CSV→`compact_npy`), no new physics, no new
   model.

---

## 7. Open questions for the lead / user

1. **Environment.** e3nn 0.6.0 + torch 2.11 are in **system python**
   (`/usr/bin/python3`), NOT `analysis/venv` (PySR/julia-only). Should the e3nn
   fitter (a) run under system python with the `feedback_cuda_ld_path`
   `LD_LIBRARY_PATH`, or (b) get its own torch venv pinned in a requirements file
   so the fitter env is reproducible and gitignored like the PySR venv? I lean (b)
   for reproducibility; flagging because the torch-vs-PySR split is currently
   implicit and undocumented.

2. **`compact_npy` schema = SDK catalog entry.** SURFACE_DESIGN says wide/array
   payloads get one `ArraySpec`-per-NPY in `python/nmr_extract/_catalog.py`. The
   per-source `compact_npy` arrays are new outputs — do they get catalog entries
   now (consistent with the discipline), or does rediscover keep its NPYs outside
   the producer SDK catalog since it is an experimental branch? (The current T2
   sidecars are written by `RecordSink` but I did not find them in
   `_catalog.py` — confirming whether rediscover NPYs are in-scope for the SDK
   catalog at all.)

3. **`sumpool_t0`/`sumpool_kernel` kept-vs-retired.** I propose keeping the
   pooling FITS and demoting only the explicit Pople-comparison arrays to labeled
   integrity checks (independent source = `bare_T0`). Alternative: retire the
   readout-vs-Pople comparison entirely and rely on PySR distillation of the
   learned `g` for the "form fell out" claim. Which does the user prefer — keep the
   labeled comparison, or lean fully on PySR?

4. **Mixed `Y2sym(r̂,n̂)` channel in e3nn.** Match the hand-roll exactly via a
   raw `1o ⊗ 1o → 2e` tensor product (parity check at R²=0.467), or use
   `FullyConnectedTensorProduct` with learnable mixing (more capable, may improve
   the number but diverges from the literal ansatz)? I propose: parity first, then
   the richer product as the "improve" step.

5. **Leave-atoms-out for T2.** FINDINGS.md notes the hand-roll only did the
   frame-split for T2; leave-atoms-out is the sharper test still TODO. Should the
   Phase-1 gate require leave-atoms-out T2 (harder, more honest) or accept the
   frame-split parity as the gate and add leave-atoms-out as a follow-on? Given the
   thinness (~7 coupled aromatic H), the user's call on whether the atom-split is
   statistically meaningful here.

---

## Summary for the lead

**Proposed boundary:** C++ stays the spine and already exports everything the
fitter needs — per-source `disp_local` / `source_normal_local` / `bond_axis_local`
/ invariant scalars / identity, and the DFT target in raw + local frames + library
T0/T1/T2 (CSV *and* `(rows,5)` NPY sidecars). Python only fits; it may never
recompute the projection, the Pople kernel, or the charge field — all are exported.
The only C++ change is the carrier (codex's fold/sink CSV→`compact_npy`), no new
physics. **e3nn plan:** rebuild `equiv_t2.py` as `o3.spherical_harmonics("2e",…)` +
tensor-product cross term + invariant radial MLP + scatter-pool, consuming the NPY
export against the library-basis T2 target via a fixture-pinned, Wigner-D-tested
5×5 change-of-basis (handling e3nn's y-z-x axis gotcha); gate on reproducing/beating
the recorded R²=0.467 / |T2| r=0.756; e3nn-jax is the only credible alternative
(rejected now for runtime churn), MACE/NequIP over-framework. **Script
dispositions:** `equiv_t2.py` rebuilt on e3nn (its `lib_T2` deleted from the model
path); `look03_coefficient.py` and `look_charge_dipole.py` rebuilt to read exported
aggregates (recomputes deleted); `sumpool_*`/`refine_kernel` keep their pooling fits
but their Pople-comparison arrays are demoted to explicitly-labeled integrity checks
justified by the independent analytic `bare_T0`; PySR distillation and the
credibility/look scripts stay as boundary consumers. **Top open questions:** which
Python env runs the e3nn fitter (system vs a pinned torch venv); whether rediscover
NPYs enter the SDK catalog; keep-vs-retire the labeled Pople comparison; exact vs
learnable e3nn cross-term; and whether the Phase-1 gate requires leave-atoms-out T2.
Full proposal: `src/rediscover/MODEL_PLACEMENT_PROPOSAL.md`. No code changed; nothing
merged.

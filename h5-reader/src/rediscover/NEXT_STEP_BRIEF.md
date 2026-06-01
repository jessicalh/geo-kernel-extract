> **SUPERSEDED (2026-06-01) — historical.** This described building the
> multi-scenario surface, which is DONE (commit `99cdc85`). Current handoff:
> `NEXT_SESSION_PROMPT.md` + `REDISCOVERY_MAP.md` + `STATE.md`. Ignore the below.

# Next-step brief — build the multi-scenario rediscover surface

Branch `h5-reader-pysr-spike` (experimental, one-shot, no integration target).
Build/verify in the lead session. Reader owns H5; additive to the reader; GUI
untouched.

## State — where we are

- The one-off ring-current extraction **works and is the ORACLE** the rebuild
  must reproduce: leave-atoms-out universal k ≈ 21 ppm·Å³, within-atom R² = 0.62
  (coupled); equivariant T2 modulation R² = 0.44, |T2| r = 0.75; lib-basis check
  4.9e-8; DFT↔H5 frame rotation ~1e-4° (aligned). McConnell scalar form also
  recovered.
- The general surface is **designed on paper** in `SURFACE_DESIGN.md`: one
  resident immutable C++ body loaded once; a data catalog (one entry per array,
  uniform `value()`); day-1 indexes (per-cloud per-frame KD-trees; the typed-atom
  `(residue,locant)->atom` index; own-ring/own-atom sets); primitive verbs
  (`at/window/neighbors/value`); relationship combinators as **iterated, curried
  closures**; a pure engine loop; a **batch-file** interface (CLI + JSON spec +
  atomic outputs + `manifest.json` + exit code; no UDP/socket control plane).
- A Phase-1 **stub draft** exists in `STUB_LANGUAGE.md` with a surfaced-issues
  list (C1–3 charges; E/T per-atom features; K1/2 cross-check kernels; I1/2
  indexes; L1/F2 frames; X3/4 schema+convention; D1 scenario-spec).
- The **supported set (9 items)** and the **picks** are settled in
  `SURFACE_DESIGN.md` ("The supported set"): charge_source is a parameter (not
  separate items); embeddings IN; polarizability KEEP; **no deferral — there is
  only today** (full stub + worked example per item; issues are findings, not
  punts).
- The **T2-frame question for the new items is RESOLVED** (`SURFACE_DESIGN.md`
  "T2 frames for the new items"): field items fit equivariantly in the lab frame
  (tumbling-safe, no per-atom frame); source-sum kernels via
  equivariant-in-lab-frame or a generic bond-local frame.

## References — read in this order
1. `src/rediscover/SURFACE_DESIGN.md` — THE design (language, set, picks, frame
   resolution, interface, indexing posits, iterated-closure shape).
2. `src/rediscover/STUB_LANGUAGE.md` — the Phase-1 stub + the surfaced issues to
   resolve.
3. `src/rediscover/{DESIGN.md, STATE.md, PLAN.md}` and
   `src/rediscover/analysis/FINDINGS.md` (note: FINDINGS lives under `analysis/`,
   not the rediscover root) — the one-off,
   the oracle numbers, the circular-vs-real credibility audit, the instantaneous
   framing (the trajectory is a geometry sampler, not a process).
4. `spec/substrate_conventions_2026-05-30.md` — SH basis (library order,
   `DecomposeLibrary`), the `(3cos²θ−1)` form, local frames per atom class,
   `charge_source` as a required enum (no default), cutoffs required+recorded,
   T2 preservation.
5. **The OBJECT MODEL** (know what is in it): `src/model/` headers — `QtProtein.h`,
   `Conformation.h`/`TrajectoryConformation.h`, `QtAtom.h` (+ typed predicates),
   `QtRing.h` (+ virtuals), `QtResidue.h`, `QtBond.h`, `QtTopology.h`, `Types.h`
   (`Mat3`/`SphericalTensor`), `DftShielding.h`. (There is NO `OBJECT_MODEL.md` in
   this tree — the headers ARE the object model; read them.)
6. **The NPYs / H5 datasets AS READER OBJECTS** (not raw array names): positions
   via `Conformation::atomPosition`; bs/mc kernels via the shielding time-series;
   `ring_nbhd` via the ring-neighbourhood time-series; `apbs_efg`/`apbs_efield`,
   `aimnet2_charge`/`_charge_response_gradient`/`_embedding`; DFT via
   `DftShieldingFrame`. The catalog (Layer 0) maps each to a reader object.
7. The one-off code `src/rediscover/*.{h,cpp}` — the grounding for the verbs.
8. The analysis scripts `src/rediscover/analysis/*.py` — the oracle + the fitters
   (`sumpool_kernel`, `equiv_t2`, `credibility2`).
9. `h5-reader/CLAUDE.md` (rules); the qt6-cpp skill (treat as a Qt **command-line**
   program — QCoreApplication, no widgets/VTK).

## The maths, both levels — and why this model

- **Goal:** rediscover whether CLASSICAL geometric kernels carry the DFT shielding
  signal, per relationship. **Correlate, not match** (R² shows the signal was
  captured; never pointwise match). Physics explanation, not prediction.
- **Scalar (T0, σ_iso):** additive kernels — ring current Pople
  `intensity·(3cos²θ−1)/r³`; McConnell bond anisotropy `(3cos²θ−1)/r³` about the
  bond axis (Δχ in the parameter); charge dipole/quadrupole multipoles; APBS
  field/EFG; Larsen H-bond. Fit by sum-pooling / SR. The one-off recovered the
  Pople form + a literature-plausible coefficient (the oracle).
- **Tensor (T2 — the angular residual, THE thesis, "T2 is sacred"):** the l=2
  part. For source-sum kernels the per-source contribution is an l=2 tensor
  (`Y2` of the source direction / dipole axis), summed — equivariant
  sum-pooling. For field items, the field-T2 vs DFT-T2 directly (equivariant,
  lab-frame, tumbling-safe). The library T2 basis (`DecomposeLibrary`, isometric
  real-SH, e3nn-style) MUST match the DFT-T2 (verified 4.9e-8).
- **Why the iterated-closure model (follow it, do NOT second-guess it):**
  `σ = Σ_sources f(geometry)` is a pure, permutation-invariant sum over a
  variable source set on an immutable resident **indexed** trajectory. Curried
  verbs (capture body+config) + `map` over the `(atom,frame)` index space + an
  inner `map` over sources + a state-carrying `fold` express this generically;
  the one-off is exactly this hand-inlined and un-curried; it proved out on ring
  current and generalizes to all 9 items. This is the right shape for our maths
  on a fully-loaded, in-memory-indexed trajectory — extend it, follow it out.

## The task

**Scope of THIS pass (explicit — it is not "all nine runnable," and not "stubs
only"):** a **runnable engine pass for `ring_current` + `mcconnell`** (build +
validate `ring_current` against the oracle), on the new spine (Catalog,
ResidentIndexes, Relationship, engine, the decided output carrier). The other
seven are **stubbed compile-and-fail-loud** (they refuse rather than emit zeros)
with their carrier/charge/frame work flagged. The DESIGN stays complete (no item
dropped); the BUILD stages on data + decision readiness. That is honest staging,
not deferral. **Agreed staging order** (both dry-runs + lead): spine (Catalog /
ResidentIndexes / Relationship / engine / carrier) → **ring_current (oracle
gate)** → mcconnell → per-atom-feature items (efg, buckingham_efield, CRG,
embedding — proves the carrier) → charge multipoles (AIMNet2 first; FF14SB/MOPAC
as data allows) → **larsen_hbond last** (most open decisions).

Resolve the `STUB_LANGUAGE.md` issues into a settled surface, keeping the 9-item
set, the picks, and no design-deferral:
- **C1–3 charges:** `charge_source` parameter; AIMNet2 per-frame (real now), MOPAC
  per-frame (writing now, lands AM), FF14SB (ADD a prmtop partial-charge read to
  the reader, additive). Fail-loud on an absent source; never auto-fallback.
- **E/T per-atom features:** make per-atom-feature a **first-class relationship
  shape** (not a degenerate `self()`); **per-relationship schemas** (drop the
  superset — X3); the record/output must carry a 256-d vector (embedding) and a
  5-component T2 (efg).
- **K1/2 cross-check:** `cross_check_kernel` becomes OPTIONAL (only ring/mc/EFG
  have one); resolve which Larsen TS is the H-bond kernel.
- **I1/2 indexes:** specify + stub the typed-atom `(residue,locant)->atom` index
  and the per-cloud day-1 KD-trees (atoms / bond-midpoints / ring-centers /
  charge-sites).
- **L1/F2 frames:** apply the resolution in `SURFACE_DESIGN.md`; address the
  latent `Types.h` T1 basis/parity inconsistency for any l=1→T1 (efield) path.
- **X4 convention:** fix the manifest example's `scipy_real`-vs-library-order
  inconsistency to match `DecomposeLibrary`.
- **input CLI (was D1 "scenario-spec"):** RESOLVED 2026-05-31 — keep the existing
  flags (`--run` / `--out` / `--case`, + optional `--charge-source` for the
  multipole cases). NO JSON spec, no config framework. Per-relationship "validation"
  is just fail-loud-on-absent-data (existing behavior), not a new ValidateScenario
  system beyond that. (Input is trivial; the per-relationship *output* schema is
  the carrier, a separate concern.)
Validate against the ORACLE: the ported ring-current relationship must reproduce
the numbers above.

## Discipline (non-negotiable)
Iterated/curried closures; named bundles in code (NOT a plugin kit; ABC
polymorphism is fine); reader owns H5 (Python consumes the emitted substrate,
never opens `trajectory.h5`); batch-file interface (no socket control plane);
additive to the reader; GUI untouched; build/verify in the lead session; the
oracle is the faithful-rebuild gate.

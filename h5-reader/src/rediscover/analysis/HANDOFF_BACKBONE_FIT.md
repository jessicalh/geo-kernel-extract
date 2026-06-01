# Handoff — broad-backbone equivariant T2 fitter (codex: run + critique)

Branch `h5-reader-pysr-spike`. NEVER MERGE. This is the codex half of the
author-with-Claude / run-with-codex split (the Claude Agent-tool sandbox denies
`import torch`/`import numpy`, so this pass AUTHORED + self-reviewed + grep-proved
the new fitter; codex runs it and critiques). The C++ model was NOT touched — only
`analysis/equiv_t2_backbone_e3nn.py` was added.

## What was authored

`analysis/equiv_t2_backbone_e3nn.py` — the broad-backbone generalization of the
validated ring fitter (`equiv_t2_e3nn.py`, frame-split T2 R²=0.466 reproduced).
Predicts the per-atom DFT **T2** across the 8 backbone atom-type strata
(N / CA / C / O / HN / HA / HA2 / HA3) from the heterogeneous broad_backbone
substrate (rings + anisotropic bonds + charge sites), on real e3nn, NO recompute.

Design (the new judgment vs the single-mechanism ring case):
- **Shared** angular machinery: `o3.spherical_harmonics("2e")` of each source's
  emitted local-frame direction, the scatter-pool (`index_add_`) to (atom,frame)
  groups, the per-atom de-meaning, and `change_of_basis.lib_to_e3nn` (the FROZEN
  C, reused — not re-derived).
- **Per-source-TYPE radial channels**: three independent invariant MLPs
  (ring / bond / charge), each fed ONLY the invariants its own kind emits
  (ring: r, cosθ, ring_intensity; bond: r, cosθ, bond_category; charge: r,
  source_q_e). A ring and a charge at the same r must scale their l=2
  contribution by different laws; a shared radial gate would force a law they do
  not obey. The heterogeneous Deep-Sets pool over all three kinds is the broad
  relationship's reason to exist.
- Per stratum: frame-split T2 R² (the GATE) + |T2| r (rotation-invariant
  magnitude) + leave-atoms-out R² (REPORTED, not gated — thin per backbone
  stratum in one protein; effective N printed alongside every number).

## The substrate it consumes (all C++-emitted; nothing recomputed)

From `--case broad_backbone` (BroadBackbone.cpp / BroadBackboneSink.cpp), in
`<out_dir>`:
- `broad_backbone_aggregated.csv` — one row per (atom,frame): `row_id`,
  `atom_index`, `atom_name`, `frame_variant` (FrameVariant ordinal), `h5_row`,
  `dft_present`, `dft_local_frame_valid`, the field/μ columns, etc.
- `broad_backbone_sources.csv` — one row per (atom,frame,source): `row_id`
  (join key), `mechanism` ∈ {ring,bond,charge}, `disp_local_{x,y,z}`, `r`,
  `cos_theta`, `ring_intensity`, `bond_category`, `source_q_e`,
  `is_self_or_bonded`.
- `broad_backbone_aggregated_target_local_T2.npy` — **(agg_rows, 5)**, the
  producer's local-frame library-basis DFT T2; the REQUIRED fit target, row-aligned
  with the aggregated CSV. The fitter FAILS LOUD if absent (no Python projection
  fallback). It is mapped library→e3nn-2e by `cob.lib_to_e3nn` (one matmul).

Strata are identified from EMITTED columns only: `frame_variant` fixes the class
(N: {4,5}, CA: {6}, C: {7}, O: {8}, HN: {1,2}, HA-class: {9}); `atom_name` splits
the HA class into HA / HA2 / HA3. No Python chemistry.

## EXACT run command (codex)

```bash
cd /shared/2026Thesis/nmr-shielding/h5-reader/src/rediscover/analysis

# LD_LIBRARY_PATH per ENV.md — torch cu130 SEGFAULTS on import without the cu13 libs.
NV=$HOME/.local/lib/python3.12/site-packages/nvidia
TORCHLIB=$HOME/.local/lib/python3.12/site-packages/torch/lib
export LD_LIBRARY_PATH="$NV/cu13/lib:$NV/cudnn/lib:$NV/nccl/lib:$NV/cusparselt/lib:$NV/nvshmem/lib:$TORCHLIB:$LD_LIBRARY_PATH"

# out_dir = /tmp/rdc-broad-backbone (codex is re-extracting the broad substrate there).
REDISCOVER_OUT=/tmp/rdc-broad-backbone \
  /usr/bin/python3 equiv_t2_backbone_e3nn.py /tmp/rdc-broad-backbone --cross both
```

If the broad substrate is in a different dir, pass it as the positional arg. The
target NPY `broad_backbone_aggregated_target_local_T2.npy` MUST be present in that
dir (re-run the extractor `--case broad_backbone` if only the CSVs are there — the
NPY sidecar is part of BroadBackboneSink's commit).

## The gate

Per the broad-backbone gate (`BROAD_BACKBONE_NEXT.md`): the fit runs across ALL 8
backbone strata (no frame/stratum blocks), produces a per-stratum frame-split T2
R² + |T2| r, and reports effective N per stratum. This is **correlate-not-match,
report-don't-oversell** — there is NO numeric R² threshold to "pass." The
deliverable is: does the equivariant heterogeneous pool COMPOSE and produce
honest per-stratum numbers with effective-N context? (The scalar σ_iso first-pass
in STATE.md found HN/N ≈ 0.45, CA ≈ 0.055 — a T2 vector fit is a different,
harder target; expect lower, report it straight.)

## What codex should check / critique

1. **Runs clean** on the real broad substrate; prints the per-stratum table + the
   final summary; no crash on empty/thin strata (guarded: strata with <2 coupled
   atoms or <4 groups are reported, not fit; LOAO returns NaN under 3 atoms).
2. **Row alignment**: the target NPY rows line up with `broad_backbone_aggregated.csv`
   rows (the fitter slices the NPY by the aggregated-CSV positional index, then
   joins sources by `row_id`). Confirm `len(NPY) == len(aggregated.csv)` (the
   fitter fails loud otherwise) and that the per-stratum group counts look right
   (≈ backbone atoms of that class × DFT frames).
3. **No physics recompute** (the law): grep the new file — it never reads the
   `dipolar` kernel column, never forms `(3cos²−1)/r³` or `q·d/r³`, never projects
   a tensor. It reads `disp_local`/invariants and lets e3nn do the angular math;
   the target is the emitted NPY mapped by the frozen `cob.lib_to_e3nn`. The only
   "recompute" in the tree is the labeled `test_change_of_basis.py`.
   Suggested grep (from `analysis/`):
   ```bash
   grep -nE 'dipolar|3 ?\* ?cos|cos.*cos.*-|r \* r \* r|q ?\* ?disp|project' equiv_t2_backbone_e3nn.py
   # expect: only the change_of_basis import + comments that say "NOT a projection"
   ```
4. **Change-of-basis reuse**: the fitter imports `change_of_basis as cob` and uses
   `cob.get_C()` / `cob.lib_to_e3nn` — it does NOT re-derive C. (Run
   `test_change_of_basis.py` first to confirm the pinned C still passes under the
   installed e3nn 0.6.0.)
5. **The flagged schema gap** (below): decide whether to extend the C++ emit. The
   fitter is correct without it (disp-only); the gap caps the angular richness.

## FLAGGED schema gap (a finding, NOT worked around)

The ring fitter oriented its l=2 contribution by the **ring normal** n̂ (the
dipole axis), reading the emitted `source_normal_local_*` columns of
`ring_current_sources.csv`. The **broad** source schema
(`BroadBackboneSink::kSourceHeader`) emits only `disp_local_*` — it does NOT write
`source_normal_local_*` (rings) or `bond_axis_local_*` (bonds), even though those
vectors are already computed and stored in the C++ `SourceSlot` by `ringAttacher`
/ `bondAttacher` (BroadBackbone.cpp lines 173, 236). So the angularly-complete
design — Y2(n̂) for rings, Y2(axis) for bonds, plus the disp⊗axis cross term the
ring fitter used — **cannot be built from the broad substrate as currently
emitted**.

The fitter therefore uses the one emitted local-frame direction, `disp_hat`, for
every source kind. This is **exact for charge** (the Coulomb field is along the
displacement) and carries the ring/bond axis ANGLE invariantly via the emitted
`cos_theta`, but not the axis VECTOR. We did NOT recompute the missing vectors in
Python (that would be the projection end-run the law forbids).

The fuller model is wired but DORMANT: `--with-axes` auto-activates Y2(axis)+cross
FOR ANY kind whose axis columns are present. **Suggested C++ follow-up (separate
codex task, additive, one sink):** add `source_normal_local_{x,y,z}` and
`bond_axis_local_{x,y,z}` to `kSourceHeader` + `writeSourceRow` in
`BroadBackboneSink.cpp` (the `SourceSlot` fields already exist — `s.source_normal_local`,
`s.bond_axis_local`). Then re-run `--case broad_backbone` and run this fitter with
`--with-axes` to light up the complete equivariant model with no Python change.
This mirrors what the ring substrate already emits, so it is a known-good shape.

## Discipline (unchanged)

Model untouched in C++; only the Python fitter added. Frozen C reused. No
recompute outside the labeled test. NEVER MERGE.

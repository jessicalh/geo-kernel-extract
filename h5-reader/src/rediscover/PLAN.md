# Rediscover — plan + stop point (2026-05-31)

A forward plan and a deliberate stop, written after the first end-to-end pass.
Living doc; reads with `STATE.md` (freshest state) + `analysis/FINDINGS.md`
(the numbers) + `GUIDANCE.md` / `DESIGN.md`.

## What we actually built (read this first)

Not "a ring-current result." An **instrument**: a typed C++ substrate that emits
a per-(atom, frame) source-neighbourhood + DFT target in one consistent frame,
plus a rediscovery pipeline (scalar sum-pooling / PySR for T0; equivariant
sum-pooling for T2) with instantaneous, leave-atoms-out honesty. We proved it on
**one (kernel, case) cell** — ring current on 1P9J — and it works: the Pople
kernel falls out (scalar k≈21 ppm·Å³, held-out-atom R²=0.62 vs independent DFT;
PySR gives the closed form), and the **tensor** carries signal too (equivariant
T2 modulation R²=0.44 held-out, |T2| r=0.75). McConnell's scalar form recovers;
its producer-kernel reconstruction has a known gap.

## The framing (Jessica, 2026-05-31)

**It will always be thin** — that is intrinsic to evaluating one case, not a flaw
to grind down by piling more frames onto 1P9J. The contribution is the
**instrument + the per-cell indicative finding**, and the fact that **aspects
generalise to other relationships**. A single (kernel, case) cell is a
proof; the value compounds across cells. Don't over-invest in any one cell.

## Generalisation — two axes

**Axis 1 — other kernel relationships, same case (1P9J).** The substrate and the
fitters are relationship-agnostic; a new relationship is a new
`RediscoveryExtraction` (stratum + source-kind + kernel form) reusing the shared
identity / DFT-target / frame / sink machinery (this is exactly what the
`RediscoveryExtraction` interface in `DESIGN.md` was built for). The H5 already
carries the per-atom kernel time-series to cross-check each, and ORCA gives the
one DFT target they all explain a piece of:
- **APBS electric-field-gradient** — `apbs_efg_time_series/t2` is already a T2
  in the H5; sources are charges. A clean, concrete next relationship: charges as
  l=0 sources at displacement vectors → predict the EFG T2 with the SAME
  equivariant machinery. Strong candidate to do next.
- **APBS E-field** — `apbs_efield_time_series/xyz` (an l=1 vector target);
  equivariant l=1 sum-pooling.
- **Charge multipoles** (dipole/quadrupole) from charges (ff14SB / AIMNet2),
  per the conventions doc's multipole definitions.
- **Larsen H-bond** — exchangeable-H stratum, donor/acceptor geometry.

Each new cell reuses: the frame check, the typed local frames, the two row
kinds, the leave-atoms-out protocol. The marginal cost is one extraction class +
a fitter run, not a rebuild.

**Axis 2 — other cases, same kernel.** 1P9J is one structure. The same substrate
runs on any calcset with DFT. Cross-case is where the "thin between-atom"
limitation lifts: many aromatic H across many environments test the universal
coefficient's transferability for real (the within-1P9J leave-atoms-out is a
hint; cross-protein is the test). This connects to the broader Stage-2 fleet —
deferred, scoped separately, not this branch's job.

## Optional sharpeners on the current cell (parked, not blocking)

- leave-ATOMS-out for T2 (sharper than the frame-split already run);
- the literature-coefficient-FIXED de-circularising check (un-fitted constant →
  DFT) — fully removes the residual circularity in the scalar claim;
- a richer e3nn gated network (vs the minimal 3-channel l=2 ansatz);
- McConnell producer-kernel reconstruction gap (R²≈0.55 — likely a fuller
  anisotropy than one bond-axis angle; check the producer's model).

## Stop point

Stopping here, deliberately. Durable artifacts:
- **Code**: `src/rediscover/*` (canonical substrate; builds + runs on 1P9J;
  `h5reader_rediscover_tests` green; T2 frame check in-tree).
- **Substrate**: `/tmp/rediscover-out-v2/` (identity-clean; machine-local).
- **Findings**: `analysis/FINDINGS.md` (numbers + the circular-vs-real audit +
  the instantaneous framing).
- **Analysis**: `analysis/*.py` (look01–03, sumpool_*, refine_kernel,
  credibility*, equiv_t2; PySR venv at `analysis/venv`, gitignored).
- **Memory**: `project_rediscover_state`, `feedback_no_parallel_h5_in_python`,
  `feedback_build_inmemory_export_dont_relitigate`.

Next session picks up from here: most likely point the instrument at a second
relationship (EFG is the natural one — its T2 is already in the H5), or run the
parked sharpeners on the ring-current cell. No blockers; nothing mid-flight.

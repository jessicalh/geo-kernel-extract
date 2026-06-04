# Rediscover — work catalog (what we want overall)

> **Current backlog trued 2026-06-04.** Old 1P9J LOAO/between positives are
> retracted and must not guide priority. The near-term work is 720-WT static
> ingest/pilot, Stage3/joint-model framing, and parked clean re-runs only if
> the lead restarts them; this file preserves the older catalog below as
> historical context.

The actionable work breakdown, complementing `REDISCOVERY_MAP.md` (the *science* roadmap —
two depths, per-cell workflow, constraints). This is the concrete unit list, structured so
each is a fast **codex-delegated cell** against the now-complete emit-catalog + the exemplar
template. Reads with `REDISCOVERY_MAP.md`, `STATE.md`, `analysis/PATTERNS.md`.

**Historical enabling-condition note:** this predates the true-LOAO retraction.
Do not treat the old 1P9J between/LOAO numbers or near-term 750/engine ordering
as current priority guidance.

## A. Done — the spine + the exemplar arc
- Functional API (engine / composed `Relationship` / curried closures) — built, byte-parity (`99cdc85`).
- `ring_current` — full arc, **Depth A**: the Pople law fell out (PySR; non-circular k≈21 / R²=0.62).
- `mcconnell` — partial (form recovers R²≈0.85).
- `broad_backbone` — **the exemplar**: B (8 strata) + A (bond McConnell recovers, atom-split 0.48)
  + de-circularised (un-fitted literature kernel predicts N 0.69 / O 0.53 / C 0.51).
  (`35f3768`, `901d1df`, `184a1ee`.)
- The e3nn fitter, the frozen change-of-basis, the literature-kernel-T2 emit — the template machinery.

## B. Engine (gates the breadth of C)
- **#29 — unify the engine fold/sink.** The broad case forced a sibling runner; unify so each
  new relationship rides ONE engine. Design DRAFTED (`ENGINE_TOTALITY_DESIGN.md` +
  `MODEL_PLACEMENT_PROPOSAL.md`); process = review (talk + independent Opus critique) → implement.
  **Do before C.**

## C. Item-2 calculators — each a per-cell arc through the unified engine (copy the exemplar)
Each: express through the engine → emit (case-by-case, spine-side) → capture (B) → distill (A)
→ de-circularise → report (correlate-not-match, effective N). Order of attack:
1. **efg** (APBS EFG → T2) — natural next; reuses the exemplar's equivariant machinery; APBS data present.
2. **buckingham_efield** (APBS E-field → l=1).
3. **charge multipoles** (dipole / quadrupole; ff14sb / aimnet2 / mopac as a parameter sweep).
4. **larsen_hbond** (exchangeable-H stratum, donor/acceptor geometry).
5. **AIMNet2 CRG** (charge-response-gradient).
6. **AIMNet2 embedding** (per-atom learned rep).

## D. Infrastructure / discipline
- **#30 — the repo tripwire guard** (backstop to `analysis/PATTERNS.md`: the doc is the push,
  the test is the catch). Build after the exemplar so it encodes the clean tree.

## E. Parked — waits on data / contributors
- **#35 — differencing system-ID** (Δσ vs Δgeometry; the response handle). PARKED: the smoothness
  gate fails at the current every-other-frame DFT sampling (Δσ ~90% noise, lag-1 ρ≈0.05) — it
  needs **finer DFT sampling** (consecutive/closely-spaced frames); and the **disambiguation**
  needs more contributors (item C). The script (`analysis/differencing_system_id.py`) + scope
  (`DIFFERENCING_EMIT_SCOPE.md`) exist and are reusable the moment finer data lands. Honest
  negative for now; the gate did its job.

## F. The data dependency
- The **1P9J 750-DFT set lands in FULL in ~3 days** (every-other-frame; the 750 took 5
  machines × 4 cores × 64 GB each, ~a month). NEAR-TERM: re-run the static rediscovery (the
  exemplar + item-C cells) on the full 750 — more frames, deeper per-stratum N. The current
  ~500-frame subset already supports items A–D.
- **Differencing (#35) is gated on SPACING, not count.** Its smoothness gate failed because
  consecutive DFT frames are far apart (every-other); the full 750 (same spacing) does NOT fix
  that — differencing needs a DENSE consecutive-frame DFT burst (a separate campaign) + more
  contributors. Transferability (G) needs the much larger, later fleet DFT.

## G. Future axes (Stage-2+)
- **Transferability across proteins** (the fleet) — the *other* ensemble axis; needs the fleet
  + the DFT (F).
- **Ensemble / coupled-factorization model** (frames × atoms × mechanisms × T2) — the
  rediscovered kernels as the mechanistic basis; the de-circularised *literature* kernels as a
  physics-grounded basis for N/O/C.
- **Symbolic distillation maturing** as the contributor set grows (more mechanisms → more to
  distill and, eventually, to disambiguate).

## The order
**B (#29 engine) → C (the calculators — the bulk; fast codex cells now that catalog + template
are done) → D (#30 guard).** E (#35) and the finer-sampling analyses wait on F (the full DFT).
G is Stage-2+.

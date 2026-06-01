# Next-session kickoff — rediscover (current as of 2026-06-01 EOD)

Continuing `rediscover` on branch `h5-reader-pysr-spike` — **NEVER MERGE**
(experimental spike; no merge / switch-to-master / PR / rebase, ever).

This supersedes the stale `NEXT_STEP_BRIEF.md` and `BROAD_BACKBONE_NEXT.md` (both
DONE / historical — ignore them).

## Read first
- **`REDISCOVERY_MAP.md`** — the science roadmap (WHERE WE'RE GOING): two depths (a
  law falls out / signal captured), the status grid, item 1 (complete the template) →
  item 2 (the other calculators), the constraints.
- **`STATE.md`** — freshest state + numbers (most recent updates at the top).
- **`analysis/PATTERNS.md`** — the controlling discipline for the Python-consumer
  surface; FRONT-LOAD it in every analysis/fitter agent brief.
- **`analysis/equiv_t2_backbone_e3nn.py`** — the EXEMPLAR (WHAT IT LOOKS LIKE), the
  template every future cell copies.
- Memory: `project_rediscover_state`, `feedback_no_python_physics_except_labeled_integrity_test`,
  `feedback_controlling_patterns_docs_steer`, `feedback_catch_debt_in_the_moment`,
  `feedback_motivated_agent_briefs`, `reference_subagent_build_agency`.

## FIRST ACTION — verify the sandbox fix (a restart took effect)
`settings.local.json` now has `permissions.defaultMode: bypassPermissions` — the fix
for subagent build/run agency (it's a permission-MODE issue, NOT sandbox; the per-call
`dangerouslyDisableSandbox` override does NOT help subagents). On this fresh session,
TEST it: spawn an Opus subagent that runs `import torch` / `cmake --build`. Works →
Opus agents have build/run agency, the author-with-Claude/build-with-codex split
collapses (route by best tool). Still denied → codex remains the build/run path.

## Where we are (the package, all on h5-reader-pysr-spike)
- Functional API built + byte-parity-proven (`99cdc85`); broad-backbone built+gated
  (`35f3768`); e3nn rebuild retiring the Python end-run (`dcc4a46`); the backbone
  equivariant-T2 EXEMPLAR (`901d1df`) — 8 strata, orientation vectors help all.
- **The science (`REDISCOVERY_MAP.md`):** ring current = a law fell out (Depth A —
  PySR Pople form, k≈21 / R²=0.62 non-circular). Backbone tensor (T2) = signal captured
  (Depth B — HN .76 / O .72 / HA .67 / N .65 / C .59 / CA .43; CA tensor 0.43 vs scalar 0.055).

## Status: the backbone exemplar is COMPLETE — full B→A→de-circularised arc
`#33` + `#34` DONE. Backbone verdict: **Depth B captured** (HN .76 / O .72 / HA .67 /
N .65 / C .59 / CA .43); **Depth A partial** (bond McConnell r⁻³ recovers, atom-split
0.48; ring null; charge partial); **de-circularised** — the UN-FITTED literature-
coefficient kernel T2 predicts DFT for N (component r 0.69), O (0.53), C (0.51), while
CA/HN/HA are weak (they need the fit). So the textbook physics, un-fitted, genuinely
carries the heavy-atom + amide-N tensors. The exemplar now demonstrates the ENTIRE arc
end-to-end (`analysis/BACKBONE_LAW_EVIDENCE.md` + the de-circularising CSV under
`/tmp/rdc-broad-backbone-axes`). Proceed per the order below.

## Then, per REDISCOVERY_MAP.md (the order)
1. **#33 completes the template** (item 1) — finish before going broad.
2. **#29 — engine totality unification** gates item 2's breadth (each new relationship
   rides ONE engine, not a sibling runner). Design DRAFTED (`ENGINE_TOTALITY_DESIGN.md`
   + `MODEL_PLACEMENT_PROPOSAL.md`); process = review (talk + an independent Opus
   critique) → implement. Do NOT go broad before this.
3. **Item 2 — the other calculators** (`#25`/`#27`): efg #4 first (APBS EFG → T2; same
   equivariant machinery; data present), then buckingham_efield / charge multipoles /
   larsen / AIMNet2 CRG + embedding. Each copies the completed exemplar + rides the
   unified engine *(any quantity a cell needs that the substrate lacks → a C++
   emit-extension, spine-side, never Python)*.
4. **#30 — the repo tripwire guard** (backstop to `analysis/PATTERNS.md`).

## Discipline (`analysis/PATTERNS.md` + the memories)
Model-is-spine (C++ emits; Python fits via e3nn on the emitted substrate); NO recompute
(projection / kernel / field) outside a labeled pinned test; equivariant = e3nn (never
hand-rolled); reuse the frozen `get_C()`; correlate-not-match (report effective N, don't
oversell thin strata); catch debt in the moment; front-load `analysis/PATTERNS.md` in
every fitter brief; **USE + SUPPORT the functional interface** (express each cell through
the engine as a composed relationship; extend it cleanly + aesthetically; don't bypass it
or proliferate sibling runners — unify #29); codex grinds / Opus + lead judge; NEVER MERGE.

## Open tasks
#33 (distillation, in flight) · #29 (engine unification, design drafted) · #25/#27
(item 2 calculators) · #30 (tripwire backstop).

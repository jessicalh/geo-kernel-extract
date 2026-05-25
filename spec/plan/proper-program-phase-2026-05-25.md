# "Proper program" phase — plan (2026-05-25)

Goal: a proper program we can build and run anywhere — clean invocation
surface, config-driven paths/DBs, and every source file documented + swept
clean. Builds on the just-landed readability sweep (`aef8f99`) and the
commit-gate removal (`e929117`).

Execution order: **2 → 5 → 3 → 4** (clean the surface + wiring first, then
document, then the catch-all sweeps last).

---

## 2. Invocation / CLI cleanup  — DONE (`ff31d9c`, `8c08ecd`)

Removal pass against the `nmr::cli` surface, back to the canonical 5 modes.
Landed:

- **`RunConfiguration` readability refactor** (`ff31d9c`): dead
  `ScanForDftPointSet` deleted (no CLI path built it); the trajectory-result
  registration spray collapsed to labelled `Produces<TRs...>` groups
  (order-preserved → bit-identical); member/accessor renames
  (`DeferredResult`, `results_to_build_`, `ResultsToBuild()`) for object-model
  clarity.
- **Coulomb retired from production + cruft flags removed** (`8c08ecd`): APBS
  and AIMNet2 are always on (not switchable); `--no-apbs` / `--no-coulomb`
  deleted; home-rolled Coulomb runs only in the `--trajectory --mopac`
  (FullFat) charge-source reconciliation probe. The five `coulomb_*` SDK specs
  stay (already optional) but are now FullFat-only.

Kept on purpose: the frame-emit window + stride controls (they compose, not
duplicate) and the FullFat MOPAC-vs-FF14SB Coulomb probe.

## 5. Configuration / portability

All DB DSNs and filesystem paths live in config (TOML / `RuntimeEnvironment`) —
**zero hardcoded paths**. Outcome: the program builds and runs anywhere given
the config.

## 3. Full Doxygen / comment pass — EVERY `.cpp` and `.h`

- Proper doxygen on every file.
- Calculator prose **starts from `intermediate/calculators_draft4.md`** (the
  calculator description set; drafts 1–4 are in `intermediate/`).
- `\ref` cross-links: each `*Result` → the `ConformationAtom` / `TrajectoryAtom`
  fields it **deposits** and the **depositor functions** (so the data-flow is
  navigable in doxygen).
- Grounded, non-AI, non-run-on comments (the readability-sweep discipline).
- Variable naming finished where the sweep didn't reach.
- `PhysicalConstants.h` rendered as a **doxygen table**.
- "Full cleanup and perfect" — comprehensive, not sampled.

## 4. Repo-wide sweeps  *(last — catch what's left)*

- **Magic numbers**: hardcoded literals → named constants / config, across
  *everything* (not just calculators).
- **Silent failures**: swallowed errors / unlogged default-returns / "return
  false with no log" across *everything* — make failures loud (per the
  fail-loud discipline).

---

## Anchors / inputs
- Calculator prose: `intermediate/calculators_draft4.md`.
- Canonical modes: CLAUDE.md "CANONICAL 5-MODE SPEC".
- Naming baseline: the 19-algorithm readability sweep (`aef8f99`).

## Out of scope (2026-05-25)
The **propka / OF3→GROMACS prep case** (use the model + published rules to
build a GROMACS-safe PDB with terminus/CYS-disulfide handling, an atom-adding
build path, and a separate prep DB of pH/pKa/ionisation) is **out of scope** —
dropped 2026-05-25. CLAUDE.md's "propka… out of scope, protonation done before
`nmr_extract`" line therefore stands as written; no canonical-spec rewrite is
pending.

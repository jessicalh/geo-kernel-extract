# Post-compact prompt for the next session

Drop this verbatim into the first user message after compact. Designed
to land the next session in the typed-tripeptide-topology refactor
with full context loaded.

---

We're continuing tripeptide-calculator work. Last session
(2026-05-10/11) landed `TripeptideBackboneShieldingResult` +
`TripeptideNeighborShieldingResult` + `TripeptideDftTable` +
`TripeptidePoseAssembler` with smoke tests green on
1UBQ_pm6dh3plus.pdb. The architectural follow-up is the **typed
tripeptide topology refactor** — design landed but not implemented.
This session is Phase 1 of that rollout.

## Read in this order before any code

1. **`spec/plan/typed-tripeptide-topology-design-2026-05-10.md`** —
   the durable design doc. Plan-agent architecture + user pose-metadata
   enrichment. Read end-to-end. The implementation order (Phase 1 / 2 / 3)
   is the work plan.

2. **`spec/plan/session_2026-05-10_part2_tripeptide_calculators.md`** —
   what landed last session, with the bugs caught + fixed and the
   discipline that's now in place.

3. **`PATTERNS.md`** — full read. Particularly: substrate boundary
   discipline, "rule application architecture for runtime atom-name
   canonicalisation" (the same shape applies to DFT-side substrate
   composition), no-adapter-classes, fail-loud on substrate gaps.

4. **`OBJECT_MODEL.md`** — particularly the `LegacyAmberTopology` +
   `ComposeAtomSemantic` sections. The DFT-side topology layer is
   modeled on that pattern.

5. **`src/SemanticEnums.h`** — the typed vocabulary
   (`AtomMechanicalIdentity`, `BackboneRole`, `Locant`, etc.) + the
   existing `AtomSemanticTable` struct. The DFT side will populate
   these same structs per row.

6. **`src/generated/LegacyAmberSemanticTables.h`** and `.cpp` (skim
   the format) — the existing build-time-generated typed table for
   protein residues. The new `TripeptideOrderingTables.{h,cpp}` will
   live alongside.

7. **`tools/topology/README.md`** and
   `tools/topology/build_semantic_tables.cpp` — the build-time
   generator. Phase 1 step 4 extends this to read the tripeptide
   orderings TOML and emit the new file.

8. **`src/TripeptidePoseAssembler.{h,cpp}`** — the current heuristic
   matcher we're replacing in Phase 2. The math (Kabsch, sidechain
   re-rotation, tensor rotation) is correct — only the matching
   logic is heuristic.

9. **`src/TripeptideDftTable.h`** — the DFT-record struct + libpq
   loader. Note `TripeptideDftRecord` has slots ready for the new
   `atom_semantics` field (will need to add).

## Auto-loaded memory entries (already in context)

These are listed in CLAUDE.md and should be present:

- `feedback_two_path_validation` — cross-substrate matching discipline.
- `feedback_residual_as_ml_feature` — Vec3 not scalar.
- `project_tripeptide_calculators_landed` — what shipped.
- `project_typed_tripeptide_topology_design` — what's queued.
- `project_larsen_neighbor_axa_reuse` (CORRECTED 2026-05-10) — paper-
  verified cap-side reading.
- `reference_gotham_assembler` (UPDATED 2026-05-10) — port complete
  with four improvements.
- `project_post_csa_calculator_queue` (UPDATED 2026-05-10) — BB +
  Neighbor LANDED; HBondHα remains queued.

## Operational state (verified end of last session)

- Local Postgres on batcave has full `tensorcs15` replica (1,870,631
  rows including SER PBE). User `jessica` has SELECT on schema
  `public`.
- `~/.nmr_tools.toml` has `[databases].tensorcs15 = "host=/var/run/postgresql
  dbname=tensorcs15 user=jessica"`.
- `nmr_shielding` library + `nmr_extract` + `structure_tests` +
  `unit_tests` all build clean. 114/114 unit tests pass. BB +
  Neighbor smoke tests pass on 1UBQ_pm6dh3plus.pdb.
- Larsen archive supplements at `/mnt/expansion/larsen_archive/`:
  `hydrogenbondnmrlogs.tar` (503 MB, not yet ingested),
  `structures/{1UBQ,1UBQ_pm6dh3plus,2OED_pm6dh3plus}.pdb` (test fixtures).

## Concrete first task — Phase 1 of the typed-topology refactor

Per `typed-tripeptide-topology-design-2026-05-10.md` Phase 1:

1. **psql dump.** Write `scripts/dump_tripeptide_orderings.py` (Python,
   uses `psycopg2`) that queries one row per `(central_residue,
   frame_type)` combination from `tensorcs15.raw_dft_calculations`
   and emits per-atom `(atom_idx, element, atom_name_if_present)` to
   JSON. Use the local DSN
   `host=/var/run/postgresql dbname=tensorcs15 user=jessica`.

2. **User reviews the dump.** Spot-check ALA, SER, ARG, LYS orderings
   for sanity. **Critical:** verify whether ORCA SER rows have the
   same atom ordering as Gaussian SER rows. If different, two SER
   ordering tables are needed.

3. **Hand-author `data/topology/tripeptide_orderings.toml`.** One
   section per `(residue, frame_type)` with atom names in DFT-row
   order. Atom names are the standard PDB names from
   `AminoAcidType::atoms`.

4. **Per-residue governing-dihedral assignment.** Walk each residue
   from CA outward and assign `GoverningDihedral` per atom based on
   which χ bond it's past. Document in
   `data/topology/governing_dihedrals.toml` or inline.

5. **Add `GoverningDihedral` enum** to `src/SemanticEnums.h`.

6. **Extend `tools/topology/build_semantic_tables.cpp`** to read the
   two TOMLs and emit `src/generated/TripeptideOrderingTables.cpp`
   with `constexpr std::array<TripeptidePositionEntry, N>` arrays
   per `(residue, frame_type)`. Also emit `kAceAtoms` / `kNmeAtoms`
   into the main residue table list.

7. **Regenerate + commit.** Run the generator, commit the regenerated
   `.cpp` files.

8. **CMake.** Append `src/generated/TripeptideOrderingTables.cpp` to
   `nmr_shielding` sources.

9. **Substrate parity test** —
   `tests/test_tripeptide_substrate_parity.cpp`. One DB query per
   `(residue, frame_type)`, assert ordering-table element pattern
   matches the DB row's element sequence.

End of Phase 1 = parity test green. Phase 2 (the runtime composer +
assembler rewrite) is the next session.

## Discipline reminders

- **T2 sacred** — preserve `Mat3` + `SphericalTensor` end-to-end. JSONB
  has them pre-decomposed; use directly.
- **No string traversals** after the load boundary. The TOML is a
  load-boundary input to the build-time generator; runtime consumes
  typed tables only.
- **Fail loud on substrate gaps.** Per `ComposeAtomSemantic`'s
  discipline. Unknown residue / missing identity = FATAL+abort.
- **Two-path validation** when matching across substrates. The DFT
  side now becomes typed; the protein-side `SemanticAt` cross-check
  remains as defence-in-depth.
- **Residual as Vec3 ML feature.** Don't gate emission on magnitude;
  the model attends to direction + magnitude.
- **Frame_type passthrough.** SER PBE rows must route distinctly from
  OPBE rows; the method tag carries through to per-atom storage.

## Alternative opening — run on 1P9J trajectory first

If you'd rather surface bugs under MD geometry variance before the
refactor, the alternative first task is:

1. Build a fixture loading 1P9J trajectory frames (the actual Stage 2
   study system; 750 frames at 20 ps stride).
2. Run BB + Neighbor calculators on each frame.
3. Check residual distributions, methods discrimination, per-atom
   assignment counts under varied MD geometry.
4. Surface any bugs first, refactor second.

The typed-topology refactor is the more architectural follow-up; the
1P9J trajectory smoke is the more empirical follow-up. User
explicitly queued the typed-topology path first per the session-wrap
discussion, but either is defensible.

## Final note

The current calculator code is correct enough to ship features
through. The typed-topology refactor is principled-clean-matching
work, not bug-fixing. If something more urgent comes up in the
session opener, treat this as queued, not blocking.

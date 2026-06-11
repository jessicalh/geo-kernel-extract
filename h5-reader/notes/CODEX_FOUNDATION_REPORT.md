# CODEX FOUNDATION REPORT

Date: 2026-06-11

## Scope

Implemented the FOUNDATION stage only. I stopped before the `per_atom_substrate`
migration and before any 720 migration. The old `per_atom_substrate` emitter was
not migrated, and the producer extractor catalog was not edited.

Commits made during this stage:

- `55b4758 Add substrate output parity harness`
- `5655ac0 Add FieldKind query foundation`
- `8d6a035 Add resident substrate indexes`
- `2be4944 Add reducer parity gates`
- `9d591db Add query-backed Parquet audit consumer`

## What Is Complete

1. Tests and old-output parity harness landed before the foundation code changes.

   The harness reads today's output rows/sidecars, validates the per-atom row
   contract, and compares old target sidecars to producer truth:
   `src/rediscover/SubstrateParity.cpp:187`, `src/rediscover/SubstrateParity.cpp:221`,
   `src/rediscover/SubstrateParity.cpp:228`, `src/rediscover/SubstrateParity.cpp:240`.
   Producer truth is read through `Catalog::value(FieldKind::OrcaTotal, ...)`:
   `src/rediscover/SubstrateParity.cpp:289`. Tests are in
   `tests/rediscover_substrate_parity_tests.cpp`.

2. FieldKind gather primitive is in place.

   `FieldRef::Producer` binds query fields to `io::FieldKind`:
   `src/rediscover/RunQuery.cpp:20`. `GatherField` gathers through
   `body.catalog.value(body, ref.kind, ...)`, not through `verbs::value` or
   `ArrayId`: `src/rediscover/RunQuery.cpp:49`.

3. Resident indexes now include the required foundation indexes.

   IUPAC-name and chemical-category indexes are defined at
   `src/rediscover/ResidentIndexes.h:87` and
   `src/rediscover/ResidentIndexes.h:116`, and populated at
   `src/rediscover/ResidentIndexes.cpp:148`. Secondary structure has explicit
   `Unknown` states and all-zero DSSP rows decode to unobserved/unknown:
   `src/rediscover/ResidentIndexes.h:25`, `src/rediscover/ResidentIndexes.h:37`,
   `src/rediscover/ResidentIndexes.h:69`. Dihedrals carry continuous radians,
   fixed bins, and presence: `src/rediscover/ResidentIndexes.h:168`,
   `src/rediscover/ResidentIndexes.h:185`,
   `src/rediscover/ResidentIndexes.h:202`,
   `src/rediscover/ResidentIndexes.h:231`.

4. Phi/psi sign is data-parity driven.

   The resident index builder computes coordinate-derived phi/psi samples and
   compares straight versus negated NPY values:
   `src/rediscover/ResidentIndexes.cpp:78`. The selected policy is applied when
   populating phi/psi, and unresolved policy leaves the angle absent rather than
   guessing: `src/rediscover/ResidentIndexes.cpp:177`,
   `src/rediscover/ResidentIndexes.cpp:182`,
   `src/rediscover/ResidentIndexes.h:212`,
   `src/rediscover/ResidentIndexes.h:225`.

5. Two-phase selectors and explicit traversal domains are implemented.

   Selectors now have separate atom and frame predicates:
   `src/rediscover/RunQuery.cpp:91`, `src/rediscover/RunQuery.cpp:99`,
   `src/rediscover/RunQuery.cpp:107`. The query runner applies atom prefilters
   before per-frame filters: `src/rediscover/RunQuery.cpp:180`,
   `src/rediscover/RunQuery.cpp:198`, `src/rediscover/RunQuery.cpp:247`.
   Traversal domains are explicit for `AllFrames`, `DftRows`, and
   `TargetPresentRows`: `src/rediscover/RunQuery.cpp:220`.
   Index-backed selectors exist for IUPAC name, chemical category, secondary
   structure, and dihedral bins/ranges: `src/rediscover/RunQuery.cpp:284`,
   `src/rediscover/RunQuery.cpp:293`, `src/rediscover/RunQuery.cpp:302`,
   `src/rediscover/RunQuery.cpp:324`.

6. Reducer parity gates are present and do not choose an authority.

   The gate runs the composed relationship reducers through existing selectors
   and attachers: `src/rediscover/ReducerParity.cpp:44`. It separately reduces
   today's per-atom pair contributions: `src/rediscover/ReducerParity.cpp:71`,
   `src/rediscover/ReducerParity.cpp:97`. It then compares both sources only
   within each relationship's stratum: `src/rediscover/ReducerParity.cpp:125`,
   `src/rediscover/ReducerParity.cpp:145`, `src/rediscover/ReducerParity.cpp:159`.
   Tests are in `tests/rediscover_reducer_parity_tests.cpp`.

7. Reader-side membrane and thin audit consumer are implemented.

   The new audit result declares its own Arrow/Parquet schema, including query
   rows, resident labels, dihedrals, producer ORCA raw values, and old
   per-atom target values: `src/rediscover/QueryAudit.cpp:210`. It loads the old
   per-atom sidecars as consumed outputs and also runs the old-output parity
   harness: `src/rediscover/QueryAudit.cpp:490`,
   `src/rediscover/QueryAudit.cpp:501`. It reads producer truth through
   `RunQuery` using a two-phase selector and `FieldRef::Producer(OrcaTotal)`:
   `src/rediscover/QueryAudit.cpp:506`,
   `src/rediscover/QueryAudit.cpp:518`. It writes one Parquet table:
   `src/rediscover/QueryAudit.cpp:543`, `src/rediscover/QueryAudit.cpp:546`.
   The CLI entry is `--case query_audit`: `src/rediscover/main_extract.cpp:523`.
   Arrow/Parquet is linked only on `h5reader_extract`:
   `CMakeLists.txt:324`, `CMakeLists.txt:444`, `CMakeLists.txt:449`.

## Validation

- `cmake --build build/linux-gcc --target h5reader_extract -j2`: passed.
- `ctest --test-dir build/linux-gcc --output-on-failure`: 17/17 passed.
- `./build/linux-gcc/h5reader_extract --help | rg -n "query_audit|per_atom_substrate|Which extraction"`: `query_audit` is advertised by the CLI.
- `git diff --check`: clean.

I also searched for existing old per-atom sidecars with:
`find /shared/2026Thesis/nmr-shielding /shared/2026Thesis/shielding-calcsets -name per_atom_substrate_rows.csv`.
No existing `per_atom_substrate_rows.csv` was present, so I did not run a live
`query_audit.parquet` emission against real old sidecars in this turn. The
consumer is intentionally fail-loud if those sidecars are absent:
`src/rediscover/QueryAudit.cpp:135`.

## Not Done

- No `per_atom_substrate` migration was started.
- No 720 migration was started.
- No `QtFieldCatalog`, generated producer field catalog, or `nmr_extract`
  extractor code was edited.
- The reducer parity gate reports drift; it does not adjudicate which reduction
  source is correct. That decision belongs before the migration fire.

## Before The Migration Fire

1. Generate or locate a current old `per_atom_substrate` output, then run
   `h5reader_extract --case query_audit --run <run> --out <old-output-dir>` and
   inspect `old_output_mismatch_rows`, `producer_mismatch_rows`, and the Parquet
   diff columns.
2. Run `AuditReducerParity` on representative real aromatic-H and HN cases.
   Treat any mismatch as a design decision, not as a blind implementation choice.
3. Keep the migration blocked unless the query audit, reducer gate, and
   old-output parity harness all pass on the selected data slice.

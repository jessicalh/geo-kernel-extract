# tools/topology — semantic tables generator

Build-time C++ generator for `src/generated/LegacyAmberSemanticTables.cpp`.

This tool reads the wwPDB Chemical Component Dictionary
(`data/ccd/components.cif`) via cifpp, runs RDKit chemistry
perception (CIP labelling, canonical rank, aromaticity,
hybridisation, bond-type), applies the per-field reconciliation
precedence table specified in
`spec/plan/bones/topology-substrate-implementation-plan-2026-05-05.md`,
and emits a typed-enum constexpr table consumed at runtime by
`LegacyAmberTopology`.

## Why a separate generator binary?

The runtime library `libnmr_shielding.a` MUST NOT link RDKit. The
string-barrier discipline lives at the linker level: chemistry-string
libraries (cifpp, RDKit) appear only in this generator's link line.
The runtime sees only the generated typed-enum table.

## Building the generator

```bash
cmake -B build-gen -DNMR_BUILD_TABLE_GENERATOR=ON
cmake --build build-gen --target build_semantic_tables
```

By default, RDKit is found in the local micromamba `mm` environment
at `/home/jessica/micromamba/envs/mm`. Override with
`-DRDKIT_ROOT=/path/to/rdkit-root` if RDKit is installed elsewhere
(must contain `include/rdkit/` and `lib/libRDKit*.so`).

cifpp is found via the standard project mechanism (`find_package(cifpp)`).

## Running the generator

```bash
./build-gen/tools/topology/build_semantic_tables \
    --ccd    data/ccd/components.cif \
    --output src/generated/LegacyAmberSemanticTables.cpp \
    --log    src/generated/LegacyAmberSemanticTables.log.txt
```

The log file is structured plain-text (section headers + key=value
records). Inspect it once per regeneration to confirm:

- Every standard residue + variant + atom has all 14 fields populated.
- ProchiralStereo entries have CIP-verified provenance.
- Source disagreements (if any) are flagged with both values stashed.

## When to regenerate

Run the generator when:

- The wwPDB CCD updates (rare; new residues, fixed entries).
- RDKit version changes in a way that affects CIP labelling
  (compare logs across versions).
- A new protonation variant or cap residue is added to the
  supported set.

The generated `.cpp` output and `.log.txt` are committed alongside
each other. The runtime build does not invoke this generator.

## What lives where

```
tools/topology/
    CMakeLists.txt                  # opt-in subproject; -DNMR_BUILD_TABLE_GENERATOR=ON
    build_semantic_tables.cpp       # main entry, orchestration
    README.md                       # this file
```

## Status

Step 1 of N (scaffolding) is the current state. Cifpp and RDKit
smoke tests run; no tables are generated yet. Subsequent commits
add: CCD entry parsing for the standard 20 residues, RDKit
perception layer, reconciliation, code emitter, process-log
discipline, coverage tests.

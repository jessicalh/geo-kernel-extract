# Evil-String Audit — `Atom::pdb_atom_name` (2026-04-28)

**Status (updated late 2026-04-28).** Inventory data is **authoritative**.
Phasing recommendation in §8 is **superseded** by
`spec/plan/openai-5.5-strong-architecture-layout.md` (the active
architecture). The shape of the fix changed from "boundary discipline
+ validators as a standalone session" to "the bug is dissolved by the
`LegacyAmberTopology` + calculator-contract architecture, which types
the implicit-AMBER assumption rather than enforcing it." The 28
consumer sites and 4 loader-side write sites listed below feed the
calculator × topology-field matrix that runs as pre-spec before any
LegacyAmberTopology implementation. Read the 5.5 spec for the active
phasing; read this document for the inventory data and gotcha list.

**Audit method.** Codebase exploration agent (2026-04-28), acting on
OpenAI external-review findings + Claude internal analysis from the
prior session. Counted, classified, and traced every read site of
`Atom::pdb_atom_name` in `src/`, `python/`, `h5-reader/`, and `ui/`.

**Headline numbers.** 28 consumer sites (not the previously-estimated
50). 4 loaders, 2 of which translate via NamingRegistry; the other 2
assign raw. NamingRegistry already has the translation rules; the
fix is largely calling them in two more places.

---

## 1. Consumer inventory by criticality

### Critical — silent corruption sites (charge / protonation / topology)

#### ChargeSource — ff14SB charge lookup
| File | Line | Function | Looks for | Failure on mismatch |
|---|---|---|---|---|
| `src/ChargeSource.cpp` | 48 | `ParamFileChargeSource::LoadCharges` (primary key) | ff14SB keys: "ALA CA", "GLY HA2", etc. | Silent 0.0 charge fallback (KNOWN_BUGS Bug 1) |
| `src/ChargeSource.cpp` | 56 | Same (residue-name fallback) | Same | Silent 0.0 charge fallback |

#### Protein backbone caching and ring detection
| File | Line | Function | Looks for | Failure on mismatch |
|---|---|---|---|---|
| `src/Protein.cpp` | 194 | `DetectAromaticRings` (name→index map) | Ring atom names from `AminoAcidType.rings` | Silent missed rings |
| `src/Protein.cpp` | 230 | `DetectAromaticRings` (HIS tautomer) | "HD1", "HE2" | Wrong HIS ring type (HID/HIE/HIP) |
| `src/Protein.cpp` | 288 | `CacheResidueBackboneIndices` | "N", "CA", "C", "O", "H", "HN", "HA", "HA2", "CB" | Backbone indices unset → cascade failures |
| `src/Protein.cpp` | 306 | `CacheResidueBackboneIndices` (chi angles) | Per-residue chi-angle atom names | Chi angles unmapped → sidechain geometry broken |

#### Protonation variant detection
| File | Line | Function | Looks for | Failure on mismatch |
|---|---|---|---|---|
| `src/ProtonationDetectionResult.cpp` | 50 | Variant detection (name→index map) | All hydrogen names per residue | Silent failed detection |
| `src/ProtonationDetectionResult.cpp` | 67 | HIS δ/ε tautomer | "HD1", "HE2" | Wrong variant (HID/HIE/HIP) → wrong charges |
| `src/ProtonationDetectionResult.cpp` | 91 | ASP protonation | "HD2" | Wrong variant (ASH vs ASP) |
| `src/ProtonationDetectionResult.cpp` | 103 | GLU protonation | "HE2" | Wrong variant (GLH vs GLU) |
| `src/ProtonationDetectionResult.cpp` | 116 | CYS disulfide | "SG" | Cannot detect disulfide |
| `src/ProtonationDetectionResult.cpp` | 123 | CYS thiol | "HG" | Wrong variant (free vs disulfide) |
| `src/ProtonationDetectionResult.cpp` | 133 | LYS amine | "HZ1", "HZ2", "HZ3" | Wrong variant (LYN vs LYS) → wrong charges |

### Cosmetic / diagnostic — display-only, not correctness

| File | Line | Use |
|---|---|---|
| `src/DsspResult.cpp` | 57 | PDB column-width formatting (size check only) |
| `src/KamlProtonator.cpp` | 42 | PDB formatting |
| `src/PropkaProtonator.cpp` | 45 | PDB formatting |
| `src/AIMNet2Result.cpp` | 169 | Error log message |
| `ui/src/RestServer.cpp` | 507 | JSON API display field |
| `ui/src/BackboneRibbonOverlay.cpp` | 131 | VTK type-display label |
| `ui/src/MainWindow.cpp` | 962, 1297 | Atom inspector display |
| `ui/src/ComputeWorker.cpp` | 305 | Documentary comment |

These read the field but do not dispatch on it. Cosmetic mismatch
when input convention differs; not a correctness bug.

### Export surface

| File | Line | Use |
|---|---|---|
| `src/TrajectoryProtein.cpp` | 206 | H5 `/atoms/pdb_atom_name` dataset (raw passthrough) |

The H5 contract is whatever convention the loader produced.
Downstream Python SDK does not currently read this; h5-reader
acknowledges the inconsistency in `notes/nmr_forensics/`.

---

## 2. Loader-side write sites

| Loader | Write site | Translation discipline |
|---|---|---|
| `src/FullSystemReader.cpp` | line 705 | `NamingRegistry::TranslateAtomName(..., Charmm, Standard)` ✓ |
| `src/GromacsEnsembleLoader.cpp` | line 337 | Same ✓ |
| `src/GromacsEnsembleLoader.cpp` | line 517 | Same ✓ (second build path) |
| `src/OrcaRunLoader.cpp` | line 204 | **Raw assignment ✗** |
| `src/PdbFileReader.cpp` | line 131 | **Raw cif++ label ✗** |

`residue_type` is correctly canonicalised across all four loaders
(via `AminoAcidType` resolution). `chain_id` is not translated;
no known issues. The bug is isolated to atom names.

---

## 3. NamingRegistry — capability vs. need

**Active rules (`src/NamingRegistry.cpp:136-308`).** ~22 rules
covering Standard ↔ CHARMM atom-name translation: backbone amide
H/HN, β-methylenes HB2/HB3 ↔ HB1/HB2, γ-methylenes HG2/HG3 ↔
HG1/HG2, plus residue-name variants HIS→HID/HIE/HIP, ASP→ASH,
GLU→GLH, CYS→CYX/CYM.

**Deferred rules (commented-out `lines 172-307`).** ~60 rules
covering CHARMM ↔ IUPAC for ILE quirks, δ/ε/α-methylenes across
ARG/LYS/PRO/GLY, plus the wildcard β-methylene over-fire on ALA.
**Not needed for the evil-string fix** — that fix moves
non-AMBER input *to AMBER*, not to IUPAC. The deferred rules are
for the future IUPACAnnotation projection at output boundary.

**What the fix needs.** Two new translation call sites:

- `OrcaRunLoader.cpp:204` — `NamingRegistry::TranslateAtomName(raw, residue, Orca, Standard)`. ORCA emits standard PDB names so the registered Orca tool context already covers it.
- `PdbFileReader.cpp:131` — needs design decision (see §6 below): assume Standard input and fail loudly on unknowns (most conservative), or detect source convention from atoms present.

NamingRegistry rejects unknown names — that's the right behaviour
(loud failure beats silent corruption) and means the boundary
becomes a strict gate.

---

## 4. Cross-protein matching surface

`src/MutationDeltaResult.cpp` uses KD-tree position matching at
0.5Å tolerance. **No atom-name dependency.** The dispatch fix is
decoupled from the future IUPACAnnotation typed-matcher migration.

`OperationRunner` compares residue counts only. No other
cross-protein atom-identity code in `src/`.

---

## 5. Export contract

| Surface | Reads `pdb_atom_name`? | Convention dependency |
|---|---|---|
| `src/ConformationResult.cpp` NPY identity | No (writes pos/element/residue_index only) | None |
| `src/TrajectoryProtein.cpp` H5 `/atoms/pdb_atom_name` | Yes (raw passthrough) | Loader-dependent |
| `python/nmr_extract/_catalog.py` SDK catalog | No entry | None |
| `python/nmr_extract/_protein.py` SDK reader | Reads H5 group but not for calculation | Cosmetic only |
| `h5-reader/notes/nmr_forensics/` | Compensates via workaround dicts | Aware of inconsistency |
| `ui/` viewer | Display only | Cosmetic only |

After the boundary fix, the H5 contract becomes "AMBER convention
guaranteed" — Python SDK and h5-reader can rely on it without
workarounds. The change is additive (downstream gains a guarantee
they did not have before).

---

## 6. Test coverage today

### Tests that exercise atom-name code paths

| Test | Coverage of bug |
|---|---|
| `test_pdb_loading.cpp` | Loads 1UBQ (reduce-protonated, AMBER input). Silent pass — convention happens to match. |
| `test_protonation_detection.cpp` | 1UBQ + protonated GMX. Silent pass — input is AMBER. |
| `test_naming_registry.cpp` | Tests translation rules explicitly. Passes; but the loaders that *don't* call it aren't tested. |
| `test_mutation_delta.cpp` | KD-tree on positions. Convention-independent. |
| `test_apbs_ff14sb.cpp` | ChargeSource on AMBER-convention 1UBQ. Silent pass. |
| `test_full_pipeline.cpp` | End-to-end on AMBER-convention input. Silent pass. |

### Silent-corruption modes with NO test coverage

1. ChargeSource with non-AMBER input (Bug 1 acute symptom)
2. ProtonationDetectionResult with CHARMM-named hydrogens
3. Backbone caching / ring detection with non-AMBER input
4. OrcaRunLoader's raw-assignment path (no test of OrcaRunLoader at all)
5. PdbFileReader with user-supplied non-AMBER PDB

**Implication.** Every existing test passes because every input is
AMBER-convention by accident. The boundary discipline is undefended.
The new session must add tests that exercise non-AMBER inputs and
verify either correct translation or loud failure.

---

## 7. Things the prior analyses missed

- **Real count is 28, not 50.** The "approximately 50" estimate from
  the 2026-04-27 audit was inflated; the actual consumer surface is
  smaller and more bounded.
- **`AminoAcidType::atoms` (`vector<AminoAcidAtom>` with `const char*` name)** is the *reference* AMBER-convention list; consumers compare *against* it. Not a consumer, not a producer of bug; it is the canonical source.
- **`AminoAcidRing::atom_names` (`vector<const char*>`)** participates in the dispatch indirectly via `DetectAromaticRings` building a name→index map. Counted in the consumer list (Protein.cpp:194,202).
- **No substring/pattern matching on atom names anywhere.** Only equality and size checks. The migration is mechanical, not semantic.
- **Reduce output is AMBER-convention.** The `--pdb` path's passthrough is correct *for reduce output* but undefended for user-supplied PDB on `--protonated-pdb`.
- **Cross-protein matching is decoupled.** MutationDeltaResult is geometry-based; the IUPACAnnotation typed-matcher migration is a separate concern from the evil-string fix.
- **Python SDK and h5-reader are downstream-only consumers** that don't read `pdb_atom_name` for calculation. The H5 contract is the only export surface that matters for the fix.
- **Source-determination problem for PdbFileReader.** User-supplied PDB could be in any convention. Most conservative path: assume Standard, let NamingRegistry fail loudly on unknowns. Alternative: detect convention from atoms present. This is a real design decision the new session must make.

---

## 8. Phasing recommendation (informs the new-session prompt)

### Phase 1 — Boundary enforcement + test coverage

1. **OrcaRunLoader.cpp:204** — add `NamingRegistry::TranslateAtomName(..., Orca, Standard)` call.
2. **PdbFileReader.cpp:131** — add translation call; design decision needed on source-context.
3. **ChargeSource.cpp:72** — replace silent 0.0 fallback with loud diagnostic logging the atom index, residue, and name.
4. **Test coverage** for all four loaders × representative dispatch sites × non-AMBER inputs. Make silent-corruption modes loud failures.
5. **Comment audit** of FullSystemReader/GromacsEnsembleLoader explaining the translation they already do.

### Phase 2 — Typed source provenance (TBD)

Three options:
- **A.** Enum `AtomSourceContext` field on `Atom`.
- **B.** Parallel array on `Protein` mapping atom_index → source.
- **C.** String invariant enforced at boundary; no source tracking needed.

Recommendation: **C** is simplest; A/B add diagnostic value but
risk recreating the IUPAC-on-Atom mistake. Decide based on whether
diagnostic surfaces need to know which loader produced each atom.

### Phase 3 — Establish AMBER invariant

- Document in `spec/CONSTITUTION.md`: post-`FinalizeConstruction`, all `Atom::pdb_atom_name` are AMBER convention.
- Mark each of the 28 consumer sites with a comment referencing the invariant.
- Optional: a validator that walks the protein post-construction and asserts the invariant.

### Phase 4 — IUPACAnnotation prerequisite met

After the above, IUPACAnnotation lands on a foundation where the
existing string surface is well-defined. The plan in memory entry
`project_proteintopology_architecture` then proceeds.

---

## 9. Canary test for "done"

```cpp
TEST(AtomNameBoundary, AllLoadersProduceAMBERConvention) {
    // Load the same protein via all four loaders with inputs in
    // different conventions where applicable (reduce output, GROMACS
    // CHARMM TPR, ORCA output, user PDB).
    // For each loaded protein:
    //   1. Walk every atom.
    //   2. Verify NamingRegistry::IsValidAMBERName(atom.pdb_atom_name,
    //      residue) is true, OR the loader returned an error result.
    //   3. Spot-check a handful of known atoms (CA, HA, HD1 on a HIS,
    //      HD2 on an ASP) for AMBER spelling.
    // Expectation: all loaders produce AMBER-convention names, OR
    // refuse loudly. No silent passes through to the consumer code.
}
```

When this test passes (and the existing silent-corruption tests
above are added and pass), the evil string is gone.

---

## Summary

| Aspect | Number / state |
|---|---|
| Consumer sites | 28 (not 50) |
| Critical-correctness sites | 13 (charge + backbone + protonation) |
| Cosmetic / display sites | 8 |
| Export-surface sites | 1 (H5) |
| Loaders writing | 4 (2 correct, 2 raw) |
| Cross-protein matching affected | 0 (geometry-based) |
| NamingRegistry rules needed but missing | 0 (Standard ↔ CHARMM is in place; just call it) |
| Test gaps (silent-corruption modes uncovered) | 5 |
| Estimated session count | 2-3 focused sessions through full test |

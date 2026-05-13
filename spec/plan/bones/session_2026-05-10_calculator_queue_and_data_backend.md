# Session 2026-05-10 — Calculator queue, 1P9J reframe, Larsen suite, tensorcs15 backend

This is the durable record of the 2026-05-10 session that closed several
calculator-queue design questions, locked the data backend for the
tripeptide DFT lookup, narrowed Stage 2 to a single-protein deep study on
1P9J, and uncovered (and corrected) a methodological misnaming around the
serine DFT regeneration. Read end-to-end before resuming the calculator
work; the discoveries are interlocked.

## Session arc

The session opened with a request for context-heavy mode covering the next
20–30 sessions of calculator work, expanded into a structured
walkthrough of the pending-calculator queue, and ended with operational
work on the data backend that unblocks the next set of calculator
landings. The arc:

1. **Context load** — re-read PATTERNS.md, OBJECT_MODEL.md (chunked),
   PLANNED_CALCULATORS_2026-04-22.md (with all 2026-05-09 amendments),
   the planned-calculator-substrate-audit, the I/O schema and diagnostics
   docs, plus McConnellResult, AIMNet2PolarisabilityResult, and
   PlanarGeometryResult as exemplars.
2. **H5 / NPY storage discussion** — five schema decisions surfaced,
   parked for later. Decision: minimal-infrastructure path — git-tag
   the current NPY format as v1, leave existing calculator outputs
   unchanged, defer zip wrapper and SDK changes, only update for science
   reasons.
3. **A.1 CSAPrincipalAxis convention pick** — three decisions confirmed
   (Haeberlen sort, σ-positive-more-shielded sign, first-nonzero-positive
   eigenvector + cross-product PAS for right-handedness). "Equivariance
   tension" clarified as no-tradeoff: keep both sorted-PAS and
   equivariant SphericalTensor representations.
4. **1P9J reframe** — Stage 2 narrowed to a single-protein deep study.
   T1E, EGF/TGFα chimera, BMRB 5801, Wingens 2003 paper. 15 ns MD at
   20 ps stride = 750 frames. The 260-DFT-pose calibration set is OUT
   of scope.
5. **A.2 Amide Tensor Geometry** — held pending Loth 2005 PDF for
   convention verification. Recommended frame: z along N-H, y normal
   to amide plane, x = y × z. Larsen-style validation framing replaced
   with modest "do our autocorrelation calculators see any of the same
   things as Wingens' Modelfree dynamics" framing.
6. **Markwick 2010 surfaced** — methodological precedent for MD-averaged
   shift prediction on dynamically heterogeneous proteins. Codex was
   right to flag it; cited explicitly in the methods narrative now.
7. **Larsen suite walkthrough** — HBondHα variant (Δσ_HαB) and tripeptide
   BB DFT assembly (σ_BB^i) identified as next-priority queue items.
   ProCS15's six-term decomposition mapped against our pipeline:
   four terms covered analogously, two real gaps.
8. **Web research on Larsen data** — found four additional ERDA archive
   files we don't have. Critically: Δσ_BB^{i±1} doesn't need new DFT
   data — Larsen 2015 Eq 3 reuses the AXA scans. Code-only addition.
9. **Gotham assembler discovery** — `TripeptideAssembler.cpp` and
   `TripeptideDatabase.cpp` exist on gotham at
   `/mnt/extfast/2026thesis/working-reference/src/External/`. Postgres-
   backed (libpq), JSONB tensor parsing, Kabsch alignment, sidechain
   re-rotation via Rodrigues. Real implementation, needs porting not
   from-scratch design.
10. **Database investigation** — tensorcs15 Postgres on gotham; three
    data sources tracked; raw_dft_calculations holds full Mat3 tensors
    pre-decomposed into irreps. chi3/chi4 latent for K/R/E/M/Q.
11. **The serine investigation** — the SER row count discrepancy
    between gotham (6,259) and batcave (0) at first looked like a
    replica gap. Followed user steer to look at original generation
    files. Found `frame_type` column as the method discriminator. Then
    found the original `frame_notes` says PBE, not r²SCAN as I'd been
    documenting. Read the README at the source location:
    "DFT functional: PBE (ORCA) vs OPBE (Gaussian) - similar". Memory
    corrected.
12. **Data backend operational work** — replicated SER subset to batcave,
    granted SELECT to user `jessica`, added [databases] section to
    `~/.nmr_tools.toml`. Verified end-to-end.

## Decisions locked

### NPY / H5 storage — minimal-infrastructure path

- **Tag git for NPY format v1.** No sentinel file, no SDK migration. Old
  thesis NPY data and new runs use the same SDK reader.
- **Existing calculator outputs stay unchanged.** Only update for
  science reasons (sign convention turning out wrong, units drift, etc.).
- **Per-frame NPY tree continues as today.** Trajectory analyses
  (autocorrelation, J(ω), memory kernel, CCR rate) emit per-TR H5
  groups via the existing per-`WriteH5Group` pattern.
- **Five schema questions parked:** lag-axis layout, irrep attributes,
  sign-convention attributes, conformation-vs-trajectory split,
  conformation-scope tensor metadata. Revisit when we have the
  trajectory chain landed and can see actual emitted shapes.
- **Sharpening waits.** Per-frame NPY already serves as the substrate
  for trajectory calculators that want to scan σ(t); no H5 reorganisation
  needed for the gain.
- **Deferred:** zip wrapper, SDK zip-aware reader, HDF5 deflate filter,
  conversion tool. Not blocking calculator work.

### A.1 CSAPrincipalAxisResult — three convention picks

1. **Sort: Haeberlen.** σ₃₃ is the principal component furthest from
   σ_iso, σ₁₁ next, σ₂₂ closest. Matches codebase irrep formalism,
   solid-state literature reporting, standard Δσ / η constructions.
   Field name carries the convention: `csa_principal_components_haeberlen`.
2. **Sign: σ-positive-more-shielded.** The σ tensor convention, not
   δ. Matches existing `*_shielding_*` storage. δ projection is
   consumer-side downstream of this calculator.
3. **Eigenvector sign-fixing: first-nonzero-positive on axes 1 and 2,
   axis 3 = axis_1 × axis_2** (right-handed PAS by construction).
   Deterministic per frame from σ tensor alone. Frame-coherent
   continuity across a trajectory deferred to a separate
   `CSAFrameCoherentTrajectoryResult` slice.
4. **Storage: both representations.** Sorted PAS (Vec3 components +
   Mat3 axes + Δσ + η) for literature-comparison consumers; existing
   equivariant Mat3 + SphericalTensor on the atom for ML. ~17 doubles
   storage cost per atom — negligible. Equivariant ML reads
   `*_shielding_contribution`; sorted PAS field name discloses convention.

### A.2 AmideTensorGeometryResult — held

- Held pending Loth 2005 PDF acquisition. The L11 Brender-Taylor-
  Ramamoorthy 2001 convention reference isn't in `references/`; the
  frame axis labelling is best verified against actual figures.
- Working assumption (NOT locked): z along N→H, y normal to amide
  plane (computed as N→C'(i-1) × N→HN), x = y × z. Right-handed.
  Three-atom plane definition (C'(i-1), N(i), HN(i)).
- Outputs three local-frame tensors per amide (σ on N, on HN, on
  C'(i-1)) plus the lab→local rotation matrix.
- NaN-fill for residues without a peptide-amide bond into them
  (residue 0, prolines).
- When Loth lands, verify axis labelling and confirm or revise.

## Key findings

### 1P9J is fundamentally different from CSA-bench proteins

- Small (36 residues), EGF-fold (mostly β + loops, limited α-helix),
  three disulfide bonds in tight fold.
- High flexibility: ensemble RMSD 2.43 ± 0.63 Å, S² down to 0.2 at
  termini, ~40% residues need R_ex term in Modelfree.
- **No published per-residue CSA tensors** — Wingens 2003 reported
  shifts (BMRB 5801), R1, R1ρ, hetNOE, S², τ_e, R_ex; not σ principal
  components.
- Citable comparisons available: BMRB shift RMSE, Yao-Bax β-sheet
  ¹⁵N CSA partial test (limited α makes the full split test
  weak), σ_N variance vs R_ex correlation (potential thesis-grade
  finding), MD-derived S² vs Modelfree S², ensemble-vs-ensemble
  consistency.
- Comparisons NOT available: per-bond Loth 2005 ubiquitin CCR table,
  per-residue Wylie 2011 solid-state CSA, Wylie 2006 ¹³C' CSA bench.
- Markwick 2010 is the methodological precedent: AMD + ensemble
  averaging cuts SHIFTX RMSD by 28% on flexible IκBα; the gain is
  largest on residues with R_ex; ¹⁵N is the biggest dynamics
  reporter (2.89 → 1.84 ppm). Our plain 15 ns MD won't reach
  μs-ms exchange — methodological caveat for R_ex residues that
  the methods text declares.

### Larsen ProCS15 six-term decomposition vs our pipeline

| Larsen term | Our coverage | Notes |
|---|---|---|
| σ_BB^i | Gotham TripeptideAssembler exists; needs port + chi3/chi4 extension | T2-preserving via raw Gaussian logs in DB |
| Δσ_BB^{i±1} | Not built | Larsen Eq 3 reuses the AXA scans — code-only addition, no new DFT |
| Δσ_HB^i | HBondResult kernel × η_elem | Analogous goal, kernel form vs Larsen-grid form |
| Δσ_HαB^i | Not built (designed) | Kernel × η form (Larsen-grid HB data not in our DB) |
| Δσ_RC^i | BS + HM + RingSus T2 kernels | Exceeds Larsen's point-dipole (full tensor vs isotropic) |
| Δσ_w^i | WaterField + HydrationShell kernels | Larsen's is a single 2.07 ppm scalar; our kernels are richer |

Plus methodology pieces: linear δ = b − aσ regression (downstream in
`learn/`), GESD outlier diagnostic (Python post-pass), missing-grid-point
handling (need to mirror Larsen's SciPy linear interpolation discipline
when porting).

### Δσ_BB^{i±1} = AXA reuse — load-bearing finding

Larsen 2015 Eq 3:

    Δσ_BB^{i−1}(residue i) = σ_BB^{i−1}(φ_{i−1}, ψ_{i−1}, χ_{i−1}, …)
                              − σ_A(φ_std, ψ_std)

with φ_std = −120°, ψ_std = 140° (the AAA reference geometry).

The neighbor sidechain effect on residue i is approximated by the
identity-vs-alanine shift on the **flanking residue's own atoms**,
evaluated in a tripeptide centered on i−1 with X = i−1's identity.
The σ_A baseline is one fixed AAA tripeptide entry at standard backbone.

**Practical consequence:** the `TripeptideNeighborShieldingResult`
(Δσ_BB^{i±1}) is a code-only addition to the same `tensorcs15` DB.
Same JSONB tensor structure, same Kabsch + tensor rotation pattern as
σ_BB^i. Different lookup key construction:

- For Δσ_BB^{i−1}: query `raw_dft_calculations` row where
  `central_residue` = (residue i−1's 1-letter code), `phi`/`psi`/`chi*`
  = i−1's actual angles. Read σ tensor at the central residue's atoms
  in that row. Subtract the AAA tripeptide row at φ=−120°, ψ=140°.
- For Δσ_BB^{i+1}: same shape, with i+1's identity and angles.

The strong assumption Larsen makes (i±1 effect on i ≈ identity-vs-Ala
shift on i−1's own atoms) is acknowledged but not blocking — Larsen
demonstrated it works.

### ERDA archive — three unfetched files (one pulled this session)

**Pulled 2026-05-10:** two files at `/mnt/expansion/larsen_archive/`,
provenance README alongside.

**`hydrogenbondnmrlogs.tar`** (503 MB, SHA256 `dbca1e85...`).

Inspecting the tar contents revealed Larsen's H-bond
parameterisation is **richer than my walkthrough had assumed** — 6
nested `.tar.bz2` archives covering 2 donor classes × 3 acceptor
classes:

- NMA / NMA (canonical amide-amide HB)
- Ac-A-NMe (Hα) / NMA (Δσ_HαB amide acceptor)
- Hα / carboxylate, Hα / hydroxyl
- amide / carboxylate, amide / hydroxyl

**Methodological implication:** Larsen distinguished acceptor types
(carbonyl O / carboxylate O / hydroxyl O) at parameterisation level,
not just element. Our existing `HBondResult` uses element-based
acceptor identification and doesn't differentiate. The Larsen-grid
form is a richer methodological substrate than our kernel × η_elem.
Worth revisiting `HBondResult` and `HBondHalphaResult` design when
we wire them against this data.

Tarball not yet ingested into `tensorcs15`. Future ingestion is its
own slice.

**`structures.tar.bz2`** (36 KB, SHA256 `ddea08df...`). Tiny because
it's just three PDB files: `1UBQ.pdb` (Larsen's ubiquitin crystal),
`1UBQ_pm6dh3plus.pdb` (PM6-D3H+ optimized — Larsen's actual NMR-input
geometry), `2OED_pm6dh3plus.pdb` (GB3, the other validation target).
**Load-bearing for porting validation:**

- Larsen 2015's published RMSD numbers were computed against these
  exact geometries. Running our ported `TripeptideBackboneShieldingResult`
  on `1UBQ_pm6dh3plus.pdb` and comparing per-atom σ_BB to Larsen's
  published values is the direct reproduction check.
- The gotham assembler's reference output on the same PDBs is the
  byte-level golden test for our port (FP tolerance, same DB lookups).
- GB3 is also L7 Yao-Bax 2010's ¹⁵N CSA bench protein — direct anchor
  for the α/β CSA split test.

This is the kind of pre-port discovery that justifies acquiring
`structures.tar.bz2` before committing to the port: it gives us a
test fixture and a published-RMSD reproduction target.

### ERDA archive — two unfetched files (original list of four reduced to two)

Public ERDA archive at
https://www.erda.dk/public/archives/YXJjaGl2ZS1TYk40VXo=/published-archive.html
holds Larsen's 2015 ProCS15 deposit (Lars Bratholm, lars.bratholm@chem.ku.dk).
We have three of seven; four unfetched:

- ⏳ **`hydrogenbondnmrlogs.tar`** — Gaussian NMR logs from H-bond scans
  on N-methylacetamide-dimer + Ac-A-NMe models. **Covers both Δσ_HB
  (amide) and Δσ_HαB (Cα-H) Larsen-grid parameterisation.** Most
  useful unfetched file. Worth as Larsen-grid cross-validation.
- ⏳ `structures.tar.bz2` — 1UBQ crystal + PM6-D3H+ optimized geometries.
- ⏳ `proteinnmrlogs.tar.bz2` — protein-level DFT validation runs.
- ⏳ `predictions.tar.bz2` — ProCS + empirical predictor outputs on
  test proteins.

Δσ_w needs no file: a single 2.07 ppm scalar.

GitHub: https://github.com/jensengroup/procs15 — CMake/C++/Phaistos
plugin, points at the same ERDA archive for data.

### tensorcs15 schema is richer than we'd been treating

Three data sources tracked in `data_sources` table:

1. `ProCSnumpyfiles.tar.bz2` (3.0 GB) — Larsen's published isotropic
   numpy lookup tables → `procs15_arrays` + `numpy_array_values`
2. `shieldings.tar` (318 MB) — parsed isotropic intermediate values →
   `parsed_isotropic_values`
3. `nmrlogs.tar` (828 MB compressed → ~30 GB uncompressed) — raw
   Gaussian NMR logs with full T2 tensors → `raw_dft_calculations`

**Per-tensor JSONB** in `raw_dft_calculations.tensors`:

```json
{
  "element": "C", "atom_idx": 1,
  "isotropic": 43.78, "anisotropy": 92.37,
  "tensor_3x3": [[88.07, -33.29, -22.59], [-25.53, -17.88, 1.48],
                 [-22.16, -10.70, 61.15]],
  "eigenvalues": [-26.69, 52.68, 105.37],
  "t2_components": [105.95, -29.41, -22.38, -4.61, 21.27]
}
```

The Mat3 is **asymmetric** (full T0 + T1 + T2 form, not symmetrized).
The irrep decomposition (isotropic, anisotropy, eigenvalues, 5-component
T2) is **already baked in at ingest** — we don't decompose at runtime,
we read `t2_components` directly from JSONB.

### Per-residue counts and chi-coverage

| Residue | n_calcs | chi3 set | chi4 set | n_tensors_max |
|---|---|---|---|---|
| R | 285,028 | 285K | 285K | 56 |
| K | 145,598 | 146K | 146K | 54 |
| E | 134,013 | 133K | 0 | 47 |
| M | 129,661 | 130K | 0 | 49 |
| Q | 130,681 | 131K | 0 | 49 |
| H/L/D/N/T/F/W/Y/I/D | ~100K each | 0 | 0 | 46–56 |
| ALA | 31,037 | 0 | 0 | 42 |
| **S (PBE)** | **6,259** | 0 | 0 | 43 |
| V/C | ~6K each | 0 | 0 | 43–48 |
| G | 344 | 0 | 0 | 39 |
| P | 258 | 0 | 0 | 46 |

**The gotham assembler queries (chi1, chi2) only and leaves chi3/chi4
data on the floor for K, R, E, M, Q.** Direct fix when porting:
extend the lookup to use chi3/chi4 columns when set.

### The serine PBE discontinuity (corrected from r²SCAN)

Larsen's unpublished T2 stash lost the ASA (Ala-Ser-Ala) Gaussian
logs. The project regenerated SER tripeptide DFT with **ORCA PBE**,
not r²SCAN as initial session notes claimed. Confirmed via:

1. README at
   `batcave:/mnt/expansion/williamsproject/experiments/exp00_ubiquitin_exploration/asa_regeneration/README.md`:
   > "DFT functional: PBE (ORCA) vs OPBE (Gaussian) - similar"
2. Every `asa_nmr.inp` first ORCA directive: `! PBE 6-31G(d,p) NMR CPCM(water)`

PBE is a pure GGA. OPBE is also a pure GGA (OPTX exchange + PBE
correlation). r²SCAN is a meta-GGA. The PBE-OPBE offset is much smaller
than r²SCAN-OPBE would have been. Per-(residue, atom_type, method)
calibration absorbs the offset.

Additional caveat: SER coordinates from **fragbuilder via tleap**
(force-field placement), while Larsen used **PM6 optimization with
frozen torsions**. Two small SER-vs-rest differences (functional and
geometry), both ride on the same `frame_type` discriminator.

### `frame_type` is the runtime DFT-method discriminator

Built into the `raw_dft_calculations` schema:

- `frame_type = 'gaussian_standard_orientation'` — Larsen's OPBE 19
  residues (1,864,372 rows)
- `frame_type = 'orca_input_orientation'` — SER PBE re-run (6,259 rows)

Calibration code reads `frame_type` per row and routes by method.
`parser_version` is a redundant discriminator. Methods text declares
explicitly.

## Operational state changes

### Database

- **batcave `tensorcs15` replica caught up with gotham** —
  pulled 6,259 SER rows via `COPY ... TO STDOUT` over SSH, applied
  via `COPY ... FROM STDIN` in a transaction. Both DBs now at
  1,870,631 rows. calc_id range for SER: 1,892,717–1,898,975 (no
  primary-key collision; SER calc_ids are above batcave's previous max).
- **GRANT applied** — user `jessica` has `CONNECT`, `USAGE`, and
  `SELECT` on schema `public`. `ALTER DEFAULT PRIVILEGES` ensures
  any future tables auto-grant.

### Configuration

`/home/jessica/.nmr_tools.toml` extended with:

```toml
[databases]
tensorcs15 = "host=/var/run/postgresql dbname=tensorcs15 user=jessica"
```

with comment block explaining `frame_type` discriminator.

### Verification

`psql "host=/var/run/postgresql dbname=tensorcs15 user=jessica" -c
"SELECT COUNT(*) FROM raw_dft_calculations"` returns 1,870,631 — full
DB readable as user `jessica`.

## Memory entries written / updated

Eight memory files, all in
`/home/jessica/.claude/projects/-shared-2026Thesis-nmr-shielding/memory/`:

1. **`project_1p9j_study_system.md`** — 1P9J replaces 260-pose; Wingens
   2003 paper acquired; 750 frames at 20 ps; supersedes Stage 2
   narrative in `project_three_stages.md`.
2. **`project_post_csa_calculator_queue.md`** — HBondHα variant and
   tripeptide BB DFT assembly are next-priority after A.1 + A.2.
   Half-sphere covered by existing proxies.
3. **`project_hbond_halpha_design.md`** — kernel × η_Hα form (NOT
   Larsen-grid; Larsen's HB grid data is not in our DB). Substrate
   gates, geometric criteria, tensor form, calibration shape.
4. **`reference_gotham_assembler.md`** — tensorcs15 Postgres on gotham
   is the authoritative tripeptide DFT source; full T2 tensors with
   pre-decomposed irreps; chi3/chi4 latent for K/R/E/M/Q;
   `frame_type` field is the method discriminator.
5. **`project_serine_pbe_discontinuity.md`** — replaces stale
   r²SCAN entry. ORCA PBE/6-31G(d,p) GIAO CPCM(water) chosen as
   closest pure-GGA approximation to Larsen's OPBE; geometry from
   fragbuilder via tleap; smaller-than-r²SCAN-would-be discontinuity;
   per-(residue, atom_type, method) calibration absorbs the offset.
6. **`project_larsen_neighbor_axa_reuse.md`** — Larsen 2015 Eq 3 reuses
   AXA scans for Δσ_BB^{i±1}; code-only addition, no new DFT data;
   φ_std = −120°, ψ_std = 140°; SER caveat applies via frame_type.
7. **`reference_erda_archive_missing_files.md`** — ERDA archive URL,
   four unfetched files, hydrogenbondnmrlogs.tar flagged as next
   download priority for HB / HαB Larsen-grid cross-validation.
8. **`project_first_ml_combination.md`** — OF3 embeddings + MACE +
   Larsen calculator suite = first ML; calculator queue is the gate;
   calculator outputs are ML feature inputs not just literature
   comparisons; consumer-side contract changes the design discipline
   for NPY shapes and metadata.

`MEMORY.md` index updated with one-line pointers to each.

Plus stale entry removed: `project_serine_rscan_discontinuity.md`
(replaced by the PBE entry above).

## Calculator queue — current state

```
✓ A.1 CSAPrincipalAxis              ready to draft (3 conventions confirmed)
⏸ A.2 AmideTensorGeometry           held pending Loth 2005 PDF
➤ TripeptideBackboneShielding       port from gotham + chi3/chi4 extension
➤ TripeptideNeighborShielding       Δσ_BB^{i±1}, AXA-reuse, code-only
➤ HBondHalphaResult                 kernel × η_Hα form (revised)
— A.3 BulkSusceptibilityAccumulator long tail
— A.4 NICSProbeEvaluator + probe API long tail
— PCS / PRE pair                    long tail (lanthanide-tag input)
— A.7–A.10 trajectory chain         long tail (autocorrelation, J(ω),
                                    memory kernel, CCR rate)
— A.11 BenchmarkBackCalculation     long tail (consumes all above)
```

Plus, on the trajectory side after the calculator queue:
- `CSAFrameCoherentTrajectoryResult` — frame-coherent realignment of
  PAS axes across a trajectory (stripped from A.1 to keep that
  calculator data-only).

Plus, beyond the calculator queue:
- Acquire `hydrogenbondnmrlogs.tar` from ERDA for Larsen-grid HB / HαB
  cross-validation reference data.
- Acquire Loth 2005 PDF to settle A.2 amide-frame convention.

## Open questions

1. **A.2 amide-frame convention** — held until Loth 2005 PDF acquired.
2. **Calibration composition for σ_BB + our kernels** — additive vs
   replace-local-kernels vs subtract-our-kernels-on-tripeptide-model.
   Defer to calibration time; calculator emits σ_BB regardless.
3. **PBE-OPBE offset characterisation** — optionally run a small
   overlap set (~100 shared geometries) with both Gaussian PBE and
   ORCA PBE to characterise the systematic offset. Useful methods-text
   material; not blocking.
4. **Five H5 schema questions** — lag-axis layout, irrep attributes,
   sign-convention attributes, conformation-vs-trajectory split,
   conformation-scope tensor metadata. Re-open when trajectory chain
   is landed and we can see emitted shapes.
5. **Δσ_HαB grid form vs kernel form** — current decision: kernel × η
   form. If `hydrogenbondnmrlogs.tar` is acquired, re-evaluate whether
   to add a Larsen-grid path as cross-check or replacement.

## Next-session jumping-off point

The data backend is unblocked. The natural first action of the next
session is the C++ porting of the gotham `TripeptideAssembler` adapted
to the current project's discipline:

1. **Create `src/TripeptideDftTable.{h,cpp}`** — libpq-backed loader.
   At session init, query `raw_dft_calculations` SELECT all columns,
   parse JSONB tensors, build a per-residue `nanoflann` kd-tree over
   `(phi, psi, chi1, chi2, chi3, chi4)` keys. Handles chi3/chi4
   extension (gotham assembler's gap). Stash `frame_type` per record.
2. **Create `src/TripeptideBackboneShieldingResult.{h,cpp}`** —
   `ConformationResult` subclass. Per residue: get φ/ψ from
   `DsspResult`, χ from `Residue::chi[k]`, look up nearest record,
   Kabsch-align tripeptide N/CA/C onto protein N/CA/C, rotate
   tensors by R, store on per-atom `tripeptide_bb_sigma_tensor` /
   `_spherical` / `_match_distance` fields. NaN-fill non-backbone
   atoms. Substrate gates via `BackboneRole`.
3. **Create `src/TripeptideNeighborShieldingResult.{h,cpp}`** —
   companion calculator using same `TripeptideDftTable`. Different
   key construction: look up flanking residue's identity + angles,
   subtract AAA standard-angle baseline, accumulate per Larsen Eq 3.
4. **Add `Session::tripeptide_table_`** as the canonical loader
   ownership site. `Session::LoadFromToml` reads
   `[databases].tensorcs15` connection string and constructs the
   table at startup.
5. **NPY emission for both calculators** — `tripeptide_bb_sigma.npy`,
   `tripeptide_bb_match_distance.npy`,
   `tripeptide_neighbor_sigma.npy`. Atom-major float64. NaN-fill for
   inapplicable atoms.
6. **Tests** — golden-value test against the gotham assembler's
   output on a single residue (validate Kabsch + tensor rotation
   produces identical numbers). Plus a smoke test on 1P9J to confirm
   end-to-end flow.

The gotham assembler at
`/mnt/extfast/2026thesis/working-reference/src/External/TripeptideAssembler.cpp`
+ `TripeptideDatabase.cpp` is the porting reference. Read first; adapt
to current project's `Protein::AtomAt(i).element`, `conf.PositionAt(i)`,
substrate-driven backbone caches, `ConformationResult` discipline, and
the project's `OperationLog` / `KernelFilterSet` patterns.

The frame-type discipline carries through: when looking up a residue
that hits an `orca_input_orientation` row, the same Kabsch +
tensor-rotation logic applies (the math is frame-agnostic; positions
and tensors are in the same frame within each row). The
`frame_type` value is stashed alongside each lookup result for
downstream calibration code.

## Things flagged for the user to decide / confirm next session

- Whether to acquire `hydrogenbondnmrlogs.tar` immediately (for
  Larsen-grid HB cross-validation) or defer until calculator chain
  lands and we want the cross-check.
- Whether to acquire Loth 2005 PDF immediately to unblock A.2 amide,
  or proceed with porting σ_BB / Δσ_BB^{i±1} first since they
  don't depend on A.2.
- Whether the next session starts with the C++ port (concrete code
  work) or with another walkthrough item (A.3 BulkSus, A.4 NICS+probe
  API, PCS/PRE pair).

The first-ML framing (OF3 embeddings + MACE + Larsen suite) means
porting the σ_BB / Δσ_BB^{i±1} calculators is the most direct
unblock for the first ML run. Recommended path: code first, walkthroughs
second.

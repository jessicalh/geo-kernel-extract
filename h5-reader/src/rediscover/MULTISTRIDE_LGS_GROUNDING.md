# Multi-Stride .LGS Grounding

Read-only grounding for a future calcset wrapper that points at the existing
1P9J ORCA DFT jobs while using the dense all-frame MOPAC extraction at:

`/shared/2026Thesis/shielding-calcsets/data/trajectories/1p9j-calibration-with-dft/extract-full-mopac-allframes-20260531`

No design decision is made here.

## 1. `.LGS` Schema And Loader

The authoritative `.LGS` schema is at:

`/shared/2026Thesis/nmr-shielding/spec/CALCSET_MANIFEST.md`

The h5-reader README points at `spec/CALCSET_MANIFEST.md`, but this branch does
not have that file under `h5-reader/spec/`; the shared top-level spec above is
present.

The current h5-reader loader lives in:

- `src/io/CalcsetManifest.h`
- `src/io/CalcsetManifest.cpp`

The current writer/helper lives in:

- `tools/lgs_write.py`

So the loader is in h5-reader, not in the producer/extractor. The writer is also
currently a h5-reader tool. The extractor writes `trajectory.h5`, NPYs, PDBs,
and `extraction_manifest.json`; it does not currently own `.LGS` generation.

Schema v1 top-level fields:

- `schema_version`: integer, must be `1`.
- `kind`: `"trajectory"`, `"single_pose"`, or `"mutant_pair"`.
- `dataset_id`: calcset id.
- `protein_id`: authoritative downstream protein id.
- `human_name`: display name.
- `trajectory`: required when `kind == "trajectory"`.
- `single_pose`: required when `kind == "single_pose"`.
- `mutant_pair`: required when `kind == "mutant_pair"`.
- `dft`: optional.
- `metadata`: informational.

`trajectory` fields:

- `md_dir`: path to MD directory.
- `topology_top`: path to `topol.top`.
- `extraction_dir`: path to extractor output directory.
- `trajectory_h5`: path to extractor `trajectory.h5`.
- `extraction_manifest`: path to extractor metadata.
- `frame_dt_ps`: positive number.
- `frame_index_basis`: string. Rediscover currently expects
  `"trr_frame_index"`.
- `reference_pdb`: optional path.

`dft` fields:

- `method`: string.
- `campaign_target_frames`: integer.
- `frame_stride`: `{ "first": int, "last": int, "step": int }`.
- `frames`: authoritative array of DFT jobs.
- `frames[].frame_index`: non-negative original trajectory frame index.
- `frames[].meta_json`: path to the ORCA job meta JSON.

Important loader behavior:

- `.LGS` may be loaded by passing either the calcset directory or the `.LGS`
  file path.
- Directory load requires exactly one `*.LGS`.
- Required declared paths must exist and have the expected file/dir kind.
- `.LGS` does not parse/validate artifact contents beyond paths.
- `protein_id` in `.LGS` is authoritative. `QtProteinLoader::LoadFromManifest`
  overrides producer-side values from H5/extraction metadata with the `.LGS`
  value.
- `dft.frame_stride` is not used for coverage or mapping. `dft.frames[]` is
  the authoritative DFT frame list.

Stride support facts:

- The schema has no first-class `extraction_stride`, no `extraction.frames[]`,
  and no explicit DFT-to-extraction-row map.
- The reader maps by original trajectory frame index:
  - H5 `/trajectory/frames/original_index` gives H5 row -> original TRR frame.
  - Per-frame NPY dirs are named `frame_NNNNNN`, where `NNNNNN` is the same
    original TRR frame.
  - `.LGS dft.frames[].frame_index` is the same original TRR frame.
- `TrajectoryFrameMap::ScanSampledRows()` scans `npys/frame_*` and maps those
  original indices back to H5 rows.
- `TrajectoryConformation::originalFrameIndex(row)` reads the H5 frame map.
- `DftShieldingStore` and rediscover `DftFrameSet` key DFT by original index.
- `rediscover::FrameMap::Build()` walks every H5 row and keeps only rows whose
  original index is present in the loaded DFT set.

Consequence: current h5-reader/rediscover does not require equal extraction and
DFT row counts when both sides share `trr_frame_index`. It can represent sparse
DFT targets over dense extraction implicitly through original frame indices.
What `.LGS` lacks is an explicit schema-level declaration of that multi-stride
relationship.

## 2. Dense MOPAC Extraction Directory

Directory:

`/shared/2026Thesis/shielding-calcsets/data/trajectories/1p9j-calibration-with-dft/extract-full-mopac-allframes-20260531`

Top-level contents:

- `atoms_category_info.npy`
- `bonds.npy`
- `command.sh`
- `extraction.log`
- `extraction.pid`
- `extraction_manifest.json`
- `npys/`
- `pdbs/`
- `residues.npy`
- `ring_membership.npy`
- `rings.npy`
- `trajectory.h5`

There is no literal `data/` subdirectory at this path. The per-frame data lives
under `npys/`.

Top-level sidecar shapes:

- `atoms_category_info.npy`: shape `(846,)`
- `bonds.npy`: shape `(862,)`
- `residues.npy`: shape `(54,)`
- `rings.npy`: shape `(16,)`
- `ring_membership.npy`: shape `(96,)`

The top-level `extraction_manifest.json` says:

- `schema_version`: `"1.0"`
- `extractor`: `"nmr_extract"`
- `extractor_version`: `"0.2.0"`
- `generated_at_utc`: `"2026-05-31T13:39:54Z"`
- `protein_id`: `"extract-full-mopac-allframes-20260531"`
- `axis_sizes.atom`: `846`
- `axis_sizes.residue`: `54`
- `axis_sizes.bond`: `862`
- `axis_sizes.ring`: `16`
- `axis_sizes.ring_membership`: `96`

The `npys/extraction_manifest.json` is a duplicate-style extractor manifest
with the same schema/axis facts, but its `protein_id` is `"npys"`.

The H5 root attribute `protein_id` is `"md"`. For downstream identity, use
`.LGS protein_id`; the existing `.LGS` stamps `1P9J_5801`.

Frame count and stride:

- `command.sh` invokes `nmr_extract --stride 1 --mopac`.
- The command comment says this is "FULL stride-1 (every TRR frame, all 1501)".
- H5 `/trajectory/frames/original_index`: 1501 entries, `0..1500`, unique diff
  `[1]`.
- H5 `/trajectory/frames/time_ps`: 1501 entries, `0.0..15000.0`, unique diff
  `[10.0]` ps.
- All H5 time-series `frame_indices` datasets checked have 1501 entries,
  `0..1500`, unique diff `[1]`.
- All H5 time-series `frame_times` datasets checked have 1501 entries,
  `0.0..15000.0`, unique diff `[10.0]` ps.
- `pdbs/` contains 1501 PDB files:
  - first: `md_f000000_t0.0.pdb`
  - last: `md_f001500_t15000.0.pdb`
- `npys/` contains 1501 frame directories:
  - first: `frame_000000`
  - last: `frame_001500`
- Every frame directory checked by distribution has 91 `.npy` files:
  `{91: 1501}`.

Per-frame NPY content:

- `frame_000000` contains 91 arrays.
- MOPAC-specific arrays include:
  - `mopac_bond_orders.npy`
  - `mopac_charges.npy`
  - `mopac_coulomb_E.npy`
  - `mopac_coulomb_efg_aromatic.npy`
  - `mopac_coulomb_efg_backbone.npy`
  - `mopac_coulomb_scalars.npy`
  - `mopac_coulomb_shielding.npy`
  - `mopac_global.npy`
  - `mopac_mc_category_T2.npy`
  - `mopac_mc_scalars.npy`
  - `mopac_mc_shielding.npy`
  - `mopac_scalars.npy`
- Example shapes from `frame_000000`:
  - `pos.npy`: `(846, 3)`, `float64`
  - `element.npy`: `(846,)`, `int32`
  - `mopac_charges.npy`: `(846,)`, `float64`
  - `mopac_bond_orders.npy`: `(896, 3)`, `float64`
  - `mopac_global.npy`: `(4,)`, `float64`
  - `mopac_scalars.npy`: `(846, 4)`, `float64`
  - `mopac_coulomb_shielding.npy`: `(846, 9)`, `float64`
  - `mopac_mc_shielding.npy`: `(846, 9)`, `float64`

Frame indexing:

- The authoritative dense extraction row map is H5
  `/trajectory/frames/original_index`.
- For this all-frame extraction, H5 row `r` maps to original TRR frame `r`.
- `npys/frame_NNNNNN` and `pdbs/md_fNNNNNN_t<ps>.pdb` use the same original
  frame number.
- No separate JSON frame-index manifest mapping NPY rows to trajectory frames
  was found. The mapping is in H5 and in per-frame directory/PDB names.

IUPAC note: topology sidecars include IUPAC/BMRB naming fields, but these are
not involved in frame-set alignment.

## 3. ORCA DFT Side

Calcset parent:

`/shared/2026Thesis/shielding-calcsets/data/trajectories/1p9j-calibration-with-dft`

Parent contents include:

- `1p9j-calibration-with-dft.LGS`
- `README.md`
- `dft/`
- `extract/`
- `extract-full-mopac-allframes-20260531/`
- `md/`
- `topol.top`

Existing `.LGS`:

`/shared/2026Thesis/shielding-calcsets/data/trajectories/1p9j-calibration-with-dft/1p9j-calibration-with-dft.LGS`

Existing `.LGS` identity:

- `schema_version`: `1`
- `kind`: `"trajectory"`
- `dataset_id`: `"1p9j-calibration-with-dft"`
- `protein_id`: `"1P9J_5801"`
- `human_name`: `"1P9J calibration trajectory (15 ns r2SCAN/def2-SVP DFTs)"`

Existing `.LGS trajectory` block:

- `md_dir`: `"md"`
- `topology_top`: `"topol.top"`
- `extraction_dir`: `"extract"`
- `trajectory_h5`: `"extract/trajectory.h5"`
- `extraction_manifest`: `"extract/extraction_manifest.json"`
- `frame_dt_ps`: `10.0`
- `frame_index_basis`: `"trr_frame_index"`

The existing `.LGS` therefore points at the sparse `extract/` directory, not
the dense all-frame MOPAC extraction.

Existing `.LGS dft` block:

- `method`: `"r2SCAN def2-SVP def2/J NMR  CPCM(Water)"`
- `campaign_target_frames`: `751`
- `frame_stride`: `{ "first": 0, "last": 1002, "step": 2 }`
- `frames[]` count: `500`

Important: `campaign_target_frames` is 751, but the current consolidated DFT
jobs and current `.LGS frames[]` contain 500 completed jobs.

Current ORCA job directory facts:

- `dft/jobs/` has 500 job directories.
- There are 500 `*_meta.json` files under `dft/jobs/*/`.
- All 500 `.LGS dft.frames[].frame_index` values match their meta JSON
  `frame_index`.
- All 500 meta JSON `frame_ps` values equal `frame_index * 10.0`.
- All 500 meta JSON `orca_exit_code` values are `0`.
- All 500 meta JSON `protein_id` values are `1P9J_5801`.
- All 500 meta JSON `atom_count` values are `846`.
- All 500 declared `files.out_primary` paths exist.

Campaign/status facts:

- `dft/_consolidation_snapshot.json` says snapshot date
  `2026-05-26T22:46:10+01:00`, `jobs: 500`, `campaign_target: 751`.
- `dft/status/jobs.csv` has 752 lines: one header plus 751 campaign rows.
- Status counts from `jobs.csv`:
  - `done`: 500
  - `running`: 4
  - `pending`: 247
- The CSV job list is the full planned even-frame set `0,2,...,1500`.
- Running frames at this snapshot: `910`, `932`, `1004`, `1006`.
- Pending frames at this snapshot: `1008,1010,...,1500`.

The current `.LGS dft.frames[]` present/completed frame set is:

- Even original TRR frame indices from `0` through `1002`,
- excluding `910` and `932`,
- no odd frames,
- no frames after `1002`.

Diffs between adjacent present DFT frame indices:

- `2` for 497 adjacent pairs.
- `4` for 2 adjacent pairs: across missing `910` and missing `932`.

So the current present DFT frame set is sparse and slightly irregular. It is
not exactly `range(0, 1004, 2)` because `910` and `932` are absent, and it is
not the full planned `range(0, 1501, 2)` because frames after `1002` are not
consolidated into `.LGS` yet.

## 4. Frame Correspondence

The all-frame MOPAC extraction and ORCA DFTs are on the same 1P9J trajectory
and use the same frame-0 origin.

Verified correspondence:

- Dense MOPAC all-frame H5 `/trajectory/frames/original_index` is
  `0,1,2,...,1500`.
- Dense MOPAC all-frame H5 `/trajectory/frames/time_ps` is
  `0.0,10.0,20.0,...,15000.0`.
- Existing sparse `extract/` H5 `/trajectory/frames/original_index` is
  `0,2,4,...,1500`.
- Existing sparse `extract/` H5 `/trajectory/frames/time_ps` is
  `0.0,20.0,40.0,...,15000.0`.
- All 751 sparse `extract/pdbs/md_fNNNNNN_t*.pdb` files are byte-identical to
  the corresponding even-frame dense MOPAC all-frame PDBs.
- All 500 present ORCA input PDBs have identical atom order, residue identity,
  element identity, and coordinates to the dense MOPAC all-frame PDB for the
  same original frame. Maximum coordinate difference across all 500 checked
  ORCA-vs-dense pairs was `0.0`.

Therefore there is no detected off-by-one, no detected time offset, and no
detected trajectory-origin mismatch.

Precise mapping:

- Dense MOPAC extraction row `r` maps to original TRR frame `r`.
- For any present ORCA DFT frame with original frame `f`, the corresponding
  dense MOPAC extraction row is `r = f`.
- Since ORCA frames are even original frames, the nominal full-campaign DFT
  index `k` maps to original frame `f = 2*k` and dense MOPAC row `r = 2*k`.
- `time_ps = 10.0 * f = 20.0 * k`.

For the current present `.LGS dft.frames[]` row index `j`, do not use
`f = 2*j` after holes. The completed-list row index is shifted by missing
frames:

- For `.LGS dft.frames[]` rows `j = 0..454`: `f = 2*j`, covering `0..908`.
- Row `455`: `f = 912`; frame `910` is absent/running.
- For rows `j = 455..464`: `f = 2*(j + 1)`, covering `912..930`.
- Row `465`: `f = 934`; frame `932` is absent/running.
- For rows `j = 465..499`: `f = 2*(j + 2)`, covering `934..1002`.

Robust mapping rule for present data:

Use `.LGS dft.frames[].frame_index`, not present-list row number and not
`dft.frame_stride`, then select dense MOPAC H5/NPY row whose original index is
the same value.

## 5. Stride-Support Gap And Factual Options

Current v1 `.LGS` can point at:

- one extraction directory and one extraction H5,
- one extractor manifest,
- an optional DFT block with explicit `frame_index -> meta_json` entries.

Current h5-reader/rediscover mapping can handle this concrete dense/sparse
case if the future `.LGS` points `trajectory.extraction_dir`,
`trajectory.trajectory_h5`, and `trajectory.extraction_manifest` at the dense
all-frame MOPAC extraction while preserving the existing `dft.frames[]` entries:

- H5 rows would be `0..1500`.
- DFT loaded set would be the 500 present original indices.
- `FrameMap::dftRows()` would select dense rows
  `0,2,4,...,908,912,...,930,934,...,1002`.

What v1 `.LGS` lacks as explicit schema:

- no `extraction_stride`;
- no explicit extraction frame list;
- no explicit DFT-to-extraction-row list;
- no explicit offset field;
- no statement that `dft.frames[].frame_index` is to be joined against
  `trajectory.h5 /trajectory/frames/original_index`;
- no way for a `.LGS`-only consumer to know dense extraction coverage without
  opening H5 or scanning `npys/frame_*`;
- no robust schema field for irregular sparse DFT target rows beyond the
  existing original-frame `dft.frames[].frame_index` values.

Factual representation options:

1. Single regular selector / stride factor.

   Example shape: DFT target rows are every `N`th extraction frame, with
   optional `first` or `offset`.

   This is compact and sufficient only when the relationship is regular and
   same-origin. It would represent the planned full campaign as `first=0,
   step=2` over dense extraction rows. It is not sufficient for the current
   present 500-job set unless it also allows gaps, because frames `910` and
   `932` are absent and the campaign stops at `1002`.

2. Explicit DFT-frame-index list into the extraction frame set.

   Example shape: each `dft.frames[]` entry carries either:

   - `extraction_row`, or
   - `extraction_original_index`, or
   - an explicit target-frame-index array for the DFT set.

   This is robust to offset, holes, irregular sampling, partial campaigns, and
   future nonuniform sampling. For this calcset the explicit list would contain
   the same values as `dft.frames[].frame_index`, because dense extraction row
   equals original frame index, but that equality is a verified fact of this
   artifact rather than a schema-level guarantee.

3. Formalize the current implicit join.

   Example shape: schema states that `dft.frames[].frame_index` is in the basis
   named by `trajectory.frame_index_basis` and must join to
   `/trajectory/frames/original_index`.

   This preserves the existing compact `.LGS` shape while making the current
   loader assumption explicit. It still requires consumers to open H5 to know
   extraction row coverage, and it still does not provide a `.LGS`-only
   extraction frame list.

No decision is made here. For this specific data, the exact safe join key is
the original TRR frame index in `trr_frame_index` basis.

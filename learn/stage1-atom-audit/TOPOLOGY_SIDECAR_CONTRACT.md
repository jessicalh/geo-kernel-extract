# Topology Sidecar Contract

Date: 2026-05-10

## 2026-05-13 Evening Addendum: Contract Landed

Commits `f2781da` + `dc50917` ship the contract. Punch-list status
from the PM addendum below:

1. Diff `atoms_category_info` schema vs codex Atom Table — **DONE**.
   6 fields added (chain_id, residue_number, insertion_code,
   parent_atom_index, ff_atom_type_string, equivalence_class).
2. Extend `ArraySpec` with native_axis / irreps / units /
   sign_convention / tensor_rank / parity — **DONE**. Populated for
   all ~108 entries. Resolves `TENTATIVE_OUTSTANDING_ISSUES.md` OI-016.
3. Emit `extraction_manifest.json` — **DONE**. Includes schema_version,
   protein_id, topology population flags, axis_sizes (atom / residue /
   bond / aromatic_ring / saturated_ring / ring / ring_membership),
   axis_alignment statements per axis.
4. Emit Bond + Ring + RingMembership + Residue table sidecars —
   **DONE**. `residues.npy` (new — codex Residue Table), `bonds.npy`,
   `rings.npy`, `ring_membership.npy`. All structured-NPY; schema
   mirrored byte-for-byte in `python/tests/_topology_fixture.py`.
5. OF3 retention — **STILL DEFERRED**. 685-fleet stopped per CLAUDE.md.
6. R-side regex-mechanism refactor — **STILL DEFERRED**. OI-120 in
   `TENTATIVE_OUTSTANDING_ISSUES.md`. Catalog now carries `mechanism`
   and `native_axis` so the R refactor reads typed columns.

Beyond the PM-addendum punch list:
- Python SDK enforces invariants in `load()`: required-file check,
  manifest-vs-actual axis sizes, bond endpoint range, ring membership
  refs, residue atom_count sum. Malformed exports fail loud.
- `learn/extract.py` STAGE1_AUDIT_OUTPUTS extended with the 5 new
  files (4 NPY + 1 JSON).
- `learn/src/secondary/loader.py` no longer silently skips on
  `FileNotFoundError` / `ValueError` — logs the reason to stderr so
  fleet runs are auditable without disappearing rows.

The body of this document remains valid as the requirements record.
The detailed AM/PM addenda below stay for archaeology.

---


Scope: next Stage 1 extraction export needed by `/learn` analysis.  This is a
contract for exposing facts already owned by the C++ topology/semantic model;
it is not a request to redesign the C++ topology layer.

## 2026-05-13 PM Addendum: Refinement Against Current Source

Read after the AM addendum below. Notes from a session of doc-cleanup +
audit (commit `d136fcb`) cross-checking the contract against current
`master`:

### Most "Atom Table" fields are already typed substrate

`AtomSemanticTable` already carries: `element`, `locant`, `branch`,
`di_index`, `backbone_role`, `prochiral`, `planar_group`,
`planar_stereo`, `pseudoatom`, `polar_h`, `ring_position`, `aromatic`,
`formal_charge`, `is_exchangeable`, `equivalence_class`. Plus
`Atom::parent_atom_index`, `Residue::terminal_state` /
`protonation_variant_index`, `LegacyAmberTopology::AtomtypeString()`.

That covers ~13 of the ~17 fields the "Atom Table Minimum" lists.
The sidecar work is mostly a PROJECTION of typed data already
present, not a redesign or new typed fields.

### Build on `atoms_category_info.npy`, don't parallel-export

`CategoryInfoProjection` already emits one structured NPY per
protein with ~31 typed fields, including the AMBER → IUPAC / BMRB
name projection (see CLAUDE.md landing list + memory entry
`project_extraction_sdk`). The right move is to **extend this file**
with any missing sidecar fields rather than create a second
atom-table-shaped export. Avoids two parallel atom tables that
would themselves need a "which is canonical" pointer.

Pre-flight task next session: diff codex's Atom Table requirements
against the actual `atoms_category_info` schema, write down the
delta. The delta is likely under 5 fields.

### Bond / Ring / RingMembership: also projection

`Bond` carries `atom_index_{a,b}`, `order`, `category`,
`is_rotatable`. `Ring` carries `atom_indices`, `type_index`,
`parent_residue_index`, `fused_partner_index`, plus aromatic /
saturated separation via `RingTopology`. `ring_geometry.npy`
exists (aromatic-only — codex's invariant requires declaring
that). Each of these can be projected without new typed fields.

### The genuinely new work

1. **Array axis registry.** Declare `native_axis` (atom / residue /
   aromatic_ring / saturated_ring / ring_contribution_pair /
   protein) per `.npy`. The catalog (`python/nmr_extract/_catalog.py`)
   has `group` but not `native_axis`. Small extension to `ArraySpec`.
2. **Feature metadata.** Irreps + units + sign convention + parity +
   tensor_rank + mechanism + calculator. This is identical in shape
   to `TENTATIVE_OUTSTANDING_ISSUES.md` **OI-016**. Same work; merge
   the targets.
3. **Mechanism / calculator / family tags.** Existing `group` field
   (`biot_savart`, `haigh_mallion`, `larsen_hbond`, etc.) already
   names the calculator. "Mechanism" is a thesis-narrative grouping
   that maps deterministically from group; declare once, don't
   regex.

### Out of scope until 685-fleet path settles

The codex doc has an "OF3 Retention Contract" requiring per-residue
`si_trunk[N_res, 384]` and similar. Per CLAUDE.md Stage 2: the 685
fleet was stopped 2026-04-30 (bad chains), OF3 is the recovery
path but no OF3-generated structures are flowing yet. Spec'ing OF3
retention before structures arrive is premature — defer the
of3_token_index requirement and the OF3 sidecar layout until OF3
output exists to test against.

### Phase / scope / "first sidecar" language

"First sidecar" / "first pass" / "first model-shard builder" is OK
as one-time scoping. If it becomes "first / second / third sidecars"
each adding fields, that's the proxy-language pattern. Aim for
one-step where the typed substrate makes it cheap; only declare
additional phases when blocked on a real dependency (e.g., the OF3
deferral above is a real dependency, not a phasing tax).

### Next-session-pickup punch list

1. Diff `atoms_category_info.npy` schema vs codex Atom Table → write
   delta (~< 5 fields).
2. Extend `ArraySpec` with `native_axis`, `irreps`, `units`,
   `sign_convention`, `tensor_rank`, `parity`. Resolves OI-016
   simultaneously.
3. Emit a small manifest JSON sibling to NPYs declaring axes + row
   counts + schema version. Codex's "Manifest" requirement.
4. Emit Bond table + Ring table + RingMembership table as additive
   NPY/JSON sidecars; project from existing typed structs.
5. Defer OF3 retention to a later session.
6. Cross-check the `/learn` regex-mechanism offenders list (in
   `PYTHON_TOPOLOGY_MODEL_REQUIREMENTS.md`) once feature metadata
   ships; the R refactor is a downstream consumer of (2).

The whole package is additive against `master` and shouldn't need
any of the listed C++ types to grow new fields.

---

## 2026-05-13 AM Addendum: Minimum Bar For The Calculator/OF3 Dataset

Tomorrow's topology export should be additive.  The current calculator payloads
do **not** need to be redesigned or re-emitted in a new format.  Add a topology
sidecar block that makes the existing arrays interpretable, then let Python
perform the heavier reshaping into model-ready shards.

The sidecar's first job is to answer one question without ambiguity:

```text
raw calculator row i == topology atom j
```

If current calculator arrays are already written in conformation atom order,
the sidecar may simply declare that atom axis and expose the topology facts for
each row in that same order.  Do not introduce a second atom ordering unless
there is a documented `raw_row -> atom_index` mapping table.

### String Projection Boundary

The naming barrier is intentional and must stay intact.  This sidecar is allowed
to export strings because Python, RefDB/BMRB joins, audit tables, and human
inspection need a string projection surface.  That does **not** mean the C++
topology or calculators should start using strings as chemistry.

Acceptable:

- render already-owned typed facts into names at the export boundary;
- export loaded/source atom names as provenance;
- export force-field atom type strings when `LegacyAmberTopology` already has
  `atomtype_string`;
- export BMRB/IUPAC/Amber-facing labels when an existing naming surface can
  provide them;
- let Python use those labels for reporting and external label joins, with
  coverage diagnostics.

Not acceptable:

- parsing atom-name strings inside calculators to decide chemistry;
- replacing `AtomSemanticTable`, `AtomMechanicalIdentity`, bonds, residues, or
  ring topology with string-keyed logic;
- synthesizing authoritative chemistry from display names in the exporter;
- making Python infer missing topology from regexes without marking that as a
  quarantined audit fallback.

In short: strings are an exported projection and provenance surface, not the
internal substrate.  If a desired label cannot be rendered from facts the
system already owns, export explicit unknown/null and let Python report the
gap.

### Current C++ Surface To Export From

Keep the first pass grounded in the fields that exist now:

- `Atom`: `element`, `pdb_atom_name`, `residue_index`, `bond_indices`,
  `parent_atom_index`;
- `Residue`: `type`, `sequence_number`, `chain_id`, `insertion_code`,
  `atom_indices`, `protonation_variant_index`,
  `protonation_state_resolved`, `terminal_state`, cached backbone indices
  (`N`, `CA`, `C`, `O`, `H`, `HA`, `CB`), and chi atom slots;
- `LegacyAmberTopology`: atom/residue counts, `BondList()`,
  `BondIndicesFor(atom)`, `HydrogenParentOf(atom)`, ring topology,
  optional invariant FF fields (`Mass`, `FfAtomTypeIndex`, `Ptype`,
  `AtomtypeString`, `Exclusions`, `FudgeQq`, `RepPow`, `Atnr`), and
  `AtomSemantic()` when populated;
- `Bond`: `atom_index_a`, `atom_index_b`, `order`, `category`,
  `is_rotatable`;
- `Ring`: `atom_indices`, `type_index`, `parent_residue_index`,
  `parent_residue_number`, `fused_partner_index`, plus aromatic/saturated
  axis separation from `RingTopology`;
- `AtomSemanticTable`: `element`, `locant`, `branch`, `di_index`,
  `backbone_role`, `prochiral`, `planar_group`, `planar_stereo`,
  `pseudoatom`, `polar_h`, `ring_position`, `aromatic`, `formal_charge`,
  `is_exchangeable`, `equivalence_class`.

Do not ask tomorrow's export to invent fields not carried by these objects.
Runtime `AtomSemanticTable` does not carry per-atom `SemanticProvenance`; that
audit trail lives in the generated semantic-table log.  The sidecar may point
to that log in the manifest, but per-atom provenance should not block the first
export.

### Non-Negotiable Minimum For The First Sidecar

The first acceptable implementation exports enough information for Python to
build a strict atom-level training artifact from:

```text
existing calculator export
+ topology sidecar
+ OF3 residue embeddings
+ RefDB/BMRB labels
```

Minimum required content:

- manifest with `schema_version`, `protein_id`, source structure identifier,
  protonation/export policy metadata, extractor version/config when available,
  and declared native axes for every array touched by the sidecar;
- atom table with one row per exported atom-axis row;
- residue table with one row per exported residue-axis row;
- bond table with atom-index endpoints and typed bond category;
- enough residue identity to map OF3 residue tokens in Python, even if the C++
  export does not know OF3 token indices directly;
- explicit null/unknown sentinels for facts the topology layer cannot yet
  provide;
- validation counts: atom rows, residue rows, bond rows, and calculator atom
  axis row count.

Do not block the first export on perfect ontology.  Block it on row-order
ambiguity, missing atom-to-residue mapping, missing bond endpoints, or hidden
string-only chemistry used as internal authority.

### Atom Table Minimum

One row per calculator atom-axis row, in the same order as atom-level
calculator arrays unless `raw_row_index` says otherwise.

Required fields:

- `atom_index`: stable integer key on this export's atom axis;
- `raw_row_index`: row in the existing calculator atom arrays, if different
  from `atom_index`;
- `residue_index`: stable integer key into the residue table;
- `chain_id`: explicit value or typed null;
- `residue_number`: biological/source residue number;
- `insertion_code`: explicit value or typed null;
- `residue_type`: three-letter or enum-coded residue identity;
- `atom_name_loaded`: source/display/provenance string, usually
  `Atom::pdb_atom_name`;
- `atom_name_canonical`: rendered label if already available; otherwise typed
  null, not a newly invented parser result;
- `element`;
- `locant`;
- `branch_outer`;
- `branch_inner`;
- `diastereotopic_index`;
- `backbone_role`;
- `prochiral`;
- `planar_group`;
- `pseudoatom_kind` or equivalent pseudoatom membership fields if already
  available from `AtomSemanticTable`;
- `ring_position` / `aromatic` if `AtomSemanticTable` is populated;
- `parent_atom_index`: for hydrogens and grouped atoms, else null;
- `terminal_flag`;
- `protonation_variant` or a pointer to protein-level protonation policy;
- `ff_atom_type_string` if `LegacyAmberTopology::AtomtypeString()` is
  populated; otherwise typed null.

The semantic fields above are required when `LegacyAmberTopology::HasAtomSemantic()`
is true.  If the load path legitimately lacks `AtomSemanticTable`, the sidecar
must declare that fact in the manifest and write typed unknowns rather than
falling back to string inference.

Nice-to-have, but not required for the first sidecar if unavailable:

- symmetry class;
- ring membership role;
- BMRB/IUPAC-facing atom names.

Those fields matter, but Python can tolerate explicit unknowns during the first
full 2300-protein production pass.  Python should not infer them from display
strings unless the inference is quarantined in an audit/report step.

### Residue Table Minimum

Required fields:

- `residue_index`;
- `chain_id`;
- `residue_number`;
- `insertion_code`;
- `residue_type`;
- `sequence_index`: zero-based sequence position in the exported protein chain
  when known;
- `prev_residue_index`;
- `next_residue_index`;
- `terminal_state`;
- `protonation_variant` when residue-specific;
- `of3_token_index` if known, otherwise enough source identity for Python to
  join to OF3 tokens deterministically.

OF3 is residue/token-level.  The model shard builder will attach
`si_trunk[of3_token_index]` to atoms through `atom.residue_index`.  Do not
duplicate OF3 vectors per atom in the C++ export.

### Bond Table Minimum

Required fields:

- `bond_index`;
- `atom_index_a`;
- `atom_index_b`;
- `bond_category`;
- `bond_order` if known, else explicit unknown;
- `is_aromatic` if known;
- `is_peptide` if known.

Optional:

- `source_provenance` if the exporter can state whether a bond came from
  geometric covalent resolution, aromatic overlay, or disulfide authority.

The first model-shard builder can derive graph edges directly from this table.
Any row whose endpoint is absent from the atom table is a hard export failure.

### Calculator Payload Handling

The sidecar must not require a rewrite of existing calculator outputs.  Python
will normalize them into a model shard with a shape like:

```text
calculator/
  names             [N_channels]
  T0                [N_atom, N_channels]
  T1                [N_atom, N_channels, 3]
  T2                [N_atom, N_channels, 5]
  mat3              [N_atom, N_channels, 3, 3]
  mask              [N_atom, N_channels]
  provenance        per channel
```

The export only needs to make the native axis and row mapping explicit enough
for Python to reorder/check arrays safely.  Silent zero fill is forbidden:
missing or not-applicable calculator values must become masks in the model
shard.

### OF3 Retention Contract

OF3 embeddings remain residue-level sidecars, not calculator outputs:

```text
of3/
  si_trunk          [N_res, 384]
  zij_trunk         [N_res, N_res, 128]  optional/lazy sidecar
  zij_row_mean      [N_res, 128]         optional cheap summary
```

The topology sidecar must carry enough residue identity for Python to map
`residue_index -> of3_token_index`.  If a residue cannot be mapped to an OF3
token, it remains valid only with an explicit OF3 mask.

### Python Model-Shard Builder Is The Heavy Layer

Expected downstream flow:

```text
raw calculator export
topology sidecar
OF3 embeddings
RefDB/BMRB labels
        |
        v
Python model-shard builder
        |
        v
atom-level graph/tensor artifact for MACE, gated tensor mixers, and baselines
```

The model shard builder may perform extensive massaging:

- canonicalize ids;
- reorder calculator arrays to topology atom order;
- attach residue and OF3 token indices;
- construct bond/spatial/residue-pair edges;
- convert per-calculator tensor fields into channel tensors;
- attach RefDB labels and masks;
- emit rejected joins and coverage reports.

The C++ sidecar should therefore prioritize faithful exported facts and
deterministic ids over premature ML-specific layout.

### First-Pass Validation Gates

Python should reject the export when any of these fail:

- atom table row count equals the declared atom axis count;
- every atom has a residue that exists in the residue table;
- every bond endpoint exists in the atom table;
- atom-level calculator arrays either share the atom axis row count or declare
  a valid mapping/projection;
- no calculator array is consumed without a declared native axis;
- unknown semantic values are explicit, not missing;
- OF3 join coverage is reported per protein;
- label join coverage is reported per atom role/element;
- tensor fields preserve both `mat3` and irreps when both are available.

This is the set bar.  More topology is welcome.  Less than this leaves the
2300-protein calculator run too hard to audit.

## Purpose

The analysis layer must stop reconstructing chemistry from names, regexes, and
array positions.  The next extraction cycle should export a topology sidecar
that lets Python/R consume the C++ semantic model directly and fail closed when
required facts are absent.

The numerical calculator payloads may remain in their existing NPY form.  The
sidecar is the semantic contract that says what each topology axis means and
how axis-level values may be projected.

## Core Rule

```text
C++ owns topology and atom semantics.
Python/R consume exported topology facts.
Python/R do not infer chemistry from strings.
```

Names are allowed as display fields and provenance fields.  They are not
allowed to be the only key for topology grouping, mechanism grouping, or axis
projection.

## Required Tables

Storage may be JSON, NPY plus JSON schema, Arrow/Parquet, or another explicit
format.  The table names below describe the contract, not mandatory filenames.

### Manifest

Required fields:

- `schema_version`
- `protein_id`
- `extractor_version`
- `extractor_config_hash` if available
- `source_structure`
- `topology_source`
- `array_axis_registry`
- `enum_vocab_refs`
- `generated_at`

The manifest should make it impossible to confuse two exports with different
axis definitions.

### Atom Table

One row per atom on the exported atom axis.

Required fields:

- `atom_index`
- `residue_index`
- `element`
- `loaded_atom_name`
- `canonical_atom_name`
- `amber_atom`
- `amber_residue`
- `iupac_atom`
- `bmrb_atom`
- `bmrb_residue`
- `atom_position`
- `nmr_class`
- `backbone_role`
- `sidechain_role`
- `parent_atom_index`
- `pseudo_atom_name`
- `symmetry_class_id`
- `prochirality`
- `terminal_flag`
- `naming_provenance`

Fields that are unavailable should be explicitly null/unknown through a typed
sentinel, not omitted or approximated by downstream code.

### Residue Table

One row per residue on the exported residue axis.

Required fields:

- `residue_index`
- `chain_id` if available, otherwise explicit null;
- `residue_number`
- `residue_type`
- `prev_residue_index`
- `next_residue_index`
- `prev_residue_type`
- `next_residue_type`
- `terminal_state`
- `protonation_variant` if known;
- `is_proline`
- `is_xpro_context`

### Ring Table

One row per topology ring, including both aromatic and saturated rings.

Required fields:

- `ring_id`
- `ring_kind`: `aromatic`, `saturated`, or typed other value;
- `native_axis`: usually `aromatic_ring` or `saturated_ring`;
- `native_axis_index`
- `ring_type`
- `source_residue_index`
- `source_residue_number`
- `fused_partner_ring_id` if present;
- `atom_count`
- optional existing calculator axis ids, when different from `ring_id`.

Conformation-dependent geometry may live here or in a separate geometry table.
If stored separately, it must still reference `ring_id`.

### Ring Membership Table

One row per ring-member atom.

Required fields:

- `ring_id`
- `atom_index`
- `ring_atom_order`
- `role_in_ring` if known;
- `is_vertex`
- `is_substituent` if known.

This is the only acceptable basis for saturated-ring pucker projection to
atoms.  If this table is absent, pucker stays on the saturated-ring axis.

### Ring Geometry Table

One row per ring geometry value, if ring geometry is exported.

Required fields:

- `ring_id`
- `geometry_axis`
- `center_x`
- `center_y`
- `center_z`
- `normal_x`
- `normal_y`
- `normal_z`
- `radius`
- `geometry_provenance`

If the existing `ring_geometry.npy` is aromatic-only, the sidecar must declare
that its axis is `aromatic_ring`.  Saturated-ring topology must not be implied
from an aromatic-only geometry table.

### Array Axis Registry

Every exported calculator or derived array must declare its native axis.

Required fields:

- `array_name`
- `source_file`
- `native_axis`
- `row_count`
- `column_count`
- `axis_table`
- `axis_id_column`
- `projection_allowed`
- `projection_method` if already projected;
- `normalization_state`
- `units` if meaningful.

Examples:

```text
atom_positions.npy          native_axis=atom
aromatic_chi2.npy           native_axis=aromatic_ring
pucker_Q.npy                native_axis=saturated_ring
ring_contributions.npy      native_axis=ring_contribution_pair
protein_scalar_features.npy native_axis=protein
```

### Feature Metadata Table

Every analysis feature produced from arrays should carry explicit metadata.

Required fields:

- `feature_id`
- `display_name`
- `source_array`
- `source_columns`
- `calculator`
- `feature_family`
- `mechanism`
- `physical_quantity`
- `tensor_rank`
- `parity`
- `native_axis`
- `target_axis_after_projection`
- `projection_method`
- `normalization_state`
- `units`
- `is_topology_inferred`
- `provenance_notes`

R should use this table for mechanism and family grouping.  Regex grouping over
feature names is discovery-code only.

## Projection Rules

Allowed projections must be explicit and table-backed:

```text
residue -> atom       through AtomTable.residue_index
ring -> member atoms  through RingMembershipTable
ring -> nearby atoms  derived geometry only, not membership
protein -> atom       broadcast, marked protein/global
```

No analysis code should silently broadcast, expand, or aggregate across axes.
Every projection should leave a metadata trail.

## Validation Invariants

Python should reject an export when any required invariant fails:

- atom table row count equals the atom axis used by atom-level arrays;
- every atom `residue_index` exists in the residue table;
- every ring-member `atom_index` exists in the atom table;
- every ring-member `ring_id` exists in the ring table;
- `native_axis_index` is unique within each ring native axis;
- aromatic ring arrays have row count equal to the aromatic-ring axis count;
- saturated-ring arrays have row count equal to the saturated-ring axis count;
- ring contribution indices reference declared ring axes;
- feature metadata exists for thesis-facing features;
- unknown or fallback semantic values remain visible in provenance columns.

An empty axis is valid only when explicitly declared.  For example, a protein
with no saturated rings may have zero pucker rows.

## Expected Payoff

After this sidecar exists, `/learn` can be rewritten as a strict semantic
consumer:

- no protein-chemistry regexes in R;
- no pucker/ring projection guesses in Python;
- no display-string category keys as primary grouping keys;
- no hidden calculator-axis assumptions;
- clean comparison against the frozen 720 discovery baseline.

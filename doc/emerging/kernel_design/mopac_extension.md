# MOPAC extension: hold the completed run whole

Objective restated: make the full understanding of a completed PM7+MOZYME
MOPAC run available downstream and emitted for reader use, without a new
architecture. The shape is still one `MopacResult`, its NPYs, and a sidecar.
Do not curate the MOPAC facts into only the features we already know how to use.

## What the current code does

`src/MopacResult.cpp` writes a `.mop` file with:

```
PM7 MOZYME 1SCF CHARGE=<net_charge> BONDS MULLIK LET GEO-OK THREADS=<threads>
```

The input body is one Cartesian atom line per conformation atom, in the same
atom order as the protein, with optimization flags set to `0`; `1SCF` makes the
run a single-point calculation, not a geometry optimization. The calculator runs
the MOPAC binary as a subprocess because the linked in-process route leaves the
Fortran `.out` incomplete when read. MOPAC writes `.out` alongside the `.mop`;
on success the current code removes `.mop`, `.out`, `.arc`, and `.den`.

The current parser reads the `.out` as text, requires normal termination, and
extracts only:

- `FINAL HEAT OF FORMATION = ...` into `heat_of_formation_`.
- `NET ATOMIC CHARGES AND DIPOLE CONTRIBUTIONS` rows into per-atom Mulliken
  charge, s population, and p population.
- `BOND ORDERS` rows matching the current regular expression, filtered to
  `bo > 0.01`.
- the `DIPOLE X Y Z TOTAL` table's `SUM` row, but only the x/y/z vector.

The current in-memory `MopacResult` exposes atom charge/s/p/derived valency,
pair and topology-bond Wiberg lookup, the unique-pair bond-order map, heat of
formation, and dipole. It stores charge/s/p/derived valency on
`ConformationAtom`, and it builds per-atom MOPAC bond-neighbour lists.

The current emit is exactly four arrays:

- `mopac_charges.npy` `(N,)`: per-atom Mulliken charge.
- `mopac_scalars.npy` `(N,4)`: charge, s population, p population, derived
  valency.
- `mopac_bond_orders.npy` `(B,3)`: unordered unique pairs, `[atom_i, atom_j,
  order]`, with `B` equal to MOPAC's unique reported-pair count, not the
  topology bond count.
- `mopac_global.npy` `(4,)`: heat of formation and dipole x/y/z.

The Python catalog mirrors those four `mopac_core` arrays, and the reader model
already treats MOPAC core as a three-axis group: atom, sparse MOPAC bond/pair,
and protein.

## Inventory: what the completed run produces

This inventory is the capture contract, not a filter. The implementation should
store every section/table MOPAC prints, even if a row type is not named here or a
future MOPAC version adds a heading.

### Run files and provenance

MOPAC's normal command-line run produces a main `.out` file and commonly an
`.arc` archive file; optional keywords can also produce `.aux`, `.den`, `.html`,
or `.pdb` style artifacts. The project already knows `.arc` and `.den` may
exist because it removes them after successful parsing.

Capture:

- raw input `.mop`, exact keyword line, title/comment lines, atom coordinate
  lines, net charge, thread count, and path/binary provenance;
- raw output artifacts exactly as written: at minimum `.out` and `.arc` when
  present, plus `.den`, `.aux`, `.html`, `.pdb`, or any other same-stem MOPAC
  product if produced;
- MOPAC version/banner/date, method, title, keyword echo, atom count, electron
  count, empirical formula, point group if printed, termination status, warnings,
  errors, SCF/MOZYME convergence records, iteration counts, CPU/timing/memory
  records, and all "normal end" diagnostics.

The raw files matter. They are the lossless source of truth when a parser misses
a MOPAC section, when a future calculator wants a field not yet indexed, or when
the reader needs to audit a suspicious value.

### Structure and labels

Capture every structural label and coordinate table MOPAC emits:

- atom ordinal, element/symbol, optional label text, and original project atom
  index mapping;
- input Cartesian coordinates and optimization flags from the `.mop`;
- any MOPAC-printed Cartesian coordinates, internal coordinates, final
  coordinates, archive coordinates, atom reordering/resequencing diagnostics,
  connectivity/topography/Lewis-structure labels, bond lengths, bond angles,
  torsions, and interatomic-distance sections;
- any PDB/archive labels MOPAC writes if PDB-like output is requested or
  generated.

For the current `1SCF` path, the final geometry should normally be the input
geometry, but the extension must not infer that as a reason to omit final
coordinate sections. Store both when both are printed.

### Energies and scalar molecular properties

Capture all graph-level numeric rows as label/value/unit records, including at
least the terms MOPAC commonly prints for final SCF results:

- final heat of formation;
- electronic energy;
- nuclear/core-core repulsion energy;
- total energy or related energy summaries when printed;
- ionization potential / frontier-orbital summary when printed;
- molecular weight, formula, area, volume, gradient norm, dielectric/solvation
  terms, or other scalar molecular properties when printed;
- the full molecular dipole table, not only the current `SUM` vector.

The dipole table is specifically full-fidelity:

- `POINT-CHG.` x/y/z/total;
- `HYBRID` x/y/z/total;
- `SUM` x/y/z/total;
- units and origin caveats. For charged systems, the dipole is origin-dependent.

### Atomic charges, electron density, and populations

The current parser captures only a small part of the MULLIK output. Capture every
row MOPAC prints under the Mulliken and charge/population sections:

- net atomic charge and electron density per atom;
- per-atom dipole-contribution fields where present;
- s, p, and any d/f or shell-resolved population columns MOPAC prints;
- AO basis labels per atom: AO ordinal, atom index, AO symmetry/type (`S`, `PX`,
  `PY`, `PZ`, etc.), zeta, principal quantum number, and any shell metadata;
- full Mulliken population analysis tables, including AO populations, overlap
  populations, and per-atom/per-AO contributions;
- if `LARGE`/related output exposes Medrano-Bochicchio-Reale population
  analysis, store its non-shared charge, active charge, total charge, valence,
  free valence, statistical promotion, and Mulliken promotion columns.

Do not collapse element-dependent missing columns into zeros without metadata.
Use `NaN` for non-applicable numeric cells in convenience NPYs and record the
reason in the sidecar column metadata; the sidecar's table records remain the
authoritative exact rows.

### Bond orders, valencies, bonding contributions

The `BONDS` keyword prints Wiberg-equivalent rotationally invariant bond orders.
The MOPAC manual states that, for MOZYME, `BONDS` prints only bonds above 0.01
and not involving hydrogen; `ALLBONDS` prints all bonds above 0.001 including
hydrogen. It also states that diagonal terms are the valencies, and that bonding
contributions of all molecular orbitals are printed before the bond matrix.

Capture:

- every printed bond-order row exactly as MOPAC prints it, including row atom,
  row element, MOPAC's parenthesized valency/diagonal value, neighbour ordinals,
  neighbour elements, bond orders, row order, and duplicate reverse rows;
- the unique-pair symmetric projection as a convenience view, but never as the
  only retained form;
- MOPAC's own per-atom valency diagonal separately from the project's derived
  valency;
- hydrogen-involving and low-order bonds when the run is configured to print
  them (`ALLBONDS`), with the print threshold recorded;
- topology-bond bridge rows, including whether each topology bond is present in
  the MOPAC print, absent because below print threshold, absent because MOZYME
  did not print it, or unmatched because MOPAC reports a pair outside topology;
- molecular-orbital / localized-orbital bonding contribution tables printed
  before the bond matrix.

The existing `mopac_bond_orders.npy` is a useful legacy projection, but it is
not the full MOPAC bond-order object: it filters, symmetrizes, loses reverse
rows and row ordinals, loses MOPAC's diagonal valencies, loses hydrogens under
`BONDS`, and stores pairs in unordered-map iteration order.

### Orbitals, matrices, and electronic-state artifacts

MOZYME uses localized molecular orbitals, and MOPAC can print or emit orbital
and matrix data via output keywords and auxiliary/density files. The extension
must preserve all such artifacts produced by the run:

- MO/LMO metadata: ordinal, energy, occupation, symmetry/label if printed,
  bonding contribution, atom/AO composition, and any local-orbital labels;
- MO/LMO coefficient matrices when printed or emitted;
- density matrix, overlap matrix, one-electron/Fock-like matrices, or other
  square/lower-triangular matrices if printed or emitted;
- `.den` contents when MOPAC writes density/LMO restart data;
- `.aux` contents when requested, because it is MOPAC's compact,
  self-descriptive machine-readable result stream.

These are large, but size is not a capture criterion. They may live in raw
sidecar artifacts plus NPY matrix blocks. If a matrix is naturally lower
triangular in MOPAC's file, store the packed original order and provide an
optional dense convenience view only as a derived projection.

## Shape of the extension

### In-memory `MopacResult`

Extend the existing `MopacResult`; do not create a new architecture. The result
becomes the held-whole MOPAC buffer:

- keep the current public query methods and current four NPYs stable;
- add a `MopacRunRecord`-like payload inside `MopacResult` containing raw
  artifacts, parsed sections, parsed tables, and typed indexes;
- expose typed accessors by natural axis:
  - atom: charge, electron density, shell/AO populations, MOPAC valency,
    project-derived valency, per-atom dipole contribution, atom labels;
  - AO: atom mapping, AO type, zeta/PQN, population, matrix row/column mapping;
  - MO/LMO: energy/occupation/label/contribution/coefficient row;
  - directed printed bond-order row: row atom plus neighbour entries exactly as
    printed;
  - symmetric atom pair: max/mean/original entries, print-threshold/provenance,
    topology bridge;
  - topology bond: current bridge plus presence state;
  - protein/global: all scalar terms and all dipole components;
  - arbitrary section/table: name, byte offsets, heading, raw lines, parsed
    cells.

The arbitrary section/table accessor is load-bearing. It prevents the typed API
from becoming an accidental editorial boundary.

### Output artifacts

Emit a sidecar directory or JSON sibling under the same per-frame/per-pose NPY
directory. Suggested shape:

```
mopac_full/
  mopac.mop
  mopac.out
  mopac.arc
  mopac.aux
  mopac.den
  mopac_full.json
```

Only files actually produced by MOPAC exist, but `mopac_full.json` records the
expected and missing artifacts explicitly. Raw text files are byte-for-byte
copies; binary files are copied byte-for-byte with size and SHA-256 in the JSON.

`mopac_full.json` is the authoritative structured sidecar:

- schema version;
- run provenance and keyword line;
- project atom mapping;
- artifact manifest with paths, byte sizes, hashes, text/binary flags;
- `sections[]`: section label, start/end byte offsets in `.out`, heading lines,
  raw line range, parser status;
- `tables[]`: table id, section id, axis, row labels, column labels, units,
  original strings, parsed numeric values, parse warnings;
- `columns[]` metadata for every convenience NPY: filename, shape, dtype, axis,
  column name, units, source table id, source column, missing-value policy,
  whether it is original or derived.

The sidecar must include original strings as well as parsed numbers. Parsed
`double` values are for computation; original strings are for fidelity and
audit.

### NPY projections

Keep existing arrays unchanged for compatibility:

- `mopac_charges.npy`;
- `mopac_scalars.npy`;
- `mopac_bond_orders.npy`;
- `mopac_global.npy`.

Add full-fidelity and convenience NPYs, each with column metadata in the sidecar:

- `mopac_global_terms.npy` `(G,)` or `(G,C)`: every graph-level numeric term,
  keyed by sidecar labels/units.
- `mopac_dipole_components.npy` `(3,4)`: rows `POINT-CHG.`, `HYBRID`, `SUM`;
  columns x/y/z/total.
- `mopac_atom_populations.npy` `(N,C_atom)`: all numeric atom-row charge,
  electron-density, shell population, dipole contribution, MOPAC valency, and
  promotion fields found in MOPAC tables.
- `mopac_ao_table.npy` structured or numeric+sidecar `(NAO,C_ao)`: AO ordinal,
  atom index, AO type id, zeta, PQN, and numeric AO population columns.
- `mopac_mulliken_overlap_sparse.npy` `(K,3+)`: sparse AO/AO or atom/atom
  overlap-population rows exactly as the table supports.
- `mopac_bond_orders_printed.npy` `(K,...)`: directed printed rows, retaining
  row ordinal, row atom, neighbour atom, elements, order, and source line.
- `mopac_bond_valencies.npy` `(N,)`: MOPAC diagonal valency, not recomputed.
- `mopac_bond_orders_unique.npy` `(U,...)`: deterministic sorted symmetric
  pairs with all original directed-row references.
- `mopac_topology_bond_orders_full.npy` `(B_topo,...)`: topology-bond bridge
  with order, present flag, source-pair id, and absence reason.
- `mopac_mo_meta.npy` `(NMO,C_mo)`: MO/LMO energy, occupation, bonding
  contribution, and labels.
- `mopac_mo_coefficients.npy` `(NMO,NAO)` when emitted by MOPAC/AUX.
- `mopac_density_packed.npy` and `mopac_overlap_packed.npy` when emitted, with
  packed-order metadata; dense forms are optional convenience projections, not
  replacements.

For variable-width/string-heavy tables, do not force bad numeric NPYs. Keep the
table in `mopac_full.json` with row/column labels and only add NPYs for numeric
blocks that benefit downstream calculators.

### Run configuration

The current keyword set prints useful but not maximal MOPAC information. For the
extension, use print keywords that expose rather than alter the PM7+MOZYME
single-point physics, with implementation-time verification that they do not
change the SCF result:

- prefer `ALLBONDS` over `BONDS` so hydrogen-involving and 0.001-to-0.01 bond
  orders are not silently unavailable;
- keep `MULLIK`;
- request MOPAC's machine-readable auxiliary output with high precision, without
  suppressing density, overlap, or MO blocks;
- preserve `.den` if MOPAC writes it, and request density/LMO output if needed
  to expose matrices that the text `.out` does not carry;
- consider `LARGE`, `LEWIS`, and `LOCAL` as disclosure keywords for population,
  Lewis/topography, and localized-orbital detail. If any such keyword changes
  the run semantics or makes MOZYME fail for proteins, the raw `.out`/`.aux`
  sidecar still captures everything the successful configured run produced, and
  the blocked keyword is recorded as an implementation caveat rather than a
  reason to drop the data class.

The sidecar is the guard against version-specific output. The parser should be
section-preserving first and table-specific second.

## Anticipated convenient forms, not capture boundaries

Known downstream features want these forms, but they do not define what is
captured:

- Charge/EFG kernels want MOPAC Mulliken charges as a parallel charge channel
  beside ff14SB, AIMNet2, APBS, and explicit-water fields. Convenient forms:
  per-atom charge, per-atom electron-density/population descriptors, and
  `1o`/`2e` field features computed by later calculators from the charge
  channel.
- McConnell/bond-anisotropy kernels want Wiberg bond orders as source weights.
  Convenient forms: deterministic symmetric pair table, topology-bond bridge,
  per-atom sorted MOPAC neighbours, bond-order present/missing flags, MOPAC
  valency, and directed raw rows for audit.
- e3nn/kernel substrate wants parallel channels, explicit axes, and irreps.
  Convenient forms: atom-axis `0e` local electronic scalars; bond/pair-axis
  source weights; AO/MO channels as optional high-dimensional local electronic
  descriptors; graph-level scalars such as energy and dipole as conditioning
  metadata, not as the only MOPAC summary.
- Reader/UI wants "what MOPAC knows" on selection. Convenient forms: selected
  atom's charge/electron density/shell populations/AO rows/MOPAC neighbours;
  selected bond's topology bond order, pair order, presence state, and raw MOPAC
  line references; selected frame's full molecular dipole decomposition and
  energy-term table.

These are anticipated views over the full sidecar, not permission to omit
anything else MOPAC prints.

## Sources read

- Project code: `src/MopacResult.{h,cpp}`, `python/nmr_extract/_catalog.py`,
  `python/nmr_extract/_tensors.py`, `h5-reader/src/model/QtMopacCoreGroup.h`,
  `h5-reader/src/model/QtResultBlocks.h`, MOPAC trajectory/Welford consumers,
  and the kernel-design notes in this directory.
- OpenMOPAC basic usage: <https://openmopac.net/guides/basics/>
- OpenMOPAC keyword list: <https://openmopac.net/Manual/allkeys.html>
- OpenMOPAC `BONDS` / `ALLBONDS`: <https://openmopac.net/Manual/bonds.html>
- OpenMOPAC `MULLIK`: <https://openmopac.net/Manual/mullik.html>
- OpenMOPAC `MOZYME`: <https://openmopac.net/Manual/mozyme.html>
- OpenMOPAC `AUX`: <https://openmopac.net/Manual/auxiliary.html>
- OpenMOPAC `OUTPUT`: <https://openmopac.net/Manual/output.html>
- Example MOPAC charge/dipole/bond-order output excerpt:
  <https://en.wikibooks.org/wiki/Computational_Chemistry/Printable_version>

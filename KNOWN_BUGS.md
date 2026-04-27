# Known bugs and architectural debt

**Scope.** Issues known to exist in the working tree at commit `130c8da`
(post-revert state, 2026-04-27). Documented honestly rather than fixed —
see "Why this document exists" below.

**Why this document exists.** A 2026-04-26 attempt to add a typed IUPAC
topology layer to `Atom` introduced architectural problems that would
have required a 10–20 session rebuild to resolve cleanly. Given finite
two-person-team capacity (the user + Claude), the rebuild's opportunity
cost (weeks of blocked thesis work, equipment costs accumulating) was
larger than its benefit. The decision was to revert the IUPAC layer and
document the bugs it surfaced — so the working tool remains productive
while the bugs are known and can be addressed when resources permit.

The architectural thinking that emerged from the attempt is preserved on
branch `architectural-thinking-2026-04-27` (origin and local). When and
if a future rebuild happens, that branch is the starting point. See
`spec/PROTEIN_TOPOLOGY_ABC_GRAPH.md` on that branch for the architecture.

---

## Bug 1 — Silent charge corruption on `--pdb` / `--protonated-pdb` / `--orca` paths

**Severity.** Physics-corrupting on affected paths. Trajectory path
unaffected.

**Where.** `src/ChargeSource.cpp:48` (and the fallback at line 56).

**What happens.** `ParamFileChargeSource::LoadCharges` builds an ff14SB
parameter-file lookup key as:

```
key = ff_resname + " " + identity.pdb_atom_name
```

The ff14SB parameter file is keyed by AMBER atom names (CA, CB, HB1,
HB2, etc.). When `pdb_atom_name` is not exactly the AMBER convention
(e.g., a CHARMM-named PDB supplied to `--protonated-pdb`, or a non-AMBER
ORCA path), the primary key misses. The fallback at line 56 tries the
canonical residue name with the same atom name; that also misses.
Control flow falls through to `ChargeSource.cpp:62-72`, which silently
assigns element-default charges (`charge = 0.0`, radius from a
hard-coded element table) with no warning, no diagnostic, no log.

**Affected load paths:**
- `--pdb` (reduce-protonated PDB → AMBER): mostly OK because reduce
  emits AMBER-convention names; can break if the input PDB had
  non-standard naming that reduce passed through.
- `--protonated-pdb` (user-attested already-protonated): broken if the
  user's input PDB used CHARMM, GROMACS-source, or other non-AMBER
  naming.
- `--orca` (tleap-built prmtop): broken if the prmtop's ATOM_NAME entries
  diverge from ff14SB-table conventions for any atom.
- `--mutant` (uses ORCA twice): same as `--orca`.

**Unaffected:**
- `--trajectory` (full-system XTC + TPR): the loader translates CHARMM →
  Standard via `NamingRegistry::TranslateAtomName` at
  `FullSystemReader.cpp:705`, so the lookup keys are correct.

**Operational guidance.** For production extractions where charge
correctness matters (calibration runs, Coulomb-derived shielding,
EFG calculations), prefer the `--trajectory` path. Or audit the
specific input proteins: the bug fires when `(ff_resname,
pdb_atom_name)` is not in `data/ff14sb_params.dat`.

**Detection.** No diagnostic surface in the current code. To detect, post
an extraction with `aimnet2_charges.npy` (or any charge-emitting NPY)
and check for atoms with `partial_charge == 0.0` whose element wouldn't
naturally produce zero.

**Architectural debt context.** This bug was identified by the audit at
`architectural-thinking-2026-04-27:spec/AUDIT_NAMING_2026-04-26/PHASE_2_CHARGES_AND_PARAMETERS.md`
during the IUPAC topology investigation. The bug is pre-existing — it
existed before the IUPAC topology landing and survives the revert.

---

## Bug 2 — String-keyed dispatch on atom names throughout `src/`

**Severity.** Architectural debt, not a runtime correctness bug in the
current code. Becomes a correctness bug in any future code that violates
the implicit "atom names are AMBER-convention" assumption.

**Where.** Approximately 50 sites across `src/` per the audit's Phase 1
inventory. Representative sites:

- `src/Protein.cpp` — `CacheResidueBackboneIndices` matches `pdb_atom_name`
  against literals "N", "CA", "C", "O", "H", "HN", "HA", "HA2", "CB".
- `src/Protein.cpp` — `DetectAromaticRings` builds a string-keyed map
  per residue and looks up ring atom names from `AminoAcidType.rings[k]`
  (currently `vector<const char*>`).
- `src/Protein.cpp` — HIS tautomer detection by string match on "HD1"
  / "HE2".
- `src/ProtonationDetectionResult.cpp` — variant detection by string
  matches on "HD1", "HE2", "HD2", "HG", "HZ1", "HZ2", "HZ3".
- `src/CovalentTopology.cpp` — bond classification keyed in places by
  atom name strings.

**Why this is debt rather than a bug.** The system works because every
load path either (a) translates to AMBER convention (trajectory path) or
(b) happens to be supplied AMBER-convention input (reduce / tleap on the
non-trajectory paths). The string-keyed dispatch is a load-time boundary
operation and the bounded universe of inputs has historically been
AMBER-convention.

**Why it would benefit from rework.** Every site is a place where future
code or future inputs in non-AMBER conventions break silently. Future
science work that involves CHARMM-source proteins via `--protonated-pdb`,
or ROSETTA outputs, or any non-AMBER convention, would land on the bug.
The architectural thinking on the preservation branch proposes a typed
topology layer with bidirectional projection-kit registry that
structurally prevents this class of bug.

**Architectural debt context.** Full inventory at
`architectural-thinking-2026-04-27:spec/AUDIT_NAMING_2026-04-26/PHASE_1_ATOM_FIELD_READS.md`.

---

## Bug 3 — Deferred `NamingRegistry` rules for CHARMM↔IUPAC gaps

**Severity.** Affects atoms produced by the `--trajectory` path on
specific positions in specific residues. Not a runtime crash; produces
wrong labels in the H5 output for those positions, and (less obviously)
mis-assigns the atoms during boundary translation if the wildcard rule
fires incorrectly.

**Where.** `src/NamingRegistry.cpp:172-307` (the commented-out rule
block).

**What's missing.** Eight categories of CHARMM↔IUPAC translation rules
were left commented out pending fleet-wide vetting. They cover:

- ILE γ-carbon, δ-methyl, γ1-methylene atom name translations.
- δ-methylene rules for ARG, LYS, PRO.
- ε-methylene rule for LYS.
- α-methylene rule for GLY.
- ALA β-methyl identity blockers (the wildcard β-methylene rule
  over-fires on ALA, producing duplicate `HB3` labels).

**Operational guidance.** For 10-protein calibration work, two typed
Python data constants in `h5-reader/notes/nmr_forensics/pack_experimental_shifts.py`
handle the affected cases in-place (per the comments in
`NamingRegistry.cpp`). Fleet-wide vetting is needed before activating
the rules in the C++ registry.

**Architectural debt context.** Documented in
`architectural-thinking-2026-04-27:spec/ChangesRequiredBeforeProductionH5Run.md`.

---

## Design choice 1 — `MutationDeltaResult` uses KD-tree position matching

**Severity.** Not a bug. A known design choice.

**Where.** `src/MutationDeltaResult.cpp` (pre-IUPAC implementation,
which is what the current tree has).

**What it does.** Cross-protein WT-to-mutant atom matching is done by
KD-tree on positions with a 0.5Å tolerance. This is the matcher that
produced the existing 723-protein calibration corpus.

**Why this is a design choice rather than a bug.** It works. The
0.5Å tolerance is tight enough to avoid false matches in practice; the
calibration results are stable. A typed-identity matcher would be more
elegant (matches by chemistry, not geometry), but the typed-identity
work is what the IUPAC topology landing was attempting and what the
preservation branch's architectural thinking targets. Until that
rebuild lands, the KD-tree matcher is the working solution.

---

## Pointer to the architectural thinking

The architectural design from the 2026-04-27 conversation is preserved
as a memory entry — `project_proteintopology_architecture` — in the
user's persistent memory store. That memory is loaded automatically
at session start regardless of branch state and is the operational
reference if the architecture is ever revisited.

In one sentence: the design is a typed-topology-as-graph container
(`TopologyProtein` ABC, with `IUPACTopologyProtein` as one concrete
subclass), 1:M with Protein, holding inhabitants at six granularities
(per-atom, per-pair, per-set, per-residue, per-sequence, per-protein);
`NamingRegistry` becomes a bidirectional projection kit per tool
authority; Atom is structurally clean (no IUPAC fields, no name strings);
and math operates on typed enum categories rather than string spellings
of atom names. See the memory entry for the full distillation.

A historical archive of the working documents (PROTEIN_TOPOLOGY_ABC_GRAPH.md,
the audit phase reports, the session-start protocol draft, the
run-state tracking template, and the 12 reverted IUPAC commits) lives
on branch `architectural-thinking-2026-04-27` on origin. The user has
indicated this branch will not be read again; the memory entry is the
durable reference.

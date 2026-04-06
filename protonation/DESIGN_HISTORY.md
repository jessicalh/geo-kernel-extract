# Protonation Design History

**Date**: 2026-04-02

A record of how the protonation design was intended, what got built
instead, and what we learned from the gap analysis. This document
exists because design decisions have reasons, and reasons that aren't
recorded become myths that the next person either follows blindly or
discards carelessly.

---

## The original vision

In structural biology, most PDB crystal structures arrive without
hydrogen atoms. Protonation is always an issue when modeling proteins.
Systems like this one, which do a one-way analysis of a protein's
physics properties, must be able to analyse the protein based on
protonation via various modeling systems — PROPKA or KaML for
example — and then recalculate the physics properties based on that
change. But since modeling happens in an instantiated set of
structures with implicit geometry and physics, you either create a
generality we do not want to support or become locked to a single
path. And while a single path is fine, it also fails to collect or
make orthogonal, as methods are added to protonate, the properties
like naming which must be tracked.

The core vision: a protein object contains geometry-invariant
properties, and conformations contain the real geometry, whether of
a PDB or MD frame poses. It was considered unlikely we would do
CpHMD in this project, so it was OK if protonation stuck to the
protein object — but there must be a pattern to support deprotonation
and reprotonation, with good C++ semantics. And ideally, even if
things turned from 1 protein instance with 1 protonation enforced
and 100 poses, to 100 protein instances with different protonations
enforced, that was OK too if we ended up handling CpHMD results.

The idea was that if you wanted to make a new protonation state, you
would use a builder. The builder would take the old protein, apply a
constructor that gave you back a hydrogen-free copy, create a
protonation state using KaML or PROPKA or whatever — even a pose's
old protonation state (classed; the fact this was not clearly defined
may have contributed to the mess) — in the state according to the how,
with the state holding receipts (not hydrogens). You would also supply
the old conformation (or conformations) and get them back with your
geometry but absolutely no calculation results, and any geometry
changes required for the hydrogens — you have to run all the
calculations again but you have your start state.

## What got built instead

Two loaders descended from the builder concept, each doing half the
job:

- **LoadProtein** (PdbFileReader): parses a PDB, builds atoms from
  whatever's in the file, no protonation step.
- **LoadOrcaRun** (OrcaRunLoader): reads an AMBER prmtop (protonated
  by tleap externally), builds the complete protein.

Both share the same final steps (create Protein, add atoms/residues,
FinalizeConstruction, create conformation) which are the builder's
body, written in two places.

The protonation infrastructure was built correctly: ProtonationState
(value type with decisions), Protonator interface (PROPKA, KaML),
NamingRegistry (translates between tool naming universes), per-variant
metadata on AminoAcidType. All tested. But nothing connects decisions
to protein construction.

## Where the design didn't survive

**1. "One ProtonationState per ProteinConformation."** The
constitution says this. The code says protonation determines the atom
list, which is per Protein. Different protonation = different atoms =
different Protein. This is correct. The constitution text is stale.

**2. Copy-and-modify.** The concept was right (keep the old protein,
produce a new one) but the mechanism was wrong. A copy constructor
that "applies" a ProtonationState would need to add/remove H atoms,
changing the atom count and invalidating every index. The real
mechanism is a builder that takes heavy-atom data and constructs
fresh. Not a copy — a rebuild.

**3. Protein is non-copyable.** Early agents deleted copy/move
entirely for pointer safety (conformations hold raw Protein*
back-pointers) instead of providing a Copy() factory method. The
builder doesn't actually need a copy constructor — it takes source
data, not a Protein reference.

## The two indexing systems

A Protein has two indexing systems. Seeing them clearly is the key
to understanding why reprotonation requires a rebuild.

**Symbolic**: residue sequence, atom names within residues, bond
categories, ring types. The chemistry. Survives reprotonation because
the heavy-atom skeleton doesn't change. Residue 17 is still HIS. Its
CG, ND1, CD2, CE1, NE2 are still there.

**Geometric**: the flat atom index array, 0 to N-1. Every ring
vertex, every bond endpoint, every ConformationAtom, every spatial
neighbour list, every calculator accumulation — all indexed into this
one flat array. Breaks on reprotonation because adding/removing H
atoms renumbers everything.

The loaders fuse these two systems at construction time: symbolic
identity (from PDB/prmtop names) gets mapped to geometric indices
(by the order atoms are added). After FinalizeConstruction, the
symbolic system is encoded IN the geometric indices — ring.atom_indices
IS the imidazole, expressed as positions in the flat array.

The builder's job is exactly this: take symbolic identity from the
old protein, take new protonation decisions, run tleap to get the
new geometric atom list (with H), and fuse them into a new flat
index space. The symbolic-to-geometric mapping happens fresh every
time. That's why it's a rebuild, not a patch.

Conformations live entirely in the geometric index space. When you
make a new Protein with a different flat array, the old conformations
are meaningless — atom 47 in the old protein is not atom 47 in the
new one. But heavy-atom positions can transfer because symbolic
identity (residue 17, atom CA) maps to a position in both the old
and new geometric index spaces. The name-within-residue lookup
bridges them.

## Findings from the gap analysis

Two rounds of agent analysis (2026-04-02) examined the builder gap
in detail. Key findings:

**Disulfide-before-protonation ordering.** Disulfide bonds must be
identified from SG-SG geometry before protonation decisions are
applied. If PROPKA predicts "deprotonate this CYS" but it's in a
disulfide, writing CYM instead of CYX produces wrong topology. The
builder must detect disulfides first, lock them as CYX, then apply
protonation decisions to remaining free cysteines.

**FinalizeConstruction circularity.** Current order: backbone cache →
ring detection → bond detection. Ring detection needs variant_index
(for HIS). Protonation detection needs bonds (for CYX disulfide
check). For the builder's "as-is" path (loading an externally
protonated structure with unknown provenance), this is circular. The
builder would need two-phase finalization: backbone → bonds →
protonation detection → rings.

**AminoAcidType atom list is single-variant.** The HIS entry lists
HIE atoms (HE2 but not HD1). Only IupacAtomIdentity uses this for
classification. No construction or calculator code depends on it.
Minor enrichment gap (HD1 on HID gets default properties), no
functional impact.

**Coordinate matching is tractable.** Match heavy atoms by
(residue_index, atom_name) between source and prmtop output. tleap
preserves heavy atom names and does not reorder them within residues.
H positions come from tleap's inpcrd. Heavy positions come from the
source. This is a dictionary lookup, not geometric matching.

**The NamingRegistry data chain is complete but uncalled.** The path
variant_index → AminoAcidType.variants[idx].registry_key →
NamingRegistry.ResolveForTool(canonical, Amber, key) → AMBER name
exists end-to-end. Nobody calls it. The builder would be the first
consumer.

**tleap failure modes.** ARN (deprotonated arginine) and TYM
(deprotonated tyrosinate) have no ff14SB parameters. PROPKA can
predict these at extreme pH. The builder needs pre-flight validation
against the force field's known residue set.

## What this means for the builder

The builder is not a rewrite. It unifies two existing loaders under
a single factory that adds the protonation step they both lack. The
core code (atom creation, residue creation, FinalizeConstruction,
conformation creation) already exists in LoadProtein and
LoadOrcaRun. The protonation infrastructure (ProtonationState,
Protonators, NamingRegistry) already exists and is tested. The
builder connects them.

The work is: extract the shared construction steps into the builder,
add protonation as a requestable step between "I have atoms" and
"call FinalizeConstruction," handle the two-phase finalization for
the as-is path, add tleap invocation with NamingRegistry translation,
add coordinate matching for heavy atom position transfer, and add
pre-flight validation for force field coverage.

All 8 calculators, all 274 tests, the Pipeline, the SampleAt
methods, and the UI surface are downstream and need no changes. They
receive complete Proteins and don't care how they were built.

## References

- `spec/BUILDER_AGENT_PROMPT.md` — the evaluation task as given
- `spec/BUILDER_ANALYSIS.md` — round 1 analysis (3400 words)
- `spec/BUILDER_ANALYSIS_CRITIQUE.md` — critique of round 1
- `spec/BUILDER_ANALYSIS_ROUND2.md` — supplementary analysis (2500 words)
- `spec/PROTEIN_BUILDER_SPEC.md` — the builder design spec
- `spec/CONSTITUTION.md` — the supreme constraint (with stale text noted)

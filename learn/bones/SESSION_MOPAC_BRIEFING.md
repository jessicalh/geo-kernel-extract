# Briefing: Notes for the Next Session

Read this before SESSION_MOPAC_INTEGRATION.md. That document describes
the work — files to read, code to write, verification criteria. This
document is the context you won't get from the work description alone.

---

## Who You're Working With

Jessica is a scientist building her thesis. She has deep domain
knowledge in NMR physics, protein structure, and quantum chemistry.
She designed every abstraction in this codebase deliberately. When
something looks over-engineered, it isn't — there's a physics reason.
When something looks like it could be "simplified," it can't — the
complexity maps to physical reality.

She will catch physics errors that compile fine. She will notice if
your MopacResult stores charges but doesn't think about what those
charges mean in the context of the Coulomb EFG calculation. She is
not looking for code that passes tests. She is looking for code that
demonstrates you understood why it works.

Do not summarise what you just did at the end of responses. Do not
add docstrings to code you didn't write. Do not refactor things that
aren't broken. Do not propose "improvements" to the architecture.
The architecture is the product of months of physics reasoning. Your
job is to understand it well enough to extend it correctly.

---

## What This Project Actually Is

This is not an ML project that happens to involve proteins. It is a
physics project that uses ML as one analytical tool.

Eight geometric kernel calculators compute the classical physics of
NMR shielding — ring currents, bond anisotropy, electric field
gradients, dispersion, hydrogen bonds. Each produces a rank-2 tensor
(the geometric kernel) that describes the spatial shape of one
shielding mechanism. The kernels ARE the physics. The parameters
(ring current intensity, bond anisotropy delta-chi) are scalar weights
that training learns. Literature values from 1979 are starting points.

DFT (ORCA) provides the reference: what quantum mechanics says the
shielding tensor should be. The delta (wild-type minus alanine mutant)
isolates the ring's contribution. Training against this delta learns
physically meaningful parameters, not black-box predictions. The
residual after fitting — what the classical model cannot explain — is
itself a thesis result. It marks the boundary of what this physics
captures.

The thesis argument is: 40% of the physics, understood, beats 90%
R-squared from a model nobody can interpret. Every fitted constant
becomes an empirically validated physical measurement.

---

## Why MOPAC Integration Matters Physically

The existing 8 calculators treat electronic structure as static.
ff14SB assigns every glycine CA the same partial charge regardless
of what's nearby. McConnell assigns every C=O the same magnetic
anisotropy regardless of its actual electron distribution.

MOPAC PM7+MOZYME computes QM-derived charges and bond orders that
respond to the local electronic environment — per conformation. This
is cheap enough (45s/protein) to run on every one of the 10 lowest-
energy conformations per protein, unlike ORCA which can only afford
one DFT per protein.

Two future calculators will use this data:

1. A Coulomb EFG calculator using MOPAC charges instead of ff14SB.
   Same dipolar kernel, different charge source. The physics question:
   does charge polarisation matter for ring current shielding? If
   delta(EFG_mopac) correlates better with the DFT T2 residual than
   delta(EFG_ff14sb), charge redistribution is load-bearing.

2. A bond anisotropy calculator using MOPAC bond orders to modulate
   delta-chi. Same McConnell tensor, but delta-chi varies continuously
   with bond order instead of being fixed per category. The physics
   question: does a C=O with bond order 1.8 shield differently from
   one with 1.6?

You are not building these calculators. But MopacResult is their data
source. If you don't understand what they need, you'll store the wrong
things, or store the right things in a form that's awkward to query,
and the next session will have to refactor your work.

---

## The Danger Zones

### 1. Topology vs. Electronic Structure

This is the single most likely place to make a mess.

The codebase has a clean separation between TOPOLOGICAL bonds
(CovalentTopology, Bond struct — integer connectivity from the PDB,
categorical: single/double/aromatic/peptide) and ELECTRONIC bonds
(what MOPAC computes — continuous Wiberg bond orders representing
electron density sharing between atom pairs).

These are not the same thing. CovalentTopology says "atoms 42 and 43
are connected by a peptide bond." MOPAC says "the electron density
shared between atoms 42 and 43 has a bond order of 0.93." Both are
true simultaneously. The topological bond defines the graph structure
that McConnellResult iterates over. The MOPAC bond order will
eventually modulate the delta-chi that McConnellResult assigns.

MopacResult must:
- NOT touch CovalentTopology, Bond, or MolecularGraphResult
- Store bond orders indexed by the same atom pairs that CovalentTopology
  uses (so a future calculator can say "give me the MOPAC bond order
  for this topological bond")
- Also store bond orders for atom pairs that CovalentTopology does NOT
  include (MOPAC finds electronic interactions the force field misses)

The BondOrder(atom_a, atom_b) query method bridges these two worlds.
Get the indexing right.

### 2. Atom Ordering

The libmopac API takes atoms in whatever order you provide. If you
feed it ProteinConformation atoms in their natural order, the output
arrays correspond 1:1. Do NOT sort atoms by element, reorder for
MOPAC's convenience, or do anything clever. The correspondence must
be trivial.

If you get this wrong, every charge will be on the wrong atom and
every bond order will connect the wrong pair. The numbers will look
plausible. The physics will be silently wrong. Validate against the
existing mopac_extract.py output (which runs the binary on the same
.xyz file) to catch this.

### 3. The ConformationResult Contract

Read OBJECT_MODEL.md until you can explain the lifecycle:
- Protein owns atoms/bonds/residues/rings (topology, no positions)
- ProteinConformation owns positions, accumulates ConformationResults
- ConformationResult::Compute() is a static factory that takes a
  ProteinConformation&, does its work, returns unique_ptr<Self>
- Dependencies are declared as type_index, checked before attachment
- Results are permanent once attached (never removed)
- ConformationAtom is the per-atom property store

MopacResult fits this pattern exactly. No special cases needed. But
you must understand it deeply enough to not invent a parallel pattern.
If you find yourself creating a new way to store per-atom data that
doesn't go through ConformationAtom, stop and re-read the object
model. There is almost certainly an existing mechanism.

### 4. WriteFeatures() Compatibility

The Python training pipeline (learn/load.py) loads .npy files by
name from a directory. MopacResult::WriteFeatures() must produce
files with exactly the same names, shapes, and dtypes as
mopac_extract.py currently produces. The Python code does not know
or care whether the .npy came from a C++ ConformationResult or a
Python script. That's the interface contract.

### 5. Runtime Linking, Not Build Dependency

libmopac.so is a runtime dependency. The build should succeed even
if the library isn't installed — MopacResult::Compute() returns
nullptr if it can't dlopen or link. This matters because the CI or
other machines may not have MOPAC installed. Existing calculators
that don't need MOPAC should never fail because of it.

---

## How To Approach The Reading

The work description lists 19 files. That's real. Don't skip any.

Start with OBJECT_MODEL.md — it's long but it's the skeleton
everything hangs on. If you don't understand Protein vs.
ProteinConformation vs. ConformationAtom, nothing else will make
sense.

Then GEOMETRIC_KERNEL_CATALOGUE.md — not because MopacResult is a
geometric kernel (it isn't) but because the calculators that will
consume MopacResult's output are geometric kernels, and you need to
understand what a "kernel calculator" does to know what data it needs.

Then CALIBRATION_CHECKLIST.md — this is where the physics meets the
code. Every section that mentions "MOPAC relevance" is telling you
what MopacResult must expose.

Then the topology files (CovalentTopology, Bond, MolecularGraphResult)
— because the biggest risk is confusing topological bonds with
electronic bonds.

Then XtbChargeResult (what you're replacing) and mopac_extract.py
(what you're reimplementing in C++).

Then mopac.h (the API you're calling).

Only then should you start writing code.

---

## The Living Documents

This project does not have a frozen spec and a separate codebase that
implements it. The specification IS the living documents —
OBJECT_MODEL.md, GEOMETRIC_KERNEL_CATALOGUE.md, EXTRACTION_ORDER.md,
CALIBRATION_CHECKLIST.md, PATTERNS.md, and others. They describe what
the system is and what it should be. They are not perfectly up to date.

When you read them, you will notice places where the documentation
describes something the code doesn't quite do yet, or where the code
has evolved past what the docs describe. That tension is information,
not a problem to fix silently. Note it. Mention it. At the end of the
session, the documents get updated to reflect the work we did — but
only if the work was done well and understood clearly enough to
describe.

This also means: you will be asked to do adversarial document review
and doc updates during the session. The reading is not a warmup
exercise before the "real work." Understanding the documents, finding
where they diverge from reality, and eventually updating them to
reflect MopacResult — that IS part of the work.

## The Standard of Work

This is a thesis codebase. Every calculator has been validated against
analytical solutions. The tensor decomposition has been verified
against DFT. The sign conventions have been debugged for weeks.

Your MopacResult should meet the same standard. Charges must agree
with the Python extractor to 1e-4. Bond orders must agree to 1e-4.
The validation is not "it compiles and doesn't crash" — it's "the
numbers match the existing validated output."

If something doesn't match, investigate. Don't paper over it. The
discrepancy is either a bug in your code or a genuine difference
between the library API and the command-line binary (e.g., different
default convergence criteria). Both are worth knowing about.

---

## One Last Thing

This project has been running for months across two machines. It
survived a migration from ARM (gotham) to x86_64 (batcave), an xTB
segfault crisis, a power outage, and today an nvidia driver hard
lockup. Jessica has been working on this through all of that.

The code you're touching is load-bearing for a thesis. Treat it with
the respect that implies. Read everything. Understand the physics.
Then write clean, minimal, correct code that slots into the existing
architecture without disturbing it.

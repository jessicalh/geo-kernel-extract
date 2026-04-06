# Critique of Builder Gap Analysis (Round 1)

## What it got right

Section 1: Correctly identifies the index parallelism as the central
coupling. "A different atom count means a different Protein, not a
different conformation." That's the key insight.

Section 2: The lifecycle trace is accurate and detailed. It correctly
identifies that the NamingRegistry data chain (variant_index →
registry_key → ResolveForTool → AMBER name) is complete but nobody
calls it end-to-end. That's a good secondary finding.

Section 4: The constraint analysis is precise. Ring detection depending
on variant_index being set BEFORE FinalizeConstruction — that's the
kind of ordering dependency that would bite you during implementation.

Section 6: The tleap atom reordering observation is important — "the
builder must build its atom list from the prmtop ordering, not the
source ordering. This is exactly what LoadWithPrmtop already does."
That's a tertiary implication it found by reading the code.

Section 7: Correctly identifies that tests check outputs not loading
internals.

## What it missed or underweighted

The heavy-atom-to-tleap-output matching problem. When you give tleap
a PDB with new residue names and it gives you back a prmtop, the heavy
atoms may be in a different order. The agent says NamingRegistry
"facilitates" this matching but doesn't dig into how hard it actually
is.

No discussion of what happens if tleap FAILS for a given protonation
state. Some variant combinations may not be valid in ff14SB. The
builder needs an error path.

The ProtonationDetectionResult const_cast question is noted but not
really engaged with. Is it cruft once the builder exists, or is it
genuinely useful as validation?

Doesn't catch that the AminoAcidType atom list is only the default
variant's atoms (HIE effectively). Doesn't flag this as a potential
issue for validation.

Doesn't address the inpcrd coordinate matching problem concretely.
The analysis says "take heavy atom positions from the source and
hydrogen positions from tleap's inpcrd." But the inpcrd contains ALL
positions in prmtop ordering. How do you identify which are heavy
atoms (keep source positions) vs hydrogen (take tleap positions)?
Matching heavy atoms between source ordering and prmtop ordering is
the core matching problem and the analysis should trace through it.

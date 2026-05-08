# Q1: Does `CovalentTopology::Resolve` actually use `rings_`?

Bounded follow-up to `ring-infrastructure-audit.md` §9.Q1. Question asked
because Bundle C wants `Protein::FinalizeConstruction` to construct rings
*after* substrate composition; the current ordering passes `rings_` into
`CovalentTopology::Resolve` at FinalizeConstruction line 232, so the
question is whether the rings parameter is a real dependency or removable
dead weight.

## Verdict

**Rings are used. Non-trivially. For exactly one chemistry decision: tagging
bonds as `BondCategory::Aromatic` / `BondOrder::Aromatic`.**

The parameter is not dead weight, but the use is narrow and isolated to a
single classification branch.

## Use sites

File: `src/CovalentTopology.cpp`

| Lines | What it does |
|-------|--------------|
| 87–90 | Builds `std::set<size_t> aromatic_atoms` by flattening every `ring->atom_indices` from the input `rings` vector. This is the only place `rings` is read. |
| 143–146 | Inside the per-bond classification loop: `else if (aromatic_atoms.count(i) > 0 && aromatic_atoms.count(j) > 0) { bond.order = BondOrder::Aromatic; bond.category = BondCategory::Aromatic; }`. Both endpoints in the aromatic set ⇒ aromatic-tagged bond. |

No other code path in `Resolve` touches `rings`. The set is built once
and consulted O(1) per bond.

The classification branch order matters: disulfide check (line 138) wins
over aromatic; aromatic wins over backbone/sidechain analysis (lines 148+).
So the rings dependency is functionally "if the bond is not S–S and both
endpoints are ring members, tag it Aromatic; else fall through to the
backbone/sidechain logic."

## What `rings` is *not* used for

- Bond detection (lines 106–122). OpenBabel `ConnectTheDots` /
  `PerceiveBondOrders` or the geometric covalent-radius fallback decides
  which pairs become bonds. Rings have no influence here.
- H parent assignment (lines 211–225). Pure nearest-bonded-heavy-atom
  geometry. No rings.
- Backbone/sidechain/peptide categorization (lines 148–197). Uses
  `is_backbone[]` and `Residue` indices; never consults rings.
- Storage. `CovalentTopology` keeps no reference to `rings` after
  `Resolve` returns; the `aromatic_atoms` set is a stack local that
  dies with the function.

## Bundle C implication

**Split. Cannot just drop the parameter.**

The aromatic categorization is a real chemistry output — downstream
calculators (any code consulting `BondCategory::Aromatic` /
`BondOrder::Aromatic`) depends on it. Removing the rings input would
demote every aromatic bond to `BackboneOther` / `SidechainOther`,
silently changing classification for HIS / PHE / TYR / TRP and any
substrate ring (e.g. Pro proline ring once Bundle C extends it).

The cleanest split, consistent with the existing structure of `Resolve`:

1. **Pre-rings pass** — bond detection + non-aromatic categorization.
   Lines 80–122 (init + detection) and lines 127–141, 148–205 (the
   non-aromatic branches of the classification loop), plus the H parent
   assignment block (lines 211–225). Default category for ring atoms
   would be the fall-through (`SidechainOther` / `BackboneOther`).
   Returns the topology in a "rings-naive" state.
2. **Post-rings pass** — a small overlay method
   (e.g. `CovalentTopology::TagAromaticBonds(const vector<unique_ptr<Ring>>&)`)
   that walks already-stored bonds, builds the same `aromatic_atoms` set
   from the now-final rings collection, and rewrites
   `category`/`order` for any bond whose both endpoints are in the set.
   Logic is line 87–90 + the test at line 143–146, lifted as-is.

This mirrors the pattern Protein already uses post-Resolve:
`OverrideDisulfides` (lines 251–348) runs after `Resolve` and applies
authoritative disulfide chemistry on top of the geometric inference.
`TagAromaticBonds` would be the same shape — overlay, not core
detection.

The split is mechanically simple because `aromatic_atoms` is built
fresh each call from `rings` (no incremental state), and the
classification branch order means demoting aromatic categorization
to a second pass does not interact with the disulfide / peptide /
sidechain branches that are kept in the first pass.

## Open questions

None for the bounded question. Bundle C drafting can proceed with the
split-pass plan.

One implementation note for whoever drafts the split: the disulfide
branch currently runs *before* the aromatic branch (line 138 vs
line 143). If the split puts disulfide tagging in pre-rings and
aromatic tagging in post-rings, an SG–SG bond between two atoms
that are also in some hypothetical aromatic ring would stay
`Disulfide` (not get demoted to `Aromatic`). That's the same behavior
as today (disulfide branch wins via `else if` ordering), and
biologically irrelevant — no aromatic ring contains an S that bonds
S — but worth a sentence in the post-rings overlay's docstring so
the precedence is documented rather than incidental.

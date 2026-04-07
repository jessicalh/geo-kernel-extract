# MOPAC Integration Handoff

**Written:** 2026-04-05, end of UI session 4.
**For:** The next session(s) doing MOPAC integration + 2 new
calculators + CalculatedRegion.

## What just happened (UI session 4)

The viewer was refactored to read the library Protein directly —
no adapter layer. An atom inspector was added (double-click →
full object model in a tree view). The viewer is ready to display
whatever new data you add to ConformationAtom or Bond.

**This means:** when you add MOPAC bond orders to Bond, or MOPAC
charges to ConformationAtom, the viewer will show them
automatically in the atom inspector. No viewer code changes needed
for the data to be *accessible*. Display (tube thickness, coloring)
is a separate step for the UI session after yours.

## What MOPAC needs to provide

Per the memory file and OBJECT_MODEL.md:

1. **Mulliken charges** — per-atom, conformation-dependent.
   Field on ConformationAtom (like `xtb_charge`).

2. **Wiberg bond orders** — per-bond, continuous (0.0–3.0).
   This is new. Bond.h has `BondOrder order` (discrete enum).
   Wiberg indices are a float — needs a new field. Options:
   - Field on Bond itself (but Bond is topology, fixed per protein)
   - Field on ProteinConformation (vector<double> parallel to bonds)
   - Per-atom BondNeighbourhood entry (like McConnell's)

   Since bond orders are conformation-dependent (they change with
   geometry), they belong on the conformation, not on Bond. A
   `vector<double> mopac_bond_orders` on ProteinConformation
   (parallel to protein.Bonds()) is the natural home. Or on
   ConformationAtom's bond_neighbours if per-atom granularity
   is wanted.

3. **Heat of formation** — per-conformation scalar.
   Field on the MopacResult itself.

4. **HOMO-LUMO gap** — per-conformation scalar.
   Field on MopacResult (xTB already stores this per-atom;
   MOPAC's is per-molecule).

## Architecture pattern

Follow XtbChargeResult exactly:
- `MopacResult : ConformationResult`
- Static `Compute(conf)` factory
- Writes to ConformationAtom fields (charges)
- Writes to conformation-level storage (bond orders, energy)
- Dependencies: GeometryResult (needs positions for input file)
- External tool: `/home/jessica/micromamba/envs/mm/bin/mopac`
- Keywords: `PM7 MOZYME 1SCF CHARGE=N BONDS MULLIK LET GEO-OK THREADS=8`
- MOZYME: linear-scaling SCF, ~45s per 889-atom protein
- Output parsing: Mulliken charges from `.out`, bond orders
  from BONDS section

## CalculatedRegion (for later, but design now)

The thesis needs to track WHERE each calculation was applied and
what it found. A CalculatedRegion describes a spatial region of the
conformation where a specific calculator was evaluated.

This intersects with the existing KernelFilterSet — filters already
decide which atom-source pairs to evaluate. CalculatedRegion would
record the result of that decision: "McConnell was evaluated for
atoms 42–67 against bonds 100–150, rejecting 12 pairs via
DipolarNearField filter."

Think about this while integrating MOPAC. The pattern might be:
- CalculatedRegion as a ConformationResult or as data on existing
  ConformationResults
- Each calculator records its evaluation footprint
- The viewer can then show "where was this calculator active?"

## The 2 new calculators

I don't know what they are yet. They'll be decided in the MOPAC
sessions. They follow the same ConformationResult pattern as the
existing 8. The viewer needs:
- 2 more physics checkboxes in the sidebar
- 2 more terms in checkedCalcT0/checkedCalcST
- 2 more sections in populateAtomInfo
- New SphericalTensor fields on ConformationAtom

## What the UI session after yours needs

See `ui/NEXT_UI_SESSION.md` — that's the message in a bottle for
the future UI session. It tells them what to read, how the
architecture works, and how MOPAC data should be displayed (bond
tube thickness, charge comparison, CalculatedRegion panels).

## Reading order for MOPAC session

1. spec/INDEX.md → reading order
2. OBJECT_MODEL.md — what exists, where new types fit
3. PATTERNS.md — how to implement ConformationResults
4. src/XtbChargeResult.h and .cpp — the template to follow
5. src/OperationRunner.cpp — where to add MOPAC to the sequence
6. src/Bond.h — where bond orders might go
7. src/ProteinConformation.h — where conformation-level data goes
8. src/ConformationAtom.h — where per-atom fields go

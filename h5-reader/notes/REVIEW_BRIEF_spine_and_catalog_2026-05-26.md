# Review brief — h5-reader topology spine, ring typing, field catalog

**For:** a second set of eyes (codex or any reviewer). **Date:** 2026-05-26.
**Tone:** this is a *fidelity* review, not a bug hunt — and a friendly one.
Thank you for reading carefully; the areas below are detailed, and a fresh
perspective is exactly what we want.

## What this is

`h5-reader/` is a standalone Qt6/VTK reader that is the **read-side mirror**
of the `nmr_shielding` C++ library (`../src/`). It re-declares the library's
typed objects — atoms, residues, bonds, the aromatic-ring class hierarchy,
the chemistry substrate — so it can interpret the library's serialised output
(per-frame NPYs + a topology sidecar) **in topology terms** rather than as
anonymous number arrays. The library and `../GEOMETRIC_KERNEL_CATALOGUE.md`
are **ground truth**; the reader must mirror them faithfully. A mismatch here
silently mis-interprets every value downstream, so it's worth a careful look.

## How to help us most

- Cite **`file:line` on both sides** (reader + library/catalogue), with
  **expected vs actual** and a one-line *why it matters*.
- **Uncertainty is welcome** — "this looks off, please check X" is genuinely
  useful. We verify every finding deterministically against the source before
  changing anything, so precise pointers help more than confident prose.
- **Peripheral observations are explicitly invited.** If something catches
  your eye outside the checklist below — naming, a sign, an off-by-one, a
  stale comment, an ordinal that smells wrong — please note it. The checklist
  is where we *think* the risk is; you may see what we didn't.
- Read-only. No changes needed from you.

## Area 1 — Ring typing

**Reader:** `h5-reader/src/model/QtRing.h` — nine subclasses
(`QtPheBenzeneRing` … `QtProPyrrolidineRing`) with virtuals
`LiteratureIntensity()` (nA/T, Giessner-Prettre 1969), `JohnsonBoveyLobeOffset()`
(Å), `NitrogenCount()`, `Aromaticity()`, `RingSizeValue()`, `TypeIndex()`,
`TypeName()`.

**Ground truth:** the library's ring-type classes in `../src/` (the `Ring`
hierarchy / `RingTopology`) and `../GEOMETRIC_KERNEL_CATALOGUE.md`.

**Please check** each value matches. Reference values from the catalogue:

| Ring | Intensity (nA/T) | Lobe offset (Å) | N-count | Ring size |
|------|------------------|-----------------|---------|-----------|
| PHE benzene | −12.0 | 0.64 | 0 | 6 |
| TYR phenol | −11.28 | 0.64 | 0 | 6 |
| TRP benzene | −12.48 | 0.64 | 0 | 6 |
| TRP pyrrole | −6.72 | 0.52 | 1 | 5 |
| TRP perimeter (indole) | −19.2 | 0.60 | 1 | 9 |
| HIS / HID / HIE imidazole | −5.16 | 0.50 | 2 | 5 |
| PRO pyrrolidine (saturated) | 0.0 | 0.0 | 1 | 5 |

Identity to confirm: `I(perimeter) == I(pyrrole) + I(benzene)`
(−19.2 == −6.72 + −12.48). Aromaticity should read Full for the six-membered
aromatics + perimeter, Reduced/Weak for pyrrole/imidazole, None for pyrrolidine.

## Area 2 — Detailed topology / chemistry substrate

**Reader:** `h5-reader/src/model/{QtAtom.h, QtResidue.h, QtBond.h, QtTopology.h,
QtAminoAcidType.*, QtSemanticEnums.h, Types.h}`.
**Ground truth:** `../src/{Types.h, SemanticEnums.h, the AtomSemanticTable,
AminoAcidType, CovalentTopology}`.

**Please check:**

1. **Enum ordinal compatibility (highest stakes).** Every enum in the reader's
   `Types.h` + `QtSemanticEnums.h` must have the **same ordinal values** as its
   library counterpart — the NPY/H5 stores raw ordinals and the loader casts
   them straight back, so any drift mis-decodes silently. Check `Element`,
   `BondOrder`, `BondCategory`, `RingTypeIndex`, `RingKind`, `AminoAcid`,
   `DsspCode`, and the substrate enums (`Locant`, `BackboneRole`,
   `PlanarGroupKind`, `PolarHKind`, `ProchiralStereo`, `PlanarStereo`,
   `DiastereotopicIndex`, `PseudoatomKind`, `RingPositionLabel`).
2. **Predicate logic.** `QtAtom.h`'s predicates (`IsBackbone`,
   `IsBackboneAmideHydrogen`, `IsAnyAlphaHydrogen` — note the GLY HA2/HA3
   special case, `IsSidechainCarboxylateOxygen`, `IsSidechainAmideOxygen`,
   `IsPolarH`, `IsInAnyRing`) should match the library's `AtomSemanticTable`
   (the header cross-references `../src/SemanticEnums.h` around line 902).
3. **AminoAcid registry.** `QtAminoAcidType` 3-letter / 1-letter codes,
   `is_aromatic` / `is_titratable` / `has_amide_h` / chi-angle count vs the
   library's `AminoAcidType`.

## Area 3 — Field catalog vs calculators

**Reader:** `h5-reader/src/io/QtFieldCatalog.gen.h`, **generated** by
`h5-reader/scripts/gen_field_catalog.py` from `../python/nmr_extract/_catalog.py`.

**Please spot-check** that the generated `FieldSpec` entries (`cols`,
`native_axis`, `irreps`) faithfully reflect `_catalog.py` **and** what the
calculators actually emit (`../src/*Result.cpp` `WriteFeatures`). Examples:
`bs_per_type_T2` cols = 40 (8 ring types × 5), `mc_category_T2` cols = 25
(5 categories × 5), EFG fields cols = 5 (T2-only, symmetric-traceless),
per-type axes ordered by `RingTypeIndex`. The generator is AST-based and
should be faithful — but a sanity read against one or two calculators is
welcome.

## Output

A list of candidate findings (reader file:line | ground-truth file:line |
expected | actual | why it matters | confidence), plus any
corner-of-the-eye notes. Thank you — this is the spine everything else
mirrors, so the careful read is appreciated.

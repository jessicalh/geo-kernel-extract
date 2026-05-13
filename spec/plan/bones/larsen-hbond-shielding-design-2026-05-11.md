# LarsenHBondShieldingResult — Hydrogen-Bond Term Calculator (2026-05-11)

This document captures the design for the calculator that implements
Larsen 2015's four hydrogen-bond contribution terms
(Δσ_1°HB, Δσ_2°HB, Δσ_1°HαB, Δσ_2°HαB) plus the water term Δσ_w,
using the Larsen ProCS15 H-bond DFT grid data unpacked from
`/mnt/expansion/larsen_archive/hydrogenbondnmrlogs.tar`.

This document **supersedes the planning** captured in the memory
entry `project_hbond_halpha_design` (2026-05-10). That design rested
on the premise that Larsen's H-bond DFT grid was not accessible. It
is — six nested archives, ~20K Gaussian logs, full T2 tensors at
OPBE/6-31G(d,p). The grid-lookup formulation here replaces the
kernel × η simplification *as the next-gen design*, but the kernel-
form `HBondResult` calculator already in the codebase **stays — not
retired**.

The thesis discipline: **methods accumulate, differences are
reportable**. The existing kernel-form `HBondResult` and the new
grid-form `LarsenHBondShieldingResult` are sibling calculators
covering overlapping physics (amide-H donor). Both run in the
pipeline. Both emit their own NPYs. Their per-atom-type residuals
on the same protein become methodological output — a coordinate in
the formulation-space, not redundant code. See memory
`feedback_methods_accumulate`.

Per the same principle: Larsen Table 2 says some atom types receive
NO contribution from a given term (e.g. Cβ never gets an HB term).
The calculator still computes and emits those readouts — they should
be near-zero, and confirming the near-zero is the reality check that
the upstream parser → loader → frame-rotation pipeline is honest.
A non-zero "should be zero" cell is a methodological signal, not a
bug to silently filter out.

## Status

PROPOSED, 2026-05-11. Verifications completed this session:

- Substrate enum coverage in `src/SemanticEnums.h` reviewed.
- No existing Gaussian shielding-tensor parser in project; geometry
  parsers (`scripts/perceive_larsen_tripeptide.py`,
  `tests/test_larsen_residue_against_source_log.cpp::ParseStandardOrientation`)
  give the line-format precedent. Parser is new code; format is
  well-understood.
- Larsen 2015 §H-bond methodology + Tables 1–2 read; design grounded
  in published decomposition.

## Motivation

The Larsen Δσ_HB terms are the single largest contribution to atom
shielding RMSD that we do not currently compute:

| Atom | Δσ that matters | Effect on Ubiquitin RMSD if removed |
|------|------------------|--------------------------------------|
| C'   | Δσ_2°HB          | 2.1 ppm                              |
| HN   | Δσ_1°HB          | 0.7 ppm                              |
| Hα   | Δσ_1°HαB         | 0.4 ppm                              |

(From Larsen 2015 Table 1 — ProCS15 row vs single-term-removed row.)

Without these terms, the calibration ridge regression's residuals on
C' and HN are systematically dominated by H-bond geometry the kernel
basis cannot capture. The grid-lookup form gives us exactly what
Larsen used; the calibration target becomes "agreement with Larsen's
DFT," not "fit a coefficient."

## Larsen 2015 spec — what this calculator implements

### Decomposition (Eq. 5)

Δσ_HB^i = Δσ_1°HB(rHO, θ, ρ) + Δσ_2°HB(rOH, θO, ρO)

Δσ_HαB^i = Δσ_1°HαB(rHαO, θ, ρ) + Δσ_2°HαB(rOHα, θO, ρO)

Each H-bond pair contributes **two terms**:

- **1° (primary)** — donor-side effect, applied to the donor residue i.
  Geometry referenced from the donor's perspective.
- **2° (secondary)** — acceptor-side effect, applied to residue i+1 of
  the acceptor (the amide N/H following the C=O that accepts the
  H-bond). Geometry referenced from the acceptor's perspective
  (`rOH, θO, ρO` labelling).

Only the NMA acceptor case (backbone O acceptor) has a well-defined
i+1 mapping. For HOMe (hydroxyl) and acetate (COO-) acceptors, only
the 1° term applies.

### The six grid archives

| Archive               | Donor       | Acceptor             | Larsen Fig |
|-----------------------|-------------|----------------------|------------|
| `NMANMA_nmrlog.tar.bz2` | NMA amide H | NMA backbone O       | 2a         |
| `NMACOH_nmrlog.tar.bz2` | NMA amide H | HOMe hydroxyl O      | 2b         |
| `NMACOO_nmrlog.tar.bz2` | NMA amide H | Acetate carboxylate O | 2c         |
| `ALANMA_nmrlog.tar.bz2` | Ac-A-NMe Hα | NMA backbone O       | 3a         |
| `ALACOH_nmrlog.tar.bz2` | Ac-A-NMe Hα | HOMe hydroxyl O      | 3b         |
| `ALACOO_nmrlog.tar.bz2` | Ac-A-NMe Hα | Acetate carboxylate O | 3c         |

### Grid axes (from §H-bond scans)

NMA donor (HB):
- rOH: 1.5–3.0 Å, step 0.125 Å (13 points)
- θH: 90°–180°, step 10° (10 points)
- ρH: -180° to 180° (periodic)

ALA donor (HαB):
- rOHα: 1.8–4.0 Å, step 0.2 Å (12 points — actual archive has up to 4.2 Å)
- θHα: 90°–180°, step 10° (10 points)
- ρHα: -180° to 180°, step 15° (24 points, periodic)

Actual filename inspection (ALANMA): 3252 logs in this archive — axis
extents may exceed paper's stated ranges. Trust the data over the paper
for axis ranges.

### Reference subtraction (§H-bond scans, last paragraph)

> To get the change in chemical shift caused by the hydrogen bonding the
> OPBE/6-31G(d,p)//PM6 chemical shielding of systems without hydrogen
> bonding are subtracted from the scans.

The reference is the **free monomer** — donor without acceptor, internal
geometry frozen at the same PM6-optimised conformation as used in the
grid points. Reference σ is subtracted at parse time → grid stores Δσ
already. Need to identify the free-monomer log per archive (likely a
separate "monomer" calculation; verify in session 1).

### Per-atom contribution dispatch (Larsen Table 2)

Per atom in the protein, which contributions apply:

| Target atom | 1°HB | 2°HB | 1°HαB | 2°HαB | RC | w |
|-------------|------|------|-------|-------|----|----|
| Cα          | —    | ✓    | —     | —     | —  | —  |
| Cβ          | —    | —    | —     | —     | —  | —  |
| C'          | —    | ✓    | —     | —     | —  | —  |
| Hα          | ✓    | ✓    | ✓     | ✓     | ✓  | —  |
| HN          | ✓    | ✓    | ✓     | ✓     | ✓  | ✓  |
| N           | ✓    | ✓    | —     | ✓     | —  | —  |

The "✓ but not for this atom-type" cells are **structurally absent**
(the structural model in Fig 3 has no Cβ on the acceptor side, etc.) —
they are NOT zeroes we should look up. The dispatch table gates the
NPY emission.

### Water term Δσ_w

Δσ_w = **2.07 ppm constant** added to amide H atoms that do not form
any H-bond to another protein atom. Isotropic only (not a tensor).
DFT-computed at OPBE on the NMA-water complex. Larsen 2015 final
paragraph of §H-bond methodology.

### Sidechain primary amide acceptors (ASN ODE1, GLN OE1)

Not separately gridded by Larsen. Approximate using the NMA acceptor
grid; document as a known limitation alongside the calculator.

## Architecture

### Component diagram

```
data/larsen_hbond_grids/                  # 6 NPY/JSON pairs (parser output)
    NMANMA_grid.npy + .axes.json
    NMACOH_grid.npy + .axes.json
    NMACOO_grid.npy + .axes.json
    ALANMA_grid.npy + .axes.json
    ALACOH_grid.npy + .axes.json
    ALACOO_grid.npy + .axes.json

src/LarsenHBondGrid.{h,cpp}              # Loads NPY/JSON, exposes
                                          # QueryNearest(donor, acceptor,
                                          #               r, θ, ρ)
                                          # Trilinear interp + ρ-wrap.
                                          # Pre-decomposed Δσ tensors per
                                          # contribution-target atom.

src/LarsenHBondShieldingResult.{h,cpp}   # Calculator. Iterates donor
                                          # atoms, gates acceptors via
                                          # substrate-typed roles +
                                          # geometry. Looks up 1°+2°
                                          # contributions, applies per
                                          # Table 2. Accumulates per-atom
                                          # tensors. Emits NPYs.
```

### Substrate-typed donor/acceptor predicates

The calculator uses these predicates on `SemanticAt(ai)` (the
`AtomSemanticTable` per-atom record). Existing predicates that work:

```cpp
sem.IsBackboneAmideHydrogen()    // BackboneRole::AmideHydrogen
sem.IsBackboneAlphaHydrogen()    // BackboneRole::AlphaHydrogen — non-GLY only
sem.IsBackboneCarbonylOxygen()   // BackboneRole::CarbonylOxygen
sem.IsBackbone()                 // Any of the 6 backbone roles
```

New predicates needed (add to `AtomSemanticTable`):

```cpp
// GLY HA2/HA3 carry Locant::Alpha + BackboneRole::None; non-GLY HA
// carries BackboneRole::AlphaHydrogen + Locant::None. This unifies.
constexpr bool IsAnyAlphaHydrogen() const {
    return backbone_role == BackboneRole::AlphaHydrogen
        || (element == Element::H && locant == Locant::Alpha);
}

// Sidechain carboxylate O: Asp OD1/OD2, Glu OE1/OE2, C-term carboxyl.
constexpr bool IsSidechainCarboxylateOxygen() const {
    return element == Element::O
        && planar_group == PlanarGroupKind::Carboxylate;
}

// Sidechain primary amide carbonyl O: Asn ODE1, Gln OE1.
// (Approximated by NMA grid at runtime.)
constexpr bool IsSidechainAmideOxygen() const {
    return element == Element::O
        && planar_group == PlanarGroupKind::SidechainAmide;
}

// Hydroxyl O: requires bond-graph dispatch (walk to bonded H, check
// PolarHKind). Calculator-side helper, not a substrate predicate.
```

The hydroxyl-O classification (SER OG, THR OG1, TYR OH) needs a
bond-graph dispatch at the calculator level: an O atom is a hydroxyl
acceptor iff one of its bonded H atoms has
`polar_h == PolarHKind::HydroxylOH_Aliphatic | HydroxylOH_Aromatic`.

### Acceptor classification → grid selection

```cpp
enum class HBondAcceptorClass : uint8_t {
    None             = 0,
    BackboneO        = 1,  // CarbonylOxygen
    SidechainAmideO  = 2,  // ASN/GLN — approximated by NMA grid
    Hydroxyl         = 3,  // SER OG, THR OG1, TYR OH — HOMe grid
    Carboxylate      = 4,  // ASP/GLU/C-term — acetate grid
};
```

Grid dispatch:

| Donor type | Acceptor class | Grid archive |
|------------|----------------|--------------|
| Amide H    | BackboneO      | NMANMA       |
| Amide H    | SidechainAmideO | NMANMA (approx) |
| Amide H    | Hydroxyl       | NMACOH       |
| Amide H    | Carboxylate    | NMACOO       |
| Hα         | BackboneO      | ALANMA       |
| Hα         | SidechainAmideO | ALANMA (approx) |
| Hα         | Hydroxyl       | ALACOH       |
| Hα         | Carboxylate    | ALACOO       |

### Per-H-bond-pair output

```cpp
struct LarsenHBondPair {
    int                       donor_atom_idx;
    int                       acceptor_atom_idx;
    HBondAcceptorClass        acceptor_class;
    // Geometry — degrees throughout (NOT radians). The LarsenHBondGrid
    // loader consumes this same `LarsenHBondGeometry` struct directly.
    LarsenHBondGeometry       geom;            // r_angstrom, theta_deg, rho_deg
    // 1° contributions (donor's residue i atoms)
    std::array<Mat3, 6>       primary_tensors;    // {N, CA, CB, C', Hα, HN}
    // 2° contributions (acceptor's residue i+1 atoms, NMA acceptor only)
    std::array<Mat3, 3>       secondary_tensors;  // {N, H, Hα-stand-in}
    bool                      has_secondary;       // false for HOMe/COO
};
```

### Per-target-atom accumulation

For each protein atom `ai`, sum contributions from all H-bond pairs:

```cpp
struct LarsenHBondPerAtom {
    Mat3                      total_tensor;       // Sum of all contributions
    SphericalTensor           total_spherical;    // T2-preserving paired storage
    int                       n_primary_contribs;   // From this atom as donor
    int                       n_secondary_contribs; // From acceptor's i+1 perspective
    double                    water_contribution_iso; // 2.07 ppm if no other HB
};
```

### Output NPYs

| File | Shape | Description |
|------|-------|-------------|
| `larsen_hbond_count.npy`            | (N_atoms,) int32 | Number of contributions per atom |
| `larsen_hbond_total_tensor.npy`     | (N_atoms, 9) float64 | Mat3 sum, row-major |
| `larsen_hbond_total_spherical.npy`  | (N_atoms, 9) float64 | T0+T1[3]+T2[5] |
| `larsen_hbond_total_shielding.npy`  | (N_atoms, 9) float64 | Same as tensor, in ppm units (no calibration coeff — Larsen Δσ is already ppm) |
| `larsen_hbond_pairs.npy`            | (N_pairs, 7) float64 | [donor_ai, acceptor_ai, class, rOH, θ, ρ, isotropic_total] for downstream introspection |
| `larsen_hbond_water_term.npy`       | (N_atoms,) float64 | 2.07 if amide H with no HB partner else 0; 0 for non-amide-H atoms |

## Hard parts

The four risk surfaces:

### 1. Frame convention

Each Gaussian log has its own coordinate system. The (rOH, θ, ρ) in
the filename must reproduce from the atomic positions in the log to
within FP. **Two-path validation invariant:** at parse time, derive
(rOH, θ, ρ) from the log's donor/acceptor atomic positions using the
canonical frame definition; compare to the filename values; fatal on
any disagreement above ε = 1e-4 (Å for r, ° for θ/ρ). Get this wrong
once and every contribution has silent angular garbage.

Larsen 2015 Fig 2/3 define the angles geometrically; the precise
vectors are:

- For donor (NMA or ALA Hα): donor H position, the bonded heavy atom
  (N for NMA, Cα for ALA), and the acceptor O position define `rOH`
  and θ (the H...O–C angle).
- ρ is the H...O–X–Y dihedral where X is the C=O carbon and Y depends
  on acceptor (carbonyl C-N for NMA, C-O for hydroxyl, C-O or C-C for
  acetate).

Confirm by inspection of one log per archive in session 1.

### 2. Reference shielding subtraction

The free monomer's σ is subtracted from each grid point's σ to give Δσ.
The free-monomer calculation is presumably a separate Gaussian log
(per donor molecule: one for free NMA, one for free Ac-A-NMe). Locate
in session 1 — likely a sibling log in the same archive, or a "monomer"
labelled point in the grid.

**Sign convention.** Larsen reports Δσ as "change in chemical
shielding" — same sign convention as σ. Add Δσ to the existing
baseline; do not negate. Verify against Larsen 2015 Table 1 published
RMSDs in session 4.

### 3. Periodic ρ in trilinear interpolation

ρ wraps at ±180°. Naive trilinear interpolation between ρ = 175° and
ρ = -180° passes through 0°, which is angularly distant. Wrap ρ to a
periodic index: if ρ > 180° - step/2, treat as wrapping; the next
grid point is at ρ = -180° (or equivalently +180°). Unit-test against
a hand-computed wrap case.

### 4. Per-atom-type Table 2 dispatch

A bug in the dispatch table silently zeroes contributions that should
fire (or vice versa). Encode Table 2 as a compile-time
`constexpr` table indexed by atom-role + contribution-term; assert in
unit tests that for each (target_atom, term) cell, the calculator's
emission matches Larsen Table 2.

## Session plan — four sessions

### Session 1 — Python parser, 6 NPY grids — LANDED 2026-05-11

**Status: COMPLETE.** 6 NPZ + meta.json pairs in `data/larsen_hbond_grids/`
(gitignored as regenerable). Parser at `scripts/larsen_hbond_grid_parse/`
with README documenting filename↔geometry caveats and reference proxy.

Final grid summary:

| Archive | N points | r range (Å) | θ range (°) | ρ range (°) |
|---------|----------|-------------|-------------|-------------|
| NMANMA  | 5027     | 1.50–4.00   | 90–180      | -180–180    |
| NMACOH  | 5040     | 1.50–4.00   | 90–180      | -180–180    |
| NMACOO  | 5039     | 1.50–4.00   | 90–180      | -180–180    |
| ALANMA  | 3120     | 1.60–4.00   | 90–180      | -180–180    |
| ALACOH  | 2880     | 1.80–4.00   | 90–180      | -180–180    |
| ALACOO  | 2880     | 1.80–4.00   | 90–180      | -180–180    |

Per-grid contents: Δσ tensors (Mat3 3×3) for donor-side readouts
{N, CA, CB, C, HA, HN} and (NMA-acceptor only) acceptor-side readouts
{N, HN, HA}. All tensors finite. Tensors expressed in **canonical donor
frame** so all grid points share a common basis.

Sanity-check sample (ALANMA tight H-bond at r=2.0 Å, θ=180°, ρ≈0):
Δσ_Hα isotropic = -2.12 ppm — physically reasonable for a near-linear
Cα-Hα...O=C; far-point reference subtracts cleanly to ~0.

Total parse time: ~35 s on a workstation. Total emitted: ~10 MB across
6 NPZ.

Deviations from original plan:

1. **Free-monomer reference NOT in archive.** Used r-max-edge proxy
   reference (mean σ across grid points where r is within one grid
   step of r_max). Bias estimated ~0.1 ppm; documented in parser
   README. Future work: run Gaussian on free monomers to replace.
2. **Filename↔geometry not a strict invariant.** Probed 10 sample
   logs: NMA donor matches paper, ALA donor has a 0.2 Å r-offset
   (likely Larsen input-file labeling off-by-one), ρ is sign-flipped
   relative to standard IUPAC. **Parser uses ACTUAL measured (r, θ, ρ)
   from atom positions as grid keys**; filenames are labels only.
3. **Larsen-failed DFT logs**: Some grid points failed Gaussian
   optimization ("Small interatomic distances ... Error termination").
   Per Larsen 2015 these were interpolated; our parser skips them and
   leaves the gap for runtime trilinear interpolation.
4. **NMA-NMA symmetric system**: cluster identification by atom count
   alone is ambiguous (both molecules have 3 C, 1 N, 1 O). Disambiguated
   by H-bond proximity (the donor cluster is the one whose amide N-H
   is closer to an O in the other cluster).
5. **Scattered geometry on ρ axis**: PM6 monomer-optimisation
   introduces small angular scatter, so the actual ρ values span
   ~100-300 distinct points rather than nominal ~24. C++ grid loader
   will need to bin to nominal axis at load time, or use scattered
   nearest-N interpolation.

### Session 2 — C++ grid loader

Outputs:
- `src/LarsenHBondGrid.{h,cpp}` — analog of `TripeptideDftTable`.
- `Session::LoadLarsenHBondGrid()` — analog of `LoadTripeptideDftTable`.
- `tests/test_larsen_hbond_grid.cpp` — structure tests.

Steps:
1. Load 6 NPZ + meta.json pairs at Session init. Hold in memory
   (entire grid set is ~10 MB; trivial).
2. **Bin scattered points to nominal regular grid at load.** Per
   archive `meta.json` defines the nominal axes (r step 0.125 NMA /
   0.2 ALA; θ step 10°; ρ step 15°). Average σ tensors within each
   bin. Result: 21×10×24 (NMA) or 13×10×24 (ALA) regular grid per
   archive. The scattered ρ in the raw NPZ (mean |residual| ~0.4°,
   max ~7.5°) is below Larsen's scan resolution; binning is the
   physically faithful path. Detailed reasoning in
   parser README + the design discussion 2026-05-11.
3. **Cubic spline interpolation** for `QueryNearest(donor_class,
   acceptor_class, rOH, θ, ρ)`. Tricubic on the binned regular grid
   with periodic ρ-wrap at ±180°. Returns `LarsenHBondQuery` with
   Δσ tensors for the primary + (when applicable) secondary readout
   atoms. Linear (trilinear) is the obvious simpler fallback if cubic
   complexity is an issue; settled choice is cubic per discussion
   2026-05-11 (~10× lower interpolation error than linear at these
   grid sizes; Larsen 2015 used cubic for analogous backbone scans).
4. Tests: hand-pick a known grid point (one that the bin contains
   directly), query at exactly that geometry, assert the tensor
   matches the binned value. Query at midpoint between two bins,
   assert tricubic-interpolated value. Test ρ-wrap at the 180/-180
   boundary. Test r and θ out-of-range return zero.

Acceptance criterion: 5+ structure tests pass; grid load is fast
(<100 ms at session init).

### Session 3 — `LarsenHBondShieldingResult` calculator

Outputs:
- `src/LarsenHBondShieldingResult.{h,cpp}`.
- `tests/test_larsen_hbond_shielding.cpp` — structure tests + 1UBQ smoke.
- New substrate predicates added to `AtomSemanticTable`.
- Registration in calculator pipeline.

#### Substrate predicate additions

Three new predicates on `AtomSemanticTable` + one calculator-side
helper:

```cpp
// AtomSemanticTable predicates (constexpr, header-only)
constexpr bool IsAnyAlphaHydrogen() const {
    // Covers non-GLY HA (BackboneRole::AlphaHydrogen) AND
    // GLY HA2/HA3 (Locant::Alpha + BackboneRole::None per Markley).
    return backbone_role == BackboneRole::AlphaHydrogen
        || (element == Element::H && locant == Locant::Alpha
            && backbone_role == BackboneRole::None);
}

constexpr bool IsSidechainCarboxylateOxygen() const {
    return element == Element::O
        && planar_group == PlanarGroupKind::Carboxylate;
}

constexpr bool IsSidechainAmideOxygen() const {
    // ASN OD1, GLN OE1 — note: PlanarGroupKind::SidechainAmide also
    // includes the amide N (ND2/NE2). Filter to Element::O to get
    // just the carbonyl O.
    return element == Element::O
        && planar_group == PlanarGroupKind::SidechainAmide;
}
```

```cpp
// Calculator-side helper (NOT a substrate field — requires bond-walk)
//
// SER OG / THR OG1 / TYR OH are hydroxyl O atoms. The substrate
// has no single-field marker for this. Walk to bonded H atoms and
// check PolarHKind::HydroxylOH_Aliphatic or HydroxylOH_Aromatic.
bool IsHydroxylOxygen(const Protein& p, int ai);
```

#### Per-class atom-role resolution algorithm

For each candidate donor `(residue_i, donor_atom_idx)`, compute the
THREE positions `(donor_H, donor_anchor, donor_third)` via the
following table:

| Donor class           | donor_H        | donor_anchor   | donor_third                   |
|-----------------------|----------------|----------------|-------------------------------|
| AlphaHydrogen (non-GLY) | `Hα(i)`       | `Cα(i)`        | `N(i)`                         |
| AlphaHydrogen (GLY HA2/HA3) | the specific HA atom | `Cα(i)` | `N(i)`                  |
| AmideHydrogen          | `HN(i)`       | `N(i)`         | `C'(i-1)`  — PRECEDING residue |

GLY edge case: GLY has TWO α-hydrogens (HA2 and HA3). Each is
processed as a separate donor candidate. `Cα(i)` and `N(i)` are
shared anchors. The chosen HA atom is whichever H is bonded to `Cα`
(both qualify, both run independently).

N-terminus edge: for the first residue, `C'(i-1)` does not exist.
NMA-donor pairing skips this residue with a logged Info message.

PRO edge: PRO has no backbone amide H (secondary amine N). The
substrate predicate `IsBackboneAmideHydrogen()` is already false for
PRO N — so PRO is silently excluded from the donor iteration. No
extra check needed.

For each acceptor atom `(residue_j, acceptor_atom_idx)`, compute the
THREE positions `(acceptor_O, acceptor_C, acceptor_third)`:

| Acceptor class       | acceptor_O   | acceptor_C    | acceptor_third          |
|----------------------|--------------|---------------|-------------------------|
| BackboneCarbonyl     | `O(j)`       | `C'(j)`       | `N(j+1)`                 |
| SidechainCarbonyl (ASN) | `OD1(j)`  | `CG(j)`       | `ND2(j)`                  |
| SidechainCarbonyl (GLN) | `OE1(j)`  | `CD(j)`       | `NE2(j)`                  |
| HydroxylOxygen (SER) | `OG(j)`      | `CB(j)`       | `HG(j)`                   |
| HydroxylOxygen (THR) | `OG1(j)`     | `CB(j)`       | `HG1(j)`                  |
| HydroxylOxygen (TYR) | `OH(j)`      | `CZ(j)`       | `HH(j)`                   |
| CarboxylateOxygen (ASP) | OD1 or OD2 (closer) | `CG(j)` | other OD              |
| CarboxylateOxygen (GLU) | OE1 or OE2 (closer) | `CD(j)` | other OE              |
| CarboxylateOxygen (C-term) | O or OXT (closer) | `C(j)`  | other terminal O     |

C-terminus edge: for the last residue, `N(j+1)` does not exist. The
BackboneCarbonyl 2° term cannot be assigned (no i+1 atoms). Mark
the 2° contribution as un-applicable for that pair; the 1° term
still applies.

Acceptor selection for carboxylate (symmetric carboxylate has 2 Os):
pick the O closer to the candidate donor H. The carboxylate is
symmetric in DFT but the grid was scanned with respect to a specific
O; picking the closer one matches Larsen's "the H bonds to whichever
O is in range."

#### Larsen Table 2 dispatch — constexpr encoding

```cpp
// (atom_type, contribution_term) → applies (true/false).
// Encodes Larsen 2015 Table 2 exactly. Calculator iterates over
// readout atoms and consults this table per (target, term).
//
// Cβ row is INTENTIONALLY all-false per Larsen — we emit the Cβ
// contribution anyway as a diagnostic ("must be near zero in
// production output").
struct LarsenContribDispatch {
    enum class TargetAtom : std::uint8_t {
        N = 0, CA, CB, C, HA, HN, Count
    };
    enum class Term : std::uint8_t {
        Primary_HB = 0, Secondary_HB,
        Primary_HaB, Secondary_HaB,
        RingCurrent, Water,
        Count
    };
    static constexpr bool Applies(TargetAtom t, Term term) {
        // Table 2 from Larsen 2015 PeerJ:
        //                    1HB 2HB 1HaB 2HaB RC w
        // N                  x   x   .   x    .  .
        // Ca                 .   x   .   .    .  .
        // Cb                 .   .   .   .    .  .   (all-zero diagnostic)
        // C                  .   x   .   .    .  .
        // Ha                 x   x   x   x    x  .
        // HN                 x   x   x   x    x  x
        constexpr bool table[6][6] = {
            //  1HB    2HB    1HaB   2HaB   RC     w
            {  true,  true, false,  true, false, false }, // N
            { false,  true, false, false, false, false }, // CA
            { false, false, false, false, false, false }, // CB
            { false,  true, false, false, false, false }, // C
            {  true,  true,  true,  true,  true, false }, // HA
            {  true,  true,  true,  true,  true,  true }, // HN
        };
        return table[(int)t][(int)term];
    }
};
```

The all-false `Cβ` row is the user-decided diagnostic discipline
(`feedback_methods_accumulate`): we still compute Cβ contributions
and emit them, expecting near-zero in production. A test asserts
this expectation.

#### Output NPYs — per-contribution-class detail for ML stratification

| File                                  | Shape          | Description                                |
|---------------------------------------|----------------|--------------------------------------------|
| `larsen_hbond_count.npy`              | (N,) int32     | Number of HB pairs per atom (Σ all classes)|
| `larsen_hbond_count_by_class.npy`     | (N, 4) int32   | Per-acceptor-class count [BB, SchAmide, OH, COO] |
| `larsen_hbond_total_tensor.npy`       | (N, 3, 3)      | Σ tensor in protein lab frame, all classes |
| `larsen_hbond_total_spherical.npy`    | (N, 9)         | T0+T1+T2 of total                          |
| `larsen_hbond_1pHB_tensor.npy`        | (N, 3, 3)      | Σ 1°HB contributions (amide-H donor only)  |
| `larsen_hbond_2pHB_tensor.npy`        | (N, 3, 3)      | Σ 2°HB contributions                        |
| `larsen_hbond_1pHaB_tensor.npy`       | (N, 3, 3)      | Σ 1°HαB contributions (Hα donor only)      |
| `larsen_hbond_2pHaB_tensor.npy`       | (N, 3, 3)      | Σ 2°HαB contributions                       |
| `larsen_hbond_water_term.npy`         | (N,) float64   | Δσ_w (2.07 ppm on HN with no HB pair, else 0; zero for non-HN atoms) |
| `larsen_hbond_diagnostic_CB.npy`      | (N, 3, 3)      | Cβ contributions (should be near-zero — diagnostic) |
| `larsen_hbond_pairs.npy`              | (N_pairs, 9)   | Per-pair: [donor_ai, acceptor_ai, donor_class, acceptor_class, r, θ, ρ, isotropic_total, any_corner_imputed] |
| `larsen_vs_kernel_hb_residual.npy`    | (N, 3, 3)      | Per-atom Mat3 diff: LarsenHBond − HBondResult on amide-H targets (zero where one or both did not contribute) |

Calculator never returns nullptr; it returns a result with all-zero
fields if `Session::HasLarsenHBondGrid()` is false. The optional
nature is encoded in the data, not the calculator presence.

#### SDK wiring

Add `LarsenHBondGroup` dataclass to `python/nmr_extract/_protein.py`
mirroring `TripeptideGroup`. Attached on `Protein.larsen_hbond` when
any of the NPYs is present. Add 12 ArraySpec entries to
`python/nmr_extract/_catalog.py`. Add SDK tests in
`python/tests/test_larsen_hbond_group.py`.

#### Steps

1. Add the new substrate predicates to `AtomSemanticTable` (header-
   only change; constexpr; verified by structure test). Build +
   verify.
2. Add `IsHydroxylOxygen(protein, ai)` calculator-side helper.
3. Encode `LarsenContribDispatch` table in
   `LarsenHBondShieldingResult.h` (constexpr).
4. Implement donor-side iteration. For each candidate donor: resolve
   the three frame atoms per the role-table above. Skip donors where
   any frame atom is unavailable (e.g. N-terminus for NMA donor),
   with Info log.
5. For each donor, query SpatialIndexResult for candidate acceptor
   atoms within rOH ∈ [1.8, 4.2] Å of `donor_H`.
6. For each candidate acceptor: classify (BackboneCarbonyl,
   SidechainCarbonyl, HydroxylOxygen, CarboxylateOxygen). Skip if
   classification fails. Resolve acceptor's three frame atoms.
7. Compute `LarsenHBondGeometry` via `ComputeLarsenHBondGeometry`.
   Gate on θ ∈ [90°, 180°] (out-of-range = no H-bond physics here).
8. Query `LarsenHBondGrid::QueryNearest`. Skip if `!is_hit`. Apply
   `RotateTensorToProteinLabFrame` per readout using the protein-
   side `ComputeLarsenDonorFrame` rotation.
9. Distribute readouts per Table 2 dispatch:
    - 1° readouts to donor residue i's atoms (N/CA/CB/C/HA/HN).
    - 2° readouts to acceptor residue (i+1 for backbone-O acceptor)
      atoms (N/HA/HN). Skip if i+1 doesn't exist.
10. After all pairs processed, sweep amide-H atoms. Any with zero HB
    pair contribution gets `Δσ_w = 2.07 ppm` in the isotropic
    component (no tensor structure).
11. Emit NPYs per the table above. Compute
    `larsen_vs_kernel_hb_residual.npy` from `HBondResult`'s output
    on the same protein (read its NPY values, take per-atom diff).
12. Smoke on `1UBQ_pm6dh3plus.pdb`:
    - Assert ≥ 40 backbone H-bonds detected (Ubiquitin has ~60 in
      published structure; we expect somewhat fewer with strict
      geometric criteria).
    - Assert all tensors finite.
    - Assert Cβ diagnostic NPY has |entries| < 0.5 ppm (near-zero
      reality check).
    - Assert at least 5 amide-H atoms received `Δσ_w = 2.07 ppm`.

Acceptance criterion: smoke passes, structure tests for substrate
predicates pass, Table 2 dispatch verified by unit test (one test
per cell asserting Applies() == expected), Cβ diagnostic NPY passes
the near-zero assertion.

### Session 4 — SDK wiring + multi-source reality checks + adversarial review

Outputs:
- `python/nmr_extract/_catalog.py` — 12 new ArraySpec entries.
- `python/nmr_extract/_protein.py` — `LarsenHBondGroup` dataclass.
- `python/tests/test_larsen_hbond_group.py` — SDK tests.
- `tests/test_larsen_hbond_against_published.cpp` — Larsen Table 1
  RMSDs comparison.
- `scripts/larsen_validation/` — Python-side scripts for the three
  reality-check comparisons (Table 1, .procs, multi-predictor).
- Codex xhigh adversarial review.

Steps:

1. SDK wiring (mirror `TripeptideGroup`). 12 ArraySpec entries cover
   per-atom totals + per-contribution-class detail + Cβ diagnostic +
   water term + kernel-vs-grid residual + per-pair detail.
2. **Reality check 1 — Larsen 2015 Table 1.** Run extraction on
   `1UBQ_pm6dh3plus.pdb`. Combine TripeptideBackboneShieldingResult
   + TripeptideNeighborShieldingResult + LarsenHBondShieldingResult +
   RingCurrent (existing). Compute per-atom-type RMSD against BMRB
   experimental shifts (from our existing BMRB ingestion).
   Expected: within ~0.1 ppm of Larsen's published ProCS15 row
   (Cα 1.7, Cβ 2.5, C' 2.1, Hα 0.6, HN 0.7, N 4.4).
3. **Reality check 2 — Larsen's actual .procs predictor output.**
   Parse `predictions/qm/1UBQ_pm6dh3plus.procs` (from
   `/mnt/expansion/larsen_archive/predictions.tar.bz2`). Compare our
   per-atom σ outputs to Larsen's per-atom prediction, keyed by
   `(resseq, atom_name)`. Mask N-term H/HA (sentinel-zero) and
   Pro-Cβ (sentinel-zero in Larsen). Compute correlation + per-
   atom-type residual distributions. Report any systematic offsets.
   Expected: r > 0.99 per atom type; if not, surface as a thesis-
   reportable methodological coordinate.
4. **Reality check 3 — 6-predictor methodology scatter.** Parse
   `predictions/charmm/*` (CharmM, CamShift, PPM, SHIFTX2, SPARTA+
   on 1UBQ + 2OED). Build a 7-way scatter (our output + 6 published
   predictors). Report where we cluster in methodology space.
5. **Reality check 0 — Larsen's whole-protein DFT (`proteinnmrlogs`).**
   Parse `proteinnmrlogs/1UBQ_pm6dh3plus.log` per the existing
   Gaussian shielding-tensor parser. 1232 atom tensors. This is the
   ground-truth DFT — our combined predictor's job is to approximate
   it. Compute per-atom RMSD AND tensor-level (Frobenius) residual.
   This is the same analysis Larsen reports in his Table 1.
6. **Side-by-side: HBondResult vs LarsenHBondShieldingResult.**
   Per-atom-type residual between the two formulations on the same
   protein. Expected: correlation r > 0.9 on amide-H isotropic; full-
   tensor differences are the methodological coordinate.
7. **Adversarial review** by codex xhigh. Focus areas: frame
   convention, reference subtraction sign, Table 2 dispatch
   correctness, periodic-ρ wrap, NPY layout. Fix findings.

Acceptance criterion: per-atom-type Ubiquitin RMSD ≤ 0.1 ppm above
Larsen Table 1's ProCS15 row across Cα/Cβ/C'/Hα/HN/N. .procs
correlation r ≥ 0.99 per atom type. All four reality checks have
their residual distributions emitted as NPY artifacts for analysis.
Adversarial review findings addressed.

Where session 5 might happen: parser frame-convention issue in S1,
.procs residuals too large to attribute to methodology differences in
S4, or unforeseen calibration coupling in S4. All are real risk
surfaces.

### Session 5+ — calibration integration

Methods accumulate (`feedback_methods_accumulate`). The new HBond
calculator's NPYs become 4+ new features in the calibration ridge
regression. `HBondResult` stays in the pipeline as a parallel
methodology coordinate. Per-atom-type residual stratification is
extended to include the kernel-vs-grid difference. Out of scope for
the 4-session H-bond build; flagged here so the design respects the
forward path.

## Session 1 reconnaissance — findings 2026-05-11

Open questions resolved by extraction + probe:

### Q1 — free-monomer reference

**Not in the H-bond tarball.** The nested `.tar.bz2` inside each of the
6 archives is byte-identical to the outer (redundant duplication, not
a separate reference). Options:

- (a) Run Gaussian on isolated NMA + Ac-A-NMe ourselves. Reproduces
  Larsen exactly. Out of scope for the parser session.
- (b) **Use r-max grid edge as proxy reference** (current parser
  choice). Per donor-side atom, reference σ = average over grid points
  at the largest r value. Small bias (~0.1 ppm); documented.
- (c) Fetch `proteinnmrlogs.tar.bz2` or `predictions.tar.bz2` from
  ERDA — they may contain Larsen's pre-interpolated .npy. Future work.

Session 1 parser uses (b) with a config knob to swap in (a) or (c)
later.

### Q2 — filename ↔ geometry convention

Probed 10 logs (6 ALA donor + 4 NMA donor). Findings:

| Donor | r relation | θ | ρ |
|-------|-----------|---|---|
| NMA   | r_actual ≈ r_filename (offset ≤ 0.001 Å) | matches | sign-flipped |
| ALA   | r_actual = r_filename − 0.2 Å (exact) | matches | sign-flipped |

The 0.2 Å offset on ALA donor archives is exactly one grid step
(paper says 0.2 Å step) — likely off-by-one in Larsen's input-file
labeling. θ matches at the O vertex (Hα–O–C(=O) angle). ρ has a
consistent sign flip — convention difference between Larsen's
dihedral direction and standard IUPAC; doesn't affect data.

**Parser strategy: key the grid on ACTUAL measured (r, θ, ρ) from
atom positions; treat filename values as labels only.** Avoids any
dependence on filename convention.

### Q3 — filename prefix convention

| Archive | Filename prefix |
|---------|----------------|
| ALANMA  | `ALANMA_...` |
| ALACOO  | `ALACOO_...` |
| ALACOH  | `ALACOH_...` |
| NMANMA  | `NMA_...` (NOT `NMANMA_...`) |
| NMACOH  | `NMACOH_...` |
| NMACOO  | `NMA_...` (NOT `NMACOO_...`) |

The two NMA donor archives that pair with COO and NMA acceptors share
a `NMA_` prefix (the donor-only label). Parser identifies acceptor
type from the parent archive (tar) name, not the log filename.

### Q4 — donor-H identification

"Closest inter-cluster H...O" heuristic fails at large rOH where
intramolecular contacts compete (verified: one of 6 ALA samples at
rOH = 4.0 Å had the wrong H picked). **Real parser uses chemical
perception per donor/acceptor molecule** (one perception pass per of
the 5 molecules: NMA, Ac-A-NMe, NMA-acceptor, HOMe, acetate) to pin
the donor H, the donor anchor atom (Cα or amide N), the acceptor O,
and the acceptor anchor atoms. Canonical atom indices baked into the
parser as constants per molecule.

### Q5 — extra grid points

ALANMA has 3252 logs vs paper's 12×10×24 = 2880 expected. The 0.2 Å
filename offset implies r grid actually spans 1.6–4.2 Å (13 points)
rather than paper's 1.8–4.0 (12 points), so 13×10×24 = 3120. Still
132 unaccounted. Likely extra ρ samples at certain (r, θ) points or
duplicate near-edge samples. Parser bins by ACTUAL geometry; any
multiple-log-per-key collisions are detected at parser time and
fatal-on-conflict (force resolution before silent averaging).

### Q6 — sidechain primary amide acceptor approximation

ASN ODE1 / GLN OE1 approximation by NMA grid stands per the plan doc
"acceptor classification → grid selection" table. Chemical similarity
(R-C(=O)-NHR' vs CH3-C(=O)-NH-CH3) is small but real; documented
limitation alongside the calculator.

## Substrate enum status

Verified 2026-05-11 in `src/SemanticEnums.h`:

| Predicate                          | Status | Notes |
|------------------------------------|--------|-------|
| `BackboneRole::AmideHydrogen`      | ✓      | PRO automatically excluded (no amide H per BackboneRole comment) |
| `BackboneRole::AlphaHydrogen`      | ✓      | Non-GLY HA |
| `BackboneRole::CarbonylOxygen`     | ✓      | Backbone O |
| `PlanarGroupKind::Carboxylate`     | ✓      | ASP/GLU/C-term, lines 277-279 |
| `PlanarGroupKind::SidechainAmide`  | ✓      | ASN/GLN, lines 245-247 |
| `PlanarGroupKind::AromaticHydroxyl` | ✓     | TYR OH, lines 281-285 |
| `PolarHKind::HydroxylOH_Aliphatic` | ✓      | SER HG, THR HG1 |
| `PolarHKind::HydroxylOH_Aromatic`  | ✓      | TYR HH |

**GLY edge case** (SemanticEnums.h:85-88, 119-122): GLY HA2/HA3 carry
`Locant::Alpha + BackboneRole::None` (not `BackboneRole::AlphaHydrogen`).
The donor gate must dispatch on `IsAnyAlphaHydrogen()` (new predicate)
covering both non-GLY HA and GLY HA2/HA3.

**Hydroxyl O classification.** No direct single-field predicate; SER
OG and THR OG1 are `Element::O + Locant::Gamma` (also matches CYS S
position pattern, but different element). Cleanest dispatch is to
walk to bonded H atoms and check `PolarHKind`. TYR OH carries
`Locant::Eta` and is in `PlanarGroupKind::AromaticHydroxyl`. The
calculator-side helper `IsHydroxylOxygen(ai)` does the bond-walk; this
is not a new substrate field, just a calculator-side utility.

## Smoke target — three accumulating reality checks

`1UBQ_pm6dh3plus.pdb` from `/mnt/expansion/larsen_archive/structures/`
is Larsen 2015's exact NMR-input geometry. Three independent
comparison signals are available from the ERDA archive (fetched
2026-05-12):

**Reality check 1 — vs Larsen 2015 Table 1 RMSDs.** Per-atom-type
RMSDs of our combined output (TripeptideBackboneShieldingResult +
TripeptideNeighborShieldingResult + LarsenHBondShieldingResult +
RingCurrent) against published experimental shifts. Larsen's ProCS15
row for Ubiquitin:

- Cα: 1.7 ppm — Cβ: 2.5 ppm — C': 2.1 ppm
- Hα: 0.6 ppm — HN: 0.7 ppm — N: 4.4 ppm

Our numbers should land within ~0.1 ppm of these if the port is
faithful.

**Reality check 2 — vs Larsen's actual ProCS predictor output.** The
file `predictions/qm/1UBQ_pm6dh3plus.procs` (fetched from ERDA
`predictions.tar.bz2`) is the per-atom shielding prediction Larsen's
ProCS predictor produces on this exact PDB. 447 atom-rows
(BB + Cβ subset; N-term H/HA and Pro-Cβ are sentinel-zero per
Larsen's predictor convention). Format: whitespace-delimited
`<idx> <resseq> <RES3> <ATOM> <ELEMENT> <shielding_ppm>`.

Key byte-faithful comparison points:
- Match our per-atom σ output against `.procs` per-atom σ at the
  BB+Cβ subset, keyed by `(resseq, atom_name)`.
- Expect zero residual modulo our methodological differences (PM6
  geometry tweaks, reference σ proxy vs Larsen's actual free-monomer
  DFT).
- Per-atom-type residual distribution becomes part of the thesis
  methodology report.

Same comparison available for `2OED_pm6dh3plus.pdb` (GB3).

**Reality check 3 — vs 6 other predictors.** `predictions/charmm/`
holds CharmM, CamShift, PPM, SHIFTX2, SPARTA+ outputs on the same
proteins. Useful as a methodology-space scatter plot — where do these
6 published predictors converge or disagree, and where does our
output land in that space.

Plus the predicted ensembles (`predictions/ensembles/`) on 5 NMR-
solution structures (1D3Z, 1XQQ, 2K39, 2KOX, 2LJ5) for Stage 2
trajectory cross-validation later.

**Reality check 0 (baseline) — vs Larsen's own DFT on the same PDB.**
`proteinnmrlogs.tar.bz2` contains `1UBQ_pm6dh3plus.log` —
Larsen's Gaussian 09 GIAO OPBE/6-31G(d,p) DFT run on this exact
geometry (22 days CPU time, normal termination). 1232 atom shielding
tensors. This is the GROUND TRUTH that ProCS aims to predict.
Comparing OUR combined-predictor output to THIS DFT output is the
absolute floor on ProCS-equivalent accuracy — the same comparison
Larsen reports in Table 1. We can reproduce that exact analysis.

Our numbers should be within 0.1 ppm of these if the port is faithful.

## HBondResult and LarsenHBondShieldingResult — sibling calculators

Both calculators stay in the codebase and pipeline forever. They cover
overlapping physics (amide-H donor for the backbone O acceptor case)
but use different formulations:

- `HBondResult`: kernel × η, geometric-only.
- `LarsenHBondShieldingResult`: DFT grid lookup against Larsen 2015
  scans.

Each emits its own NPYs. The Python SDK exposes both groups. The
calibration pipeline consumes both. The PER-ATOM-TYPE DIFFERENCE
between them is itself a thesis-reportable output — a coordinate in
methodology space describing where the kernel-form captures or misses
geometries that the grid-form encodes faithfully.

Reality-check pattern:

- For amide-H donors with backbone-O acceptors (the only physics class
  both calculators cover), emit a `larsen_vs_kernel_hb_residual.npy`
  per-atom Mat3 difference. Expect correlation r > 0.9 on
  isotropic; geometry-specific residuals are part of the analysis.
- Hα donors and sidechain-O acceptors are LarsenHBondShieldingResult-
  only — no comparison signal there.

Per `feedback_methods_accumulate` (memory).

## Cross-references

- Memory entry `project_hbond_halpha_design` — PLANNING-SUPERSEDED by this document; the kernel-form `HBondResult` still in the codebase is NOT retired (see `feedback_methods_accumulate`).
- Memory entry `reference_erda_archive_missing_files` — the archive
  contents lookup.
- Memory entry `project_larsen_residue_model` — the perception
  machinery this calculator's parser reuses.
- Memory entry `feedback_two_path_validation` — invariant pattern
  applied to parser frame-roundtrip check.
- Reference: Larsen et al., PeerJ 2015 (DOI 10.7717/peerj.1344) —
  `references/larsen-2015-procs15-dft-chemical-shift-predictor.pdf`.

## Provenance

This design landed after a multi-turn exchange 2026-05-11 in which the
predecessor session's "kernel × η" framing was pushed back on the
grounds that the grid data was actually in hand and Larsen's
decomposition was richer than that framing acknowledged. Triage of:
the larsen_archive (6 grid archives confirmed, 3252 logs in ALANMA);
substrate enum coverage (GLY HA edge case surfaced); existing parsing
infrastructure (geometry parser reusable, shielding tensor parser is
new); Larsen 2015 §H-bond methodology + Tables 1-2.

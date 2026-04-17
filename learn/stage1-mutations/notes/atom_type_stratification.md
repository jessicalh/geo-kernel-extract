# Atom-Type Stratification

2026-04-15.  720 proteins, 446K atoms.  Per-atom-type ridge using
AMBER atom names from the prmtop.

Script: src/actual_physics/atom_type_calibration.py
Output: src/output/actual_physics/atom_type_calibration/

---

## The problem

All prior Stage 1 analysis pooled atoms by element only (H, C, N, O).
This hides massive variation between atom types within an element.
A backbone carbonyl carbon and a sidechain methyl carbon have
fundamentally different electronic environments and respond to
aromatic ring removal through different mechanisms.

For NMR work or any biology-facing presentation, pooling all carbons
or all nitrogens into one number is not meaningful.  The correction
is straightforward: classify by AMBER atom name from the prmtop.

---

## Atom type definitions

Carbon:
- **CA** — alpha carbon (sp3, backbone)
- **C=O** — carbonyl carbon (sp2, backbone, bonded to O)
- **CB** — beta carbon (sp3, first sidechain position)
- **C side** — all deeper sidechain carbons (CG, CD, CE, CZ, ...)

Nitrogen:
- **N bb** — backbone amide nitrogen
- **N side** — sidechain nitrogen (ND, NE, NZ, NH in HIS, ARG, LYS, ...)

Oxygen:
- **O bb** — backbone carbonyl oxygen
- **O side** — sidechain oxygen (OG, OH, OD, OE in SER, THR, TYR, ASP, GLU, ...)

Hydrogen kept as one group (R²=0.921 raw, no sub-type needed).

---

## Results: raw kernels (55 core, no normalisation)

| Atom type | n | R² |
|-----------|---------|------|
| H | 230,135 | 0.921 |
| C all | 133,488 | 0.514 |
|   CA | 29,944 | 0.541 |
|   C=O | 29,944 | 0.361 |
|   CB | 27,429 | 0.578 |
|   C side | 46,171 | **0.660** |
| N all | 39,954 | 0.210 |
|   N bb | 29,944 | 0.212 |
|   N side | 10,010 | **0.588** |
| O all | 42,429 | 0.231 |
|   O bb | 29,944 | 0.253 |
|   O side | 12,485 | 0.241 |

---

## Results: normalised + progressive scalars

| Atom type | n | Base | +Scales | +Mut | Fair |
|-----------|---------|-------|---------|------|------|
| H | 230,135 | 0.848 | 0.924 | 0.861 | 0.928 |
| C all | 133,488 | 0.471 | 0.529 | 0.512 | 0.562 |
|   CA | 29,944 | 0.539 | 0.577 | 0.597 | 0.627 |
|   C=O | 29,944 | 0.372 | 0.411 | 0.430 | 0.463 |
|   CB | 27,429 | 0.544 | 0.597 | 0.604 | 0.647 |
|   C side | 46,171 | 0.646 | 0.690 | 0.694 | **0.729** |
| N all | 39,954 | 0.245 | 0.292 | 0.345 | 0.380 |
|   N bb | 29,944 | 0.254 | 0.301 | 0.351 | 0.387 |
|   N side | 10,010 | 0.720 | 0.762 | 0.870 | **0.887** |
| O all | 42,429 | 0.274 | 0.303 | 0.358 | 0.382 |
|   O bb | 29,944 | 0.303 | 0.338 | 0.395 | 0.422 |
|   O side | 12,485 | 0.334 | 0.373 | 0.543 | **0.566** |

---

## Key findings

### 1. "Nitrogen is hard" was wrong

Pooled N R²=0.210 (raw).  But sidechain nitrogen is 0.588 raw,
0.887 fair — the second-best atom type after hydrogen.  Backbone
amide N (0.212 raw, 0.387 fair) is hard.  The pooled number was
dominated by 30K backbone atoms hiding 10K well-predicted sidechain
atoms.

### 2. Carbonyl carbon drags the carbon average down

C=O at 0.361 raw (0.463 fair) vs sidechain carbon at 0.660 raw
(0.729 fair).  The carbonyl's sp2 electronic structure with C=O pi*
excitation creates a fundamentally different shielding mechanism
(paramagnetic dominant) that the geometric kernels capture poorly.
Pooling it with CA, CB, and sidechain carbons diluted the signal.

### 3. Mutation type is atom-type-dependent

The +Mut column shows which atom types are sensitive to ring
identity (PHE vs TYR vs TRP vs HIE):

| Atom type | Base → +Mut | Delta |
|-----------|------------|-------|
| H | 0.848 → 0.861 | +0.013 |
| C side | 0.646 → 0.694 | +0.048 |
| N bb | 0.254 → 0.351 | +0.097 |
| N side | 0.720 → 0.870 | **+0.150** |
| O bb | 0.303 → 0.395 | +0.091 |
| O side | 0.334 → 0.543 | **+0.209** |

Sidechain N and O are the most sensitive to ring identity.  These
are the atoms most likely in direct contact with the aromatic ring
(histidine stacking, serine/threonine H-bonds).  Different ring
types have different charge distributions and current intensities.
Ring identity matters most where contact is closest.

Hydrogen is insensitive (+0.013) — ring current geometry alone is
sufficient regardless of ring type.

### 4. Sidechain atoms are uniformly better predicted

Every element shows sidechain > backbone.  Sidechain atoms are
closer to the removed aromatic ring and have more diverse geometric
relationships to it, giving the kernels more angular signal to work
with.

---

## Implications for the thesis

The element-pooled numbers (H=0.928, C=0.562, N=0.380, O=0.382)
remain correct as weighted averages but should never be presented
without the atom-type decomposition.  The decomposition reveals:

1. The tool is strong on sidechain atoms (0.57-0.89 fair).
2. The tool is weak on backbone atoms (0.39-0.46 fair for C=O, N bb).
3. Backbone weakness is expected: backbone geometry doesn't change
   in mechanical mutants, so backbone atoms see the perturbation
   only through long-range fields.
4. The mutation type sensitivity of sidechain N/O is Case (1995)'s
   ring-type-specific intensity factors at work in the T2 channel.

---

## Implications for Stage 2

The nonlinear signal (completeness_and_nonlinear.md) should be
re-examined per atom type.  The pooled N RF delta of +0.169 may
be concentrated in backbone N (where linear ridge is weak) while
sidechain N (already at 0.887 with ridge) may have little nonlinear
headroom.

The charge-polarisation gap for carbon (+0.197 pooled) should be
decomposed: is it concentrated in C=O (where paramagnetic mechanisms
dominate) or distributed across all carbon types?

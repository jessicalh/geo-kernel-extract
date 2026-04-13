# Dia/Para Decomposition of the Mutation Delta

2026-04-13.  Validated against C++ to < 0.003 ppm on 149 proteins.

---

## Method

Parse orca diamagnetic and paramagnetic 3x3 tensors from raw
{ID}_WT_nmr.out and {ID}_ALA_nmr.out files.  Match WT→ALA atoms
by element + 3D position (0.1 Å threshold, mechanical mutants have
identical backbone positions).  Compute delta_dia = WT_dia - ALA_dia
and delta_para = WT_para - ALA_para per matched atom.

Script: `src/actual_physics/orca_dia_para.py`

Validation: delta_dia_T2 + delta_para_T2 reproduces the C++ stored
delta_shielding T2 to max 0.003 ppm across all elements.  The
residual is floating-point rounding from orca output formatting.

## Key finding: dia and para cancel

| Element | n atoms | |delta_dia| | |delta_para| | |delta_total| | para/dia var |
|---------|---------|-----------|-------------|--------------|-------------|
| H | 46,794 | 6.44 | 6.49 | 0.78 | 1.02 |
| C | 27,953 | 7.42 | 7.67 | 1.54 | 1.06 |
| N | 8,185 | 6.96 | 7.27 | 1.58 | 1.07 |
| O | 9,265 | 6.55 | 7.15 | 2.13 | 1.20 |

The mutation perturbs both diamagnetic and paramagnetic shielding by
~7 ppm each.  They nearly cancel to ~1-2 ppm total delta.  For H the
cancellation is almost exact (0.78 from 6.5+6.5).  Oxygen shows the
most paramagnetic excess (1.20 ratio).

The para/dia variance ratio increases H→C→N→O, confirming Saito
et al. 2010 (Prog NMR Spectrosc 57:181) in the T2 delta channel.

## Kernels see the net, not the channels

| Element | BS→delta_dia | BS→delta_para | BS→delta_total |
|---------|-------------|--------------|---------------|
| H | 0.002 | 0.023 | 0.704 |
| C | 0.003 | 0.019 | 0.075 |
| N | 0.012 | 0.030 | 0.078 |
| O | 0.005 | 0.051 | 0.132 |

R² for BS or EFG against individual delta_dia or delta_para is
< 0.05 for every element.  But against delta_total: BS gets 0.70
for H, EFG gets 0.73.

The geometric kernels compute classical fields corresponding to the
total interaction, not the QM decomposition into dia/para terms.
This is physically correct — ring current is a net magnetic effect,
not separable into diamagnetic and paramagnetic contributions at the
classical level.

## Cosine alignment

BS→delta_dia cosine: ~0.51 for all elements.
BS→delta_total cosine: 0.89 for H, 0.73 for O.

The individual dia/para delta T2 vectors are poorly aligned with any
geometric kernel (~0.51 is between random 0.45 and collinear 1.0).
The total delta is well-aligned.  The cancellation that produces the
total from dia+para also produces the alignment.

## What this means for the thesis

1. The geometric kernels are correctly calibrated against the total
   shielding perturbation, which is the physically meaningful quantity.

2. The dia/para decomposition is a DFT internal bookkeeping that does
   not map to the classical kernel decomposition.  This is expected
   and not a limitation.

3. The para/dia variance ratio provides independent confirmation of
   element-dependent shielding physics (Saito 2010), now measured in
   the T2 channel of mutation deltas.

4. The near-perfect cancellation (7+7→1 ppm) means the mutation
   perturbation lives in the small residual after a large cancellation.
   This makes the measurement sensitive and explains why per-protein
   normalization is necessary — the signal is a small angular residual
   on top of a large but nearly self-cancelling perturbation.

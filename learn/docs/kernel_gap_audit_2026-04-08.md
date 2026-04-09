# Kernel Set Gap Audit — 2026-04-08

Cross-reference of the 10 output calculators (spec) against the 64-kernel
calibration pipeline (learn/src/mutation_set/kernels.py + scalars.py).

## 1. Pi-Quadrupole T2 entirely absent from kernel set

Calculator 5 computes the EFG from the ring pi-electron quadrupole moment.
The SDK exposes both `p.pi_quadrupole.per_type_T2.as_block()` (8 ring
types x 5 components) and `p.pi_quadrupole.shielding.T2` (total).  It is
pure T2 (symmetric traceless; T0=0, T1=0 from the physics).

Its distance dependence is **1/r^5**, not 1/r^3.  Every other T2 kernel
in the set is either dipolar (1/r^3) or short-range switching (dispersion,
5 A cutoff).  The pi-quadrupole fills the intermediate range: it dominates
at close approach and falls off faster than dipolar at distance.  The model
currently has **zero kernels with 1/r^5 spatial character**.

The pipeline _knows_ pi_quadrupole exists: `scalars.py:111` reads
`p.pi_quadrupole.shielding.T1` and takes its norm.  But pi-quadrupole is
pure T2 — the T1 magnitude being captured is numerical noise.  The actual
T2 tensor is discarded.

**Missing from kernels.py:**
- 8 per-type T2 kernels (same interface as BS/HM/Disp)
- 1 total T2 kernel (same interface as other calculator totals)
- Per-ring PQ T2 via `ring_contributions.pq` for the residual mechanism

**Missing from scalars.py:**
- Pi-quadrupole scalar (Buckingham A-term: `(3 cos^2 theta - 1)/r^4`)

## 2. Haigh-Mallion raw H tensor discarded

The HM calculator produces two distinct tensors:
- **H** — the raw surface integral, symmetric traceless, pure T2
- **G = n (x) (H . n)** — the shielding kernel, asymmetric, T0+T1+T2

The pipeline uses G.T2 via `per_type_T2`.  But H.T2 != G.T2 — they have
different angular structure.  H is the geometric field of the surface
integral; G filters it through the ring normal projection.

H is available per-ring via `ring_contributions.hm_H` (cols 18-26).
If the ring normal from force-field geometry is slightly wrong vs QM,
H is the more robust quantity.

## 3. Per-ring residuals incomplete

The per-ring residual mechanism (kernels 46-63) applies to BS, HM, Chi.
Two ring calculators are missing:

- **Dispersion per-ring**: the per-type sum is in the kernel set (cols
  16-23) but individual ring dispersion is in `ring_contributions` and
  not used for residuals.
- **Pi-quadrupole per-ring**: available as `ring_contributions.pq`,
  not used at all.

For multi-ring proteins (Trp has 3 ring definitions), the residual
captures angular structure the type sum collapses.  This mechanism is
only applied to 3 of the 5 ring calculators.

## 4. Ring susceptibility per-type — confirmed not available

The SDK exposes `ring_susceptibility` as a bare `ShieldingTensor` with no
`per_type_T2`.  The pipeline correctly uses the total.  Not a gap.

## Scale

| Missing kernel group          | Count | Physics                                    |
|-------------------------------|-------|--------------------------------------------|
| PQ per-type T2                | +8    | 1/r^5 angular structure per ring type       |
| PQ total T2                   | +1    | 1/r^5 total                                |
| PQ per-ring residuals (K=6)   | +6    | Individual ring quadrupole residuals        |
| HM H per-ring (K=6)          | +6    | Raw surface integral, pure T2              |
| Disp per-ring residuals (K=6) | +6    | Individual ring dispersion residuals        |
| **Total**                     | **+27** | Nearly half the current 64-kernel count   |

The pi-quadrupole alone (15 kernels) is the worst omission: an entirely
different distance dependence the model cannot currently represent.
Anything it explains is forced into dipolar kernels or left as residual.

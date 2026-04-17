# EFG Signal Is Geometric, Not Charge-Dependent

## The finding

Three independent physical models get the same R^2 for hydrogen
shielding T2 in the mechanical mutant calibration:

| Model | R^2 | Mechanism |
|---|---|---|
| BS (8 per-type) | 0.772 | Magnetic dipole (ring current) |
| HM (8 per-type) | 0.784 | Surface integral (ring current) |
| MopacEFG_aro | 0.766 | Electric field gradient (MOPAC charges) |
| efg group (all) | 0.770 | All EFG variants |

Angular correlation with DFT target T2 direction:

| Kernel | cos_H |
|---|---|
| EFG_aro (topology charges) | 0.9296 |
| MopacEFG_aro (MOPAC charges) | 0.9323 |

Difference: 0.003. The DIRECTION of the EFG tensor is determined
by the geometric T-tensor (3 r_a r_b - r^2 delta_ab)/r^5, not by
the charge magnitudes. MOPAC vs topology charges change the
magnitude, not the angular structure.

## Why the old analysis was misleading

Forward selection is greedy. MopacEFG_aro entered first at 0.766.
Topology EFG_aro entered at position 18, adding +0.0003 — because
MopacEFG_aro already captured the angular signal. This made MOPAC
look essential. It isn't.

If forward selection were run WITHOUT any MOPAC kernels, topology
EFG_aro would enter first and score ~0.76 for H.

## What this means for the ensemble analysis

The geometry-only extraction computes topology Coulomb EFG. This
carries the same directional information as MOPAC EFG (cos 0.93).
We do not need MOPAC charges for the ensemble investigation.

The R^2 = 0.77 from EFG is really R^2 = 0.77 from RING GEOMETRY
seen through the electric field mechanism. BS, HM, and EFG are
three views of the same geometric relationship between ring atoms
and field atoms. They agree because the geometry is what matters.

## SPECULATIVE implication

This means the entire 0.93 R^2 for H (ring_current group) is
geometric. The charges contribute <1% additional angular information.
For the ensemble analysis, every kernel we have is as good as we
need it to be. No charge calculation needed.

The interesting question becomes: in the ensemble, when the
GEOMETRY changes across frames, all three views (BS, HM, EFG)
should move together. Measuring their co-movement tests whether
the geometric interpretation holds under conformational variation.
When they DIVERGE (different frames where BS and EFG disagree),
that flags a geometry where the three models aren't equivalent —
which is a near-field non-dipolar effect worth investigating.

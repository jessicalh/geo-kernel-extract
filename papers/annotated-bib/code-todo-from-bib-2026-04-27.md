# Code-side notes surfaced by the bib reading

These are implementation observations that came up while writing the annotated bibliography. They live here in scholarship territory until a code session moves them to `spec/` (likely `spec/PLANNED_CALCULATORS_TIME_SERIES_2026-04-24.md` or a new spec doc). Listed in the order they surfaced.

## 2026-04-27 — Ring-normal-stable BS/HM kernel (Sahakyan & Vendruscolo 2013)

The current `BiotSavartResult` and `HaighMallionResult` calculators compute the ring's normal vector from three fixed ring atoms. Sahakyan & Vendruscolo 2013 (J. Phys. Chem. B 117, 1989–1998) flag this as a real instability for MD trajectories:

> "...the ring normal computed from three ring atoms is unstable to out-of-plane fluctuations of the atomic positions in molecular dynamics..."

They propose averaging two independently-computed normals (each from a different three-atom subset of the ring) as a fix. For our Stage-2 trajectory work, where every frame is fed through the kernels, this stability issue could bias T2 predictions in regions where ring atoms wobble out of plane.

**Concrete to-do**: a ring-normal-stable variant of the BS and HM kernels — either (a) two-subset averaging per their suggestion, or (b) least-squares plane fitting through all ring atoms, or (c) test both against full-DFT shielding on out-of-plane-fluctuating ring frames to see which is closer.

**Stage**: Stage-2 implementation; surfaces as a kernel-level fix to the existing `BiotSavartResult` / `HaighMallionResult` calculators rather than a new calculator family.

**Reference for the spec doc**: Sahakyan & Vendruscolo 2013, *J. Phys. Chem. B* 117(7), 1989–1998. doi:10.1021/jp3057306. PDF in `references/sahakyan-vendruscolo-2013-ring-current-electric-field-contributions.pdf`.

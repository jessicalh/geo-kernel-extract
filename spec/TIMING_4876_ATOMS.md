# Calculator Timing: 1Q8K_10023 (4876 atoms, 300 residues)

Measured 2026-04-12 on batcave. CHARMM36m charges from TPR.
No MOPAC. 6 fleet poses (boltzmann minimum + 5 dihedral deviations).

## Per-calculator (pose 1, representative)

| Calculator | Time (ms) | Notes |
|---|---|---|
| GeometryResult | 8 | |
| SpatialIndexResult | 1722 | 2.99M neighbour pairs, 15A cutoff |
| EnrichmentResult | 1 | |
| DsspResult | ~180 | (gap between Enrichment end and ChargeAssignment start) |
| ChargeAssignmentResult | 0 | CHARMM36m preloaded from TPR |
| BiotSavartResult | 613 | 19057 atom-ring pairs, 26 rings |
| HaighMallionResult | 1344 | 19057 atom-ring pairs |
| **McConnellResult** | **5418** | **1.17M atom-bond pairs, 4887 bonds** |
| RingSusceptibilityResult | 91 | |
| PiQuadrupoleResult | 89 | |
| DispersionResult | 104 | |
| **CoulombResult** | **24781** | **O(N*k), 20A cutoff, dominates pipeline** |
| **HBondResult** | **7129** | **Unexpectedly expensive at this N** |
| SasaResult | 511 | Shrake-Rupley |
| **Total** | **41997** | |

## With APBS (pose 1)

| Calculator | Time (ms) |
|---|---|
| ApbsFieldResult | 4587 |
| CoulombResult | 24492 |
| **Total** | **46220** |

APBS adds 4.2s (10%) to a 42s pipeline. Fixed 161^3 grid — cost is
nearly independent of atom count.

## Across all 6 poses (with APBS, no MOPAC)

| Pose | APBS | Coulomb | McConnell | HBond | Total |
|---|---|---|---|---|---|
| 1 | 4587 | 24492 | 5332 | 7065 | 46220 |
| 2 | 4681 | 25629 | 5697 | 7266 | 48161 |
| 3 | 3898 | 25576 | 5577 | 7268 | 47139 |
| 4 | 3168 | 22829 | 5121 | 7061 | 42736 |
| 5 | 4207 | 25502 | 5459 | 6945 | 46810 |
| 6 | 3638 | 26124 | 5479 | 7296 | 47294 |

## Decision: APBS replaces Coulomb

- CoulombResult: 25s, O(N*k), vacuum partial charges, home-rolled
- ApbsFieldResult: 4s, O(1) grid solve, solvated Poisson-Boltzmann, well-cited
- APBS is 6x faster AND better physics on this molecule
- Coulomb retained for special comparison tests only

## Cost projections

At 46s/pose (with APBS, no Coulomb, no MOPAC): ~17s/pose.
At 17s/pose × 6 poses = ~100s per protein.
At 17s/pose × 685 proteins × 6 poses = ~20 hours on one machine.
With 3 machines and 40 days: easily feasible even with more frames.

## Scaling suspects

- **McConnell** (5.4s): 1.17M atom-bond pairs, effectively O(N*B)
- **HBond** (7.1s): needs investigation — may have O(N^2) loop
- **SpatialIndex** (1.7s): KD-tree build, scales well

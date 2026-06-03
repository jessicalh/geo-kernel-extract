# Partition Filter Design

Status: durable summary of the Build1/Build1Fix architecture (`a3531040655955a4d393e314c5a69de51b1c54bf`, fixed by `d9ba53dc0a7929eb4d4ce8185f6d9e6ff868c374`).

## Purpose

The partition filter is the C++ input-side layer that lets clean mechanism examples emerge from the all-atom substrate without selecting on DFT targets, residuals, or fitted coefficients.

It exists to separate two jobs:

- **Law layer:** dominance-isolated, mechanism-clean exemplars used to recover an expected form and coefficient.
- **Model layer:** mixed all-atom/per-type fits that explain variance but may shrink, pool, or interact terms.

Dominance is the gate between them: use it to decide whether a row/window is clean enough for a law read; do not let the coarse model wash decide the law.

## C++ Boundary

Python receives emitted scalars, bins, sidecars, and manifests; it fits, slices, plots, and audits only.

C++ owns source association, geometry, exclusions, reductions, isolation, bin ids, and case selection because those are physics/model-state reads over the resident trajectory.

The analysis script must not open `trajectory.h5`, ORCA outputs, older run dirs, or per-source dumps to reconstruct these values.

## Resident Indexes

`ResidentIndexes` pairs `RunData` with day-one immutable indexes:

- `TypedAtomIndex` for chemically typed scopes and habitats.
- `SpatialIndexSet` for per-frame source clouds: atoms, bond midpoints, ring centers, charge sites, and all-bond midpoints.
- `RingGeometryCache` for stable ring centers/normals.
- `TemporalIndex` for frame windows.

These indexes keep partitioning as queries over resident model state, not as a Python join/rebuild problem.

## Isolation Primitives

Build1 emits compact per-(atom, frame) scalars:

- `nearest_*_r` and `gap_to_2nd_{ring,charge,bond}_r` for source separation.
- `dominant_fraction_{ring,charge,mc}` for whether one mechanism dominates the local sum.
- mechanism magnitudes and modulation-by-atom features for driver exercise.
- partition bin ids for distance, gap, magnitude, and modulation axes.

Build1Fix made the charge path use the same self/same-residue exclusion as the charge contribution source set and made DFT selection reads throw during `CaseHunter`.

## CaseHunter

`CaseHunter` searches typed habitats over frame windows using only input-side predicates: isolation, driver motion, quiet competing mechanisms, and support.

It reports candidate strips plus `dft_recovery_R2` only after selection; that metric is a measurement label, not a selection input.

The anti-circular guard is real: any DFT target read during selection raises an error.

## Current Result

The live substrate is `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build1`.

Build2 showed the first partition result but exposed pooled-fit anti-prediction for H strata; Build3 fixed the fit architecture by adding per-type fits and dominance response curves.

Next architecture step, if re-emitting, is to move dominance quantile bin ids into C++ beside the other partition bins.

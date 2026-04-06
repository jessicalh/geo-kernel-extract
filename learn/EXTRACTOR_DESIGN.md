# Parameterised Geometry Extractor — Design

## TOML Parameter Schema

```toml
# =============================================================
# Calculator parameters
# =============================================================
# Default values match current hardcoded constants.
# Override any value to run a sweep experiment.

[calculators]

  [calculators.biot_savart]
  cutoff = 15.0                  # Å — ring search radius
  # Per-ring-type intensity (nA) and lobe offset (Å).
  # These are the learnable physical parameters.
  [calculators.biot_savart.PHE]
  intensity = -12.00
  lobe_offset = 0.64
  [calculators.biot_savart.TYR]
  intensity = -11.28
  lobe_offset = 0.64
  [calculators.biot_savart.TRP6]
  intensity = -12.48
  lobe_offset = 0.64
  [calculators.biot_savart.TRP5]
  intensity = -6.72
  lobe_offset = 0.52
  [calculators.biot_savart.TRP9]
  intensity = -19.20
  lobe_offset = 0.60
  [calculators.biot_savart.HIE]
  intensity = -5.16
  lobe_offset = 0.50
  [calculators.biot_savart.HID]
  intensity = -5.16
  lobe_offset = 0.50
  [calculators.biot_savart.HIP]
  intensity = -5.16
  lobe_offset = 0.50

  [calculators.haigh_mallion]
  cutoff = 15.0                  # Å — same ring search radius
  subdiv_threshold_L1 = 2.0     # Å — adaptive quadrature level 0→1
  subdiv_threshold_L2 = 1.0     # Å — adaptive quadrature level 1→2

  [calculators.mcconnell]
  cutoff = 10.0                  # Å — bond midpoint search radius

  [calculators.dispersion]
  r_cut = 5.0                   # Å — hard cutoff
  r_switch = 4.3                # Å — CHARMM switching onset

  [calculators.coulomb]
  sanity_limit = 100.0          # V/Å — E-field magnitude clamp

  [calculators.hbond]
  max_dist = 50.0               # Å — maximum N...O distance
  count_radius = 3.5            # Å — proximity counting

  # Pi-quadrupole and ring susceptibility share the ring cutoff.
  # They have no independent distance constants.

# =============================================================
# Filters
# =============================================================

[filters]

  [filters.dipolar_near_field]
  extent_factor = 0.5            # reject if distance < factor × source_extent
  # factor = 0.5 means: skip atoms closer than half the source size.
  # Physics: multipole expansion invalid inside the source distribution.

  [filters.sequential_exclusion]
  min_residue_separation = 2     # skip atom-source pairs within N residues

  [filters.ring_bonded_exclusion]
  enabled = true                 # skip atoms that are ring vertices or bonded to them

# =============================================================
# Cutoff shape (NEW)
# =============================================================
# Applies to all ring calculators (BS, HM, PQ, Disp, RingSusc).
# Default: spherical (a = b = cutoff).
#
# Ellipsoidal in ring-frame cylindrical coordinates:
#   (ρ/a)² + (z/b)² < 1
# where ρ = radial distance in ring plane, z = height above ring.
#
# When a = b, this is a sphere. When b > a, the cutoff extends
# further along the ring normal than in-plane — physically
# motivated by the (3cos²θ - 1)/r³ angular factor.
#
# The spatial index query still uses a spherical radius of
# max(a, b) for the initial candidate set, then the ellipsoidal
# test refines.

[cutoff_shape]
mode = "sphere"                  # "sphere" | "sextic" | "ellipsoid"
# sphere:    r < cutoff.  Standard approach.
# sextic:    r⁶ < C⁶ × (3z² + r²).  Exact iso-|B| for point dipole.
#            Aspect ratio locked at 4^(1/6) = 1.260.  One parameter C.
# ellipsoid: (ρ/a)² + (z/b)² < 1.  Two free parameters.
#            Recovers sphere at a=b, approximates sextic at b/a=1.26.

# Parameters (mode-dependent):
C = 15.0                         # Å — sextic scale parameter
a = 12.0                        # Å — ellipsoid in-plane extent
b = 15.1                        # Å — ellipsoid axial extent (b/a=1.26 = dipolar)
# Smooth taper (optional, 0 = hard boundary):
taper_width = 0.0               # Å — CHARMM-style switching within this shell

# =============================================================
# Output control
# =============================================================

[output]
per_ring = true                  # write per-ring contributions (not just per-type sums)
per_bond = false                 # write per-bond contributions (McConnell)
```

## Per-Ring Output Format

When `output.per_ring = true`, the extractor writes an additional
file per conformation:

### `ring_contributions.npy` — sparse (atom, ring) pair array

Shape: `(P, C)` where P = number of (atom, ring) pairs within
the cutoff, C = number of columns.

Each row is one atom-ring interaction:

```
Column layout (C = 48):
  [0]     atom_index          int → float64
  [1]     ring_index          int → float64
  [2]     ring_type           RingTypeIndex (0-7) → float64
  [3]     distance            Å
  [4]     rho                 Å (radial, in ring plane)
  [5]     z                   Å (signed height above ring)
  [6]     theta               rad (cone angle from normal)
  [7]     mcconnell_factor    (3cos²θ - 1) / r³
  [8]     exp_decay           exp(-distance / 4.0)
  --- BS kernel (9 = SphericalTensor) ---
  [9:18]  bs_G                [T0, T1[3], T2[5]]
  --- HM kernel ---
  [18:27] hm_G                [T0, T1[3], T2[5]]
  --- PiQuadrupole ---
  [27:36] pq_G                [T0, T1[3], T2[5]]
  --- RingSusceptibility ---
  [36:45] chi_G               [T0, T1[3], T2[5]]
  --- Dispersion ---
  [45]    disp_scalar         1/r⁶ (with switching)
  [46]    disp_contacts       number of vertex contacts
  [47]    gaussian_density    Gaussian overlap density
```

### Why sparse?

Different atoms have 0–20 nearby rings. A dense (N_atoms × max_rings × C)
array wastes memory and makes batching awkward. The sparse format has
exactly one row per evaluated kernel pair. Python loads it as:

```python
data = np.load("ring_contributions.npy")  # (P, 48)
atom_idx = data[:, 0].astype(int)
ring_idx = data[:, 1].astype(int)
rho = data[:, 4]
z = data[:, 5]
bs_t2 = data[:, 13:18]   # T2 components of BS kernel
```

Grouping by atom or by ring is a one-line `np.unique` or pandas groupby.

### `ring_geometry.npy` — ring reference data

Shape: `(R, 10)` where R = number of rings in the protein.

```
  [0]     ring_index
  [1]     ring_type           RingTypeIndex (0-7)
  [2]     residue_index
  [3:6]   center              (x, y, z) in Å
  [6:9]   normal              (nx, ny, nz) unit vector
  [9]     radius              Å
```

This allows Python to reconstruct ring-frame coordinates,
verify the C++ (ρ, z) computation, and compute new geometric
features without re-running the extractor.

## Sweep Experiments

### Distance cutoff sweep

```bash
for cutoff in 6 8 10 12 15 18 22 28; do
    sed "s/cutoff = 15.0/cutoff = ${cutoff}/" base.toml > sweep_${cutoff}.toml
    nmr_extract --orca $DIR --output features_${cutoff}/ --config sweep_${cutoff}.toml
done
```

### Cutoff shape investigation

#### Theoretical basis

For a magnetic dipole the iso-|B| surface is NOT an ellipsoid.
It is a sextic algebraic surface:

    (x² + y² + z²)³  =  C² (x² + y² + 4z²)

The polar-to-equatorial aspect ratio is exactly 4^(1/6) = 1.260.
This is a known result from classical electrodynamics — the
dipolar field strength falls off faster in-plane than along the
axis because of the (3cos²θ - 1) angular factor.

No mature computational field (MD, crystallography, FMM, QM)
uses non-spherical interaction cutoffs. This is novel. The
justification is not convention but physics: the iso-field
contour IS anisotropic, and a spherical cutoff at 15 Å either
wastes computation on weak in-plane contributions or misses
strong axial contributions at the same distance.

#### Experimental design

The experiment compares three cutoff shapes against the DFT
meter, carrying all three through to the final NMR prediction:

**Shape 1: Sphere (standard approach, control)**
    r < R.  One parameter R.
    Sweep: R = 6, 8, 10, 12, 15, 18, 22, 28 Å.
    This is what every code in the literature uses.

**Shape 2: Sextic iso-|B| surface (theoretical prediction)**
    r(θ) < C × (3cos²θ + 1)^(1/6).  One parameter C.
    Sweep: C = 6, 8, 10, 12, 15, 18, 22, 28 Å.
    Aspect ratio locked at 1.260 — pure theory, no fitting.
    This is the null hypothesis: if the dipolar model is correct,
    this shape should outperform the sphere.

**Shape 3: Ellipsoid (exploratory)**
    (ρ/a)² + (z/b)² < 1.  Two parameters (a, b).
    Sweep: a = 8, 10, 12, 15; b/a = 1.0, 1.15, 1.26, 1.4, 1.6.
    4 × 5 = 20 configurations.
    Recovers the sphere at b/a = 1.0 and approximates the sextic
    at b/a = 1.26. Deviations from 1.26 indicate physics beyond
    the point-dipole model.

All three shapes are run through the SAME downstream pipeline:
extract → train Level A → train Level C → feed to NMR predictor.
Each produces a final prediction error. The comparison is:

- Sphere vs sextic: does the theoretical shape help?
- Sextic vs ellipsoid: does relaxing the aspect ratio help further?
- Optimal b/a vs 1.260: does the data agree with dipolar theory?

The result is 3 pages of thesis regardless of outcome:
- If sextic beats sphere → angular cutoffs are justified by data.
- If ellipsoid beats sextic → the effective interaction range
  is not purely dipolar (quadrupolar, distributed source, etc.).
- If sphere wins → the cutoff shape doesn't matter at this
  distance range, and a sphere is adequate. Also a finding.

#### Implementation

The sextic filter is cheap: given (ρ, z, r) already computed
in KernelEvaluationContext, the test is:

    r < C × (3 × (z/r)² + 1)^(1/6)

or equivalently:

    r⁶ < C⁶ × (3z² + r²)        [avoids the sixth root]

One comparison, no transcendental functions.

### Near-field extent sweep

```bash
for factor in 0.2 0.3 0.4 0.5 0.6 0.8 1.0 1.2 1.5; do
    # generate toml with extent_factor=$factor
    ...
done
```

9 configurations. The curve should show:
- Too small (0.2): including atoms in source distribution → divergence
- Too large (1.5): excluding real contributions → signal loss
- Minimum somewhere in 0.3–0.8 range

## DistributedRingSusceptibility Calculator

Tests the point-dipole approximation. Instead of one McConnell
dipole at ring center with direction = ring normal, place one
dipole at each ring vertex with direction = ring normal.

```
Existing RingSusc:
    σ(atom) = Σ_rings  M(atom ← ring_center, ring_normal) / r³

Distributed:
    σ(atom) = Σ_rings Σ_vertices  M(atom ← vertex_pos, ring_normal) / r³  ÷ n_vertices
```

The two converge at large distance (point source = distributed source
far away). They diverge within ~2× ring radius. The DIFFERENCE
between them measures the point-dipole error:

```
error(atom) = RingSusc(atom) - DistributedRingSusc(atom)
```

If this difference correlates with the DFT delta T2 residual
(after fitting), the point-dipole approximation is the bottleneck
and a distributed model improves the result. If the difference
does NOT correlate with the residual, the error is elsewhere.

Output: same format as RingSusc (shielding SphericalTensor),
plus the per-ring sparse contribution format if per_ring = true.

## Implementation Order

1. **TOML loading** — extend RuntimeEnvironment to read calculator
   and filter sections. All constants get a TOML key. Absent keys
   use current defaults. No behaviour change without a config file.

2. **Per-ring output** — add WriteRingContributions() that walks
   ConformationAtom::ring_neighbours and writes the sparse array.
   Data already exists in memory; this is pure serialisation.

3. **Ellipsoidal filter** — new KernelEvaluationFilter subclass.
   Takes (a, b) from TOML. Computes (ρ, z) from KernelEvaluationContext
   (which already carries distance and ring geometry). Rejects if
   (ρ/a)² + (z/b)² ≥ 1.

4. **DistributedRingSusceptibility** — new ConformationResult.
   Same structure as RingSusceptibilityResult. Iterates ring vertices
   instead of using ring center. ~100 lines of new code.

5. **Sweep runner** — Python script that generates TOML files,
   runs extractions, trains models, plots curves.

# Extraction SDK — python/nmr_extract/

## What it is

A Python reader for the NPY files that the C++ extractor produces.
Typed objects with `e3nn.o3.Irreps`, `torch.Tensor`, and `numpy`
access.  75 registered arrays covering identity, all classical
calculators (BiotSavart, HaighMallion, McConnell, PiQuadrupole,
Dispersion, RingSusceptibility, HBond, Coulomb, APBS), DSSP + chi
angles, SASA + surface normal, AIMNet2 charges + EFG + embedding,
EEQ charges + coordination number, explicit-water field and
hydration shell + hydration geometry (trajectory), GROMACS energy
(trajectory), MOPAC core + Coulomb + McConnell, Orca DFT,
mutation delta, per-ring sparse contributions, and ring geometry.

## Why it exists

1. The C++ engine is the single source of truth for what gets
   computed.  The Python reader is the single source of truth
   for how those results are consumed.  Both live in this repo.

2. Multiple consumers need the same data: equivariant calibration
   (learn/), upstream GNN (nmr-training/), future work.  They
   install this package as a dependency — they do not copy or
   modify the reader.

3. Trust boundary: when an AI agent works on a consumer project,
   it cannot edit the reader.  The format contract is enforced
   by the fact that the reader code is in a different repo from
   the consumer code.

## What it covers

- **Identity**: positions, elements, residue types
- **Ring current kernels**: Biot-Savart, Haigh-Mallion, pi-quadrupole,
  dispersion, ring susceptibility — each with per-type T0/T2
  decompositions by ring type
- **Per-ring sparse contributions**: (P, 59) table, one row per
  (atom, ring) pair, carrying geometry + all 5 kernel
  SphericalTensors + dispersion scalars + azimuthal angle
- **Ring geometry**: (R, 10) reference table with centers, normals,
  radii
- **Bond anisotropy**: McConnell + MOPAC-McConnell per bond category
- **Electrostatics**: Coulomb + MOPAC-Coulomb E-field + EFG; APBS
  solvated E-field + EFG; AIMNet2 Coulomb EFG (backbone + aromatic)
- **AIMNet2**: Hirshfeld charges, 256-dim AIM embedding, Coulomb EFG
- **H-bond, DSSP** (8-class SS, H-bond energies, chi1-4 dihedrals),
  **SASA** (per-atom Shrake-Rupley)
- **MOPAC electronic structure**: core charges + orbitals + bond orders
- **Orca DFT reference**: total/diamagnetic/paramagnetic shielding
- **Mutation delta**: WT-ALA shielding delta, match metadata,
  removed ring geometry
- **Trajectory-only**: explicit water E-field + EFG (total + first shell),
  hydration shell geometry, GROMACS energy terms per frame

## Tensor types

Every SphericalTensor (shielding, EFG) is 9 components:
`[T0, T1[3], T2[5]]`.  T2 m-ordering matches e3nn: -2,-1,0,+1,+2.

```python
st.irreps       # Irreps("1x0e+1x1o+1x2e")
st.T2            # ndarray (N, 5)
st.torch()       # torch.Tensor (N, 9)
```

Normalization is isometric (Frobenius-preserving), not
orthonormal-on-sphere.  Per-protein normalization in the
calibration pipeline absorbs this.

## Key files

| File | Role |
|------|------|
| `_catalog.py` | Format contract: stem → ArraySpec |
| `_tensors.py` | SphericalTensor, VectorField, scalar types |
| `_ring.py` | RingContributions (P, 59), RingGeometry (R, 10) |
| `_protein.py` | Protein dataclass + `load()` |
| `_types.py` | RingType, BondCategory enums |

## For the detailed API

See **python/API.md** — every field, every column index, every
property, with code examples.

## When to update

Update the SDK when:
- A new `WriteFeatures` method adds an NPY file → add to `_catalog.py`
- A column layout changes (e.g. ring_contributions grew from 48→59) →
  update `_ring.py` column indices + `_catalog.py` cols
- A new calculator is added → add wrapper class, catalog entries,
  group dataclass, wire into `load()`

Always run `python -m pytest python/tests/` after changes.

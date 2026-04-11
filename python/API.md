# nmr_extract API

Reader for NMR shielding tensor extraction data.  Loads NPY files
produced by the C++ extractor into typed Python objects with e3nn
Irreps, torch tensors, and numpy arrays.

## Install

```bash
pip install -e /path/to/nmr-shielding/python
```

Requires numpy, torch, e3nn.

## Load a protein

```python
from nmr_extract import load

p = load("path/to/extraction/directory")
```

Returns a `Protein` with all extracted features.  Errors if required
files are missing.  Warns on unrecognised files.

## Protein fields

```
p.protein_id            str
p.n_atoms               int
p.pos                   VectorField (N, 3)
p.element               ndarray (N,) int32 — atomic number
p.residue_type          ndarray (N,) int32
p.residue_index         ndarray (N,) int32
```

### Ring current calculators

```
p.biot_savart           BiotSavartGroup
  .shielding            ShieldingTensor (N, 9)
  .per_type_T0          PerRingTypeT0 (N, 8)
  .per_type_T2          PerRingTypeT2 (N, 40)
  .total_B              VectorField (N, 3)
  .ring_counts          RingCounts (N, 4)

p.haigh_mallion         RingKernelGroup
p.pi_quadrupole         RingKernelGroup
p.dispersion            RingKernelGroup
  .shielding            ShieldingTensor (N, 9)
  .per_type_T0          PerRingTypeT0 (N, 8)
  .per_type_T2          PerRingTypeT2 (N, 40)

p.ring_susceptibility   ShieldingTensor (N, 9)
```

### Per-ring sparse data

```
p.ring_contributions    RingContributions (P, 59)
p.ring_geometry         RingGeometry (R, 10)
```

P = number of (atom, ring) pairs evaluated.  R = number of rings.

### Bond calculators

```
p.mcconnell             McConnellGroup
  .shielding            ShieldingTensor (N, 9)
  .category_T2          PerBondCategoryT2 (N, 25)
  .scalars              McConnellScalars (N, 6)

p.coulomb               CoulombGroup
  .shielding            ShieldingTensor (N, 9)
  .E                    VectorField (N, 3)
  .efg_backbone         EFGTensor (N, 9)
  .efg_aromatic         EFGTensor (N, 9)
  .scalars              CoulombScalars (N, 4)

p.hbond                 HBondGroup
  .shielding            ShieldingTensor (N, 9)
  .scalars              HBondScalars (N, 3)

p.dssp                  DsspScalars (N, 5)
```

### Optional groups (None if not extracted)

```
p.mopac                 MopacGroup | None
  .core.charges         ndarray (N,)
  .core.scalars         MopacScalars (N, 4)
  .core.bond_orders     BondOrders (B, 3)
  .core.global_         MopacGlobal (4,)    — [hof, dipole_x, dipole_y, dipole_z]
  .coulomb              MopacCoulombGroup — same as CoulombGroup
  .mcconnell            MopacMcConnellGroup — same as McConnellGroup

p.apbs                  APBSGroup | None
  .E                    VectorField (N, 3)
  .efg                  EFGTensor (N, 9)

p.orca                  OrcaGroup | None
  .total                ShieldingTensor (N, 9)
  .diamagnetic          ShieldingTensor (N, 9)
  .paramagnetic         ShieldingTensor (N, 9)

p.delta                 DeltaGroup | None
  .shielding            ShieldingTensor (N, 9)
  .scalars              DeltaScalars (N, 6)
  .apbs                 DeltaAPBS (N, 12) | None
  .ring_proximity       DeltaRingProximity (N, R*6)

p.aimnet2               AIMNet2Group | None
  .charges              AIMNet2Charges (N,) — Hirshfeld charges
  .aim                  AIMNet2AimEmbedding (N, 256) — electronic structure embedding
  .efg                  EFGTensor (N, 9) — Coulomb EFG from AIMNet2 charges
  .efg_aromatic         EFGTensor (N, 9) — aromatic decomposition
  .efg_backbone         EFGTensor (N, 9) — backbone decomposition
  .charge_sensitivity   AIMNet2ChargeSensitivity (N,) — intrinsic polarisability proxy
```

## Tensor types

### SphericalTensor (and ShieldingTensor, EFGTensor)

9-component packing: [T0, T1[3], T2[5]].

```python
st = p.biot_savart.shielding

st.data                 # ndarray (N, 9)
st.torch()              # torch.Tensor (N, 9)
st.irreps               # Irreps("1x0e+1x1o+1x2e")

st.T0                   # ndarray (N, 1)
st.T1                   # ndarray (N, 3)
st.T2                   # ndarray (N, 5) — m = -2, -1, 0, +1, +2
st.isotropic            # ndarray (N,) — T0 squeezed
st.T2_magnitude         # ndarray (N,) — L2 norm of T2
```

T2 component ordering matches e3nn.  No permutation needed.

### VectorField

```python
v = p.coulomb.E

v.data                  # ndarray (N, 3)
v.torch()               # torch.Tensor (N, 3)
v.irreps                # Irreps("1x1o")
v.x, v.y, v.z           # ndarray (N,)
v.magnitude             # ndarray (N,)
```

### PerRingTypeT2

8 ring types x 5 T2 components = 40.

```python
t2 = p.biot_savart.per_type_T2

t2.irreps               # Irreps("8x2e")
t2.for_type(RingType.PHE)   # ndarray (N, 5)
t2.as_block()           # ndarray (N, 8, 5)
t2.total                # ndarray (N, 5) — sum over types
```

### PerBondCategoryT2

5 bond categories x 5 T2 components = 25.

```python
mc = p.mcconnell.category_T2

mc.irreps               # Irreps("5x2e")
mc.for_category(BondCategory.CO_nearest)  # ndarray (N, 5)
mc.as_block()           # ndarray (N, 5, 5)
```

## RingContributions

Sparse (P, 59) table — one row per (atom, ring) pair.

```python
rc = p.ring_contributions

rc.n_pairs              # int
rc.atom_index           # ndarray (P,) int
rc.ring_index           # ndarray (P,) int
rc.ring_type            # ndarray (P,) int — RingType enum values
rc.distance             # ndarray (P,) — Angstroms
rc.rho, rc.z, rc.theta  # cylindrical coordinates in ring frame
rc.cos_phi, rc.sin_phi  # azimuthal angle in ring plane (relative to vertex 0)

# Physics kernels — each is a SphericalTensor (P, 9)
rc.bs                   # Biot-Savart shielding kernel G
rc.hm_H                 # Haigh-Mallion raw integral H (pure T2)
rc.hm                   # Haigh-Mallion shielding kernel G (intensity * H)
rc.pq                   # Pi-quadrupole
rc.chi                  # Ring susceptibility
rc.disp_scalar          # ndarray (P,) — 1/r^6

# T2 from any kernel
rc.bs.T2                # ndarray (P, 5)
rc.bs.torch()           # torch.Tensor (P, 9)

# Filter
rc.for_atom(42)         # RingContributions — rows for atom 42
rc.for_ring_type(RingType.PHE)  # rows for PHE rings only
```

### Column layout (59 columns)

```
[0]     atom_index
[1]     ring_index
[2]     ring_type           0-7
[3]     distance            Angstroms
[4]     rho                 Angstroms
[5]     z                   Angstroms (signed)
[6]     theta               radians
[7]     mcconnell_factor    (3cos^2 theta - 1) / r^3
[8]     exp_decay           exp(-distance / 4.0)
[9:18]  bs_G                SphericalTensor — BS shielding kernel
[18:27] hm_H                SphericalTensor — HM raw integral (pure T2)
[27:36] hm_G                SphericalTensor — HM shielding kernel (intensity * H)
[36:45] pq_G                SphericalTensor
[45:54] chi_G               SphericalTensor
[54]    disp_scalar         1/r^6
[55]    disp_contacts       vertex contact count
[56]    gaussian_density    (placeholder)
[57]    cos_phi             azimuthal angle cosine (relative to vertex 0)
[58]    sin_phi             azimuthal angle sine (relative to vertex 0)
```

## RingGeometry

Reference table for ring identity and spatial properties.

```python
rg = p.ring_geometry

rg.n_rings              # int
rg.ring_index           # ndarray (R,) int
rg.ring_type            # ndarray (R,) int
rg.residue_index        # ndarray (R,) int
rg.center               # ndarray (R, 3)
rg.normal               # ndarray (R, 3) — unit vectors
rg.radius               # ndarray (R,)
```

## Enums

```python
from nmr_extract import RingType, BondCategory

RingType.PHE            # 0
RingType.TYR            # 1
RingType.TRP_benzene    # 2
RingType.TRP_pyrrole    # 3
RingType.TRP_perimeter  # 4
RingType.HIS            # 5
RingType.HID            # 6
RingType.HIE            # 7

BondCategory.backbone_total   # 0
BondCategory.sidechain_total  # 1
BondCategory.aromatic_total   # 2
BondCategory.CO_nearest       # 3
BondCategory.CN_nearest       # 4
```

## Catalog

The format contract.  Maps every NPY filename to its metadata.

```python
from nmr_extract import CATALOG

for stem, spec in CATALOG.items():
    print(f"{stem:30s} {spec.group:20s} required={spec.required}")
```

## HM: two representations

Haigh-Mallion provides both forms per ring:

- `rc.hm_H` — raw surface integral H (pure T2, geometry only)
- `rc.hm` — shielding kernel G = intensity x H (T0+T1+T2)

`rc.hm` sums match `p.haigh_mallion.per_type_T2`.  `rc.hm_H` gives
the unscaled geometric kernel for analysis.

# Categorical Analysis Pass 1: EquivariantExtractor v1

Every operation, decision, comparison, and extraction criterion in the v1
EquivariantExtractor, enumerated line by line.

Source files:
- `src/features/EquivariantExtractor.h`
- `src/features/EquivariantExtractor.cpp` (1740 lines)

---

## Anonymous-namespace helper functions (lines 13-88)

### `is_backbone_atom(name)` (line 15)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 1 | String equality test: name == "N" | atom_name string | Categorical: is this the backbone nitrogen? | bool (OR'd into result) | n/a (helper) | none |
| 2 | String equality test: name == "CA" | atom_name string | Categorical: is this the alpha carbon? | bool (OR'd into result) | n/a | none |
| 3 | String equality test: name == "C" | atom_name string | Categorical: is this the backbone carbonyl carbon? | bool (OR'd into result) | n/a | none |
| 4 | String equality test: name == "O" | atom_name string | Categorical: is this the backbone carbonyl oxygen? | bool (OR'd into result) | n/a | none |
| 5 | String equality test: name == "H" | atom_name string | Categorical: is this the backbone amide hydrogen? | bool (OR'd into result) | n/a | none |
| 6 | String equality test: name == "HA" | atom_name string | Categorical: is this the alpha hydrogen? | bool (OR'd into result) | n/a | none |
| 7 | String equality test: name == "HN" | atom_name string | Categorical: is this the alternate amide hydrogen name? | bool (OR'd into result) | n/a | none |

Note: The backbone atom set is {N, CA, C, O, H, HA, HN}. This is used in AtomSite::build and hbond_features. The set in `aromatic_field_features` is different (broader, includes CB, HB, HB2, HB3, HA2, HA3).

### `check_amide_H(name)` (line 20)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 8 | String equality test: name == "H" | atom_name string | Categorical: PDB standard amide H | bool (OR'd) | n/a | none |
| 9 | String equality test: name == "HN" | atom_name string | Categorical: alternate amide H name | bool (OR'd) | n/a | none |

### `check_alpha_H(name)` (line 24)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 10 | String equality test: name == "HA" | atom_name string | Categorical: standard alpha H | bool (OR'd) | n/a | none |
| 11 | String equality test: name == "HA2" | atom_name string | Categorical: glycine alpha H variant 1 | bool (OR'd) | n/a | none |
| 12 | String equality test: name == "HA3" | atom_name string | Categorical: glycine alpha H variant 2 | bool (OR'd) | n/a | none |

### `check_methyl(name)` (lines 28-40)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 13 | Length check: name.size() < 3 | atom_name string | Decision: name too short to be methyl H | early return false | n/a | none |
| 14 | Character test: name[0] != 'H' | atom_name string | Decision: must start with H | early return false | n/a | none |
| 15 | Read second character: name[1] | atom_name string | Extract parent-letter | char `second` | n/a | none |
| 16 | Read third character: name[2] | atom_name string | Extract position indicator | char `third` | n/a | none |
| 17 | Character test: second in {B, G, D, E} | char `second` | Decision: parent must be beta/gamma/delta/epsilon carbon | early return false if not | n/a | none |
| 18 | Digit test: std::isdigit(third) | char `third` | Decision: third character must be a digit (methyl numbering) | bool result | n/a | none |

Comment notes previous version incorrectly flagged HD1/HD2 (His ring H), HE1/HE2 (Trp ring H), HE (Tyr OH H), HD21/HD22 (Asn amide H) as methyl.

### `ring_has_nitrogen(ring)` (line 42)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 19 | Iterate ring.atom_names | ring.atom_names vector | Loop over all atom names in ring | iteration | n/a | none |
| 20 | Character test: name[0] == 'N' for each name | atom name first character | Decision: is this atom a nitrogen? (checks only first character, not full element) | bool | n/a | none |

### `ring_type_idx(ring)` (line 48)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 21 | Delegate to RingCurrentParams::type_index | ring.residue_name, ring.ring_type | Maps (residue_name, ring_type) to integer 0-6. PHE=0, TYR=1, TRP+6-ring=2, TRP+5-ring=3, HIS/HIP=4, HID=5, HIE=6 | int type index | n/a | none |

Note: RingCurrentParams::type_index returns -1 for unrecognized residue names. TRP without matching ring_type defaults to 2 (6-ring).

### `to_ring_cylindrical(v, center, normal, point)` (line 53)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 22 | Compute displacement d = point - center | point (Vec3), center (Vec3) | none | Vec3 d | n/a | none |
| 23 | Project onto normal: z_comp = d.dot(normal) | d, normal | none | double z_comp | n/a | op 22 |
| 24 | Compute in-plane component: in_plane = d - z_comp * normal | d, z_comp, normal | none | Vec3 in_plane | n/a | ops 22-23 |
| 25 | Compute radial distance: rho = in_plane.norm() | in_plane | none | double rho | n/a | op 24 |
| 26 | Degenerate check: rho < 1e-10 | rho | Decision: atom on ring axis (no radial direction defined) | if true: return (v.dot(normal), 0, 0) | n/a | op 25 |
| 27 | Compute radial unit vector: e_rho = in_plane / rho | in_plane, rho | none | Vec3 e_rho | n/a | ops 24-25 |
| 28 | Compute azimuthal unit vector: e_phi = normal.cross(e_rho) | normal, e_rho | none | Vec3 e_phi | n/a | op 27 |
| 29 | Project v onto cylindrical basis: (v.dot(normal), v.dot(e_rho), v.dot(e_phi)) | v, normal, e_rho, e_phi | none | Vec3 (B_n, B_rho, B_phi) | n/a | ops 27-28 |

Magic number: 1e-10 threshold for rho degeneracy.

### `is_sp2_parent(parent_name, residue_name)` (line 68)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 30 | String test: parent_name == "N" | parent atom_name | Decision: backbone N is always sp2 | early return true | n/a | none |
| 31 | String test: parent_name == "C" | parent atom_name | Decision: backbone C (carbonyl) is always sp2 | early return true | n/a | none |
| 32 | Call is_aromatic_residue(residue_name) | residue_name | Decision: is the residue aromatic? (PHE/TYR/TRP/HIS/HID/HIE/HIP) | bool | n/a | none |
| 33 | If aromatic: string tests for 13 atom names (CG, CD1, CD2, CE1, CE2, CE3, CZ, CZ2, CZ3, CH2, NE1, ND1, NE2) | parent_name | Decision: is parent one of the known ring/adjacent sp2 carbons or nitrogens? | bool | n/a | op 32 |

Note: CA is NOT in this list, so alpha H is never flagged as sp2-parented. CB also not included.

### `safe_unit(v)` (line 83)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 34 | Compute norm: n = v.norm() | Vec3 v | none | double n | n/a | none |
| 35 | Degenerate check: n > 1e-10 | n | Decision: is vector non-zero? | Vec3 (v/n) or Vec3::Zero() | n/a | op 34 |

Magic number: 1e-10 threshold for zero-vector.

---

## AtomSite::build (lines 95-201)

Factory method constructing the observation-point context for one atom.

### Distance-sorted ring list (lines 102-124)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 36 | Reserve vector capacity: rings.size() | rings vector | none | memory allocation | n/a | none |
| 37 | For each ring i: compute distance (atom.position - rings[i].center).norm() | atom.position, rings[i].center | none | double d | n/a | none |
| 38 | Push {index=i, distance=d} into rings_by_distance | i, d | none | RingNeighbor appended | n/a | op 37 |
| 39 | Sort rings_by_distance ascending by distance | rings_by_distance | none | sorted vector | n/a | op 38 |
| 40 | For each RingNeighbor: if distance <= 3.0, increment n_rings_within_3A | rn.distance | Threshold: 3.0 Angstroms | int counter | n/a | op 39 |
| 41 | For each RingNeighbor: if distance <= 5.0, increment n_rings_within_5A | rn.distance | Threshold: 5.0 Angstroms | int counter | n/a | op 39 |
| 42 | For each RingNeighbor: if distance <= 8.0, increment n_rings_within_8A | rn.distance | Threshold: 8.0 Angstroms | int counter | n/a | op 39 |
| 43 | For each RingNeighbor: if distance <= 12.0, increment n_rings_within_12A | rn.distance | Threshold: 12.0 Angstroms | int counter | n/a | op 39 |
| 44 | Check: rings_by_distance not empty | rings_by_distance.empty() | Decision: any rings at all? | controls mean computation | n/a | op 39 |
| 45 | Sum all ring distances | rn.distance for all rings | none | double sum | n/a | op 39 |
| 46 | Divide sum by count to get mean_ring_distance | sum, count | none | double mean_ring_distance | n/a | op 45 |

Note: If rings_by_distance is empty, mean_ring_distance stays at its default: 99.0 (sentinel).

### Atom classification (lines 127-132)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 47 | Element test: atom.element == "H" | atom.element | Categorical: is this atom hydrogen? | bool is_hydrogen | n/a | none |
| 48 | Call is_backbone_atom(atom.atom_name) | atom.atom_name | Categorical: is atom in backbone set? (ops 1-7) | bool is_backbone | n/a | none |
| 49 | Call check_amide_H(atom.atom_name) | atom.atom_name | Categorical: is atom amide H? (ops 8-9) | bool is_amide_H | n/a | none |
| 50 | Call check_alpha_H(atom.atom_name) | atom.atom_name | Categorical: is atom alpha H? (ops 10-12) | bool is_alpha_H | n/a | none |
| 51 | Call check_methyl(atom.atom_name) | atom.atom_name | Categorical: is atom methyl H? (ops 13-18) | bool is_methyl | n/a | none |
| 52 | Call is_aromatic_residue(atom.residue_name) | atom.residue_name | Categorical: is residue in {PHE, TYR, TRP, HIS, HID, HIE, HIP}? | bool is_on_aromatic_residue | n/a | none |

### Parent heavy atom search (lines 134-149)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 53 | Guard: only run if is_hydrogen | is_hydrogen | Decision: skip for heavy atoms | conditional block | n/a | op 47 |
| 54 | Initialize best_d = 1.8 | hardcoded constant | Magic number: 1.8 Angstrom H-X bond length cutoff | double best_d | n/a | none |
| 55 | For each pdb_atom pa: skip if pa.element == "H" | pa.element | Decision: only consider heavy atoms as parents | skip | n/a | none |
| 56 | Compute distance: d = (pa.position - atom.position).norm() | pa.position, atom.position | none | double d | n/a | none |
| 57 | Update nearest: if d < best_d, update best_d and parent pointer | d, best_d | Decision: is this the closest heavy atom within 1.8 A? | AtomInfo* parent | n/a | op 56 |
| 58 | Guard: only proceed if parent was found (non-null) | parent pointer | Decision: did we find a parent? | conditional block | n/a | op 57 |
| 59 | Compute bond direction: safe_unit(atom.position - parent->position) | atom.position, parent->position | none (uses safe_unit op 34-35) | Vec3 bond_direction | n/a | op 57 |
| 60 | Compute sp2 flag: is_sp2_parent(parent->atom_name, atom.residue_name) | parent->atom_name, atom.residue_name | Categorical: is parent sp2 hybridized? (ops 30-33) | bool parent_is_sp2 | n/a | op 57 |

Note: The parent search iterates ALL pdb_atoms. It is a brute-force O(N) search. The bond direction points from parent TO hydrogen (atom.position - parent->position).

### Nearest H-bond acceptor (lines 152-162)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 61 | For each pdb_atom pa: check same-residue exclusion (pa.residue_number == atom.residue_number AND pa.chain_id == atom.chain_id) | pa.residue_number, pa.chain_id, atom.residue_number, atom.chain_id | Decision: skip atoms on the same residue | skip | n/a | none |
| 62 | Element filter: pa.element must be "O" or "N" | pa.element | Decision: only O and N are acceptors | skip if neither | n/a | none |
| 63 | Compute distance: d = (pa.position - atom.position).norm() | pa.position, atom.position | none | double d | n/a | none |
| 64 | Update nearest: if d < nearest_acceptor_dist | d, nearest_acceptor_dist | Decision: is this the closest acceptor found? | double nearest_acceptor_dist, Vec3 dir_to_nearest_acceptor | n/a | op 63 |
| 65 | Compute direction: safe_unit(pa.position - atom.position) | pa.position, atom.position | none | Vec3 dir_to_nearest_acceptor | n/a | op 64 |

Note: This finds the nearest O or N on a DIFFERENT residue, regardless of whether the current atom is a donor. The acceptor search makes no distinction between N-as-donor and N-as-acceptor. Default: nearest_acceptor_dist = 99.0 (sentinel), dir = Zero.

### Nearest backbone carbonyl C=O (lines 165-181)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 66 | For each pdb_atom pa: filter by atom_name == "O" | pa.atom_name | Decision: only consider atoms named "O" (backbone oxygen) | skip if not | n/a | none |
| 67 | Compute distance: d = (pa.position - atom.position).norm() | pa.position, atom.position | none | double d | n/a | none |
| 68 | Update nearest: if d < nearest_co_dist | d, nearest_co_dist | Decision: closest backbone O so far? | double nearest_co_dist, Vec3 co_O_pos | n/a | op 67 |
| 69 | Store O position: co_O_pos = pa.position | pa.position | none | Vec3 co_O_pos | n/a | op 68 |
| 70 | Inner loop: for each pdb_atom ca, find backbone C in same residue | ca.residue_number == pa.residue_number, ca.chain_id == pa.chain_id, ca.atom_name == "C", ca.element == "C" | Decision: is this the carbonyl carbon of the same residue? | Vec3 co_C_pos | n/a | op 68 |

Note: The search for backbone C uses BOTH atom_name == "C" AND element == "C". This double check distinguishes the backbone carbonyl carbon from any other atom. But the O search only checks atom_name == "O" without checking if it's backbone (no element check, no is_backbone check). Any atom named "O" in any residue will match, including sidechain O atoms in Thr (OG1), Ser (OG), Tyr (OH). However, PDB convention names backbone O as "O" and sidechain oxygens with different names (OD1, OE1, OG, etc.), so in practice this works for standard PDB files. The atom_name "O" is reliably backbone only.

Default: nearest_co_dist = 99.0 (sentinel).

### Packing density (lines 184-188)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 71 | For each pdb_atom pa: skip if pa.element == "H" | pa.element | Decision: only count heavy atoms | skip | n/a | none |
| 72 | Compute squared distance: (pa.position - atom.position).squaredNorm() | pa.position, atom.position | none | double squared_dist | n/a | none |
| 73 | Threshold check: squared_dist < 64.0 (i.e., distance < 8.0 A) | squared_dist | Decision: within 8 Angstrom shell? | increment n_heavy_within_8A | n/a | op 72 |

Magic number: 64.0 = 8.0^2 (squared distance threshold for 8 A).

Note: This count INCLUDES the atom itself (no self-exclusion).

### Nearest ring atom distance (lines 191-198)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 74 | Guard: rings_by_distance not empty | rings_by_distance.empty() | Decision: any rings? | conditional block | n/a | op 39 |
| 75 | Get nearest ring: rings[rings_by_distance[0].index] | rings_by_distance[0].index | none | RingInfo& nearest | n/a | op 39 |
| 76 | For each vertex rv in nearest ring: compute (rv - atom.position).norm() | rv, atom.position | none | double d | n/a | op 75 |
| 77 | Update minimum: if d < min_ring_atom_dist | d, min_ring_atom_dist | Decision: closest ring atom vertex so far? | double min_ring_atom_dist | n/a | op 76 |

Note: Only checks vertices of the NEAREST ring (by center distance). Default: min_ring_atom_dist = 99.0 (sentinel).

---

## ring_current_features (lines 220-450)

### Per-ring accumulation loop (lines 237-264)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 78 | Initialize 7-element arrays sum_G_T0 and sum_G_T2 to zero | hardcoded NT=7 | none | arrays | n/a | none |
| 79 | Initialize B_total = Vec3::Zero() | none | none | Vec3 | n/a | none |
| 80 | Initialize G_T2_exp_sum to zero, G_iso_sum_exp to 0 | none | none | accumulators | n/a | none |
| 81 | For each ring neighbor rn (sorted by distance): look up ring | rn.index, rings vector | none | RingInfo& ring | n/a | AtomSite.rings_by_distance |
| 82 | Compute ring type index: ring_type_idx(ring) | ring.residue_name, ring.ring_type | Categorical: maps to 0-6 (op 21) | int tidx | n/a | none |
| 83 | Compute JB B-field: johnson_bovey_field(ring.vertices, ring.normal, jb_offset, 1.0, pos) | ring geometry, atom position, type-specific JB lobe offset | Physics: Johnson-Bovey double-loop B-field at atom position, with unit current (1.0) | Vec3 B | L1e (accumulated into B_total) | ops 37, 82 |
| 84 | Accumulate B_total += B | B, B_total | none | Vec3 B_total | n/a | op 83 |
| 85 | Compute geometric factor tensor: ring_geometric_factor_jb(ring.vertices, ring.normal, jb_offset, pos) | ring geometry, atom position, JB offset | Physics: full 3x3 shielding geometric factor from JB model | Mat3 G | n/a | op 82 |
| 86 | Decompose tensor: decompose_tensor(G) | Mat3 G | Spherical decomposition: T0 (scalar) + T2 (5 components) | SphericalDecomposition sd | n/a | op 85 |
| 87 | Type-index bounds check: tidx >= 0 && tidx < NT | tidx | Decision: valid ring type? | conditional accumulation | n/a | op 82 |
| 88 | Accumulate sum_G_T0[tidx] += sd.T0 | sd.T0, tidx | none | per-type scalar sum | n/a | ops 86-87 |
| 89 | Accumulate sum_G_T2[tidx][m] += sd.T2[m] for m in 0..4 | sd.T2, tidx | none | per-type tensor sum | n/a | ops 86-87 |
| 90 | Compute exponential decay: decay = exp(-rn.distance / 4.0) | rn.distance | none (decay length = 4.0 Angstroms) | double decay | n/a | none |
| 91 | Accumulate G_iso_sum_exp += sd.T0 * decay | sd.T0, decay | none | double | n/a | ops 86, 90 |
| 92 | Accumulate G_T2_exp_sum[m] += sd.T2[m] * decay for m in 0..4 | sd.T2, decay | none | array | n/a | ops 86, 90 |
| 93 | Threshold check: if rn.distance <= 8.0, push sd.T0 into G_iso_near | rn.distance | Threshold: 8.0 Angstroms for variance computation | vector append | n/a | op 86 |

Magic numbers: 4.0 (exponential decay length in Angstroms), 8.0 (radius for G_iso variance).

### L=0 type-summed isotropic factors [7 features] (lines 269-270)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 94 | Push sum_G_T0[0] through sum_G_T0[6] (PHE, TYR, TRP6, TRP5, HIS, HID, HIE) | sum_G_T0 array | none | 7 L0 scalars: jb_sum_G_T0_{type} | L0 | ops 88 |

### L=0 B-field magnitude [1 feature] (line 273)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 95 | Compute B_total.norm() | B_total | none | L0 scalar: jb_B_total_mag | L0 | op 84 |

### L=0 ring environment counts [4 features] (lines 276-279)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 96 | Cast n_rings_within_3A to double | site.n_rings_within_3A | none | L0 scalar: ring_n_within_3A | L0 | op 40 |
| 97 | Cast n_rings_within_5A to double | site.n_rings_within_5A | none | L0 scalar: ring_n_within_5A | L0 | op 41 |
| 98 | Cast n_rings_within_8A to double | site.n_rings_within_8A | none | L0 scalar: ring_n_within_8A | L0 | op 42 |
| 99 | Cast n_rings_within_12A to double | site.n_rings_within_12A | none | L0 scalar: ring_n_within_12A | L0 | op 43 |

### L=0 ring distances [4 features] (lines 282-289)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 100 | Lambda rd(k): if rings_by_distance.size() > k, return distance, else 99.0 | rings_by_distance, k | Decision: does k-th ring exist? Sentinel: 99.0 | double | n/a | op 39 |
| 101 | Push rd(0) | 1st ring distance or 99.0 | none | L0 scalar: ring_dist_1st | L0 | op 100 |
| 102 | Push rd(1) | 2nd ring distance or 99.0 | none | L0 scalar: ring_dist_2nd | L0 | op 100 |
| 103 | Push rd(2) | 3rd ring distance or 99.0 | none | L0 scalar: ring_dist_3rd | L0 | op 100 |
| 104 | Push site.mean_ring_distance | mean or 99.0 (default) | none | L0 scalar: ring_mean_dist | L0 | op 46 |

### L=0 nearest ring (#1) ring-frame scalars [20 features] (lines 292-330)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 105 | Guard: rings_by_distance not empty | rings_by_distance.empty() | Decision: any rings? True branch: compute 20 features. False branch: push 20 zeros | conditional | n/a | op 39 |
| 106 | Get ring1: rings[rings_by_distance[0].index] | sorted ring list | none | RingInfo& ring1 | n/a | op 39 |
| 107 | Get r1: rings_by_distance[0].distance | sorted ring list | none | double r1 | n/a | op 39 |
| 108 | Compute displacement: d = pos - ring1.center | atom position, ring center | none | Vec3 d | n/a | none |
| 109 | Compute axial component: z = d.dot(ring1.normal) | d, ring1.normal | none | double z (signed height above/below ring plane) | n/a | op 108 |
| 110 | Compute in-plane component: ip = d - z * ring1.normal | d, z, ring1.normal | none | Vec3 ip | n/a | ops 108-109 |
| 111 | Compute radial distance: rho = ip.norm() | ip | none | double rho | n/a | op 110 |
| 112 | Compute polar angle: theta = atan2(rho, z) | rho, z | none | double theta (angle from ring normal) | n/a | ops 109, 111 |
| 113 | Compute r^3: r3 = r1 * r1 * r1 | r1 | none | double r3 | n/a | op 107 |
| 114 | Push r1 | r1 | none | L0 scalar: near1_r | L0 | op 107 |
| 115 | Push rho | rho | none | L0 scalar: near1_rho | L0 | op 111 |
| 116 | Push z | z (signed) | none | L0 scalar: near1_z | L0 | op 109 |
| 117 | Push theta | theta | none | L0 scalar: near1_theta | L0 | op 112 |
| 118 | Compute and push 1/r3 | r3 | none | L0 scalar: near1_inv_r3 | L0 | op 113 |
| 119 | Compute cos_theta = cos(theta) | theta | none | double cos_theta | n/a | op 112 |
| 120 | Compute and push McConnell factor: (3*cos_theta^2 - 1) / r3 | cos_theta, r3 | none | L0 scalar: near1_mcconnell | L0 | ops 113, 119 |
| 121 | Compute B1: johnson_bovey_field for ring1 at pos | ring1 geometry, JB offset for ring1's type, atom position | Physics: JB B-field from nearest ring only | Vec3 B1 | n/a | op 106 |
| 122 | Transform B1 to cylindrical: to_ring_cylindrical(B1, ring1.center, ring1.normal, pos) | B1, ring1 geometry, atom position | Coordinate transform (ops 22-29) | Vec3 B1_cyl = (B_n, B_rho, B_phi) | n/a | op 121 |
| 123 | Push B1_cyl(0) | B_n component | none | L0 scalar: near1_Bn | L0 | op 122 |
| 124 | Push B1_cyl(1) | B_rho component | none | L0 scalar: near1_Brho | L0 | op 122 |
| 125 | Push B1_cyl(2) | B_phi component | none | L0 scalar: near1_Bphi | L0 | op 122 |
| 126 | Compute G1: ring_geometric_factor_jb for ring1 at pos | ring1 geometry, JB offset, atom position | Physics: geometric factor tensor | Mat3 G1 | n/a | op 106 |
| 127 | Decompose: sd1 = decompose_tensor(G1) | G1 | Spherical decomposition | SphericalDecomposition sd1 | n/a | op 126 |
| 128 | Push sd1.T0 | T0 of nearest ring | none | L0 scalar: near1_G_iso | L0 | op 127 |
| 129 | Push ring1.radius | effective ring radius | none | L0 scalar: near1_radius | L0 | op 106 |
| 130 | Check ring_has_nitrogen(ring1) and push 1.0 or 0.0 | ring1.atom_names | Categorical: does nearest ring contain nitrogen? (op 19-20) | L0 scalar: near1_has_N | L0 | op 106 |
| 131 | Compute ring type one-hot: for t in 0..6, push (t == tidx1) ? 1.0 : 0.0 | ring_type_idx(ring1) | Categorical: 7-class one-hot of nearest ring type | 7 L0 scalars: near1_is_{PHE,TYR,TRP6,TRP5,HIS,HID,HIE} | L0 | op 106 |
| 132 | Compute and push exp(-r1 / 4.0) | r1 | none (decay length = 4.0 A) | L0 scalar: near1_exp_decay | L0 | op 107 |

Magic numbers: 4.0 (decay length). False branch: pushes 20 zeros via a size-tracking while loop (line 328-329). The expected count comment says "5+1+3+1+1+1+7+1 = 20" which matches 12 named scalars + 7 one-hot + 1 decay = 20.

### L=0 second nearest ring (#2) key scalars [4 features] (lines 333-352)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 133 | Guard: rings_by_distance.size() >= 2 | rings count | Decision: is there a second ring? | conditional | n/a | op 39 |
| 134 | Get ring2: rings[rings_by_distance[1].index] | sorted ring list | none | RingInfo& ring2 | n/a | op 39 |
| 135 | Get r2: rings_by_distance[1].distance | sorted ring list | none | double r2 | n/a | op 39 |
| 136 | Clamp r2: if r2 < MIN_DISTANCE (0.1), set r2 = MIN_DISTANCE | r2, constants::MIN_DISTANCE | Decision: prevent division by zero | double r2 (clamped) | n/a | op 135 |
| 137 | Compute d2 = pos - ring2.center | atom position, ring2 center | none | Vec3 d2 | n/a | none |
| 138 | Compute z2 = d2.dot(ring2.normal) | d2, ring2.normal | none | double z2 | n/a | op 137 |
| 139 | Compute r2_3 = r2^3 | r2 | none | double r2_3 | n/a | op 136 |
| 140 | Compute cos2 = z2 / r2 | z2, r2 | none | double cos2 | n/a | ops 136, 138 |
| 141 | Compute G2: ring_geometric_factor_jb for ring2 at pos | ring2 geometry, JB offset, atom position | Physics: geometric factor tensor for 2nd ring | Mat3 G2 | n/a | op 134 |
| 142 | Decompose: sd2 = decompose_tensor(G2) | G2 | Spherical decomposition | SphericalDecomposition sd2 | n/a | op 141 |
| 143 | Push r2 | 2nd ring distance | none | L0 scalar: near2_r | L0 | op 136 |
| 144 | Push z2 | 2nd ring axial height | none | L0 scalar: near2_z | L0 | op 138 |
| 145 | Push sd2.T0 | 2nd ring isotropic geometric factor | none | L0 scalar: near2_G_iso | L0 | op 142 |
| 146 | Compute and push (3*cos2^2 - 1) / r2_3 | cos2, r2_3 | none | L0 scalar: near2_mcconnell | L0 | ops 139-140 |

False branch: pushes 4 zeros.

### L=0 multi-ring geometry [5 features] (lines 355-381)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 147 | Guard: rings_by_distance.size() >= 2 | rings count | Decision: are there at least 2 rings? | conditional | n/a | op 39 |
| 148 | Compute n1_dot_n2 = r1.normal.dot(r2.normal) | normals of 2 nearest rings | none | double (cosine of angle between ring planes) | n/a | op 39 |
| 149 | Compute n1_cross_n2 = r1.normal.cross(r2.normal).norm() | normals of 2 nearest rings | none | double (sine of angle between ring planes) | n/a | op 39 |
| 150 | Compute d2_d1_ratio: if dist[0] > 1e-10, ratio = dist[1]/dist[0], else 0 | 1st and 2nd ring distances | Decision: avoid division by zero (threshold 1e-10) | double | n/a | op 39 |
| 151 | Push n1_dot_n2 | dot product | none | L0 scalar: ring_n1_dot_n2 | L0 | op 148 |
| 152 | Push n1_cross_n2 | cross product norm | none | L0 scalar: ring_n1_cross_n2 | L0 | op 149 |
| 153 | Push G_iso_sum_exp | exponentially weighted sum of isotropic factors | none | L0 scalar: ring_G_iso_sum_exp | L0 | op 91 |
| 154 | Compute variance of G_iso within 8 A: if G_iso_near.size() >= 2, compute sample variance (with N denominator, not N-1) | G_iso_near vector | Decision: need at least 2 values for variance | double G_iso_var | n/a | op 93 |
| 155 | Push G_iso_var | variance | none | L0 scalar: ring_G_iso_var_8A | L0 | op 154 |
| 156 | Push d2_d1_ratio | distance ratio | none | L0 scalar: ring_d2_d1_ratio | L0 | op 150 |

Note: Variance uses population formula (divide by N), not sample variance (N-1).

### L=1 B-fields (pseudovectors, 1e) [3 features] (lines 387-401)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 157 | Push B_total | total JB B-field over all rings | none | L1 vector: jb_B_total | L1e | op 84 |
| 158 | Guard: rings_by_distance not empty for B_near1 | rings count | Decision | conditional | n/a | op 39 |
| 159 | Compute B_near1: johnson_bovey_field for nearest ring | ring[0] geometry, JB offset, atom position, unit current | Physics: B-field from nearest ring alone | Vec3 B_near1 | n/a | op 39 |
| 160 | Guard: rings_by_distance.size() >= 2 for B_near2 | rings count | Decision | conditional | n/a | op 39 |
| 161 | Compute B_near2: johnson_bovey_field for 2nd nearest ring | ring[1] geometry, JB offset, atom position, unit current | Physics: B-field from 2nd nearest ring alone | Vec3 B_near2 | n/a | op 39 |
| 162 | Push B_near1 (or Zero if no ring) | B_near1 | none | L1 vector: jb_B_near1 | L1e | op 159 |
| 163 | Push B_near2 (or Zero if < 2 rings) | B_near2 | none | L1 vector: jb_B_near2 | L1e | op 161 |

### L=1 ring normals (pseudovectors, 1e) [2 features] (lines 404-410)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 164 | Guard: rings not empty for normal1 | rings count | Decision | conditional | n/a | op 39 |
| 165 | Get normal1 = nearest ring normal | rings[0].normal | none | Vec3 normal1 | n/a | op 39 |
| 166 | Guard: >= 2 rings for normal2 | rings count | Decision | conditional | n/a | op 39 |
| 167 | Get normal2 = 2nd nearest ring normal | rings[1].normal | none | Vec3 normal2 | n/a | op 39 |
| 168 | Push normal1 (or Zero) | normal1 | none | L1 vector: ring_normal_1 | L1e | op 165 |
| 169 | Push normal2 (or Zero) | normal2 | none | L1 vector: ring_normal_2 | L1e | op 167 |

### L=1 displacement directions (polar vectors, 1o) [2 features] (lines 413-419)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 170 | Guard: rings not empty for dir1 | rings count | Decision | conditional | n/a | op 39 |
| 171 | Compute dir1 = safe_unit(ring1.center - pos) | ring1 center, atom position | Direction FROM atom TO ring center | Vec3 dir1 | n/a | op 39 |
| 172 | Guard: >= 2 rings for dir2 | rings count | Decision | conditional | n/a | op 39 |
| 173 | Compute dir2 = safe_unit(ring2.center - pos) | ring2 center, atom position | Direction FROM atom TO ring center | Vec3 dir2 | n/a | op 39 |
| 174 | Push dir1 (or Zero) | dir1 | none | L1 vector: dir_to_ring_1 | L1o | op 171 |
| 175 | Push dir2 (or Zero) | dir2 | none | L1 vector: dir_to_ring_2 | L1o | op 173 |

### L=2 type-summed anisotropic geometric factors [7 features] (lines 424-425)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 176 | Push sum_G_T2[0] through sum_G_T2[6] | per-type accumulated T2 tensors | none | 7 L2 tensors: jb_sum_G_T2_{PHE,TYR,TRP6,TRP5,HIS,HID,HIE} | L2e | op 89 |

### L=2 nearest ring G_T2 [3 features] (lines 428-448)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 177 | Guard: rings not empty | rings count | Decision | conditional | n/a | op 39 |
| 178 | Compute G for nearest ring and decompose to T2 | ring[0] geometry, JB offset, atom position | Physics: JB geometric factor tensor for nearest ring | array<5> G_T2_near1 | n/a | op 39 |
| 179 | Guard: >= 2 rings | rings count | Decision | conditional | n/a | op 39 |
| 180 | Compute G for 2nd nearest ring and decompose to T2 | ring[1] geometry, JB offset, atom position | Physics: JB geometric factor tensor for 2nd ring | array<5> G_T2_near2 | n/a | op 39 |
| 181 | Push G_T2_near1 (or zeros) | T2 of nearest ring | none | L2 tensor: jb_G_T2_near1 | L2e | op 178 |
| 182 | Push G_T2_near2 (or zeros) | T2 of 2nd ring | none | L2 tensor: jb_G_T2_near2 | L2e | op 180 |
| 183 | Push G_T2_exp_sum | exponentially weighted sum of T2 over all rings | none | L2 tensor: jb_G_T2_exp_sum | L2e | op 92 |

### Return value (line 449)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 184 | Return B_total | accumulated B-field | none | Vec3 (consumed by structural_features) | n/a | op 84 |

**ring_current_features totals**: 47 L0 scalars, 7 L1 vectors, 10 L2 tensors.

---

## bond_anisotropy_features (lines 463-544)

### Backbone peptide bond loop (lines 478-491)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 185 | For each PeptideBond pb: compute mcconnell_bond_factor(pb.C_pos, pb.O_pos, pos) | C position, O position, atom position | Physics: (3cos^2 theta - 1)/R^3 from C=O midpoint | double f_co | n/a | ProteinContext.peptide_bonds |
| 186 | Accumulate co_sum += f_co | f_co, co_sum | none | double co_sum | n/a | op 185 |
| 187 | Compute midpoint mid = 0.5 * (C + O) | pb.C_pos, pb.O_pos | none | Vec3 mid | n/a | none |
| 188 | Compute distance d = (pos - mid).norm() | atom position, midpoint | none | double d | n/a | op 187 |
| 189 | Update nearest C=O: if d < best_co_dist | d, best_co_dist | Decision: is this the closest C=O bond midpoint? | updates best_co_dist, co_nearest, nearest_co_mid | n/a | op 188 |
| 190 | Compute mcconnell_bond_factor(pb.C_pos, pb.N_next_pos, pos) | C position, N_next position, atom position | Physics: McConnell factor for C-N bond | double f_cn | n/a | none |
| 191 | Accumulate cn_sum += f_cn | f_cn, cn_sum | none | double cn_sum | n/a | op 190 |

### Sidechain C=O loop (lines 494-495)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 192 | For each SidechainCO sc: compute mcconnell_bond_factor(sc.C_pos, sc.O_pos, pos) | sidechain C, O positions, atom position | Physics: McConnell factor for sidechain C=O | double | n/a | ProteinContext.sidechain_cos |
| 193 | Accumulate side_co += factor | factor, side_co | none | double side_co | n/a | op 192 |

### Nearest C-N bond search (lines 498-507)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 194 | For each PeptideBond pb: compute C-N midpoint | pb.C_pos, pb.N_next_pos | none | Vec3 mid | n/a | none |
| 195 | Compute distance d = (pos - mid).norm() | atom position, midpoint | none | double d | n/a | op 194 |
| 196 | Update nearest: if d < best_cn_dist | d, best_cn_dist | Decision: closest C-N midpoint? | updates best_cn_dist, nearest_cn_mid | n/a | op 195 |

### Nearest C=O dipolar T2 (lines 510-518)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 197 | Guard: best_co_dist < 1e20 AND best_co_dist > MIN_DISTANCE (0.1) | best_co_dist | Decision: was a C=O found and not too close? | conditional | n/a | op 189 |
| 198 | Compute displacement: d = pos - nearest_co_mid | atom position, nearest C=O midpoint | none | Vec3 d | n/a | op 189 |
| 199 | Compute R = d.norm(), R3 = R^3, R5 = R^5 | d | none | doubles | n/a | op 198 |
| 200 | Compute dipolar kernel: T_ab = 3*d_a*d_b/R5 - delta_ab/R3 | d, R3, R5 | Physics: rank-2 dipolar tensor | Mat3 T2_CO_near | n/a | ops 198-199 |

### Nearest C-N dipolar T2 (lines 521-529)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 201 | Guard: best_cn_dist < 1e20 AND best_cn_dist > MIN_DISTANCE | best_cn_dist | Decision: was a C-N found and not too close? | conditional | n/a | op 196 |
| 202 | Compute displacement and dipolar kernel for nearest C-N | atom position, nearest C-N midpoint | Physics: rank-2 dipolar tensor | Mat3 T2_CN_near | n/a | op 196 |

### L=0 outputs [5 features] (lines 532-536)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 203 | Push co_sum | sum of C=O McConnell factors over all peptide bonds | none | L0 scalar: mcconnell_co_sum | L0 | op 186 |
| 204 | Push co_nearest | McConnell factor of nearest C=O | none | L0 scalar: mcconnell_co_nearest | L0 | op 189 |
| 205 | Push cn_sum | sum of C-N McConnell factors over all peptide bonds | none | L0 scalar: mcconnell_cn_sum | L0 | op 191 |
| 206 | Push side_co | sum of sidechain C=O McConnell factors | none | L0 scalar: mcconnell_side_co_sum | L0 | op 193 |
| 207 | Push best_cn_dist (or sentinel 99.0 via NO_DATA_SENTINEL) | best_cn_dist | Decision: if best_cn_dist >= 1e20, use sentinel | L0 scalar: mcconnell_cn_nearest_dist | L0 | op 196 |

Magic number: 1e20 used as "not found" threshold (compared against initial 1e30). Sentinel: constants::NO_DATA_SENTINEL = 99.0.

### L=1 output [1 feature] (line 539)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 208 | Push safe_unit(nearest_co_mid - pos) | nearest C=O midpoint, atom position | Direction FROM atom TO nearest C=O midpoint | L1 vector: mcconnell_dir_nearest_co | L1o | op 189 |

### L=2 outputs [2 features] (lines 542-543)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 209 | Push decompose_tensor(T2_CO_near).T2 | nearest C=O dipolar tensor | Decompose to T2 components | L2 tensor: peptide_T2_CO_nearest | L2e | op 200 |
| 210 | Push decompose_tensor(T2_CN_near).T2 | nearest C-N dipolar tensor | Decompose to T2 components | L2 tensor: peptide_T2_CN_nearest | L2e | op 202 |

**bond_anisotropy_features totals**: 5 L0, 1 L1, 2 L2.

---

## coulomb_field_features (lines 558-628)

### Electric field and EFG accumulation (lines 575-599)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 211 | Initialize E_total, E_backbone to Zero | none | none | Vec3 accumulators | n/a | none |
| 212 | Initialize EFG_total, EFG_backbone to Zero | none | none | Mat3 accumulators | n/a | none |
| 213 | For each pdb_atom i: read charge q = ctx.charges[i] | ProteinContext.charges[i] | none | double q | n/a | ProteinContext.charges |
| 214 | Skip if |q| < 1e-10 | q | Decision: skip uncharged atoms | skip | n/a | op 213 |
| 215 | Compute displacement: dr = pos - pdb_atoms[i].position | atom position, other atom position | none | Vec3 dr | n/a | none |
| 216 | Compute distance: r = dr.norm() | dr | none | double r | n/a | op 215 |
| 217 | Skip if r < MIN_DISTANCE (0.1) | r | Decision: avoid self-interaction / singularity | skip | n/a | op 216 |
| 218 | Compute r^3 and r^5 | r | none | doubles r3, r5 | n/a | op 216 |
| 219 | Compute E-field contribution: contrib = q * dr / r3 | q, dr, r3 | Physics: Coulomb field from point charge | Vec3 contrib | n/a | ops 213, 215, 218 |
| 220 | Accumulate E_total += contrib | contrib | none | Vec3 E_total | n/a | op 219 |
| 221 | Check backbone flag: ctx.is_backbone[i] | ProteinContext.is_backbone[i] | Decision: is source atom backbone? | conditional accumulation | n/a | none |
| 222 | If backbone: accumulate E_backbone += contrib | contrib | none | Vec3 E_backbone | n/a | ops 219, 221 |
| 223 | Compute EFG tensor: V(a,b) = q * (3*dr_a*dr_b/r5 - delta_ab/r3) | q, dr, r3, r5 | Physics: electric field gradient (dipolar kernel) | Mat3 V | n/a | ops 213, 215, 218 |
| 224 | Accumulate EFG_total += V | V | none | Mat3 EFG_total | n/a | op 223 |
| 225 | If backbone: accumulate EFG_backbone += V | V | none | Mat3 EFG_backbone | n/a | ops 223, 221 |

Magic numbers: 1e-10 (charge threshold), MIN_DISTANCE = 0.1 (singularity cutoff).

### L=0 outputs [4 features] (lines 601-619)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 226 | Compute E_mag = E_total.norm() | E_total | none | double E_mag | n/a | op 220 |
| 227 | Push E_mag | total E-field magnitude | none | L0 scalar: coulomb_E_mag | L0 | op 226 |
| 228 | Guard: site.parent exists AND E_mag > 1e-15 | parent pointer, E_mag | Decision: only compute bond projection for H atoms with non-zero field | conditional | n/a | ops 57, 226 |
| 229 | Compute E_bond_proj = E_total.dot(site.bond_direction) | E_total, bond_direction | Physics: Buckingham E-bond projection (E projected onto parent->H direction) | double E_bond_proj | n/a | ops 59, 220 |
| 230 | Push E_bond_proj (or 0.0 if guard fails) | E_bond_proj | none | L0 scalar: coulomb_E_bond_proj | L0 | op 229 |
| 231 | Guard: rings_by_distance not empty | rings | Decision: need a ring normal to project onto | conditional | n/a | op 39 |
| 232 | Compute E_ring_proj = E_total.dot(nearest ring normal) | E_total, ring[0].normal | Physics: E projected onto nearest ring normal | double E_ring_proj | n/a | ops 39, 220 |
| 233 | Push E_ring_proj (or 0.0) | E_ring_proj | none | L0 scalar: coulomb_E_ring_proj | L0 | op 232 |
| 234 | Compute backbone fraction: if E_mag > 1e-15, E_backbone.norm() / E_mag, else 0 | E_backbone, E_mag | Decision: avoid division by zero (threshold 1e-15) | double | n/a | ops 222, 226 |
| 235 | Push backbone fraction | fraction | none | L0 scalar: coulomb_E_backbone_frac | L0 | op 234 |

Magic number: 1e-15 used as E-field zero threshold.

### L=1 outputs [2 features] (lines 622-623)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 236 | Push E_total | total Coulomb E-field vector | none | L1 vector: coulomb_E_total | L1o | op 220 |
| 237 | Push E_backbone | backbone-only E-field vector | none | L1 vector: coulomb_E_backbone | L1o | op 222 |

### L=2 outputs [2 features] (lines 626-627)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 238 | Push decompose_tensor(EFG_total).T2 | total EFG tensor | Decompose to T2 | L2 tensor: coulomb_EFG_T2_total | L2e | op 224 |
| 239 | Push decompose_tensor(EFG_backbone).T2 | backbone EFG tensor | Decompose to T2 | L2 tensor: coulomb_EFG_T2_backbone | L2e | op 225 |

**coulomb_field_features totals**: 4 L0, 2 L1, 2 L2.

---

## atom_identity_features (lines 635-651)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 240 | Push (element == "H" ? 1.0 : 0.0) | atom.element | Categorical: is hydrogen? | L0 scalar: atom_is_H | L0 | none |
| 241 | Push (element == "C" ? 1.0 : 0.0) | atom.element | Categorical: is carbon? | L0 scalar: atom_is_C | L0 | none |
| 242 | Push (element == "N" ? 1.0 : 0.0) | atom.element | Categorical: is nitrogen? | L0 scalar: atom_is_N | L0 | none |
| 243 | Push (element == "O" ? 1.0 : 0.0) | atom.element | Categorical: is oxygen? | L0 scalar: atom_is_O | L0 | none |
| 244 | Push (element == "S" ? 1.0 : 0.0) | atom.element | Categorical: is sulfur? | L0 scalar: atom_is_S | L0 | none |
| 245 | Push (double)atom.Z | atom.Z (atomic number) | none | L0 scalar: atom_Z | L0 | none |
| 246 | Push (is_backbone ? 1.0 : 0.0) | site.is_backbone | Categorical: op 48 | L0 scalar: atom_is_backbone | L0 | op 48 |
| 247 | Push (is_amide_H ? 1.0 : 0.0) | site.is_amide_H | Categorical: op 49 | L0 scalar: atom_is_amide_H | L0 | op 49 |
| 248 | Push (is_alpha_H ? 1.0 : 0.0) | site.is_alpha_H | Categorical: op 50 | L0 scalar: atom_is_alpha_H | L0 | op 50 |
| 249 | Push (is_methyl ? 1.0 : 0.0) | site.is_methyl | Categorical: op 51 | L0 scalar: atom_is_methyl | L0 | op 51 |

**atom_identity_features totals**: 10 L0, 0 L1, 0 L2.

---

## structural_features (lines 658-705)

### L=0 outputs [10 features] (lines 667-700)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 250 | Push nearest_acceptor_dist | site.nearest_acceptor_dist | none (default 99.0 sentinel) | L0 scalar: struct_acceptor_dist | L0 | op 64 |
| 251 | Push nearest_co_dist | site.nearest_co_dist | none (default 99.0 sentinel) | L0 scalar: struct_co_dist | L0 | op 68 |
| 252 | Compute co_vec = co_O_pos - co_C_pos | site.co_O_pos, site.co_C_pos | none | Vec3 co_vec | n/a | op 70 |
| 253 | Guard: co_vec.norm() > 0.1 AND nearest_co_dist < 90.0 | co_vec norm, nearest_co_dist | Decision: is there a valid C=O nearby? Threshold: 0.1 for vector norm, 90.0 for distance | conditional | n/a | ops 252, 68 |
| 254 | Compute co_mcconnell = mcconnell_bond_factor(co_C_pos, co_O_pos, pos) | site.co_C_pos, site.co_O_pos, atom position | Physics: McConnell factor from nearest C=O | double co_mcconnell | n/a | op 253 |
| 255 | Push co_mcconnell (or 0.0) | co_mcconnell | none | L0 scalar: struct_co_mcconnell | L0 | op 254 |
| 256 | Push (double)n_heavy_within_8A | site.n_heavy_within_8A | none | L0 scalar: struct_n_heavy_8A | L0 | op 73 |
| 257 | Push min_ring_atom_dist | site.min_ring_atom_dist | none (default 99.0 sentinel) | L0 scalar: struct_min_ring_atom_dist | L0 | op 77 |
| 258 | Compute ring_center_dist: if rings empty, 99.0; else rings_by_distance[0].distance | rings_by_distance | Decision: any rings? | double ring_center_dist | n/a | op 39 |
| 259 | Compute ratio: if ring_center_dist > 0.1, min_ring_atom_dist / ring_center_dist, else 1.0 | min_ring_atom_dist, ring_center_dist | Decision: avoid division by zero (threshold 0.1) | double | n/a | ops 257-258 |
| 260 | Push ratio | min_ring_atom_dist / ring_center_dist | none | L0 scalar: struct_min_center_ratio | L0 | op 259 |
| 261 | Guard: site.parent exists AND B_total.norm() > 1e-15 | parent, B_total | Decision: only for H atoms with non-zero B-field | conditional | n/a | ops 57, 184 |
| 262 | Compute B_bond_cos = bond_direction.dot(B_total.normalized()) | bond_direction, B_total | Physics: cosine between parent->H bond and total B-field | double B_bond_cos | n/a | ops 59, 184 |
| 263 | Push B_bond_cos (or 0.0) | B_bond_cos | none | L0 scalar: struct_B_bond_cos | L0 | op 262 |
| 264 | Guard: site.parent exists AND rings not empty | parent, rings_by_distance | Decision: H atom with a nearby ring | conditional | n/a | ops 57, 39 |
| 265 | Compute ch_pi_cos = |bond_direction.dot(nearest ring normal)| | bond_direction, ring[0].normal | Physics: |cos| of angle between C-H bond and ring normal (CH-pi interaction indicator) | double ch_pi_cos | n/a | ops 59, 39 |
| 266 | Push ch_pi_cos (or 0.0) | ch_pi_cos | none | L0 scalar: struct_ch_pi_cos | L0 | op 265 |
| 267 | Push (is_on_aromatic_residue ? 1.0 : 0.0) | site.is_on_aromatic_residue | Categorical: op 52 | L0 scalar: struct_is_aromatic_res | L0 | op 52 |
| 268 | Push (parent_is_sp2 ? 1.0 : 0.0) | site.parent_is_sp2 | Categorical: op 60 | L0 scalar: struct_parent_sp2 | L0 | op 60 |

Magic numbers: 0.1 (C=O vector norm threshold), 90.0 (C=O distance threshold), 0.1 (ring center distance threshold), 1e-15 (B-field zero threshold).

### L=1 outputs [2 features] (lines 703-704)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 269 | Push dir_to_nearest_acceptor | site.dir_to_nearest_acceptor | none (Zero if none found) | L1 vector: struct_dir_acceptor | L1o | op 65 |
| 270 | Push bond_direction | site.bond_direction | none (Zero for heavy atoms) | L1 vector: struct_bond_direction | L1o | op 59 |

**structural_features totals**: 10 L0, 2 L1, 0 L2.

---

## residue_encoding (lines 712-726)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 271 | Read residue_name from atom | site.atom->residue_name | none | string rn | n/a | none |
| 272 | Normalize histidine variants: if rn in {HID, HIE, HIP}, set rn = "HIS" | rn | Decision: collapse all histidine protonation states to HIS for one-hot | string rn (modified) | n/a | op 271 |
| 273 | For each of 20 standard amino acids: push (rn == aa[i]) ? 1.0 : 0.0 | rn, aa20 list | Categorical: 20-class one-hot. Order: ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL | 20 L0 scalars: res_{AA} | L0 | op 272 |

Note: Non-standard residues (e.g., selenomethionine MSE) will produce an all-zero one-hot vector.

**residue_encoding totals**: 20 L0, 0 L1, 0 L2.

---

## haigh_mallion_features (lines 741-788)

### Per-ring accumulation (lines 755-773)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 274 | Initialize 7-element arrays hm_T0_by_type and hm_T2_by_type to zero | NT=7 | none | arrays | n/a | none |
| 275 | Initialize hm_B_total = Zero, hm_sd_nearest = default | none | none | accumulators | n/a | none |
| 276 | For each ring neighbor rn: get ring and compute type index | rn.index, ring.residue_name, ring.ring_type | Categorical: type_index mapping (same as op 82) | int tidx | n/a | AtomSite.rings_by_distance |
| 277 | Compute HM geometric factor tensor: haigh_mallion_tensor_factor(ring, pos) | ring geometry, atom position | Physics: Haigh-Mallion surface integral (7-point Gaussian quadrature on fan triangulation) | Mat3 G_HM | n/a | none |
| 278 | Decompose: sd = decompose_tensor(G_HM) | G_HM | Spherical decomposition | SphericalDecomposition sd | n/a | op 277 |
| 279 | Type-index bounds check and accumulate hm_T0_by_type[tidx] += sd.T0 | tidx, sd.T0 | Decision: valid type index? | per-type scalar | n/a | ops 276, 278 |
| 280 | Accumulate hm_T2_by_type[tidx][m] += sd.T2[m] for m in 0..4 | tidx, sd.T2 | none | per-type tensor | n/a | ops 276, 278 |
| 281 | Compute HM B-field: haigh_mallion_B_field(ring, pos) and accumulate | ring geometry, atom position | Physics: HM B-field from each ring | Vec3 hm_B_total | n/a | none |
| 282 | Check if this is the first (nearest) ring: &rn == &rings_by_distance[0] | pointer comparison | Decision: is this the nearest ring? (pointer address check) | conditional | n/a | none |
| 283 | If nearest: store hm_sd_nearest = sd | sd | none | SphericalDecomposition | n/a | op 278 |

### L=0 outputs [9 features] (lines 776-779)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 284 | Push hm_T0_by_type[0..6] | per-type HM isotropic factors | none | 7 L0 scalars: hm_sum_T0_{PHE..HIE} | L0 | op 279 |
| 285 | Push hm_sd_nearest.T0 | HM isotropic factor of nearest ring | none | L0 scalar: hm_T0_nearest | L0 | op 283 |
| 286 | Push hm_B_total.norm() | total HM B-field magnitude | none | L0 scalar: hm_B_mag | L0 | op 281 |

### L=1 output [1 feature] (line 782)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 287 | Push hm_B_total | total HM B-field vector | none | L1 vector: hm_B_total | L1e | op 281 |

### L=2 outputs [8 features] (lines 785-787)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 288 | Push hm_T2_by_type[0..6] | per-type HM T2 tensors | none | 7 L2 tensors: hm_sum_T2_{PHE..HIE} | L2e | op 280 |
| 289 | Push hm_sd_nearest.T2 | nearest ring HM T2 | none | L2 tensor: hm_T2_nearest | L2e | op 283 |

**haigh_mallion_features totals**: 9 L0, 1 L1, 8 L2.

---

## graph_features (lines 803-837)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 290 | Get atom index: idx = site.atom->index | atom.index | none | int idx | n/a | none |
| 291 | Bounds check: idx < 0 OR idx >= graph.n_atoms() | idx, graph size | Decision: is atom in graph? Fallback: push 10 zeros | conditional | n/a | op 290, requires graph != nullptr |
| 292 | Get graph environment: graph.environment(idx) | idx | none | GraphEnvironment& env | n/a | op 291 |
| 293 | Push dist_to_ring_atom (or 99.0 sentinel if -1) | env.dist_to_ring_atom | Decision: if dist < 0, use sentinel 99.0 | L0 scalar: graph_dist_ring | L0 | op 292 |
| 294 | Push dist_to_nitrogen (or 99.0 sentinel if -1) | env.dist_to_nitrogen | Decision: if dist < 0, use sentinel 99.0 | L0 scalar: graph_dist_N | L0 | op 292 |
| 295 | Push dist_to_oxygen (or 99.0 sentinel if -1) | env.dist_to_oxygen | Decision: if dist < 0, use sentinel 99.0 | L0 scalar: graph_dist_O | L0 | op 292 |
| 296 | Push eneg_sum_1 | env.eneg_sum_1 | none | L0 scalar: graph_eneg_sum_1 | L0 | op 292 |
| 297 | Push eneg_sum_2 | env.eneg_sum_2 | none | L0 scalar: graph_eneg_sum_2 | L0 | op 292 |
| 298 | Push (hybridization == sp ? 1.0 : 0.0) | env.hybridization | Categorical: sp? | L0 scalar: graph_is_sp | L0 | op 292 |
| 299 | Push (hybridization == sp2 ? 1.0 : 0.0) | env.hybridization | Categorical: sp2? | L0 scalar: graph_is_sp2 | L0 | op 292 |
| 300 | Push (hybridization == sp3 ? 1.0 : 0.0) | env.hybridization | Categorical: sp3? | L0 scalar: graph_is_sp3 | L0 | op 292 |
| 301 | Push (double)n_pi_bonds_3 | env.n_pi_bonds_3 | none | L0 scalar: graph_n_pi_bonds_3 | L0 | op 292 |
| 302 | Push (is_conjugated ? 1.0 : 0.0) | env.is_conjugated | Categorical: conjugated bond chain? | L0 scalar: graph_is_conjugated | L0 | op 292 |

**graph_features totals**: 10 L0, 0 L1, 0 L2.

---

## full_bond_anisotropy_features (lines 853-885)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 303 | Get atom index: idx = site.atom->index | atom.index | none | int idx | n/a | none |
| 304 | Bounds check: idx < 0 OR idx >= graph.n_atoms() | idx, graph size | Decision: is atom in graph? Fallback: 3 L0 zeros, 1 L1 zero, 3 L2 zeros | conditional | n/a | requires graph != nullptr |
| 305 | Call compute_bond_anisotropy(site.atom->position, idx, graph, pdb_atoms) | atom position, atom index, graph, all atoms | Physics: McConnell tensor for all bonds, grouped by category, excluding self-bonds | BondAnisotropyResult | n/a | op 303 |
| 306 | Push result.f_backbone | backbone McConnell factor sum | none | L0 scalar: bond_aniso_f_backbone | L0 | op 305 |
| 307 | Push result.f_sidechain | sidechain McConnell factor sum | none | L0 scalar: bond_aniso_f_sidechain | L0 | op 305 |
| 308 | Push result.f_aromatic | aromatic ring bond McConnell factor sum | none | L0 scalar: bond_aniso_f_aromatic | L0 | op 305 |
| 309 | Push result.dir_nearest_midpoint | direction to nearest bond midpoint (any category) | none | L1 vector: bond_aniso_dir_nearest | L1o | op 305 |
| 310 | Push result.T2_backbone | backbone dipolar kernel T2 | none | L2 tensor: bond_aniso_T2_backbone | L2e | op 305 |
| 311 | Push result.T2_sidechain | sidechain dipolar kernel T2 | none | L2 tensor: bond_aniso_T2_sidechain | L2e | op 305 |
| 312 | Push result.T2_aromatic | aromatic dipolar kernel T2 | none | L2 tensor: bond_aniso_T2_aromatic | L2e | op 305 |

**full_bond_anisotropy_features totals**: 3 L0, 1 L1, 3 L2.

---

## apbs_field_features (lines 892-928)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 313 | Copy apbs_E to local E | apbs_E (const Vec3&) | none | Vec3 E | n/a | requires apbs_efield != nullptr |
| 314 | For each component i in 0..2: sanitize NaN/Inf to 0.0 | E(i) | Decision: if NaN or Inf, clamp to zero | Vec3 E (sanitized) | n/a | op 313 |
| 315 | Compute E_mag = E.norm() | E | none | double E_mag | n/a | op 314 |
| 316 | Push E_mag | APBS E-field magnitude | none | L0 scalar: apbs_E_mag | L0 | op 315 |
| 317 | Guard: site.parent exists AND E_mag > 1e-15 | parent, E_mag | Decision: H atom with non-zero field | conditional | n/a | ops 57, 315 |
| 318 | Compute E_bond = E.dot(site.bond_direction) | E, bond_direction | Physics: APBS E projected onto parent->H bond | double E_bond | n/a | ops 59, 314 |
| 319 | Push E_bond (or 0.0) | E_bond | none | L0 scalar: apbs_E_bond_proj | L0 | op 318 |
| 320 | Push E | APBS E-field vector | none | L1 vector: apbs_E_total | L1o | op 314 |
| 321 | Guard: E_mag > 1e-15 for normalized outer product | E_mag | Decision: avoid division by zero | conditional | n/a | op 315 |
| 322 | Compute normalized outer product: EE(a,b) = E(a)*E(b) / E_mag^2 | E, E_mag | Physics: unit-scale directional tensor (EFG proxy) | Mat3 EE | n/a | ops 314, 315 |
| 323 | Push decompose_tensor(EE).T2 | EE tensor | Decompose to T2 | L2 tensor: apbs_EFG_T2 | L2e | op 322 |

Magic number: 1e-15 (E-field zero threshold).

**apbs_field_features totals**: 2 L0, 1 L1, 1 L2.

---

## hbond_features (lines 935-1019)

### Helper functions (lines 937-947)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 324 | is_hbond_donor(name): check if name in {N, ND2, NE2, NZ, NH1, NH2} | atom_name | Categorical: H-bond donor nitrogens | bool | n/a | none |
| 325 | is_hbond_acceptor(name): check if name in {O, OD1, OD2, OE1, OE2, OG, OG1, OH} | atom_name | Categorical: H-bond acceptor oxygens | bool | n/a | none |

### Main method (lines 950-1019)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 326 | Classify current atom: atom_is_donor = is_hbond_donor(atom_name) | site.atom->atom_name | Categorical: is current atom a potential H-bond donor? | bool | n/a | op 324 |
| 327 | Classify current atom: atom_is_acceptor = is_hbond_acceptor(atom_name) | site.atom->atom_name | Categorical: is current atom a potential H-bond acceptor? | bool | n/a | op 325 |
| 328 | Initialize best_dist = 99.0 (NO_DATA_SENTINEL) | sentinel constant | none | double | n/a | none |
| 329 | For each pdb_atom pa: skip if same residue AND same chain | pa.residue_number, pa.chain_id vs site.atom | Decision: exclude same-residue atoms | skip | n/a | none |
| 330 | Check partner compatibility: if atom_is_donor AND pa is acceptor, OR if atom_is_acceptor AND pa is donor | donor/acceptor flags | Decision: is this a valid D...A or A...D pair? | bool partner_ok | n/a | ops 326-327, 324-325 |
| 331 | Skip if !partner_ok | partner_ok | Decision | skip | n/a | op 330 |
| 332 | Skip sequential backbone N...O: if atom is "N" and pa is "O" and |residue_number difference| < 2 | atom_name, pa.atom_name, residue_numbers | Decision: exclude i,i+1 backbone N-O pairs (not real H-bonds) | skip | n/a | none |
| 333 | Skip sequential backbone O...N: if atom is "O" and pa is "N" and |residue_number difference| < 2 | atom_name, pa.atom_name, residue_numbers | Decision: same exclusion reversed | skip | n/a | none |
| 334 | Compute distance: d = (pa.position - pos).norm() | pa.position, atom position | none | double d | n/a | none |
| 335 | If d < 3.5: increment count_within_3_5 | d | Threshold: 3.5 Angstroms for H-bond counting | int counter | n/a | op 334 |
| 336 | Update nearest: if d < best_dist AND d > MIN_DISTANCE (0.1) | d, best_dist | Decision: closest valid H-bond partner | updates best_dist, best_dir, best_is_backbone | n/a | op 334 |
| 337 | Compute direction: best_dir = safe_unit(pa.position - pos) | pa.position, atom position | Direction FROM atom TO partner | Vec3 best_dir | n/a | op 336 |
| 338 | Classify H-bond backbone: both atom and partner have backbone atom names? | is_backbone_atom check for both | Decision: is this a backbone-backbone H-bond? | bool best_is_backbone | n/a | ops 1-7 |

Magic numbers: 3.5 A (H-bond counting radius), MIN_DISTANCE = 0.1 A (singularity cutoff).

### L=0 outputs [6 features] (lines 994-1000)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 339 | Push best_dist | nearest H-bond partner distance (or 99.0 sentinel) | none | L0 scalar: hbond_dist | L0 | op 336 |
| 340 | Compute inv_d3: if best_dist < 50.0, 1/(best_dist^3), else 0 | best_dist | Decision: threshold 50.0 for "found" | L0 scalar: hbond_inv_d3 | L0 | op 339 |
| 341 | Push (atom_is_donor ? 1.0 : 0.0) | atom_is_donor | Categorical | L0 scalar: hbond_is_donor | L0 | op 326 |
| 342 | Push (atom_is_acceptor ? 1.0 : 0.0) | atom_is_acceptor | Categorical | L0 scalar: hbond_is_acceptor | L0 | op 327 |
| 343 | Push (best_is_backbone ? 1.0 : 0.0) | best_is_backbone | Categorical | L0 scalar: hbond_is_backbone | L0 | op 338 |
| 344 | Push (double)count_within_3_5 | count | none | L0 scalar: hbond_count_3.5A | L0 | op 335 |

Magic number: 50.0 used as threshold for "H-bond partner exists" in inv_d3 computation.

### L=1 output [1 feature] (line 1003)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 345 | Push best_dir | direction to nearest H-bond partner (or Zero) | none | L1 vector: hbond_dir | L1o | op 337 |

### L=2 output [1 feature] (lines 1006-1018)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 346 | Guard: best_dist < 50.0 AND best_dist > MIN_DISTANCE | best_dist | Decision: valid H-bond distance? | conditional | n/a | op 336 |
| 347 | Reconstruct unnormalized displacement: d = best_dir * best_dist | best_dir, best_dist | none | Vec3 d | n/a | ops 336-337 |
| 348 | Compute dipolar kernel: T(a,b) = 3*d(a)*d(b)/d^5 - delta_ab/d^3 | d, best_dist | Physics: rank-2 dipolar tensor | Mat3 T | n/a | op 347 |
| 349 | Push decompose_tensor(T).T2 (or 5 zeros if guard fails) | T | Decompose to T2 | L2 tensor: hbond_T2 | L2e | op 348 |

**hbond_features totals**: 6 L0, 1 L1, 1 L2.

---

## dispersion_tensor_features (lines 1037-1080)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 350 | Define R_CUT = 5.0, R_MIN = 1.5 | hardcoded constants | none | doubles | n/a | none |
| 351 | Initialize T_disp = Mat3::Zero(), disp_scalar = 0, n_contacts = 0 | none | none | accumulators | n/a | none |
| 352 | For each ring neighbor rn: early exit if rn.distance > R_CUT + 3.0 = 8.0 | rn.distance | Decision: ring center too far for any atom to be within R_CUT (heuristic: ring radius ~ 3 A) | break | n/a | AtomSite.rings_by_distance |
| 353 | For each ring vertex rv: compute displacement d = pos - rv | atom position, vertex position | none | Vec3 d | n/a | none |
| 354 | Compute R = d.norm() | d | none | double R | n/a | op 353 |
| 355 | Skip if R < R_MIN (1.5) or R > R_CUT (5.0) | R | Decision: distance filter (exclude covalent-range and far contacts) | skip | n/a | op 354 |
| 356 | Compute R^2, R^3, R^6, R^8 | R | none | doubles | n/a | op 354 |
| 357 | Accumulate London-weighted dipolar tensor: T_disp(a,b) += 3*d_a*d_b/R^8 - delta_ab/R^6 | d, R^6, R^8 | Physics: dispersion-scale (1/R^6, 1/R^8) dipolar kernel | Mat3 T_disp | n/a | ops 353, 356 |
| 358 | Accumulate disp_scalar += 1/R^6 | R^6 | Physics: London dispersion scalar | double | n/a | op 356 |
| 359 | Increment n_contacts | counter | none | int | n/a | op 355 |

Magic numbers: R_CUT = 5.0 A, R_MIN = 1.5 A, R_CUT + 3.0 = 8.0 A (early exit heuristic).

### L=0 outputs [2 features] (lines 1075-1076)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 360 | Push disp_scalar | sum of 1/R^6 over ring atom contacts | none | L0 scalar: disp_sum_inv_r6 | L0 | op 358 |
| 361 | Push (double)n_contacts | number of ring atom contacts within range | none | L0 scalar: disp_n_contacts | L0 | op 359 |

### L=2 output [1 feature] (line 1079)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 362 | Push decompose_tensor(T_disp).T2 | accumulated dispersion tensor | Decompose to T2 | L2 tensor: disp_T2_contact | L2e | op 357 |

**dispersion_tensor_features totals**: 2 L0, 0 L1, 1 L2.

---

## aromatic_field_features (lines 1092-1155)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 363 | Initialize E_arom, EFG_arom to Zero, n_sidechain = 0 | none | none | accumulators | n/a | none |
| 364 | Define backbone_plus_cb set: {N, CA, C, O, H, HA, HA2, HA3, CB, HB, HB2, HB3} | hardcoded set (static) | none | string set | n/a | none |

Note: This backbone set is BROADER than is_backbone_atom (ops 1-7): includes HA2, HA3, CB, HB, HB2, HB3. This defines what "sidechain beyond CB" means.

| 365 | For each pdb_atom pa: check is_aromatic_res | pa.residue_name | Categorical: is residue PHE, TYR, TRP, HIS, HID, HIE, or HIP? (hardcoded inline, NOT calling is_aromatic_residue) | bool is_aromatic_res | n/a | none |

Note: This duplicates the is_aromatic_residue() logic from PDBParser.h but as inline string comparisons. The set is identical: {PHE, TYR, TRP, HIS, HID, HIE, HIP}.

| 366 | Skip if !is_aromatic_res | is_aromatic_res | Decision | skip | n/a | op 365 |
| 367 | Check backbone_plus_cb membership: if pa.atom_name in the set, skip | pa.atom_name | Decision: exclude backbone and CB atoms (keep only sidechain beyond CB) | skip | n/a | op 364 |
| 368 | Read charge: q = ctx.charges[i] | ProteinContext.charges[i] | none | double q | n/a | none |
| 369 | Compute displacement: dr = pos - pa.position | atom position, pa position | none | Vec3 dr | n/a | none |
| 370 | Compute distance: r = dr.norm() | dr | none | double r | n/a | op 369 |
| 371 | Skip if r < MIN_DISTANCE (0.1) | r | Decision: singularity guard | skip | n/a | op 370 |
| 372 | Compute r^3, r^5 | r | none | doubles | n/a | op 370 |
| 373 | Accumulate E_arom += q * dr / r^3 | q, dr, r^3 | Physics: Coulomb E-field from aromatic sidechain atom | Vec3 | n/a | ops 368-372 |
| 374 | Accumulate EFG_arom(a,b) += q * (3*dr_a*dr_b/r^5 - delta_ab/r^3) | q, dr, r^3, r^5 | Physics: EFG tensor from aromatic sidechain atom | Mat3 | n/a | ops 368-372 |
| 375 | Increment n_sidechain | counter | none | int | n/a | op 367 |
| 376 | Compute E_mag = E_arom.norm() | E_arom | none | double | n/a | op 373 |
| 377 | Guard: site.parent AND E_mag > 1e-15 | parent, E_mag | Decision: H atom with non-zero aromatic field | conditional | n/a | ops 57, 376 |
| 378 | Compute E_bond_proj = E_arom.dot(site.bond_direction) | E_arom, bond_direction | Physics: aromatic field projected onto bond | double | n/a | ops 59, 373 |

### L=0 outputs [3 features] (lines 1146-1148)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 379 | Push E_mag | aromatic sidechain E-field magnitude | none | L0 scalar: arom_E_mag | L0 | op 376 |
| 380 | Push E_bond_proj (or 0.0) | bond projection | none | L0 scalar: arom_E_bond_proj | L0 | op 378 |
| 381 | Push (double)n_sidechain | count of aromatic sidechain atoms contributing | none | L0 scalar: arom_n_sidechain_atoms | L0 | op 375 |

### L=1 output [1 feature] (line 1151)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 382 | Push E_arom | aromatic sidechain E-field vector | none | L1 vector: arom_E_total | L1o | op 373 |

### L=2 output [1 feature] (line 1154)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 383 | Push decompose_tensor(EFG_arom).T2 | aromatic sidechain EFG tensor | Decompose to T2 | L2 tensor: arom_EFG_T2 | L2e | op 374 |

**aromatic_field_features totals**: 3 L0, 1 L1, 1 L2.

---

## per_type_nearest_features (lines 1166-1206)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 384 | Define N_COLLAPSED = 4 collapsed types: PHE=0, TYR=1, TRP=2, HIS-family=3 | hardcoded | none | int constant | n/a | none |
| 385 | Initialize nearest_dist array to 1e30 | 4 elements | none | array | n/a | none |
| 386 | Initialize nearest_T2 array to zeros | 4 x 5 elements | none | array | n/a | none |
| 387 | For each ring neighbor rn: determine collapsed type ct | ring.residue_name | Categorical: PHE->0, TYR->1, TRP->2, everything else (HIS/HID/HIE/HIP)->3 | int ct | n/a | none |

Note: This collapses the 7-type system to 4 types. The collapse is done by string comparison on residue_name, NOT by using ring_type_idx. TRP 5-ring and 6-ring are both ct=2. HIS/HID/HIE/HIP are all ct=3. Any unrecognized residue also goes to ct=3 (the else branch).

| 388 | Update nearest: if rn.distance < nearest_dist[ct] | rn.distance, nearest_dist[ct] | Decision: is this the closest ring of this collapsed type? | update nearest_dist and nearest_T2 | n/a | op 387 |
| 389 | If updated: compute full ring_type_idx (7-type) for JB offset lookup | ring.residue_name, ring.ring_type | Maps to JB offset via ring_type_idx (7 types) | int tidx | n/a | op 388 |
| 390 | Compute JB geometric factor and decompose to T2 | ring geometry, JB offset, atom position | Physics: JB geometric factor for this specific nearest ring | array<5> T2 | n/a | op 389 |

### L=0 outputs [4 features] (lines 1199-1201)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 391 | For each collapsed type t in 0..3: push nearest_dist[t] (or sentinel 99.0 if > 1e20) | nearest_dist[t] | Decision: if no ring of this type found, use NO_DATA_SENTINEL | 4 L0 scalars: near_dist_{PHE,TYR,TRP,HIS} | L0 | op 388 |

### L=2 outputs [4 features] (lines 1204-1205)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 392 | For each collapsed type t in 0..3: push nearest_T2[t] | T2 of nearest ring per type | none | 4 L2 tensors: jb_G_T2_near_{PHE,TYR,TRP,HIS} | L2e | op 390 |

**per_type_nearest_features totals**: 4 L0, 0 L1, 4 L2.

---

## quadrupole_field_features (lines 1233-1292)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 393 | Initialize quad_scalar_sum = 0, quad_EFG_total = Zero | none | none | accumulators | n/a | none |
| 394 | For each ring neighbor rn: early exit if rn.distance > 15.0 | rn.distance | Decision: 15.0 Angstrom cutoff for quadrupole field | break | n/a | AtomSite.rings_by_distance |
| 395 | Compute displacement: d = pos - ring.center | atom position, ring center | none | Vec3 d | n/a | none |
| 396 | Compute distance: r = d.norm() | d | none | double r | n/a | op 395 |
| 397 | Skip if r < MIN_DISTANCE (0.1) | r | Decision: singularity guard | skip | n/a | op 396 |
| 398 | Get ring normal: n = ring.normal | ring.normal | none | Vec3 n | n/a | none |
| 399 | Compute cos_theta = d.dot(n) / r | d, n, r | none | double cos_theta | n/a | ops 395-396, 398 |
| 400 | Compute r^4 | r | none | double r4 | n/a | op 396 |
| 401 | Compute quadrupole scalar: (3*cos_theta^2 - 1) / r^4 | cos_theta, r4 | Physics: angular x distance factor for axial quadrupole field | double quad_factor | n/a | ops 399-400 |
| 402 | Accumulate quad_scalar_sum += quad_factor | quad_factor | none | double | n/a | op 401 |
| 403 | Compute r^5, r^7 | r | none | doubles | n/a | op 396 |
| 404 | Compute dn = d.dot(n) (= r * cos_theta) | d, n | none | double dn | n/a | ops 395, 398 |
| 405 | Compute quadrupole EFG tensor components for each (a,b): t1 = 15*dn^2*d_a*d_b/r^7, t2 = -3*dn*(n_a*d_b + n_b*d_a)/r^5, t3 = -3*d_a*d_b/r^5, t4 = delta_ab * (3*dn^2/r^5 - 1/r^3) | d, n, dn, r^5, r^7 | Physics: second derivative of quadrupole potential | Mat3 contribution | n/a | ops 395, 398, 403-404 |
| 406 | Accumulate quad_EFG_total += contribution | contribution | none | Mat3 | n/a | op 405 |

Magic number: 15.0 A cutoff for quadrupole contributions.

### L=0 output [1 feature] (line 1288)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 407 | Push quad_scalar_sum | sum of quadrupole angular factors over all rings | none | L0 scalar: quad_scalar_sum | L0 | op 402 |

### L=2 output [1 feature] (line 1291)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 408 | Push decompose_tensor(quad_EFG_total).T2 | total quadrupole EFG tensor | Decompose to T2 | L2 tensor: quad_EFG_T2 | L2e | op 406 |

**quadrupole_field_features totals**: 1 L0, 0 L1, 1 L2.

---

## ring_chi_aniso_features (lines 1322-1369)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 409 | Initialize chi_scalar_sum = 0, chi_nearest = 0, best_dist = 1e30, chi_tensor_total = Zero | none | none | accumulators | n/a | none |
| 410 | For each ring neighbor rn: early exit if rn.distance > 15.0 | rn.distance | Decision: 15.0 Angstrom cutoff | break | n/a | AtomSite.rings_by_distance |
| 411 | Compute displacement: d = pos - ring.center | atom position, ring center | none | Vec3 d | n/a | none |
| 412 | Compute distance: r = d.norm() | d | none | double r | n/a | op 411 |
| 413 | Skip if r < MIN_DISTANCE (0.1) | r | Decision: singularity guard | skip | n/a | op 412 |
| 414 | Get ring normal: n = ring.normal | ring.normal | none | Vec3 n | n/a | none |
| 415 | Compute cos_theta = d.dot(n) / r | d, n, r | none | double cos_theta | n/a | ops 411-412, 414 |
| 416 | Compute r^3 | r | none | double r3 | n/a | op 412 |
| 417 | Compute McConnell factor: (3*cos_theta^2 - 1) / r^3 | cos_theta, r3 | Physics: whole-ring McConnell with ring center as reference | double mcconnell | n/a | ops 415-416 |
| 418 | Accumulate chi_scalar_sum += mcconnell | mcconnell | none | double | n/a | op 417 |
| 419 | Update nearest: if r < best_dist, set chi_nearest = mcconnell | r, best_dist | Decision: is this the nearest ring? | double chi_nearest | n/a | ops 412, 417 |
| 420 | Compute r^5 | r | none | double r5 | n/a | op 412 |
| 421 | Accumulate ring-centre dipolar tensor: chi_tensor_total(a,b) += 3*d_a*d_b/r^5 - delta_ab/r^3 | d, r3, r5 | Physics: dipolar kernel from ring centre | Mat3 | n/a | ops 411, 416, 420 |

Magic number: 15.0 A cutoff (same as quadrupole).

### L=0 outputs [2 features] (lines 1364-1365)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 422 | Push chi_scalar_sum | sum of ring-McConnell over all rings | none | L0 scalar: chi_aniso_sum | L0 | op 418 |
| 423 | Push chi_nearest | McConnell factor for nearest ring | none | L0 scalar: chi_aniso_nearest | L0 | op 419 |

### L=2 output [1 feature] (line 1368)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 424 | Push decompose_tensor(chi_tensor_total).T2 | total ring-centre dipolar tensor | Decompose to T2 | L2 tensor: chi_aniso_T2 | L2e | op 421 |

**ring_chi_aniso_features totals**: 2 L0, 0 L1, 1 L2.

---

## extract_atom (lines 1376-1456) -- Top-level orchestrator

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 425 | Create IrrepFeatures out | none | none | empty container | n/a | none |
| 426 | Build AtomSite: AtomSite::build(atom, rings, pdb_atoms, ctx) | all inputs | Computes all of ops 36-77 | AtomSite site | n/a | none |
| 427 | Call ring_current_features(site, rings, out) -> B_total | site, rings | ops 78-184 | B_total, L0/L1/L2 appended | all | op 426 |
| 428 | Call bond_anisotropy_features(site, ctx, out) | site, ctx | ops 185-210 | L0/L1/L2 appended | all | op 426 |
| 429 | Call coulomb_field_features(site, pdb_atoms, ctx, rings, out) | site, pdb_atoms, ctx, rings | ops 211-239 | L0/L1/L2 appended | all | op 426 |
| 430 | Call atom_identity_features(site, out) | site | ops 240-249 | L0 appended | L0 | op 426 |
| 431 | Call structural_features(site, rings, B_total, out) | site, rings, B_total | ops 250-270 | L0/L1 appended | L0, L1 | ops 426-427 |
| 432 | Call residue_encoding(site, out) | site | ops 271-273 | L0 appended | L0 | op 426 |
| 433 | Call haigh_mallion_features(site, rings, out) | site, rings | ops 274-289 | L0/L1/L2 appended | all | op 426 |
| 434 | Guard: graph != nullptr | graph pointer | Decision: are molecular graph features available? | conditional | n/a | none |
| 435 | If graph: call graph_features(site, *graph, out) | site, graph | ops 290-302 | L0 appended | L0 | ops 426, 434 |
| 436 | If graph: call full_bond_anisotropy_features(site, *graph, pdb_atoms, out) | site, graph, pdb_atoms | ops 303-312 | L0/L1/L2 appended | all | ops 426, 434 |
| 437 | If NO graph: push 10 L0 zeros (graph_features placeholder) | none | Decision: graph absent, zero-fill | 10 zeros | L0 | op 434 |
| 438 | If NO graph: push 3 L0 zeros + 1 L1 zero + 3 L2 zeros (full_bond_aniso placeholder) | none | Decision: graph absent, zero-fill | zeros | all | op 434 |
| 439 | Guard: apbs_efield != nullptr | efield pointer | Decision: is APBS data available? | conditional | n/a | none |
| 440 | If APBS: call apbs_field_features(site, *apbs_efield, out) | site, efield | ops 313-323 | L0/L1/L2 appended | all | ops 426, 439 |
| 441 | If NO APBS: push 2 L0 zeros + 1 L1 zero + 1 L2 zero | none | Decision: APBS absent, zero-fill | zeros | all | op 439 |
| 442 | Call hbond_features(site, pdb_atoms, out) | site, pdb_atoms | ops 324-349 | L0/L1/L2 appended | all | op 426 |
| 443 | Call dispersion_tensor_features(site, rings, out) | site, rings | ops 350-362 | L0/L2 appended | L0, L2 | op 426 |
| 444 | Call aromatic_field_features(site, pdb_atoms, ctx, rings, out) | site, pdb_atoms, ctx, rings | ops 363-383 | L0/L1/L2 appended | all | op 426 |
| 445 | Call per_type_nearest_features(site, rings, out) | site, rings | ops 384-392 | L0/L2 appended | L0, L2 | op 426 |
| 446 | Call quadrupole_field_features(site, rings, out) | site, rings | ops 393-408 | L0/L2 appended | L0, L2 | op 426 |
| 447 | Call ring_chi_aniso_features(site, rings, out) | site, rings | ops 409-424 | L0/L2 appended | L0, L2 | op 426 |
| 448 | Assert: out.L0.size() == L0_names().size() | out.L0, L0_names | Validation: feature count must match registry | assert | n/a | all above |
| 449 | Assert: out.L1.size() == L1_names().size() | out.L1, L1_names | Validation: feature count must match registry | assert | n/a | all above |
| 450 | Assert: out.L2.size() == L2_names().size() | out.L2, L2_names | Validation: feature count must match registry | assert | n/a | all above |
| 451 | Set out.element = atom.element | atom.element | none | metadata | n/a | none |
| 452 | Set out.protein_id = protein_id | protein_id parameter | none | metadata | n/a | none |
| 453 | Set out.atom_name = atom.atom_name | atom.atom_name | none | metadata | n/a | none |
| 454 | Set out.residue_name = atom.residue_name | atom.residue_name | none | metadata | n/a | none |
| 455 | Set out.residue_number = atom.residue_number | atom.residue_number | none | metadata | n/a | none |
| 456 | Set out.position = atom.position | atom.position | none | metadata | n/a | none |
| 457 | Set out.target_T0 = delta.T0 | DeltaTensor.T0 | none | target | n/a | none |
| 458 | Set out.target_T2 = delta.T2 | DeltaTensor.T2 | none | target | n/a | none |

---

## extract_atom_no_target (lines 1458-1477)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 459 | Create dummy DeltaTensor with zero matrix and ala_index = -1 | atom fields | none | DeltaTensor | n/a | none |
| 460 | Copy atom metadata into dummy (element, atom_name, residue_name, residue_number, position) | atom fields | none | DeltaTensor fields | n/a | none |
| 461 | Delegate to extract_atom(atom, dummy, rings, pdb_atoms, ctx, "", graph, apbs_efield) | all inputs + dummy | Produces features with zero targets | IrrepFeatures | all | ops 425-458 |

---

## extract_all (lines 1479-1531)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 462 | For each DeltaTensor d: distance filter | d.dist_nearest_ring, min_distance, max_distance | Decision: skip atoms outside [min_distance, max_distance] range | skip | n/a | none |
| 463 | Construct AtomInfo `matched` from DeltaTensor fields | d.element, d.position, etc. | none | AtomInfo | n/a | none |
| 464 | Compute Z: element_to_z(d.element) | d.element | Categorical: element string to atomic number | int Z | n/a | none |
| 465 | Check: matched.atom_name.empty() | atom_name | Decision: if atom_name missing, do spatial lookup | conditional | n/a | none |
| 466 | If atom_name empty: iterate pdb_atoms, find closest by position with matching element, threshold 1.0 A | d.position, pa.position, d.element | Decision: spatial matching with 1.0 A threshold and element filter | AtomInfo matched | n/a | op 465 |
| 467 | Guard: apbs_efield != nullptr AND matched.index valid | efield pointer, atom index | Decision: is per-atom APBS data available for this atom? | conditional | n/a | none |
| 468 | Guard APBS value: !isnan(e.norm()) AND e.norm() < 100.0 | APBS E-field vector | Decision: skip NaN or unreasonably large (> 100 V/A) fields | conditional | n/a | op 467 |
| 469 | Call extract_atom(matched, d, rings, pdb_atoms, ctx, protein_id, graph, atom_efield) | all inputs | Produce features for this atom | IrrepFeatures | all | op 466 |

Magic numbers: 1.0 A (spatial matching threshold), 100.0 V/A (APBS sanity threshold).

---

## Feature name registries (lines 1542-1740)

### L0_names() (lines 1542-1631)

Returns 121 L0 feature names. Count breakdown:
- ring_current_features: 7 (type G_T0) + 1 (B_mag) + 4 (ring counts) + 4 (ring dists) + 20 (near1) + 4 (near2) + 5 (multi-ring) = 45
- bond_anisotropy_features: 5
- coulomb_field_features: 4
- atom_identity_features: 10
- structural_features: 10
- residue_encoding: 20
- haigh_mallion_features: 9
- graph_features: 10
- full_bond_anisotropy_features: 3
- apbs_field_features: 2
- hbond_features: 6
- dispersion_tensor_features: 2
- aromatic_field_features: 3
- per_type_nearest_features: 4
- quadrupole_field_features: 1
- ring_chi_aniso_features: 2
- **Total L0: 136** (from summing above: 45+5+4+10+10+20+9+10+3+2+6+2+3+4+1+2)

### L1_names() (lines 1634-1654)

Returns 17 L1 feature names:
- jb_B_total, jb_B_near1, jb_B_near2 (3, 1e)
- ring_normal_1, ring_normal_2 (2, 1e)
- dir_to_ring_1, dir_to_ring_2 (2, 1o)
- mcconnell_dir_nearest_co (1, 1o)
- coulomb_E_total, coulomb_E_backbone (2, 1o)
- struct_dir_acceptor, struct_bond_direction (2, 1o)
- hm_B_total (1, 1e)
- bond_aniso_dir_nearest (1, 1o)
- apbs_E_total (1, 1o)
- hbond_dir (1, 1o)
- arom_E_total (1, 1o)
- **Total L1: 17** (6 pseudovectors (1e) + 11 polar vectors (1o))

### L1_parities() (lines 1657-1667)

Returns 17 parity strings: "e","e","e","e","e","o","o","o","o","o","o","o","e","o","o","o","o"

### L2_names() (lines 1669-1714)

Returns 35 L2 feature names:
- ring_current_features: 7 (per-type G_T2) + 3 (near1, near2, exp_sum) = 10
- bond_anisotropy_features: 2 (CO_nearest, CN_nearest)
- coulomb_field_features: 2 (EFG total, EFG backbone)
- haigh_mallion_features: 7 (per-type) + 1 (nearest) = 8
- full_bond_anisotropy_features: 3 (backbone, sidechain, aromatic)
- apbs_field_features: 1
- hbond_features: 1
- dispersion_tensor_features: 1
- aromatic_field_features: 1
- per_type_nearest_features: 4
- quadrupole_field_features: 1
- ring_chi_aniso_features: 1
- **Total L2: 35** (all even parity, 2e)

### L2_parities() (lines 1717-1722)

All L2 features are even parity.

### e3nn_irreps_string() (lines 1724-1740)

| # | What | Inputs | Category/Decision | Output | Irrep | Dependencies |
|---|------|--------|--------------------|--------|-------|--------------|
| 470 | Count L0, L1e, L1o, L2 features from name/parity registries | name and parity lists | none | string like "136x0e+6x1e+11x1o+35x2e" | n/a | ops 448-450 |

---

## Summary of all magic numbers and sentinel values

| Value | Where used | Meaning |
|-------|-----------|---------|
| 99.0 | AtomSite defaults, rd() lambda, ring_center_dist, graph_dist sentinel | NO_DATA_SENTINEL: "no data available" |
| 1.8 | AtomSite::build parent search | Maximum H-X bond length (Angstroms) |
| 1e-10 | safe_unit, to_ring_cylindrical, d2_d1_ratio | Near-zero threshold for vector norms |
| 1e-15 | Coulomb E_bond_proj, APBS E_bond, structural B_bond_cos | Near-zero threshold for field magnitudes |
| 1e20 | bond_anisotropy best_co_dist/best_cn_dist guard | "Not found" threshold (initial values are 1e30) |
| 1e30 | bond_anisotropy, per_type_nearest, ring_chi_aniso initial best_dist | "Not found" initial value |
| 0.1 | Constants::MIN_DISTANCE | Singularity guard (Angstroms) |
| 3.0 | Ring count threshold | Ring within 3 Angstroms |
| 5.0 | Ring count threshold, dispersion R_CUT | Ring within 5 A / dispersion cutoff |
| 8.0 | Ring count threshold, packing density (64.0 = 8^2), G_iso variance | Ring within 8 A / heavy atom packing |
| 12.0 | Ring count threshold | Ring within 12 A |
| 15.0 | Quadrupole and chi_aniso cutoff | Ring contribution cutoff |
| 4.0 | Exponential decay length | exp(-r/4) weighting |
| 0.1 | co_vec.norm() threshold, ring_center_dist threshold | Minimum vector magnitude for valid C=O |
| 90.0 | nearest_co_dist threshold | Maximum C=O distance for McConnell |
| 50.0 | hbond inv_d3 threshold | Maximum H-bond distance for 1/d^3 |
| 3.5 | hbond count radius | H-bond partner counting shell |
| 1.5 | dispersion R_MIN | Minimum distance for dispersion (covalent exclusion) |
| 8.0 | dispersion early exit (R_CUT + 3.0) | Ring-center distance for early break |
| 1.0 | extract_all spatial matching threshold | Position matching tolerance (Angstroms) |
| 100.0 | APBS E-field sanity | Maximum E-field magnitude (V/A) |
| 2 | hbond sequential exclusion threshold | Residue number difference for N...O exclusion |

---

## Summary of all categorical classification decisions

| Operation(s) | Classification | Where used | Note |
|--------------|---------------|-----------|------|
| ops 1-7 | Backbone atom identity | AtomSite::build, hbond_features | Set: {N, CA, C, O, H, HA, HN} |
| op 364 | Backbone + CB atom identity | aromatic_field_features | BROADER set: adds {HA2, HA3, CB, HB, HB2, HB3} |
| ops 8-9 | Amide H identity | AtomSite::build | Set: {H, HN} |
| ops 10-12 | Alpha H identity | AtomSite::build | Set: {HA, HA2, HA3} |
| ops 13-18 | Methyl H identity | AtomSite::build | Pattern: H + {B,G,D,E} + digit |
| ops 19-20 | Ring nitrogen detection | ring_current_features | First-character check only |
| op 21 | Ring type (7-class) | Multiple methods | PHE/TYR/TRP6/TRP5/HIS/HID/HIE |
| ops 30-33 | sp2 parent hybridization | AtomSite::build | N, C (backbone) + 13 ring atom names |
| op 52 | Aromatic residue | AtomSite::build | Delegates to is_aromatic_residue |
| op 62 | H-bond acceptor element | AtomSite::build | O or N (any non-same-residue) |
| ops 240-244 | Element one-hot (5-class) | atom_identity_features | H, C, N, O, S |
| ops 271-273 | Residue one-hot (20-class) | residue_encoding | Standard amino acids, HID/HIE/HIP -> HIS |
| op 131 | Nearest ring type one-hot (7-class) | ring_current_features | Per RingCurrentParams types |
| ops 298-300 | Hybridization one-hot (3-class) | graph_features | sp, sp2, sp3 |
| op 324 | H-bond donor N atoms | hbond_features | Set: {N, ND2, NE2, NZ, NH1, NH2} |
| op 325 | H-bond acceptor O atoms | hbond_features | Set: {O, OD1, OD2, OE1, OE2, OG, OG1, OH} |
| op 365 | Aromatic residue (inline) | aromatic_field_features | Same as is_aromatic_residue but inlined |
| op 387 | Collapsed ring type (4-class) | per_type_nearest_features | PHE/TYR/TRP/HIS-family |
| ops 332-333 | Sequential backbone exclusion | hbond_features | |resnum diff| < 2 for N...O or O...N |

---

## Summary of all physics computations called

| Physics Model | Method(s) | External function | What it computes |
|---------------|-----------|-------------------|------------------|
| Johnson-Bovey B-field | ring_current_features | johnson_bovey_field() | B(r) from double-loop pi-electron model |
| JB geometric factor | ring_current_features, per_type_nearest | ring_geometric_factor_jb() | Full 3x3 shielding geometric factor G |
| McConnell bond factor | bond_anisotropy, structural | mcconnell_bond_factor() | (3cos^2 theta - 1) / R^3 from bond midpoint |
| Coulomb E-field | coulomb_field, aromatic_field | inline q*dr/r^3 sum | E(r) = sum q_i (r-r_i)/|r-r_i|^3 |
| EFG tensor | coulomb_field, aromatic_field | inline dipolar kernel sum | V_ab = q (3dr_a dr_b/r^5 - delta_ab/r^3) |
| Haigh-Mallion tensor | haigh_mallion_features | haigh_mallion_tensor_factor() | Surface-integral geometric factor |
| HM B-field | haigh_mallion_features | haigh_mallion_B_field() | HM B-field from ring surface integral |
| Bond anisotropy | full_bond_anisotropy | compute_bond_anisotropy() | McConnell tensor for all graph bonds by category |
| Dispersion kernel | dispersion_tensor | inline 1/R^6 and 1/R^8 sum | London-weighted dipolar tensor |
| Quadrupole EFG | quadrupole_field | inline second derivative of axial quadrupole potential | Pi-electron quadrupole gradient tensor |
| Ring chi aniso | ring_chi_aniso | inline (3cos^2 theta - 1)/r^3 from ring center | Whole-ring McConnell factor |
| H-bond dipolar tensor | hbond_features | inline dipolar kernel | 3d_a d_b / R^5 - delta_ab / R^3 |
| APBS proxy EFG | apbs_field | inline E outer product | E_a E_b / |E|^2 decomposed to T2 |
| Spherical decomposition | all methods producing T2 | decompose_tensor() | T0 (trace/3) + T2 (5 STF components) |

---

## Total feature counts (from registry)

| Irrep | Count | Components per feature | Total scalars |
|-------|-------|----------------------|---------------|
| L0 (0e) | 136 | 1 | 136 |
| L1e | 6 | 3 | 18 |
| L1o | 11 | 3 | 33 |
| L2e | 35 | 5 | 175 |
| **Total** | **188 features** | | **362 scalars** |

e3nn irreps string: `136x0e+6x1e+11x1o+35x2e`

---

## Structural dependency graph

```
AtomSite::build (ops 36-77)
  |
  +-- ring_current_features (ops 78-184) --> returns B_total
  |
  +-- bond_anisotropy_features (ops 185-210) -- needs ProteinContext
  |
  +-- coulomb_field_features (ops 211-239) -- needs ProteinContext, rings
  |
  +-- atom_identity_features (ops 240-249) -- AtomSite only
  |
  +-- structural_features (ops 250-270) -- needs B_total from ring_current
  |
  +-- residue_encoding (ops 271-273) -- AtomSite only
  |
  +-- haigh_mallion_features (ops 274-289) -- needs rings
  |
  +-- graph_features (ops 290-302) -- needs MolecularGraph (optional)
  |
  +-- full_bond_anisotropy_features (ops 303-312) -- needs MolecularGraph (optional)
  |
  +-- apbs_field_features (ops 313-323) -- needs APBS E-field (optional)
  |
  +-- hbond_features (ops 324-349) -- needs pdb_atoms
  |
  +-- dispersion_tensor_features (ops 350-362) -- needs rings
  |
  +-- aromatic_field_features (ops 363-383) -- needs pdb_atoms, ctx, rings
  |
  +-- per_type_nearest_features (ops 384-392) -- needs rings
  |
  +-- quadrupole_field_features (ops 393-408) -- needs rings
  |
  +-- ring_chi_aniso_features (ops 409-424) -- needs rings
```

Only one inter-method dependency: structural_features depends on B_total from ring_current_features. All other methods depend only on AtomSite and their direct input data.

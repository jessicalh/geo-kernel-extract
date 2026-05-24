# haigh-mallion — codex review

- **Targets:** src/HaighMallionResult.{h,cpp}
- **Model:** codex `gpt-5.5`, reasoning `xhigh`
- **Date:** 2026-05-24
- **Tokens:** 22395
- **Brief:** `../_brief.md`

---

Scope: reviewed the supplied `HaighMallionResult.h/.cpp` only. Contract-dependent items are marked “check”.

**1. Algorithm Correctness**
`src/HaighMallionResult.h:15` — documented `G` formula omits the minus sign used in `cpp:287-294`; I don’t know which convention is authoritative — check BiotSavart/literature sign and make them agree.

`src/HaighMallionResult.h:25` — `G = n (x) V` reads transposed relative to the component formula/code (`V_a n_b`) — use `G_ab = -V_a n_b`.

`src/HaighMallionResult.cpp:117` — subdivision tests only vertex distances, so close interior/edge approaches can be under-refined — check whether ring geometry/thresholds guarantee coverage.

`src/HaighMallionResult.cpp:310` — `d = atom_pos - geom.center` is center-to-atom, then stored as `direction_to_center`; if the field name is literal, sign is reversed — check `RingNeighbourhood` convention.

`src/HaighMallionResult.cpp:316` — `abs(z)` folds above/below ring into one `theta`; correct only if `theta` is intentionally unsigned — check stored-angle convention.

`src/HaighMallionResult.cpp:375` — sampler near-field guard uses center distance/radius, not distance to ring surface; it can skip finite off-plane points and still allow near-edge cases — check sampler domain.

No definite issue found in the quadrature constants, triangle area scaling, kernel component formula, or fan loop bounds.

**2. Variable Naming Clarity**
`src/HaighMallionResult.cpp:81` — `r` hides that this is the target/field point — `field_point`.

`src/HaighMallionResult.cpp:90` — `rS` is cryptic for quadrature surface point — `surface_point`.

`src/HaighMallionResult.cpp:285` — `V` could read as potential/velocity, not `H.n` field — `effective_field`.

`src/HaighMallionResult.cpp:310` — `d` hides frame/direction — `center_to_atom`.

`src/HaighMallionResult.cpp:315` — scalar `rho` reuses the separation-vector symbol for ring-plane radius — `in_plane_radius`.

`src/HaighMallionResult.cpp:316` — `theta` hides the unsigned/folded convention — `unsigned_theta`.

**3. Grouping Of Math Operations**
`src/HaighMallionResult.cpp:291` — outer-product construction is a meaningful tensor step but appears as raw component loops — label as `outer product`.

`src/HaighMallionResult.cpp:313` — ring-frame coordinate derivation mixes signed height, plane projection, radius, and folded angle — label as `ring-frame coordinates`.

`src/HaighMallionResult.cpp:325` — storage block mixes raw tensors, decompositions, and contractions with long inline comments — split labels: `raw integral`, `shielding kernel`.

**4. Method / Function Naming**
`src/HaighMallionResult.cpp:44` — `Gauss7` omits triangle/rule context — `TriangleGauss7Rule`.

`src/HaighMallionResult.cpp:79` — `AccumulateTensor` does not say which tensor/domain — `AccumulateTriangleH`.

`src/HaighMallionResult.cpp:125` — `AccumulateAdaptive` hides that it is adaptive triangle-H accumulation — `AccumulateAdaptiveTriangleH`.

`src/HaighMallionResult.cpp:158` — `SurfaceIntegral` does not say it returns raw `H`, not shielding `G` — `ComputeSurfaceIntegralH`.

`src/HaighMallionResult.cpp:401` — `PackST_HM` is abbreviated — `PackHMSphericalTensor`.

**5. Comments**
`src/HaighMallionResult.cpp:32` — misleading area normalization wording; weights sum to 1 in this implementation — `area-normalized weights`.

`src/HaighMallionResult.cpp:108` — “close to a triangle” overstates the vertex-distance test — `vertex-distance refinement`.

`src/HaighMallionResult.cpp:244` — TODO is stale/unclear in this math path — remove or replace with `source ring record`.

`src/HaighMallionResult.cpp:281` — verbose step comments restate nearby code — `surface integral`.

`src/HaighMallionResult.cpp:284` — verbose contraction comment — `normal contraction`.

`src/HaighMallionResult.cpp:287` — good sign note, but too long across three lines — `shielding sign`.

`src/HaighMallionResult.cpp:326` — inline comments overclaim exact `T0/T1` zeros after floating-point quadrature — `raw H integral`.

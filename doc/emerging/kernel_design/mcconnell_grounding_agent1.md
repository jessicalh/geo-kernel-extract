# McConnell grounding — Agent 1 (codex xhigh; web + local corpus; write-nothing), 2026-06-06

Grounding-path-1 of the two Jessica named (code/web/e3nn reality-check; the grinder corpus-walk is path-2, parked). Brief: /tmp/mcconnell_grounding_brief.txt.

---

**Verdict**

The core McConnell replacement should be built on the point-dipole susceptibility propagator, but I would tighten the current synthesis in one important way: the literature-standard object is not “the old `M` made cleaner.” It is the dipole propagator acting on a source susceptibility tensor, decomposed into declared irreps. If you emit only the symmetric-traceless part, call that the primary `2e` McConnell kernel. If you emit the full shielding response, it is generally `0e ⊕ 1e ⊕ 2e`, not a single `2e`.

The clean pair expression should be:

```text
D_is = (3 rhat_is rhat_is^T - I) / r_is^3

Q_s = Δχ_ax,c (u_s u_s^T - I/3)
    + Δχ_rh,c/2 (v_s v_s^T - w_s w_s^T)   where rhombic data exist

σ_is^McC = κ · D_is · [ b_s · Q_s ]
```

with sign/prefactor absorbed into calibration if needed. `b_s` is the MOPAC source-strength modifier, preferably a calibrated function of Wiberg bond order, not a hard-coded claim that susceptibility is exactly linear in bond order.

**Evidence**

The current code really does compute the older mixed tensor. [McConnellResult.cpp](/shared/2026Thesis/nmr-shielding/src/McConnellResult.cpp:35) documents and implements `M = 9 cosθ d⊗b - 3 b⊗b - (3 d⊗d - I)`, stores `K` separately, then accumulates full `M` and decomposes it after symmetrizing category sums. [MopacMcConnellResult.cpp](/shared/2026Thesis/nmr-shielding/src/MopacMcConnellResult.cpp:95) only multiplies that same old object by MOPAC bond order. So the MOPAC seed is useful, but it is not yet the proper expression.

The field-standard anchor is the PCS/McConnell point-dipole form. Case’s biomolecular review gives the remote-group shielding tensor as susceptibility times the dipole field propagator. Suturina and Kuprov write the same point tensor and its isotropic PCS average, then rewrite the scalar through second-rank spherical harmonics with five irreducible susceptibility components; their text also says isotropic and rank-1 susceptibility components do not enter PCS ([local chunk](/shared/2026Thesis/nmr-shielding/references-text/suturina-2016-pseudocontact-shifts-from-mobile-spin-labels-text-2.txt:53), [RSC article](https://pubs.rsc.org/en/content/articlehtml/2016/CP/C6CP05437D)). Ben Mahmoud et al. explicitly decompose NMR rank-2 tensors as `T(0) ⊕ T(1) ⊕ T(2)` for equivariant ML, with even parity for all three components ([local chunk](/shared/2026Thesis/nmr-shielding/references-text/ben-mahmoud-2024-gnn-solid-state-nmr-spherical-tensors-text-1.txt:191), [arXiv](https://arxiv.org/abs/2412.15063)).

**Irreps And Parity**

Use e3nn `0e`, `1e`, `2e` for shielding/susceptibility tensor pieces. Do not use `1o` for the antisymmetric piece of a rank-2 shielding tensor. e3nn’s own representation rules say irreps carry degree and parity, and tensor-product parity multiplies; two polar vectors (`1o × 1o`) give even branches, including an even pseudovector ([e3nn docs](https://docs.e3nn.org/en/stable/api/o3/o3_irreps.html), [local e3nn chunk](/shared/2026Thesis/nmr-shielding/references-text/geiger-smidt-2022-e3nn-euclidean-neural-networks-text-4.txt:71)). So the prompt’s “outer-product-of-polar-vectors gives `1o`” warning is not correct under e3nn convention.

If the thesis wants the McConnell responsibility kernel to be pure and simple, emit the `2e` symmetric-traceless projection as the primary feature. But the scalar PCS is a `0e` coupling of `D` and `χ`; it is not the trace of a pure `2e` object. That point in the one-agent synthesis is too compressed.

**MOPAC Composition**

MOPAC bond order belongs in the source model. OpenMOPAC states its `BONDS/ALLBONDS` output is rotationally invariant and equivalent to Wiberg bond order, with ethane/ethylene/acetylene roughly 1/2/3; `ALLBONDS` lowers the threshold and includes H bonds ([OpenMOPAC BONDS](https://openmopac.net/Manual/bonds.html), [ALLBONDS](https://openmopac.net/Manual/allbonds.html)). Your local capture also has the right fuller shape: atom populations, electron density, AO populations, valencies, directed printed bond rows, topology bridge, and full accessors ([MopacResult.h](/shared/2026Thesis/nmr-shielding/src/MopacResult.h:90), [mopac_extension.md](/shared/2026Thesis/nmr-shielding/doc/emerging/kernel_design/mopac_extension.md:141)).

The honest calibrated combination is:

```text
b_s model:
  b_s = 1 + λ_c (BO_s / BO_ref,c - 1)
or equivalently emit two basis channels:
  Σ D_is Q_s                 fixed category source
  Σ BO_s D_is Q_s            Wiberg-weighted source
and learn β_c0, β_c1.
```

That makes MOPAC count without pretending the literature proves exact linear `Δχ ∝ Wiberg`. This is consistent with broader chemistry practice where AM1-Wiberg fractional bond orders interpolate valence parameters, but it should be described as a measured QM source-strength modulation plus calibration, not a first-principles susceptibility calculation ([OpenFF SMIRNOFF](https://openforcefield.github.io/standards/standards/smirnoff/)).

MOPAC charges/populations should mostly stay `0e` conditioning scalars or feed Charge/EFG. They should not be folded into McConnell as invented “population-to-Δχ” corrections unless you later find a specific susceptibility model.

**Complete Emit**

Per non-aromatic source category:

```text
mc_<cat>_0e_fixed, mc_<cat>_1e_fixed, mc_<cat>_2e_fixed
mc_<cat>_0e_bo,    mc_<cat>_1e_bo,    mc_<cat>_2e_bo
```

If you choose the stricter “primary McConnell is pure T2” thesis surface, keep:

```text
mc_<cat>_2e_fixed
mc_<cat>_2e_bo
mc_<cat>_pcs0_fixed
mc_<cat>_pcs0_bo
```

and label the scalar as PCS/McConnell scalar, not as a projection of a pure `2e`.

Aromatic McConnell must be zero while Biot-Savart/Haigh-Mallion are active. The spine already states this partition ([PIPELINE_SPINE.md](/shared/2026Thesis/nmr-shielding/doc/emerging/kernel_design/PIPELINE_SPINE.md:327)); the current code still accumulates aromatic McConnell channels, including in the MOPAC path ([MopacMcConnellResult.cpp](/shared/2026Thesis/nmr-shielding/src/MopacMcConnellResult.cpp:221)).

**Limits**

Keep the point-dipole ceiling explicit. Suturina shows distributed-source corrections can exceed 1 ppm out to ~10 Å in mobile spin-label PCS examples and that near the source extra multipoles are needed ([local chunk](/shared/2026Thesis/nmr-shielding/references-text/suturina-2016-pseudocontact-shifts-from-mobile-spin-labels-text-4.txt:37)). Martin and Brown show McConnell’s alkene shielding-cone model can have the wrong sign within about 3 Å because it omits electric field, orbital, and dispersion terms ([local summary](/shared/2026Thesis/nmr-shielding/references-meta/martin-2000-graphical-model-proton-deshielding-summary.txt:5)). Do not hide that with cleaner equivariance.

**What Changes**

The plan is right to replace current `M`. It is right to use spherical irreps, calibrated Δχ, aromatic zeroing, and MOPAC bond orders. It is thin where it calls the contracted PCS form a single `2e` and where it treats trace/asymmetry as automatically unphysical. The defensible thesis version is: clean source physics, correct even-parity irreps, MOPAC bond order as calibrated source strength, MOPAC populations as scalar conditioning, and the point-dipole approximation disclosed as a limit.

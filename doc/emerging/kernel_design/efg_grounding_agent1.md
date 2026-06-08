Historical -- not current truth; see doc/emerging/kernel_design/efg_spec.md and doc/emerging/CONTROLLING_SPEC.md.
# EFG / charge grounding — first-stage (codex xhigh; read-only; web + corpus), 2026-06-06

Pre-staged for the EFG work-through (after McConnell). Brief: /tmp/efg_grounding_brief.txt. Companion to mcconnell_grounding_agent1.md.

---

## Cross-kernel: ring ↔ EFG overlap (2026-06-06)

Ring (BS/HM) and charge/EFG do **not** overlap at the *mechanism* level — ring is **magnetic**
(induced ring current → secondary magnetic field); EFG is **electric** (static field + gradient from
charges). Neither uses the other's source, so there is **no double-count to zero** (contrast
ring↔McConnell, which *is* the same magnetic physics → aromatic-zero). They **coexist** as
complementary channels.

They **do** correlate at the *geometric-descriptor* level near aromatic groups: the ring `2e` and the
aromatic part of the EFG `2e` are both low-order multipole moments of the *same* neighbourhood (the
dual-nature — a kernel is a physical field AND a geometric descriptor). That is **correlation, not
double-counting** — the law study (NEW step 1) splits shared vs unique variance by ablation and
*reports* the collinearity; it is NOT zeroed.

Real overlaps live **within** families: ring↔McConnell (magnetic — aromatic-zero), **EFG↔π-quadrupole**
(electric — partition-contested, the electric double-count risk; see the π-Quadrupole section), the
EFG source-fork (MOPAC/FF14SB/APBS), and MOPAC-populations↔AIMNet2.

---

**Parity / Irreps**

The controlling spec is right: charge/EFG should be `1o ⊕ 2e` [CONTROLLING_SPEC.md](/shared/2026Thesis/nmr-shielding/doc/emerging/CONTROLLING_SPEC.md:131).

The electric field is a true polar vector. Under inversion, position changes sign and charge does not, so the field changes sign. In e3nn language this is `1o`. e3nn’s irreps use angular order `l` plus inversion parity, with vectors as `l=1` and odd/even parity tracked separately; the local e3nn paper explains this convention around irreps and parity [geiger-smidt-2022 text](/shared/2026Thesis/nmr-shielding/references-text/geiger-smidt-2022-e3nn-euclidean-neural-networks-text-2.txt:59). Official e3nn docs use the same `1o`, `2e` notation: https://docs.e3nn.org/en/stable/api/o3/o3_irreps.html.

The EFG is the gradient of the electric field, equivalently the Hessian of a scalar electrostatic potential. Gradient is odd and field is odd, so the product is even. After symmetrization and traceless projection it is a pure `2e` object. This matches modern NMR tensor practice: EFG is treated as a symmetric traceless rank-2 tensor, naturally represented as an `l=2` spherical tensor [ben-mahmoud-2024 text](/shared/2026Thesis/nmr-shielding/references-text/ben-mahmoud-2024-gnn-solid-state-nmr-spherical-tensors-text-3.txt:33).

This differs from the McConnell shielding tensor case. Shielding itself is an even response tensor, so its decomposition is `0e/1e/2e`. That “never `1o`” rule must not be imported here. The charge field is not a shielding tensor; it is an electrostatic input to shielding.

**Literature Grounding**

Buckingham is the right anchor. The held summary says Buckingham expands shielding in the local electric field, with a linear field term often dominating the quadratic term for polar groups [buckingham-1960 summary](/shared/2026Thesis/nmr-shielding/references-meta/buckingham-1960-chemical-shifts-polar-groups-summary.txt:7). Case’s ring-current calibration also includes a Buckingham-style electrostatic polarization term using AMBER partial-charge fields [case-1995 summary](/shared/2026Thesis/nmr-shielding/references-meta/case-1995-ring-current-calibration-summary.txt:5). Sitkoff’s peptide DFT work supports electric-field polarization terms and shows local charge geometry matters, especially near lone pairs and polar sites [sitkoff-1997 summary](/shared/2026Thesis/nmr-shielding/references-meta/sitkoff-1997-density-functional-proton-chemical-shifts-model-peptides-summary.txt:5). Sahakyan and Vendruscolo show ring-current and electric-field terms can be calibrated separately against DFT, with electric-field contributions especially important for heavier nuclei [sahakyan-2013 summary](/shared/2026Thesis/nmr-shielding/references-meta/sahakyan-vendruscolo-2013-ring-current-electric-field-contributions-summary.txt:5).

For APBS, the relevant claim is not “better charges”; it is continuum solvent screening. APBS solves the Poisson-Boltzmann electrostatic problem using dielectric, ionic-accessibility, and solute charge information [baker-2001 summary](/shared/2026Thesis/nmr-shielding/references-meta/baker-2001-apbs-electrostatics-nanosystems-summary.txt:5).

**Source Fork: Map, Don’t Decide**

The central ambiguity is real and should stay open.

FF14SB charges are stable, cheap, MD-consistent, and historically connected to the Case/Sitkoff style of field correction. Their weakness is that they are fixed force-field charges: conformation affects geometry, not charge redistribution. They also currently use cutoff machinery, which makes them operationally different from MOPAC.

MOPAC Mulliken charges are conformation-dependent and may capture electronic redistribution missing from FF14SB. Their weakness is that Mulliken charges are method-dependent population analyses, not observables, and the current implementation changes both charge source and pair range at once.

APBS is best treated as a screened/reaction-field channel: “FF14SB-like charge model plus continuum PB environment,” not a peer charge assignment method. Its strength is solvent physics; its current weakness is PB parameter validity, especially radii.

So the defensible design is parallel source channels plus reconciliation diagnostics. A “better-replaces-cruft” outcome is possible if one source dominates DFT and experiment after ablation. A “butterfly-friends” outcome is also plausible: agreement between MOPAC, FF14SB, and APBS-derived fields can be treated as validation of a robust electrostatic signal. The existing MOPAC-vs-FF14SB reconciliation cosine is already pointing in this direction.

**π-Quadrupole Partition**

This must be flagged, not resolved. The project’s own π-quadrupole note states the risk clearly: if charge/EFG already includes aromatic atom partial charges, then a point π-quadrupole can double-count the same far-field electrostatic multipole unless the partition is explicit [pi_quadrupole.md](/shared/2026Thesis/nmr-shielding/doc/emerging/kernel_design/pi_quadrupole.md:123). The current `PiQuadrupoleResult` computes a rank-2 EFG-like quadrupole tensor and scalar, not a full `1o ⊕ 2e` quadrupole field. Treat π-quad as a participant whose partition against charge/EFG is contested.

**Defensible Generous Emit Shape**

For each source channel, emit:

`E_1o`: `(N, 3)`, V/Angstrom, polar vector.  
`EFG_2e`: `(N, 5)`, V/Angstrom², real spherical symmetric traceless tensor.  
`0e` diagnostics: `|E|`, `|E|^2`, bond projection `E·b`, EFG Frobenius norm, eigenvalue/asymmetry summaries, source-partition magnitudes, cutoff/source-count diagnostics, and source agreement cosines.

Source channels should remain separate:

`ff14sb_vacuum_total`, plus backbone/sidechain/aromatic where available.  
`mopac_vacuum_total`, plus the same partitions.  
`apbs_pb_total`.  
`apbs_minus_ff14sb_vacuum` as a solvent/reaction-field delta when both are available.

I would not collapse source channels into one averaged field before the law study. Fido can take everything; the tensor model needs clean irreps and source tags; the law study needs ablatable source comparisons.

**Feed The Three New Steps**

Step 1, law/correlation study: use scalar and angular comparisons per source. Scalar-only baselines should include `|E|`, `E·bond`, `|EFG|`, EFG eigenvalue/asymmetry summaries, and source agreement metrics. Angular models should receive full `1o` and `2e` objects. This directly tests “with angular part vs without angular part” against DFT.

Step 2, equivariant shielding-tensor model: feed source-tagged `1o ⊕ 2e` channels plus `0e` diagnostics. The target shielding tensor remains even, so the model may consume `1o` field inputs but should emit shielding irreps as even tensor components, primarily `0e ⊕ 2e` unless antisymmetric shielding is deliberately modeled.

Step 3, experimental-shift fido predictor: include everything ablatable: raw field/EFG channels, current NPY scalar exports, APBS deltas, MOPAC/FF14SB agreement, source partitions, π-quad participant features, AIMNet2 charge/EFG features if available, and diagnostics marking cutoff/radii/source validity.

**Thin Spots**

The biggest thin spots are naming and comparability. Current `*_shielding.npy` files contain bare EFG tensors, not calibrated shielding. FF14SB and MOPAC differ in both charge model and cutoff/all-pair treatment. Sidechain EFG is computed but not cleanly exported like total/backbone/aromatic. APBS has explicit parameter-validity warnings and a clamp behavior the code itself calls a physics/blessing question. Static NPY outputs do not carry parity/unit metadata as cleanly as the newer H5 APBS trajectory outputs.

The honest scientific conclusion: the kernel form is defensible, the parity is clear, and the source fork is not decidable from theory alone. Keep the emit generous, keep the sources separate, and let DFT plus experimental ablations decide whether one source replaces the rest or whether agreement between independent electrostatic channels is itself useful signal.
tokens used
273,837
**Bottom Line**

The charge/EFG calculator is genuinely an electrostatic field kernel, not just a mislabeled shielding feature. The defensible forward object is:

`charge/EFG = E_1o ⊕ V_2e`

where `E` is the polar electric field and `V` is the symmetric traceless electric-field gradient. The existing code mostly computes that physics, but the surface is uneven: historical “shielding” filenames actually contain bare EFG-like tensors, source handling is not normalized, and APBS is a screened PB field rather than just another charge set.

**What The Code Computes**

`CoulombResult` computes vacuum Coulomb electric field and EFG from FF14SB partial charges:

`E_a = k_e sum_j q_j r_a / r^3`  
`V_ab = k_e sum_j q_j (3 r_a r_b / r^5 - delta_ab / r^3)`

That formula is explicit in [CoulombResult.cpp](/shared/2026Thesis/nmr-shielding/src/CoulombResult.cpp:30). It uses a cutoff/spatial index and filters self or near-self sources, then decomposes total, backbone, sidechain, and aromatic source contributions. It stores `E_total`, `EFG_total`, selected source partitions, `|E|`, bond projection, and an aromatic field magnitude. The `coulomb_shielding.npy` export is a historical name: it packs the decomposed EFG as a 9-component spherical tensor, not a calibrated shielding tensor [CoulombResult.cpp](/shared/2026Thesis/nmr-shielding/src/CoulombResult.cpp:310).

`MopacCoulombResult` computes the same Coulomb kernel, but from MOPAC Mulliken charges rather than FF14SB charges [MopacCoulombResult.cpp](/shared/2026Thesis/nmr-shielding/src/MopacCoulombResult.cpp:29). Important difference: it is full all-pairs, not the same cutoff implementation as FF14SB [MopacCoulombResult.cpp](/shared/2026Thesis/nmr-shielding/src/MopacCoulombResult.cpp:41). That means “MOPAC vs FF14SB” currently compares both charge model and range treatment unless this is made explicit or harmonized.

`ApbsFieldResult` is different in kind. It samples an APBS Poisson-Boltzmann potential grid, computes `E = -grad(phi)` by finite differences, computes the gradient of `E`, symmetrizes it, projects it traceless, then converts to V/Angstrom and V/Angstrom² [ApbsFieldResult.cpp](/shared/2026Thesis/nmr-shielding/src/ApbsFieldResult.cpp:72). This is a solvated continuum electrostatic field, not an independent Mulliken-like charge source. The code itself warns that placeholder PB radii make APBS “populated/finite but NOT physically validated” when authoritative radii are missing [ApbsFieldResult.cpp](/shared/2026Thesis/nmr-shielding/src/ApbsFieldResult.cpp:125).

The tensor machinery is appropriate: `SphericalTensor` packs `T0`, `T1`, and `T2`; `Mat3::Decompose()` gives trace scalar, antisymmetric vector, and symmetric traceless rank-2 components in a real spherical basis [Types.cpp](/shared/2026Thesis/nmr-shielding/src/Types.cpp:25). For a physical EFG from external electrostatics, `T0` and `T1` should be structural zeros after projection; the useful object is `T2`.

**Parity / Irreps**

The controlling spec is right: charge/EFG should be `1o ⊕ 2e` [CONTROLLING_SPEC.md](/shared/2026Thesis/nmr-shielding/doc/emerging/CONTROLLING_SPEC.md:131).

The electric field is a true polar vector. Under inversion, position changes sign and charge does not, so the field changes sign. In e3nn language this is `1o`. e3nn’s irreps use angular order `l` plus inversion parity, with vectors as `l=1` and odd/even parity tracked separately; the local e3nn paper explains this convention around irreps and parity [geiger-smidt-2022 text](/shared/2026Thesis/nmr-shielding/references-text/geiger-smidt-2022-e3nn-euclidean-neural-networks-text-2.txt:59). Official e3nn docs use the same `1o`, `2e` notation: https://docs.e3nn.org/en/stable/api/o3/o3_irreps.html.

The EFG is the gradient of the electric field, equivalently the Hessian of a scalar electrostatic potential. Gradient is odd and field is odd, so the product is even. After symmetrization and traceless projection it is a pure `2e` object. This matches modern NMR tensor practice: EFG is treated as a symmetric traceless rank-2 tensor, naturally represented as an `l=2` spherical tensor [ben-mahmoud-2024 text](/shared/2026Thesis/nmr-shielding/references-text/ben-mahmoud-2024-gnn-solid-state-nmr-spherical-tensors-text-3.txt:33).

This differs from the McConnell shielding tensor case. Shielding itself is an even response tensor, so its decomposition is `0e/1e/2e`. That “never `1o`” rule must not be imported here. The charge field is not a shielding tensor; it is an electrostatic input to shielding.

**Literature Grounding**

Buckingham is the right anchor. The held summary says Buckingham expands shielding in the local electric field, with a linear field term often dominating the quadratic term for polar groups [buckingham-1960 summary](/shared/2026Thesis/nmr-shielding/references-meta/buckingham-1960-chemical-shifts-polar-groups-summary.txt:7). Case’s ring-current calibration also includes a Buckingham-style electrostatic polarization term using AMBER partial-charge fields [case-1995 summary](/shared/2026Thesis/nmr-shielding/references-meta/case-1995-ring-current-calibration-summary.txt:5). Sitkoff’s peptide DFT work supports electric-field polarization terms and shows local charge geometry matters, especially near lone pairs and polar sites [sitkoff-1997 summary](/shared/2026Thesis/nmr-shielding/references-meta/sitkoff-1997-density-functional-proton-chemical-shifts-model-peptides-summary.txt:5). Sahakyan and Vendruscolo show ring-current and electric-field terms can be calibrated separately against DFT, with electric-field contributions especially important for heavier nuclei [sahakyan-2013 summary](/shared/2026Thesis/nmr-shielding/references-meta/sahakyan-vendruscolo-2013-ring-current-electric-field-contributions-summary.txt:5).

For APBS, the relevant claim is not “better charges”; it is continuum solvent screening. APBS solves the Poisson-Boltzmann electrostatic problem using dielectric, ionic-accessibility, and solute charge information [baker-2001 summary](/shared/2026Thesis/nmr-shielding/references-meta/baker-2001-apbs-electrostatics-nanosystems-summary.txt:5).

**Source Fork: Map, Don’t Decide**

The central ambiguity is real and should stay open.

FF14SB charges are stable, cheap, MD-consistent, and historically connected to the Case/Sitkoff style of field correction. Their weakness is that they are fixed force-field charges: conformation affects geometry, not charge redistribution. They also currently use cutoff machinery, which makes them operationally different from MOPAC.

MOPAC Mulliken charges are conformation-dependent and may capture electronic redistribution missing from FF14SB. Their weakness is that Mulliken charges are method-dependent population analyses, not observables, and the current implementation changes both charge source and pair range at once.

APBS is best treated as a screened/reaction-field channel: “FF14SB-like charge model plus continuum PB environment,” not a peer charge assignment method. Its strength is solvent physics; its current weakness is PB parameter validity, especially radii.

So the defensible design is parallel source channels plus reconciliation diagnostics. A “better-replaces-cruft” outcome is possible if one source dominates DFT and experiment after ablation. A “butterfly-friends” outcome is also plausible: agreement between MOPAC, FF14SB, and APBS-derived fields can be treated as validation of a robust electrostatic signal. The existing MOPAC-vs-FF14SB reconciliation cosine is already pointing in this direction.

**π-Quadrupole Partition**

This must be flagged, not resolved. The project’s own π-quadrupole note states the risk clearly: if charge/EFG already includes aromatic atom partial charges, then a point π-quadrupole can double-count the same far-field electrostatic multipole unless the partition is explicit [pi_quadrupole.md](/shared/2026Thesis/nmr-shielding/doc/emerging/kernel_design/pi_quadrupole.md:123). The current `PiQuadrupoleResult` computes a rank-2 EFG-like quadrupole tensor and scalar, not a full `1o ⊕ 2e` quadrupole field. Treat π-quad as a participant whose partition against charge/EFG is contested.

**Defensible Generous Emit Shape**

For each source channel, emit:

`E_1o`: `(N, 3)`, V/Angstrom, polar vector.  
`EFG_2e`: `(N, 5)`, V/Angstrom², real spherical symmetric traceless tensor.  
`0e` diagnostics: `|E|`, `|E|^2`, bond projection `E·b`, EFG Frobenius norm, eigenvalue/asymmetry summaries, source-partition magnitudes, cutoff/source-count diagnostics, and source agreement cosines.

Source channels should remain separate:

`ff14sb_vacuum_total`, plus backbone/sidechain/aromatic where available.  
`mopac_vacuum_total`, plus the same partitions.  
`apbs_pb_total`.  
`apbs_minus_ff14sb_vacuum` as a solvent/reaction-field delta when both are available.

I would not collapse source channels into one averaged field before the law study. Fido can take everything; the tensor model needs clean irreps and source tags; the law study needs ablatable source comparisons.

**Feed The Three New Steps**

Step 1, law/correlation study: use scalar and angular comparisons per source. Scalar-only baselines should include `|E|`, `E·bond`, `|EFG|`, EFG eigenvalue/asymmetry summaries, and source agreement metrics. Angular models should receive full `1o` and `2e` objects. This directly tests “with angular part vs without angular part” against DFT.

Step 2, equivariant shielding-tensor model: feed source-tagged `1o ⊕ 2e` channels plus `0e` diagnostics. The target shielding tensor remains even, so the model may consume `1o` field inputs but should emit shielding irreps as even tensor components, primarily `0e ⊕ 2e` unless antisymmetric shielding is deliberately modeled.

Step 3, experimental-shift fido predictor: include everything ablatable: raw field/EFG channels, current NPY scalar exports, APBS deltas, MOPAC/FF14SB agreement, source partitions, π-quad participant features, AIMNet2 charge/EFG features if available, and diagnostics marking cutoff/radii/source validity.

**Thin Spots**

The biggest thin spots are naming and comparability. Current `*_shielding.npy` files contain bare EFG tensors, not calibrated shielding. FF14SB and MOPAC differ in both charge model and cutoff/all-pair treatment. Sidechain EFG is computed but not cleanly exported like total/backbone/aromatic. APBS has explicit parameter-validity warnings and a clamp behavior the code itself calls a physics/blessing question. Static NPY outputs do not carry parity/unit metadata as cleanly as the newer H5 APBS trajectory outputs.

The honest scientific conclusion: the kernel form is defensible, the parity is clear, and the source fork is not decidable from theory alone. Keep the emit generous, keep the sources separate, and let DFT plus experimental ablations decide whether one source replaces the rest or whether agreement between independent electrostatic channels is itself useful signal.

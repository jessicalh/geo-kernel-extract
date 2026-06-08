# Charge / EFG — shareable spec

*The electrostatic field effect on nuclear shielding: the linear field (`1o`) and its gradient (`2e`).*

**Provenance.** Second of the Three to reach spec (McConnell was the template). Folds the first-stage
grounding (`efg_grounding_agent1.md`), the spec-vs-code reality-check (`efg_reality_check.md`), and the
2026-06-07 codex confirm pass. **Register: confirm-and-refine, not rescue.** The point-charge `q/r³` term
is our most reliable result; the job here is to state the standard, defensible form the field uses, check
our shape against it, and name honestly where the physics is clean and where it is not. Nothing below is
novel — it is the standard machinery, assembled and cited.

---

## 1. The wonder, and the clean maths

A charged group a few ångström from a nucleus — a backbone carbonyl, a Lys⁺, an ordered water — sits in
the same space as that nucleus's electron cloud, and tugs on it. The cloud is what sets the shielding, so
the tug shifts the shielding. This is a real *through-space* effect, distinct from the through-bond φ/ψ
dependence and from ring currents.

The standard framing is **Buckingham's (1960)**: expand the shielding tensor as a Taylor series in the
electric field **E** and its gradient **∇E** at the nucleus,

  σ_αβ(E) = σ⁰_αβ + σ_αβ,γ · E_γ + ½ σ⁽²⁾_αβ,γδ · E_γE_δ + σ⁽∇⁾_αβ,γδ · (∂E_γ/∂r_δ) + …

The **linear** coefficient σ_αβ,γ contracts the field **E** (`l=1`); the **quadratic** coefficient σ⁽²⁾
contracts **E⊗E**; the *distinct* **field-gradient** coefficient σ⁽∇⁾ contracts the **EFG ∇E** (`l=2`).
For a proton in an X–H bond the famous reduced form is **Δσ = −A·E_z − B·E²** — the field component *along
the bond* dominates (Buckingham 1960, *Can. J. Chem.* —
`references-meta/buckingham-1960-chemical-shifts-polar-groups-summary.txt`). The modern protein form makes
the EFG term explicit: **An 2014** (on GB3) fits **δ_HN(E) = δ_HN(0) − A·E − B·E² − C·EFG**, with A the
dipole shielding polarizability and C the EFG/quadrupole-shielding response — the held empirical anchor for
*both* the field and the EFG terms (`references-meta/an-2014-protein-apparent-dielectric-constant-summary.txt`).
(This electric-field expansion is the anchor we kept *here* rather than in McConnell, whose physics is
magnetic, not electric.)

Both objects come from a single function — the Coulomb potential `1/ρ` of a unit source — by one or two
spatial derivatives:

- **First derivative → the field.** `∂_a(1/ρ) = −ρ̂_a/ρ²`. A vector, falling as `1/ρ²`. An **`l=1`**
  quantity, and **odd under inversion** (position flips sign, charge doesn't, the field flips) → **`1o`**.
- **Second derivative → the EFG.** `∂_b∂_a(1/ρ) = (3 ρ̂_a ρ̂_b − δ_ab)/ρ³`. Symmetric, and **traceless**
  because `∇²(1/ρ) = 0` (Laplace, away from the source). So for charges *external* to the nucleus the EFG
  is **exactly a 5-component symmetric-traceless `l=2` object — `2e`, no `l=0`, no `l=1`, by
  construction, not approximation.** This is the single cleanest claim in the whole feature family.

The parity is the mirror of McConnell's, and the trap to avoid is importing McConnell's rule: shielding is
an even response tensor (`0e⊕1e⊕2e`, never `1o`); the **charge field is not a shielding tensor** — it is an
electrostatic *input* to shielding, and its field genuinely **is** odd (`1o`). (e3nn parity convention:
Geiger–Smidt 2022 — `references-text/geiger-smidt-2022-e3nn-euclidean-neural-networks-text-2.txt`. EFG as a
symmetric-traceless `l=2` spherical tensor is modern NMR-tensor practice: Ben-Mahmoud 2024 —
`references-text/ben-mahmoud-2024-gnn-solid-state-nmr-spherical-tensors-text-3.txt`.)

**So: charge/EFG = `E (1o) ⊕ V (2e)`.** Both multipole orders earn their place — the Buckingham expansion
couples to both, and they are not derivable from one another (different `l`, different distance weighting:
`E ~ q/r²`, `V ~ q/r³`).

---

## 2. What we compute

Per target nucleus `i`, per charge source, summed over a neighbourhood:

  E_a = k_e · Σ_j q_j · d_a / r³  
  V_ab = k_e · Σ_j q_j · (3 d_a d_b − δ_ab r²) / r⁵

with an **explicit traceless projection** after accumulation (to kill floating-point trace drift), then a
once-and-documented Cartesian→real-`l=2` conversion. The EFG's `T0` (trace) and `T1` (antisymmetric) are
**structural zeros** by §1 — only the 5-component `T2` carries signal (`_tensors.py:EFGTensor`, `1x2e`,
5-wide; the pre-2026-05-18 9-wide packing is rejected).

This is **target-based** physics: we sample the field *at* the nucleus. There is therefore **no
source-tensor to ride the molecular frame** (the contrast with ring/McConnell, whose unit group-tensor
rotates with the group — see the surgery guide, rule 1). The EFG **is** the physical tensor already; no
scale surgery, no `Δχ` analogue.

**Source channels — kept separate, the fork left emergent** (map, don't decide):

| source | E (`1o`) | EFG (`2e`) | range | note |
|---|---|---|---|---|
| FF14SB (`CoulombResult`) | yes | yes | **cutoff 20 Å** | the FullFat (`--mopac`) reconciliation source — PerFrame *skips* vacuum Coulomb; we run everything FullFat now, so it's present in practice (not unconditionally always-on in code). Cutoff left as-is. |
| MOPAC (`MopacCoulombResult`) | yes | yes | **all-pairs** | conformation-responsive Mulliken charges |
| APBS (`ApbsFieldResult`) | yes | yes | PB grid | screened/continuum field, not a peer charge set |
| water (`WaterFieldResult`) | yes | yes | explicit | first/second-shell separable; rides *with* charge/EFG |
| AIMNet2 (`AIMNet2Result`) | **no** | yes | cutoff | EFG + the learned embedding; **not** a field source today |

The FF14SB↔MOPAC range difference (cutoff vs all-pairs) is a caveat of the MOPAC-vs-FF14SB reconciliation
diagnostic, **not an EFG deliverable to harmonize**. The winning/default source is decided later by
ablation, not here.

---

## 3. The structural biology

The field term is a **leading backbone contribution**, not a correction. De Dios, Pearson & Oldfield
(1993, *Science*) established the **charge-field-perturbation** approach and showed that peptide ¹⁵N
shielding needs electrostatics alongside torsion, H-bonding, and neighbour identity — the amide nuclei are
genuinely field-sensitive (`references-meta/de-dios-1993-secondary-tertiary-structural-effects-summary.txt`).
**Long-range electrostatics** — the kind a coherent helix dipole produces — reach the backbone from beyond
the first shell, which is why a **generous cutoff (or all-pairs) matters for the field**: `E ~ q/r²` falls
off slowly. The EFG (`~ q/r³`) falls off faster, so a cutoff captures **most** of *it* — the physics that
argues for range lives in the `1o`, not the `2e`.

Polar and charged side chains (Asp/Glu carboxylates, Lys/Arg, His), backbone carbonyls, and ordered first-
shell water are the sources that matter; the empirical predictors (SHIFTS, SPARTA and successors) all
carry an electrostatic/field term beside ring-current and H-bond terms. (Case 1995 ring-current
calibration includes a Buckingham-style polarization term; Sitkoff 1997 peptide DFT; Sahakyan–Vendruscolo
2013 separate ring-current and field terms against DFT, field especially important for heavier nuclei —
all in `references-meta/`.)

---

## 4. What's clean, and what's not

**Clean (high confidence):**
- **The EFG as `2e`.** Symmetric and traceless by Laplace — exactly 5 components, no approximation. The
  *geometry* of the `q/r³` field gradient is unambiguous given the charges.
- **The field as `1o` and the `q/r²` law.** The directional, distance-weighted sum is exact for a given
  charge set. Our reliable point-charge result lives here — because *this* part is geometry, not the part
  where the uncertainty sits.

**Less clean (state the honest range, don't paper over it):**
- **The charge source.** The field's geometry is exact; the *charges* are a model choice. FF14SB (fixed
  RESP), MOPAC (responsive Mulliken), AIMNet2 (responsive Hirshfeld), QM/RESP give meaningfully different
  fields near polar/charged groups, and Mulliken/Hirshfeld populations are not observables. **No single
  right answer at this scale** — emit more than one, keep them separate, let the fit and the ablation
  weigh them. This *is* the emergent source fork.
- **Solvent screening / dielectric.** The biggest honest uncertainty. Fitting NMR chemical-shift
  perturbations to a single apparent dielectric gives **εₐ ≈ 8.6** (An 2014, on GB3), while a
  Poisson–Boltzmann protein interior wants **εₚ ≈ 2–4** — different numbers for "the" dielectric, because
  no single constant represents the protein/water boundary that does the screening. The honest range on
  the screened field's magnitude is **a factor of ~2–4** (An 2014's MD shows bulk water screening the
  direct field ≈4×). So the unscreened/effective-ε Coulomb field **and** the APBS PB field go in as
  **separate features**, and the gap between them is itself informative — never collapse to one "corrected"
  field. (An 2014 — `references-meta/an-2014-protein-apparent-dielectric-constant-summary.txt`; APBS =
  Poisson–Boltzmann continuum screening, Baker 2001 —
  `references-meta/baker-2001-apbs-electrostatics-nanosystems-summary.txt`.)
- **The quadratic `B·E²` term.** Real but second-order; whether it carries signal beyond the linear term
  at our field strengths is empirical. Cheap to expose as `0e` invariants (`|E|²`, EFG asymmetry η); no
  strong prior it dominates.
- **π-quadrupole partition.** Aromatic atom partial charges already enter the EFG; a point π-quadrupole
  can re-encode the same far-field multipole. The reality-check verdict: **not a simple double-count** —
  `PiQuadrupoleResult` builds a *distinct* geometric quadrupole basis that may carry π-anisotropy the
  point charges miss. **Decide by ablation, not a code-grounds drop.** π-quad is caged (a participant),
  so this is its ablation question, not an EFG obligation.

---

## 5. The danger-zone picture

Picture two nested regions around a backbone ¹⁵N:

- **The unambiguous core** — given a fixed charge set, the field `E` and gradient `V` at the nucleus are
  exact functions of geometry. Rotate the protein, both rotate correctly; the `2e` EFG is traceless to
  machine precision after projection. Nothing here is a judgment call.
- **The honest envelope around it** — two things blur the *magnitude*, not the geometry: (a) which charges
  you believe (FF14SB vs MOPAC vs AIMNet2 vs QM — divergent near the carboxylate, the guanidinium, the
  imidazole), and (b) how much the protein/water boundary screens the field on its way in (a factor of
  ~2–3, position- and direction-dependent, that no scalar dielectric captures). The picture to keep in
  the thesis is **a sharp inner kernel inside a fuzzy magnitude shell** — and the discipline that follows
  is: emit every charge source and both Coulomb-and-PB screenings as *separate* channels, and let the data
  draw the shell, rather than pre-choosing one field and hiding the width.

---

## 6. Fields, irreps, and the model

**The irreps, per quantity (loud, because the trap is real):**
- **E → `1o`** — one `l=1`, **odd** parity. The full vector in the molecular frame, *unprojected* (the
  bond projection `E·b̂` is the right physical scalar for a hand-built predictor, but inside e3nn it
  pre-imposes a frame and throws away equivariant content; let the network form the invariant by tensor
  product). **Do NOT import McConnell's all-even rule** — the field is genuinely odd.
- **EFG → `2e`** — one `l=2`, even, 5 real-spherical components. The clean object from §1.
- Optional `0e` invariants (`|E|²`, EFG η and largest-eigenvalue magnitude) — cheap, give the network easy
  access to the Buckingham `B·E²`/quadrupolar quantities, a sane hedge if message-passing is shallow.

**Per-Step altitude (the dragon):** EFG emits its `1o` field and `2e` gradient, but it is not a Step-1
input. Step 2 may use it only if ablation earns it. The EFG kernel can still be reported as a
**characterized descriptor↔DFT relationship** with confounds, not optimized as a model *input*.

**The surgery (confirm-and-refine, no new physics — the producer already computes the physics):**
1. **Naming.** `coulomb_shielding.npy` / `mopac_coulomb_shielding.npy` hold a **bare EFG** decomposition,
   not calibrated shielding. Rename to an EFG name (`*_efg`), so the label stops lying.
2. **Switch on the dropped channels, separably.** Sidechain EFG (computed, projected, *not written*) and
   the per-source partitioned **E** (only total emitted today) become separable emits — re-extract is
   <1 day, and you cannot conjure a channel you never wrote. Stay lean (<15 GB); drop downstream.
3. **Parity metadata.** The catalog marks E-fields `parity="odd"` but `irreps="1e"` — should be `1o`
   (`_catalog.py`). And the BS/HM trajectory-writer parity bug (`0e+1o+2e`) is a *ring* item, flagged in
   the ring spec — but the EFG catalog parity is fixed here.
4. **AIMNet2.** Contributes its **EFG (`2e`) and its embedding, not an E-field (`1o`)**. The missing `1o`
   is **not a gap** — the field is covered by FF14SB/MOPAC/APBS/water, AIMNet2's distinctive value is the
   embedding, and adding the (trivially-computable) `1o` is an **open item under the emergent fork**, not
   a fix. The old `charge_efg.md` "AIMNet2-primary" line is stale.
5. **SDK consumers.** New EFG emits need matching SDK entries: `CoulombGroup` / `MopacCoulombGroup`
   currently lack sidechain-EFG / partitioned-E fields, and `EFGTensor` rejects the old 9-wide packing —
   add the fields with the emits.

**Build-checks (mirror McConnell's):** rotation equivariance (`2e` by the `l=2` Wigner-D in the
real-spherical convention; `1o` as an odd vector); EFG `T0`/`T1` are structural zeros after projection;
traceless projection holds to machine precision; per-source channels stay separate; near-field (<3 Å)
accepted/rejected counts reported, not hidden.

**What this spec deliberately does NOT do:** pick the winning charge source (emergent); harmonize the
FF14SB/MOPAC range (not an EFG deliverable); touch FF14SB's status (stays); resolve the π-quad partition
(an ablation, and π-quad is caged); claim the tensor predicts the solution shift (the dragon — the `2e` is
DFT-angular evidence, never a solution-shift predictor).

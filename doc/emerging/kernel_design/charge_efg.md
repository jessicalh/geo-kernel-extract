Historical -- not current truth; see doc/emerging/kernel_design/efg_spec.md and doc/emerging/CONTROLLING_SPEC.md.

# Kernel design: the electrostatic contribution — electric field (l=1) and EFG (l=2)

Scope: one input feature for the e3nn equivariant shielding predictor — the
electrostatic effect of the surrounding charge distribution on the nuclear
magnetic shielding at a target nucleus. This is the **field-effect** family:
the linear electric field (a vector, l=1) and the electric-field gradient
(EFG, a rank-2 traceless symmetric tensor, l=2), framed by the Buckingham
field expansion.

Register: **confirm and refine.** This kernel is comparatively solid in our
hands — the point-charge q/r³ term is our most reliable result. The job here
is not rescue from scratch; it is to state the standard, defensible form the
field uses, check our shape against it, and name honestly where the physics
is clean and where it is not. "Comparatively solid" is not a licence to skip
scrutiny.

Terms are defined on first use. Nothing below is novel; it is the standard
machinery, assembled.

---

## 1. What the electrostatic effect on shielding is

The nuclear magnetic shielding σ at a nucleus is set by how the surrounding
electron cloud responds to an applied magnetic field. Any external
electric field perturbs that electron cloud, and so perturbs the shielding.
A charge sitting a few ångström from a nucleus — a backbone carbonyl, a
charged side chain, an ordered water — produces an electric field and an
electric-field gradient at that nucleus, and both shift its shielding. This
is a real, through-space physical effect, distinct from the through-bond
conformational dependence (φ/ψ) and from ring currents.

### The Buckingham field expansion

The standard framing is Buckingham's (1960): expand the shielding tensor of a
nucleus as a Taylor series in the electric field **E** and its gradient **∇E**
evaluated at the nucleus,

  σ_αβ(E) = σ_αβ⁰ + σ_αβ,γ E_γ + ½ σ_αβ,γδ E_γ E_δ + σ_αβ,γδ (∂E_γ/∂r_δ) + …

(summation over repeated Cartesian indices γ, δ). The coefficients are
intrinsic molecular response properties of the nucleus and its bonding
environment:

- **σ_αβ,γ** — the *dipole shielding polarizability*, a rank-3 tensor. It
  contracts with the **linear electric field E** (l=1). This is the leading,
  usually dominant, term.
- **σ_αβ,γδ** (the E_γE_δ term) — the *quadratic* field response, a rank-4
  tensor contracting with the outer product **E⊗E**. Smaller; matters in
  strong fields or for very polarizable bonds.
- the **field-gradient term** — contracts the EFG (l=2) with a corresponding
  response tensor.

For a proton in an X–H bond the famous reduced form is
**Δσ = −A·E_z − B·E²**, where E_z is the field component **along the bond**,
A and B are bond constants, and the linear A·E_z term dominates. This is the
"linear electric field effect" (LEFE) on the proton shift. The projection of
**E onto a bond axis** is the physically meaningful scalar invariant of the
l=1 part for that nucleus — it is what de Dios/Oldfield-style protein shift
predictors actually fit.

Sources:
[Buckingham 1960, *Can. J. Chem.* — chemical shift of molecules with polar groups](https://cdnsciencepub.com/doi/10.1139/v60-040);
[Nuclear shielding surface / field-expansion review (Springer)](https://link.springer.com/chapter/10.1007/978-94-011-1652-7_5);
[EFG effects on nuclear magnetic shieldings (academia.edu PDF)](https://www.academia.edu/32048535/Electric_field_gradient_effects_on_nuclear_magnetic_shieldings).

### How serious people compute it for proteins

The canonical protein line is **de Dios, Pearson & Oldfield (1993)** and the
"charge-field perturbation" / shielding-polarizability program that followed.
The recipe:

1. Take the protein's point charges (a force-field or QM-derived set).
2. Sum the **Coulomb field E** and the **EFG** at each backbone nucleus from
   all surrounding charges.
3. Multiply by precomputed *ab initio* shielding polarizabilities (the
   σ_αβ,γ tensors), or fit empirical coefficients, to get the electrostatic
   contribution to the shift.

They established (a) that the **peptide ¹⁵N and ¹HN nuclei are highly
polarizable** and therefore very sensitive to the local electric field, so
this term is a first-order effect for the backbone, and (b) that long-range
features such as the **α-helix macrodipole** measurably deshield ¹⁵N.
Empirical predictors (SHIFTS, SPARTA, and successors) carry an electrostatic
/ field term alongside ring-current and H-bond terms.

Sources:
[de Dios, Pearson, Oldfield 1993, *Science* — secondary/tertiary structure effects on protein shifts, ab initio charge-field approach](https://www.science.org/doi/10.1126/science.8502992);
[Calculation of dipole-shielding polarizabilities σ_αβγ for backbone nuclei (PubMed)](https://pubmed.ncbi.nlm.nih.gov/12903999/);
[DFT shielding of backbone ¹⁵N, helix-dipole electrostatics (Springer)](https://link.springer.com/article/10.1007/s10858-009-9358-3).

### The field itself: point charges, dielectric, Poisson–Boltzmann

Two standard ways to get **E** and **∇E** at the nucleus from the charges:

- **Direct point-charge Coulomb sum with an effective dielectric** ε_eff.
  E = k Σ_j q_j (r̂_j)/(ε_eff r_j²); the EFG is the second derivative of the
  same potential. Simple, fast, transparent. The catch is ε_eff: fitting NMR
  chemical-shift perturbations to charge changes in plastocyanin gives an
  optimal **ε_eff ≈ 6.5** for ¹⁵N — a single bulk number that crudely
  stands in for both protein and solvent screening.
- **Finite-difference Poisson–Boltzmann (PB)** — treat protein interior and
  solvent as distinct dielectric regions (protein ε_p, water ε_w ≈ 78) and
  solve for the field. The same NMR-CSP fit gives **protein ε_p ≈ 2–5,
  optimum ≈ 3**. PB does better than a single ε_eff precisely because
  **bulk-water screening strongly attenuates field propagation inside the
  protein**, and PB represents that boundary explicitly where a single
  dielectric constant cannot.

Source:
[Protein dielectric constants from NMR chemical-shift perturbations — Coulomb ε_eff≈6.5 vs PB ε_p≈3 (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4480918/);
[Probing electric fields in proteins by NMR (PubMed)](https://pubmed.ncbi.nlm.nih.gov/18214953/).

### Field is l=1, EFG is l=2 — and why the EFG is the clean one

The two orders are physically and mathematically distinct, and e3nn wants
them kept distinct:

- **E** is a vector → an **l=1** spherical-tensor (3 components, odd parity).
- **∇E** (the EFG, V_αβ = ∂²φ/∂r_α∂r_β) is a rank-2 Cartesian tensor. Its
  irreducible decomposition is l=0 (trace) ⊕ l=1 (antisymmetric part) ⊕ l=2
  (symmetric traceless part). But for a field produced by charges **external
  to the nucleus**, Laplace's equation ∇²φ = 0 forces the **trace to vanish**,
  and a gradient of a gradient is symmetric, so the antisymmetric l=1 part
  vanishes too. The EFG from external charge is therefore **exactly a 5-component,
  symmetric, traceless l=2 spherical tensor** — no l=0, no l=1 contamination,
  by construction. This is the cleanest statement available: the EFG is a
  pure l=2 object as a matter of electrostatics, not approximation.

Sources:
[Electric field gradient — symmetric, traceless by Laplace's equation, 5 independent components (Wikipedia)](https://en.wikipedia.org/wiki/Electric_field_gradient);
[Irreducible spherical-tensor decomposition of a rank-2 Cartesian tensor into k=0,1,2 (grad QM notes)](https://etneil.github.io/grad_qm_lec_notes/wigner_eckart.html).

---

## 2. Featurizing for e3nn — what the irreps should be

e3nn is an equivariant network: its features live in **irreducible
representations (irreps)** of O(3), labelled by degree l (dimension 2l+1) and
parity. Equivariance means a rotation of the input rotates the features by
the matching Wigner-D matrix; the network never sees a hand-imposed local
frame, and you must not give it Cartesian components dressed up as scalars.
The standard NMR-tensor ML practice (Venetos et al. on ²⁹Si shift tensors;
the zeolite/solid-state spherical-tensor GNNs) is to **carry quantities in
their spherical-tensor irreps and let Clebsch–Gordan tensor products combine
them**, decomposing any rank-2 tensor into its l=0/l=1/l=2 parts rather than
feeding nine Cartesian numbers.

Sources:
[Venetos et al. 2023 — equivariant GNN for full ²⁹Si shift tensors, spherical-tensor target (JPCA / PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10026072/);
[GNN prediction of solid-state NMR parameters via spherical-tensor decomposition (arXiv 2412.15063)](https://arxiv.org/abs/2412.15063).

The electrostatic feature is **at once a physical field and a lossy geometric
descriptor of the neighbourhood**, and that dual nature is exactly what makes
it a good e3nn input: the same q/r³ sum that estimates the physical field also
encodes the directional, distance-weighted arrangement of nearby charged
groups. Designed honestly in irrep terms:

- **E-field → `1o`** (one l=1, odd-parity irrep). The full vector at the
  nucleus, in the global/molecular frame, raw — no projection onto a local
  bond axis. (The bond projection E·b̂ is the right *physical* scalar for a
  hand-built predictor, but inside e3nn it is the wrong move: it pre-imposes
  a frame and throws away equivariant content. Let the network form the
  invariant from the `1o` field and the bond geometry via a tensor product.)
- **EFG → `2e`** (one l=2, even-parity irrep, 5 components). The symmetric
  traceless tensor mapped to its real l=2 spherical components. This is the
  clean object from §1.
- **Recommendation: feed both `1o ⊕ 2e`.** The electrostatic effect spans
  more than one multipole order — the field is l=1, the gradient is l=2 — and
  the Buckingham expansion couples to **both**. Feeding only one throws away a
  physically real, non-redundant order. They are not derivable from each other
  by the network (different l, different distance weighting: E ~ q/r², EFG ~
  q/r³), so both earn their place. Optionally add the **scalar invariants**
  |E|² (`0e`) and the EFG's rotational invariants (η, the asymmetry; and the
  largest-eigenvalue magnitude) as explicit `0e` channels — they cost nothing,
  give the network easy access to the quantities the Buckingham B·E² and
  quadrupolar terms actually use, and are a sane hedge if the higher-l message
  passing is shallow.

Per-source decomposition is a defensible refinement, not a requirement:
emitting separate `1o`/`2e` fields for **backbone**, **side-chain/charged**,
**aromatic**, and **solvent** charge subsets lets the network weight by
chemical origin (a carbonyl dipole field and a Lys⁺ monopole field have
different physics). It is more channels, not more assumptions; keep it if the
fit uses it, drop it if it doesn't.

---

## 3. The defensible, sane recommendation

**Compute, per target nucleus, two equivariant features and feed both:**

1. **Electric field `E`** as a `1o` irrep — the Coulomb vector field from
   surrounding charges within a neighbourhood cutoff, in the molecular frame,
   unprojected.
2. **EFG `V`** as a `2e` irrep — the symmetric traceless field-gradient
   tensor, in real l=2 spherical components.
   Plus, optionally, **`0e` invariants** |E|² and the EFG asymmetry/magnitude.

This is the standard Buckingham field expansion, sorted by multipole order
into the irreps e3nn consumes, with the network left to form the rotational
invariants (E·b̂, E:V, etc.) by tensor product rather than us pre-projecting.

**Charge source — most defensible choice, with a backup.**
- **Primary: AIMNet2 charges.** They are environment-responsive (they
  redistribute with conformation and local field), which is the physically
  correct behaviour for a through-space field term and matches the project's
  existing AIMNet2 dependency. This is the defensible default *for the field*.
- **Backup / reference: ff14SB fixed RESP charges.** Fixed, transferable,
  and the charges the MD that generated the geometry actually used — so the
  field is self-consistent with the trajectory. Cheap to emit alongside as a
  second channel; the fitter chooses.
- **Gold-standard anchor: QM/RESP** on a cluster, where affordable, as a
  positive control to check the cheaper charge fields against — not as a
  production feature.
  Emitting more than one charge source is honest: the charge set is the least
  clean input (see §4), so providing the alternative and letting the model
  weigh them beats betting the feature on one choice.

Source on charge sets:
[AIMNet2 — environment-responsive partial charges (PMC)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12057637/);
[Charge-derivation methods vs ff14SB for amino acids (ACS Omega)](https://pubs.acs.org/doi/10.1021/acsomega.8b00438).

**Point-charge vs continuum.**
- **Compute the explicit point-charge Coulomb field/EFG as the primary
  feature.** It is local, exact for a given charge set, fully equivariant,
  cheap with a KD-tree neighbourhood, and it is the term the field's
  shift-predictors are built on. This is the q/r³ "cookie."
- **Carry the APBS Poisson–Boltzmann field as a parallel feature**, because
  PB is the established way to represent **bulk-water screening**, which a raw
  vacuum sum or single-ε_eff sum gets wrong inside a protein. APBS already
  solves the field on a grid; sampling E and ∇E at each nucleus gives a
  second, screened `1o`/`2e` pair. The honest framing: vacuum/effective-ε
  Coulomb and PB are **two estimates of the same physical field under
  different screening assumptions**, and the gap between them is itself
  informative. Feed both; do not collapse to one.
- **Solvent (explicit water) field** is a third, legitimate channel: the
  field from ordered/first-shell waters that the continuum smears out. Keep
  it as its own `1o`/`2e` source rather than folding it into the protein sum.

**Use the fast linear algebra where it buys defensibility, not for its own
sake:** KD-tree neighbourhoods so the sum is honest about its cutoff; direct
EFG-tensor assembly V_αβ = k Σ_j q_j (3 d_α d_β − δ_αβ r²)/r⁵ with an explicit
traceless projection after accumulation to kill floating-point trace drift;
Cartesian→real-l=2 conversion done once with a fixed, documented basis;
spherical-harmonic sampling of the PB grid. None of this is novel — it is
careful bookkeeping that makes the irreps clean.

---

## 4. Where the physics is clean, and where it is not

**Clean (high confidence):**
- **The EFG as an l=2 object.** For charges external to the nucleus, the EFG
  is symmetric and traceless by Laplace's equation — exactly 5 components,
  exactly `2e`, no approximation. This is the single cleanest claim in the
  whole feature family. The *geometry* of the q/r³ field gradient is
  unambiguous given the charges.
- **The field as l=1 and the q/r² law.** The directional, distance-weighted
  Coulomb sum is exact for a given charge set. Our point-charge q/r³ result is
  reliable because *this* part — the geometry — is not where the uncertainty
  lives.

**Less clean (state the honest range, don't paper over it):**
- **The charge source.** The geometry of the field is exact; the *charges*
  are a model choice. ff14SB (fixed RESP), AIMNet2 (responsive), and QM/RESP
  give meaningfully different fields, especially near polar/charged groups.
  There is no single right answer at this scale; the defensible move is to
  emit more than one and let the model weigh them. **I don't know which is
  best for shielding specifically — that is an empirical question for the fit,
  and the honest range is "ff14SB to AIMNet2 to QM, with QM as the anchor."**
- **Solvent screening / dielectric.** The biggest honest uncertainty. The
  same NMR data are fit by Coulomb with **ε_eff ≈ 6.5** *or* PB with
  **ε_p ≈ 3** — two different numbers for "the" dielectric, because a single
  constant cannot represent the protein/water boundary that actually does the
  screening. Inside a real protein the effective screening is **position- and
  direction-dependent**, somewhere between vacuum and bulk water, and no
  scalar captures it. PB is the better-grounded representation; it is still a
  continuum approximation that smears explicit ordered water. The honest range
  on the screened field's magnitude is **a factor of ~2–3** depending on the
  screening model — which is exactly why both the unscreened/effective-ε
  Coulomb field and the PB field should go in as separate features rather than
  one "corrected" field.
- **The quadratic (B·E²) and explicit field-gradient response terms.** Real
  but second-order; whether they carry signal beyond the linear term at our
  field strengths is empirical. Cheap to expose as `0e` invariants; no strong
  prior that they dominate.

---

## 5. What we currently compute vs the recommendation (neutral)

We already compute, per atom, both orders of this feature: a vacuum
point-charge **electric field E** (`CoulombResult`, "vacuum Coulomb E-field")
and a point-charge **EFG tensor V** that is explicitly assembled with the
dipolar kernel and **traceless-projected** after accumulation, then decomposed
to a spherical (l=2) tensor — and the code already documents E as the rank-1 /
Buckingham-A,B leg and V as the rank-2 / γ leg, kept as two separate kernels
rather than one unified tensor. We compute an **APBS Poisson–Boltzmann** field
sampled on a grid (the continuum/screened estimate), a **WaterFieldResult**
from explicit solvent (with first/second-shell separation), a **MOPAC/AIMNet2
charge** path alongside the ff14SB charges, and rediscover-side feature cells
(`BuckinghamEfield`, `EfgFeature`, `ChargeDipoleNeighborhood`) that already
project/emit these quantities. In short: the pieces the recommendation names —
both multipole orders, point-charge and PB and explicit-solvent fields,
multiple charge sources, the traceless l=2 EFG — already exist in the producer;
the recommendation is about which of them to feed e3nn and **in what form**
(raw `1o` field and `2e` EFG irreps in the molecular frame, both orders
together, multiple charge/screening channels left for the fit to weigh), as
opposed to pre-projected scalars in an imposed local frame. This paragraph is
descriptive only; it carries no verdict on the current code.

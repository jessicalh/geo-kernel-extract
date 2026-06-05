# Maxwell → Shielding — a step-by-step derivation
### Prospectus (the plan; the document itself is built next, with Mathematica)

## The idea, in one breath
The physics-architecture tower (`../thesis-overview/fig_shielding_to_maxwell.pdf`) *asserts* a
chain: **Maxwell → 1/ρ → D_ab → fork (magnetic / electrical) → the calculators → shielding.**
This document turns each arrow of that tower into an actual derivation — so the classical
calculators stop being formulas we trust and become *the only thing they could have been*: named
projections of one kernel that is the near-field shadow of Maxwell. The figure is the map; this
document walks every rung.

## Why we are doing it
- **Grounding.** It answers the deepest "why these calculators?" — not by appeal to the literature,
  but by derivation. Each calculator is a specific instance of one descent from first principles.
- **It is the bedrock the categorical-grounding chapter reaches for.** That chapter grounds the
  calculators in history and atom-type; this grounds them in Maxwell.
- **It is honest.** If the chain holds at every level (it does, on inspection), the thesis's central
  claim — *one geometry behind all the calculators* — is a theorem we worked, not a slogan.

## How it gets made
- **Mathematica** (on the Windows machine, where it lives) is the engine and the check: it does the
  symbolic algebra and *verifies* each step — computes `∂_a∂_b(1/ρ)`, confirms it equals
  `(3ρ̂ρ̂−δ)/ρ³`, gets the eigenvalues, shows the magnetic-dipole field tensor and the
  electric-field-gradient are the *same* `D_ab`, runs the multipole expansion, finds the magic-angle
  root, plots the angular function.
- **Claude** structures the document, narrates each step in the thesis voice (every symbol grounded,
  no hand-waving), ties each derivation to the actual calculator code, holds the "why," and catches
  convention/sign errors before they propagate.
- **Jessica** drives — the physics judgment, the voice, what depth is enough, what stays.

## The chain — one section per rung (each: *from → how → to → meaning → which calculator*)

1. **Maxwell's source equations → the potentials.**
   `∇·E = ρ/ε₀` (Gauss) and `∇×B = μ₀J` (Ampère) → `E = −∇φ`, `B = ∇×A` → Poisson/Laplace for φ, A.

2. **→ the 1/ρ Green's function.**
   Solve `∇²G = −4πδ(ρ)` → `G = 1/ρ`. The single object both φ and A are built on. *(This is the
   waist: where electric and magnetic become one geometry.)*

3. **→ D_ab.**
   The field of a localized source = derivatives of `1/ρ`. Second derivative →
   `D_ab = ∂_a∂_b(1/ρ) = (3ρ̂_aρ̂_b − δ_ab)/ρ³`. **Show explicitly that the magnetic dipole field
   `B_a = D_ab m_b` and the electric field gradient `∂_a E_b = D_ab q` are the same tensor.**

4. **→ the geometry of D_ab.**
   Traceless, symmetric; eigenvalues `+2/ρ³` (along ρ̂), `−1/ρ³` (×2); the angular law
   `(3cos²θ−1)/ρ³`; the magic-angle node at `cos²θ = 1/3` (54.7°); radial `1/ρ³`.

5. **→ the fork: the two doorways and their calculators.**
   - **Magnetic** — `D_ab` contracted with an induced moment `m` (ring current χ; bond Δχ):
     → **McConnell, Johnson–Bovey, Haigh–Mallion, RingSusceptibility, BiotSavart, bond anisotropy.**
   - **Electrical** — `D_ab` contracted with a charge `q` (the EFG / field):
     → **Coulomb, APBS, MOPAC, EEQ, WaterField, π-quadrupole.**
   For each: what physical quantity it computes, and *which step of this chain it implements*.

6. **→ back to shielding.**
   Ramsey: the shielding tensor is the response of the electron currents to `B₀`. The through-space
   pieces enter σ by two routes — magnetic mechanisms **add a field at the nucleus**; electric
   mechanisms **perturb the electronic response**. Each calculator computes its piece of σ; the sum
   `σ_ab(i) ≈ Σ_s I_s D_ab(r_is)` is the closing of the loop.

7. **The calculator level — what we are looking at, and why.**
   The synthesis: a table mapping every calculator to its rung, its source, and its intensity `I_s`,
   so a reader sees the whole library as one derivation instantiated many times.

## The eventual document's shape
- **LaTeX** (`maxwell-to-shielding.tex`), one section per rung above; each section = *starting point
  → derivation (key Mathematica steps shown) → result → physical meaning → calculator(s) grounded →
  link to the next rung.*
- The tower figure as the frontispiece/map; each section fills one rung of it.
- **Mathematica notebooks** as the working + verification layer (likely an appendix or a sibling
  `notebooks/` folder); the LaTeX shows the *result and the meaning*, the notebooks hold the *grind
  and the proof*.

## Open decisions (for when we sit down on the Windows machine)
- **Units/conventions:** SI vs Gaussian; the sign/normalisation choices that bite (flag them once,
  up front).
- **Depth:** full derivations, or key steps with the algebra in the notebooks?
- **Scope of calculators:** all of them, or one representative per branch with the rest tabulated?
- **Format split:** how much lives in LaTeX vs the Mathematica notebooks; how the notebooks export.
- **Frame convention:** keep raw molecular-frame throughout (consistent with the equivariant work).

## Handoff
The derivation + Mathematica work happens on the **Windows machine** (Mathematica's home), with
Claude there. This prospectus and the three figures (`fig_shielding_to_maxwell`,
`fig_governing_equation`, `fig_shielding_mag_elec`) travel with it as the map and the plan. This
subdirectory holds the prospectus now, the LaTeX as it's built, and the exported notebook work.

*(Directory name `maxwell-to-shielding/` is a working choice — easily renamed.)*

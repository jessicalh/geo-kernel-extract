# Equivariant methodology and field survey for tensor-resolved NMR shielding

Date: 2026-06-06

This note grounds the proposed `e3nn`-style shielding-tensor model and gives an honest field survey for the current thesis plan in `doc/emerging/CONTROLLING_SPEC.md`. I have treated that controlling spec as the current design and ignored older "three stage" material where it conflicts with it.

The plan being surveyed is:

1. Compute physically motivated per-atom or per-neighbour tensor features from geometry: ring-current fields, Haigh-Mallion ring-current scalars, McConnell neighbour magnetic-anisotropy tensors, and charge/electric-field-gradient features.
2. Express these features in typed equivariant irreducible representations, chiefly `0e + 1e + 2e`.
3. Feed them, together with geometry, into an equivariant graph model that predicts the NMR magnetic shielding tensor.
4. Calibrate or train against DFT shielding/EFG data.
5. Apply the predictor per frame over MD trajectories and average, while stating that finite MD does not fully sample the experimental motional/ensemble average.

The short verdict is: the linear-algebra and equivariant-ML methodology are recognised ground; tensorial NMR prediction with equivariant models is real but still a small, materials-heavy subfield; the full pipeline of hand-computed NMR-physics tensor kernels into an equivariant tensor predictor and then MD averaging is genuinely unusual. I found close neighbours, but no paper that does this whole combination in the biomolecular shielding setting.

## Executive verdict

| Question | Honest classification | Evidence |
|---|---:|---|
| Equivariant atomistic ML on molecular geometry | Recognised ground | Tensor Field Networks, Cormorant, `e3nn`, NequIP, MACE, Allegro, SEGNN, SE(3)-Transformer. |
| Equivariant prediction of NMR shielding/EFG tensors | Quiet but real | Venetos et al. predict full 29Si shielding tensors in silicates; Ben Mahmoud et al. predict shielding and EFG tensors in silica; recent materials work extends this to EFGs and tensorial NMR metrics. |
| Physical or computed tensors as equivariant input features | Architecturally recognised, NMR-specific use unusual | SEGNN explicitly allows physical vector/tensor node and edge attributes. Protein shift predictors use physical scalar terms. iShiftML uses cheap QM shielding tensors as ML features. I did not find the specific pattern "ring-current/McConnell/EFG tensors as `0e + 1e + 2e` input channels to an e3nn shielding tensor model." |
| MD-trajectory averaged scalar chemical shifts | Recognised ground | Li and Bruschweiler, Robustelli et al., de Gortari et al., PPM, SPARTA+/SHIFTX2 workflows. |
| MD-trajectory averaged equivariant tensor prediction | Emerging, mostly materials | Ben Mahmoud et al. apply tensorial NMR ML to cristobalite dynamics; Zakary and Lantto use Allegro+MatTen to predict 129Xe shielding tensors from MLMD snapshots and average isotropic shifts; Schmiedmayer et al. do MLMD plus equivariant EFG prediction for MAPbI3. |
| Whole proposed biomolecular pipeline | Genuinely unusual | The pieces exist, but I found no direct published analogue combining hand physical tensor kernels, e3nn tensor learning, DFT calibration, and MD-ensemble averaging for biomolecular shielding tensors. |

The safest thesis phrasing is therefore not "we follow a standard published pipeline". It is closer to:

> We combine established equivariant tensor-learning machinery with established NMR physical descriptors and established trajectory averaging. The particular integration of precomputed NMR-physics tensor kernels as typed equivariant features for a DFT-calibrated shielding tensor predictor appears to be unusual and, to our knowledge from the surveyed literature, not a standard published recipe.

That is a defensible MSc-methods novelty if it is framed as integration and tested with ablations. It becomes a red flag only if the thesis implies that the whole pipeline is already routine, if the DFT reference is treated as experimental truth, or if MD averaging is presented as fully sampled experimental motional averaging.

## Part A: equivariant methodology

### What equivariance means here

For an atomistic structure with positions `r_i`, atom types, and optional typed physical features, an equivariant model should respond predictably to rotations, translations, inversion/reflection, and atom permutations.

For a group element `g` in `E(3)` or `O(3)`, a model `F` is equivariant when

```text
F(D_in(g) x) = D_out(g) F(x).
```

If the output is a scalar energy, `D_out(g)` is the identity. If the output is a vector, it rotates as a vector. If the output is a rank-2 shielding tensor, it transforms as

```text
T -> R T R^T
```

under a rotation `R`. Translational invariance is usually enforced by using relative positions `r_ij = r_j - r_i`; permutation equivariance is enforced by shared message functions and symmetric aggregation over neighbours.

The modern equivariant atomistic-ML literature uses this principle explicitly. Tensor Field Networks introduced rotation- and translation-equivariant convolutions over 3D point clouds using spherical harmonics and Clebsch-Gordan tensor products [Thomas et al. 2018](https://arxiv.org/abs/1802.08219). Cormorant built covariant molecular neural networks with similar representation-theoretic machinery [Anderson et al. 2019](https://arxiv.org/abs/1906.04015). The `e3nn` library generalises these operations for Euclidean neural networks and typed geometric tensors [Geiger and Smidt 2022](https://arxiv.org/abs/2207.09453). NequIP, MACE, Allegro, and SEGNN are prominent descendants or neighbours in atomistic modelling [Batzner et al. 2022](https://www.nature.com/articles/s41467-022-29939-5), [Batatia et al. 2022](https://arxiv.org/abs/2206.07697), [Musaelian et al. 2023](https://www.nature.com/articles/s41467-023-36329-y), [Brandstetter et al. 2022](https://arxiv.org/abs/2110.02905).

### Irreps, parity, and Wigner-D matrices

For 3D rotations, the irreducible representations are indexed by angular degree `ell = 0, 1, 2, ...`. An `ell` irrep has dimension `2 ell + 1`.

Common cases:

| irrep | dimension | informal object |
|---|---:|---|
| `0` | 1 | scalar |
| `1` | 3 | vector-like object |
| `2` | 5 | traceless symmetric rank-2 object / quadrupole-like object |

For `O(3)` rather than just `SO(3)`, each irrep also has a parity label. In `e3nn` notation this is usually `e` for even and `o` for odd. Inversion acts as `+1` on even features and `-1` on odd features. A displacement vector is `1o`: it changes sign under inversion. A rank-2 tensor formed from two polar vectors has even parity, because `1o x 1o` gives even parity.

For a rotation `R`, a feature of type `ell` transforms by the Wigner-D matrix:

```text
x_m' = sum_n D^ell_{m n}(R) x_n,    m,n = -ell,...,+ell.
```

In a real-valued implementation such as `e3nn`, the basis is usually a real spherical-harmonic basis rather than the complex `Y_ell^m` basis. The representation is equivalent; the component ordering and normalisation differ.

For a general `O(3)` element written as a rotation plus optional inversion, the parity label supplies the additional sign. Thus a `1o` vector changes sign under inversion, while a `1e` axial vector does not.

### Rank-2 Cartesian tensor decomposition: `0e + 1e + 2e`

Let `T` be a real Cartesian rank-2 tensor. For NMR shielding, this is a 3 by 3 tensor. Under rotation,

```text
T -> T' = R T R^T.
```

Decompose `T` into three orthogonal Cartesian pieces:

```text
T = T_iso I + A + Q
```

where

```text
T_iso = tr(T) / 3
A     = (T - T^T) / 2
S     = (T + T^T) / 2
Q     = S - T_iso I
```

The trace `T_iso` is invariant under `T -> R T R^T`, because `tr(R T R^T) = tr(T)`. It is the `ell = 0` even scalar component, `0e`.

The antisymmetric part `A` has three independent components and can be encoded as an axial vector:

```text
a = [
  (T_yz - T_zy) / 2,
  (T_zx - T_xz) / 2,
  (T_xy - T_yx) / 2
].
```

This vector-like object transforms by the `ell = 1` rotation matrix. Because it arises from the product of two polar directions, it is even under inversion: `1e`. In most diamagnetic shielding applications the antisymmetric part is often small or not the experimentally dominant term, but for a full tensor target it is the correct slot for it.

The remaining symmetric traceless tensor `Q` has five independent components:

```text
Q = Q^T,     tr(Q) = 0.
```

One convenient real `ell = 2` encoding, matching the style already used in the project docs, is:

```text
q = [
  sqrt(2)       Q_xy,
  sqrt(2)       Q_yz,
  sqrt(3 / 2)   Q_zz,
  sqrt(2)       Q_xz,
  (Q_xx - Q_yy) / sqrt(2)
].
```

Normalisation conventions can vary; what matters is that the encoding is fixed consistently for training and reconstruction.

Why this is the `ell = 2` piece: associate each traceless symmetric `Q` with the homogeneous quadratic polynomial

```text
p_Q(x) = x^T Q x.
```

The Laplacian is

```text
Delta p_Q = 2 tr(Q) = 0,
```

so `p_Q` is a harmonic polynomial of degree 2. The five-dimensional space of degree-2 harmonic polynomials is exactly the `ell = 2` irreducible representation of `SO(3)`. Under rotation,

```text
p_{R Q R^T}(x) = p_Q(R^T x),
```

which is the degree-2 Wigner-D action in a Cartesian disguise. This is the concrete reason a symmetric traceless tensor belongs in a `2e` channel.

Thus a completely general rank-2 polar-polar response tensor decomposes as:

```text
rank-2 tensor = 0e + 1e + 2e.
```

For an electric field gradient in a charge-free local region, the EFG tensor is symmetric and traceless, so the natural target is only `2e`. For shielding, the full tensor can in principle contain all three pieces; many practical NMR observables emphasise the scalar trace and symmetric anisotropy.

### Tensor products and Clebsch-Gordan coupling

Equivariant networks need a way to combine typed features without breaking symmetry. The key rule is the Clebsch-Gordan decomposition:

```text
(ell_1, p_1) tensor (ell_2, p_2)
  -> direct sum over ell = |ell_1 - ell_2|, ..., ell_1 + ell_2
     of (ell, p_1 p_2).
```

In components, using a complex spherical basis for clarity,

```text
z^{ell_3}_{m_3}
  = sum_{m_1,m_2}
      C^{ell_3 m_3}_{ell_1 m_1, ell_2 m_2}
      x^{ell_1}_{m_1}
      y^{ell_2}_{m_2},
```

where the `C` coefficients are Clebsch-Gordan coefficients. Real `e3nn` bases use equivalent real coupling matrices.

The concrete case that matters for molecular geometry is two polar vectors:

```text
1o tensor 1o -> 0e + 1e + 2e.
```

For two ordinary vectors `u` and `v`, these pieces are recognisable:

```text
0e:  s = u dot v

1e:  a = u cross v

2e:  Q_uv = (u v^T + v u^T) / 2 - (u dot v / 3) I
```

This is exactly the algebra behind many physical kernels. A point-dipole or Hessian-like propagator uses `r_hat r_hat` minus a trace term, which is a `2e` object. A full non-symmetric rank-2 construction can also produce the `1e` antisymmetric slot.

### Steerable convolutions and equivariant message passing

In an equivariant message-passing network, node features are not just scalar vectors. They are direct sums of irreps:

```text
f_i = f_i^{0e} + f_i^{1o} + f_i^{1e} + f_i^{2e} + ...
```

For an edge `i <- j`, the relative direction `r_hat_ij` is expanded in spherical harmonics:

```text
Y_L(r_hat_ij)
```

which itself transforms as an irrep of degree `L` with parity `(-1)^L`. A typical steerable message has the form:

```text
m_ij^{ell_out, p_out}
  = sum_paths W_path(|r_ij|)
      [ f_j^{ell_in, p_in} tensor Y_L(r_hat_ij) ]_{ell_out, p_out}.
```

The radial weights `W_path(|r_ij|)` are learned scalar functions of distance. The bracket means "project the tensor product onto the requested output irrep using Clebsch-Gordan coefficients." Messages are summed over neighbours:

```text
M_i = sum_{j in N(i)} m_ij.
```

Because each ingredient transforms predictably, the sum also transforms predictably. Nonlinearities are usually applied to scalar channels directly and to higher-order channels through gates or norm-based mechanisms so that equivariance is preserved.

This is enough linear algebra for the proposed thesis purpose: the model is equivariant because every operation is built from interatomic relative geometry, spherical harmonics, radial scalar weights, Clebsch-Gordan projections, and symmetric neighbour aggregation.

### Canonical frameworks and what each contributes

Tensor Field Networks are the early direct reference for spherical-harmonic filters and tensor products on 3D point clouds [Thomas et al. 2018](https://arxiv.org/abs/1802.08219).

Cormorant demonstrates covariant molecular neural networks using irreducible representations and Clebsch-Gordan coupling for molecular properties [Anderson et al. 2019](https://arxiv.org/abs/1906.04015).

SE(3)-Transformer adds equivariant attention over 3D graphs [Fuchs et al. 2020](https://arxiv.org/abs/2006.10503). It is not essential to the present plan, but it is part of the same architectural family.

`e3nn` is the practical library anchor: it represents features as direct sums of `O(3)` irreps and supplies tensor products, spherical harmonics, gates, and layers for Euclidean neural networks [Geiger and Smidt 2022](https://arxiv.org/abs/2207.09453).

NequIP uses E(3)-equivariant message passing to learn interatomic potentials with strong data efficiency. It is an important precedent for atomistic models whose internal features are geometric tensors even when the final output is scalar energy [Batzner et al. 2022](https://www.nature.com/articles/s41467-022-29939-5).

MACE uses higher-order equivariant message passing inspired by many-body/atomic-cluster ideas and is now a standard reference for accurate equivariant force fields [Batatia et al. 2022](https://arxiv.org/abs/2206.07697).

Allegro is a strictly local E(3)-equivariant architecture designed for scalable atomistic dynamics, used in large MLMD settings [Musaelian et al. 2023](https://www.nature.com/articles/s41467-023-36329-y).

SEGNN is especially relevant to the present thesis because it explicitly allows geometric and physical node/edge attributes beyond invariant scalars, including vectors and tensors [Brandstetter et al. 2022](https://arxiv.org/abs/2110.02905). That is the clearest architecture-level precedent for "feed physical equivariant quantities into the network." It is not, by itself, a precedent for the specific NMR shielding kernels.

## Part B: field survey

### Search basis and negative evidence

I searched the local paper corpus through `references-meta/INDEX.md` and `references-text/`, using variants of:

- `equivariant`, `e3nn`, `NequIP`, `MatTen`, `tensor`, `spherical tensor`
- `NMR shielding tensor`, `chemical shift tensor`, `EFG`, `electric field gradient`
- `ring current`, `magnetic anisotropy`, `McConnell`, `electric field`, `physical feature`
- `MD averaged chemical shift`, `trajectory chemical shift`, `ensemble shift`

I also searched the web for current literature through 2026-06-06, including recent 2025 publications. Negative evidence is never proof of absence, but the pattern was consistent: the exact full pipeline did not appear in the corpus or broad web search, while several partial neighbours did.

### Equivariant prediction of NMR shielding and EFG tensors

This is no longer empty territory, but it is not a large biomolecular field.

The clearest early full-shielding-tensor neighbour is Venetos, Wen, and Persson, "Machine Learning Full NMR Chemical Shift Tensors of Silicon Oxides with Equivariant Graph Neural Networks" [J. Phys. Chem. A 2023](https://pubs.acs.org/doi/10.1021/acs.jpca.2c07530); the open PMC copy summarises that their equivariant GNN learns full 29Si NMR shielding/shift tensors in silicates and outperforms invariant eigenvalue approaches [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC10026072/). This is a strong precedent for tensor-resolved NMR prediction, but it is a materials/silicate model, not a protein model, and it learns from geometry rather than from hand-computed NMR physical kernels.

The held anchor, Ben Mahmoud, Rosset, Yates, and Deringer, is directly on point for spherical-tensor decomposition. Their 2024 preprint and 2025 JCP article systematically study graph-neural-network prediction of solid-state NMR tensor quantities: anisotropic magnetic shielding and the electric field gradient, with downstream chemical shifts, quadrupolar couplings, tensor orientations, spectra, large-cell transfer, and cristobalite dynamics [arXiv:2412.15063](https://arxiv.org/abs/2412.15063), [JCP 163, 024118, DOI 10.1063/5.0274240](https://www.materials.ox.ac.uk/publication/2245400/dimensions). This is the closest methodological anchor for "spherical tensor decomposition + equivariant GNN + NMR tensor targets."

ShiftML3 extends deep learning of organic-solid shieldings to full shielding tensors and anisotropy at large scale [Kellner et al. 2025](https://arxiv.org/abs/2506.13146), with an online description noting full chemical shielding tensor prediction over a large organic-crystal training set [ShiftML3 predictor](https://shiftml.dev.materialscloud.cscs.ch/). From the local summary, ShiftML3 is not best described as a strict e3nn irrep pipeline; it is still important because it shows that full tensor shieldings, not only isotropic shifts, are becoming a practical ML target.

Recent EFG-specific materials work also matters. Schmiedmayer et al. combine a machine-learned force field with a second symmetry-preserving model for electric field gradients to predict quadrupolar coupling constants in the MAPbI3 phase transition [arXiv:2507.19435](https://arxiv.org/abs/2507.19435). This is not shielding, but EFG is a rank-2 tensor NMR observable and it is close to the dynamics/tensor side of the thesis.

The conclusion for tensorial NMR is: there is a real, citable subfield, but it is small and concentrated in solid-state/materials examples. For biomolecular chemical shielding tensors, especially protein or structural-biology settings, I did not find a comparably populated literature using e3nn-style tensor prediction.

### Physical or computed tensors as equivariant input features

This is the part where the proposed thesis is most unusual.

There is strong general precedent for feeding physical equivariant quantities into an equivariant GNN. SEGNN explicitly generalises message passing so node and edge attributes can include covariant physical information such as position, force, velocity, spin, vectors, and tensors [Brandstetter et al. 2022](https://arxiv.org/abs/2110.02905). This supports the representation-theoretic legitimacy of supplying `0e`, `1o`, `1e`, and `2e` physical fields as model inputs.

There is also strong NMR precedent for physical descriptors, but mostly in scalar or semiempirical models rather than equivariant tensor channels. Case's review of biomolecular chemical shifts discusses physical contributions such as local electronic structure, aromatic/ring-current effects, magnetic anisotropy, hydrogen bonding, electric fields, and conformational dependence [Case 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3877577/). SPARTA uses torsion/sequence matches plus neighbour, ring-current, and hydrogen-bond corrections [Shen and Bax 2007](https://colab.ws/articles/10.1007%2Fs10858-007-9166-6). SPARTA+ and SHIFTX2 are scalar shift predictors with machine-learning and structural descriptors [Shen and Bax 2010](https://pubmed.ncbi.nlm.nih.gov/20628786/), [Han et al. 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3085061/). UCBShift 2.0 is a recent high-performing scalar protein chemical-shift model [Ptaszek et al. 2024](https://pubs.acs.org/doi/10.1021/jacs.4c10474).

iShiftML is a useful intermediate case: Li et al. propose a feature representation informed by cheap low-level QM calculations of atomic shielding tensors, then train to high-level shielding targets [arXiv:2306.08269](https://arxiv.org/abs/2306.08269). That is a precedent for "computed shielding-like tensors as features", but it is not the same as typed equivariant tensor features in an e3nn graph model and is not the ring-current/McConnell/EFG-kernel construction proposed here.

What I did not find: papers that compute NMR-inspired physical field tensors such as ring-current tensors, McConnell magnetic-anisotropy tensors, or EFG-like charge tensors and feed them explicitly as `0e + 1e + 2e` feature channels to an equivariant GNN predicting shielding tensors.

So the honest classification is mixed:

- Recognised at the architecture level: yes, physical steerable attributes are a known pattern.
- Recognised in NMR scalar predictors: yes, physical descriptors are long-established.
- Recognised as this exact NMR tensor-kernel/e3nn-feature pattern: I found little to none. This part is genuinely unusual.

### MD trajectory averaging and dynamics

For scalar chemical shifts, trajectory or ensemble averaging is established. Li and Bruschweiler used ensemble-averaged back-calculated shifts to assess MD force fields and found improvements over individual snapshots [J. Phys. Chem. Lett. 2010](https://pubs.acs.org/doi/10.1021/jz9001345). Robustelli, Stafford, and Palmer evaluated dynamically averaged backbone shifts from unbiased MD simulations and compared them with experimental shifts [JACS 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3324661/). De Gortari et al. studied time averaging of NMR chemical shifts in the MLF peptide using DFT and MD/conformer sampling [JACS 2010](https://pubs.acs.org/doi/abs/10.1021/ja9062629). PPM was explicitly designed for assessing protein conformational ensembles [Li and Bruschweiler 2012](https://ohiostate.elsevierpure.com/en/publications/ppm-a-side-chain-and-backbone-chemical-shift-predictor-for-the-as/). MDTraj also exposes wrappers for SPARTA+, PPM, and SHIFTX2 over trajectories, which reflects how normal the scalar workflow has become [MDTraj NMR example](https://www.mdtraj.org/1.3.0/examples/nmr.html).

For tensor-resolved equivariant NMR prediction over trajectories, the literature is much thinner but not empty.

Ben Mahmoud et al. apply their tensorial NMR ML models to the dynamics of the alpha-to-beta inversion in cristobalite [arXiv:2412.15063](https://arxiv.org/abs/2412.15063), and the 2025 JCP record explicitly lists static and dynamic behaviour of complex materials as part of the scope [Oxford/JCP record](https://www.materials.ox.ac.uk/publication/2245400/dimensions).

Zakary and Lantto are a very close dynamics neighbour for the workflow shape. They use E(3)-equivariant GNNs to build both an ML interatomic potential and an NMR-ML model for 129Xe in porous liquids. The article states that Allegro enables MLMD and MatTen predicts the 129Xe magnetic shielding tensor from MLMD snapshots, which is then converted to chemical shift tensors and isotropic shifts for comparison to experiment [J. Phys. Chem. Lett. 2025](https://oulurepo.oulu.fi/bitstream/handle/10024/59168/nbnfioulu-202511136704.pdf?isAllowed=y&sequence=1), [repository record](https://oulurepo.oulu.fi/handle/10024/59168). This is important: it means "equivariant tensor prediction per MD snapshot and then averaging an NMR shift observable" has now been published in at least one non-biomolecular system.

Schmiedmayer et al. provide a similar pattern for EFGs: ML force-field trajectories plus a second rotationally/translationally symmetric EFG model, with finite-temperature quadrupolar coupling predictions [arXiv:2507.19435](https://arxiv.org/abs/2507.19435).

Thus the dynamics part should not be claimed as wholly unprecedented. The honest wording is: scalar shift averaging over MD is standard; tensorial equivariant per-frame NMR prediction is emerging in materials and porous-liquid work; applying it to biomolecular shielding tensors with explicit NMR physical tensor kernels remains unusual.

### Closest published neighbours to the whole pipeline

No single paper I found matches the whole pipeline. The closest neighbours are:

| Neighbour | What matches | What is missing relative to this thesis |
|---|---|---|
| Ben Mahmoud et al. 2024/2025 | Spherical-tensor decomposition; shielding and EFG tensors; equivariant graph models; DFT labels; large-cell and dynamics tests. | No hand-computed ring-current/McConnell/charge physical kernels as input features; silica/materials rather than biomolecular structural biology. |
| Venetos et al. 2023 | Full shielding tensor prediction with an equivariant GNN; DFT labels; tensor target rather than scalar-only shift. | Static silicate setting; no MD averaging; no physical tensor kernels. |
| Zakary and Lantto 2025 | Allegro MLMD plus MatTen NMR-ML; shielding tensor predicted from MLMD snapshots; isotropic shift averaged/compared to experiment. | Porous-liquid xenon rather than protein; no ring-current/McConnell/EFG physical kernel inputs; the reported experimental comparison is isotropic 129Xe shift. |
| Schmiedmayer et al. 2025 | MLMD plus equivariant EFG prediction; finite-temperature NMR quadrupolar observable. | EFG/quadrupolar coupling rather than shielding; materials phase transition; no hand NMR kernels. |
| SEGNN | Physical vector/tensor features as steerable attributes. | General architecture paper, not NMR and not shielding. |
| SPARTA, SPARTA+, SHIFTX2, PPM, UCBShift | Physical/empirical chemical-shift descriptors; scalar shifts from structures; trajectory or ensemble use in some workflows. | Not tensor-resolved, not equivariant tensor-channel learning, mostly scalar empirical prediction. |
| iShiftML | Cheap computed shielding tensors used as ML-informed features for high-level shift prediction. | Not e3nn-style equivariant tensor message passing; not biomolecular MD; not ring-current/McConnell kernels. |

This places the thesis between two existing traditions:

1. Equivariant tensor learning for materials NMR.
2. Biomolecular chemical-shift prediction with physical descriptors and ensemble averaging.

The proposed contribution is the bridge: typed physical NMR tensor kernels inside an equivariant shielding-tensor predictor for biomolecular structures and MD frames. That bridge is not a crowded room.

## What is defensible, and what should not be overclaimed

Defensible claims:

- The decomposition of full rank-2 shielding tensors into `0e + 1e + 2e` is standard linear algebra and matches the e3nn/O(3) representation framework.
- Feeding physical vectors/tensors as steerable attributes is supported by the SEGNN/equivariant-message-passing literature.
- Predicting NMR shielding and EFG tensors with equivariant GNNs is now published, especially in solid-state/materials settings.
- Averaging per-frame NMR predictions over MD trajectories is standard for scalar shifts and emerging for tensorial equivariant NMR workflows in materials/porous liquids.
- The proposed combination is a modest methods novelty if evaluated carefully.

Claims to avoid:

- "No one has predicted NMR tensors with equivariant ML." This is false.
- "No one has used equivariant ML over MD snapshots for NMR." This is now too strong because of Ben Mahmoud et al., Zakary and Lantto, and EFG dynamics work.
- "Physical tensor kernels are a standard input representation for NMR equivariant GNNs." I did not find support for this.
- "MD-averaged prediction equals the experimental chemical shift." The experimental shift is a motional/ensemble average; finite MD is an approximation whose sampling limits must be disclosed.
- "DFT calibration validates against experiment." It validates against the DFT reference. Experimental agreement requires separate referencing, uncertainty, and sampling discussion.

## Method risks and necessary controls

The unusual part is not automatically a problem, but it changes what the thesis must prove.

The central ablation should be:

```text
raw geometry only
vs raw geometry + scalar physical features
vs raw geometry + typed physical tensor kernels
```

If tensor kernels do not improve DFT-tensor prediction or experimental scalar shifts after fair baselines, the negative result is still meaningful: the physical decomposition was symmetry-correct but not useful under the available data/model regime.

Other controls that matter:

- Split train/test by molecule, system, or trajectory block, not by highly correlated MD frames alone.
- Report errors separately for `0e`, `1e`, and `2e` components, plus reconstructed tensor invariants such as isotropic shielding and anisotropy.
- Include rotation tests: rotate structures and confirm predicted tensors rotate as `R T R^T`.
- Include kernel ablations: ring current, McConnell, charge/EFG, and combinations.
- State the reference level: GIPAW/DFT or cluster DFT is the training target, with its own systematic error.
- State the MD sampling limit explicitly: finite MD does not fully sample the experimental ensemble.
- Avoid interpreting learned corrections as unique physical mechanisms unless separate diagnostics support that.

## Bottom line for the thesis

This project is not out on a limb mathematically. The equivariant machinery is established, and full tensor NMR prediction is now a real, citable direction.

But the exact scientific recipe is not standard. The combination of NMR-inspired physical field tensors, typed as e3nn irreps, used as inputs to a DFT-calibrated equivariant shielding tensor model and then averaged over biomolecular MD frames appears genuinely unusual. It is best presented as a careful hybrid method: standard symmetry machinery plus standard NMR physics ideas, combined in a way that seems under-published.

That is a valid MSc position if the thesis is honest about novelty and limits. The examiner should see that the method is not arbitrary, but they should not be led to think this is a populated, routine workflow in structural-biology NMR.

## References

- Anderson, B., Hy, T. S., and Kondor, R. "Cormorant: Covariant Molecular Neural Networks." NeurIPS 2019. [arXiv:1906.04015](https://arxiv.org/abs/1906.04015).
- Batatia, I. et al. "MACE: Higher Order Equivariant Message Passing Neural Networks for Fast and Accurate Force Fields." NeurIPS 2022. [arXiv:2206.07697](https://arxiv.org/abs/2206.07697).
- Batzner, S. et al. "E(3)-equivariant graph neural networks for data-efficient and accurate interatomic potentials." Nature Communications 2022. [DOI:10.1038/s41467-022-29939-5](https://www.nature.com/articles/s41467-022-29939-5).
- Ben Mahmoud, C., Rosset, L. A. M., Yates, J. R., and Deringer, V. L. "Graph-neural-network predictions of solid-state NMR parameters from spherical tensor decomposition." 2024 preprint. [arXiv:2412.15063](https://arxiv.org/abs/2412.15063). Later JCP article: "Graph-neural-network predictions of solid-state NMR parameters in silica from spherical tensor decomposition." J. Chem. Phys. 163, 024118, 2025. [DOI:10.1063/5.0274240](https://www.materials.ox.ac.uk/publication/2245400/dimensions).
- Brandstetter, J., Hesselink, R., van der Pol, E., Bekkers, E. J., and Welling, M. "Geometric and Physical Quantities Improve E(3) Equivariant Message Passing." ICLR 2022. [arXiv:2110.02905](https://arxiv.org/abs/2110.02905).
- Case, D. A. "Chemical shifts in biomolecules." Current Opinion in Structural Biology 2013. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC3877577/).
- de Gortari, I. et al. "Time Averaging of NMR Chemical Shifts in the MLF Peptide in the Solid State." JACS 2010. [DOI:10.1021/ja9062629](https://pubs.acs.org/doi/abs/10.1021/ja9062629).
- Fuchs, F. B. et al. "SE(3)-Transformers: 3D Roto-Translation Equivariant Attention Networks." NeurIPS 2020. [arXiv:2006.10503](https://arxiv.org/abs/2006.10503).
- Geiger, M. and Smidt, T. "e3nn: Euclidean Neural Networks." 2022. [arXiv:2207.09453](https://arxiv.org/abs/2207.09453).
- Han, B., Liu, Y., Ginzinger, S. W., and Wishart, D. S. "SHIFTX2: significantly improved protein chemical shift prediction." J. Biomol. NMR 2011. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC3085061/).
- Kellner, M. et al. "A deep learning model for chemical shieldings in molecular organic solids including anisotropy." 2025. [arXiv:2506.13146](https://arxiv.org/abs/2506.13146); [ShiftML3 predictor](https://shiftml.dev.materialscloud.cscs.ch/).
- Li, D.-W. and Bruschweiler, R. "Certification of Molecular Dynamics Trajectories with NMR Chemical Shifts." J. Phys. Chem. Lett. 2010. [DOI:10.1021/jz9001345](https://pubs.acs.org/doi/10.1021/jz9001345).
- Li, D.-W. and Bruschweiler, R. "PPM: a side-chain and backbone chemical shift predictor for the assessment of protein conformational ensembles." J. Biomol. NMR 2012. [Record](https://ohiostate.elsevierpure.com/en/publications/ppm-a-side-chain-and-backbone-chemical-shift-predictor-for-the-as/).
- Li, J. et al. "Highly Accurate Prediction of NMR Chemical Shifts from Low-Level Quantum Mechanics Calculations Using Machine Learning." 2023. [arXiv:2306.08269](https://arxiv.org/abs/2306.08269).
- Musaelian, A. et al. "Learning local equivariant representations for large-scale atomistic dynamics." Nature Communications 2023. [DOI:10.1038/s41467-023-36329-y](https://www.nature.com/articles/s41467-023-36329-y).
- Ptaszek, A. L. et al. "UCBShift 2.0: Bridging the Gap from Backbone to Side Chain Protein Chemical Shift Prediction for Protein Structures." JACS 2024. [DOI:10.1021/jacs.4c10474](https://pubs.acs.org/doi/10.1021/jacs.4c10474).
- Robustelli, P., Stafford, K. A., and Palmer, A. G. III. "Interpreting Protein Structural Dynamics from NMR Chemical Shifts." JACS 2012. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC3324661/).
- Schmiedmayer, B. et al. "Equivariant machine learning of Electric Field Gradients - Predicting the quadrupolar coupling constant in the MAPbI3 phase transition." 2025. [arXiv:2507.19435](https://arxiv.org/abs/2507.19435).
- Shen, Y. and Bax, A. "Protein backbone chemical shifts predicted from searching a database for torsion angle and sequence homology." J. Biomol. NMR 2007. [DOI record](https://colab.ws/articles/10.1007%2Fs10858-007-9166-6).
- Shen, Y. and Bax, A. "SPARTA+: a modest improvement in empirical NMR chemical shift prediction by means of an artificial neural network." J. Biomol. NMR 2010. [PubMed](https://pubmed.ncbi.nlm.nih.gov/20628786/).
- Thomas, N. et al. "Tensor Field Networks: Rotation- and Translation-Equivariant Neural Networks for 3D Point Clouds." 2018. [arXiv:1802.08219](https://arxiv.org/abs/1802.08219).
- Venetos, M. C., Wen, M., and Persson, K. A. "Machine Learning Full NMR Chemical Shift Tensors of Silicon Oxides with Equivariant Graph Neural Networks." J. Phys. Chem. A 2023. [DOI:10.1021/acs.jpca.2c07530](https://pubs.acs.org/doi/10.1021/acs.jpca.2c07530), [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC10026072/).
- Zakary, O. and Lantto, P. "Equivariant Neural Networks Reveal How Host-Guest Interactions Shape 129Xe NMR in Porous Liquids." J. Phys. Chem. Lett. 2025. [DOI/repository record](https://oulurepo.oulu.fi/handle/10024/59168), [PDF](https://oulurepo.oulu.fi/bitstream/handle/10024/59168/nbnfioulu-202511136704.pdf?isAllowed=y&sequence=1).

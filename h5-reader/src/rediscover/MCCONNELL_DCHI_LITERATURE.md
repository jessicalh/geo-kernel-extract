# McConnell Delta-Chi Literature Values

Status: research note only. Final parameter selection should stay with the
lead.

Local convention note: `../doc/calculators/mcconnell.tex` was present and read;
`src/rediscover/GEOMETRIC_KERNEL_CATALOGUE.md` was not present in this checkout.
The local convention file says the base calculator stores unscaled Angstrom^-3
geometry. Its category scalar sums use
`f = (3 c^2 - 1) / r_A^3` for peptide C=O, peptide C-N, sidechain C=O, and
aromatic bonds. No susceptibility constant is applied by the producer.

## Unit System And Prefactor

All table values below are in cgs molar susceptibility anisotropy units:

`q = Delta chi / (10^-6 cm^3 mol^-1)`.

Many NMR papers report `Delta chi` in `10^-30 cm^3 molecule^-1` or equivalently
`10^-6 Angstrom^3 molecule^-1`. Convert to the requested units by:

`q_molar = 0.602214076 * q_molecule`.

For the local category scalar

`f = (3 cos^2 theta - 1) / r_A^3`

the literature chemical-shift-sign McConnell term is:

`delta_ppm = +0.5535130224 * q_molar * f`.

For the DFT shielding target, use the opposite sign:

`sigma_ppm = -0.5535130224 * q_molar * f`.

Equivalently, if a downstream feature uses shielding-signed geometry
`M_sigma = (1 - 3 cos^2 theta) / r_A^3`, then:

`sigma_ppm = +0.5535130224 * q_molar * M_sigma`.

Derivation:

`chi_SI = 4*pi*10^-12 * q_molar m^3 mol^-1`.

`m_ind = chi_SI * B0 / (mu0 * N_A)`.

`B_dip/B0 = (mu0/4*pi) * m_ind * r_m^-3`, so the `mu0` and `4*pi` factors
cancel the cgs-to-SI susceptibility conversion. With `r_m^-3 = 10^30 r_A^-3`,
`ppm = 10^6`, and the axial traceless factor `1/3`:

`prefactor = 10^24 / (3 N_A) = 0.5535130224`.

No bulk molar volume should be inserted when the input is molar susceptibility;
molar volume is only needed when starting from volume susceptibility. If the
producer stores the raw dipole tensor `(I - 3 u u^T)/r_A^3` and the code forms
the traceless susceptibility tensor as `q*(s s^T - I/3)`, the tensor-product
prefactor before that `1/3` contraction is `10^24/N_A = 1.660539067`.

Sign: the values below are literature `Delta chi` values, not sign-flipped
parameters. The local scalar `f = 3 cos^2 theta - 1` is chemical-shift-signed
for the usual McConnell scalar expression. To emit DFT shielding `sigma`, apply
the negative prefactor above. This is the same shielding-vs-field sign choice
as the ring path's `-n.B`.

## Source Anchors

- McConnell introduced the long-range dipolar shielding model: H. M. McConnell,
  J. Chem. Phys. 27, 226 (1957), DOI
  [10.1063/1.1743676](https://doi.org/10.1063/1.1743676).
- Case's biomolecular-shift review states the tensor form
  `sigma = chi/r^3 - 3 r r^T chi/r^5` and notes that isotropic susceptibility
  cancels, leaving susceptibility anisotropy:
  [Chemical shifts in biomolecules](https://pmc.ncbi.nlm.nih.gov/articles/PMC3877577/).
- Hooper and Kaiser explicitly found multiple carbonyl susceptibility component
  sets from amide spectra, in `10^-6 cm^3 mole^-1`, depending on assumptions:
  [Can. J. Chem. 43, 2363 (1965)](https://ouci.dntb.gov.ua/en/works/4rdPe6N9/),
  DOI [10.1139/v65-318](https://doi.org/10.1139/v65-318).
- Abraham and Ainger report ketone carbonyl anisotropies and a combined
  anisotropy/electric-field model:
  [J. Chem. Soc., Perkin Trans. 2, 1999, 441](https://pubs.rsc.org/en/Content/ArticleLanding/1999/P2/A808908F),
  DOI [10.1039/A808908F](https://doi.org/10.1039/A808908F).
- Abraham, Mobli, and Smith compare carbonyl values from Zurcher, ApSimon,
  Schneider, Williamson, and Abraham, and report conjugated/aromatic carbonyl
  values:
  [Part 19 PDF](https://www.modgraph.co.uk/Downloads/Chemical%20shifts%20of%20aromatic%20aldehydes%20and%20ketones.pdf).
- Abraham, Griffiths, and Perez report amide carbonyl anisotropies for
  aliphatic and aromatic amides, plus a separate nitrogen anisotropy:
  [DOI 10.1002/mrc.3920](https://doi.org/10.1002/mrc.3920), abstract mirror
  [Zendy](https://zendy.io/title/10.1002/mrc.3920).
- Kohlbacher/Lenhof describe a protein shift model using Williamson-Asakura
  C=O and C-N peptide-bond parameters:
  [thesis PDF](https://publikationen.sulb.uni-saarland.de/bitstream/20.500.11880/25764/1/OliverKohlbacher_ProfDrHansPeterLenhof.pdf).
- Aromatic whole-ring susceptibility anisotropy is usually a ring-current
  quantity, not an independent local bond term. A recent secondary web source
  quotes benzene-ring values near `49.5-54.7 * 10^-6 cm^3 mol^-1`:
  [ACS Omega 2024](https://pubs.acs.org/doi/10.1021/acsomega.4c07152).

## Candidate Values

`C=O` and `C-N` entries below are the paired amide/carbonyl-axis constants used
in McConnell/Williamson-style models. When the source reports values per
molecule, the table shows converted molar cgs values.

| Category | Source/model | Candidate q values in 10^-6 cm^3 mol^-1 | Citation support | Confidence |
|---|---:|---:|---|---|
| Peptide C=O, peptide C-N | Williamson-Asakura peptide/protein model | C=O `+2.41`; C-N `-5.42` | Secondary web-cited from Abraham Part 19 comparison (`4, -9` in 10^-30 cm^3 molecule^-1); Kohlbacher confirms WA C=O/C-N use in protein shifts | High for protein lineage; medium for primary traceability |
| Peptide C=O, peptide C-N | Abraham aromatic/conjugated carbonyl | C=O `+3.83`; C-N/perpendicular axis `-7.15` | Web-cited Part 19 PDF (`6.36, -11.88` molecule units); paper states these are close to Williamson peptide values | Medium-high |
| Peptide C=O, peptide C-N | Abraham aliphatic amides | C=O `+6.34`; C-N/perpendicular axis `-14.25` | Web-cited 2013 amide abstract (`10.53, -23.67` molecule units) | Medium for neutral amides; lower for backbone peptides |
| Peptide C=O, peptide C-N | ApSimon aliphatic carbonyl | C=O `+12.65`; C-N/perpendicular axis `-3.61` | Secondary web-cited from Abraham Part 19 comparison (`21, -6` molecule units) | Low-medium for proteins; useful disagreement bound |
| Peptide C=O, peptide C-N | Schneider aliphatic carbonyl | C=O `+14.45`; C-N/perpendicular axis `-7.23` | Secondary web-cited from Abraham Part 19 comparison (`24, -12` molecule units) | Low-medium for proteins; useful disagreement bound |
| Sidechain C=O | Abraham aliphatic amide | C=O `+6.34` | Web-cited 2013 amide abstract; best direct neutral amide sidechain analogue for Asn/Gln | Medium |
| Sidechain C=O | Williamson peptide-like | C=O `+2.41` | Web-cited secondary peptide value; safer if sidechain category is mixed with peptide-like amides | Medium-low |
| Sidechain C=O | Abraham ketone/aldehyde | C=O `+10.30`; perpendicular `+1.93` | Web-cited RSC abstract (`17.1, 3.2` molecule units) | Low for protein sidechains; ketone model not carboxylate/amide |
| Sidechain C=O | Abraham aromatic/conjugated carbonyl | C=O `+3.83` | Web-cited Part 19; relevant only for conjugated/aromatic carbonyls | Low-medium |
| Aromatic ring bond | Drop when RING is active | `0` | Convention-supported: RING already uses Pople/Johnson-Bovey-style ring current intensities | High |
| Aromatic ring bond | Whole benzene ring anisotropy converted to six equivalent in-plane C-C axial bonds | about `+16.5` to `+18.2` per aromatic C-C bond, inferred from whole-ring `Delta chi_ax = -49.5` to `-54.7` | Secondary web-cited; inference assumes six identical in-plane axial bond tensors, where `Delta chi_ring = -3 Delta chi_bond` | Low; standalone-only |

Topology note: if `sidechain C=O` mixes neutral amide carbonyls (Asn/Gln) with
carboxylate/carboxylic acid carbonyls (Asp/Glu), the category is chemically
too broad. This report did not find a tight-pass web-primary carboxylate value.
Do not solve that by IUPAC name parsing; split by topology/atom typing if the
lead wants separate values.

## Recommended Literature Set To Try

Use this as the first literature-scaled report set, not as a final lock-in:

| Typed category | Recommended q in 10^-6 cm^3 mol^-1 | Why |
|---|---:|---|
| Peptide C=O | `+2.41` | Williamson-Asakura is the protein-predictor lineage and is close to Abraham's conjugated/aromatic carbonyl re-fit. |
| Peptide C-N | `-5.42` | Same WA peptide lineage; keeps the amide-plane anisotropy paired with the C=O value. |
| Sidechain C=O | `+6.34` | Best direct neutral aliphatic amide analogue for Asn/Gln sidechain carbonyls from Abraham 2013. Mark low-confidence if Asp/Glu carboxylates are included. |
| Aromatic | `0` with RING active | Avoid double-counting the same aromatic pi-current physics already emitted by RING. |

Alternative conservative set if the lead wants one source family only:

| Typed category | q in 10^-6 cm^3 mol^-1 |
|---|---:|
| Peptide C=O | `+2.41` |
| Peptide C-N | `-5.42` |
| Sidechain C=O | `+2.41` |
| Aromatic | `0` |

This single-family set is less chemically specific for sidechains, but it avoids
mixing Abraham sidechain-amide values with Williamson peptide values.

## Aromatic Double-Count Flag

Do not enable an aromatic McConnell `Delta chi` bond term at the same time as
the explicit RING/Johnson-Bovey/Biot-Savart path unless the term is proven to
represent only a residual sigma-bond anisotropy with the pi ring current
removed.

Reason: the large aromatic susceptibility anisotropy is the macroscopic version
of the induced pi ring current. The RING path already scales explicit
Biot-Savart kernels by literature currents (`PHE -12.0`, `TYR -11.28`,
`TRP-6 -12.48`, `TRP-5 -6.72`, `HIS -5.16 nA/T`). Adding a benzene-derived
McConnell aromatic bond value would count the same physical current a second
time.

For a McConnell-only ablation with RING disabled, a possible low-confidence
standalone aromatic value is the inferred per-C-C-bond `+18.2` above. It should
not be used in the production RING+McConnell feature set without a reconciliation
test against the existing ring current emit.

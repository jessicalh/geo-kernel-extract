# Broad-backbone T2 law-distillation evidence

Run date: 2026-06-01. Branch: `h5-reader-pysr-spike`.

This is evidence only. It does not declare that a law fell out.

## Inputs and artifacts

- Substrate: `/tmp/rdc-broad-backbone-axes`
- Evidence dir: `/tmp/rdc-backbone-law-evidence`
- Main script: `backbone_distill_evidence.py`
- PySR script: `backbone_pysr_distill.py`
- Model checkpoints: `/tmp/rdc-backbone-law-evidence/models/backbone_<stratum>.pt`
- Tables: `frame_split_metrics.csv`, `radial_fit_summary.csv`,
  `radial_gate_scan.csv`, `atom_split_validation.csv`,
  `pysr_summary_gate_chosen.csv`, `pysr_summary_gate_norm.csv`
- Plots: `/tmp/rdc-backbone-law-evidence/figures/radial_{ring,bond,charge}.png`

## Backbone model reproduction

The e3nn axes model was retrained for 4000 epochs per stratum with the frozen
change-of-basis path. Frame-split T2 R2 / |T2| r:

| stratum | atoms | T2 R2 | \|T2\| r |
|---|---:|---:|---:|
| N | 54 | 0.647 | 0.764 |
| CA | 54 | 0.427 | 0.594 |
| C | 54 | 0.585 | 0.804 |
| O | 54 | 0.713 | 0.836 |
| HN | 52 | 0.757 | 0.864 |
| HA | 50 | 0.674 | 0.785 |
| HA2 | 4 | 0.838 | 0.934 |
| HA3 | 4 | 0.915 | 0.966 |

HA2/HA3 are thin strata (4 atoms).

## Learned radial readout

Radial fits use sampled learned gates from each trained model. Ring/bond use the
axis gate as the analytic-compatible signed gate when axes are present; charge
uses the displacement gate. `radial_gate_scan.csv` also scans all gates and the
gate norm.

### Ring

- Chosen signed gate: `log(|gate/intensity|) ~ log(r)` median power -2.45,
  median log-space R2 0.14 across strata.
- Chosen signed gate vs `intensity * r^-3`: median R2 0.061, range 0.0002-0.466.
- All-gate scan: median power -2.56; median log-space R2 0.18; best
  `intensity * r^-3` R2 0.519 on one gate/stratum.

Evidence status: weak/partial radial match, not a robust Pople-like readout in
this broad-backbone model.

### Bond

- Chosen signed gate: `log(|gate|) ~ log(r) + category offsets` median power
  -2.57, median log-space R2 0.53.
- Chosen signed gate vs category-specific `r^-3`: median R2 0.557, range
  0.157-0.992.
- All-gate scan: median `r^-3` comparator R2 0.652; strongest strata/gates:
  C axis 0.992, O axis 0.987, HN/N axis about 0.96/0.89.

Evidence status: the strongest radial evidence is bond-like; it is uneven across
target strata.

### Charge

- Chosen signed gate: `log(|gate/q|) ~ log(r)` median power -2.70, median
  log-space R2 0.17.
- Chosen signed gate vs `q * r^-2`: median R2 0.269, range 0.008-0.920.
- Strongest match is N (`q * r^-2` R2 0.920); CA/HA/O are moderate; C/HN/HA2
  are weak.

Evidence status: mixed. One stratum matches a field-like radial strongly, but
the pooled/readout evidence is not uniformly Buckingham-like.

## PySR closed forms

PySR was run in `analysis/venv` on sampled learned per-source readouts.

Signed selected gate (`pysr_summary_gate_chosen.csv`):

| source | expression | R2 |
|---|---|---:|
| ring | `(ring_intensity * 0.073672436) - -0.9611622` | 0.054 |
| bond pooled | `(cos_theta / ((r - cos_theta) - 2.4455519)) / (square(bond_category / cos_theta) - (cos_theta - 2.4899824))` | 0.121 |
| charge | `1.3059723 / (square(r) - 4.9755645)` | 0.334 |

Gate norm with per-bond-category fits (`pysr_summary_gate_norm.csv`):

| source | expression | R2 |
|---|---|---:|
| ring | `9.850933 - r` | 0.095 |
| bond pooled | `7.089535 - r` | 0.221 |
| bond cat0 | `7.1976647 - r` | 0.251 |
| bond cat1 | `square(square(square(1.2887528 - square((r * 0.32364473) + -1.0542474))))` | 0.434 |
| bond cat3 | `square(5.023032 / r)` | 0.269 |
| bond cat4 | `square((cos_theta / r) * 8.781364)` | 0.251 |
| charge | `-1.0386896 / (2.1492085 - r)` | 0.202 |

Evidence status: PySR did not return a compact high-R2 textbook-form expression
for the pooled broad-backbone readouts. The best symbolic evidence is weak to
moderate, and some higher-R2 expressions are not physically clean.

## Atom-split validation

Validation used source-kind-only e3nn models trained on one 50/50 atom split per
stratum for 300 epochs. This tests learned source-kind signal transfer to atoms
held out from the fit; it is not an exhaustive leave-one-atom-out sweep.

| source | held-out atoms median | T2 R2 median | R2 range | \|T2\| r median |
|---|---:|---:|---:|---:|
| ring | 26.5 | -0.041 | -4.982 to 0.043 | 0.146 |
| bond | 26.5 | 0.478 | 0.359 to 0.600 | 0.682 |
| charge | 26.5 | 0.380 | 0.148 to 0.469 | 0.622 |

HA2/HA3 validation splits hold out only 2 atoms each and should be treated as
thin. Ring source-kind-only validation is negative in this setup; bond and
charge transfer positively on the single split.

## Literature-coefficient-fixed gate

Not run as a T2 fixed-coefficient predictor. The current emitted substrate does
not contain per-source or aggregated fixed-literature T2 tensors for ring, bond,
or charge. Constructing those in this Python analysis would require building the
textbook angular T2/source field in Python rather than reading an emitted column.
That is flagged as a substrate gap under `analysis/PATTERNS.md`.

The radial-only comparator fits above use requested `r^-3` / `r^-2` readouts
against learned gates; they are not a no-fit fixed-coefficient DFT predictor.

## Discipline check

- No C++ was edited by this pass.
- No merge was performed and nothing was committed.
- `cob.get_C()` was called and the run printed max `|C.T C - I| = 1.110e-16`.
- e3nn is used through the existing `equiv_t2_backbone_e3nn.py` path.
- PySR was run through `analysis/venv`.
- Grep for forbidden recompute patterns in the new scripts found no
  `(3cos...)`, no `q*disp`, no H5 access, and no tensor-target construction.
  The new script contains only the explicitly requested radial-only comparator
  terms `np.power(r, -3.0)` / `np.power(r, -2.0)`.


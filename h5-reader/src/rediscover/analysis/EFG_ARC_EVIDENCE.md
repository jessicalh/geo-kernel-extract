# APBS-EFG T2 Arc Evidence

Run date: 2026-06-01. Branch: `h5-reader-pysr-spike`.

This is corrected capstone evidence after removing the EFG lab-frame rotation
confound. It is evidence only; it does not claim that a closed-form law has
fallen out.

## Corrected Inputs And Discipline

- Corrected substrate: `/tmp/rdc-efg-localframe-capstone`
- Local-frame audit: `/tmp/rdc-efg-localframe-audit-capstone`
- Evidence dir: `/tmp/rdc-efg-arc-evidence-localframe-capstone`
- Main scripts: `efg_localframe_audit.py`, `equiv_t2_efg_e3nn.py`,
  `efg_distill_evidence.py`
- Tables: `efg_localframe_audit.csv`, `efg_distill_summary.csv`,
  `efg_gate_samples.csv`, `efg_distill_run.json`
- Read only emitted sidecars: `efg_aggregated.csv`, `efg_feature_T2.npy`,
  `efg_target_T2.npy`, `efg_feature_lab_T2.npy`, `efg_target_lab_T2.npy`
- Reused frozen `change_of_basis.get_C()`.

The corrected C++ emitter writes local-frame EFG feature/target sidecars for the
main model and lab-frame sidecars only for audit. Python does not read H5, does
not recompute APBS/DFT EFG tensors, and does not define a second physical model.

The corrected extraction emitted 423000 total EFG rows, 163000 finite
local-frame backbone rows, and 163000 valid local frames. Frame validity is
100% for the six backbone strata used below.

## Rotation-Confound Audit

The high frame-to-frame autocorrelation in the DFT EFG target was a lab-frame
orientation signature. After rotating feature and target tensors into the
backbone local frame before decomposition, the target lag-1 correlation drops
from roughly 0.75-0.86 to roughly 0.05-0.23.

| stratum | lab target lag1 | local target lag1 | lab gamma R2 | local gamma R2 |
|---|---:|---:|---:|---:|
| C | 0.863 | 0.055 | 0.115 | 0.0146 |
| CA | 0.752 | 0.090 | 0.0007 | 0.0031 |
| HA | 0.789 | 0.233 | 0.0325 | 0.0001 |
| HN | 0.772 | 0.158 | 0.0725 | 0.0144 |
| N | 0.861 | 0.077 | 0.0108 | 0.0234 |
| O | 0.803 | 0.054 | 0.315 | 0.0004 |

This is the gating check for the corrected rerun: the old O and C EFG signals
were lab-frame rotation artifacts. The local-frame substrate removes that
confound before the model/evidence pass.

## Corrected Depth-A Readout

Model form:

```text
DFT_T2 ~= g(|EFG|) * EFG_T2
```

The corrected local-frame substrate leaves only weak/null held-out predictive
signal. The fitted constant-g linear form is reported as a diagnostic, not as an
independently fixed literature coefficient.

Blocked frame split, purge=1. `cross_split_lag1_pairs` is zero in the corrected
blocked run.

| stratum | atoms | N_eff_lag1 | constant R2 | nonlinear R2 | nonlinear gain | constant `|T2| r` | nonlinear `|T2| r` |
|---|---:|---:|---:|---:|---:|---:|---:|
| N | 54 | 25219 | 0.0209 | 0.0229 | 0.0020 | 0.082 | 0.059 |
| CA | 54 | 23495 | 0.0031 | 0.0045 | 0.0014 | 0.033 | 0.042 |
| C | 54 | 26025 | 0.0146 | 0.0165 | 0.0019 | 0.042 | 0.064 |
| O | 54 | 25708 | -0.0006 | -0.0004 | 0.0001 | 0.008 | 0.004 |
| HN | 52 | 21504 | 0.0122 | 0.0243 | 0.0121 | 0.202 | 0.180 |
| HA | 58 | 21121 | -0.0043 | -0.0047 | -0.0004 | 0.144 | 0.080 |

Random-split diagnostics are not used as the headline. In the variance
decomposition rerun, random EFG produced 170 cross-split lag-1 pairs, while the
blocked corrected run produced zero.

## Local-Frame Vs Lab-Frame Fit Comparison

The same fitter on the old lab-frame substrate reproduces why the fix matters:
the apparent lab-frame C/O performance disappears in the corrected local-frame
rerun.

| stratum | lab R2 | lab `|T2| r` | local R2 | local `|T2| r` |
|---|---:|---:|---:|---:|
| N | -0.012 | 0.130 | 0.023 | 0.059 |
| CA | -0.013 | -0.050 | 0.005 | 0.042 |
| C | 0.112 | -0.105 | 0.016 | 0.064 |
| O | 0.342 | 0.136 | -0.000 | 0.002 |
| HN | 0.063 | 0.072 | 0.024 | 0.180 |
| HA | 0.003 | 0.050 | -0.005 | 0.076 |

## De-Circularisation

No defensible fixed literature coefficient was identified for a universal
`sigma_T2 = gamma * APBS_EFG_T2` predictor across these backbone strata.
Buckingham-style shielding response coefficients are nucleus- and
environment-dependent, and the explicit EFG-shielding literature does not
provide a transferable scalar coefficient for this emitted APBS rank-2 tensor
test.

Therefore the de-circularised result is form/evidence only. The corrected local
frame shows that APBS-EFG is not the current explanatory axis for the DFT T2
target once the orientation confound is removed.

## Caveats

- T2 headline metrics use `|T2| r`; scalar gamma R2 is only a diagnostic for a
  fitted one-coefficient projection.
- T2 variance shares in the companion variance table are pooled-Frobenius shares
  over emitted rank-2 components.
- The DFT substrate is r2SCAN/def2-SVP with CPCM(Water). Fixed-charge/APBS EFG
  comparisons carry a solvation-treatment and basis-set mismatch.
- F1 is resolved elsewhere: McConnell/ring current use the canonical `K * chi`
  traceless-chi PCS convention and were not edited.

## Verification

Commands represented by this doc:

```text
build/linux-gcc/h5reader_extract --case efg --out /tmp/rdc-efg-localframe-capstone
python3 src/rediscover/analysis/efg_localframe_audit.py /tmp/rdc-efg-localframe-capstone --out-dir /tmp/rdc-efg-localframe-audit-capstone
python3 src/rediscover/analysis/equiv_t2_efg_e3nn.py /tmp/rdc-efg-localframe-capstone --epochs 4000 --hidden 32 --lr 3e-3 --split blocked --purge-frames 1
python3 src/rediscover/analysis/efg_distill_evidence.py /tmp/rdc-efg-localframe-capstone --evidence-dir /tmp/rdc-efg-arc-evidence-localframe-capstone --epochs 4000 --hidden 32 --lr 3e-3 --split blocked --purge-frames 1
```

No merge was performed.

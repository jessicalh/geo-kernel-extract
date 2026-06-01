# APBS-EFG T2 Arc Evidence

Run date: 2026-06-01. Branch: `h5-reader-pysr-spike`.

This is evidence only. The "law fell out" verdict is reserved for the lead +
Jessica.

## Inputs and discipline

- Substrate: `/tmp/rdc-efg`
- Evidence dir: `/tmp/rdc-efg-arc-evidence`
- Main scripts: `equiv_t2_efg_e3nn.py`, `efg_distill_evidence.py`
- Tables: `efg_distill_summary.csv`, `efg_gate_samples.csv`
- Read only emitted sidecars: `efg_aggregated.csv`, `efg_feature_T2.npy`,
  `efg_target_T2.npy`
- Reused frozen `change_of_basis.get_C()`; printed
  `max |C.T C - I| = 1.110e-16`
- No H5 read, no EFG/projection recompute, no C++ edit in this pass.

## Depth A readout

Model form:

```text
DFT_T2 ~= g(|EFG|) * EFG_T2
```

The trained MLP gate is not numerically flat when sampled over the emitted
`|EFG|` range. However, the |EFG|-dependent part adds little held-out predictive
power: the fitted constant-g predictor accounts for essentially all of the
nonlinear scalar-gate R2 in the signal-bearing strata. I therefore treat the
supported form as the linear Buckingham-like form, with a fitted coefficient:

```text
sigma_T2 ~= gamma_stratum * EFG_T2
```

The coefficient is not independently fixed by literature in this run.

Frame split is the same deterministic frame split as the capture fitter.
`N_eff_lag1` is a lag-1 frame-equivalent count from per-atom target-|T2|
autocorrelation; `atoms` is the cross-environment unit count.

| stratum | gate form readout | atoms | N_eff_lag1 | fitted constant g | constant R2 | nonlinear R2 | nonlinear gain | absT2 r constant | absT2 r nonlinear |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| O | learned gate drifts, constant form predictive | 54 | 5179.97 | +152.415 | +0.316 | +0.327 | +0.0109 | +0.204 | +0.311 |
| C | learned gate drifts, constant form predictive | 54 | 3679.80 | -20.038 | +0.115 | +0.120 | +0.0057 | +0.044 | +0.134 |
| HN | learned gate drifts, constant form predictive | 52 | 9033.43 | +1.9467 | +0.073 | +0.077 | +0.0036 | +0.199 | +0.198 |
| HA | weak signal; constant form predictive | 58 | 10488.27 | +2.0402 | +0.030 | +0.031 | +0.0010 | +0.003 | +0.005 |
| N | weak/null signal; gate drift not interpretable | 54 | 5645.30 | -9.6606 | +0.010 | +0.014 | +0.0040 | +0.205 | +0.133 |
| CA | null signal; gate drift not interpretable | 54 | 10342.75 | -0.7064 | +0.000 | +0.002 | +0.0021 | +0.006 | -0.051 |

Gate sample diagnostics over the central emitted `|EFG|` range:

| stratum | EFG p05-p95 | gate p05-p95 | gate-vs-EFG r | constant share of nonlinear R2 |
|---|---:|---:|---:|---:|
| O | 1.120-2.246 | +122.4 to +197.6 | -0.99 | 96.7% |
| C | 0.974-2.351 | -27.37 to -17.00 | +0.75 | 95.3% |
| HN | 0.784-1.746 | +1.444 to +2.604 | +0.82 | 95.3% |
| HA | 0.324-0.834 | +1.633 to +2.906 | +0.31 | 96.7% |
| N | 0.711-1.819 | -27.27 to -4.79 | +0.58 | 72.0% of tiny R2 |
| CA | 0.827-1.582 | -2.17 to -0.067 | +0.61 | 9.4% of null R2 |

Interpretation by stratum:

- O: strongest signal. Nonlinear gate improves R2 by only +0.011 over the
  fitted constant-g law, so the defensible form is linear despite a visibly
  drifting learned gate.
- C: same pattern at lower signal; constant-g captures 95% of nonlinear R2.
- HN: signal is weaker but the constant and nonlinear models are effectively
  tied.
- HA: weak but positive; nonlinear dependence is not useful.
- N and CA: weak/null. The learned gate drift is not a reliable law readout
  because the target signal is near zero.

PySR was not run: the readout is one-dimensional and the constant-g comparator
already answers the Depth-A question more directly than symbolic regression.

## De-circularisation

I did not find a defensible fixed literature coefficient that can be plugged in
as an unfitted `sigma_T2 = gamma * EFG_T2` predictor for these backbone strata
and this APBS-EFG T2 substrate.

What the literature supports:

- Buckingham's 1960 paper gives the classical electrostatic-shielding framework
  and a proton X-H electric-field coefficient, but not a universal rank-2
  APBS-EFG-to-shielding-tensor coefficient for backbone O/C/N/H/CA/HA
  environments.
- Protein CSP work uses Buckingham's equation with electric-field and
  field-gradient response terms. Those terms are nucleus-dependent response
  properties, and in practice values are drawn from separate prior work or
  fitted/approximated for specific reporter nuclei such as backbone HN and N.
- The explicit EFG-shielding theory literature describes shielding
  polarizability tensors for uniform electric-field gradients, but also flags
  origin/gauge dependence of the tensor components. That is not a clean scalar
  `gamma_stratum` for the emitted APBS EFG T2 vectors.
- Finite-field shielding-polarizability papers compute molecule/nucleus-specific
  tensors for small molecules. They do not supply a transferable protein
  backbone coefficient for this APBS-EFG tensor test.

Sources checked: the local references corpus
`/shared/2026Thesis/nmr-shielding/references/` and primary literature/web
records for Buckingham 1960 (`doi:10.1139/v60-040`), Batchelor 1975
(`doi:10.1021/ja00845a022`), McDowell & Buckingham 1992
(`doi:10.1039/FT9928803281`), Lazzeretti 2003
(`doi:10.1016/S0166-1280(03)00264-1`), Kukic et al. 2013
(`doi:10.1021/ja406995j`), Facelli 2011, and Sahakyan & Vendruscolo 2013
(`doi:10.1021/jp3057306`).

Therefore the de-circularised result is not a fixed-coefficient DFT prediction.
It is form recovery only: the fitted constant-g linear form carries the EFG
signal, but the coefficient remains fitted, not fixed.

## Verification

Commands run:

```text
python3 -m py_compile src/rediscover/analysis/equiv_t2_efg_e3nn.py src/rediscover/analysis/efg_distill_evidence.py
python3 src/rediscover/analysis/efg_distill_evidence.py /tmp/rdc-efg --evidence-dir /tmp/rdc-efg-arc-evidence --epochs 4000 --hidden 32 --lr 3e-3
python3 src/rediscover/analysis/equiv_t2_efg_e3nn.py /tmp/rdc-efg --epochs 4000 --hidden 32 --lr 3e-3
```

The original fitter reproduced the capture metrics:

```text
N  R2=+0.014  |T2| r=+0.133
CA R2=+0.002  |T2| r=-0.051
C  R2=+0.120  |T2| r=+0.134
O  R2=+0.327  |T2| r=+0.311
HN R2=+0.077  |T2| r=+0.198
HA R2=+0.031  |T2| r=+0.005
```

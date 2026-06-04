# 1P9J interpolation build result - 2026-06-04

Built the source-e3nn v0 advisor graph from a fresh bounded ring-current producer
re-emission. The extractor was not modified. The composed-vs-procedural oracle
passed byte/numeric parity for 141,000 source rows, 20,500 aggregate rows, and all
ring-current CSV/NPY sidecars. Through-space training consumed 110,500 emitted
source rows, not the 558,360 all-atom aggregate rows.

Clean held-out protocol: blocked/purged frame split, train-only atom centering,
train-only feature normalization, `cross_split_lag1_pairs = 0`, frozen
`change_of_basis.get_C()` orthogonality `1.11e-16`, 5-component T2 target kept
intact. Best model was `learnable` cross. Honest held-out number: T2 component
R2 `+0.4776`; `|T2|` r `+0.8188`; T0 modulation R2 `+0.4537`.

Graph:
`/tmp/rediscover-runs/2026-06-04-interp-1p9j-e3nn/interp_1p9j_advisor_graph.png`
and `.pdf` beside it. Panel A is held-out centered T2 component scatter. Panel B
is held-out `|T2|` modulation recovery. Panel C uses T0 modulation as the main
axis and demotes restored `sigma_iso` to a baseline-dominated inset. Panel D
states the within-axis / geometry-sampler / correlate-not-match /
direction-not-destination caveat.

Artifacts:
`interp_1p9j_predictions.csv`, `interp_1p9j_metrics.json`, and
`interp_1p9j_run_audit.json` are in the same run directory. The clean number is
the gate baseline; the old `0.466/0.757` is pre-protocol-fix leaky-suspect and
was not used. Ungrounded/not claimed: no transferability, no between-axis
recovery, no BMRB validation, no dynamics/process model; the graph is the bounded
ring-current source-e3nn branch, not the full broad/all-atom Stage-3 chewer.

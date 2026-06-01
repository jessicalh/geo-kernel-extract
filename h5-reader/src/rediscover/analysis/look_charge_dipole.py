#!/usr/bin/env python3
"""look_charge_dipole — within-atom signal check for FF14SB charge dipoles.

Consumes the C++-EMITTED charge_dipole_aggregated.csv only: the producer ships
the charge dipole mu = (mu_x, mu_y, mu_z, mu_norm) as a reducer output. Python
fits on those emitted columns; it does NOT recompute any field.

REMOVED (feedback_no_python_physics_except_labeled_integrity_test, 2026-06-01 +
lead decision): the `E = Sigma q*disp/r^3` field recompute that previously built
field_{x,y,z} from the per-source source_q_e / disp_local / r columns. The charge
ELECTRIC FIELD / EFG is not a Python quantity to rebuild: it is the C++
`buckingham_efield` relationship (APBS-derived), a fail-loud stub on this branch
until its data flows. When that relationship emits, read its field columns here;
never recompute q*d/r^3. The emitted charge MULTIPOLE on this branch is mu, so
this script fits mu only.
"""
import os
import sys

import numpy as np
import pandas as pd


out_dir = sys.argv[1] if len(sys.argv) > 1 else os.environ.get(
    "REDISCOVER_OUT", "/tmp/rediscover-chargedipole"
)
agg_path = f"{out_dir}/charge_dipole_aggregated.csv"


def finite_mask(*cols):
    m = np.ones(len(cols[0]), dtype=bool)
    for c in cols:
        m &= np.isfinite(np.asarray(c, dtype=float))
    return m


def corr(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = finite_mask(x, y)
    if ok.sum() < 3:
        return np.nan
    x = x[ok]
    y = y[ok]
    if np.std(x) <= 0 or np.std(y) <= 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def r2_score(y, yhat):
    y = np.asarray(y, dtype=float)
    yhat = np.asarray(yhat, dtype=float)
    ok = finite_mask(y, yhat)
    if ok.sum() < 3:
        return np.nan
    den = float(np.sum(y[ok] ** 2))
    if den <= 0:
        return np.nan
    return float(1.0 - np.sum((y[ok] - yhat[ok]) ** 2) / den)


def fit_all(df, cols):
    X = df[[f"{c}_w" for c in cols]].to_numpy(float)
    y = df["sigma_w"].to_numpy(float)
    ok = finite_mask(y, *[X[:, i] for i in range(X.shape[1])])
    if ok.sum() <= X.shape[1]:
        return np.full(X.shape[1], np.nan), np.nan
    coef, *_ = np.linalg.lstsq(X[ok], y[ok], rcond=None)
    return coef, r2_score(y[ok], X[ok] @ coef)


def leave_atom_out(df, cols):
    pred = np.full(len(df), np.nan)
    Xall = df[[f"{c}_w" for c in cols]].to_numpy(float)
    yall = df["sigma_w"].to_numpy(float)
    atoms = np.sort(df.atom_index.unique())
    usable = 0
    for atom in atoms:
        te = df.atom_index.to_numpy() == atom
        tr = ~te
        X = Xall[tr]
        y = yall[tr]
        ok = finite_mask(y, *[X[:, i] for i in range(X.shape[1])])
        if ok.sum() <= X.shape[1]:
            continue
        coef, *_ = np.linalg.lstsq(X[ok], y[ok], rcond=None)
        pred[te] = Xall[te] @ coef
        usable += 1
    return r2_score(yall, pred), usable


def lag1(x):
    x = np.asarray(x, dtype=float)
    ok = np.isfinite(x)
    x = x[ok]
    if len(x) < 4 or np.std(x) <= 0:
        return np.nan
    return corr(x[:-1], x[1:])


def ar1_effective_rows(df, feature):
    total = 0.0
    used = 0
    for _, g in df.sort_values("h5_row").groupby("atom_index"):
        x = g[f"{feature}_w"].to_numpy(float)
        y = g["sigma_w"].to_numpy(float)
        n = len(g)
        rx = lag1(x)
        ry = lag1(y)
        if not np.isfinite(rx) or not np.isfinite(ry):
            continue
        prod = float(np.clip(rx * ry, -0.95, 0.95))
        ne = n * (1.0 - prod) / (1.0 + prod)
        total += float(np.clip(ne, 1.0, n))
        used += 1
    return total, used


df = pd.read_csv(agg_path)
df = df[df.dft_present == 1].copy()

# Fit only on the C++-EMITTED charge dipole mu. (No field recompute — see header.)
feature_sets = {
    "mu_norm": ["mu_norm"],
    "mu_xyz": ["mu_x", "mu_y", "mu_z"],
}
scalar_features = ["mu_norm", "mu_x", "mu_y", "mu_z"]

df["sigma_w"] = df.dft_sigma_iso - df.groupby("atom_index").dft_sigma_iso.transform("mean")
for c in scalar_features:
    df[f"{c}_w"] = df[c] - df.groupby("atom_index")[c].transform("mean")

atoms = np.sort(df.atom_index.unique())
frames_per_atom = df.groupby("atom_index").size()
print(
    f"rows={len(df)}  HN atoms={len(atoms)}  "
    f"median_frames_per_atom={frames_per_atom.median():.0f}"
)
print(
    "effective_N: leave-atoms-out unit is atom; "
    f"raw rows={len(df)} should not be read as independent samples"
)
for feat in ["mu_norm"]:
    ne, used = ar1_effective_rows(df, feat)
    print(f"  AR1-adjusted row N for pooled corr {feat}: {ne:.1f} across {used} atoms")

print("\n== pooled within-atom scalar correlations ==")
pooled = {}
for c in scalar_features:
    pooled[c] = corr(df[f"{c}_w"], df.sigma_w)
    print(f"  {c:>10s}: r={pooled[c]:+.4f}")

print("\n== shared-coefficient fits on within-atom modulation ==")
for name, cols in feature_sets.items():
    coef, in_r2 = fit_all(df, cols)
    loo_r2, usable_atoms = leave_atom_out(df, cols)
    coef_s = ", ".join(f"{v:+.4g}" for v in coef)
    print(
        f"  {name:>10s}: in-sample R2={in_r2:+.4f}  "
        f"leave-one-atom-out R2={loo_r2:+.4f}  atoms_predicted={usable_atoms}  coef=[{coef_s}]"
    )

print("\n== per-atom best scalar within-atom correlation ==")
per_rows = []
for atom, g in df.groupby("atom_index"):
    vals = {c: corr(g[f"{c}_w"], g.sigma_w) for c in scalar_features}
    finite = {k: v for k, v in vals.items() if np.isfinite(v)}
    if finite:
        best = max(finite.items(), key=lambda kv: abs(kv[1]))
        best_name, best_r = best
    else:
        best_name, best_r = "none", np.nan
    per_rows.append(
        {
            "atom_index": int(atom),
            "atom_name": str(g.atom_name.iloc[0]),
            "n": int(len(g)),
            "sigma_sd": float(g.sigma_w.std(ddof=0)),
            "best_feature": best_name,
            "best_r": best_r,
            "r_mu_norm": vals["mu_norm"],
        }
    )

per = pd.DataFrame(per_rows).sort_values("atom_index")
for row in per.itertuples(index=False):
    print(
        f"  atom={row.atom_index:4d} {row.atom_name:>5s} n={row.n:3d} "
        f"sigma_sd={row.sigma_sd:8.4f} best={row.best_feature:>10s} r={row.best_r:+.4f} "
        f"r(|mu|)={row.r_mu_norm:+.4f}"
    )

best_set = max(feature_sets, key=lambda k: leave_atom_out(df, feature_sets[k])[0])
best_loo, _ = leave_atom_out(df, feature_sets[best_set])
print("\n== read ==")
if np.isfinite(best_loo) and best_loo > 0.02:
    print(f"  strongest universal readout is {best_set} with leave-one-atom-out R2={best_loo:+.4f}")
else:
    print(f"  no convincing universal predictor in this first look; best leave-one-atom-out R2={best_loo:+.4f} ({best_set})")

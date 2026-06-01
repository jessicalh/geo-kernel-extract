#!/usr/bin/env python3
"""pysr_distill — find the CLOSED-FORM symbolic ring-current kernel.

Ring current is a sum over a variable source set, and no (atom,frame) group has
a single source (~7 rings each), so a pooling model is needed to assign the
per-source contribution. sumpool_kernel.py did that: it learned g(r,cos,intensity)
whose sum reconstructs the producer kernel (test R2=0.97), and the read-out
matches the Pople form (angular R2=0.99, axial 1/r^3 R2=0.98).

PySR now distils that learned per-source function g into a symbolic equation.
If it returns ~ intensity*(3cos^2-1)/r^3, the ring-current law has fallen out in
closed form. Far-field only (r>=4 A), where the point-dipole form is valid.
"""
import numpy as np
from pysr import PySRRegressor

d = np.load("/tmp/rediscover-out/sumpool_kernel_readout.npz")
gsrc, r, ct, inten = d["gsrc"], d["r"], d["ct"], d["inten"]
dipolar = d["dipolar"]            # C++-emitted (3cos^2-1)/r^3 per source
far = (r >= 4.0) & np.isfinite(gsrc)
X = np.column_stack([r[far], ct[far], inten[far]])
y = gsrc[far]

# subsample for tractable evolutionary search
rng = np.random.default_rng(0)
idx = rng.choice(len(y), size=min(6000, len(y)), replace=False)
X, y = X[idx], y[idx]
print(f"PySR distillation: {len(y)} far-field samples; vars = r, ct(=cos theta), inten")

model = PySRRegressor(
    niterations=60,
    populations=24,
    population_size=40,
    maxsize=22,
    binary_operators=["+", "-", "*", "/"],
    unary_operators=["square"],
    variable_names=["r", "ct", "inten"],
    elementwise_loss="L2DistLoss()",
    constraints={"/": (-1, 9)},
    timeout_in_seconds=600,
    progress=False,
    deterministic=True,
    parallelism="serial",
    random_state=0,
    verbosity=0,
    temp_equation_file=True,
)
model.fit(X, y)

print("\n== Pareto front (complexity vs accuracy) ==")
eq = model.equations_
for _, row in eq.iterrows():
    print(f"  cx={int(row['complexity']):2d}  loss={row['loss']:.3e}  {row['equation']}")

best = model.get_best()
print(f"\nbest: {best['equation']}")
pred = model.predict(X)
r2 = 1 - ((y-pred)**2).sum()/((y-y.mean())**2).sum()
print(f"best R2 on far-field samples: {r2:+.4f}")

# explicit comparison: how well does the producer's EMITTED kernel fit the
# learned g? (read the emitted `dipolar` column, do NOT recompute (3ct^2-1)/r^3)
pople = (inten * dipolar)[far]          # intensity * C++-emitted dipolar
s, b = np.polyfit(pople, y, 1)
r2p = 1 - ((y-(s*pople+b))**2).sum()/((y-y.mean())**2).sum()
print(f"emitted intensity*dipolar: R2={r2p:+.4f}  coeff={s:+.4g}  "
      f"(comparator is the C++-emitted kernel column, not a numpy recompute)")

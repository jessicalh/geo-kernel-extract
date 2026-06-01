#!/usr/bin/env python3
"""refine_kernel — far-field coefficient + a figure of the recovered kernel.

Loads the sum-pooling read-out (the learned per-source g on every real source,
from sumpool_kernel.py) and, in the regime where the point-dipole form is valid
(r >= 4 A), fits the single coefficient A in

    g  ~  A * intensity * dipolar          (dipolar = C++-emitted (3cos^2-1)/r^3)

NO Python physics recompute: `dipolar` is read from the C++ producer's emitted
`dipolar_3cos2m1_over_r3` column (saved into the readout npz by sumpool_kernel),
not rebuilt from r and cos_theta. This script is a labeled coefficient/figure
read-out against the producer's own emitted kernel; the symbolic "is it the
Pople form" claim lives in pysr_distill.py.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = np.load("/tmp/rediscover-out/sumpool_kernel_readout.npz")
gsrc, r, ct, inten = d["gsrc"], d["r"], d["ct"], d["inten"]
dipolar = d["dipolar"]            # C++-emitted (3cos^2-1)/r^3 per source

def fit(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    s, b = np.polyfit(x[ok], y[ok], 1)
    yh = s*x[ok]+b
    return s, b, 1-((y[ok]-yh)**2).sum()/((y[ok]-y[ok].mean())**2).sum(), np.corrcoef(x[ok], y[ok])[0,1]

far = r >= 4.0
print(f"far-field sources (r>=4 A): {far.sum()}/{len(r)}")
kern_far = (inten*dipolar)[far]          # intensity * emitted dipolar
s, b, R2, rr = fit(kern_far, gsrc[far])
print(f"  g ~ A * intensity*dipolar(emitted)   A={s:+.5g}  intercept={b:+.4g}  "
      f"R2={R2:+.4f}  r={rr:+.4f}")
s2, b2, R22, rr2 = fit(dipolar[far], gsrc[far])
print(f"  g ~ A * dipolar(emitted, unweighted) A={s2:+.5g}  R2={R22:+.4f}  r={rr2:+.4f}")

# ---- figure (built from the learned read-out and the EMITTED dipolar) ----
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))

# (1) angular law: bin learned g by cos theta in a narrow r band (4-5 A) so r and
# intensity are roughly fixed; overlay the emitted-dipolar angular shape by
# binning the EMITTED dipolar the same way (read, not recomputed).
band = (r >= 4.0) & (r < 5.0)
cb = ct[band]; gb = gsrc[band]; db = dipolar[band]
edges = np.linspace(-1, 1, 16); mids = 0.5*(edges[:-1]+edges[1:])
binned = np.array([gb[(cb >= edges[i]) & (cb < edges[i+1])].mean() for i in range(len(mids))])
binned_dip = np.array([db[(cb >= edges[i]) & (cb < edges[i+1])].mean() for i in range(len(mids))])
ok = np.isfinite(binned) & np.isfinite(binned_dip)
sa, ba, R2a, _ = fit(binned_dip[ok], binned[ok])
ax[0].plot(mids, binned, "o", ms=5, label="binned g (r=4-5 A)")
order = np.argsort(binned_dip[ok])
ax[0].plot(binned_dip[ok][order]*sa+ba, "-", lw=2, label=f"a*emitted_dipolar+b\nR2={R2a:.3f}")
ax[0].set_xlabel("cos theta"); ax[0].set_ylabel("mean learned g")
ax[0].set_title("Angular law (r=4-5 A band)"); ax[0].legend(fontsize=8)

# (2) far-field scatter g vs intensity*emitted dipolar
ax[1].scatter(kern_far, gsrc[far], s=2, alpha=0.15)
xx = np.linspace(np.nanpercentile(kern_far,1), np.nanpercentile(kern_far,99), 50)
ax[1].plot(xx, s*xx+b, "r-", lw=2, label=f"A={s:.3g}\nR2={R2:.3f}")
ax[1].set_xlabel("intensity*dipolar (emitted)"); ax[1].set_ylabel("learned g")
ax[1].set_title("Far-field (r>=4 A): kernel recovery"); ax[1].legend(fontsize=8)

# (3) radial law axial vs in-plane (show real sources)
axial = (ct > 0.7) & far
inplane = (np.abs(ct) < 0.2) & far
for m, lab, c in [(axial, "axial cos>0.7", "C0"), (inplane, "in-plane |cos|<0.2", "C3")]:
    ax[2].scatter(r[m], gsrc[m], s=3, alpha=0.2, color=c, label=lab)
ax[2].set_xlabel("r (A)"); ax[2].set_ylabel("learned g")
ax[2].set_title("Radial falloff"); ax[2].legend(fontsize=8); ax[2].set_xlim(3, 12)

fig.tight_layout()
fig.savefig("/tmp/rediscover-out/ring_kernel_recovered.png", dpi=110)
print("\n  figure -> /tmp/rediscover-out/ring_kernel_recovered.png")

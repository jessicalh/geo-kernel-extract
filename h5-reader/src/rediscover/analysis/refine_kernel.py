#!/usr/bin/env python3
"""refine_kernel — far-field coefficient + a figure of the recovered kernel.

Loads the sum-pooling read-out (g on every real source) and, in the regime
where the point-dipole Pople form is valid (r >= 4 A), fits
    g  ~  A * intensity * (3cos^2-1)/r^3
to recover the single coefficient A. Then draws g overlaid on the Pople kernel.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = np.load("/tmp/rediscover-out/sumpool_kernel_readout.npz")
gsrc, r, ct, inten = d["gsrc"], d["r"], d["ct"], d["inten"]

def fit(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    s, b = np.polyfit(x[ok], y[ok], 1)
    yh = s*x[ok]+b
    return s, b, 1-((y[ok]-yh)**2).sum()/((y[ok]-y[ok].mean())**2).sum(), np.corrcoef(x[ok], y[ok])[0,1]

far = r >= 4.0
print(f"far-field sources (r>=4 A): {far.sum()}/{len(r)}")
kern_far = (inten*(3*ct**2-1)/r**3)[far]
s, b, R2, rr = fit(kern_far, gsrc[far])
print(f"  g ~ A * intensity*(3cos^2-1)/r^3   A={s:+.5g}  intercept={b:+.4g}  "
      f"R2={R2:+.4f}  r={rr:+.4f}")
# pure geometry far-field (fixing out intensity by using only one ring type would be cleaner;
# here intensity varies, so the intensity-weighted form is the honest comparison)
s2, b2, R22, rr2 = fit(((3*ct**2-1)/r**3)[far], gsrc[far])
print(f"  g ~ A * (3cos^2-1)/r^3 (unweighted) A={s2:+.5g}  R2={R22:+.4f}  r={rr2:+.4f}")

# ---- figure (built from the correct per-source read-out arrays) ----
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))

# (1) angular law: bin learned g by cos theta in a narrow r band (4-5 A) so r and
# intensity are roughly fixed; overlay the (3cos^2-1) shape.
band = (r >= 4.0) & (r < 5.0)
cb = ct[band]; gb = gsrc[band]
edges = np.linspace(-1, 1, 16); mids = 0.5*(edges[:-1]+edges[1:])
binned = [gb[(cb >= edges[i]) & (cb < edges[i+1])].mean() for i in range(len(mids))]
binned = np.array(binned)
ok = np.isfinite(binned)
ang = 3*mids**2 - 1
sa, ba, R2a, _ = fit(ang[ok], binned[ok])
ax[0].plot(mids, binned, "o", ms=5, label="binned g (r=4-5 A)")
xs = np.linspace(-1, 1, 60)
ax[0].plot(xs, sa*(3*xs**2-1)+ba, "-", lw=2, label=f"a(3cos^2-1)+b\nR2={R2a:.3f}")
ax[0].set_xlabel("cos theta"); ax[0].set_ylabel("mean learned g")
ax[0].set_title("Angular law (r=4-5 A band)"); ax[0].legend(fontsize=8)

# (2) far-field scatter g vs intensity*(3cos^2-1)/r^3
ax[1].scatter(kern_far, gsrc[far], s=2, alpha=0.15)
xx = np.linspace(np.nanpercentile(kern_far,1), np.nanpercentile(kern_far,99), 50)
ax[1].plot(xx, s*xx+b, "r-", lw=2, label=f"A={s:.3g}\nR2={R2:.3f}")
ax[1].set_xlabel("intensity*(3cos^2-1)/r^3"); ax[1].set_ylabel("learned g")
ax[1].set_title("Far-field (r>=4 A): kernel recovery"); ax[1].legend(fontsize=8)

# (3) radial law axial vs in-plane (from a fresh grid would need the model; show real sources)
axial = (ct > 0.7) & far
inplane = (np.abs(ct) < 0.2) & far
for m, lab, c in [(axial, "axial cos>0.7", "C0"), (inplane, "in-plane |cos|<0.2", "C3")]:
    ax[2].scatter(r[m], gsrc[m], s=3, alpha=0.2, color=c, label=lab)
ax[2].set_xlabel("r (A)"); ax[2].set_ylabel("learned g")
ax[2].set_title("Radial falloff"); ax[2].legend(fontsize=8); ax[2].set_xlim(3, 12)

fig.tight_layout()
fig.savefig("/tmp/rediscover-out/ring_kernel_recovered.png", dpi=110)
print("\n  figure -> /tmp/rediscover-out/ring_kernel_recovered.png")

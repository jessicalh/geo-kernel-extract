#!/usr/bin/env python3
"""make_rama_figs.py — redraw coherent Ramachandran shielding maps from the atlas
dihedral_sigma_bins table (the atlas's own PNGs render blank / as scattered specks;
the binned data underneath is rich). One filled phi/psi heatmap of sigma_mean per
IUPAC context (residue x atom), proper degree axes + ppm colorbar.

Usage:
  make_rama_figs.py --out <dir> [--view pooled] [--only GLY:N,VAL:N,...] [--limit N]
"""
import argparse, csv, collections, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TBL = os.environ.get("ATLAS_TABLES", "/shared/2026Thesis/nmr-shielding-emit-build/output/contribution_atlas_20260612T000000Z/tables") + "/dihedral_sigma_bins.csv"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--view", default="pooled")          # pooled | between_720 | within_1p9j
    ap.add_argument("--only", default="")                # RES:ATOM,RES:ATOM
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--min-n", type=int, default=3)       # drop bins with fewer rows (despeckle)
    args = ap.parse_args()
    only = set(x.strip().upper() for x in args.only.split(",") if x.strip())
    os.makedirs(args.out, exist_ok=True)

    rows = []
    with open(TBL) as f:
        for d in csv.DictReader(f):
            if (d["granularity"] == "per_residue_type" and d["dihedral"] == "phi_psi"
                    and d["dataset_scope"] == args.view):
                rows.append(d)

    xb = [int(float(d["x_bin"])) for d in rows]
    yb = [int(float(d["y_bin"])) for d in rows]
    nb = max(max(xb), max(yb)) + 1                        # square binning over [-180,180]
    edges = np.linspace(-180, 180, nb + 1)

    byctx = collections.defaultdict(list)
    for d in rows:
        atom = d["nucleus_id"].split("|")[1]
        byctx[(d["residue_scope"], atom)].append(d)

    made = 0
    for (res, atom) in sorted(byctx):
        if only and f"{res}:{atom}" not in only:
            continue
        if args.limit and made >= args.limit:
            break
        ds = byctx[(res, atom)]
        grid = np.full((nb, nb), np.nan)
        nfilled = 0
        for d in ds:
            try:
                if int(float(d["n_rows"])) < args.min_n:
                    continue
                x = int(float(d["x_bin"])); y = int(float(d["y_bin"]))
                grid[y, x] = float(d["sigma_mean"]); nfilled += 1
            except Exception:
                pass
        if nfilled < 8:
            continue
        vals = grid[~np.isnan(grid)]
        lo, hi = np.nanpercentile(vals, 2.5), np.nanpercentile(vals, 97.5)
        ctr = np.nanmedian(vals); half = max(hi - ctr, ctr - lo, 1e-6)

        fig, ax = plt.subplots(figsize=(4.7, 4.0))
        cmap = plt.cm.RdBu_r.copy(); cmap.set_bad("#f2f2f2")
        pcm = ax.pcolormesh(edges, edges, np.ma.masked_invalid(grid), cmap=cmap,
                            vmin=ctr - half, vmax=ctr + half, shading="flat")
        ax.set_xlim(-180, 180); ax.set_ylim(-180, 180); ax.set_aspect("equal")
        ax.set_xticks([-180, -90, 0, 90, 180]); ax.set_yticks([-180, -90, 0, 90, 180])
        ax.set_xlabel(r"$\varphi$ (deg)"); ax.set_ylabel(r"$\psi$ (deg)")
        ax.axhline(0, color="k", lw=0.3, alpha=0.3); ax.axvline(0, color="k", lw=0.3, alpha=0.3)
        ax.set_title(f"{res}:{atom}   $\\sigma$ vs $\\varphi/\\psi$   ({args.view}, {nfilled} bins)",
                     fontsize=9)
        cb = fig.colorbar(pcm, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label(r"$\sigma$ (ppm)", fontsize=8); cb.ax.tick_params(labelsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(args.out, f"{res}_{atom}.png"), dpi=130)
        plt.close(fig)
        made += 1
    print(f"wrote {made} rama figures to {args.out}  (view={args.view}, nb={nb})")

if __name__ == "__main__":
    main()

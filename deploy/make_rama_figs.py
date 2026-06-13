#!/usr/bin/env python3
"""make_rama_figs.py — coherent shielding maps from the atlas dihedral_sigma_bins
table (the atlas's own PNGs render blank / specky; the binned data is rich).

sigma here is the DFT **absolute isotropic shielding** <sigma_iso> (orca_total_T0,
r2SCAN) in ppm — NOT chemical shift, NOT a contribution. Figures are labelled per
nucleus so that is unambiguous (the "is that ppm real?" question: yes, absolute
sigma, which is why a 13C map reads 0-160 rather than the 0-60 of delta).

Two figure families per IUPAC context (residue x atom):

  A. phi/psi heatmap of <sigma_iso>, colour SHARED PER IUPAC ATOM NAME (all CA on
     one scale, all O on another) so strength is comparable across residues and
     the atoms of a residue stop looking alike. phi/psi is a BACKBONE coordinate,
     identical for every atom of a residue, so the occupancy outline is shared
     across a residue's atoms — only sigma differs.

  B. sidechain atoms only: sigma-vs-chi1/chi2 marginals — the axis on which a
     sidechain actually moves. phi/psi cannot differentiate sidechain dynamics.

SUPPORT IS HONEST (the "n=" work): every cell is shown down to n=1, but per-bin
support n is mapped to OPACITY — well-sampled cells are solid, sparse ones fade to
the grey background. So nothing is hidden by a hard count cutoff, and nothing
sparse masquerades as solid. (--n-full is the count at which a cell reads fully
opaque; sqrt scaling with a small floor keeps n=1 faintly visible.) On the chi
panels, support is the marker SIZE. This is physics with some confounding DFT
water, so the aggregate matters but the per-bin support has to be legible too.

Usage: make_rama_figs.py --out <dir> [--view pooled] [--only RES:ATOM,...]
                         [--limit N] [--min-n 1] [--n-full 30]
"""
import argparse, csv, collections, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

TBL = os.environ.get(
    "ATLAS_TABLES",
    "/shared/2026Thesis/nmr-shielding-emit-build/output/contribution_atlas_20260612T000000Z/tables",
) + "/dihedral_sigma_bins.csv"

BACKBONE_ATOMS = {"N", "CA", "C", "O", "H", "HA", "HA2", "HA3",
                  "H1", "H2", "H3", "OXT", "OT1", "OT2"}
NUCLEUS = {"C": r"$^{13}$C", "N": r"$^{15}$N", "O": r"$^{17}$O",
           "H": r"$^{1}$H", "S": r"$^{33}$S"}
ALPHA_FLOOR = 0.12  # n=1 stays faintly visible rather than invisible


def nucleus_of(atom):
    return NUCLEUS.get(atom[:1], atom[:1])


def n_weighted_cells(rows):
    """Collapse ss3/ss8-split rows to an n-weighted mean per bin cell.
    Returns {(x, y): (mean, total_n)}; y is None for 1D (chi) dihedrals."""
    acc = {}
    for d in rows:
        try:
            nr = int(float(d["n_rows"]))
            x = int(float(d["x_bin"]))
            yraw = d.get("y_bin", "")
            y = int(float(yraw)) if yraw not in ("", None) else None
            sm = float(d["sigma_mean"])
        except Exception:
            continue
        a = acc.setdefault((x, y), [0.0, 0])
        a[0] += sm * nr
        a[1] += nr
    return {k: (s / n, n) for k, (s, n) in acc.items() if n > 0}


def support_alpha(n, n_full):
    return float(np.clip(np.sqrt(max(n, 0) / float(n_full)), ALPHA_FLOOR, 1.0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--view", default="pooled")
    ap.add_argument("--only", default="")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--min-n", type=int, default=1)      # show all support; honesty via opacity
    ap.add_argument("--n-full", type=int, default=30)    # bin count that reads fully opaque
    args = ap.parse_args()
    only = set(x.strip().upper() for x in args.only.split(",") if x.strip())
    os.makedirs(args.out, exist_ok=True)
    nfull = max(1, args.n_full)

    # ---- load ----
    pp = collections.defaultdict(list)
    chi = collections.defaultdict(lambda: collections.defaultdict(list))
    nb = nbins_chi = 0
    with open(TBL) as f:
        for d in csv.DictReader(f):
            if d["dataset_scope"] != args.view or d["granularity"] != "per_residue_type":
                continue
            nid = d["nucleus_id"].split("|")
            if len(nid) < 2:
                continue
            atom = nid[1]
            dih = d["dihedral"]
            if dih == "phi_psi":
                pp[(d["residue_scope"], atom)].append(d)
                try:
                    nb = max(nb, int(float(d["x_bin"])) + 1, int(float(d["y_bin"])) + 1)
                except Exception:
                    pass
            elif dih in ("chi1", "chi2"):
                chi[(d["residue_scope"], atom)][dih].append(d)
                try:
                    nbins_chi = max(nbins_chi, int(float(d["x_bin"])) + 1)
                except Exception:
                    pass

    # ---- pass 1: phi/psi cells + shared per-atom-name colour scale ----
    pp_cells = {}
    sig_by_atom = collections.defaultdict(list)
    for key, rows in pp.items():
        cells = {k: v for k, v in n_weighted_cells(rows).items()
                 if k[1] is not None and v[1] >= args.min_n}
        pp_cells[key] = cells
        sig_by_atom[key[1]].extend(m for (m, n) in cells.values())

    scale = {}
    for atom, vals in sig_by_atom.items():
        v = np.asarray(vals, float); v = v[np.isfinite(v)]
        if v.size == 0:
            continue
        ctr = float(np.median(v))
        half = float(max(np.percentile(v, 98) - ctr, ctr - np.percentile(v, 2), 1e-6))
        scale[atom] = (ctr, half)

    # ---- pass 2: draw ----
    made_pp = made_chi = 0
    for (res, atom) in sorted(set(pp_cells) | set(chi)):
        if only and f"{res}:{atom}" not in only:
            continue
        if args.limit and made_pp >= args.limit:
            break
        nuc = nucleus_of(atom)

        # A. phi/psi heatmap: colour = sigma (shared per-atom scale), alpha = support n
        cells = pp_cells.get((res, atom), {})
        if len(cells) >= 8:
            ctr, half = scale.get(atom, (0.0, 1.0))
            norm = Normalize(vmin=ctr - half, vmax=ctr + half)
            cmap = plt.cm.RdBu_r
            rgba = np.zeros((nb, nb, 4))
            ntot = 0
            for (x, y), (m, n) in cells.items():
                rgba[y, x] = cmap(norm(m))
                rgba[y, x, 3] = support_alpha(n, nfull)
                ntot += n
            fig, ax = plt.subplots(figsize=(4.7, 4.1))
            ax.set_facecolor("#f2f2f2")
            ax.imshow(rgba, origin="lower", extent=[-180, 180, -180, 180],
                      aspect="equal", interpolation="nearest")
            ax.set_xlim(-180, 180); ax.set_ylim(-180, 180)
            ax.set_xticks([-180, -90, 0, 90, 180]); ax.set_yticks([-180, -90, 0, 90, 180])
            ax.set_xlabel(r"$\varphi$ (deg)"); ax.set_ylabel(r"$\psi$ (deg)")
            ax.axhline(0, color="k", lw=0.3, alpha=0.3); ax.axvline(0, color="k", lw=0.3, alpha=0.3)
            ax.set_title(f"{res}:{atom} ({nuc})  abs. isotropic shielding vs $\\varphi/\\psi$\n"
                         f"{args.view} · {len(cells)} bins, n={ntot} · shared {atom} scale · "
                         f"opacity = support (n$\\geq${nfull} solid)", fontsize=7.5)
            sm = ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
            cb = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
            cb.set_label(r"$\langle\sigma_{\mathrm{iso}}\rangle$  (ppm, absolute)", fontsize=8)
            cb.ax.tick_params(labelsize=7)
            fig.tight_layout()
            fig.savefig(os.path.join(args.out, f"{res}_{atom}.png"), dpi=130)
            plt.close(fig)
            made_pp += 1

        # B. sidechain sigma-vs-chi marginals: marker SIZE = support n
        if atom not in BACKBONE_ATOMS and (res, atom) in chi and nbins_chi > 0:
            series = []
            for dih in ("chi1", "chi2"):
                rows = chi[(res, atom)].get(dih)
                if not rows:
                    continue
                cc = {k: v for k, v in n_weighted_cells(rows).items() if v[1] >= args.min_n}
                pts = sorted((x, m, n) for (x, _y), (m, n) in cc.items())
                if len(pts) >= 4:
                    ang = [-180 + (x + 0.5) * 360.0 / nbins_chi for x, _, _ in pts]
                    series.append((dih, ang, [m for _, m, _ in pts], [n for _, _, n in pts]))
            if series:
                fig, ax = plt.subplots(figsize=(4.7, 3.5))
                for dih, ang, sm_, nn in series:
                    sizes = np.clip(np.array(nn, float) / nfull * 45.0 + 6.0, 6.0, 80.0)
                    line, = ax.plot(ang, sm_, lw=0.8, alpha=0.5, label=dih)
                    ax.scatter(ang, sm_, s=sizes, color=line.get_color(), alpha=0.9, edgecolors="none")
                if atom in scale:
                    ax.axhline(scale[atom][0], color="k", lw=0.5, alpha=0.35, ls="--",
                               label=f"all-{atom} median")
                ax.set_xlim(-180, 180); ax.set_xticks([-180, -90, 0, 90, 180])
                ax.set_xlabel(r"$\chi$ (deg)")
                ax.set_ylabel(r"$\langle\sigma_{\mathrm{iso}}\rangle$  (ppm)")
                ax.set_title(f"{res}:{atom} ({nuc})  abs. isotropic shielding vs sidechain $\\chi$\n"
                             f"{args.view} · marker size = support (n)", fontsize=7.5)
                ax.legend(fontsize=7, frameon=False)
                ax.grid(True, lw=0.3, alpha=0.3)
                fig.tight_layout()
                fig.savefig(os.path.join(args.out, f"{res}_{atom}_chi.png"), dpi=130)
                plt.close(fig)
                made_chi += 1

    print(f"wrote {made_pp} phi/psi + {made_chi} sidechain-chi figures to {args.out} "
          f"(view={args.view}, nb={nb}, nbins_chi={nbins_chi}, n_full={nfull})")


if __name__ == "__main__":
    main()

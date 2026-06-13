#!/usr/bin/env python3
"""make_iupac_report.py — per-IUPAC LaTeX report from the contribution atlas.

Per IUPAC nucleus (residue x atom), per view, ONE merged table:
  Contributor (plain name) | Mechanism | R | R^2 | redundancy group
with a credibility header (sigma spread in ppm, n_proteins, effective independent
dimensions among the listed contributors) so no R^2 is read out of context. The
redundancy group + effective-dim count ARE the linear independence (collinear
contributors share a group; from compute_li.py's sigma-prediction correlations).
Then a redrawn sigma-vs-(phi,psi) Ramachandran map. A grounded glossary closes it.

Usage: make_iupac_report.py --out <dir> --figdir <rama dir> --li <li_data.json> [--only RES,..] [--limit N]
"""
import argparse, csv, collections, json, os, glob

T = os.environ.get("ATLAS_TABLES", "/shared/2026Thesis/nmr-shielding-emit-build/output/contribution_atlas_20260612T000000Z/tables")
AA = ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS",
      "MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]
BB = ["N","CA","C","O","H","HA"]

ACR = {"efg":"EFG","mopac":"MOPAC","aimnet2":"AIMNet2","apbs":"APBS","bs":"BS","hm":"HM",
       "dssp":"DSSP","eeq":"EEQ","sasa":"SASA","ff":"FF","mc":"McConnell","co":"C=O",
       "bo":"bond-order","t0":"T0","t1":"T1","t2":"T2","aim":"AIM","sd":"SD"}

# (prefix, tool + what it observes) -- checked in order, first match wins
TOOL = [
 ("mopac_coulomb", "MOPAC PM7 (semiempirical QM): Coulomb field/energy from the converged charge density"),
 ("mopac", "MOPAC PM7 (semiempirical QM): wavefunction-derived charges / energetics"),
 ("aimnet2", "AIMNet2 (ML interatomic potential): predicted charges, charge response, or embedding"),
 ("apbs", "APBS (Poisson-Boltzmann continuum electrostatics): electric field at the atom"),
 ("bs", "Biot-Savart ring-current kernel (ours): magnetic field from aromatic ring currents"),
 ("hm", "Haigh-Mallion ring-current kernel (ours): ring-current geometry factor"),
 ("mc", "McConnell bond-anisotropy kernel (ours): anisotropic susceptibility of nearby bonds"),
 ("ff_coulomb", "Force-field point charges: Coulomb field / energy"),
 ("ff", "Force-field (MM) parameters: per-atom partial charge / Poisson-Boltzmann radius"),
 ("coulomb", "Force-field point charges: Coulomb field / energy"),
 ("eeq", "EEQ charge model: electronegativity-equilibrated partial charges / coordination number"),
 ("dssp", "DSSP: secondary-structure assignment and backbone H-bond geometry"),
 ("gromacs", "GROMACS: per-frame simulation energy"),
 ("bonded", "GROMACS bonded terms: per-atom bonded-energy share"),
 ("atom_sasa", "Shrake-Rupley: solvent-accessible surface area"),
 ("water", "Explicit solvent: local water field / polarization"),
 ("disp", "Dispersion (D4-type): per-type dispersion term"),
 ("pq", "pi-quadrupole kernel (ours): aromatic quadrupole geometry factor (3cos^2theta-1)/r^4 at the site"),
 ("larsen", "Larsen ProCS15-grid H-bond term (ours): per-class backbone N-H / C=O hydrogen-bond shielding from the ProCS15 lookup"),
 ("hbond", "H-bond kernel (ours): donor-H..acceptor geometry, the angular/distance scalar (3cos^2theta-1)/r^3 at the site"),
 ("omega", "Geometry: backbone omega dihedral / deviation"),
 ("phi", "Geometry: backbone phi"), ("psi", "Geometry: backbone psi"),
 ("pyramid", "Geometry: pyramidalisation"), ("planar", "Geometry: planarity deviation"),
 ("ring_geometry", "Geometry: aromatic-ring geometry"),
 ("residue", "Topology: residue identity / index"), ("element", "Topology: element"),
]
PHYS = {
 "electrostatic_efg": "A local electric field and its gradient polarize the electron density at the nucleus, shifting the paramagnetic shielding.",
 "ring_current": "Circulating pi-electron currents in aromatic rings make a secondary magnetic field that adds to or opposes B0 at nearby nuclei (the ring-current shielding/deshielding cone).",
 "bond_anisotropy": "Nearby bonds have anisotropic magnetic susceptibility (McConnell); their orientation relative to the nucleus shifts sigma.",
 "charges": "Local partial charge / electron density sets the diamagnetic-paramagnetic balance of the shielding.",
 "secondary_structure": "Backbone H-bonding and the phi/psi geometry of the secondary-structure element correlate with sigma.",
 "solvation": "Solvent exposure and the solvent reaction field perturb the electron density and sigma.",
 "geometry": "Local geometric distortion (planarity, dihedral, pyramidalisation) modulates the electronic structure and sigma.",
 "gromacs_runtime": "Whole-frame energetics -- an indirect, weak proxy for the structural state.",
 "ring_efg": "Electric field / gradient contributed specifically by aromatic-ring atoms.",
 "ring_dispersion": "Dispersion interaction with aromatic rings (a weak van der Waals proxy).",
 "hbond_kernel": "Hydrogen-bond geometry / strength at the site, which shifts sigma (especially for H, N, O).",
 "hbond_grid": "Backbone N-H / C=O hydrogen-bond geometry mapped through the ProCS15 empirical shielding surface (a per-class lookup), which shifts sigma at backbone N, H, and C=O.",
 "topology": "An identity or positional label (residue, element, locant) -- a descriptor, not a physical field.",
}

def esc(s):
    s = str(s)
    for a, b in [("\\", r"\textbackslash{}"), ("_", r"\_"), ("#", r"\#"), ("%", r"\%"),
                 ("&", r"\&"), ("$", r"\$"), ("{", r"\{"), ("}", r"\}"),
                 ("~", r"\textasciitilde{}"), ("^", r"\textasciicircum{}")]:
        s = s.replace(a, b)
    return s

def pretty(cid):
    return " ".join(ACR.get(t.lower(), t) for t in cid.split("_"))

def tool_of(cid):
    for pre, desc in TOOL:
        if cid.startswith(pre):
            return desc
    return "--"

def fnum(x):
    try: return float(x)
    except: return None

def groups_from_corr(corr, thr=0.85):
    n = len(corr)
    adj = [[i] for i in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if abs(corr[i][j]) >= thr:
                adj[i].append(j); adj[j].append(i)
    seen = {}; gid = 0; labels = [None] * n
    for i in range(n):
        if i in seen: continue
        stack = [i]; comp = []
        while stack:
            u = stack.pop()
            if u in seen: continue
            seen[u] = True; comp.append(u)
            for v in adj[u]:
                if v not in seen: stack.append(v)
        lab = chr(ord("A") + gid); gid += 1
        for u in comp: labels[u] = lab
    return labels

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--figdir", required=True)
    ap.add_argument("--li", required=True)
    ap.add_argument("--only", default="")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()
    only = set(x.strip().upper() for x in args.only.split(",") if x.strip())

    li = json.load(open(args.li))
    # score: (res,atom,view) -> {cid: (R,R2,neff)} ; and cid->mechanism
    score = collections.defaultdict(dict); mech = {}
    for d in csv.DictReader(open(T + "/contributor_score_by_nucleus.csv")):
        if d["granularity"] != "per_residue_type": continue
        cid = d["contributor_id"]; mech[cid] = d["mechanism"]
        score[(d["residue_scope"], d["nucleus_id"].split("|")[1], d["dataset_scope"])][cid] = (
            fnum(d["R_oof"]), fnum(d["R2_oof"]), d["n_eff"])

    ctxs = []
    for p in glob.glob(os.path.join(args.figdir, "*.png")):
        base = os.path.basename(p)[:-4]
        if "_" not in base: continue
        res, atom = base.split("_", 1)
        if only and res not in only: continue
        ctxs.append((res, atom, p))
    byres = collections.defaultdict(list)
    for res, atom, p in ctxs:
        rank = (0, BB.index(atom)) if atom in BB else (1, atom)
        byres[res].append((rank, atom, p))

    used_cids = set()

    def view_table(res, atom, view, label):
        rec = li.get(f"{res}|{atom}|{view}")
        if not rec or not rec.get("contributors"):
            return []
        cids = rec["contributors"]
        labels = groups_from_corr(rec["corr"]) if rec.get("corr") and len(rec["corr"]) == len(cids) else None
        sc = score.get((res, atom, view), {})
        hdr = []
        ss = rec.get("sigma_sd"); npr = rec.get("n_proteins"); eff = rec.get("eff_dims")
        bits = []
        if ss is not None: bits.append(r"$\sigma$ spread %.1f ppm" % ss)
        if view == "between_720" and npr: bits.append(r"%d proteins" % npr)
        if view == "within_1p9j": bits.append(r"%d atom-frames" % rec.get("n_atoms", 0))
        if eff is not None: bits.append(r"%.1f eff.\ dims of %d" % (eff, len(cids)))
        hdr.append(r"\textbf{%s} \small{(%s)}\\[1pt]" % (label, " \\textbullet\\ ".join(bits)))
        usegrp = labels is not None
        cols = "llrr" + ("c" if usegrp else "")
        head = r"Contributor & Mechanism & $R$ & $R^2$" + (r" & grp" if usegrp else "") + r"\\"
        hdr.append(r"\begin{tabular}{%s}\toprule %s\midrule" % (cols, head))
        for i, cid in enumerate(cids):
            used_cids.add(cid)
            R, R2, ne = sc.get(cid, (None, None, ""))
            row = r"%s & %s & %s & %s" % (
                esc(pretty(cid)), esc(mech.get(cid, "")),
                ("%.3f" % R) if R is not None else "--",
                ("%.3f" % R2) if R2 is not None else "--")
            if usegrp: row += r" & %s" % labels[i]
            hdr.append(row + r"\\")
        hdr.append(r"\bottomrule\end{tabular}\\[5pt]")
        return hdr

    L = [r"""\documentclass[10pt]{article}
\usepackage[a4paper,margin=2cm]{geometry}
\usepackage{graphicx,booktabs,longtable,float,caption}
\usepackage{helvet}\renewcommand{\familydefault}{\sfdefault}
\usepackage[T1]{fontenc}
\usepackage{hyperref}\hypersetup{colorlinks=true,linkcolor=blue}
\captionsetup{font=small}\setlength{\parindent}{0pt}
\title{\textbf{NMR Shielding Contribution Atlas}}
\date{2026-06-12}\author{}
\begin{document}\maketitle
"""]

    n = 0; stop = False
    for res in sorted(byres, key=lambda r: AA.index(r) if r in AA else 99):
        L.append(r"\section*{%s}" % esc(res))
        L.append(r"\addcontentsline{toc}{section}{%s}" % esc(res))
        for _, atom, fig in sorted(byres[res], key=lambda t: t[0]):
            if args.limit and n >= args.limit: stop = True; break
            L.append(r"\subsection*{%s\,:\,%s}" % (esc(res), esc(atom)))
            L += view_table(res, atom, "between_720", "Between proteins (720)")
            L += view_table(res, atom, "within_1p9j", "Within protein (1P9J)")
            if os.path.exists(fig):
                L.append(r"\begin{figure}[H]\centering\includegraphics[width=0.58\textwidth]{\detokenize{%s}}" % fig)
                L.append(r"\caption{$\sigma$ vs.\ $(\varphi,\psi)$ --- %s:%s (pooled).}\end{figure}" % (esc(res), esc(atom)))
            n += 1
        if stop: break

    # glossary
    L.append(r"\clearpage\section*{Glossary of input fields}")
    L.append(r"\small Contributors that appear in this atlas: the field, the tool that observed it "
             r"and what that tool measures, and the physical mechanism relating it to the shielding $\sigma$. "
             r"Linear-independence groups are empirical (sample collinearity), not claims of physical independence.\\[4pt]")
    L.append(r"\begin{longtable}{p{3.5cm}p{6.0cm}p{6.0cm}}\toprule "
             r"Field & Tool \& what it observes & Physics link to $\sigma$\\\midrule\endhead")
    prev_tool = prev_phys = None
    for cid in sorted(used_cids, key=lambda c: (mech.get(c, "zzz"), c)):
        idtex = cid.replace("_", r"\_\allowbreak{}")   # break long names at underscores
        tool = tool_of(cid); phys = PHYS.get(mech.get(cid, ""), "--")
        # collapse a cell to "as above" when identical to the row directly above
        tool_cell = r"\textit{as above}" if tool == prev_tool else esc(tool)
        phys_cell = r"\textit{as above}" if phys == prev_phys else esc(phys)
        L.append(r"{\scriptsize\texttt{%s}}\newline %s & %s & %s\\[2pt]" % (
            idtex, esc(pretty(cid)), tool_cell, phys_cell))
        prev_tool, prev_phys = tool, phys
    L.append(r"\bottomrule\end{longtable}")

    L.append(r"\end{document}")
    os.makedirs(args.out, exist_ok=True)
    tex = os.path.join(args.out, "iupac_atlas.tex")
    with open(tex, "w") as f:
        f.write("\n".join(L))
    print(f"wrote {tex}  ({n} IUPAC contexts, {len(used_cids)} glossary fields)")

if __name__ == "__main__":
    main()

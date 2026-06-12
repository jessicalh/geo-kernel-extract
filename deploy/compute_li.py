#!/usr/bin/env python3
"""compute_li.py — per-IUPAC linear independence + credibility controls, from the substrate.

For each (residue, atom, view): select the context's atoms (DFT target present),
fit each top contributor's OLS prediction of sigma_T0, and report:
  * the correlation among those sigma-predictions (the real collinearity),
  * effective independent dimensions (participation ratio of the corr eigenvalues),
  * credibility controls: sigma spread (ppm), n_proteins, n_atoms.
Writes li_data.json keyed "RES|ATOM|view".
"""
import numpy as np, csv, collections, json, os

E = "/shared/2026Thesis/nmr-shielding-emit-build/output/consolidated_fp_emit_20260612T000000Z"
T = "/shared/2026Thesis/nmr-shielding-emit-build/output/contribution_atlas_20260612T000000Z/tables"
CACHE = "/tmp/li_rows_cache.npz"
OUT = "/shared/2026Thesis/atlas_report_20260612/li_data.json"

def fnum(x):
    try: return float(x)
    except: return None

def load_rows():
    if os.path.exists(CACHE):
        z = np.load(CACHE, allow_pickle=True)
        return z["res"], z["iup"], z["dsid"], z["prot"], z["tgt"]
    res=[]; iup=[]; dsid=[]; prot=[]; tgt=[]
    for d in csv.DictReader(open(E + "/rows.csv")):
        res.append(d["residue_type"]); iup.append(d["iupac_name"]); dsid.append(d["dataset_id"])
        prot.append(d["protein_id"]); tgt.append(d["target_total_present"] == "1")
    res=np.array(res); iup=np.array(iup); dsid=np.array(dsid); prot=np.array(prot); tgt=np.array(tgt, bool)
    np.savez(CACHE, res=res, iup=iup, dsid=dsid, prot=prot, tgt=tgt)
    return res, iup, dsid, prot, tgt

def main():
    res, iup, dsid, prot, tgt = load_rows()
    sigma = np.asarray(np.load(E + "/targets/orca_total_T0.npy"), float)
    dsl = dsid.astype(str)
    views = {"between_720": np.char.find(dsl, "720") >= 0,
             "within_1p9j": (np.char.find(dsl, "1p9j") >= 0) | (np.char.find(dsl, "1P9J") >= 0)}

    top = collections.defaultdict(list)
    for d in csv.DictReader(open(T + "/contributor_score_by_nucleus.csv")):
        if d["granularity"] != "per_residue_type": continue
        r2, q = fnum(d["R2_oof"]), fnum(d["q_value"])
        if r2 and r2 > 0 and (q is None or q < 0.05):
            top[(d["residue_scope"], d["nucleus_id"].split("|")[1], d["dataset_scope"])].append((r2, d["contributor_id"]))

    def yhat(idx, cid):
        p = E + f"/measures/{cid}.npy"
        if not os.path.exists(p): return None
        X = np.asarray(np.load(p, mmap_mode='r')[idx], float)
        if X.ndim == 1: X = X[:, None]
        y = sigma[idx]
        ok = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
        if ok.sum() < 20: return None
        A = np.hstack([np.ones((int(ok.sum()), 1)), X[ok]])
        beta, *_ = np.linalg.lstsq(A, y[ok], rcond=None)
        yh = np.full(len(idx), np.nan); yh[ok] = A @ beta
        return yh

    out = {}
    ctxs = sorted(set((r, a) for (r, a, v) in top))
    for (rt, atom) in ctxs:
        for view, mask in views.items():
            idx = np.where((res == rt) & (iup == atom) & tgt & mask)[0]
            if len(idx) < 20: continue
            y = sigma[idx]
            rec = {"n_atoms": int(len(idx)), "n_proteins": int(len(set(prot[idx].tolist()))),
                   "sigma_sd": float(np.nanstd(y))}
            cids = [c for _, c in sorted(top[(rt, atom, view)], reverse=True)][:6]
            Y = []; labs = []
            for c in cids:
                yh = yhat(idx, c)
                if yh is not None: Y.append(yh); labs.append(c)
            rec["contributors"] = labs
            if len(Y) >= 2:
                Y = np.array(Y); ok = np.all(np.isfinite(Y), axis=0)
                if ok.sum() >= 20:
                    C = np.corrcoef(Y[:, ok])
                    ev = np.linalg.eigvalsh(C); ev = ev[ev > 1e-9]
                    rec["eff_dims"] = round(float((ev.sum() ** 2) / (ev ** 2).sum()), 2)
                    rec["corr"] = [[round(float(C[i, j]), 2) for j in range(len(labs))] for i in range(len(labs))]
            out[f"{rt}|{atom}|{view}"] = rec
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(out, open(OUT, "w"))
    print(f"wrote {OUT} with {len(out)} context-views")

if __name__ == "__main__":
    main()

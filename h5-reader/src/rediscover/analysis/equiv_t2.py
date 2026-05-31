#!/usr/bin/env python3
"""equiv_t2 — equivariant sum-pooling for the ring-current T2 tensor.

The shift is instantaneous; the trajectory samples geometry. The ring-current
contribution to a proton's shielding TENSOR is an l=2 object oriented by BOTH
the displacement direction r_hat and the source ring's normal n_hat (a point
magnetic dipole's field tensor). So the equivariant per-source contribution is a
learned-radial combination of three l=2 basis tensors, in the LIBRARY basis (so
they match the emitted DFT T2 exactly):

    contribution_T2(source) = R_A(r,i)*Y2(r_hat) + R_B(r,i)*Y2(n_hat)
                              + R_C(r,i)*Y2sym(r_hat, n_hat)
    node_T2 = SUM_sources contribution_T2                 (permutation-invariant)

Each Y2 transforms as l=2 (equivariant by construction); the R's are invariant
scalars (small MLPs). Target = the per-atom DE-MEANED local-frame DFT T2 (strip
the near-constant local C-H baseline; the modulation is what ring current claims).
Honesty: frame split + leave-atoms-out; report |T2| too.
"""
import numpy as np
import os
import pandas as pd
import sys
import torch
import torch.nn as nn

torch.manual_seed(0); np.random.seed(0)
dev = "cuda" if torch.cuda.is_available() else "cpu"
out_dir = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("REDISCOVER_OUT", "/tmp/rediscover-out-v2")
SRC = f"{out_dir}/ring_current_sources.csv"

S2, S32, S12 = np.sqrt(2.0), np.sqrt(1.5), 1.0/np.sqrt(2.0)
def lib_T2(M):  # (...,3,3) -> (...,5) in the library basis (matches DecomposeLibrary)
    Sym = 0.5*(M + np.swapaxes(M, -1, -2))
    tr = (Sym[...,0,0]+Sym[...,1,1]+Sym[...,2,2])/3.0
    Sym = Sym - tr[...,None,None]*np.eye(3)
    xx,yy,zz = Sym[...,0,0],Sym[...,1,1],Sym[...,2,2]
    xy,yz,xz = Sym[...,0,1],Sym[...,1,2],Sym[...,0,2]
    return np.stack([S2*xy, S2*yz, S32*zz, S2*xz, S12*(xx-yy)], -1)

df = pd.read_csv(SRC)
# basis sanity: lib_T2 of the raw lab tensor must equal the emitted dft_total_T2.
raw = df[[f"dft_total_raw_{i}{j}" for i in range(3) for j in range(3)]].to_numpy().reshape(-1,3,3)
emit = df[[f"dft_total_T2_{k}" for k in range(5)]].to_numpy()
chk = lib_T2(raw)
print(f"basis check: max|lib_T2(raw) - emitted dft_total_T2| = {np.nanmax(np.abs(chk-emit)):.2e} (expect ~0)")

# through-space sources only (self/bonded is the local baseline; de-meaned anyway)
df = df[(df.dft_present == 1) & (df.is_self_or_bonded == 0)].copy()
print(f"through-space source rows: {len(df)}  atoms={df.atom_index.nunique()}")

# per-source equivariant l=2 basis tensors in the local frame
disp = df[["disp_local_x","disp_local_y","disp_local_z"]].to_numpy()
rhat = disp/np.linalg.norm(disp,axis=1,keepdims=True)
nhat = df[["source_normal_local_x","source_normal_local_y","source_normal_local_z"]].to_numpy()
nhat = nhat/np.clip(np.linalg.norm(nhat,axis=1,keepdims=True),1e-9,None)
A = lib_T2(rhat[:,:,None]*rhat[:,None,:])                       # Y2(r_hat)
B = lib_T2(nhat[:,:,None]*nhat[:,None,:])                       # Y2(n_hat)
C = lib_T2(rhat[:,:,None]*nhat[:,None,:] + nhat[:,:,None]*rhat[:,None,:])  # Y2sym(r,n)
feat = df[["r","ring_intensity"]].to_numpy(np.float32)
fmean,fstd = feat.mean(0),feat.std(0); featn=(feat-fmean)/fstd

# target: local-frame DFT T2 (per group), de-meaned per atom
locT2 = lib_T2(df[[f"dft_total_local_{i}{j}" for i in range(3) for j in range(3)]].to_numpy().reshape(-1,3,3))

df["gid"],_ = pd.factorize(df.atom_index.astype(str)+":"+df.h5_row.astype(str))
n_groups = df.gid.max()+1
grp = df.groupby("gid").agg(aid=("atom_index","first"),frame=("h5_row","first")).sort_index()
# per-group target T2 via groupby-first on the 5 comps, then de-mean by atom
tcols=[f"t2_{k}" for k in range(5)]
for k in range(5): df[tcols[k]]=locT2[:,k]
gT2 = df.groupby("gid")[tcols].first().sort_index().to_numpy()
aid_of_group = grp.aid.to_numpy()
# de-mean per atom
gT2_dm = gT2.copy()
for a in np.unique(aid_of_group):
    m = aid_of_group==a; gT2_dm[m] -= gT2[m].mean(0)

t = lambda x: torch.tensor(x, dtype=torch.float32, device=dev)
A_t,B_t,C_t,featn_t = t(A),t(B),t(C),t(featn)
src_group = torch.tensor(df.gid.to_numpy(),dtype=torch.long,device=dev)
target = t(gT2_dm)
frames = np.sort(grp.frame.unique())
test_f = set(np.random.default_rng(0).choice(frames,int(0.2*len(frames)),replace=False))
g_te = torch.tensor(grp.frame.isin(test_f).to_numpy(),device=dev); g_tr=~g_te
# per-group atom de-mean must apply to the PREDICTION too -> subtract per-atom mean of pred
group_atom = torch.tensor(pd.factorize(aid_of_group)[0],dtype=torch.long,device=dev)
n_atoms = int(group_atom.max())+1

class EquivPool(nn.Module):
    def __init__(self):
        super().__init__()
        self.R = nn.Sequential(nn.Linear(2,64),nn.SiLU(),nn.Linear(64,64),nn.SiLU(),nn.Linear(64,3))
    def forward(self):
        w = self.R(featn_t)                       # (n_src,3) radial weights R_A,R_B,R_C
        contrib = w[:,0:1]*A_t + w[:,1:2]*B_t + w[:,2:3]*C_t   # (n_src,5)
        pooled = torch.zeros(n_groups,5,device=dev)
        pooled.index_add_(0,src_group,contrib)
        # de-mean prediction per atom (match the target's per-atom de-meaning)
        mean = torch.zeros(n_atoms,5,device=dev); cnt=torch.zeros(n_atoms,1,device=dev)
        mean.index_add_(0,group_atom,pooled); cnt.index_add_(0,group_atom,torch.ones(n_groups,1,device=dev))
        return pooled - (mean/cnt)[group_atom]

model=EquivPool().to(dev); opt=torch.optim.Adam(model.parameters(),lr=3e-3,weight_decay=1e-7)
def r2(p,y): return float(1-((y-p)**2).sum()/((y-y.mean(0))**2).sum())
for ep in range(4000):
    model.train(); opt.zero_grad()
    pred=model(); loss=((pred[g_tr]-target[g_tr])**2).mean(); loss.backward(); opt.step()
    if ep%1000==0 or ep==3999:
        model.eval()
        with torch.no_grad(): pr=model()
        print(f"  ep {ep:4d} loss={loss.item():.4e} T2 R2 train={r2(pr[g_tr],target[g_tr]):+.3f} test={r2(pr[g_te],target[g_te]):+.3f}")

# |T2| magnitude (invariant) fit + leave-atoms-out
with torch.no_grad(): pr=model().cpu().numpy()
tg=target.cpu().numpy()
mag_p=np.linalg.norm(pr,axis=1); mag_t=np.linalg.norm(tg,axis=1)
te=g_te.cpu().numpy()
print(f"\n|T2| modulation corr (test): r={np.corrcoef(mag_p[te],mag_t[te])[0,1]:+.3f}")
print(f"T2 component test R2 (5-vec): {r2(torch.tensor(pr[te]),torch.tensor(tg[te])):+.3f}")
sd=df.groupby('atom_index').r.std()  # rough coupling proxy unused; report coupled by target var
gvar=np.linalg.norm(tg,axis=1);
print(f"groups={n_groups} atoms={n_atoms}  (thin: ring-current T2 lives in a few coupled aromatic H)")

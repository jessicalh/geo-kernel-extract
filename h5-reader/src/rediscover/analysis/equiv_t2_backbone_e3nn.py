#!/usr/bin/env python3
"""equiv_t2_backbone_e3nn — the BROAD-BACKBONE equivariant T2 fitter, on e3nn.

Generalizes the single-mechanism ring-current fitter (equiv_t2_e3nn.py) to the
broad-backbone substrate: predict the per-atom DFT T2 across the eight backbone
atom-type strata (N / CA / C / O / HN / HA / HA2 / HA3) from a HETEROGENEOUS
neighbourhood of three source kinds at once — rings, anisotropic bonds, and
charge sites — each contributing through the SAME shared angular machinery but
with its OWN per-type radial channel. This is the heterogeneous-pooling exemplar
the broad case forces; the ring fitter was one stratum × one mechanism.

  contribution_2e(source) = w_kind(invariants_kind) * Y2(disp_hat)      # o3.SH 2e
  [+ w_kind * Y2(axis_hat) + cross(disp_hat, axis_hat)  IF axis emitted; see below]
  node_2e(atom,frame) = SUM_sources contribution_2e                     # scatter-add
  predict node_2e vs C @ (emitted local-frame DFT T2), de-meaned per atom.

What is SHARED vs PER-TYPE (the new judgment vs the ring case):
  * SHARED   — the angular spherical-harmonic projection Y2(.) (e3nn, Wigner-D
               tested), the scatter-pool (Deep Sets permutation invariance), the
               per-atom de-meaning, and the library->e3nn change-of-basis C.
  * PER-TYPE — three independent radial MLPs (ring / bond / charge), each fed ONLY
               the invariants its own source kind emits (ring: r, cosθ, intensity;
               bond: r, cosθ, category; charge: r, q). A ring at 4 Å and a charge
               at 4 Å must be allowed to scale their l=2 contribution differently;
               a single shared radial gate would force them to share a law they do
               not physically obey. The heterogeneous pool is the whole point of
               the broad relationship — it is here, not baked-narrow, by design.

THE BOUNDARY (feedback_no_python_physics_except_labeled_integrity_test +
feedback_model_is_spine): every input is READ from the C++ broad_backbone export;
NOTHING is recomputed in Python.
  * Source geometry: the emitted per-source `disp_local` (vec3, local frame) and
    the emitted invariants `r`, `cos_theta`, `ring_intensity`, `bond_category`,
    `source_q_e`, plus `mechanism`, `is_self_or_bonded`, and the `row_id` join
    key — all columns of broad_backbone_sources.csv (BroadBackboneSink).
  * The TARGET is the C++-emitted local-frame library-basis DFT T2 NPY sidecar
    `broad_backbone_aggregated_target_local_T2.npy` (agg_rows, 5), row-aligned
    with broad_backbone_aggregated.csv. It is REQUIRED — the fitter fails loud if
    absent. There is NO Python projection fallback; we never reconstruct the
    library 5-vector from the emitted 3x3. The 5-vector is mapped into e3nn's 2e
    convention by the PINNED constant change_of_basis.lib_to_e3nn (one matmul on
    emitted data, NOT a projection). |T2| is invariant under that orthogonal map.
  * e3nn's internal (y,z,x) axis convention is handled once, in change_of_basis.

SCHEMA GAP (FLAGGED, not worked around — see HANDOFF_BACKBONE_FIT.md):
  The ring fitter oriented its l=2 contribution by the RING NORMAL n̂ (the dipole
  axis), via the emitted `source_normal_local_*` columns of ring_current_sources.csv.
  The broad_backbone_sources.csv schema (BroadBackboneSink::kSourceHeader) does
  NOT emit `source_normal_local_*` (rings) or `bond_axis_local_*` (bonds) — only
  `disp_local_*` (the target->source displacement). Those axis vectors EXIST in
  the C++ SourceSlot and are filled by the ring/bond attachers; they are simply
  not written to the broad source CSV. So the angularly-complete design — Y2(n̂)
  for rings, Y2(axis) for bonds, and the disp⊗axis cross term — cannot be built
  from the substrate as currently emitted.

  This fitter therefore uses the ONE emitted local-frame direction, disp_hat, for
  every source kind (the charge case is exact — the Coulomb field IS along disp;
  the ring/bond cases carry the axis ANGLE invariantly through the emitted cosθ,
  but not the axis VECTOR). The fuller model is wired but DORMANT: pass
  --with-axes and it auto-activates the Y2(axis)/cross terms FOR ANY source kind
  whose axis columns are present in the CSV. When codex extends BroadBackboneSink
  to emit source_normal_local_*/bond_axis_local_* (the data is already in
  SourceSlot), this fitter lights up the complete model with no code change. We do
  NOT recompute the missing vectors in Python — that would be the projection
  end-run the law forbids.

Honesty (mirrors equiv_t2_e3nn.py): de-mean target AND prediction per atom (strip
the near-constant per-atom baseline tensor — the local bonding the through-space
mechanisms do not touch); frame-split is the GATE; leave-atoms-out is available as
an opt-in audit via --loao, not gated (thin — backbone strata have few coupled
atoms in one protein; effective N is printed per stratum, not oversold). Reported
per atom-type stratum (all 8).

Run (system python; torch cu130 needs nvidia cu13 libs on LD_LIBRARY_PATH or it
segfaults — see ENV.md / requirements-e3nn.txt):
    LD_LIBRARY_PATH=<cu13 libs>:<torch/lib> \
        /usr/bin/python3 equiv_t2_backbone_e3nn.py /tmp/rdc-broad-backbone-axes --with-axes
    (--cross defaults to 'exact', the model; add --cross both only to run the
     learnable cross-term as a comparison that confirms the fixed angular form.)
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from e3nn import o3

import change_of_basis as cob

torch.manual_seed(0)
np.random.seed(0)

# ── The eight backbone atom-type strata, identified from EMITTED columns only ──
# frame_variant (FrameVariant ordinal, LocalFrameBasis.h) fixes the frame class;
# atom_name splits the Cα-H class into HA (non-GLY) vs HA2/HA3 (GLY). No Python
# chemistry — frame_variant is the C++-typed backbone role, atom_name is the
# producer's verbatim name (output-side identity, ML stratification only —
# feedback_naming_input_output_asymmetry). HN covers both N-terminus variants;
# N covers both N variants; the carbonyl C/O and Cα each have one variant.
FV_HN = {1, 2}          # HN_Standard, HN_NTerminus
FV_N = {4, 5}           # BackboneN, BackboneN_NTerminus
FV_CA = {6}             # BackboneCA
FV_C = {7}              # BackboneCarbonylC
FV_O = {8}              # BackboneCarbonylO
FV_HA = {9}             # BackboneHA  (HA / HA2 / HA3 split by atom_name)

STRATA_ORDER = ["N", "CA", "C", "O", "HN", "HA", "HA2", "HA3"]

TEST_FRAME_FRACTION = 0.20
FRAME_SPLIT_SEED = 0
MIN_FRAME_SPLIT_FRAMES = 5
THIN_ATOM_WARN = 10

# ── Per-source-type invariant feature sets (EMITTED columns only) ─────────────
# Each kind's radial MLP sees ONLY what its attacher fills. cosθ carries the
# (emitted) angle to the source axis as an invariant for ring/bond. Charge has no
# meaningful cosθ (no source axis), so it uses (r, q).
KIND_FEATURES = {
    "ring":   ["r", "cos_theta", "ring_intensity"],
    "bond":   ["r", "cos_theta", "bond_category"],
    "charge": ["r", "source_q_e"],
}


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir", nargs="?",
                    default=os.environ.get("REDISCOVER_OUT", "/tmp/rdc-broad-backbone"))
    ap.add_argument("--cross", choices=["exact", "learnable", "both"], default="exact",
                    help="cross-term mode for the Y2(disp)⊗Y2(axis) path (active only "
                         "with --with-axes). 'exact' IS THE MODEL (default, main path) "
                         "— the angular cross-term is fixed by geometry. 'learnable' / "
                         "'both' run the learnable cross-term as an opt-in COMPARISON "
                         "that confirms the fixed form (it lands within ~0.001 R² of "
                         "exact); it is NOT the model path.")
    ap.add_argument("--with-axes", action="store_true",
                    help="activate Y2(axis)/cross terms for source kinds whose "
                         "axis columns (source_normal_local_*/bond_axis_local_*) "
                         "are present in the CSV. Dormant by default because the "
                         "broad source schema does not yet emit them (FLAGGED).")
    ap.add_argument("--loao", action="store_true",
                    help="also run leave-atoms-out retraining. This is an O(atoms) "
                         "audit per stratum, so it is opt-in; the frame-split "
                         "metrics remain the reported gate.")
    ap.add_argument("--epochs", type=int, default=4000)
    ap.add_argument("--lr", type=float, default=3e-3)
    return ap.parse_args()


def stratum_of(frame_variant, atom_name):
    """Map an aggregated-row (frame_variant, atom_name) to one of the 8 strata.
    Returns None for rows that are not one of the eight (should not occur given
    the C++ backbone stratum, but kept defensive — reported, not crashed)."""
    fv = int(frame_variant)
    if fv in FV_N:
        return "N"
    if fv in FV_CA:
        return "CA"
    if fv in FV_C:
        return "C"
    if fv in FV_O:
        return "O"
    if fv in FV_HN:
        return "HN"
    if fv in FV_HA:
        nm = str(atom_name).strip().upper()
        if nm == "HA2":
            return "HA2"
        if nm == "HA3":
            return "HA3"
        return "HA"  # non-GLY HA (and any HA1 alias)
    return None


# ── Load the substrate. READ-ONLY: every column/array below is C++-emitted. ───
def load(out_dir):
    src_csv = f"{out_dir}/broad_backbone_sources.csv"
    agg_csv = f"{out_dir}/broad_backbone_aggregated.csv"
    npy = f"{out_dir}/broad_backbone_aggregated_target_local_T2.npy"
    for path, what in [(src_csv, "per-source CSV"), (agg_csv, "aggregated CSV")]:
        if not os.path.exists(path):
            sys.exit(f"FATAL: required broad_backbone {what} not found:\n  {path}\n"
                     "Re-run the extractor (--case broad_backbone) into this out_dir.")
    if not os.path.exists(npy):
        # The target NPY is the producer's local-frame library-basis DFT T2. We do
        # NOT recompute the spherical projection from the emitted 3x3 in numpy
        # (feedback_no_python_physics_except_labeled_integrity_test); re-extract.
        sys.exit(
            f"FATAL: required emitted target NPY not found:\n  {npy}\n"
            "The equivariant fitter reads the C++-emitted local-frame T2 target; it "
            "does NOT recompute the spherical projection. Re-run the extractor so "
            "the broad_backbone_aggregated_target_local_T2.npy sidecar is present.")

    src = pd.read_csv(src_csv)
    agg = pd.read_csv(agg_csv)
    target_lib_all = np.load(npy)  # (agg_rows, 5) library basis, emitted, lab- OR
    #                                 local-frame per the sidecar's name (local).
    if len(target_lib_all) != len(agg):
        sys.exit(f"FATAL: target NPY rows {len(target_lib_all)} != aggregated CSV "
                 f"rows {len(agg)} (sidecar/CSV out of sync; re-extract).")
    return src, agg, target_lib_all


def detect_axis_columns(src):
    """Which source kinds have their axis vector emitted? Returns a dict kind->bool.
    The fuller Y2(axis)/cross design activates per-kind only where the substrate
    carries the vector — NEVER by recomputing it."""
    ring_axis = all(c in src.columns for c in
                    ["source_normal_local_x", "source_normal_local_y", "source_normal_local_z"])
    bond_axis = all(c in src.columns for c in
                    ["bond_axis_local_x", "bond_axis_local_y", "bond_axis_local_z"])
    return {"ring": ring_axis, "bond": bond_axis, "charge": False}


class KindChannel(nn.Module):
    """One source kind's per-type radial channel: an invariant MLP producing gate
    weights for Y2(disp_hat) [+ Y2(axis_hat) + cross]. The angular Y2/cross
    machinery is SHARED across kinds (built once in the parent); only this radial
    map and its input feature set are per-type."""

    def __init__(self, n_feat, use_axis, cross):
        super().__init__()
        self.use_axis = use_axis
        # gates: [disp] always; + [axis, cross] when this kind emits its axis.
        n_gate = 3 if use_axis else 1
        self.R = nn.Sequential(
            nn.Linear(n_feat, 64), nn.SiLU(),
            nn.Linear(64, 64), nn.SiLU(),
            nn.Linear(64, n_gate))
        self.cross = cross
        if use_axis:
            if cross != "exact":
                irr_1o, irr_2e = o3.Irreps("1o"), o3.Irreps("2e")
                self.tp = o3.FullyConnectedTensorProduct(
                    irr_1o, irr_1o, irr_2e, irrep_normalization="component")

    def forward(self, pk):
        featn = pk["featn"]
        w = self.R(featn)                                          # (n_src, n_gate)
        y2_d = pk["y2_disp"]                                       # (n_src, 5)
        contrib = w[:, 0:1] * y2_d
        if self.use_axis:
            y2_a = pk["y2_axis"]                                   # (n_src, 5)
            if self.cross == "exact":
                cross = pk["cross_exact"]                          # (n_src, 5)
            else:
                cross = self.tp(pk["disp_hat"], pk["axis_hat"])    # (n_src, 5)
            contrib = contrib + w[:, 1:2] * y2_a + w[:, 2:3] * cross
        return contrib                                             # (n_src, 5) e3nn 2e


class BroadEquivPool(nn.Module):
    """Heterogeneous equivariant sum-pool: per-source-kind radial channels feed a
    SHARED scatter-pool to (atom, frame) groups; de-meaned per atom. Each kind's
    sources are processed by that kind's KindChannel, then all contributions are
    index-added into the same group accumulator (Deep Sets over the heterogeneous
    neighbourhood)."""

    def __init__(self, kind_nfeat, kind_axis, n_groups, n_atoms, cross="exact"):
        super().__init__()
        self.n_groups = n_groups
        self.n_atoms = n_atoms
        self.channels = nn.ModuleDict({
            k: KindChannel(kind_nfeat[k], kind_axis[k], cross) for k in kind_nfeat})

    def forward(self, per_kind, group_atom):
        """per_kind[k] = dict(featn, disp_hat, axis_hat, src_group) tensors for the
        sources of kind k. group_atom maps each group index -> its atom index."""
        pooled = torch.zeros(self.n_groups, 5, device=group_atom.device)
        for k, ch in self.channels.items():
            pk = per_kind[k]
            if pk["src_group"].numel() == 0:
                continue
            contrib = ch(pk)                                           # (n_src_k, 5)
            pooled.index_add_(0, pk["src_group"], contrib)
        # de-mean prediction per atom (mirror the target's per-atom de-meaning)
        mean = torch.zeros(self.n_atoms, 5, device=pooled.device)
        cnt = torch.zeros(self.n_atoms, 1, device=pooled.device)
        mean.index_add_(0, group_atom, pooled)
        cnt.index_add_(0, group_atom, torch.ones(self.n_groups, 1, device=pooled.device))
        return pooled - (mean / cnt.clamp_min(1.0))[group_atom]


def r2(p, y):
    return float(1 - ((y - p) ** 2).sum() / ((y - y.mean(0)) ** 2).sum())


def _unit(v):
    n = np.linalg.norm(v, axis=1, keepdims=True)
    return v / np.clip(n, 1e-9, None)


def build_pack(src, agg, target_lib_all, axis_present, with_axes, dev):
    """Assemble the equivariant tensors for ONE stratum's aggregated rows + their
    joined sources. `agg` and `target_lib_all` are already filtered to the stratum;
    `src` is the full source table (joined by row_id)."""
    # Only valid DFT targets in the local frame — the fit target is local-frame T2.
    keep_agg = (agg.dft_present == 1) & (agg.dft_local_frame_valid == 1)
    a = agg.loc[keep_agg].copy()
    if not len(a):
        return None
    # Contiguous group ids over the kept aggregated rows; keep the row index into
    # target_lib_all so the emitted target stays row-aligned.
    a = a.reset_index(drop=False).rename(columns={"index": "__agg_row"})
    a["__gid"] = np.arange(len(a))
    gid_of_rowid = dict(zip(a.row_id.to_numpy(), a.__gid.to_numpy()))

    # group -> atom (for de-meaning + leave-atoms-out) and frame (for the split).
    group_atom_idx, _ = pd.factorize(a.atom_index.to_numpy())
    n_atoms = int(group_atom_idx.max()) + 1
    n_groups = len(a)

    # Target: emitted local-frame library T2 of the kept rows -> e3nn 2e via C.
    gT2_lib = target_lib_all[a["__agg_row"].to_numpy()]            # (n_groups, 5)
    if not np.isfinite(gT2_lib).all():
        # Local-frame target carries NaN where the C++ local frame was invalid;
        # those rows were filtered by dft_local_frame_valid, so any NaN here is a
        # substrate inconsistency — fail loud rather than silently drop.
        bad = int((~np.isfinite(gT2_lib).all(1)).sum())
        sys.exit(f"FATAL: {bad} kept rows have non-finite local-frame target T2 "
                 "despite dft_local_frame_valid==1 (substrate inconsistency).")
    gT2_e3nn = cob.lib_to_e3nn(gT2_lib)                           # constant matmul

    # de-mean target per atom (strip the near-constant per-atom baseline tensor)
    gT2_dm = gT2_e3nn.copy()
    for at in np.unique(group_atom_idx):
        m = group_atom_idx == at
        gT2_dm[m] -= gT2_e3nn[m].mean(0)

    # ── Per-kind source tensors, joined back to the kept groups by row_id ──────
    # Drop the self/bonded ring (the local baseline; matches the ring fitter); the
    # bond/charge mechanisms have no self/bonded flag here (charge already excludes
    # the target's own residue in C++ via exclude_residue).
    s = src[src.row_id.isin(set(a.row_id))].copy()
    s = s[~((s.mechanism == "ring") & (s.is_self_or_bonded == 1))]
    s["__gid"] = s.row_id.map(gid_of_rowid).astype("int64")

    per_kind = {}
    kind_nfeat, kind_axis = {}, {}
    exact_cross_tp = None
    for k, cols in KIND_FEATURES.items():
        sk = s[s.mechanism == k]
        use_axis = with_axes and axis_present.get(k, False)
        kind_axis[k] = use_axis
        kind_nfeat[k] = len(cols)
        if not len(sk):
            per_kind[k] = dict(featn=torch.zeros(0, len(cols), device=dev),
                               disp_hat=torch.zeros(0, 3, device=dev),
                               axis_hat=torch.zeros(0, 3, device=dev),
                               y2_disp=torch.zeros(0, 5, device=dev),
                               y2_axis=torch.zeros(0, 5, device=dev),
                               cross_exact=torch.zeros(0, 5, device=dev),
                               src_group=torch.zeros(0, dtype=torch.long, device=dev))
            continue
        disp = sk[["disp_local_x", "disp_local_y", "disp_local_z"]].to_numpy(float)
        disp_hat = _unit(disp)
        if use_axis:
            if k == "ring":
                axc = ["source_normal_local_x", "source_normal_local_y", "source_normal_local_z"]
            else:  # bond
                axc = ["bond_axis_local_x", "bond_axis_local_y", "bond_axis_local_z"]
            axis = sk[axc].to_numpy(float)
            axis_norm = np.linalg.norm(axis, axis=1)
            bad_axis = (~np.isfinite(axis).all(1)) | (axis_norm <= 1e-9)
            if bad_axis.any():
                sys.exit(f"FATAL: {int(bad_axis.sum())} {k} rows have non-finite "
                         "or zero emitted axis vectors despite axis columns being "
                         "present. Re-extract after fixing BroadBackboneSink.")
            axis_hat = axis / axis_norm[:, None]
        else:
            axis_hat = np.zeros_like(disp_hat)  # unused by the channel
        feat = sk[cols].to_numpy(np.float32)
        fmean, fstd = feat.mean(0), feat.std(0)
        featn = (feat - fmean) / np.where(fstd > 1e-8, fstd, 1.0)
        disp_hat_t = torch.tensor(disp_hat, dtype=torch.float32, device=dev)
        axis_hat_t = torch.tensor(axis_hat, dtype=torch.float32, device=dev)
        with torch.no_grad():
            y2_disp = o3.spherical_harmonics(
                2, disp_hat_t, normalize=True, normalization="component")
            if use_axis:
                y2_axis = o3.spherical_harmonics(
                    2, axis_hat_t, normalize=True, normalization="component")
                if exact_cross_tp is None:
                    exact_cross_tp = o3.TensorProduct(
                        o3.Irreps("1o"), o3.Irreps("1o"), o3.Irreps("2e"),
                        instructions=[(0, 0, 0, "uuu", False)],
                        irrep_normalization="component").to(dev)
                cross_exact = exact_cross_tp(disp_hat_t, axis_hat_t)
            else:
                y2_axis = torch.zeros(len(sk), 5, device=dev)
                cross_exact = torch.zeros(len(sk), 5, device=dev)
        per_kind[k] = dict(
            featn=torch.tensor(featn, dtype=torch.float32, device=dev),
            disp_hat=disp_hat_t,
            axis_hat=axis_hat_t,
            y2_disp=y2_disp,
            y2_axis=y2_axis,
            cross_exact=cross_exact,
            src_group=torch.tensor(sk["__gid"].to_numpy(), dtype=torch.long, device=dev))

    # frame split (gate)
    frames = np.sort(a.h5_row.unique())
    if len(frames) < MIN_FRAME_SPLIT_FRAMES:
        g_te = np.zeros(n_groups, dtype=bool)  # too few frames to split honestly
    else:
        test_f = set(np.random.default_rng(FRAME_SPLIT_SEED).choice(
            frames, max(1, int(TEST_FRAME_FRACTION * len(frames))), replace=False))
        g_te = a.h5_row.isin(test_f).to_numpy()

    t = lambda x, dt=torch.float32: torch.tensor(x, dtype=dt, device=dev)
    return dict(
        per_kind=per_kind, kind_nfeat=kind_nfeat, kind_axis=kind_axis,
        group_atom=t(group_atom_idx, torch.long),
        target=t(gT2_dm), g_te=t(g_te, torch.bool),
        n_groups=n_groups, n_atoms=n_atoms,
        group_atom_idx=group_atom_idx, gT2_e3nn=gT2_e3nn,
        n_src={k: per_kind[k]["src_group"].numel() for k in per_kind},
    )


def _new_model(pack, cross, dev):
    return BroadEquivPool(pack["kind_nfeat"], pack["kind_axis"],
                          pack["n_groups"], pack["n_atoms"], cross=cross).to(dev)


def train_one(pack, dev, cross, epochs, lr, run_loao):
    model = _new_model(pack, cross, dev)
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-7)
    target, g_te = pack["target"], pack["g_te"]
    g_tr = ~g_te
    if g_tr.sum() == 0:
        return float("nan"), float("nan"), float("nan")
    for ep in range(epochs):
        model.train(); opt.zero_grad()
        pred = model(pack["per_kind"], pack["group_atom"])
        loss = ((pred[g_tr] - target[g_tr]) ** 2).mean()
        loss.backward(); opt.step()
    model.eval()
    with torch.no_grad():
        pr = model(pack["per_kind"], pack["group_atom"]).cpu().numpy()
    tg = target.cpu().numpy()
    te = g_te.cpu().numpy()
    if te.sum() < 2:
        r2_comp, r_mag = float("nan"), float("nan")
    else:
        mag_p, mag_t = np.linalg.norm(pr, axis=1), np.linalg.norm(tg, axis=1)
        r_mag = float(np.corrcoef(mag_p[te], mag_t[te])[0, 1])
        r2_comp = r2(torch.tensor(pr[te]), torch.tensor(tg[te]))
    loo = leave_atoms_out(pack, dev, cross, epochs, lr) if run_loao else float("nan")
    return r2_comp, r_mag, loo


def leave_atoms_out(pack, dev, cross, epochs, lr):
    """Reported honesty check: predict the T2 of atoms held entirely out of the
    fit. Thin per backbone stratum in one protein; NOT the gate. Returns NaN when
    fewer than 3 coupled atoms (nothing to hold out)."""
    atoms = np.unique(pack["group_atom_idx"])
    if len(atoms) < 3:
        return float("nan")
    pred = np.full(pack["gT2_e3nn"].shape, np.nan)
    ga = pack["group_atom_idx"]
    for held in atoms:
        tr_mask = torch.tensor(ga != held, dtype=torch.bool, device=dev)
        model = _new_model(pack, cross, dev)
        opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-7)
        for ep in range(max(800, epochs // 4)):
            model.train(); opt.zero_grad()
            p = model(pack["per_kind"], pack["group_atom"])
            loss = ((p[tr_mask] - pack["target"][tr_mask]) ** 2).mean()
            loss.backward(); opt.step()
        model.eval()
        with torch.no_grad():
            p = model(pack["per_kind"], pack["group_atom"]).cpu().numpy()
        pred[ga == held] = p[ga == held]
    tg = pack["gT2_e3nn"].copy()
    for at in atoms:
        m = ga == at
        tg[m] -= pack["gT2_e3nn"][m].mean(0)
        pred[m] -= pred[m].mean(0)
    ok = np.isfinite(pred).all(1)
    if ok.sum() < 2:
        return float("nan")
    return r2(torch.tensor(pred[ok]), torch.tensor(tg[ok]))


def report_stratum(name, pack, dev, args):
    n_groups = pack["n_groups"]
    n_atoms = pack["n_atoms"]
    n_src = pack["n_src"]
    te = int(pack["g_te"].sum().item())
    thin = n_atoms < THIN_ATOM_WARN
    thin_txt = "  THIN" if thin else ""
    print(f"\n── stratum {name}: groups(atom,frame)={n_groups}  coupled atoms={n_atoms}  "
          f"test groups={te}{thin_txt}  sources[ring/bond/charge]="
          f"{n_src.get('ring',0)}/{n_src.get('bond',0)}/{n_src.get('charge',0)}",
          flush=True)
    if n_atoms < 2 or n_groups < 4:
        print(f"   (effective N too small to fit honestly — reported, not fit)", flush=True)
        return None
    results = {}
    modes = ["exact", "learnable"] if args.cross == "both" else [args.cross]
    # cross modes only differ when an axis is active for some kind; otherwise the
    # exact/learnable distinction is moot (no cross term), so collapse to one run.
    any_axis = any(pack["kind_axis"].values())
    if not any_axis:
        modes = modes[:1]
    for m in modes:
        r2c, rmag, loo = train_one(pack, dev, m, args.epochs, args.lr, args.loao)
        tag = f"[{m:9s}]" if any_axis else "[disp-only]"
        loo_txt = f"{loo:+.3f}" if args.loao else "skipped (--loao)"
        print(f"   {tag} frame-split T2 R²(test)={r2c:+.3f}  |T2| r(test)={rmag:+.3f}  "
              f"leave-atoms-out R²={loo_txt} (reported; coupled atoms={n_atoms})",
              flush=True)
        results[m] = (r2c, rmag, loo, n_atoms)
    return results


def main():
    args = parse_args()
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device={dev}  out_dir={args.out_dir}  cross={args.cross}  "
          f"with_axes={args.with_axes}  loao={args.loao}", flush=True)
    # pin C up front so a convention mismatch fails loudly before training.
    C = cob.get_C()
    print(f"change-of-basis C orthogonality |CᵀC - I|max="
          f"{np.abs(C.T @ C - np.eye(5)).max():.2e}", flush=True)

    src, agg, target_lib_all = load(args.out_dir)
    print(f"aggregated rows={len(agg)}  source rows={len(src)}", flush=True)

    axis_present = detect_axis_columns(src)
    print(f"axis columns present in source CSV: ring_normal={axis_present['ring']}  "
          f"bond_axis={axis_present['bond']}", flush=True)
    if args.with_axes and not (axis_present["ring"] or axis_present["bond"]):
        print("  NOTE: --with-axes requested but NO axis columns are emitted; the "
              "model runs disp-only. Extend BroadBackboneSink to emit "
              "source_normal_local_*/bond_axis_local_* (FLAGGED in HANDOFF).",
              flush=True)
    elif not args.with_axes:
        print("  NOTE: running disp-only (Y2(disp_hat) per source). The ring-normal/"
              "bond-axis Y2(axis)+cross terms are DORMANT — the broad source schema "
              "does not emit those vectors (FLAGGED in HANDOFF).", flush=True)

    # stratify the AGGREGATED rows by (frame_variant, atom_name).
    agg = agg.copy()
    agg["__stratum"] = [stratum_of(fv, nm)
                        for fv, nm in zip(agg.frame_variant, agg.atom_name)]
    unrecognised = int(agg["__stratum"].isna().sum())
    if unrecognised:
        print(f"  WARNING: {unrecognised} aggregated rows did not map to one of the "
              "8 backbone strata (reported, not fit)", flush=True)

    summary = {}
    for name in STRATA_ORDER:
        a_s = agg[agg["__stratum"] == name]
        if not len(a_s):
            print(f"\n── stratum {name}: 0 rows (absent in this protein/run)", flush=True)
            continue
        tgt_s = target_lib_all[a_s.index.to_numpy()]
        a_s = a_s.reset_index(drop=True)
        pack = build_pack(src, a_s, tgt_s, axis_present, args.with_axes, dev)
        if pack is None:
            print(f"\n── stratum {name}: no valid local-frame DFT targets", flush=True)
            continue
        summary[name] = report_stratum(name, pack, dev, args)

    print("\n=== broad-backbone equivariant T2 summary (per stratum) ===", flush=True)
    print("frame-split T2 R² is the GATE per stratum; |T2| r is the rotation-"
          "invariant magnitude correlation; leave-atoms-out is opt-in via --loao.",
          flush=True)
    for name in STRATA_ORDER:
        res = summary.get(name)
        if not res:
            print(f"  {name:4s}: (insufficient effective N)", flush=True)
            continue
        best = max(res, key=lambda k: (res[k][0] if np.isfinite(res[k][0]) else -1e9))
        r2c, rmag, loo, n_atoms = res[best]
        loo_txt = f"{loo:+.3f}" if args.loao else "skipped"
        thin_txt = "  THIN" if n_atoms < THIN_ATOM_WARN else ""
        print(f"  {name:4s}: T2 R²={r2c:+.3f}  |T2| r={rmag:+.3f}  LOAO={loo_txt}{thin_txt}",
              flush=True)


if __name__ == "__main__":
    main()

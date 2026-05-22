#!/usr/bin/env python3
"""
test_inspector.py — integration test for the viewer's Atom Inspector
via the REST API.

The viewer's atom_dump + list_atoms commands return the same typed
per-atom data the Inspector tree displays, as machine-readable JSON.
This script uses them to verify that every calculator family the
viewer surfaces is reading non-default values across enough atoms to
be confident the wiring is correct.

Discipline (per feedback_log_overages_dont_assert):
  - This script LOGS findings; it does not assert hard pass/fail.
  - Pipeline produces imperfect output for legitimate reasons
    (chi-fallback misses, perception edges, Tripeptide DB coverage,
    Larsen grid skips). We report counts and let the operator judge.
  - Only structural invariants (e.g. backbone N must have
    graph_dist_N == 0) are reported as ERROR.

Usage:
  # Launch the viewer first:
  scripts/run_with_cuda_env.sh build/nmr-viewer --pdb \\
      tests/data/external/1UBQ.pdb &

  # Then run this script (it talks to localhost:9147):
  python3 ui/tests/test_inspector.py

  # Optional flags:
  python3 ui/tests/test_inspector.py --host 127.0.0.1 --port 9147
  python3 ui/tests/test_inspector.py --verbose         # dump all atom records
"""

from __future__ import annotations

import argparse
import json
import socket
import sys
from typing import Any


# ─── REST client ────────────────────────────────────────────────────

def send_cmd(host: str, port: int, cmd: dict[str, Any]) -> dict[str, Any]:
    """Send one JSON command, return the parsed response. Raises on bad reply.

    The viewer's RestServer (ui/src/RestServer.cpp) writes one JSON
    object followed by '\\n' per command and keeps the connection open
    for more commands — it does NOT close after replying. We read
    line-buffered until we get the first newline-terminated reply.
    Using a per-call short-lived socket is the simplest model; the
    test does ~1300 small calls and localhost TCP handles it fine.

    Scaling note: this opens one TCP connection per command. Fine up to
    ~5000 atoms; above that, TIME_WAIT socket pile-up starts mattering
    (Linux default net.ipv4.tcp_max_tw_buckets = 262144 is plenty for
    one run, but back-to-back runs across many large proteins can
    accumulate). Switch to a persistent socket + delimiter-driven
    framing when that boundary actually matters."""
    with socket.create_connection((host, port), timeout=30.0) as s:
        s.sendall((json.dumps(cmd) + "\n").encode("utf-8"))
        rf = s.makefile("r", encoding="utf-8")
        line = rf.readline()
    raw = line.strip()
    if not raw:
        raise RuntimeError(f"empty reply to {cmd}")
    obj = json.loads(raw)
    if not obj.get("ok"):
        raise RuntimeError(f"REST error on {cmd}: {obj.get('error', '?')}")
    return obj["result"]


# ─── Reporter (PASS / NOTE / WARN / ERROR), no asserts ──────────────

class Reporter:
    def __init__(self) -> None:
        self.counts = {"PASS": 0, "NOTE": 0, "WARN": 0, "ERROR": 0}
        self.first_failures: list[str] = []

    def _record(self, level: str, msg: str) -> None:
        self.counts[level] = self.counts.get(level, 0) + 1
        prefix = {
            "PASS":  "\033[32mPASS\033[0m ",
            "NOTE":  "\033[36mNOTE\033[0m ",
            "WARN":  "\033[33mWARN\033[0m ",
            "ERROR": "\033[31mERROR\033[0m",
        }.get(level, level)
        print(f"  {prefix}  {msg}")
        if level in ("WARN", "ERROR") and len(self.first_failures) < 20:
            self.first_failures.append(f"[{level}] {msg}")

    def passes(self, msg: str) -> None:  self._record("PASS",  msg)
    def note(self,   msg: str) -> None:  self._record("NOTE",  msg)
    def warns(self,  msg: str) -> None:  self._record("WARN",  msg)
    def errors(self, msg: str) -> None:  self._record("ERROR", msg)

    def summary(self) -> None:
        print()
        print("═" * 70)
        print(f"Summary: {self.counts['PASS']} pass · "
              f"{self.counts['NOTE']} note · "
              f"{self.counts['WARN']} warn · "
              f"{self.counts['ERROR']} error")
        if self.first_failures:
            print("First failures:")
            for line in self.first_failures[:10]:
                print(f"  {line}")


# ─── Per-section checks ─────────────────────────────────────────────

def section(title: str) -> None:
    print()
    print(f"── {title} {'─' * (66 - len(title))}")


def check_status_and_listing(host: str, port: int, r: Reporter) -> tuple[dict, list[dict]]:
    section("status + list_atoms")
    status = send_cmd(host, port, {"cmd": "status"})
    print(f"  protein={status.get('protein', '?')}  "
          f"n_atoms={status.get('n_atoms')}  "
          f"n_residues={status.get('n_residues')}  "
          f"n_rings={status.get('n_rings')}  "
          f"computing={status.get('computing')}")
    if status.get("computing", True):
        r.errors("viewer still computing — wait longer and re-run")
        sys.exit(1)
    if status.get("n_atoms", 0) == 0:
        r.errors("no protein loaded in the viewer")
        sys.exit(1)
    r.passes(f"status responded; protein has {status['n_atoms']} atoms")

    listing = send_cmd(host, port, {"cmd": "list_atoms"})
    atoms = listing["atoms"]
    if len(atoms) == status["n_atoms"]:
        r.passes(f"list_atoms returned all {len(atoms)} atoms")
    else:
        r.errors(f"list_atoms returned {len(atoms)}, "
                 f"expected {status['n_atoms']}")
    return status, atoms


def fetch_all_dumps(host: str, port: int, atoms: list[dict],
                     verbose: bool = False) -> list[dict]:
    """One pass: atom_dump for every atom. All later checks read from
    this in-memory list rather than re-hitting the REST endpoint."""
    section(f"fetching atom_dump × {len(atoms)} (one pass)")
    dumps: list[dict] = []
    for i, a in enumerate(atoms):
        d = send_cmd(host, port, {"cmd": "atom_dump", "atom": a["index"]})
        dumps.append(d)
        if verbose and (i % 100 == 0):
            print(f"    fetched {i}/{len(atoms)}")
    print(f"    done — {len(dumps)} atom dumps cached")
    return dumps


def check_graph_topology(r: Reporter, dumps: list[dict]) -> None:
    section("graph topology invariants (MolecularGraphResult)")
    bn = [d for d in dumps if d["identity"]["role"] == "BackboneN"]
    bo = [d for d in dumps if d["identity"]["role"] == "BackboneO"]
    arc = [d for d in dumps if d["identity"]["role"] == "AromaticC"]

    bad_N = sum(1 for d in bn if d["graph_topology"]["graph_dist_N"] != 0)
    if bad_N == 0 and bn:
        r.passes(f"every BackboneN ({len(bn)}) has graph_dist_N == 0")
    elif bn:
        r.errors(f"{bad_N}/{len(bn)} BackboneN atoms have graph_dist_N != 0")

    bad_O = sum(1 for d in bo if d["graph_topology"]["graph_dist_O"] != 0)
    if bad_O == 0 and bo:
        r.passes(f"every BackboneO ({len(bo)}) has graph_dist_O == 0")
    elif bo:
        r.errors(f"{bad_O}/{len(bo)} BackboneO atoms have graph_dist_O != 0")

    bad_ring = sum(1 for d in arc if d["graph_topology"]["graph_dist_ring"] != 0)
    if arc and bad_ring == 0:
        r.passes(f"every AromaticC ({len(arc)}) has graph_dist_ring == 0")
    elif arc:
        r.errors(f"{bad_ring}/{len(arc)} AromaticC atoms have "
                 f"graph_dist_ring != 0")
    else:
        r.note("no AromaticC atoms in this protein")

    # Conjugated flag should fire on aromatic ring atoms and peptide-bond
    # atoms (C=O, C-N). At minimum, every aromatic C should be flagged.
    bad_conj = sum(1 for d in arc if not d["graph_topology"]["is_conjugated"])
    if arc and bad_conj == 0:
        r.passes(f"every AromaticC has is_conjugated == true")
    elif arc:
        r.warns(f"{bad_conj}/{len(arc)} AromaticC atoms have "
                f"is_conjugated == false (possible perception gap)")


def check_amber_substrate(r: Reporter, dumps: list[dict]) -> None:
    section("AMBER substrate (LegacyAmberTopology)")
    missing = sum(1 for d in dumps if "amber_substrate" not in d)
    if missing == 0:
        r.passes(f"AMBER substrate populated on all {len(dumps)} atoms")
    else:
        r.errors(f"{missing}/{len(dumps)} atoms lack amber_substrate "
                 f"(HasAtomSemantic gate)")

    # Aromatic ring atoms should have in_any_ring == true.
    arc = [d for d in dumps if d["identity"]["role"] == "AromaticC"]
    bad = sum(1 for d in arc
              if not d.get("amber_substrate", {})
                       .get("ring_position", {})
                       .get("in_any_ring", False))
    if arc and bad == 0:
        r.passes(f"every AromaticC has substrate ring_position.in_any_ring "
                 f"== true ({len(arc)} atoms)")
    elif arc:
        r.errors(f"{bad}/{len(arc)} AromaticC atoms have in_any_ring == false")

    # Inventory of planar_group counts (eye check).
    pg_counts: dict[str, int] = {}
    for d in dumps:
        pg = d.get("amber_substrate", {}).get("planar_group", "—")
        pg_counts[pg] = pg_counts.get(pg, 0) + 1
    r.note(f"planar_group distribution: " +
           ", ".join(f"{k}={v}" for k, v in sorted(pg_counts.items())))

    # PolarH distribution (only H atoms are polar)
    polar_counts: dict[str, int] = {}
    for d in dumps:
        ph = d.get("amber_substrate", {}).get("polar_h", "—")
        if ph != "NotPolar" and ph != "—":
            polar_counts[ph] = polar_counts.get(ph, 0) + 1
    if polar_counts:
        r.note(f"polar_h distribution (H atoms only): " +
               ", ".join(f"{k}={v}" for k, v in sorted(polar_counts.items())))


def check_protonation_variants(r: Reporter, dumps: list[dict]) -> None:
    section("Residue protonation variants")
    # One variant per (residue_index) — take first atom seen.
    per_residue: dict[int, tuple[str, str]] = {}  # ridx → (rtype, vname)
    for d in dumps:
        ident = d["identity"]
        ridx = ident["residue_index"]
        if ridx in per_residue:
            continue
        rt = ident["residue_type"]
        vn = ident.get("protonation_variant_name", "—")
        per_residue[ridx] = (rt, vn)

    his = [(rt, vn) for (rt, vn) in per_residue.values() if rt == "HIS"]
    if his:
        r.passes(f"HIS protonation variants: "
                 f"{', '.join(sorted(set(vn for (_, vn) in his)))} "
                 f"({len(his)} HIS residue(s))")
    else:
        r.note("no HIS residues in this protein")

    cys = [(rt, vn) for (rt, vn) in per_residue.values() if rt == "CYS"]
    if cys:
        r.note(f"CYS variants: "
               f"{', '.join(sorted(set(vn for (_, vn) in cys)))} "
               f"({len(cys)} CYS residue(s))")

    all_variants = sorted({vn for (_, vn) in per_residue.values() if vn != "—"})
    if all_variants:
        r.passes(f"variant names resolved: {', '.join(all_variants)}")


def check_aimnet2_ranges(r: Reporter, dumps: list[dict]) -> None:
    section("AIMNet2 charge + EFG plausibility")
    extreme = [(d["index"], d["identity"]["pdb_atom_name"],
                d["charges"]["aimnet2_hirshfeld"])
               for d in dumps
               if abs(d["charges"]["aimnet2_hirshfeld"]) > 1.2]
    if not extreme:
        r.passes(f"AIMNet2 Hirshfeld charges all |q| ≤ 1.2 across "
                 f"{len(dumps)} atoms")
    else:
        r.warns(f"{len(extreme)} atoms with |q_aimnet2| > 1.2; "
                f"first: {extreme[0]}")

    zero_emb = sum(1 for d in dumps
                   if d["aimnet2_embedding"]["l2_squared"] <= 0.0)
    if zero_emb == 0:
        r.passes(f"AIMNet2 embedding L2² > 0 on all {len(dumps)} atoms")
    else:
        r.errors(f"{zero_emb}/{len(dumps)} atoms have aimnet2 "
                 f"embedding L2² == 0")

    zero_grad = sum(1 for d in dumps
                    if d["aimnet2_polarisability"]["gradient_magnitude"] == 0.0)
    if zero_grad < len(dumps) * 0.05:
        r.passes(f"AIMNet2 polarisability non-zero on "
                 f"{len(dumps) - zero_grad}/{len(dumps)} atoms")
    else:
        r.warns(f"AIMNet2 polarisability == 0 on "
                f"{zero_grad}/{len(dumps)} atoms")


def check_tripeptide_coverage(r: Reporter, dumps: list[dict]) -> None:
    section("Tripeptide DFT coverage (ProCS15)")
    n_bb = sum(1 for d in dumps if d.get("tripeptide_dft", {})
                                    .get("backbone", {})
                                    .get("has_match"))
    n_nb = sum(1 for d in dumps if d.get("tripeptide_dft", {})
                                    .get("neighbor", {})
                                    .get("has_match"))
    methods: dict[str, int] = {}
    for d in dumps:
        bb = d.get("tripeptide_dft", {}).get("backbone", {})
        if bb.get("has_match"):
            m = bb.get("method_name", "—")
            methods[m] = methods.get(m, 0) + 1
    r.note(f"tripeptide BB: {n_bb}/{len(dumps)} atoms matched")
    r.note(f"tripeptide neighbor: {n_nb}/{len(dumps)} atoms matched")
    if methods:
        r.note("BB method distribution: " +
               ", ".join(f"{k}={v}" for k, v in methods.items()))
    if n_bb == 0:
        r.warns("zero tripeptide BB matches — DSN may not be configured")
    else:
        r.passes(f"tripeptide BB DFT readback active "
                 f"({n_bb} atoms with σ_BB)")


def check_larsen_hbond_coverage(r: Reporter, dumps: list[dict]) -> None:
    section("Larsen H-bond coverage (grid lookup)")
    n_pairs = sum(1 for d in dumps
                  if d.get("larsen_hbond", {}).get("n_pairs", 0) > 0)
    pair_total = sum(d.get("larsen_hbond", {}).get("n_pairs", 0)
                     for d in dumps)
    n_water = sum(1 for d in dumps
                  if abs(d.get("larsen_hbond", {})
                          .get("water_term_ppm", 0.0)) > 1e-9)
    n_imputed = sum(1 for d in dumps
                    if d.get("larsen_hbond", {})
                        .get("any_corner_imputed", False))

    r.note(f"Larsen H-bond: {n_pairs} atoms with ≥1 pair, "
           f"{pair_total} total pair contributions")
    r.note(f"Larsen water term: {n_water} amide H atoms (Δσ_w)")
    if n_imputed:
        r.note(f"Larsen corner imputed: {n_imputed} atoms had at "
               f"least one nearest-neighbour-filled grid bin")
    if n_pairs == 0:
        r.warns("zero Larsen H-bond pairs — grid dir may be unconfigured")
    else:
        r.passes(f"Larsen H-bond grid lookup active "
                 f"({n_pairs} atoms with contribution)")

    # Cross-check: water term atoms should be amide-H per Larsen Table 2
    # (only HN atoms get Δσ_w, and only when they have no H-bond pair).
    water_role_violation = 0
    for d in dumps:
        if abs(d.get("larsen_hbond", {})
                .get("water_term_ppm", 0.0)) > 1e-9:
            if d["identity"]["role"] != "AmideH":
                water_role_violation += 1
    if n_water > 0 and water_role_violation == 0:
        r.passes(f"all {n_water} Larsen water-term atoms are AmideH role")
    elif n_water > 0:
        r.errors(f"{water_role_violation}/{n_water} water-term atoms "
                 f"are not AmideH role (Larsen Table 2 says HN only)")


def check_bonds(r: Reporter, dumps: list[dict]) -> None:
    section("Bonds (covalent / mcconnell / mopac)")
    missing = sum(1 for d in dumps if "bonds" not in d)
    if missing:
        r.errors(f"{missing}/{len(dumps)} atom dumps lack 'bonds' block")
        return
    r.passes(f"every atom has bonds block ({len(dumps)} atoms)")

    # Every atom should have at least one covalent bond (except isolated
    # ions/waters — not present in 1UBQ).
    no_covalent = sum(1 for d in dumps if not d["bonds"]["covalent"])
    if no_covalent == 0:
        r.passes(f"every atom has ≥1 covalent bond")
    else:
        r.warns(f"{no_covalent}/{len(dumps)} atoms have no covalent bonds")

    # Covalent bond entries must have matching category between this
    # atom's view and the topology — same field exposed twice in
    # mcconnell entries for the cross-check.
    mismatches = 0
    for d in dumps:
        for b in d["bonds"]["mcconnell"]:
            if b["category"] != b["topology_category"]:
                mismatches += 1
    if mismatches == 0:
        r.passes("mcconnell.bond_neighbours category agrees with topology")
    else:
        r.errors(f"{mismatches} mcconnell entries with category != topology_category")

    # Total covalent bonds should be (sum over atoms of bond_indices.size())
    # = 2 × BondCount (each bond touches 2 atoms). Library reports
    # BondCount in status.
    total_covalent_endpoints = sum(len(d["bonds"]["covalent"]) for d in dumps)
    r.note(f"total covalent endpoints across atoms: "
           f"{total_covalent_endpoints} (should be 2 × n_bonds)")

    # MOPAC bonds — empty when skip_mopac=true. Just report.
    n_with_mopac = sum(1 for d in dumps if d["bonds"]["mopac"])
    if n_with_mopac == 0:
        r.note("mopac bond orders: empty everywhere (skip_mopac=true expected)")
    else:
        r.note(f"mopac bond orders present on {n_with_mopac} atoms")


def check_shielding_signal(r: Reporter, dumps: list[dict]) -> None:
    section("Per-calculator shielding signal")
    # Per-kernel classification:
    #   "T0"        — kernel carries a trace (T0 ≠ 0 expected when active)
    #   "T2only"    — structurally traceless tensor (PiQuad EFG,
    #                 AIMNet2 EFG, Coulomb EFG, Dispersion 3d⊗d/r⁸-I/r⁶):
    #                 T0 ≡ 0 by construction; |T2| is the signal.
    #   "skipped"   — intentionally not run in this load path
    #                 (viewer sets skip_coulomb=true; APBS replaces vacuum
    #                  Coulomb EFG).
    kernel_kind = {
        "bs":                  "T0",
        "hm":                  "T0",
        "mc":                  "T0",
        "ringchi":             "T0",
        "hbond_kernel":        "T0",
        "hbond_larsen":        "T0",
        "tripeptide_bb":       "T0",
        "tripeptide_neighbor": "T0",
        "mopac_coulomb":       "skipped",   # ComputeWorker sets skip_mopac=true
        "mopac_mc":            "skipped",   # ditto
        "coulomb":             "skipped",   # ComputeWorker sets skip_coulomb=true
        "dispersion":          "T2only",    # K_ab = S(r)(3d_a d_b/r⁸ - δ/r⁶), traceless
        "piquad":              "T2only",    # EFG-class, traceless
        "aimnet2_efg":         "T2only",    # EFG-class, traceless
    }

    def t2_norm2(st: dict[str, Any]) -> float:
        t2 = st.get("T2", [0.0]*5)
        return sum(c*c for c in t2)

    for k, kind in kernel_kind.items():
        sc_field = [d["shielding_contributions"].get(k, {}) for d in dumps]
        if kind == "skipped":
            # Should be all zeros. If anything non-zero, that's surprising.
            nonzero_T0 = sum(1 for st in sc_field
                             if abs(st.get("T0", 0.0)) > 1e-12)
            nonzero_T2 = sum(1 for st in sc_field if t2_norm2(st) > 1e-24)
            if nonzero_T0 == 0 and nonzero_T2 == 0:
                r.note(f"{k}: skipped (skip_* flag in ComputeWorker)")
            else:
                r.warns(f"{k}: marked skipped but found "
                        f"T0≠0 on {nonzero_T0}, |T2|>0 on {nonzero_T2} atoms")
        elif kind == "T2only":
            # T0 must be 0; |T2| should be non-zero for atoms in range.
            n_T2 = sum(1 for st in sc_field if t2_norm2(st) > 1e-24)
            n_T0 = sum(1 for st in sc_field
                       if abs(st.get("T0", 0.0)) > 1e-12)
            if n_T0 > 0:
                r.warns(f"{k}: structurally-traceless kernel but T0≠0 on "
                        f"{n_T0} atoms — possible decomposition bug")
            if n_T2 == 0:
                r.warns(f"{k}: |T2| is zero on every atom — calculator "
                        f"may not be writing aggregate")
            else:
                r.passes(f"{k}: {n_T2}/{len(dumps)} atoms with |T2| > 0 "
                         f"(traceless kernel; T0 ≡ 0 as expected)")
        else:  # T0
            n_T0 = sum(1 for st in sc_field
                       if abs(st.get("T0", 0.0)) > 1e-12)
            if n_T0 == 0:
                r.warns(f"{k}: zero atoms with non-zero T0")
            else:
                r.passes(f"{k}: {n_T0}/{len(dumps)} atoms with non-zero T0")


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--host", default="127.0.0.1")
    p.add_argument("--port", type=int, default=9147)
    p.add_argument("--verbose", action="store_true",
                   help="dump full atom record per atom checked")
    args = p.parse_args()

    print(f"Connecting to {args.host}:{args.port} …")
    try:
        send_cmd(args.host, args.port, {"cmd": "status"})
    except Exception as e:
        print(f"ERROR: cannot reach viewer REST: {e}", file=sys.stderr)
        print("  Launch the viewer first:", file=sys.stderr)
        print("    scripts/run_with_cuda_env.sh build/nmr-viewer --pdb "
              "tests/data/external/1UBQ.pdb", file=sys.stderr)
        return 2

    r = Reporter()
    status, atoms = check_status_and_listing(args.host, args.port, r)
    dumps = fetch_all_dumps(args.host, args.port, atoms, args.verbose)
    check_graph_topology(r, dumps)
    check_amber_substrate(r, dumps)
    check_protonation_variants(r, dumps)
    check_aimnet2_ranges(r, dumps)
    check_tripeptide_coverage(r, dumps)
    check_larsen_hbond_coverage(r, dumps)
    check_bonds(r, dumps)
    check_shielding_signal(r, dumps)
    r.summary()
    # Always exit 0 — per feedback_log_overages_dont_assert, this is a
    # report tool, not a strict pass/fail gate. Errors are visible in
    # the summary; the operator decides.
    return 0


if __name__ == "__main__":
    sys.exit(main())

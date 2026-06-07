# McConnell rhombic component — spec + scoping

*Spec cadence sibling of `mcconnell_spec.md` / `ring_spec.md` / `efg_spec.md`. Scopes adding the
**rhombic** (second-component) bond susceptibility anisotropy to the **one** McConnell calculator, for
the two non-cylindrical bond types where the physics demands it: the **carbonyl C=O** and the **C=C
double bond**. Additive by construction — **drops no existing output**. Status: scoped, value-gated on
acquisition. 2026-06-07.*

---

## 1. What, and why it is the one honest gap in McConnell's tensor

A bond with cylindrical symmetry — C≡C, C≡N, C–C, S–S — has **one** magnetic susceptibility anisotropy
`Δχ` along the bond axis. Its source shape is uniaxial, and that is exactly what the producer builds
today:

```
source_shape = (û ⊗ û − I/3)          // McConnellResult.cpp:230 — strictly AXIAL
```

A **flat** group — the carbonyl **C=O** and the **C=C** double bond — is not rod-symmetric: `χ` has
**three inequivalent principal values** (an axial part `Δχ_∥` *and* a rhombic part `Δχ_⊥`), because the
deshielding sits at the bond *ends* with shielding above the plane (`bond_anisotropy.md:73–85, 101–114`,
Abraham Part 13 / ApSimon / RSC Adv. 2014). The current code carries only one of the two components for
these bonds. The metadata says so honestly: `rhombic_status: "absent_no_primary_table"`
(`McConnellResult.cpp:187`), and the spec stance is *"carry a rhombic term only where a real value
exists, and mark its absence"* (`mcconnell_spec.md:105–106`).

This is the single place McConnell's emitted tensor is poorer than the physics it models. Closing it puts
**real `2e` structure exactly where it physically lives**, on a featured kernel — the opposite of
grinding McConnell down.

## 2. The construction (rides the existing per-bond machinery)

Extend the source tensor, **not** the propagator. The dipole propagator `D(r) = (3n⊗n − I)/r³` and the
contraction `response = D(r)·Q̂` (`McConnellResult.cpp:229–231`) are unchanged. Only `Q̂` gains a rhombic
term for the two flat categories:

```
Q̂_axial(û)        = (û ⊗ û − I/3)                          // today, all categories
Q̂_rhombic(û, m̂)   = Q̂_axial(û) + κ · (ê_in ⊗ ê_in − ê_out ⊗ ê_out)   // C=O, C=C only
```

where `m̂` is the **sp²-plane normal** of the bond's flat group, `ê_out = m̂` (out-of-plane), `ê_in = m̂ ×
û` (in-plane perpendicular), and `κ` is the **unit rhombicity shape weight** (dimensionless; the
`(ê_in⊗ê_in − ê_out⊗ê_out)` block is already traceless, so `Q̂_rhombic` stays traceless and is a clean
`2e` object to hand e3nn). The physical rhombic magnitude `Δχ_rh` is **not baked in** — it rides as the
learned/calibrated scale, the same `emit-unit-kernel / scale-rides-separately` pattern the rest of
McConnell and ring already follow.

**Why this is buildable without the held value.** Because `Δχ` is calibrated downstream (literature values
scatter 2–5×; `bond_anisotropy.md:128–132`), the producer emits the **geometric rhombic shape** in `Å⁻³`,
unscaled. The held/cited `Δχ_rh` value is needed only for the *literature-scaled physical-hypothesis*
channel (direction-4), not for the geometric emit. So the build and the citation decouple cleanly.

## 3. The new geometry — the sp²-plane normal (the real added work)

The axial term needs only `û` (the bond axis, already in `conf.bond_directions[bi]`). The rhombic term
needs a **second axis per flat bond**: the plane normal `m̂`. It is derived from topology the producer
already holds:

- **C=O (peptide & sidechain carbonyl):** the sp² plane of the carbonyl carbon — normal to the plane
  through `C`, `O`, and the carbon's third substituent (amide `N` / `Cα` for backbone; the two C
  neighbours for a ketone/Asn/Gln sidechain). `m̂ = normalize((r_O − r_C) × (r_X − r_C))`.
- **C=C:** the alkene plane through the two sp² carbons and one substituent each. `m̂` = normal to that
  plane.

This is a bounded, parameter-free geometric computation per rhombic bond, riding the existing
connectivity — **no new calculator, no new dependency.** Degenerate cases (collinear substituents,
missing third atom) fall back to `κ = 0` (axial only) and are flagged in metadata, never invented.

## 4. The value / citation status — the gate (have-it-to-cite-it)

The scholarship audit (held-corpus read, 2026-06-07) is unambiguous:

| Bond | What we HOLD on disk | Class | Rhombic? |
|---|---|---|---|
| **C=O** | peptide-group **axial** molar Δχ = **−5.36 ×10⁻⁶ cm³ mol⁻¹** (Pauling 1979, p. 2293) | HELD-PRIMARY | **axial only** |
| **C=O** | ester **8.8**, carboxyl **4.5 ×10⁻⁶** (Worcester 1978 p. 5475, quoting Lonsdale) | HELD-PRIMARY (value as reported); Lonsdale itself HELD-SECONDARY | axial |
| **C=O** | **Hooper & Kaiser 1965** amide rhombic sets (acetamide/formamide; EF-corrected & EF-neglected) | **HELD-PRIMARY** ✓ | **rhombic** (3 principal values) |
| **C=O** | **Abraham & Ainger 1999** ketone rhombic χ_∥−χ_⊥=10.30, χ_out−χ_⊥=1.93 ×10⁻⁶ cm³ mol⁻¹ | **HELD-PRIMARY** ✓ | **rhombic** |
| **C=O** | WA `+2.41/−5.42`, ApSimon `+12.65/−3.61`, Schneider `+14.45/−7.23` (two-component) | POINTER-ONLY (superseded by the held primaries above) | rhombic |
| **C=C** | **nothing** — Martin 2000 gives a GIAO shielding *surface* (Δσ ppm), not a Δχ, and abandons the cone model | HELD paper, no Δχ | **none, axial or rhombic** |

**Stated loudly (updated 2026-06-07 — C=O rhombic LANDED):** the **C=O two-component value is now
HELD-PRIMARY** — Hooper & Kaiser 1965 and Abraham & Ainger 1999 were fetched, ingested, and
provenance-verified (right papers). **C=C is still not held.** Held C=O values, three principal values
`(χ_out-of-plane, χ_along C=O, χ_in-plane⊥)` in `10⁻⁶ cm³ mol⁻¹`:

| source | set (assumption) | (χ_out, χ_∥, χ_⊥) | cite |
|---|---|---:|---|
| Hooper & Kaiser 1965 | acetamide, **EF-corrected** | `(−5.4, +4.0, −14)` | p. 2366 Table III |
| Hooper & Kaiser 1965 | acetamide, **EF-neglected** | `(+18, −19, −14)` | p. 2367 Table V |
| Hooper & Kaiser 1965 | formamide, EF-corrected | `(−13, +3.2, −5.0)` | p. 2366 Table IV |
| Abraham & Ainger 1999 | ketone C=O (relative) | `χ_∥−χ_⊥=10.30`, `χ_out−χ_⊥=1.93` | pp. 446–447 |

**The design call buried in the data (→ §7):** Hooper reports **EF-corrected** (pure magnetic anisotropy,
electric-field subtracted) and **EF-neglected** (effective χ that *absorbs* the carbonyl's electric-field
effect) sets — and they differ in **sign/pattern**, not just magnitude. Because we model the E-field
**separately** (charge/EFG), the physically-clean, **no-double-count** choice is **EF-corrected** — the
exact same logic that zeros aromatic McConnell against ring. Hooper found the EF-corrected fit *poorer*
(crude EF model), but for us EF-corrected avoids counting the carbonyl field twice. Sign A/B is unresolved
in Hooper; since the scale is **learned**, the fit absorbs sign, so a representative EF-corrected set is
pinned and calibration decides magnitude.

So the gate is **OPEN for C=O**: rhombic shape buildable now (§2–3) **and** a held, citable `Δχ_rh` exists
(EF-corrected set, pending the §7 double-count call). **C=C** remains: no held Δχ — ApSimon (acquire), the
held Martin **shielding-surface** (different object), or the **model-compound calc** (§5b).

## 5. Acquisition gate (route through codex download + scholarship vet)

Tiered by which bond each fixes. Acquire → ingest (`scripts/references/ingest_pdf.sh`) → summarise →
then the pointer rows become citable.

- **Tier 1 — C=O rhombic:** **Hooper & Kaiser**, Can. J. Chem. 43, 2363 (1965), DOI 10.1139/v65-318
  (multiple carbonyl component sets from amide spectra — the most on-target primary); **Abraham, Mobli &
  Smith** Part 19 (the comparison table all current pointers descend from — one acquisition upgrades ~5
  rows to HELD-SECONDARY); **Abraham & Ainger**, Perkin Trans. 2 1999, 441, DOI 10.1039/A808908F
  (ketone axial+perpendicular pair).
- **Tier 2 — peptide C=O / C–N (the production category):** **Williamson–Asakura** primary (the
  protein-predictor lineage; pin the exact citation on acquisition); **Abraham, Griffiths & Perez**, DOI
  10.1002/mrc.3920 (amides); **Lonsdale** Proc. R. Soc. A 171, 541 (1939) (upgrades the held
  ester/carboxyl numbers to HELD-PRIMARY).
- **Tier 3 — C=C:** **ApSimon** alkene anisotropy primary (the original of Martin's lineage) for a Δχ
  term; or the **Martin/Allen** algorithm primaries if the C=C term is recast as a shielding surface.

### 5b. Computed-proxy note — the ORCA route is closed (2026-06-07)

The ORCA-magnetizability cross-check is **ruled out**: our existing ORCA runs were not configured to
emit magnetizability, and there is no budget to re-run the protein campaigns for it. So the cited value is
**literature** (§5) — *or* the **one credible computed proxy**: a tiny standalone GIAO/CSGT
**magnetizability calc on model compounds** (formamide / N-methylacetamide → peptide C=O; acetone → ketone;
ethylene → C=C). That yields the full susceptibility tensor (three principal values → axial `Δχ_∥` +
rhombic `Δχ_⊥`) **directly, non-circular, held**, dodging the acquisition gate — and it rests on the same
model-compound→protein-bond transferability the literature values already assume (so it is at least as
defensible, and ours). It is a **seconds-scale calc on 3 small molecules**, not a protein campaign — a
different category from the "one more ORCA run" budget. **Not** credible as proxies: MOPAC polarizability
(electric, not magnetic) and the DFT CSA tensors (shielding `σ`, not susceptibility `χ`) — wrong physical
object.

## 6. Scope — the drop-nothing guarantee

- **Additive only.** New rhombic channels and the `m̂` geometry are *added*; every existing
  `mc_<category>_<channel>.npy`, the axial source shape, and all current per-atom fields are **byte-for-byte
  unchanged** for the axial categories. The flat categories gain the rhombic block; their existing axial
  emit is untouched (rhombic is a separate channel/array, not a mutation of the axial one).
- **One calculator.** Folds into `McConnellResult` (`bond_anisotropy.md:179` — extend the single
  calculator, do not stand up a second). Aromatic stays unconditionally zeroed (ring owns that physics).
- **Re-bless touches only the new arrays.** No change to existing baselines beyond the added channels;
  consistent with [[feedback_dont_overhold_builds_bless]].
- **Parity.** The rhombic block is a genuine `2e` susceptibility shape → even, `0e+1e+2e` like the rest of
  McConnell (consistent with the #1 parity sweep).

## 7. Open design questions (Jessica)

1. **Which categories get rhombic?** C=O backbone (PeptideCO) + sidechain (SidechainCO) + C=C — or
   backbone-only first? (Sidechain carboxylate vs amide is chemically broad — `MCCONNELL_DCHI_LITERATURE.md`
   flags it; split by topology, not IUPAC name.)
2. **EF-corrected vs EF-neglected Δχ — the double-count call (§4).** Hooper gives both; the EF-*neglected*
   set absorbs the carbonyl electric-field effect we *already* model in charge/EFG. Recommend
   **EF-corrected** (pure magnetic anisotropy, no double-count — same logic as aromatic↔ring), despite
   Hooper's poorer fit there. Confirm.
3. **C=C construction fork.** No held C=C Δχ exists. Acquire ApSimon (Δχ term, uniform with C=O), or use
   the **held Martin 2000 shielding surface** (different object, but on disk today)? These are not
   interchangeable — one is a susceptibility, one is a GIAO σ-increment field.
4. **Build-now vs value-gate.** Emit the rhombic *shape* now (additive, scale-learned, defensible as
   geometry), and let the cited `Δχ_rh` land when acquired — or hold the whole thing until Tier-1 is on
   disk? (**C=O Tier-1 is now on disk** — Hooper + Abraham held; the shape is safe to emit; the value is a
   separate, later, citable scale.)
5. **Second-axis robustness.** Confirm the sp²-plane derivation per category against the real topology
   (which third atom defines the plane for each carbonyl/alkene class).
6. **Where in the landing plan** (§8).

## 8. Placement relative to the landing state

The worktree currently holds **#1 parity hygiene, uncommitted** (build-confirmed); **#2 complete-emit /
#3 earn-"Haigh–Mallion" / #4 two-path** are queued and not started. The rhombic is naturally a **#5** —
and it pairs with #3/#4 as "the McConnell + ring tensor-richness pass." It is **value-gated**: the
geometric-shape build can proceed in the #2–#5 window, but the *cited* rhombic scale waits on Tier-1
acquisition. Recommendation: keep #1–#4 as the committed line; schedule the rhombic **shape** build after
#2 (it touches the same emit surface), fire the **acquisition** in parallel (codex), and wire the cited
value when it lands.

## 9. Unit landmines (carry into the build)

- Three systems all written "×10⁻⁶": molar cgs `10⁻⁶ cm³ mol⁻¹` (Pauling/Worcester/the note table — the
  **target**); per-molecule `10⁻³⁰ cm³ molecule⁻¹` ≡ `10⁻⁶ Å³ molecule⁻¹` (Abraham/ApSimon tables —
  where the pointers live); convert `q_molar = 0.602214076 · q_molecule`.
- Pauling/Worcester are **already molar cgs** — do not re-apply the 0.602 factor.
- McConnell prefactor (`MCCONNELL_DCHI_LITERATURE.md`): `δ_ppm = +0.5535130224 · q_molar · f`, flip sign
  for DFT shielding σ; tensor-form prefactor before the 1/3 contraction is `10²⁴/N_A = 1.660539067`.
- **CSA-vs-Δχ trap (load-bearing):** Cornilescu's ¹³C′ tensor and Baranac-Stojanović's σ components are
  **shielding (ppm)** tensors, not susceptibilities — never feed them through the Δχ prefactor.

## 10. Citations

- **Held (citable today):** Pauling, *Diamagnetic anisotropy of the peptide group*, PNAS 1979 (axial
  peptide Δχ); Worcester, *Structural origins of diamagnetic anisotropy in proteins*, PNAS 1978 (ester /
  carboxyl axial, quoting Lonsdale). `references/pauling-1979-…`, `references/worcester-1978-…`.
- **Pointers (acquire — §5):** Hooper & Kaiser 1965; Abraham Part 13 / Part 19 / mrc.3920; ApSimon;
  Williamson–Asakura; Schneider; Zürcher; Lonsdale 1939; Martin/Allen (C=C surface).
- **Scholarship-hygiene follow-up:** `MCCONNELL_DCHI_LITERATURE.md` lists the unheld primaries as "Source
  Anchors" as if they back the table — relabel those POINTER-ONLY and add the two genuinely-held axial
  anchors (Pauling 1979, Worcester 1978) it currently omits. [scholarship workstream]

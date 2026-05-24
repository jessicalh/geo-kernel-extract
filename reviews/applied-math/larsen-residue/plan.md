# Fix plan — LarsenResidue.{h,cpp}

## 1. Summary

`LarsenResidue.h/.cpp` implement bond-graph + K=3 Weisfeiler-Lehman
perception of atom identity inside a Larsen ProCS15 tripeptide DFT
record. This is **perception, not a tensor producer**: the file assigns
typed `AtomMechanicalIdentity` to each DFT atom and copies the
pre-decomposed shielding tensor through unchanged (`EmitPiece` lines
946-949 just copy `shielding_tensor` / `isotropic` / `anisotropy` /
`t2_components` from the DB row). It computes **no signs and no physical
numbers** of its own. The only literals are geometric bond cutoffs
(Å), two FNV-1a hash constants, a 31-bit mask, and amide-detection
distance windows — none of these flow into a shielding value; they only
decide graph topology.

The file mostly tells a coherent story. The header narrates the 5
pieces + 6-step algorithm before any declaration; the .cpp is
banner-sectioned in exactly that order (bond graph → amide ID →
components → canonical chemistry → WL match → emit → assemble). Both
reviews agree the story holds and neither found a definite correctness
bug. The friction is concentrated in three spots: the WL signature core
(symbol-soup locals, a mid-expression magic mask, two near-identical
`*WLSignatures` functions), the `ParentHeavyForH` ladder, and the HIS
variant block. The header carries heavy history/prototype prose mixed
with the live contract.

**This pass will:** improve internal local-variable names in the WL and
amide code, add 2-4 word signposts in the WL / stamping / HIS / emit
blocks, name the bond-cutoff and amide-distance literals and the WL
mask, drop one dead type alias, and trim a runtime FATAL string of its
test-suite prose. **This pass will not:** touch any number, sign,
cutoff value, hash constant, or the algorithm; rename any output /
serialized name; rename the cross-file slot fields (`N_idx` … `O_idx`)
without flagging the carry-through cost (see §3); add abstractions.

**Boundary trace (done up front, governs every rename below):**
- No name in `LarsenResidue.{h,cpp}` is a serialized contract.
  `python/nmr_extract/_catalog.py` keys (`tripeptide_bb_*`,
  `tripeptide_neighbor_*`, `larsen_hbond_*`) are emitted by
  `TripeptidePoseAssembler` / the `*ShieldingResult` classes, **not** by
  this file. Nothing here crosses into NPY/H5/SDK names.
- The **public** symbols consumed across `src/`:
  `PerceiveLarsenTripeptide`, `LarsenTripeptide` + its fields
  (`ace/n_cap/central/c_cap/nme`, `central_variant_index`, `TotalAtoms`,
  `FindByDftIdx`), `LarsenResidue` + its slot fields
  (`N_idx/H_idx/CA_idx/HA_idx/CB_idx/C_idx/O_idx`), `PerAtom` +
  `canonical_assignment_ambiguous`, `LookupByIdentity`,
  `HasAllRequiredSlots`. Consumers: `TripeptidePoseAssembler.cpp`,
  `TripeptideNeighborShieldingResult.cpp`, `TripeptideDftTable.cpp`.
  These are a header contract — renames carry cross-file cost and are
  declined or flagged below.
- Everything in the anonymous namespace of the .cpp (the entire
  perception engine: `BuildBondGraph`, `FindPeptideAmides`,
  `CanonicalPiece`, `MatchPiece`, `EmitPiece`, the WL functions, all
  locals) is **file-local** — free to rename for clarity.
- The `canonical_assignment_ambiguous` contract **is honored** by the
  consumer: `TripeptidePoseAssembler.cpp:479` reads
  `perc.canonical_assignment_ambiguous` into `relaxed` and drops
  BranchAddress/DiastereotopicIndex exactly as the header (h:99-111)
  and `MatchPiece` (cpp:692-710) promise. Producer and consumer agree.

---

## 2. Review-finding ledger

Every finding from both reviews, with disposition. (No
`codex-correctness.md` present.)

### codex.md

| # | Finding | Disposition |
|---|---------|-------------|
| C1 | h:21 header opens with prototype/history prose before the data model; compress, move Python/prototype history to design doc | **Adopted** → §3 E1 (trim, keep pointer to design doc) |
| C2 | h:99 ambiguity-flag comment packs too much (WL classes, automorphism, BranchAddress, DI, matcher policy); reduce | **Declined** — Claude (CL-h99) judges every sentence earns its place because it defines the whole downstream matching contract that lives outside this file; I agree (the consumer at PoseAssembler:453-478 restates exactly this). Trimming would orphan the contract. Keep; optional one-line tightening of the parenthetical only (§3 E2). |
| C3 | cpp:357 `ParentHeavyForH` long order-dependent if-chain; add group labels (`// backbone H`, `// ring H`, etc.) | **Adopted (partial)** → §3 E3. One header line stating prefix-match + order-matters; light group labels. Not full fragmentation. |
| C4 | cpp:534 WL signature model hard to reconcile: comment says neighbor degree, code stores masked prior-round hashes; rename + restate | **Adopted** → §3 E4 + Usage note U1. The `Signature` type genuinely mixes round-0 (degree) and round-k (folded prior hash) meaning of its `int`/pair fields. |
| C5 | cpp:1138 canonical-cache lambdas interrupt the narrative; add `// canonical templates` signpost | **Adopted** → §3 E9 |
| C6 | cpp:1205 HIS branch fuses policy/variant-order/match/logging; split with labels | **Adopted** → §3 E10 (signpost labels only, no restructuring) |
| C7 | h:129 `N_idx`/`H_idx` etc. omit "local atom index for role slot"; rename to `nitrogen_local_idx` etc. | **Declined (rename)** — cross-file contract consumed in `TripeptidePoseAssembler.cpp:247-261, 387-401`. Carry-through cost outweighs gain; the header comment at h:129 already says "Role-pinned slot cache" + "local index." See §4 / §3 E2-note. Flagged, not done. |
| C8 | cpp:128 `[c, nj]` / `n_heavy_other` hide chemistry; use `carbonyl_c` / `amide_n` / `noncarbonyl_heavy_neighbors` | **Adopted** → §3 E5 |
| C9 | cpp:172 `p`/`pi`/`ai` premature compression; use `piece`/`piece_idx`/`atom_idx` | **Adopted** → §3 E6 |
| C10 | cpp:576 `cur`/`next`/`nbr_sig`/`s` symbol soup; use `previous_round_hashes`/`next_round_hashes`/`neighbor_labels`/`round_signature` | **Adopted** → §3 E7 (file-local; applied to BOTH WL functions for parallelism) |
| C11 | cpp:1220 `v`/`tried`/`any_match` obscure HIS state; use `variant_name`/`attempted_variants`/`matched_variant` | **Adopted** → §3 E10 (note: `variant_name` collides with `CanonicalPiece::variant_name` semantics; use `variant_code` to stay distinct) |
| C12 | cpp:48 bond cutoffs embedded literals; name physical thresholds (`kHydrogenHeavyCutoffAngstrom`, `kSulfurSulfurCutoffAngstrom`) | **Adopted** → §3 E11 + Usage note U2 |
| C13 | cpp:103 amide magic numbers `1.32`/`1.20`/`1.50`; name as C=O and peptide C-N windows | **Adopted** → §3 E12 + Usage note U3 |
| C14 | cpp:561 WL should read seed → refine → group; name `kWlLabelMask` + prior-label step | **Adopted** → §3 E4/E7/E8 |
| C15 | cpp:930 `EmitPiece` fuses payload copy / identity / ambiguity / slot fill; add `// atom payload` + `// role slots` | **Adopted** → §3 E13 |
| C16 | cpp:1240 HIS matching fuses count filter + graph match; label `// atom-count prefilter` + `// graph match` | **Adopted** → §3 E10 |
| C17 | cpp:67 `Distance` omits units → `DistanceAngstrom` | **Adopted** → §3 E14 (file-local free function) |
| C18 | cpp:97 `FindPeptideAmides` returns C-N bond pairs not groups → `FindPeptideAmideBonds` | **Adopted** → §3 E15 (file-local; rename callsite at cpp:1091) |
| C19 | cpp:543 `HashSignature` too generic → `HashWlSignature` | **Adopted** → §3 E16 (file-local) |
| C20 | cpp:651 `MatchPiece` → `MatchPieceByWLSignature` | **Declined** — Claude (CL-651) judges `MatchPiece` fine because the comment block at cpp:627-645 makes the WL-and-index contract explicit; the longer name is verbose at its 3 callsites. Weighed call: keep `MatchPiece`, the signpost already carries it. |
| C21 | cpp:901 `EmitPiece` vague → `PopulateMatchedPiece` / `EmitMatchedResidue` | **Declined** — `EmitPiece` is the consistent verb paired with the `*Piece` family (`MatchPiece`, `OrderPiecesAlongAmideChain`, `CanonicalPiece`). The §1 narrative banner ("Per-piece emission") plus the doc comment carry intent. Renaming one of a coherent family adds churn without net clarity. |
| C22 | h:21 verbose process comment → `// perception pipeline` | **Duplicate** of C1. |
| C23 | h:99 dense block → `// WL-equivalent atom` | **Duplicate** of C2 (declined). |
| C24 | cpp:1 file banner mostly process/history → `// Larsen tripeptide perception` | **Declined** — the banner's provenance (Python POC validated against 20 DB combos + AAA log byte-for-byte) is load-bearing project context a reviewer wants; it is not interfering with reading code. Trim only redundant repetition with the header (§3 E17, light). |
| C25 | cpp:512 long WL rationale interrupts impl; keep `// WL refinement`, move detail to docs | **Declined** — Claude (CL-508) explicitly calls this block "excellent grounding; the CD1/CD2 example is exactly what a chemist needs; keep." Two reviewers disagree; I side with keep — this is the one place the *why* of K=3 lives, and the brief warns against over-trimming. |
| C26 | cpp:536 comment stale (says neighbor degree, code uses prior-round hash) → `// hashed neighbor labels` | **Adopted** → §3 E4 (this is the comment-conforms-to-code fix; see U1) |
| C27 | cpp:789 20+ lines stamping prose too much; reduce to `// typed identity stamping` + `// methyl-H collapse` + `// topology lookup` | **Adopted (partial)** → §3 E8. Add the two inner signposts; keep the contract paragraph (the methyl-collapse mirror-of-ComposeAtomSemantic rule is non-obvious and consumed nowhere else), trim only restated history. |
| C28 | cpp:857 fatal string contains test/process prose; keep identifiers + mismatch reason, move test claims out of runtime message | **Adopted** → §3 E18 (matches Claude CL-860) |
| C29 | cpp:1206 HIS policy comment longer than branch needs → `// HIS variant policy` | **Adopted (partial)** → §3 E10. Keep the hint-vs-no-hint contract (it documents a silent-misassignment hazard); add the short signpost above the loop. |
| C30 | cpp:583 (Correctness) check whether 31-bit WL hash truncation is intentional/documented; collision tradeoff unexplained | **Adopted as comment** → §3 E4 + Usage note U4. Traced: intentional and consistent on both sides (cpp:583 and cpp:614). Not a bug; document the fold. |

### claude.md

| # | Finding | Disposition |
|---|---------|-------------|
| L1 | cpp:565-624 `CanonicalWLSignatures`/`PerceivedWLSignatures` same algorithm twice; add signpost on each stating deliberate parallel impls | **Adopted** → §3 E7-note. Add `// canonical side, dense vector-indexed` / `// perceived side, sparse dft-indexed` headers. (Brief forbids refactor/merge.) |
| L2 | cpp:578-584 `cur[j] & 0x7fffffff` mid-loop; name or comment the 31-bit fold | **Duplicate** of C30/C14 → §3 E4. |
| L3 | cpp:543-558 `HashSignature` bare FNV magic; 2-word `// FNV-1a 64-bit` label on seed line | **Adopted** → §3 E16-note (the inner comment already says "FNV-1a 64-bit"; move/echo it to the seed-constant line). Minor. |
| L4 | cpp:360-400 `ParentHeavyForH` prefix-vs-exact ladder; existing per-line comments adequate, consider one header comment | **Duplicate** of C3 → §3 E3. |
| L5 | cpp:1147-1158 `cached_residue` bare `21` + `idx=20` fallback; comment `// 20 standard AAs + index 20 = Unknown` | **Adopted** → §3 E9-note. Traced: `AminoAcid` enum is 0-19 standard + 20=Unknown (distinct from `kBaseVariantIdx=255`). Comment it. |
| L6 | cpp:1205-1270 HIS variant block densest control flow; 3-word signpost before loop | **Duplicate** of C6/C29 → §3 E10. |
| L7 | cpp:128 `nj` reads as "n-sub-j" (reused 181, 1101); use `n_idx`/`amide_n` | **Duplicate** of C8 → §3 E5. Note all three sites (128, 181, 1101). |
| L8 | cpp:254-265 `AddBond` `A`/`B` (sets) vs `a`/`b` (params) one-char distinction; rename sets `adj_a`/`adj_b` | **Adopted** → §3 E19 |
| L9 | cpp:540 `SigKey` alias declared but unused; drop if dead | **Adopted** → §3 E20. Confirmed dead by grep (only declaration site; code uses `std::vector<std::uint64_t>` / `std::map<int,std::uint64_t>` directly). Remove. |
| L10 | cpp:534 `Signature` tuple `int` is degree, pair is (nbr-elem, nbr-prior-hash-folded); alias could note `// (element, degree, neighbour multiset)` | **Duplicate** of C4 → §3 E4. |
| L11 | cpp:705 `ambiguous` well-named; good | **No-op** (positive). |
| L12 | h:130-136 `N_idx`/`CA_idx` slot names clear and domain-correct; no change | **Adopted (as decline of C7)** — confirms §4 verdict to keep slot names. |
| L13 | cpp:573-591 WL round body fuses gather→sort→hash; blank-line separation + named sequence | **Adopted** → §3 E7 (named intermediates + blank lines achieve this). |
| L14 | cpp:97-138 `FindPeptideAmides` well-staged (candidates → filter); good | **No-op** (positive). |
| L15 | cpp:812-894 `StampChainIdentitiesViaTable` clear 4-step; methyl-collapse block labeled; good | **No-op** (positive; the C27 inner signposts reinforce). |
| L16 | cpp:684-711 cardinality check + within-class assignment two stages; blank line + `// validate` / `// assign` | **Adopted** → §3 E21 |
| L17 | cpp:48 `BondCutoffAngstrom` exemplary naming | **No-op** (positive). |
| L18 | cpp:169 `OrderPiecesAlongAmideChain` good | **No-op** (positive). |
| L19 | cpp:651 `MatchPiece` fine, contract explicit at 627-645 | **Adopted (as decline of C20)** — confirms keep. |
| L20 | cpp:236 `FinalizeAdjacency` vague; `BuildIndexAdjacency` says it on the tin | **Adopted** → §3 E22 (file-local; rename + callsites at 302, 322, 439, 497). |
| L21 | cpp:565/597 `CanonicalWLSignatures`/`PerceivedWLSignatures` clear and parallel; good | **No-op** (positive). |
| L22 | cpp:508-532 WL "why multi-round" block excellent; keep | **Adopted (as decline of C25)** — confirms keep. |
| L23 | cpp:860-875 FATAL message embeds ~15 lines process/history; trim to actionable line, move rationale to comment | **Duplicate** of C28 → §3 E18. |
| L24 | cpp:989 `// heuristic; refinement deferred` good terse signpost | **No-op** (positive; keep). |
| L25 | h:99-110 ambiguity comment long but every sentence earns place; acceptable | **Adopted (as decline of C2)** — confirms keep. |
| L26 | cpp:49 `// Cutoffs match the Python POC.` + 328/359 history prose; reassures but doesn't help math; keep one at file top, redundant per-function | **Adopted (partial)** → §3 E11/E8-note. Keep file-top provenance; the per-function "matches Python POC" lines stay only where they're the sole provenance, drop where redundant. Light. |
| L27 | cpp:1138-1140 `// Thread-safe per C++11 ...` correct and useful; keep | **No-op** (positive; keep). |
| L28 | cpp:578-583 (Correctness) 31-bit fold raises collision prob vs full 64-bit; confirm intended; consistent with perceived side at 614 | **Duplicate** of C30 → Usage note U4 (traced: intentional, consistent, harmless at ~50 atoms). |
| L29 | cpp:112/121 (Correctness) amide cutoffs `<1.32`, `>1.20 && <1.50` geometry-dependent; check they hold for DFT-optimized geometries (r²SCAN/PBE) in DB, not just POC inputs | **Question** → §6 Q1. Cannot confirm from this file; needs DB geometry check. Not a bug claim. |
| L30 | cpp:705-708 (Correctness) multi-atom within-class binding by container order is arbitrary; correct only if every downstream consumer honors the flag | **Resolved (coherent)** → Usage note U5. Traced consumer: `TripeptidePoseAssembler.cpp:479` honors it. Not a bug. |
| L31 | cpp:1241/1273 (Correctness) HIS path checks size before `EmitPiece`, non-HIS relies on EmitPiece's internal check; harmless redundancy | **Resolved (coherent)** → Usage note U6. Confirmed redundant, not a bug; no change. |
| L32 | cpp:1014-1039 (Correctness) ACE `HasAllRequiredSlots` returns `atoms.size()==6` only; comment claims "Required slots: C, O" but C/O slots not verified | **Adopted as comment fix** → §3 E23 + Usage note U7. The comment overstates; code conforms-to-reality is `size==6` only. Either tighten the comment or add the slot checks. Recommend comment fix (no number/behavior change); flag the stronger option as a question (§6 Q2). |
| L33 | cpp:970-971 (Correctness) CB slot `C && Beta`; GLY exempt at 1035; consistent, no bug | **No-op** (positive; confirmed coherent). |
| L34 | cpp:1055 `FindByDftIdx` linear-scans 5 pieces; fine for ~50 atoms; noted only | **Declined** — not a readability or correctness issue; O(N) over ~50 atoms is fine and the brief forbids refactors. No change. |

---

## 3. Edits that don't move numbers

All edits below are comment fixes, signposts, named literals, named
intermediates, or **file-local** renames. None changes a value, sign,
cutoff, hash constant, the 31-bit mask, or control flow. No
output/serialized name touched.

**Named literals (file-local constants, identical values):**

- E11 — `LarsenResidue.cpp:48-64` (`BondCutoffAngstrom`) — introduce
  named `constexpr double` for each cutoff: e.g.
  `kHydrogenSulfurCutoffA = 1.50`, `kHydrogenHeavyCutoffA = 1.30`,
  `kSulfurSulfurCutoffA = 2.30`, `kSulfurHeavyCutoffA = 2.10`,
  `kOxygenOxygenCutoffA = 1.60`, `kHeavyHeavyCutoffA = 1.85`. Values
  unchanged. Keep the file-top "match Python POC" provenance; drop the
  per-function restatement at line 49 (covered by E17). (C12, L26)
- E12 — `LarsenResidue.cpp:112,121` (`FindPeptideAmides`) — name the
  amide-geometry windows: `kCarbonylCODoubleBondMaxA = 1.32`,
  `kPeptideCNMinA = 1.20`, `kPeptideCNMaxA = 1.50`. Values unchanged.
  (C13)
- E4 — `LarsenResidue.cpp:583,614` — name the 31-bit fold:
  `constexpr std::uint64_t kWlPriorHashFoldMask = 0x7fffffffull;` and
  use at both WL sites. Add comment at the seed/first-use:
  `// fold prior-round hash to 31 bits so it pairs with Element in the`
  `// neighbour-label multiset; both WL sides fold identically`.
  (C14, C30, L2, L10)

**Comment fixes (comment conforms to code):**

- E1 — `LarsenResidue.h:21-31` — trim the prototype/Python-history
  paragraph to ~2 lines + the existing design-doc pointer; keep the
  6-step list (it is the live contract). (C1, C22)
- E4 — `LarsenResidue.cpp:534-540` — `Signature` alias comment: state
  that the `int` is **degree** and the pair is **(neighbour element,
  neighbour prior-round folded hash)**; at round 0 the multiset is
  empty and `int` is 0 (seed). The current comment at 536-539 says
  "(neighbour-element, neighbour-degree)" which is stale for rounds
  ≥1 — fix to "(neighbour-element, neighbour folded prior-round hash)".
  (C4, C26, L10)
- E17 — `LarsenResidue.cpp:1,49,328,359` — drop redundant
  "matches/mirrors the Python POC" restatements where the file-top
  banner already says it; keep one provenance statement at file top.
  Do **not** gut the file banner (C24 declined). (L26)
- E18 — `LarsenResidue.cpp:860-875` — trim the FATAL `fprintf` body to
  the actionable line ("canonical chemistry in LarsenResidue.cpp does
  not match the generated topology table; update one to match the
  other; positive coverage:
  `LarsenResiduePerceptionTest.AllCombinationsPerceiveCleanly`") plus
  the identifiers/identity dump (keep those — they are the debug
  payload). Move the "why this is unreachable from the standard-20
  substrate" paragraph to a code comment above the block. (C28, L23)
- E23 — `LarsenResidue.cpp:1017-1019` — fix the `Kind::AceCap` comment:
  it claims "Required slots: C, O" but the check is `atoms.size()==6`
  only. Reword to state what the code actually guarantees ("ACE: 6
  atoms total; the carbonyl C/O roles are stamped by
  StampCapIdentities and not re-verified here"). See §6 Q2 for the
  stronger alternative. (L32)

**Signposts (2-4 words, no fragmentation):**

- E3 — `LarsenResidue.cpp:360` — one header comment on `ParentHeavyForH`:
  `// Prefix-match H name → heavy parent; order matters (specific`
  `// before general). H/HA backbone, then β/γ/δ/ε/ζ/η by shell.`
  Keep existing per-line residue comments. (C3, L4)
- E8 — `LarsenResidue.cpp:831,852` — inside `StampChainIdentitiesViaTable`
  add `// methyl-H collapse` and `// generated-table lookup (authority)`
  signposts; trim restated history but keep the contract paragraph and
  the ComposeAtomSemantic line reference. (C27, L15)
- E9 — `LarsenResidue.cpp:1141` — `// canonical templates (cached once
  per kind)` signpost above the cached lambdas. At cpp:1156 add
  `// AminoAcid 0-19 standard, 20 = Unknown fallback`. (C5, L5)
- E10 — `LarsenResidue.cpp:1217,1235,1240,1251` — signposts in the HIS
  block: `// variant try-order: strict hint, else HID/HIE/HIP`,
  `// atom-count prefilter`, `// graph match`, `// no variant matched`.
  Trim the 1206-1216 prose to keep the silent-misassignment hazard note
  + the new signpost. No restructuring. (C6, C16, C29, L6)
- E13 — `LarsenResidue.cpp:935,955` — in `EmitPiece` per-atom loop:
  `// atom payload (copy DB row through)` and `// role slots (typed
  dispatch, no name compare)`. (C15)
- E21 — `LarsenResidue.cpp:684,702` — in `MatchPiece`: `// validate:
  signature sets + per-class cardinality` and `// assign: bind
  perceived → canonical within each class`. (L16)
- E7-note / L1 — `LarsenResidue.cpp:565,597` — one line on each
  `*WLSignatures`: `// canonical side: dense, vector-indexed by node`
  and `// perceived side: sparse, keyed by dft atom index — same
  refinement as CanonicalWLSignatures`. (L1)
- E2 — `LarsenResidue.h:99-111` — keep the block; optionally shorten the
  parenthetical example list to one line. Low priority. (C2 partial)

**Named intermediates (file-local locals):**

- E5 — `LarsenResidue.cpp:128-134,181,1101` — rename `nj`→`amide_n`,
  `c`→`carbonyl_c`, `n_heavy_other`→`noncarbonyl_heavy_neighbours`.
  All three structured-binding sites (128, 181, 1101). (C8, L7)
- E6 — `LarsenResidue.cpp:172-176` (`OrderPiecesAlongAmideChain`) —
  `p`→`piece`, `pi`→`piece_idx`, `ai`→`atom_idx`. (C9)
- E7 / L13 — `LarsenResidue.cpp:568-591` and `602-622` (both WL
  functions, for parallelism) — `cur`→`prev_round_hash`,
  `next`→`next_round_hash`, `nbr_sig`→`neighbour_labels`,
  `s`→`round_signature`. Add blank lines between the gather / sort /
  hash steps so the "symmetrise then hash" shape is visible. (C10, L13)
- E10 — `LarsenResidue.cpp:1220-1247` — `v`→`variant_code`,
  `tried`→`attempted_variants`, `any_match`→`matched`. (C11)

**File-local function renames (no header, no cross-file):**

- E14 — `LarsenResidue.cpp:67` `Distance`→`DistanceAngstrom`; callsites
  at 80, 112, 120. (C17)
- E15 — `LarsenResidue.cpp:97` `FindPeptideAmides`→`FindPeptideAmideBonds`;
  callsite at 1091. (C18)
- E16 — `LarsenResidue.cpp:543` `HashSignature`→`HashWlSignature`;
  callsites at 571, 588, 605, 619. Echo the "FNV-1a 64-bit" label onto
  the seed-constant line at 550. (C19, L3)
- E22 — `LarsenResidue.cpp:236` `FinalizeAdjacency`→`BuildIndexAdjacency`;
  callsites at 302, 322, 439, 497. (L20)
- E19 — `LarsenResidue.cpp:255-264` (`AddBond`) — adjacency sets
  `A`→`adj_a`, `B`→`adj_b` (keep params `a`/`b`). (L8)

**Deletion:**

- E20 — `LarsenResidue.cpp:540` — remove the unused `SigKey` type alias
  (confirmed dead; nothing references it). (L9)

---

## 4. Usage notes (the reasons, traced)

This file produces no signs or physical numbers; the "values" the
reviews flagged are graph-decision literals and a hash. Each traced:

- **U1 — the WL `Signature` `int` field (degree vs folded hash).**
  *Reason:* round 0 seeds with `(element, 0, {})` — the `int` is a
  placeholder 0, no neighbours folded yet (cpp:570, 604). Rounds 1-3
  build `(element, degree, sorted multiset of (neighbour-element,
  neighbour-folded-prior-hash))` (cpp:586-588, 617-619). The `int` in
  the *pair* is the neighbour's folded prior-round hash, **not** its
  degree — the original comment at 536-539 was stale (it predates the
  prior-hash refinement). Comment fixed to match code (E4).
  *Consumed:* only inside `MatchPiece` (cpp:667-710) to group nodes
  into signature classes. Never serialized.

- **U2 — bond cutoffs (Å).** *Reason:* element-pair covalent-distance
  thresholds that define the perceived bond graph; "match the Python
  POC" that was validated against all 20 DB (tripeptide, frame_type)
  combinations + the AAA Gaussian log byte-for-byte (file banner).
  *Consumed:* only `BuildBondGraph` (cpp:80) → graph topology →
  amide cut / components / WL. No tensor or shielding value depends on
  them. Producer-only; no cross-consumer to reconcile.

- **U3 — amide-detection windows (1.32 / 1.20 / 1.50 Å).** *Reason:*
  geometric discriminators for "this C=O + C-N pair is a peptide
  amide" (cpp:112, 121). Same POC provenance. *Consumed:* only
  `FindPeptideAmideBonds` → the 4-amide count gate and the
  cut-into-5-pieces step (cpp:1092, 1101). See §6 Q1 for the one open
  question (do these windows hold for the r²SCAN/PBE DB geometries).

- **U4 — 31-bit fold (`& 0x7fffffff`).** *Reason:* the WL refinement
  pairs each neighbour's prior-round hash with its element into a
  `std::pair<Element,int>`; the `int` truncation folds the 64-bit hash
  into the pair's `int` slot. *Traced both sides:* canonical at
  cpp:583 and perceived at cpp:614 fold **identically** — so the two
  signature spaces are comparable, which is what `MatchPiece` requires.
  At ~50 atoms/piece the collision risk is negligible. **Coherent
  (expected)**; not a bug. Documented via E4, not changed.

- **U5 — within-class binding by container order (cpp:705-708).**
  *Reason:* for multi-atom (graph-automorphic) WL classes, the
  perceived→canonical assignment is arbitrary by design, and
  `canonical_assignment_ambiguous=true` is set to signal it.
  *Consumer traced:* `TripeptidePoseAssembler.cpp:479` reads the flag
  into `relaxed`, and `candidate_protein_atoms(..., relaxed)` drops
  BranchAddress + DiastereotopicIndex (PoseAssembler:442, 480-481),
  then resolves by **nearest-spatial** (PoseAssembler:493-499) — exactly
  the contract h:99-111 and MatchPiece cpp:692-701 promise. **Producer
  and consumer agree.** Coherent; not a bug (L30 resolved).

- **U6 — HIS size-check redundancy (cpp:1241 vs EmitPiece cpp:907).**
  The HIS path prefilters on `piece.size()==canon.atoms.size()` before
  calling `EmitPiece`, which re-checks the same at 907. The prefilter
  lets the loop `continue` to the next variant cheaply (so a wrong-size
  variant doesn't run a full WL match). Harmless, intentional.
  Coherent; no change (L31 resolved).

- **U7 — ACE `HasAllRequiredSlots`.** *Reason:* ACE has no
  backbone N/CA/H/HA/CB slots in the chain sense; its carbonyl C/O
  identities are stamped by `StampCapIdentities` (cpp:778-786) and
  `EmitPiece`'s typed dispatch (cpp:968-969) sets `C_idx`/`O_idx`.
  The `atoms.size()==6` gate is the only post-emit check for ACE; the
  comment claiming "Required slots: C, O" overstates what's verified.
  Code-conforms fix is to correct the comment (E23). Stronger option
  (actually verifying C_idx/O_idx ≥ 0) is a behavior change → §6 Q2.

---

## 5. Bug-by-exhaustion candidates

**None.** No finding survived tracing as a producer-side bug. The two
correctness items the reviews raised about consumer-honored conventions
(L30 ambiguity flag, L31 size redundancy) traced as **coherent** with
the consumer in `TripeptidePoseAssembler.cpp`. The 31-bit fold (C30/L28)
traced as intentional and consistent. The remaining correctness flags
are a comment overstatement (L32 → comment fix E23, not a number) and a
DB-geometry question (L29 → Q1) that this file cannot answer.

---

## 6. Questions & Ambiguities

- **Q1 (from L29) — do the amide-detection distance windows
  (C=O < 1.32 Å, peptide C-N in 1.20-1.50 Å) hold for the
  DFT-optimized geometries actually in `tensorcs15`?** The values were
  validated against the Python POC's inputs and the AAA log; the file
  banner asserts all-20-combo coverage by
  `LarsenResiduePerceptionTest.AllCombinationsPerceiveCleanly`, which
  would catch a window that excluded a real amide (it would fail the
  4-amide gate). So the *test* is the empirical answer — but I did not
  re-run it against the live r²SCAN/PBE rows in this pass. If the
  fix-implementer wants certainty before touching the comment, run that
  test. Not a bug claim; a verification ask. No change proposed to the
  values regardless (governing prior: numbers earned by iteration).

- **Q2 (from L32) — should ACE `HasAllRequiredSlots` verify `C_idx ≥ 0
  && O_idx ≥ 0` in addition to `atoms.size()==6`?** This would be a
  behavior change (could turn a currently-passing malformed-ACE into a
  decline), so it is out of scope for a readability pass under the
  brief. The default plan (E23) fixes the comment to match the
  weaker-than-claimed check. Flagging for the human: if the intent was
  always "C and O slots required," the slot check is the faithful
  implementation and the comment was the spec — but that is the human's
  call, not mine to make in a no-numbers pass.

- **Q3 — review disagreement on the two big comment blocks.** codex
  wants C25 (cpp:508-532 WL rationale) and C2 (h:99 ambiguity) trimmed;
  Claude wants both kept verbatim as load-bearing. I sided with keep on
  both (the brief warns over-trimming is its own unclarity, and these
  are the only places the *why* of K=3 and the downstream matching
  contract live). If the human prefers codex's leaner header, E2/C25
  can be revisited — but I'd keep the CD1/CD2 worked example regardless.

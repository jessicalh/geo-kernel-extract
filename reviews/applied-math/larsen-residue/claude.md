# larsen-residue — claude review (readability focus)

- **Targets:** src/LarsenResidue.{h,cpp}
- **Model:** Claude (general-purpose agent, opus)
- **Date:** 2026-05-24
- **Brief:** `../_brief-claude.md`

---

## Verdict

**LarsenResidue.h** — Reads cleanly. The header narrates the whole perception pipeline (5 pieces, 6-step algorithm, failure modes) before any declaration, and every field carries a units/semantics comment. A chemist who has never seen the codebase can follow the through-line top to bottom. The only friction is that the per-atom comment block on `canonical_assignment_ambiguous` is doing heavy conceptual lifting that the reader must absorb before the algorithm makes sense — but it is doing it correctly.

**LarsenResidue.cpp** — Mostly coherent, organized into clearly-banner-sectioned stages (bond perception → amide ID → components → canonical chemistry → WL matching → emission → assembly) that mirror the header's 6 steps. The story holds. Where it breaks is in the WL signature core: a few load-bearing magic constants (the `0x7fffffff` truncation, the FNV seeds) and the dual near-identical `*WLSignatures` functions force the reader to diff two blocks to confirm they compute the same thing. The `ParentHeavyForH` table and the HIS `try_order` block are long but honest. Correctness looks sound modulo two conventions I flag below.

---

## 1. Coherent story / readability (primary)

- LarsenResidue.cpp:565-624 — `CanonicalWLSignatures` and `PerceivedWLSignatures` are the same algorithm written twice (one over `std::vector`+`adj_by_idx`, one over `std::map`+`AdjacencySet`); a reader must diff them line-by-line to confirm they refine identically — add a one-line signpost on each (`// canonical side, vector-indexed` / `// perceived side, sparse dft-indexed`) stating they are deliberate parallel implementations of one WL refinement.
- LarsenResidue.cpp:578-584 — `cur[j] & 0x7fffffff` mid-loop: as written, a reader must guess why a 64-bit hash is truncated to 31 bits before re-hashing — name it or comment `// fold prior hash into 31 bits to fit the pair<Element,int> discriminator`.
- LarsenResidue.cpp:543-558 — `HashSignature` uses bare FNV-1a magic constants (`0xcbf29ce484222325`, `0x100000001b3`); fine, but a 2-word label `// FNV-1a 64-bit` on the seed line saves the reader a lookup (the inner comment exists but the constants read as noise without it adjacent).
- LarsenResidue.cpp:360-400 — `ParentHeavyForH` is a 40-line ladder of `rfind(...,0)==0` prefix tests; the logic is sound but the reader must hold "is this a prefix match or exact?" per line — the existing per-line residue comments (`// Thr HG1`, `// Asn`) carry it adequately, so this is acceptable, but consider a single header comment `// prefix-match H name → heavy parent; order matters (specific before general)`.
- LarsenResidue.cpp:1147-1158 — `cached_residue` builds a 21-entry table with a bare `21` and `idx=20` fallback; reader cannot tell 20==Unknown sentinel without cross-referencing the enum — comment `// 20 standard AAs + index 20 = Unknown fallback`.
- LarsenResidue.cpp:1205-1270 — the HIS variant block is the densest control flow in the file (hint vs. no-hint, try_order, tried-string accumulation, matched_variant_idx); the story is followable but long — a 3-word signpost before the loop (`// strict hint or first-fit fallback`) would orient the reader faster than the prose comment at 1206-1216, which mixes contract and rationale.

## 2. Naming carries meaning

- LarsenResidue.cpp:128 — `nj` for the amide-N index (and reused at 181, 1101) reads as "n-sub-j"; `n_idx` or `amide_n` is clearer for a chemist scanning carbonyl-C/amide-N pairs.
- LarsenResidue.cpp:254-265 — `AddBond` locals `A`/`B` (capital, the adjacency sets) vs. params `a`/`b` (the names) is a one-char distinction doing real work; rename the sets `adj_a`/`adj_b`.
- LarsenResidue.cpp:540 — type alias `SigKey = std::vector<std::uint64_t>` is declared but appears unused in the file (the code uses `std::vector<std::uint64_t>` / `std::map<int,std::uint64_t>` directly) — if dead, drop it; if intended, use it.
- LarsenResidue.cpp:534 — `Signature` tuple is `(Element, int, vector<pair<Element,int>>)`; the `int` is degree and the pair is (neighbour-element, neighbour-prior-hash-folded) — a domain reader can't infer that from the alias; the comment at 536-539 helps but the alias itself could note `// (element, degree, neighbour multiset)`.
- LarsenResidue.cpp:705 — `ambiguous` is well-named; good.
- LarsenResidue.h:130-136 — `N_idx`/`CA_idx`/etc. slot names are clear and domain-correct; no change.

## 3. Visible math structure (grouping)

- LarsenResidue.cpp:573-591 — the WL round body fuses three steps (gather neighbour sigs → sort multiset → hash) into one loop; sound, but the `nbr_sig` build, the `std::sort`, and the `Signature s{...}` would read as a named sequence with blank-line separation — currently the de-facto "symmetrize then hash" structure is implicit.
- LarsenResidue.cpp:97-138 — `FindPeptideAmides` is well-staged (candidates → filter) with the two-phase split visible; good grouping.
- LarsenResidue.cpp:812-894 — `StampChainIdentitiesViaTable` has a clear 4-step shape (parse → build tentative → methyl collapse → table lookup) matching its doc comment; the methyl-collapse inner block (834-850) is the one nested stretch a reader must work through, but the comment at 831-833 labels it. Good.
- LarsenResidue.cpp:684-711 — the cardinality check (685-690) and the within-class assignment (702-710) are two distinct stages run back-to-back; a blank line + 2-word label (`// validate`, `// assign`) would make the "verify isomorphism feasible, then bind" shape pop.

## 4. Function / method naming

- LarsenResidue.cpp:48 — `BondCutoffAngstrom` — name states quantity and units; exemplary.
- LarsenResidue.cpp:169 — `OrderPiecesAlongAmideChain` — says what it returns; good.
- LarsenResidue.cpp:651 — `MatchPiece` returns node indices (not names); the name is fine and the header/comment make the index-not-name contract explicit at 627-645. Good.
- LarsenResidue.cpp:236 — `FinalizeAdjacency` — slightly vague ("finalize" = translate name-keyed to index-keyed); the comment carries it, but `BuildIndexAdjacency` would say it on the tin.
- LarsenResidue.cpp:565/597 — `CanonicalWLSignatures` / `PerceivedWLSignatures` — clear and parallel; good.

## 5. Comments as signposts

- LarsenResidue.cpp:508-532 — the WL "why multi-round" block is excellent grounding (the CD1/CD2 scrambling example is exactly what a chemist needs); keep.
- LarsenResidue.cpp:860-875 — the FATAL `fprintf` message embeds ~15 lines of process/history prose (which test covers it, how to fix, why it's unreachable); this is verbose for a runtime abort string — trim to the actionable line ("canonical chemistry doesn't match generated table; see AllCombinationsPerceiveCleanly") and move the rationale to a code comment above the block.
- LarsenResidue.cpp:989 — `// heuristic; refinement deferred` on `b.order = 1` — good terse signpost.
- LarsenResidue.cpp:99-110 (header) — the `canonical_assignment_ambiguous` comment is long but every sentence earns its place (it defines the entire downstream matching contract); acceptable, though the parenthetical example list could be one line shorter.
- LarsenResidue.cpp:49 — `// Cutoffs match the Python POC.` and similar Python-POC provenance comments (328, 359) are history/process prose; they reassure but don't help a reader understand the math — fine to keep one at file top, redundant per-function.
- LarsenResidue.cpp:1138-1140 — `// Thread-safe per C++11 function-local static init rules.` is correct and useful here (concurrency is non-obvious); keep.

## 6. Correctness (secondary)

- LarsenResidue.cpp:578-583 — folding the prior-round hash to 31 bits (`& 0x7fffffff`) before re-hashing slightly raises collision probability vs. carrying the full 64-bit hash; for ~50 atoms this is almost certainly harmless, but it is a deliberate narrowing of the WL discriminator — confirm intended (the parallel perceived side at 614 does the same, so at least it's consistent).
- LarsenResidue.cpp:112 / 121 — amide detection uses hard distance cutoffs (C=O `< 1.32`, C–N `> 1.20 && < 1.50`); these are convention/geometry-dependent magic numbers — check they hold for the DFT-optimized geometries (r²SCAN/PBE) actually in the DB, not just the Python POC's inputs.
- LarsenResidue.cpp:705-708 — for multi-atom (ambiguous) WL classes the within-class perceived→canonical binding is by container order (`perc_ids[k]`↔`canon_ids[k]`), which is arbitrary; the code documents that `canonical_ambiguous=true` defers to nearest-spatial downstream — correct only if every downstream consumer honors the flag and drops BranchAddress/DiastereotopicIndex before trusting the binding. The contract is stated but enforcement lives outside this file; check the assembler side.
- LarsenResidue.cpp:1241 / 1273 — HIS and non-HIS central paths gate on `piece.size() == canon.atoms.size()` before/inside `EmitPiece` (which re-checks at 907); the HIS path checks size before calling, the non-HIS path relies on EmitPiece's internal check — harmless redundancy, not a bug.
- LarsenResidue.cpp:1014-1039 — `HasAllRequiredSlots` for `Kind::AceCap` returns `atoms.size() == 6` only (no slot check), while NME checks `N_idx/H_idx` plus size; ACE's carbonyl C/O slots are *not* verified despite the comment at 1017-1018 saying "Required slots: C, O" — if a malformed ACE perceived with 6 atoms but mis-stamped roles, this would pass; low risk given upstream EmitPiece, but the check is weaker than its own comment claims.
- LarsenResidue.cpp:970-971 — CB slot is `Element::C && Locant::Beta`; for GLY (no CB) this correctly never fires, and `HasAllRequiredSlots` exempts GLY at 1035 — consistent, no bug.
- LarsenResidue.cpp:1055 — `FindByDftIdx` linear-scans all 5 pieces every call; fine for ~50 atoms, not a correctness issue, noted only because the header advertises it as a lookup.

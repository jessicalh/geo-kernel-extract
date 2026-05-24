# larsen-residue — codex review (readability focus)

- **Targets:** src/LarsenResidue.{h,cpp}
- **Model:** codex `gpt-5.5` `xhigh`
- **Date:** 2026-05-24
- **Tokens:** 33,726
- **Brief:** `../_brief.md`

---

**Verdict**

`src/LarsenResidue.h` tells a mostly coherent API story: typed residue pieces, cached slots, and a single perception entry point. The main readability problem is over-documentation: current contracts are mixed with history, prototype notes, and downstream matcher behavior, so the simple model is harder to see than it should be.

`src/LarsenResidue.cpp` has a real through-line by section order: bond graph → amides → pieces → canonical graph → WL match → emit → perceive. It breaks in the WL/stamping/HIS blocks, where dense symbols and long process comments force the reader to reverse-engineer intent from implementation machinery. I do not see a definite correctness bug from the inlined text.

**1. Coherent Story / Readability**

- `src/LarsenResidue.h:21` — the header opens with long prototype/history prose before the data model; as written, a reader must separate current contract from archaeology — compress to “current perception pipeline” and move Python/prototype history to the design doc.
- `src/LarsenResidue.h:99` — the ambiguity flag comment packs WL classes, graph automorphism, `BranchAddress`, `DiastereotopicIndex`, and matcher policy into one block — reduce to the meaning of the flag plus one sentence on relaxed spatial resolution.
- `src/LarsenResidue.cpp:357` — `ParentHeavyForH` is a long order-dependent if-chain; as written, a reader must scan protein atom-name conventions linearly — group with labels like `// backbone H`, `// beta/gamma H`, `// ring H`, `// polar H`.
- `src/LarsenResidue.cpp:534` — the WL signature model is hard to reconcile with the code: comments describe neighbor degree, but the implementation stores masked prior-round hashes — rename and restate the representation directly.
- `src/LarsenResidue.cpp:1138` — canonical-cache lambdas interrupt the main perception narrative just after piece ordering — add a short `// canonical templates` signpost or move this machinery out of the visible pipeline path.
- `src/LarsenResidue.cpp:1205` — the HIS branch combines policy, variant-order construction, matching, and failure logging in one nested block — split visually with labels: `// variant order`, `// try variants`, `// no match`.

**2. Naming Carries Meaning**

- `src/LarsenResidue.h:129` — `N_idx`, `H_idx`, etc. omit that these are local atom indices for role slots — consider `nitrogen_local_idx`, `amide_h_local_idx`, `carbonyl_c_local_idx`.
- `src/LarsenResidue.cpp:128` — `[c, nj]` and `n_heavy_other` hide the chemistry — use `carbonyl_c`, `amide_n`, `noncarbonyl_heavy_neighbors`.
- `src/LarsenResidue.cpp:172` — `p`, `pi`, and `ai` in chain ordering are premature compression — use `piece`, `piece_idx`, `atom_idx`.
- `src/LarsenResidue.cpp:576` — `cur`, `next`, `nbr_sig`, and `s` make the WL step read like symbol soup — use `previous_round_hashes`, `next_round_hashes`, `neighbor_labels`, `round_signature`.
- `src/LarsenResidue.cpp:1220` — `v`, `tried`, and `any_match` obscure the HIS policy state — use `variant_name`, `attempted_variants`, `matched_variant`.

**3. Visible Math Structure**

- `src/LarsenResidue.cpp:48` — bond cutoffs are embedded literals with only “matches Python POC” as explanation — name the physical thresholds, e.g. `kHydrogenHeavyCutoffAngstrom`, `kSulfurSulfurCutoffAngstrom`.
- `src/LarsenResidue.cpp:103` — amide detection has distinct physical tests but magic numbers `1.32`, `1.20`, `1.50` hide them — use named local constants for carbonyl C=O and peptide C-N windows.
- `src/LarsenResidue.cpp:561` — WL computation should visibly read seed → refine rounds → group classes; currently the masked prior-hash detail appears mid-expression — name `kWlLabelMask` and the prior-label step.
- `src/LarsenResidue.cpp:930` — `EmitPiece` fuses atom payload copying, identity assignment, ambiguity flagging, and slot filling in one loop — add signposts `// atom payload` and `// role slots`.
- `src/LarsenResidue.cpp:1240` — HIS matching fuses atom-count filtering and graph matching — label these as `// atom-count prefilter` and `// graph match`.

**4. Function / Method Naming**

- `src/LarsenResidue.cpp:67` — `Distance` omits units — `DistanceAngstrom`.
- `src/LarsenResidue.cpp:97` — `FindPeptideAmides` returns C-N bond pairs, not amide groups — `FindPeptideAmideBonds`.
- `src/LarsenResidue.cpp:543` — `HashSignature` is too generic inside a file with several graph concepts — `HashWlSignature`.
- `src/LarsenResidue.cpp:651` — `MatchPiece` does a specific WL canonical-node match — `MatchPieceByWLSignature`.
- `src/LarsenResidue.cpp:901` — `EmitPiece` is vague for “populate a LarsenResidue from matched canonical atoms” — `PopulateMatchedPiece` or `EmitMatchedResidue`.

**5. Comments As Signposts**

- `src/LarsenResidue.h:21` — verbose process comment — replace in-code signpost with `// perception pipeline`.
- `src/LarsenResidue.h:99` — dense explanatory block — replace local signpost with `// WL-equivalent atom`.
- `src/LarsenResidue.cpp:1` — file banner is mostly process/history — replace with `// Larsen tripeptide perception`.
- `src/LarsenResidue.cpp:512` — long WL rationale interrupts the implementation — keep only `// WL refinement` here and move the detailed rationale to docs.
- `src/LarsenResidue.cpp:536` — comment appears stale: it says neighbor degree while code uses prior-round hash labels — replace with `// hashed neighbor labels`.
- `src/LarsenResidue.cpp:789` — 20+ lines of stamping prose are too much for the local code — reduce to `// typed identity stamping`, then local labels `// methyl-H collapse` and `// topology lookup`.
- `src/LarsenResidue.cpp:857` — fatal string contains testing/process prose that a code reader must wade through — keep the identifiers and mismatch reason; move test-suite claims out of the runtime message.
- `src/LarsenResidue.cpp:1206` — HIS policy comment is much longer than the branch needs — replace with `// HIS variant policy`.

**6. Correctness**

No definite correctness bug found from the inlined text. The only thing I would check, not call a bug from this review alone, is whether the 31-bit WL hash truncation at `src/LarsenResidue.cpp:583` is intentional and documented; the current comments do not explain the collision tradeoff.

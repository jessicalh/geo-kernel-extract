# tripeptide-pose-assembler — codex review (readability focus)

- **Targets:** src/TripeptidePoseAssembler.{h,cpp}
- **Model:** codex `gpt-5.5` `xhigh`
- **Date:** 2026-05-24
- **Tokens:** 15,501
- **Brief:** `../_brief.md`

---

**Verdict**

`TripeptidePoseAssembler.h` tells the intended validation story, but it does so with too much policy/history prose and one major contradiction around residual rejection. A domain reader can understand the contract, but only after reconciling comments that say residuals are never gates with later comments saying cap-path residuals are gates.

`TripeptidePoseAssembler.cpp` is mostly ordered sensibly: Kabsch setup, cap assembly, typed central assembly, public dispatch. The main readability failure is stale or overlong commentary that obscures the actual algorithm, especially the obsolete central-assembly block and long incident-history comments inside the per-atom matching loop.

**1. Coherent Story / Readability**

`src/TripeptidePoseAssembler.h:37` — says residuals are “NOT a rejection gate” globally, but `:135-149` later says cap-path residuals are excluded and central-path residuals are kept. As written, a reader must discover a hidden path split and mentally override the top-level contract. Suggestion: qualify the first block as “central path residual feature” or remove it and rely on the path-specific contract.

`src/TripeptidePoseAssembler.h:121` — “Aggregate diagnostics (across atoms that passed validation)” conflicts with `n_above_threshold` being incremented for rejected cap atoms and retained central atoms. Suggestion: “Aggregate diagnostics for this assembly attempt.”

`src/TripeptidePoseAssembler.h:123` — `residual > threshold (rejected)` is false for central assembly. Suggestion: “cap rejected; central counted only.”

`src/TripeptidePoseAssembler.cpp:315` — the “Central assembly” comment describes element+nearest-distance matching and sidechain re-rotation, but the actual code below says typed identity, no element heuristic, no sidechain re-rotation. This is the worst story break. As written, a reader must decide whether the code or comment is authoritative. Suggestion: delete this obsolete block or replace with `// typed central assembly`.

`src/TripeptidePoseAssembler.cpp:414` — the aligned-position comment explains historical sidechain re-rotation, ML rationale, and design history in one large paragraph. The actual step is simple. Suggestion: reduce to `// backbone-frame positions` plus one sentence: “Residuals are retained as features; no sidechain-only reorientation.”

`src/TripeptidePoseAssembler.cpp:453` — long dispatch comment mixes WL classes, Markley/CIP convention, examples, and a prior bug. The important logic is strict vs relaxed identity matching. As written, a reader must parse project history before seeing the branch condition. Suggestion: split into `// strict identity` and `// relaxed symmetry class` labels, with examples only if needed.

**2. Naming Carries Meaning**

`src/TripeptidePoseAssembler.cpp:36` — `src` / `dst` are conventional but not domain-specific. Suggestion: `dft_backbone` and `protein_backbone`.

`src/TripeptidePoseAssembler.cpp:37` — `r` hides that it is the alignment result. Suggestion: `fit`.

`src/TripeptidePoseAssembler.cpp:45` — `H` is standard in Kabsch but opaque to non-implementers. Suggestion: `covariance` or `cross_covariance`.

`src/TripeptidePoseAssembler.cpp:54` — `D` is mathematically standard but not self-documenting. Suggestion: `reflection_fix`.

`src/TripeptidePoseAssembler.cpp:437` — lambda parameter `pid` looks like protein id, but it is the perceived/DFT identity. Suggestion: `perceived_identity`.

`src/TripeptidePoseAssembler.cpp:480` — `cand` is compressed. Suggestion: `identity_matches`.

`src/TripeptidePoseAssembler.cpp:507` — `adat` is symbol soup. Suggestion: `aligned_atom`.

`src/TripeptidePoseAssembler.cpp:451` — `perc` is understandable only after context. Suggestion: `perceived_atom`.

**3. Visible Math Structure / Grouping**

`src/TripeptidePoseAssembler.cpp:36-63` — Kabsch computation is compact but readable; the stages could be more visible. Suggestion: add terse labels: `// centroids`, `// cross-covariance`, `// reflection guard`, `// fitted RMSD`.

`src/TripeptidePoseAssembler.cpp:450-525` — central per-atom loop fuses identity dispatch, candidate search, nearest assignment, diagnostics, tensor rotation, and output construction. The computation is coherent but dense. Suggestion: add block labels only: `// identity candidates`, `// nearest pose match`, `// residual diagnostic`, `// rotate tensor`.

`src/TripeptidePoseAssembler.cpp:170-173` and `:521-524` — tensor-frame step is visible enough, but a signpost would help domain readers. Suggestion: `// rotate shielding tensor`.

`src/TripeptidePoseAssembler.cpp:616-625` — residual aggregation is clear.

**4. Function / Method Naming**

`src/TripeptidePoseAssembler.cpp:67` — `ApplyKabsch` is clear enough, but it returns a transformed position, not a generic application. Suggestion: `AlignPosition`.

`src/TripeptidePoseAssembler.cpp:72` — `RotateTensor` is clear but omits frame direction. If convention matters, suggestion: `RotateTensorToProteinFrame`.

`src/TripeptidePoseAssembler.cpp:91` — `SubstrateRoleMatches` is clear.

`src/TripeptidePoseAssembler.cpp:127` — `EmitAlignedAtom` underplays that it may reject the atom on role/residual checks. Suggestion: `ValidateAndEmitCapAtom`.

`src/TripeptidePoseAssembler.cpp:185` — `LarsenLocalToRecIdx` is understandable but compressed. Suggestion: `LarsenLocalAtomToRecordIndex`.

`src/TripeptidePoseAssembler.cpp:335` — `IdentityCompatible` is clear enough.

`src/TripeptidePoseAssembler.cpp:348` — `ProteinIdentityAt` is clear.

`src/TripeptidePoseAssembler.cpp:361` — `AssembleCentralTyped` is clear.

**5. Comments As Signposts**

`src/TripeptidePoseAssembler.h:2-48` — too much history/process prose for a header contract. Suggestion: keep a short contract and move history elsewhere. Terse labels: `// typed identity contract`, `// residual semantics`.

`src/TripeptidePoseAssembler.cpp:50-52` — history/process detail in Kabsch comment. Suggestion: `// reflection guard`.

`src/TripeptidePoseAssembler.cpp:95-97` — comment is useful but wordy. Suggestion: `// semantics unavailable`.

`src/TripeptidePoseAssembler.cpp:217-221` — verbose and slightly duplicates the code path. Suggestion: `// require typed substrate`.

`src/TripeptidePoseAssembler.cpp:230-235` — good substance, too long. Suggestion: `// perceived cap slots`.

`src/TripeptidePoseAssembler.cpp:298-302` — too much invariant history for a loop guard. Suggestion: `// absent optional slot`.

`src/TripeptidePoseAssembler.cpp:528-534` — useful operational note, but too long for inline algorithm reading. Suggestion: `// identity-miss summary`.

`src/TripeptidePoseAssembler.cpp:585-587` — decorative banner adds little. Suggestion: `// public API`.

**6. Correctness**

`src/TripeptidePoseAssembler.h:122-123` — documentation correctness issue: comments say diagnostics are across passed atoms and threshold means rejected, but code counts central over-threshold atoms without rejecting and cap over-threshold atoms before emission. Suggestion: fix comments to match the two path semantics.

`src/TripeptidePoseAssembler.cpp:315-322` — stale comment describes an algorithm the code no longer runs. This is not a computation bug, but it is dangerous documentation because it misstates matching and frame behavior. Suggestion: remove or replace.

No definite computational bug found from the inlined text.

# Fix plan — TripeptidePoseAssembler (Kabsch/SVD backbone superposition + R·σ·Rᵀ tensor frame rotation)

## 1. Summary

`TripeptidePoseAssembler.{h,cpp}` mostly tells a coherent story. The math
core is sound and well-staged: `KabschAlign` runs centroids → centered
columns → cross-covariance `H` → SVD → reflection guard → rotation → fitted
RMSD as a clean linear sequence, and the two consumed quantities (the rotated
tensor `R σ Rᵀ` and the position residual `aligned − protein`) are each a
single labeled line. The validation contract — typed-identity matching across
both substrates, residual-as-ML-feature-not-gate, perception-or-nothing — is
documented at the top of the header and is consistent with how the two
downstream calculators (`TripeptideBackboneShieldingResult`,
`TripeptideNeighborShieldingResult`) and the SDK catalog actually read the
outputs.

Two genuine readability defects exist, both flagged by both reviews:

1. **A stale comment banner** at `cpp:315–322` describes a retired
   element+nearest-distance / sidechain-re-rotation central algorithm that the
   code beneath it no longer runs. It directly contradicts the accurate banner
   immediately below it (`cpp:324–333`). This is the worst story break in the
   file.
2. **History prose burying the algorithm** in the central per-atom match loop:
   a ~25-line essay (`cpp:453–478`) sits inside a ~75-line loop whose actual
   computation (dispatch on ambiguity flag → gather candidates → nearest
   aligned → emit) is ~15 lines.

**This pass will touch:** the stale comment block, over-long history/provenance
comments, a handful of opaque internal local names (`P/Q`, `perc`, `adat`, the
shadowed local `rec`), and a few terse signposts inside the two dense blocks.

**This pass will NOT touch:** any number, sign, or algorithm; the Kabsch
reflection guard (verified earned, see §4); the `R σ Rᵀ` rotation; the
residual sign convention; any output name (NPY/SDK contract names are stable);
struct field names that cross into `ConformationAtom` / catalog / consumers.

---

## 2. Review-finding ledger

Every finding from `codex.md` and `claude.md`, with disposition. (No
`codex-correctness.md` present.)

### codex.md

| # | Finding | Disposition |
|---|---------|-------------|
| C1 | `.h:37` "NOT a rejection gate" globally conflicts with `.h:135–149` cap-path-is-a-gate | **Adopted** → §3 E1. Qualify the top block to point at the path-specific contract. |
| C2 | `.h:121` "across atoms that passed validation" conflicts with counters | **Adopted** → §3 E2. |
| C3 | `.h:123` "residual > threshold (rejected)" false for central path | **Adopted** → §3 E3. |
| C4 | `.cpp:315` stale central-assembly comment (element+nearest + Rsc) contradicts code | **Adopted** → §3 E4 (the central fix). Delete the stale banner. |
| C5 | `.cpp:414` aligned-position comment is one large history paragraph | **Adopted** → §3 E5. Compress. |
| C6 | `.cpp:453` dispatch comment mixes WL/CIP/examples/prior-bug | **Adopted** → §3 E6. Compress to a 2-line signpost. |
| C7 | `.cpp:36` rename `src`/`dst` → `dft_backbone`/`protein_backbone` | **Declined.** These are `KabschAlign(const Vec3 src[3], const Vec3 dst[3])` parameters — `src`/`dst` is the correct general contract of a reusable aligner; the DFT-vs-protein meaning lives at the two call sites (`cpp:266–276`, `398–408`), which is the right place for it. Renaming the generic aligner's params to domain terms would mislabel the abstraction. |
| C8 | `.cpp:37` rename `r` → `fit` | **Adopted** → §3 E7 (rename local `r` to `result` inside `KabschAlign`; `fit` is fine too, `result` matches the local naming used elsewhere). Minor. |
| C9 | `.cpp:45` rename `H` → `covariance`/`cross_covariance` | **Adopted** → §3 E8. `cross_covariance` with the comment retained. |
| C10 | `.cpp:54` rename `D` → `reflection_fix` | **Declined.** `D` is the canonical Kabsch diagonal and is anchored by the adjacent reflection-guard comment and the identical sibling in `RmsdTrackingTrajectoryResult`; keeping `D` preserves cross-file recognisability of the shared idiom. A signpost (E8) carries the meaning without breaking the math-notation match. Weighed call, not cost-driven. |
| C11 | `.cpp:437` lambda param `pid` reads like protein-id but is perceived identity | **Adopted** → §3 E9. Rename `pid` → `perceived_identity`. |
| C12 | `.cpp:480` `cand` → `identity_matches` | **Adopted** → §3 E9 (rename to `candidates`). |
| C13 | `.cpp:507` `adat` → `aligned_atom` | **Adopted** → §3 E9. |
| C14 | `.cpp:451` `perc` → `perceived_atom` | **Adopted** → §3 E9 (`perceived`). |
| C15 | `.cpp:36–63` add terse Kabsch stage labels (`// centroids`, etc.) | **Adopted, restrained** → §3 E10. Add `// centroids`, `// cross-covariance`, `// reflection guard`, `// fitted RMSD` — four 2-word labels, no fragmentation. |
| C16 | `.cpp:450–525` add block labels in central loop | **Adopted, restrained** → §3 E11. `// identity candidates`, `// nearest aligned pose`, `// emit (mirrors EmitAlignedAtom)`. |
| C17 | `.cpp:170–173` / `:521–524` add `// rotate shielding tensor` signpost | **Adopted** → §3 E11/E12. |
| C18 | `.cpp:616–625` residual aggregation clear | **Declined** (no-op; confirmation only). |
| C19 | `.cpp:67` `ApplyKabsch` → `AlignPosition` | **Declined.** `ApplyKabsch` accurately names "apply the stored Kabsch transform to a point"; it is used at both cap and central sites and reads clearly. `AlignPosition` loses the "which transform" specificity. Internal-name churn without readability gain. |
| C20 | `.cpp:72` `RotateTensor` → `RotateTensorToProteinFrame` | **Declined as a rename**, **adopted as a comment.** The function is a pure `R σ Rᵀ` primitive and is frame-agnostic by design; the protein-frame meaning comes from which `R` is passed. A one-line comment at the def (E12) states the frame contract instead, matching claude §6's verified "source→target" analysis. |
| C21 | `.cpp:127` `EmitAlignedAtom` → `ValidateAndEmitCapAtom` | **Declined.** Marginal. The function is cap-only by call structure and the substrate/residual checks are visible in its body; the longer name asserts "cap" in a helper that is conceptually the per-atom emit. Not worth the change; the body already reads as validate-then-emit. (If the human disagrees, it is a safe internal rename — noted in §6.) |
| C22 | `.cpp:185` `LarsenLocalToRecIdx` → `LarsenLocalAtomToRecordIndex` | **Adopted** → §3 E13. Internal file-local static; cheap, clearer. |
| C23 | `.cpp:335/348/361` `IdentityCompatible`/`ProteinIdentityAt`/`AssembleCentralTyped` clear | **Declined** (no-op; confirmation). |
| C24 | `.h:2–48` header history prose too long for a contract | **Adopted, restrained** → §3 E14. Trim the `History:` paragraph (`.h:44–48`) to one durable line + the spec pointer; keep the discipline block. |
| C25 | `.cpp:50–52` history/provenance in Kabsch comment | **Adopted** → §3 E15. Keep the reflection-guard rationale, drop the dated cross-reference line. |
| C26 | `.cpp:95–97` substrate-unavailable comment wordy | **Adopted, light** → §3 E16. Trim, keep the "permit + let residual catch it" fact. |
| C27 | `.cpp:217–221` verbose fail-loud comment | **Adopted, light** → §3 E16. Trim to the durable fact. |
| C28 | `.cpp:230–235` perceived-cap-slots comment too long | **Adopted, light** → §3 E16. |
| C29 | `.cpp:298–302` invariant-history on loop guard | **Adopted, light** → §3 E16. Keep "defensive guard against rec.atoms[-1]"; drop the history. |
| C30 | `.cpp:528–534` CYX/HG note too long | **Declined.** Both reviews note this elsewhere; claude §5 explicitly says it "earns its length" because it explains the dominant production warning at a non-hot site. Keep. (Duplicate of claude L52.) |
| C31 | `.cpp:585–587` decorative `Public API` banner | **Declined.** The box banners are navigation aids that claude §5 recommends keeping; this one marks the public/anonymous boundary. Harmless, consistent with the file's other banners. |
| C32 | `.h:122–123` doc-correctness: counters vs "passed/rejected" wording | **Duplicate** of C2/C3 → §3 E2/E3. |
| C33 | `.cpp:315–322` stale comment is dangerous documentation | **Duplicate** of C4 → §3 E4. |
| C34 | "No definite computational bug found" | **Noted.** Concurs with this plan's §4/§5 verdict. |

### claude.md

| # | Finding | Disposition |
|---|---------|-------------|
| L1 | `.cpp:315–322` stale comment contradicts code; delete it, 324–333 is accurate | **Duplicate** of C4 → §3 E4. |
| L2 | `.cpp:450–526` 15-line algorithm buried in 75; move 453–478 to 2-line signpost | **Duplicate** of C6 → §3 E6/E11. |
| L3 | `.h:44–48`, `.cpp:217–221`, `:233–235`, `:473–478` retired-design history prose | **Duplicate** of C24/C27/C28/C6 → §3 E14/E16/E6. |
| L4 | `.cpp:67–69, 72–74` `ApplyKabsch`/`RotateTensor` crisp — no change | **Declined** (confirmation; supports declining C19). |
| L5 | `.cpp:40` rename `P,Q` → `src_centered`/`dst_centered`, keep `H` | **Adopted** → §3 E8. (Note: this is the opposite of C9 on `H` — see §6 Q1; resolved by keeping `H`→`cross_covariance` since the comment already calls it covariance, and renaming `P,Q`.) |
| L6 | `.cpp:139` local `rec` (`AlignedDftAtom`) shadows the `rec` param meaning used file-wide | **Adopted** → §3 E17. Rename the local to `out_atom`. Real readability hazard (same name, two types). |
| L7 | `.cpp:451, 507` `perc`/`adat` cryptic | **Duplicate** of C13/C14 → §3 E9. |
| L8 | `.cpp:57` `sumSq` clear — no change | **Declined** (confirmation). |
| L9 | `K,R,sigma,H,U,V,D` conventional, comments anchor them | **Adopted as principle** — supports declining C10 (`D`) and keeping `U,V,K`. |
| L10 | `.cpp:36–64` `KabschAlign` well-staged, shape obvious | **Adopted** → §3 E10 (terse labels only, no restructure). |
| L11 | `.cpp:493–498` nearest-spatial scan clean | **Declined** (confirmation). |
| L12 | `.cpp:507–525` central emit duplicates `EmitAlignedAtom` build inline; add a signpost | **Adopted** → §3 E11. One-line signpost `// emit (mirrors EmitAlignedAtom; central path keeps outliers)`. Not a refactor. |
| L13 | `.cpp:185` `LarsenLocalToRecIdx` good | **Note:** conflicts with C22's expand suggestion. **Adopted C22** (slight expansion) — both agree it is clear; the expansion is a marginal improvement, low cost. |
| L14 | `.cpp:335` `IdentityCompatible` reads well | **Declined** (confirmation). |
| L15 | `.cpp:67` `ApplyKabsch` good | **Declined** (confirmation; supports declining C19). |
| L16 | function names uniformly accurate | **Declined** (confirmation). |
| L17 | `.cpp:50–52` reflection-guard signpost good, drop the dated provenance line | **Duplicate** of C25 → §3 E15. |
| L18 | `.cpp:414–425` 12-line "no Rsc" comment accurate but verbose; compress | **Duplicate** of C5 → §3 E5. |
| L19 | `.cpp:528–534` CYX/HG note earns its length — keep | **Adopted** (supports declining C30). |
| L20 | `.cpp:297–302` defensive-guard comment honest — fine | **Partially conflicts** with C29. **Adopted a light trim** (E16): keep the defensive-guard fact, drop the HasAllRequiredSlots history. |
| L21 | `.cpp:77–83, 123–125, 179–181` box banners are good navigation — keep | **Adopted** (supports declining C31; keep banners). |
| L22 | `.cpp:53–56` reflection guard `det(V·Uᵀ)` + `D=diag(1,1,sign)` + `R=VDUᵀ` is canonical — no bug | **Adopted as verification** → §4. |
| L23 | `.cpp:72–74` `R σ Rᵀ` correct for source→target rotation; verified `src`=DFT, `dst`=protein at both call sites | **Adopted as verification** → §4. |
| L24 | `.cpp:113–115` `BackboneHA` slot matches `AlphaHydrogen` role OR `locant==Alpha` — looser than its name; safe given call constraints, latent looseness | **Declined as an edit; logged** → §6 Q3. Per the governing prior this is "I don't see why this is loose" not a traced bug; it is only ever dispatched for `res.HA` (cpp:292), so no α-carbon reaches it. Comment-only clarification optional; flagged as a question, not changed. |
| L25 | `.cpp:503–505` central increments `n_above_threshold` without rejecting — matches contract | **Adopted as verification** → §4. |
| L26 | `.cpp:614` `ok && !aligned_atoms.empty()` semantics correct | **Adopted as verification** → §4. |
| L27 | `.cpp:386, 236, 244` optional-derefs all guarded | **Adopted as verification** → §4. |

---

## 3. Edits that don't move numbers

Comment fixes, signposts, named intermediates, and internal renames. No
output names, no struct fields crossing the catalog/consumer boundary, no
numbers.

- **E1** — `TripeptidePoseAssembler.h:37–42` — Re-point the "residual is NOT a
  rejection gate" block so it does not read as a global claim. Change to state
  that this is the **central-path** semantics, and that the cap path *does*
  excise on residual; cross-reference the per-path contract at `:135–149`.
  (Addresses C1.)
- **E2** — `TripeptidePoseAssembler.h:121` — `// Aggregate diagnostics (across
  atoms that passed validation).` → `// Aggregate diagnostics for this
  assembly attempt (both emitted and rejected atoms contribute to counters).`
  (Addresses C2/C32.)
- **E3** — `TripeptidePoseAssembler.h:123` — `// residual > threshold
  (rejected)` → `// residual > threshold: cap path rejects; central path
  counts only.` (Addresses C3/C32.)
- **E4** — `TripeptidePoseAssembler.cpp:315–322` — **Delete the stale banner**
  (the element+nearest-distance / gotham-Rsc description). The accurate banner
  at `:324–333` remains as the sole description of central assembly. (Addresses
  C4/C33/L1 — the central fix.)
- **E5** — `TripeptidePoseAssembler.cpp:414–425` — Compress the 12-line
  no-Rsc paragraph to: `// Aligned positions via the single BB Kabsch K — no
  sidechain re-rotation. χ1-grid coarseness is left in residual_vec as an ML
  feature (feedback_residual_as_ml_feature) rather than rotated away.`
  (Addresses C5/L18.)
- **E6** — `TripeptidePoseAssembler.cpp:453–478` — Replace the 25-line essay
  with a 2-line signpost: `// Strict identity for chemistry-distinct branches
  (e.g. ILE CG1/CG2); relaxed (drop BranchAddress/DiastereotopicIndex, nearest-
  spatial tiebreak) for graph-automorphic pairs (PHE CD1/CD2, ARG NH1/NH2,
  methyl Hs) that no WL round can split. See larsen-residue-design-2026-05-11.md.`
  (Addresses C6/L2.)
- **E7** — `TripeptidePoseAssembler.cpp:37` — rename local `KabschResult r;` →
  `KabschResult result;` (and its uses through `:63`). (Addresses C8.)
- **E8** — `TripeptidePoseAssembler.cpp:40, 45` — rename `P, Q` →
  `src_centered, dst_centered`; rename `const Mat3 H` → `cross_covariance`
  (comment notes `H = src_centeredᵀ-form covariance`). Keep `U`, `V`, `D` as
  canonical Kabsch notation. (Addresses C9/L5.)
- **E9** — `TripeptidePoseAssembler.cpp:437, 451, 480, 507` — internal local
  renames: lambda param `pid` → `perceived_identity`; loop local `perc` →
  `perceived`; `cand` → `candidates`; `adat` → `aligned_atom`. (Addresses
  C11/C12/C13/C14/L7.)
- **E10** — `TripeptidePoseAssembler.cpp:36–63` — add four terse stage labels
  inside `KabschAlign`: `// centroids`, `// cross-covariance`, `// reflection
  guard`, `// fitted RMSD`. No restructure. (Addresses C15/L10.)
- **E11** — `TripeptidePoseAssembler.cpp:450–526` — add three block labels in
  the central loop: `// identity candidates` (at the dispatch+gather),
  `// nearest aligned pose` (at the best_atom scan), `// emit (mirrors
  EmitAlignedAtom; central path keeps outliers)` (at the AlignedDftAtom build).
  (Addresses C16/L12.)
- **E12** — `TripeptidePoseAssembler.cpp:72` — add a one-line comment on
  `RotateTensor`: `// R σ Rᵀ: rotate a rank-2 Cartesian tensor by R (source→
  target frame; here DFT→protein via K.rotation).` (Addresses C17/C20.)
- **E13** — `TripeptidePoseAssembler.cpp:185` — rename file-local static
  `LarsenLocalToRecIdx` → `LarsenLocalAtomToRecordIndex`. File-local only; no
  cross-file carry-through. (Addresses C22; L13 agrees it is already clear, so
  this is marginal.)
- **E14** — `TripeptidePoseAssembler.h:44–48` — trim the `History:` paragraph
  to one line + the spec pointer: `// The pre-2026-05-11 element-pattern +
  threshold-rejection design is retired; see
  spec/plan/larsen-residue-design-2026-05-11.md.` Keep the discipline block
  above. (Addresses C24/L3.)
- **E15** — `TripeptidePoseAssembler.cpp:50–52` — keep the reflection-guard
  rationale, drop the dated provenance sentence: `// Canonical Kabsch:
  sign(det(V·Uᵀ)) for the reflection guard, not the raw determinant — the
  product is orthogonal (|det|=1), so only the sign matters.` (Addresses
  C25/L17.)
- **E16** — light trims, each keeping the durable fact and dropping history:
  - `cpp:95–97` → `// Substrate not populated (stub fixture): permit; the
    residual gate catches mismappings.`
  - `cpp:217–221` → `// Cross-substrate matching needs typed AtomSemantic;
    without it SubstrateRoleMatches is a no-op, so decline rather than emit
    untyped data.`
  - `cpp:230–235` → `// Cap slots come from typed LarsenResidue perception; no
    heuristic fallback (perception is the single source of DFT-side identity).`
  - `cpp:298–302` → `// Defensive guard against rec.atoms[-1]: slot is -1 only
    if a future invariant relaxation drops a required cap slot.`
  (Addresses C26/C27/C28/C29/L3/L20.)
- **E17** — `TripeptidePoseAssembler.cpp:139` — rename the local
  `AlignedDftAtom rec;` → `out_atom` (and its uses through `:175`) to stop it
  shadowing the file-wide `rec` (`TripeptideDftRecord`). Function-local only.
  (Addresses L6 — a real same-name-two-types hazard.)

---

## 4. Usage notes (the sign/value reasons discovered)

### 4a. Kabsch reflection guard — `sign(det(V·Uᵀ))`, not raw det. COHERENT (expected).

**Reason.** `KabschAlign` (cpp:46–56) builds `R = V · diag(1,1,s) · Uᵀ` with
`s = sign(det(V·Uᵀ))`. The guard exists to prevent the optimal-alignment SVD
from returning an improper rotation (a reflection, `det = −1`) when the point
sets are near-coplanar — three backbone atoms N/CA/C are exactly the
borderline case. Using the *sign* of the determinant rather than the raw value
is correct precisely because `V·Uᵀ` is orthogonal, so its determinant is
±1 up to floating error; the magnitude carries no information and feeding the
raw determinant into the diagonal would corrupt the rotation. This is textbook
Kabsch (Kabsch 1976; the Umeyama 1991 form).

**Where consumed / agreement.** The identical guard, with a fuller
explanatory comment, lives in `src/RmsdTrackingTrajectoryResult.cpp:71–83`
("…the product is theoretically orthogonal with |det| = 1, so the reflection
guard only needs the sign. Passing the raw det…"). Both producers agree. The
rotation `K.rotation` is consumed three ways, all consistent:
- copied to `out.rotation` (cpp:277, 409) → surfaced as `AssembledTripeptide::
  rotation` for diagnostics;
- applied to positions via `ApplyKabsch` (cpp:67, 144, 428);
- applied to the tensor via `RotateTensor` (cpp:72, 171, 522).

**Verdict:** earned convention, do not flip. The only edits are comment
trims (E10/E15). This matches the brief's explicit instruction
("verify, don't flip") and codex C34 / claude L22.

### 4b. Tensor rotation `R σ Rᵀ` — source(DFT)→target(protein) frame. COHERENT (expected).

**Reason.** `RotateTensor(σ, R) = R σ Rᵀ` (cpp:72–74) is the standard
similarity transform for a rank-2 Cartesian tensor under a rotation. Because
`K.rotation` maps the *source* (DFT/Larsen) frame onto the *target* (protein)
frame — verified by the call structure `KabschAlign(src, dst)` with
`src = perceived/cap N/CA/C` and `dst = protein N/CA/C` at both call sites
(cpp:266–276 cap, cpp:398–408 central) — the rotated tensor
`shielding_tensor_aligned` lands in the protein frame, which is what the
downstream calculators need.

**Where consumed / agreement.**
- `TripeptideBackboneShieldingResult.cpp:258` copies
  `shielding_tensor_aligned` → `tripeptide_bb_shielding_tensor` (σ_BB^i,
  catalog `tripeptide_bb_shielding`, "Mat3 ppm").
- `TripeptideNeighborShieldingResult.cpp:305–317` reads the AXA-side and
  AAA-reference tensors **both already rotated by their own Kabsch into the
  protein frame**, then forms the Larsen Eq 3 difference
  `Δσ = σ(AXA) − σ(AAA)` (catalog `tripeptide_neighbor_shielding`). Subtracting
  two same-frame tensors is the correct frame algebra; the reference
  subtraction is the whole point of the neighbor term.
Producer and both consumers agree on the frame. No edit beyond the E12
signpost.

### 4c. Position residual `aligned − protein`. COHERENT (expected).

**Reason / sign.** `residual_vec = aligned_position − conf.PositionAt(atom)`
(cpp:91 doc, cpp:145–146 cap, cpp:512 central). The sign convention "DFT-
aligned minus protein" is the displacement *of the DFT pose from the protein
pose*, captured as a Vec3 ML feature (direction + magnitude), never as a
rejection gate on the central path (`feedback_residual_as_ml_feature`).

**Where consumed / agreement.** The SDK catalog
(`python/nmr_extract/_catalog.py:417`) documents `tripeptide_bb_residual_vec`
as exactly `"aligned_dft - protein position; Vec3 ML feature"`, and the
neighbor variants (`:425`, `:427`) inherit the same sense.
`TripeptideBackboneShieldingResult.cpp:349–351` writes the three components
in (x,y,z) order with no sign flip. Producer and SDK agree byte-for-byte on
sign and component order. No edit.

### 4d. `n_above_threshold` two-path semantics. COHERENT (expected) — comment only.

**Reason.** The same field carries two meanings the header conflated:
- Cap path (`EmitAlignedAtom`, cpp:165–168): residual > threshold ⇒ atom is
  **excluded** *and* counted.
- Central path (`AssembleCentralTyped`, cpp:503–505): residual > threshold ⇒
  atom is **kept** and counted (diagnostic only; residual is the ML feature).
This is intentional (claude L25 confirms; the header `:131–152` even spells
the split out lower down). The defect is purely that the field's one-line
header comment (`:121`, `:123`) reads as if it described a single regime. Fixed
by E2/E3 — comment conforms to code per the brief.

---

## 5. Bug-by-exhaustion candidates

**None.**

Both reviews independently report no computational bug ("No definite
computational bug found" — codex C34; claude §6 verifies the reflection guard,
the tensor rotation frame, the two-path counter, the `ok` semantics, and every
optional-deref as correct). The sign/value tracing in §4 found a reason for
every convention raised, and producer/consumer agree at the catalog and both
consuming calculators. No quantity could not be reconciled, so the exhaustion
bar for declaring a bug is not met anywhere. Every change in §3 is a comment,
signpost, or internal-local rename — zero touch the numbers.

---

## 6. Questions & Ambiguities

- **Q1 (review disagreement on `H`).** codex C9 wants `H` →
  `covariance`/`cross_covariance`; claude L5 explicitly wants `H` *kept* and
  only `P,Q` renamed. Resolved in this plan toward renaming both (E8:
  `src_centered`/`dst_centered` and `cross_covariance`) because the adjacent
  code already computes and comments it as the covariance and `H` is opaque to
  a first-time chemist reader. If the human prefers to preserve the
  Kabsch-textbook `H` for cross-file recognisability against
  `RmsdTrackingTrajectoryResult` (which also uses `H`), keeping `H` is an
  equally defensible call — flagging so the human picks one consistently across
  both files.

- **Q2 (`D` rename, declined — confirm).** I declined codex C10 (`D` →
  `reflection_fix`) to keep the math-notation match with the sibling Kabsch in
  `RmsdTrackingTrajectoryResult.cpp:82` and because the reflection-guard
  comment anchors it. If the human would rather optimise for the
  non-implementer reader over cross-file symbol consistency, `reflection_fix`
  is safe (file-local). Please confirm the preference.

- **Q3 (`BackboneHA` slot looseness, claude L24).** The predicate
  `sem.backbone_role == AlphaHydrogen || sem.locant == Alpha` (cpp:113–115) is
  looser than the slot name: an α-locant non-hydrogen (CA itself) would satisf
  the `locant == Alpha` clause. It is safe today because the slot is only
  dispatched for `res.HA` (cpp:292), so no α-carbon reaches it. Per the
  governing prior this is "I don't see why it's loose," not a traced bug, so I
  am **not** changing the predicate. Open question for the human: is the
  `|| locant == Alpha` clause an intentional tolerance for a known
  HA-role-not-stamped case (in which case a one-word comment would help), or a
  latent over-match that should be tightened to `&& element == Hydrogen`? I
  could not determine which from the substrate-stamping code alone — leaving it
  as a question rather than guessing.

- **Q4 (`EmitAlignedAtom` rename, declined).** codex C21 wants
  `ValidateAndEmitCapAtom`. Declined as marginal, but it is a safe file-local
  rename if the human wants the "cap" + "validate" assertion in the name. No
  consumer carry-through (anonymous-namespace static).

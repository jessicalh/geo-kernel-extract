# Applied-mathematics algorithm reviews

Adversarial reviews of the numerical methods inventoried in
`spec/APPLIED_MATHEMATICS.md`. Each algorithm gets a review per model. Two
models were run — **codex** and **Claude** — each readability-first per its
own tuned brief. The Qwen sweep was **skipped** (2026-05-24): two models were
judged sufficient.

## Layout — retrievable per algorithm, across models

```
reviews/applied-math/
  _brief.md                 shared dispatch brief (all models, all algorithms)
  <algorithm>/
    codex.md                codex gpt-5.5 xhigh — readability-first (canonical)
    codex-correctness.md    codex gpt-5.5 xhigh — earlier correctness-first pass
    claude.md               Claude (agent) — readability
    synthesis.md            reconciliation across models (optional, added last)
  INDEX.md                  this file
```

One directory = one algorithm; one file = one model. To compare an algorithm
across models, read its directory; to compare a model across algorithms, read
that filename down the tree.

## Dispatch convention (keep identical across models)

Same mechanism each time: `brief + line-numbered target file(s)` → the model.
- **codex:** `codex exec - < <prompt>` where the prompt is `_brief.md` followed
  by each target file fenced and `cat -n`-numbered (so findings cite
  `file:line`). gpt-5.5, reasoning `xhigh`, project trusted.
- **Claude:** an agent reads its tuned brief (`_brief-claude.md`) + the named
  source files and writes its review directly to `<algorithm>/claude.md`.

## Status

Legend: codex-R = codex readability pass (canonical), codex-C = codex
correctness pass (earlier, archived alongside).

| Algorithm | Targets | codex-R | codex-C | claude | synthesis |
|---|---|---|---|---|---|
| haigh-mallion | `HaighMallionResult.{h,cpp}` | ✅ | ✅ | ✅ | ⬜ |
| types-spherical-tensor | `Types.{h,cpp}` | ✅ | ✅ | ✅ | ⬜ |
| ring-svd-plane-fit | `Ring.{h,cpp}` | ✅ | ✅ | ✅ | ⬜ |
| eeq-result | `EeqResult.{h,cpp}` | ✅ | ✅ | ✅ | ⬜ |
| planar-geometry | `PlanarGeometryResult.{h,cpp}` | ✅ | ✅ | ✅ | ⬜ |
| biot-savart | `BiotSavartResult.{h,cpp}` | ✅ | ⬜ | ✅ | ⬜ |
| mcconnell | `McConnellResult` + `MopacMcConnellResult` | ✅ | ⬜ | ✅ | ⬜ |
| coulomb | `CoulombResult` + `MopacCoulombResult` | ✅ | ⬜ | ✅ | ⬜ |
| dispersion | `DispersionResult.{h,cpp}` | ✅ | ⬜ | ✅ | ⬜ |
| pi-quadrupole | `PiQuadrupoleResult.{h,cpp}` | ✅ | ⬜ | ✅ | ⬜ |
| ring-susceptibility | `RingSusceptibilityResult.{h,cpp}` | ✅ | ⬜ | ✅ | ⬜ |
| hbond | `HBondResult.{h,cpp}` | ✅ | ⬜ | ✅ | ⬜ |
| sasa | `SasaResult.{h,cpp}` | ✅ | ⬜ | ✅ | ⬜ |
| larsen-residue | `LarsenResidue.{h,cpp}` | ✅ | ⬜ | ✅ | ⬜ |
| tripeptide-pose-assembler | `TripeptidePoseAssembler.{h,cpp}` | ✅ | ⬜ | ✅ | ⬜ |
| larsen-hbond-grid | `LarsenHBondGrid.{h,cpp}` | ✅ | ⬜ | ✅ | ⬜ |
| apbs-field | `ApbsFieldResult.{h,cpp}` | ✅ | ⬜ | ✅ | ⬜ |
| aimnet2-charge-response | `AIMNet2ChargeResponseGradientResult.{h,cpp}` | ✅ | ⬜ | ✅ | ⬜ |
| bonded-energy | `BondedEnergyResult.{h,cpp}` | ✅ | ⬜ | ✅ | ⬜ |

`codex.md` = readability pass (canonical brief); `codex-correctness.md` =
earlier correctness-first pass (5 algorithms only).

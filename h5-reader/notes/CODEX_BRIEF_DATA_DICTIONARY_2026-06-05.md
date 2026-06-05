# Codex — data dictionary DRAFT: every reader metric → meaningful name + explanation + units

## (Standard qt6-cpp + VTK-docs preamble is prepended above. This is a READ-ONLY research+author task — no code, no git. Verify from the code; cite file:line.)

## Why
The reader's metric selection currently keys off raw descriptor ids / storage paths — infra
gibberish to a chemist. The lead wants selection to lead with a **meaningful NAME + one-line
EXPLANATION + UNITS**, with the origin/infra detail pushed to "the part of the list nobody
sees." Step one is the LIST: a human-meaningful dictionary of everything the reader can
surface. You are drafting that list. It becomes the glossary + the selection UI later; the
storage format and UI wiring are a SEPARATE decision — do not wire anything, just produce the
content.

## Read (authoritative list of what the reader surfaces)
- `src/model/TrajectorySignalCatalog.{h,cpp}` — the descriptors: labels, concept keys,
  channels (T0/T1/T2, x/y/z, magnitude, component), display modes, storage paths, the
  dense-H5 / NPY / ORCA-DFT field families. THIS is the authoritative set of metrics.
- `src/model/DisplayModeCapability.h` and the channel helpers — how a metric's channels/modes
  are structured.
- For PHYSICS MEANING you MAY read, READ-ONLY, the calculator exposition docs at
  `../doc/calculators/` (the per-calculator explainers) and the catalog's own comments.
  Do NOT modify anything outside `h5-reader/`, and do NOT read the `nmr_shielding` library
  source (the root `src/`). Ground meanings in those docs + standard NMR/MD domain knowledge.

## Produce — `notes/DATA_DICTIONARY_DRAFT_2026-06-05.md`
A grouped markdown table. Group by physics/data category (ring current, EFG, bond/CSA
anisotropy, electrostatics/APBS, AIMNet2, GROMACS energy, geometry/dihedral, DSSP/topology,
H-bonds/hydration, kernels, ORCA-DFT reference, …). Per metric, columns:
- **Name** — what a chemist would call it (human, meaningful).
- **What it is** — one line, the physical meaning.
- **Units**.
- **Rank/shape** — scalar / vector / rank-2 tensor; mark T0/T1/T2 where it decomposes
  (T2 is the thesis argument — never collapse it; surface it).
- **Origin / infra** (the de-emphasized column): descriptor id, concept key, storage path,
  H5 dataset or NPY field name, producing calculator/source. Cite `TrajectorySignalCatalog`
  `file:line` for each entry.

## Rules
- Do NOT invent confident physics. Where a meaning or unit is uncertain, write the best draft
  and tag it **[VERIFY]** for human/doc refinement. A flagged honest gap beats a confident
  wrong meaning (project ethos).
- This is a DRAFT scaffold — well-drafted helps the lead shape the final wording; she'll
  refine the meanings. Cover EVERYTHING the catalog can surface; note anything you couldn't
  locate a meaning for.
- No code changes, no git, no UI wiring. Output the markdown only. Note the total metric count
  and how many you tagged [VERIFY].

# sasa — claude review (readability focus)

- **Targets:** src/SasaResult.{h,cpp}
- **Model:** Claude (general-purpose agent, opus)
- **Date:** 2026-05-24
- **Brief:** `../_brief-claude.md`

---

## Verdict

**SasaResult.h** — Clean. The header comment lays out the Shrake-Rupley method, the two TOML parameters, and both output NPYs with shape/units/frame, so a chemist sees the whole story before reading the implementation. The class surface is small and the names are honest.

**SasaResult.cpp** — Reads well. The Compute loop follows the textbook through-line a chemist already knows (build lattice → for each atom expand by probe → test occlusion → exposed-fraction × sphere area → average normal), and the two non-obvious choices (search radius, buried-atom normal) are commented. A couple of magic numbers and one slightly costly inner-loop pattern are the only friction. Nothing reads as symbol soup.

## 1. Coherent story / readability
- SasaResult.cpp:69-86 — clean: the point-loop reads exactly as Shrake-Rupley occlusion testing; no decoding needed.
- SasaResult.cpp:94 — occlusion threshold pulled from config inline inside the per-atom loop, breaking the read of the normal block — fine as-is but a reader briefly wonders why a config lookup sits mid-arithmetic; no change needed beyond awareness.

## 2. Naming carries meaning
- SasaResult.cpp:59 — `r_i` is the probe-expanded radius (r_vdW + r_probe), not the vdW radius; the name reads as bare radius — consider `r_expanded_i` or a one-line `// vdW + probe`.
- SasaResult.cpp:63 — `max_vdw` holds a bare Bondi radius while `r_i` already includes the probe, making the `search_radius` sum on line 64 read inconsistently — fine numerically; a `// bare vdW, probe added separately` note would settle it.
- SasaResult.cpp:73-77 — `r_j` here includes the probe (neighbour expanded radius), unlike `max_vdw` on line 63; same `r_*` prefix, two conventions — worth a half-line note.

## 3. Visible math structure (grouping)
- SasaResult.cpp:60 — `sphere_area = 4πr²` computed up top but consumed only at line 89; reads slightly out of sequence but harmless. Clean otherwise.
- SasaResult.cpp:62-65 — the search-radius derivation is well isolated and commented; good grouping.

## 4. Function / method naming
- SasaResult.cpp:15 — `BondiRadius` is a one-line passthrough to `BondiVdwRadius`; the local alias adds an indirection a reader must chase — the name is clear, but the wrapper itself is the only thing to question (not a rename).
- Otherwise clean: `FibonacciSphere`, `Compute`, `WriteFeatures`, `AtomSASA`, `AllSASA` all say what they return.

## 5. Comments as signposts
- SasaResult.cpp:63 — `// largest Bondi radius in table` is good but assumes S is the max; if the table changes this silently under-searches — terse guard note like `// assumes S = max Bondi` would flag the dependency.
- SasaResult.cpp:54, :100 — both GeometryChoice comments restate the code (`// record the parameters used` / `// Record a single GeometryChoice summarising…`) — drop or shorten to `// record params`.
- SasaResult.cpp:91-92, :141 — concise and grounded; keep.
- SasaResult.cpp:123 — `// Rebuild from atoms on demand (rare path…)` is a useful signpost; keep.

## 6. Correctness (secondary)
- SasaResult.cpp:63 — `max_vdw` recomputed every atom iteration from a constant `Element::S`; hoist outside the `i` loop — no result change, avoids N redundant table lookups (minor).
- SasaResult.cpp:73-76 — `BondiRadius(...)` and `r_j*r_j` recomputed for every test point p × neighbour j though they depend only on j; correct, but the per-atom neighbour radii could be precomputed once per atom — flagging as accumulation/cost, not a bug.
- SasaResult.cpp:47 — `sasa_n_points` cast double→int via `static_cast`; truncates rather than rounds (e.g. 91.9 → 91). Check the TOML always stores an exact integer.
- SasaResult.cpp:77 — occlusion uses strict `<` (point inside neighbour expanded sphere). Convention looks standard for Shrake-Rupley; no issue, just confirm boundary points (dist == r_j) should count as exposed.
- SasaResult.cpp:94 — buried-atom guard correctly zeroes the normal; consistent with header line 19. Good.

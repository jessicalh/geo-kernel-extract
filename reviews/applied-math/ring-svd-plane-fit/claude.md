# ring-svd-plane-fit — claude review (readability focus)

- **Targets:** src/Ring.{h,cpp}
- **Model:** Claude (general-purpose agent, opus)
- **Date:** 2026-05-24
- **Brief:** `../_brief-claude.md`

---

**Verdict — Ring.h:** This reads cleanly as a coherent story. The header is a flat, regular type hierarchy where each ring class is a table of named const properties; a chemist can scan it top to bottom and read off intensity, nitrogen count, and aromaticity per residue without decoding. The only friction is that the physical meaning of `JBLobeOffset` and the sign convention of the literature intensities are never stated, but the structure itself is exemplary.

**Verdict — Ring.cpp:** Equally coherent. `ComputeGeometry` is a textbook centroid → SVD-normal → orientation-fix → radius sequence, each stage already signposted; `CreateRing` is a plain dispatch with a well-justified fail-loud tail. A reader following the math once will not get lost. Findings are minor.

---

**1. Coherent story / readability**
- Ring.h — clean. The 8-types-in-3-categories docblock (lines 9–12) sets up the hierarchy before the reader sees it; the through-line holds.
- Ring.cpp — clean. The four geometry stages run in the obvious physical order and each is already labeled.

**2. Naming carries meaning**
- Ring.h:64 — `JBLobeOffset()` — "JB" is unexpanded (Johnson–Bovey?) and the unit/meaning is opaque to a domain expert who doesn't know this codebase — add a 2–4 word comment naming it (e.g. `// Johnson-Bovey lobe z-offset, Å`).
- Ring.h:42 — `total_G_T0_diagnostic` — "G" and "T0" are unglossed here — a one-word trailing note (`// isotropic ring-current kernel trace`) would fix it.
- Ring.h:43 — `mutual_B_from` keyed by `size_t` — key identity (which ring index?) isn't stated — trailing `// keyed by source ring index`.

**3. Visible math structure (grouping)**
- Ring.cpp — clean. Centroid, normal, orientation guard, and radius are visually separated into four blocks; intermediates (`edge01`, `edge02`) are named.

**4. Function / method naming**
- Ring.h:75 / Ring.cpp:8 — `ComputeGeometry` returns a `RingGeometry` (center/normal/radius/vertices) — name is accurate and the return type carries the units story. Clean.
- Ring.h:63,95 — `LiteratureIntensity()` vs `Intensity()` — the distinction (calibrated config value vs published reference) is clear from context. Clean.

**5. Comments as signposts**
- Ring.cpp:65–70 — the fail-loud rationale is 6 lines of process/justification prose where a 2-line signpost suffices — trim to `// Fail loud: an out-of-range enum cast is a programmer error, not recoverable.`
- Ring.h:178–185 — the ProPyrrolidine docblock is long but earns it (justifies the literal 0.0 as physics not calibration); could shed the "Joule & Mills 2010 ch. 7" sentence to a trailing `// see Joule & Mills 2010 ch.7` but acceptable as-is.
- Ring.cpp:23–25 — SVD-normal comment is correct and well-placed. Clean.
- Ring.cpp:34–35 — "right-hand rule" orientation note is exactly the right signpost. Clean.

**6. Correctness (secondary)**
- Ring.cpp:32 — `svd.matrixV().col(2)` hardcodes column index 2 as the smallest singular value — correct only because `coords` is N×3 so V is always 3×3; safe here, but a `// col 2 = smallest σ (V is 3×3)` would document the assumption.
- Ring.cpp:38 — orientation guard uses vertices 0,1,2; if the first three ring atoms were ever collinear `edge01.cross(edge02)` ≈ 0 and the sign flip is undefined — for real aromatic/saturated rings this cannot happen, so not a bug, but check that callers never pass a degenerate 3-atom set.
- Ring.cpp:21,46 — centroid and radius divide by `vertices.size()` with no zero guard, but the `atom_indices.empty()` early return at line 10 covers it. Safe.
- Ring.cpp:15 — `positions[idx]` is unchecked; relies on caller supplying indices in range — consistent with the project's fail-loud-on-programmer-error stance, no in-loop guard needed.

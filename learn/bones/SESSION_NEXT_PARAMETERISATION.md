# Session prompt: CalculationArea — integration and Ring.h refactor

## What happened last session (2026-04-06 afternoon)

Wrote CalculationArea.h and CalculationArea.cpp. These are COMPLETE —
14 concrete types covering all 39 choice points in the catalogue.
Committed as 54cdd41 on master.

Also wrote a Ring.h refactor PROPOSAL (not implemented) at
spec/RING_OBJECT_MODEL_PROPOSAL.md. Do NOT implement it without
discussing with the user first.

The full session record with design reasoning and quoted conversation
is at spec/PARAMETERISATION_CAMPAIGN.md. Read sections 1-7 for context.
Sections 8-9 (next steps) are partially stale — the header is done,
ignore "Phase 1: complete the header."

## What is DONE (do not redo)

1. CalculationArea ABC with Name(), Description(), Domain()
2. DomainKind enum (13 values)
3. All 14 concrete types with class definitions AND implementations:
   RadialThreshold, ShellBoundary, SourceRelativeExclusion,
   RingBondedExclusion, SelfSourceExclusion, SequenceGate,
   SwitchingFunction, DecayFunction, RingCurrent, LobeOffset,
   NumericalAccuracy, ValueClamp, ValueGate, SentinelValue,
   WholeDomain
4. Physics comments on every type (the comments ARE the deliverable —
   they document provenance and justification for every choice point)
5. Catalogue CSV at spec/CALCULATION_AREA_CATALOGUE.csv (39 rows)

## What is NOT done

1. No calculator includes CalculationArea.h. It is not wired in.
2. No CMakeLists.txt change. The files will not compile in the
   current build — that is expected, they are not integrated yet.
3. No named instances created. The types exist but nobody constructs
   a RingHorizon or MultipoleInnerBoundary object yet.
4. Ring.h refactor is a PROPOSAL only. The existing 11-class hierarchy
   still exists and all calculators use it.
5. No binary comparison test (extract before/after, diff .npy files).
6. How CalculationArea objects enter the pipeline is an OPEN design
   question — the user has not decided this yet. Do not assume a
   factory function or any specific integration pattern.

## Stale documents (read with caution)

- spec/SESSION_NEXT_CALCULATION_AREAS.md — this was the PREVIOUS
  session's briefing. It says "write the first 3 concrete types."
  We wrote 14. It lists the reading list (still valid) and design
  constraints (still valid) but its deliverables are superseded.

- spec/PARAMETERISATION_CAMPAIGN.md sections 8-9 — Phase 1 is done.
  Phase 3 was corrected to "TBD." The rest of the document (sections
  1-7, 10-11) is accurate and contains the design reasoning.

## What to read (if you need physics context)

Only read these if the session requires understanding the calculators.
Do not re-read them just because they are listed — check what the
user is asking for first.

1. GEOMETRIC_KERNEL_CATALOGUE.md — kernel physics
2. OBJECT_MODEL.md (first 250 lines) — ConformationResult pattern
3. src/KernelEvaluationFilter.h — existing filter framework
4. src/BiotSavartResult.cpp — example calculator using 3 area patterns

## Design rules (settled, do not re-derive)

1. ABC is minimal: Name(), Description(), Domain(). Nothing else.
2. Type = pattern, not calculator.
3. At most 1 abstract layer between ABC and concrete types.
4. Areas are CONSTANTS — they do not Accept()/Reject(). Filters evaluate.
5. Ring parameters are in the hierarchy ("rings need to count").
6. Any Ring.h refactor must produce binary identical extraction output.
7. Physics comments on each type are the primary deliverable — this is
   about making 15k lines of working code readable as science.
8. The user designed the object model. Do not reorganise it. Ask first.

## Known design tensions (discussed, not resolved)

1. Ring identity fragments across CalculationArea types (RingHorizon
   is a RadialThreshold, PheRingCurrent is a RingCurrent, etc.).
   This is the right tradeoff for indexing by pattern. Ring.h stays
   the home for "what IS a ring."

2. Ring current intensity and lobe offset are coupled (change d, need
   different I). "A ball of worms but one we need to open just not
   today."

3. Whether intensity stays as a literature constant, becomes learnable,
   or gets replaced by a calculated estimate is an OPEN question.

4. The Dispersion switching function (onset 4.3A, cutoff 5.0A) was
   merged into one SwitchingFunction instance. The catalogue had them
   as two separate entries.

## What the user likely wants next

Ask them. Possibilities include:
- Wire CalculationArea into the build (CMakeLists, compile test)
- Discuss how areas get instantiated and consumed
- Ring.h refactor (constexpr table — spec/RING_OBJECT_MODEL_PROPOSAL.md)
- Binary comparison test infrastructure
- Continue with other project work entirely

Do NOT start coding until you have talked to the user about what
they want this session.

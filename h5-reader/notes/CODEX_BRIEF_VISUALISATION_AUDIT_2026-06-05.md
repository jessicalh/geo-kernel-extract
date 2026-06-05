# Codex — visualisation class-system audit + must-have proposal

## (Standard qt6-cpp + VTK-docs preamble is prepended above. This is a READ-ONLY audit+proposal — no code, no git. Verify from the code; cite file:line.)

## Why
The lead's read of the reader: "We only show strips right now. None of the other visualisation
tools exist. And yet they are spattered all over the UI as options, and I doubt there is even a
class system to organise them." She wants to **decide what we MUST have (strips + 0–3 more),
support those, and remove all the hollow bits.** Before any cutting, ground three things:
what actually renders, whether there's a real class system, and a clean proposal for her to
decide against. ("Understand before cutting" is the lead's rule.)

## Read (verify from code; cross-ref but don't trust the prior inventory)
- `src/model/DisplayModeCapability.h` (the capability table), `src/app/DashboardDisplayController.{h,cpp}`
  (mode→panel dispatch), `src/app/StripStackWidget.*`, the `*Panel` classes, the scene/reveal
  overlays (`SceneRevealOverlay`, glyph/tensor overlays), `src/app/SignalDisplayDialog.*`
  (where modes are offered), the catalog's `displayModeIds`.
- `notes/DISPLAY_MODE_INVENTORY_2026-06-05.md` — the prior inventory (which modes render).
  Use it as a map, but re-verify the class-structure question from the code.

## Answer (each, with file:line)
1. **Class system — is there one?** Is there an actual typed hierarchy organising visualisation
   types (a base "visualisation/panel" with typed subclasses answering for themselves), or is
   it string mode-ids + switch/if dispatch + ad-hoc panels? Map how a metric+mode actually
   becomes a rendered thing today. Judge it against the project's OO discipline (typed objects,
   virtual methods, no string dispatch on identity).
2. **What actually renders.** Confirm exactly what draws pixels: strips (which strip modes),
   which `static.*` panels, which scene overlays. Everything else is hollow/advertised-dead.
3. **Clean class-system proposal.** Sketch a minimal typed visualisation hierarchy
   (capability-driven, no string dispatch) that organises the must-haves AND makes hollow modes
   structurally un-offerable — you can't advertise a visualisation that has no class behind it.
   Keep it small; this is a desk-ready basic screen, not a framework.
4. **Must-have proposal (for the lead to DECIDE, not decide for her).** The reader's job is
   vetting H5 fields before they become thesis claims — surfacing complex data (per-metric
   colouring, tensor glyphs preserving T2, strip charts, maybe one structural/sequence view).
   Propose strips (have) + **0–3** candidate must-haves, each with: what it shows, why it's
   must-have for vetting, what data it needs, rough build cost (reuse vs new). Keep T2/tensor
   visible — the rank-2 tensor is the thesis.
5. **Remove list.** Every hollow/dead mode + option to cut, and WHERE each is "spattered"
   (which dialog/menu/toolbar offers it) so the removal is mechanical later.

## Output — `notes/VISUALISATION_AUDIT_AND_MUSTHAVES_2026-06-05.md`
Read-only proposal. No code, no git. Cite `file:line` throughout. End with a crisp
"DECISIONS OWED TO THE LEAD" section (must-have set; keep/cut boundary; class-system go/no-go).

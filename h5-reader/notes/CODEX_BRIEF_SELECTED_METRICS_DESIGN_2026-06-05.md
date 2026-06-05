# Codex — selected-metrics tracking + panel-visibility DESIGN review (READ-ONLY)

Working root `/shared/2026Thesis/nmr-shielding/h5-reader` (branch h5-reader-pysr-spike). **READ-ONLY: map +
design. Do NOT edit, build, git, or launch.**

## FIRST — load the qt6-cpp skill, and respect the portable project
Read the qt6-cpp skill before reviewing: `/home/jessica/.claude/skills/qt6-cpp/SKILL.md` +
`references/model-view.md` + `references/architecture.md`. The Qt model/view + signal-slot idioms are the lens
for the selected-metrics design. **This is a PORTABLE Qt6 project (Windows / macOS / Linux)** — treat the
Windows and macOS source paths as first-class; do NOT disparage, skip, or "hold your nose" at the cross-platform
code. The design must hold on all three.

## The problem (the lead — sharper than a bug)
The reader "mostly works but is too funky to hand off," and the metric/panel area is the worst of it: it feels
like something is **failing silently and the UI state is not reflecting reality.** Concretely: the lead selects
a field, goes to the panel selector/button, and sees **no panels enabled — the selector is "constantly greyed
out"** instead of reflecting what's selected. (A recent fix that greys non-rendering mode checkboxes made this
feel worse.) This is not "find the one missing call" — it's "is the underlying model even coherent."

## What to figure out
1. **Is there a clear, SINGLE internal list of selected metrics — a source of truth?** Where does it live
   (`DashboardSignalModel`? `DashboardPanelModel`? the strip tracks? scattered across them)? Is it one coherent
   list or duplicated/derived in several places? Map it precisely.
2. **What events should mutate that list, and do they?** Lead's model: adding a measure ADDS to it; closing a
   panel SHRINKS it. Trace add/remove end-to-end — is the list the authority, or does state diverge?
3. **How should panel/dock VISIBILITY relate to the list?** Lead's model: the panel is UP by default whenever
   there are selected metrics, hidden ONLY if the user explicitly disabled it. Map the current visibility logic
   (docks start hidden; the toggle path; the dead-mode greying) and exactly where it diverges from this.
4. **Why is the panel selector "constantly greyed out"?** Trace the enable/disable state of the panel selector,
   the metric-mode checkboxes, and the dock toggles. Pin where the DISPLAYED state diverges from the actual
   selected-metrics list — the silent failure.
5. **Are the relevant internals good?** Assess `DashboardSignalModel` / `DashboardPanelModel` /
   `DashboardDisplayController` / `DashboardStripDock` / `StripStackWidget` / `SignalDisplayDialog` as a model
   layer: single source of truth? complete signals? coherent state? Flag the rot, in qt model/view terms.

## The design target (propose it, grounded in the qt model/view idiom)
ONE internal list of selected metrics is the source of truth; add/remove events mutate it; panel/dock
visibility DERIVES from it (visible when non-empty unless the user disabled it); and the selector + checkboxes
REFLECT the list state accurately — never greyed-out-disconnected-from-reality. Propose the model, the events,
the visibility rule, and the internals fixes. State plainly whether the current code can be steered to this or
needs a model rebuild.

## ALSO — relate it to DATA AVAILABILITY and the existing UI DATA-FLOW model (the bar is "well done")
The selected-metrics design does NOT stand alone. It must compose with two existing systems — center these:

**A. Data availability (`TrajectoryFieldAvailability`).** A loaded run already classifies each field
{Absent / NoFramePayload / AllMissing / AllZeroStructural / AllZeroObserved / Available}. Availability is the
GATE on what is even selectable. Map how (or whether) the metric-selection path uses availability today.
Critically: is the current dead-mode greying **conflating "this display MODE has no renderer" with "this DATA
is unavailable"**? Those are different axes; untangle them. Design the relationship cleanly — availability
gates what the picker OFFERS; the selected-metrics list is the user's chosen subset of AVAILABLE metrics; the
selector's enabled/greyed state derives from availability + selection, never from an ad-hoc test disconnected
from "is there data."

**B. The existing UI data-flow model (model-is-spine).** H5 / snapshot → the typed model (QtProtein /
Conformation / QtFrame) → views (inspector / dashboard / scene); the dashboard signal→panel→strip flow is a
branch of it. Map where the selected-metrics state and the panel rendering sit in that flow. The design must
FIT model-is-spine: the single selected-metrics list lives in the model layer, views DERIVE from it, values
come through the existing snapshot/time-series path — do NOT bolt a parallel state machine onto the widgets.

**The unified picture to propose:** availability (what data exists) + selection (what the user chose to show)
as related MODEL state; panel/dock visibility DERIVES from selection (non-empty → up unless disabled); views +
selector REFLECT both — the UI cannot display state disconnected from reality. Show how the three
(availability, selection, data-flow) compose into ONE coherent model. This is the bar.

## Output
Write `notes/SELECTED_METRICS_DESIGN_2026-06-05.md`: (1) the CURRENT tracking mapped — the list (or its
absence) and the events, cited `file:line`; (2) the silent-failure / state-divergence root cause; (3) the
PROPOSED clean model (list + events + visibility + selector-reflects-reality) in qt model/view terms; (4) the
internals verdict + a fix plan (steerable vs rebuild). READ-ONLY — design, do not implement.

## Boundaries
Read-only; qt skill loaded first; portable project (Windows/macOS source is first-class); `h5-reader/src`
scope; cite `file:line`; concrete, not "it's complex."

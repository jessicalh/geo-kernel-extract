# Planning Artifacts

This directory holds bounded migration plans, architecture-layout notes,
and matrix-style working documents. These are execution artifacts used
to organise multi-session work and round-robin review.

Use this directory for:

- calculator/topology migration matrices
- staged implementation plans
- round-robin review records
- known untestable-window notes
- rollback-point records
- per-stage acceptance criteria

Keep root-level `spec/` documents for stable project contracts. Planning
documents may graduate into root-level specs only after the design has
landed and the plan no longer needs to change.

Current planning packet:

- `spec/plan/current-topology-anchor-2026-04-29.md`
- `spec/plan/amber-terminal-charge-generation-2026-04-29.md`
- `spec/plan/legacy-amber-implementation-brief-2026-04-29.md`
- `spec/plan/pre-iupac-cruft-map-2026-04-29.md`
- **`spec/plan/amber-implementation-plan-2026-04-29.md`** —
  TODAY'S CENTRAL CAPTURE-OF-DECISIONS. Six implementation steps
  GREEN in the working tree (62/62 tests). Locks O2/O3/O4,
  the crystal projection rule, the substrate-vs-conformation
  split, the drift policy, and the post-slice sequencing
  (PHASE 0 → PHASE 1 (N1.A-G) → PHASE 2 → OpenBabel exit → N4).
- `spec/plan/openai-5.5-strong-architecture-layout.md`
- `spec/plan/legacy-amber-topology-matrix-prespec-2026-04-28.md`
  (matrix template; superseded as a deliverable by the central
  plan above, but the row-template structure remains useful)

Current implementation entry point:

- `spec/plan/current-topology-anchor-2026-04-29.md` (anchor)
- `spec/plan/amber-implementation-plan-2026-04-29.md` (central plan)

Latest OpenAI response:

- `spec/plan/openai-round-4-response-2026-04-28.md`

Review records in this directory, including `opus-round-*.md`, are
historical agent-to-agent responses. They may contain superseded names or
proposals. Treat them as review history, not as active spec text.

Current baseline:

```text
The planning packet is the baseline for the next round-robin review.
Prior half-finished reviewer edits are not active instructions.
New sessions should read the packet from the top, review the model, and
propose concrete edits or questions against the current files.
```

Current topology anchor:

```text
ProteinTopology is the small ABC.
LegacyAmberTopology is the concrete current calculator topology.
CovalentTopology is the bond graph component owned by LegacyAmberTopology.
Protein::BondTopology() is a compatibility/read path to that component.
Protein::Topology(), if present, is only an old covalent-topology synonym.
```

Round-robin rule:

```text
Each stage gets a fresh review by the active Opus sessions, OpenAI 5.5,
and the user. Agreement per stage precedes implementation.
```

Review protocol:

```text
Respond to the current packet, not to stale session momentum.
Keep claims tied to file sections, object names, or matrix rows.
Ask questions when project intent changes the model.
Do not implement source changes during planning review.
Do not replace the ABC topology model with an untyped mapping model.
```

# Brief: the object model and runtime architecture

You're writing the document that explains how this extractor is built: how it holds a
protein, how it manages its buffers, what calls what during a run, and why the object
model is shaped the way it is. It is the outward-facing companion to the internal
`OBJECT_MODEL.md` — the version you hand someone who needs to understand the
architecture without living in the source.

Take your time and do this well. This is the heart of the engineering argument, and it
goes to readers who will judge the design.

## Who reads this

Someone with a year of calculus-based physics and linear algebra (Gilbert Strang on
the shelf) and twenty years of C and C++. They know ownership, lifetime, RAII, virtual
dispatch, value-versus-reference semantics, cache behaviour, and call graphs cold — you
do not explain those. What you explain is *this system's* choices: why a protein is
split into identity/topology and per-conformation geometry, how the buffers are owned
and reused, what the per-frame call sequence is, and why the named-operation pattern
looks the way it does. Assume a reader who will push back on a hand-wave; give them the
real structure, named and precise.

## The standard

Real grounding, no fake-plain filler, no superlatives or absolutes. Name the real
types and methods — they are the subject of this document, not jargon to hide. When you
introduce a type, say what it owns and what question it answers about itself. Where you
describe a call sequence, it must match the code's actual order. An analogy only where
it maps the architecture exactly (the project's own "named operations are rooms, not
wrappers" is one such picture — use it only if the code bears it out), never a
decorative one.

## The code is the final authority; the design docs are intent

Read `OBJECT_MODEL.md` and `PATTERNS.md` (repo root) and `spec/CONSTITUTION.md` in
full — they are the design intent and the reasoning, they are long, and you read them
rather than skim. Then read the source, which is the final authority:
`src/Protein.{h,cpp}`, `src/ProteinConformation.{h,cpp}`, `src/ProteinTopology.h`,
`src/ConformationResult.{h,cpp}`, `src/DenseBuffer.h`, `src/OperationRunner.{h,cpp}`,
`src/OperationLog.{h,cpp}`, `src/Session.{h,cpp}`, `src/Trajectory.{h,cpp}`,
`src/TrajectoryProtein.{h,cpp}`, `src/RunConfiguration.{h,cpp}`. Identify the real
per-conformation and per-frame geometry types from the headers rather than assuming a
name. Where the design docs and the code disagree, follow the code and note the
divergence — the docs are large and may have drifted.

## Load-bearing invariants to get right (the reader will check these)

- A `Protein` is identity and topology only. Per-atom geometry lives on the
  conformation; per-frame geometry on the frame type. State this and show where the
  boundary is enforced.
- The tensor output (the T0/T1/T2 rank-2 decomposition) is preserved end to end; the
  data model does not collapse it to a scalar.
- Objects answer questions about themselves through typed methods rather than string
  dispatch on identity.

Make these true by pointing at the code, not by assertion.

## What the document covers

1. The shape of the system in one page: the few primary objects and how a run flows
   through them, from input reader to emitted output.
2. The primary object model: `Protein` (identity + topology), the per-conformation
   geometry object, the per-frame geometry object, and the result object that carries
   computed kernels. What each owns, what each answers, where the geometry boundary
   sits.
3. Buffer management: how `DenseBuffer` and the per-conformation/per-frame buffers are
   allocated (the buffers-from-the-constructor choice), owned, reused across frames,
   and read back. The lifetime story, in real terms.
4. The swimlanes — what is called when: session/run setup, the calculator/operation
   dispatch (`OperationRunner`, `OperationLog`), and the per-frame loop, with the
   actual call sequence and how short the inner loop is. A diagram is welcome if it
   matches the code.
5. The complete protein object model, documented: the type hierarchy and relationships
   (atoms, residues, rings, topology, conformation, frame, result), enough that a
   reader could navigate the source from this document.

## Method

Write the draft, then re-read it against the source: does every ownership claim, every
call-order claim, every type relationship match the code? Fix the drift. Then the copy
pass — no superlatives, no smeared ideas, no undefined term, no claim the code does not
bear.

## LaTeX and PDF

`article`, 11pt, `geometry`, `graphicx`, `amsmath` if needed, `\texttt` for type and
method names, a `listings`/`verbatim` block for any code or call-sequence sketch, and a
`\tableofcontents` (the document is long enough to want one). Expect this document to run
well beyond five pages — the calculator notes were short by design, this one is not. Do
not compress the object model to fit a few pages; completeness and clarity are the goal,
and length follows the material. Do not pad either: every section carries its weight.
Compile with `pdflatex` run twice (for the table of contents) to a clean build.

## Diagrams — use Mermaid

This document should carry diagrams, and they earn their place: a sequence/swimlane
diagram of the per-frame call flow (who calls what, in order) and an object-model
relationship diagram (Protein / topology / conformation / frame / result, with the
ownership edges). Add more where they clarify. The diagrams are part of the verification
surface — each must match the code's actual structure, not an idealized one.

Render Mermaid with the project's known-good pipeline (do not fight puppeteer); work from
`doc/calculators/`:
1. Write each diagram as a `.mmd` file.
2. Create `doc/calculators/.puppeteerrc.json` containing
   `{"args": ["--no-sandbox", "--disable-setuid-sandbox"]}` (required on this Ubuntu).
3. Render to vector PDF: `mmdc -i NAME.mmd -o NAME.pdf -p .puppeteerrc.json --pdfFit`
   (`--pdfFit` crops the page to the diagram). If a diagram gives mmdc trouble as PDF,
   fall back to `mmdc -i NAME.mmd -o NAME.png -p .puppeteerrc.json -s 3`.
4. Include with `\includegraphics[width=0.75\textwidth]{NAME.pdf}` inside a `center`.

Keep each diagram to roughly 5–8 nodes so it stays legible; split a dense flow into
several focused diagrams rather than one tangle. Mermaid label text avoids the LaTeX
Unicode problem (it is rendered to an image), but keep node labels to real type and
method names. Confirm the rendered files exist and the final `pdflatex` build embeds them
with no missing-file or undefined-reference errors.

## Done when

The primary object model, the buffer-management story, the call-sequence swimlanes, and
the complete protein model are each grounded in the source; the load-bearing invariants
are shown true against the code; it reads clearly for the Strang-plus-twenty-years-C++
reader; it compiles clean.

# doc-discrepancies/

Bottom-up doc audit for the truthfulness campaign (see memory
`project_docs_single_source_of_truth`). While fixing `src/` comments (or any
code work), confirm code + comment against the now-truthful doc set and log
every **code-vs-doc mismatch** here as a residual DOC error.

Convention (one shared file races — never use one):
- Each agent writes its OWN file: `<UTC-timestamp>-<short-slug>-disc.md`
  (e.g. `20260526T214530Z-mcconnell-disc.md`). Never append to another file.
- One entry per line: `<doc>:<line> | <doc claim> | <code truth> | <src file:line> | DOC|COMMENT`
- A later aggregation pass dedupes all `*-disc.md` into one final doc-cleanup
  sweep, then clears them.

`*-disc.md` are gitignored (transient); this README is the durable convention.

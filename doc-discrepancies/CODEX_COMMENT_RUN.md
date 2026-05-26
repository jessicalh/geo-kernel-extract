# Autonomous comment-truthfulness run (codex-owned)

You (codex) own the `src/` comment-truthfulness pass from here. Run it to
completion, in chunks, committing as you go. This is COMMENT-ONLY work on a
FROZEN library — you change comments, never code. Work the repo at
`/shared/2026Thesis/nmr-shielding`.

## Read first
`doc-discrepancies/COMMENT_ARC.md` — the full doctrine, the 16-group arc, and the
per-group file membership + progress checkboxes. That file is authoritative; the
summary below is the load-bearing subset. Also skim a few already-finished files
to calibrate (e.g. `git show` the diffs of commits whose subject starts
`docs(src comments):` — Protein.h, ConformationResult.h, errors.h, GeometryChoice.h).

## The priority (why you are doing this)
Most comments get retired into a future single source of truth, so this is not
precious polishing. Two moves:
- **Cut aggressively:** banners (`// ====`, `// ----`), section/group labels,
  type/hierarchy maps, restatements of adjacent code, derivable usage examples,
  AI/marketing filler, ALL changelog (dates, commit hashes, "was X / now Y",
  "Bundle/Slice/phase" labels, "replaces", "future consumers", "reserved for
  future"). When in doubt, cut.
- **Preserve carefully — the irreplaceable core:** the non-obvious *why*, the
  gotcha, the precondition/contract, the unit, the citation — anything that
  exists ONLY in the comment and a reader could NOT reconstruct from the code.
  These are the point. Verify them against the code; do not destroy them.

## The keep/cut test (master rule)
Keep a comment ONLY if a careful reader doing a deep-dive of the code could not
reconstruct it from the code itself. If they'd figure it out anyway → cut.

## Hard rules
1. **COMMENT-ONLY. Never change code** — not a statement, signature, initializer,
   include, whitespace, or a lint "fix". After editing each file, run the
   self-check below; if it reports CODE CHANGED, `git checkout -- <file>` and redo
   comment-only.
2. **SOURCE IS TRUTH.** Verify each comment against the `.cpp` implementation
   (header comments lie too). A claim that is completely out of date may be
   *removed* — you need not re-ground it.
3. **Serve THIS code.** No domain textbook, no pipeline/architecture tour ("X
   does this, then Y's job…"). Recast such story into a crisp, verified
   precondition or statement of this code's own effect. Do not assume a doc holds
   a real constraint — keep the constraint, phrased locally.
4. **Citations/factoids are claims.** Verify. If you cannot verify a citation
   (can't confirm the source/value), preface the comment `From notes(check):` —
   never present an unverified reference as authoritative. Canonical, on-topic
   references (e.g. Johnson & Bovey 1958, Caldeweyher 2019) may stay as-is.
5. **Finish each file in one pass.** No half-done files; if finishing needs you to
   read a calculator / TrajectoryAtom / the .cpp, go read it now.
6. **Don't push. Don't touch** `src/generated/` (generated), or anything outside
   `src/`.

## Method per file
read the file + its implementation → edit comments to the rules above →
self-check (below) → commit with a clear message (`docs(src comments): <file or
group> — <what changed>`) → tick the file/group in COMMENT_ARC.md. Commit per file
or per small coherent group. Re-read the diff critically before committing — you
are your own reviewer now.

## Comment-only self-check (run before every commit)
```bash
# usage: pass the files you edited
for f in "$@"; do
  if diff <(git show "HEAD:$f" 2>/dev/null | sed 's|//.*||; s/[[:space:]]*$//' | grep -v '^[[:space:]]*$') \
          <(sed 's|//.*||; s/[[:space:]]*$//' "$f" | grep -v '^[[:space:]]*$') >/dev/null; then
    echo "OK comment-only: $f"
  else
    echo "CODE CHANGED — reverting: $f"; git checkout -- "$f"
  fi
done
```
(Strips `//` line comments + trailing whitespace + blank lines, then diffs the
code tokens vs HEAD. Empty diff = comment-only. The codebase is `//`-comment
heavy; if you touch a `/* */` block, reason about comment-only manually too.)

## Order / resume
Work the arc in `COMMENT_ARC.md` top to bottom. **Start where it is unfinished:**
G1's last two files — `src/Protein.cpp` (full pass; only Protein.h is done) and
`src/ConformationAtom.h` per-calculator field blocks (the `// === X (set by
YResult) ===` sections + their field comments; the header is done) — then G2
(Types, SemanticEnums, AminoAcidType, PhysicalConstants — the citation rule's real
test), then G3…G16. Enumerate a group's files with `find src -maxdepth 2 \(
-name '*.h' -o -name '*.cpp' \)` against the group's membership rule. To see
what's already done: `git log --oneline | grep 'docs(src comments)'`. Tick boxes
in COMMENT_ARC.md as you complete files/groups, so the next run resumes cleanly.

## Log doc discrepancies (optional, cheap)
If a comment reveals that a `.md` doc (OBJECT_MODEL, ARCHITECTURE, CATALOGUE,
API.md, …) is wrong about the code, append one line to
`doc-discrepancies/<UTC-timestamp>-codex-disc.md`:
`doc:line | doc claim | code truth | src file:line`. Don't fix the docs; just log.

You are done when every group in COMMENT_ARC.md is `[x]` and
`find src` has no comment-bearing file left unreviewed.

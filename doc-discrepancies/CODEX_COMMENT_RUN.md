# Autonomous comment-truthfulness run (codex-owned)

You (codex) own the `src/` comment-truthfulness pass from here. Run it to
completion, in chunks, committing as you go. This is COMMENT-ONLY work on a
FROZEN library — you change comments, never code. Work the repo at
`/shared/2026Thesis/nmr-shielding`.

## Stance — calm, unhurried, truth-seeking
Take your time. There is no deadline and no token pressure — this run may take
days, and that is fine. Better one file made truthful than ten rushed.

Go check the truth of things freely. When a comment makes a claim you cannot
immediately confirm, that is not something to cut your way out of — go read: the
`.cpp`, the header, the callers, the spec docs, a cited paper, whatever settles
it. Investigating is the work, not a detour. You have the access and the time.

Don't freak out. Uncertainty is a cue to slow down and look, never to slash in a
panic. A large file or a thick claim is worked methodically, not feared. You are
a careful reader making true statements truer — not a deadline-driven trimmer.

## Read first
`doc-discrepancies/COMMENT_ARC.md` — use it for the 16-group arc, the per-group
file membership, and the progress checkboxes. Its older doctrine (cutting cruft,
killing darlings, "all changelog goes", the deep-dive keep/cut test) is SUPERSEDED
by this file — we are no longer editing on style. The only doctrine now is "The
one job" below: make untrue comments true.

## The one job: make every comment TRUE
That is the entire task. For each comment, verify it against the code — read the
`.cpp`, the callers, a cited source; take the time to be sure. Then exactly one of:
- **True → leave it untouched.** Banners (`// ====`, `// ----`), section headers,
  changelog notes, restatements of adjacent code, decorative asides, long essays
  (NamingRegistry included) — all stay, exactly as written. We are not cutting,
  shortening, de-cluttering, or restyling anything.
- **Untrue → correct it, or delete it.** If the code does not do what the comment
  says, rewrite it into a brief, true, plain statement; or, if it is not worth
  re-grounding, delete it. Correcting is preferred to deleting.

Untruth is the only trigger. If a comment is true, you do nothing to it — however
redundant, decorative, verbose, or stale-looking it is.

## Hard rules
1. **COMMENT-ONLY. Never change code** — not a statement, signature, initializer,
   include, whitespace, or a lint "fix". After editing each file, run the
   self-check below; if it reports CODE CHANGED, `git checkout -- <file>` and redo
   comment-only.
2. **SOURCE IS TRUTH.** Verify each comment against the `.cpp` implementation
   (header comments lie too). A claim that is completely out of date may be
   *removed* — you need not re-ground it.
3. **Don't judge style — only truth.** We are not shortening, de-cluttering,
   removing banners or changelog, de-textbooking, or killing darlings. A true
   comment stays whatever its shape — verbose, decorative, a section banner, an
   essay. The ONLY reason to touch a comment is that it is untrue.
4. **Citations/factoids are claims.** Verify (research is allowed). If you cannot
   confirm the source/value, reword the comment to `probable source: <ref>` —
   never present an unconfirmed reference as authoritative. A confirmed, on-topic
   reference (e.g. Johnson & Bovey 1958, Caldeweyher 2019) stays as a real citation.
5. **Finish each file in one pass.** No half-done files; if finishing needs you to
   read a calculator / TrajectoryAtom / the .cpp, go read it now.
6. **Don't push. Don't touch** `src/generated/` (generated), or anything outside
   `src/`.
7. **NO META, no secret history.** When you fix or delete an untruth, leave NO
   trace of the edit: no "corrected", "was/now", "previously", "formerly X", no
   date, no parenthetical preserving the old wrong version, no scaffolding added
   to justify that a claim is true. The result must read as if it had always been
   brief and true. Truth silently replaces untruth.
8. **Stage by explicit path only.** The working tree holds unrelated uncommitted
   WIP outside `src/` AND the restored (reverted) `src/` files. Commit ONLY the
   files you edited, by path (`git add src/Foo.cpp`). NEVER `git add -A`, `git add
   .`, or `git commit -a`. Never stage anything outside `src/`.

## Method per file
read the file + its implementation → correct or delete ONLY the untrue comments,
leaving every true comment untouched → self-check (below) → commit with a clear
message (`docs(src comments): <file> — fix untrue comments`) → tick the file/group
in COMMENT_ARC.md. Commit per file or per small coherent group. Re-read the diff
before committing: every changed line should be an untruth you fixed or removed,
nothing else.

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
-name '*.h' -o -name '*.cpp' \)` against the group's membership rule. The
COMMENT_ARC checkboxes are the source of truth for progress; the working tree was
reset to the pre-run state, so trust the checkboxes and the file contents, not the
log. WARNING: the git log holds a prior `docs(src comments): G1…G7` run that was
REVERTED — do NOT treat those commit subjects as "done". Tick boxes in
COMMENT_ARC.md as you complete files/groups, so the next run resumes cleanly.

## Log doc discrepancies (optional, cheap)
If a comment reveals that a `.md` doc (OBJECT_MODEL, ARCHITECTURE, CATALOGUE,
API.md, …) is wrong about the code, append one line to
`doc-discrepancies/<UTC-timestamp>-codex-disc.md`:
`doc:line | doc claim | code truth | src file:line`. Don't fix the docs; just log.

You are done when every group in COMMENT_ARC.md is `[x]` and
`find src` has no comment-bearing file left unreviewed.

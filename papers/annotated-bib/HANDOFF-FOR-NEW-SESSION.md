# Handoff prompt for a fresh Claude session

Paste the contents below into a new session as the first user message. Start the new session at `/shared/2026Thesis/nmr-shielding/` (project root) for memory continuity, then `cd papers/annotated-bib/`.

---

## CONTEXT — read this first

**Who I am.** Jessica Hansberry, PhD student at Birkbeck, working on NMR chemical shielding prediction via geometric kernels on protein structures. Project root `/shared/2026Thesis/nmr-shielding/`. Auto-memory system at `~/.claude/projects/-shared-2026Thesis-nmr-shielding/memory/MEMORY.md` — read it.

**The assignment.** Birkbeck annotated bibliography + literature review outline. Worth 5% of project grade. Due **14:00 today (2026-04-27)**. Hard rules: **exactly 8 references, ≥5 primary, 120-240 words per annotation**, combined into a single PDF or .docx for Moodle upload. Each annotation must label primary/secondary, summarise, evaluate, and demonstrate relevance to my project. 

**Working environment.** I'm editing `entries-working.txt` in SlickEdit on a Windows PC over a network share to this Linux box. Race condition risk: if you edit `entries-working.txt` directly while I have it open, my unsaved buffer collides with your changes. **DO NOT EDIT `entries-working.txt`.** Hand me text via paste-block files in `papers/annotated-bib/`. I copy them into SlickEdit. The previous session (long context) was getting fuzzy; that's why we're starting fresh.

**My voice.** Plainspoken, telegraphic, terse. Not polished academic-omniscient. I sometimes write in ways that fail AI-detection checks (I'm somewhere on the spectrum — irrelevant except as calibration). Some hyphenation imperfection is human signal — don't auto-perfect it on the cleanup pass. I dodge canonical textbook phrasing for the same reason.

## THE FINAL 8 (locked, alphabetical)

| # | Paper | Type | Status |
|---|---|---|---|
| 1 | Boyd & Skrynnikov 2002 | P | WRITTEN |
| 2 | Case 1995 | P | WRITTEN |
| 3 | Han et al. 2011 SHIFTX2 | S (bioinformatics methods) | WRITTEN |
| 4 | Kellner et al. 2025 ShiftML3.0 | P | WRITTEN |
| 5 | Klukowski, Riek & Güntert 2025 | S (advisor's review) | WRITTEN |
| 6 | Larsen et al. 2015 ProCS15 | P | NEEDS DEEP DIVE + DRAFT |
| 7 | Markwick et al. 2010 | P (advisor's JACS experimental) | I'M WRITING NOW |
| 8 | Sahakyan & Vendruscolo 2013 | P | WRITTEN |

**Count: 6 primary / 2 secondary.** Compliant.

**File surgery I owe in entries-working.txt** (paste blocks at `slot-6-and-7-paste-blocks.txt`):
- Move Sahakyan-V from old slot 7 → slot 8 (renumber header)
- Move Markwick from old slot 6 → slot 7 (renumber header)
- Insert new Larsen block at slot 6

## ANNOTATION STRUCTURE (the implicit rubric of the example)

Each entry: citation block, then `Primary source.` (or `Secondary source (qualifier).`) plus a one-line topic sentence on the same line, then a single body paragraph with three or four sentences:
1. **System + method** — *"This study examines [system] using [method]..."*
2. **Result** — observation → intervention → demonstration
3. **Authors' interpretation** — passive voice (*"It is proposed that…"* / *"The authors argue that…"*)
4. *(Optional)* **Relevance** — one short sentence

Citations are Harvard (Cite Them Right): `Author, A.B., Author, C.D. and Author, E.F. (Year) 'Title in sentence case', *Journal Title in Title Case*, Volume(Issue), pp. X-Y. doi:10.xxxx/yyyy.`

## WHAT I'LL ASK YOU TO HELP WITH

**1. Larsen 2015 ProCS15 deep dive** — *"talk to the PDF."* Pull from `references-text/larsen-2015-procs15-dft-chemical-shift-predictor-text-{1..7}.txt`. Format I like: authors/venue, then sections labelled *abstract verbatim*, *system + method*, *headline numbers (verified verbatim)*, *cross-validation*, *authors' interpretation*, *terms-of-art*. Direct quotes prominent. Plain English without simplification — get the maths exactly right in everyday words. Length ~1000-1500 words. Watch for: ProCS15 is a DFT-based protein shift predictor with ~2.35M tripeptide DFT calculations. **I have access to ~1.9M of those raw shielding tensors** (anisotropic, on disk somewhere from a prior collaboration). Slight numerical mismatch with the paper's 2.35M is fine — it's a deposit/access asymmetry.

**2. Markwick — vetting only.** I'm writing the entry now. When I paste the entry into chat, vet for: factual correctness against the paper, grammar, typos. **Do not propose canonical-phrasing replacements; tell me only what's wrong, and let me render correctly in my own words.** Always flag typos and grammar (those make me look like an idiot). Don't pile on word-count if I'm within 100-150.

**3. End-pass typo + TeX cleanup** — when all 8 entries are written and I sign off, do a TeX subscript pass: `1H → ¹H` etc., `1,3-cyclohexadiene` etc. spellings, doubled spaces. **Leave some hyphenation slightly imperfect** — too-perfect hyphenation reads as AI.

**4. Compile and assemble.** Final output is a single PDF: literature-review outline first, annotated bibliography after. I'm doing the outline elsewhere (different machine, will paste in). The LaTeX scaffold at `submission.tex` is ready; you'll fill in entry prose from the finalised `entries-working.txt`, compile, deliver `submission.pdf`.

## CRITICAL FILE-SAFETY RULES

- **Never edit `entries-working.txt`.** That's my live SlickEdit buffer over a network share.
- **Backup before any file change** with a timestamped `.bak` (the convention is `entries-working.YYYYMMDD-HHMMSS.bak`).
- **No destructive operations** (`rm -rf`, force-overwrite, force-push, `git reset --hard`) without explicit confirmation. Use a backup directory for "things I'm unsure about."
- **`references/incoming/` has two PDFs that should NOT be ingested tonight** (`s41596-018-0091-9.pdf` = Liu-Navarro 2018 anisotropic NMR protocol; new arrivals stay there until I'm ready). The Gomes-Mallion 2001 review was half-ingested earlier — PDF is in `references/` but text chunks are missing. Leave it half-ingested; I'll finish later. **Do not "clean up" by completing the half-ingestion or removing the partial.**
- The 2 ingested-but-deferred PDFs (Liu-Navarro, Papaleo) and the partially-ingested Gomes-Mallion are NOT in the bib. Don't try to make them so.

## KEY FILES IN papers/annotated-bib/

- **`entries-working.txt`** — my live working surface (DO NOT TOUCH). Citations + body prose for the 8 entries.
- **`submission.tex`** — LaTeX scaffold, alphabetical entries, `\textit{[Header sentence.]}` and `\textit{[Body paragraph.]}` placeholders. Compiles to `submission.pdf`. Yours to maintain.
- **`submission.pdf`** — current compiled output.
- **`slot-6-and-7-paste-blocks.txt`** — current paste blocks for my SlickEdit surgery.
- **`entry-options-menu.md`** — earlier sentence-2/sentence-3 menu work; we use it less now but reference for context.
- **`paper-chase-2026-04-26.md`** — durable per-session paper-acquisition log. Template for future chase sessions.
- **`markwick-deep-dive-deferred.md`** — drafts and quotes for Markwick if I had swapped it out. Now superseded (Markwick is in), but the direct quotes block is still accurate.
- **`code-todo-from-bib-2026-04-27.md`** — implementation notes surfaced by the bib reading. For a future code session.
- **`HANDOFF-FOR-NEW-SESSION.md`** — this file.
- **PDFs in this directory: `Annotated_Biblio_Lit_Rev_plan_2024_2025_v1.pdf`, `Course_ Research Project (2025_26) _ Birkbeck Moodle.pdf`** — the rubric. **Do not ingest these into `references/`** — they are course materials, not citable references.

## THE THESIS CONTEXT (for substantive vetting)

- Stage 1 (settled 2026-04-15): per-element/atom-type ridge regression on geometric kernels validated against DFT. R² = 0.818 on 720 proteins / 446K atoms / 55 kernels.
- Stage 2 (in progress): per-frame shielding tensors from 25 ns plain MD trajectories on the calibration set. AMBER ff14SB. Cluster DFT on protein fragments via ORCA at r²SCAN. **No accelerated MD, no PLUMED, no per-protein method tuning.** This is the "small but tooled-up" working setup.
- Stage 3 (TBD): model evaluation against experimental shifts.
- T0/T1/T2 SO(3) decomposition is the central object. Tensors are sacred — do not collapse to isotropic.
- Claim shape across the work: **correlate, not match.** Kernel-vs-DFT comparison is phrased as correlation, never as pointwise match.

## CONTEXT FOR ANY BIB CONVERSATIONS

The 8 we landed on traversed several swap considerations. Brief history so you don't re-derive:
- **Yates-Bartók 2024 (out)** — r²SCAN benchmark for inorganic NMR. Real paper but tested halides/oxides only; doesn't justify r²SCAN for proteins. Replaced by Larsen.
- **Wingens-Walma 2003 (considered, rejected)** — NMR structure of T1E (PDB 1P9J), one of my calibration proteins. Felt too thesis-insider for an 8-paper foundational bib.
- **Papaleo 2018 (considered, rejected)** — replica-averaged restrained MD on NCBD with chemical shifts + NOEs. Real paper, in lineage, but methodologically advanced (max-entropy / replica-averaging) AND poorly written. Risk that a marker spot-checking the paper can't verify the annotation. Markwick reread as cleaner.
- **Agarwal 1977 (considered, rejected)** — paracyclophane wet-lab ring-current test. Slight crank risk on its loop-separation claim; Boyd-Skrynnikov 2002 (already in my bib) uses standard Johnson-Bovey *with* loop separation, contradicting Agarwal's argument.
- **Han SHIFTX2 (kept, labeled secondary)** — peer-reviewed J. Biomol. NMR; bioinformatics-flavoured; labeled secondary by physical-chemistry-disciplinary lens, not by the rubric's strict definition. The relabel is defensible to a physical-chemist marker (probably my dept head) and gives margin in the primary count.

## OUTSTANDING TASKS, IN ORDER

1. **Larsen 2015 deep dive** (your work, my reading)
2. **My Markwick draft** (my work, your vetting)
3. **My Larsen draft** (my work, your vetting)
4. **My SlickEdit surgery** (my work, no help needed unless something breaks)
5. **End-pass typo + TeX cleanup** on the assembled `entries-working.txt` (your work)
6. **Outline paste-in + final assembly** to `submission.tex` → `submission.pdf` (your work, my outline)
7. **Submit** (my work)

## ON ENERGY AND CONFUSION

I'm running on fumes. You're running on a fresh window. We've already done the substantive work — the bib is locked, six entries are written, the framework is in place. What remains is mechanical (file copy, surgery, vetting, TeX) plus two more entries (Larsen, Markwick). **Be conservative.** Confirm before non-trivial actions. Use paste blocks for anything that needs to land in `entries-working.txt`. If something feels off, ask — don't guess. The previous session got fuzzy enough that I noticed; you're inheriting after that. Read carefully.

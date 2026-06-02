# Brief: the nmr_extract use guide

You're writing the operator's guide for `nmr_extract`: how to run it, what inputs
each mode needs, and the rules those inputs must satisfy. Someone reaches for this
when they have structures in hand and need to run the tool correctly the first time.
It has to be right — every flag, every file requirement, every rejection rule vetted
against the code, because a wrong instruction here wastes a long extraction run.

Take your time. There is no rush and no token budget to protect on your end.

## Who reads this

Someone who will run `nmr_extract` from a Linux command line and needs to get the
inputs exactly right. Comfortable in a shell; competent with structural-biology file
formats (PDB, GROMACS `.tpr`/`.trr`, AMBER `prmtop`) at a working level; not assumed
to know this codebase or to have read its design docs. Write so they can prepare
inputs and run each mode without guessing. Define a convention the first time it
matters, then use it.

## The standard

Real-plain grounding, no fake-plain filler. No superlatives or absolutes. State each
rule as a rule: what is required, what is rejected, what the default is. Where a value
has units or a specific threshold, give the exact value from the code. Do not invent
flags, defaults, or file names — every one must be in the source.

## The code is the final authority

The CLI surface lives in `src/Cli/`: `ModeSpec.h` (mode and flag definitions),
`Parse.cpp`, `PrintUsage.cpp`, `CommonOptions.h`, and `Validate.cpp`/`Validate.h`
(the input-validation rules and their rejection messages). `src/RunConfiguration.{h,cpp}`
is the resolved run shape. The readers are `src/FullSystemReader.{h,cpp}` (the
`--trajectory DIR` file derivation), `src/PdbFileReader.{h,cpp}`, and
`src/ReduceProtonation.{h,cpp}`. These files are the final authority for which flags
exist, what inputs are required, and what is rejected. Where any prose disagrees with
them, the code is right.

## The canonical mode set — do not exceed it

`CLAUDE.md` at the repo root carries the canonical 5-mode specification: the only five
supported modes; MOPAC an orthogonal `--mopac`/`--no-mopac` toggle; APBS and AIMNet2
always on and not switchable; propka, kaml, `--analysis`, and `--fleet` out of scope;
`--no-apbs`/`--no-coulomb` removed. Document exactly those five modes and no more. If
you find a flag in the code the canonical spec does not list, flag the discrepancy
rather than presenting it as a feature.

## What the guide covers

1. What `nmr_extract` does, in two or three sentences, and the always-on truth (APBS
   and AIMNet2 run in every mode; MOPAC is the toggle).
2. The five modes, in order, each with the flag(s) that select it, the exact input
   files it requires and where it expects them, and what it emits: (1) `--pdb`,
   (2) `--protonated-pdb`, (3) `--orca --root NAME`, (4) `--mutant --wt NAME --ala NAME`,
   (5) `--trajectory DIR`.
3. The GROMACS trajectory layout for mode 5: the required directory contents and how
   the reader derives each file by name (no globbing), the hard requirement that the
   run carry coordinates in a `.trr` (quote the `.xtc`-rejected validation check and
   its condition), the independent NPY stride `m` and PDB stride `n`, and the
   `--mopac` FullFat variant.
4. The mutant pair layout for mode 4: what the WT and ALA poses each are (mode-3
   poses) and how they are named and located.
5. PDB rules for modes 1–2: bare versus protonated, what `reduce` does in mode 1,
   ff14SB charge assignment, and any residue/atom expectations the readers enforce.
6. The MOPAC toggle and what it changes (the run shape per mode).
7. Out of scope: the flags and use cases that are explicitly not supported, so a
   reader does not go looking for them.

## Method

Write the draft, then re-read it against the code: is every flag real, every required
file correct, every rejection rule quoted accurately, every default the code's actual
default? Fix what drifted. Then the copy pass — no superlatives, no idea-smearing, no
unexplained term.

## LaTeX and PDF

`article`, 11pt, `geometry`, `amsmath` only if needed, and a monospace treatment for
commands and paths (`\texttt`, or a `verbatim`/`listings` block). Show example
invocations as command blocks. Length follows completeness — be complete and exact, do
not pad; there is no fixed page cap, but every page must earn its place. Compile with
`pdflatex` to a clean build: no overfull boxes, no undefined references.

## Done when

Every mode, flag, required file, and rejection rule matches `src/Cli/` and the readers;
the documented mode set matches the canonical 5-mode spec; example invocations are
correct; it compiles clean; a competent operator could prepare inputs and run any mode
from it without guessing.

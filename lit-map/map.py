#!/usr/bin/env python3
"""lit-map aggregator — turn a crawl's jsonl into a dead-zone ELIMINATION WORKLIST.

The model is corpus-blind; this corpus-aware step decides have-vs-missing. A real dead zone
is "a statement we want to say but cannot source"; at the area level we approximate it as the
foundational refs our corpus leans on but does not hold. Crucially, every paper in our jsonl
IS held, so the papers citing a missing ref are its secondary-citation candidates, already on
the shelf. We surface, per missing ref: demand + those candidates, to drive the elimination
(acquire primary / secondary-cite after judging value / drop the claim).

Usage:  map.py <slug>
Reads   lit-map/results/<slug>.jsonl + corpus inventory (references-meta/*-summary.txt)
Writes  lit-map/results/<slug>-MAP.md and prints it.
"""
import json, os, re, sys, glob
from collections import Counter, defaultdict

ROOT = "/shared/2026Thesis/nmr-shielding"
slug = sys.argv[1]
jsonl = f"{ROOT}/lit-map/results/{slug}.jsonl"
out_md = f"{ROOT}/lit-map/results/{slug}-MAP.md"

corpus = [os.path.basename(p)[:-len("-summary.txt")].lower()
          for p in glob.glob(f"{ROOT}/references-meta/*-summary.txt")]

YEAR = re.compile(r'\b(1[89]\d{2}|20\d{2})\b')
PARTICLES = {"van", "von", "de", "der", "den", "del", "di", "la", "le", "el", "and"}

def cite_key(cite):
    ym = YEAR.search(cite)
    if not ym:
        return None
    toks = [t for t in re.findall(r"[A-Za-z][A-Za-z'\-]+", cite[:ym.start()])
            if t.lower() not in PARTICLES]
    if not toks:
        return None
    surname = toks[0].lower().strip("-'")          # "Pople-1956" -> "pople", not "pople-"
    return (surname, ym.group(0)) if surname else None

def held(key):
    s, y = key
    return any(b.startswith(s) and y in b for b in corpus)

def cites_of(content):
    for line in content.splitlines():
        if line.startswith("CITES:"):
            return line.split(":", 1)[1].strip()
    return None

modes = Counter()
support = defaultdict(list)
refs = defaultdict(lambda: {"papers": set(), "raw": set(), "held": False})
anchors = defaultdict(set)   # held paper -> {missing ref keys it cites}

records = 0
for ln in open(jsonl):
    ln = ln.strip()
    if not ln:
        continue
    try:
        r = json.loads(ln)          # tolerate a partial trailing line on a live file
    except json.JSONDecodeError:
        continue
    records += 1
    eng = r.get("engages", "")
    modes[eng] += 1
    if eng == "no":
        continue
    support[eng].append(r["basename"])
    cites = cites_of(r.get("content", ""))
    if not cites or cites.lower() == "none":
        continue
    for c in re.split(r'[;]', cites):
        c = c.strip()
        if not c or c.lower() == "none":
            continue
        k = cite_key(c)
        if not k:
            continue
        refs[k]["papers"].add(r["basename"])
        refs[k]["raw"].add(c)
        refs[k]["held"] = held(k)
        if not refs[k]["held"]:
            anchors[r["basename"]].add(k)

dead = {k: v for k, v in refs.items() if not v["held"]}
have = {k: v for k, v in refs.items() if v["held"]}
engaged = sum(v for m, v in modes.items() if m != "no")

L = [f"# Lit-map dead-zone worklist — {slug}\n",
     f"Records {records} · engaged {engaged} · corpus {len(corpus)} · dead-zone refs {len(dead)}\n",
     "## Depth — how much of the corpus bears on this area"]
for m in ("deals-with", "touches-intro", "touches-discussion", "yes", "no"):
    if modes.get(m):
        L.append(f"- {m}: {modes[m]}")

L.append("\n## Support set — held papers to cite")
for m in ("deals-with", "touches-intro", "touches-discussion", "yes"):
    if support.get(m):
        L.append(f"\n**{m}** ({len(support[m])})")
        L += [f"- {b}" for b in sorted(support[m])]

L.append("\n## DEAD ZONE worklist — cited but NOT held (acquire OR secondary-cite)")
L.append("_demand = # held papers leaning on it; listed papers = secondary-citation candidates_\n")
for k in sorted(dead, key=lambda k: (-len(dead[k]["papers"]), k)):
    d = dead[k]
    L.append(f"- **{k[0]} {k[1]}**  demand {len(d['papers'])}  ({' / '.join(sorted(d['raw']))})")
    L.append(f"    - secondary via: {', '.join(sorted(d['papers']))}")

L.append("\n## Best secondary anchors — held papers covering the most dead zones")
for b in sorted(anchors, key=lambda b: -len(anchors[b]))[:15]:
    if anchors[b]:
        L.append(f"- {b} — covers {len(anchors[b])}: " +
                 ", ".join(f"{s} {y}" for s, y in sorted(anchors[b])))

L.append("\n## Held foundational refs (on disk)")
for k in sorted(have, key=lambda k: -len(have[k]["papers"])):
    L.append(f"- {k[0]} {k[1]} ×{len(have[k]['papers'])}")

report = "\n".join(L) + "\n"
open(out_md, "w").write(report)
print(report)

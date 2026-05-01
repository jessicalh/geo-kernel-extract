#!/usr/bin/env python3
"""qwen_summary.py — Draft 2 summary via a local OpenAI-compatible LLM.

Reads chunked text from references-text/<basename>-text-N.txt,
sends them in one shot to a long-context model, writes the result to
references-meta/<basename>-summary-qwen-test.txt by default.

Test-only path. Canonical summaries (`<basename>-summary.txt`) are
untouched until output is vetted.
"""

import argparse
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path("/shared/2026Thesis/nmr-shielding")

SYSTEM_PROMPT = """\
You are summarising a scholarly paper for an NMR shielding-tensor PhD thesis. The thesis builds geometric kernels — ring-current Biot-Savart and Haigh-Mallion, electric-field-gradient, bond magnetic anisotropy, hydrogen-bond — on protein structures, calibrated against DFT shielding values. Rank-2 tensor (T2) preservation is the central novelty; the system outputs kernels, and calibration turns kernels into shielding.

The summary you write will be the durable record of this paper for the project. The thesis author, committee, and future readers will rely on it weeks and months from now to make decisions about citations, methods, and arguments — they will not re-read the paper. The summary IS what is remembered. So accuracy and specificity matter; a vague summary is worse than no summary, and a fabricated detail is worse than an honest gap.

A good summary preserves the wealth of true detail in the paper:
- Authors, year, title, journal, volume, page range. DOI when present in the text. Reference count and any funding or credit notes that are stated.
- Specific numerical findings — r and r² values, RMSEs, intensity factors per residue, coefficient values, dataset sizes, distance ranges, basis sets, functionals, force fields. Quote the numbers the paper actually gives.
- The paper's claims, distinct from what the surrounding field claims.
- Specific authors and works the paper builds on or points toward, named from the paper's own bibliography.
- The paper's stance — research paper, Perspective, review, calibration study, methodological note. State it honestly.

A good summary avoids:
- Inventing methodology to fill a slot. If the chunks do not state how something was done (e.g. error estimation), do not name a technique. The reader will trust your characterisation, so silence beats fabrication.
- Manufacturing thesis connections that are not there. Most papers are useful background; few are direct ancestors. Do not stretch.
- Excited or evaluative language ("groundbreaking", "elegant", "comprehensive"). The register is dry plain English; define terms of art briefly inline.
- Editorialising — "our approach differs" framing belongs elsewhere, not in a summary.

Two project conventions to honour:
- Kernel-vs-DFT comparisons are framed as correlation: the r² shows the signal was captured. They are not framed as match (pointwise agreement). Default vocabulary: "correlate", "captures the signal", "r² of …". The thesis's Stage 1 R² = 0.818 is a correlation claim, not an accuracy claim.
- Per-element complexity is the substance of the thesis (Stage 1: H ≈ 20 effective dims, C ≈ 6, N ≈ 3, O ≈ 12). Do not collapse a tensor to a scalar except where the source paper itself does.

Format.

Header paragraph: citation (authors, year, title, journal, volume, pages), DOI when visible, page count, reference count if stated, and any load-bearing access, funding, or credit notes.

Then five paragraphs, no section titles, paragraphs flow:
1. What the paper is — scope, stance, audience, what it argues or measures.
2. Methods or central results (research) / central theorems or organising principles (review).
3. Findings — magnitudes, what was observed, what was not.
4. Relations to the wider literature — the works this paper builds on and points toward, named from its own bibliography.
5. Implications for our work, considered widely and honestly. Do not stretch beyond what the paper supports.

Length scales with paper size: 4–6 pages → 200–300 words body; 6–12 pages → 300–450; 10+ page dense review or Perspective → 450–500.

Take the time to think carefully about the paper before composing. You are a careful reader writing for a careful reader. The numbers, the names, the distinctions — those are the substance.
"""

VERIFY_PROMPT = """\
Please give your draft summary one more careful read against the paper text, with the accuracy of the science in mind. The thesis author and committee will trust this summary as the durable record of the paper, so a wrong number or a mis-attributed method becomes a real error in the work that follows. Getting it right matters.

Three things in particular deserve a check.

Statistics — every number in the draft: r and r² values, RMSEs, dataset sizes, distance ranges, intensity factors per residue, coefficient values, page and reference counts. For each, find it in the chunks and confirm it matches. If a number is off or you cannot find it, correct it.

Methodology — the procedures and tools attributed to the paper: functionals, basis sets, force fields, fitting procedures, error-estimation techniques, dataset construction. Confirm the paper says these are what the authors did, not just what they cite or compare against. If anything is mis-attributed, correct it.

Names — every author, year, work, residue, ring, molecule, or chemical species named in the draft. A small substitution can lose a real distinction. Confirm each name against the paper.

If the draft is already accurate, return it unchanged. If something needs correction, return a revised summary in the same header-plus-five-paragraphs format, with the corrections applied. Do not rewrite for style — only adjust where the science requires it. Reply with the final summary only.
"""


def discover_model(endpoint_base: str) -> str:
    with urllib.request.urlopen(f"{endpoint_base}/models", timeout=10) as r:
        data = json.load(r)
    return data["data"][0]["id"]


def read_chunks(basename: str) -> str:
    chunk_dir = ROOT / "references-text"
    paths = sorted(
        chunk_dir.glob(f"{basename}-text-*.txt"),
        key=lambda p: int(p.stem.rsplit("-", 1)[1]),
    )
    if not paths:
        sys.exit(f"no chunks found for basename: {basename}")
    parts = []
    for p in paths:
        parts.append(f"=== {p.name} ===\n{p.read_text()}")
    return "\n\n".join(parts)


def strip_thinking(text: str) -> tuple[str, str]:
    """Return (body_without_think, think_block_or_empty)."""
    m = re.search(r"<think>(.*?)</think>\s*", text, flags=re.DOTALL)
    if not m:
        return text, ""
    return (text[: m.start()] + text[m.end() :]).lstrip(), m.group(1).strip()


def chat_complete(endpoint: str, payload: dict) -> tuple[dict, float]:
    """POST a chat completion. Return (parsed_json, wall_seconds)."""
    req = urllib.request.Request(
        f"{endpoint}/chat/completions",
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
    )
    t0 = time.time()
    with urllib.request.urlopen(req, timeout=1800) as r:
        result = json.load(r)
    return result, time.time() - t0


def extract_message(result: dict) -> tuple[str, str, str]:
    """Pull (content, reasoning_content, finish_reason) from a chat completion result.

    Promotes any embedded <think>...</think> block in content into reasoning_content
    when the server didn't separate them on its own.
    """
    msg = result["choices"][0]["message"]
    content = msg.get("content") or ""
    reasoning = msg.get("reasoning_content") or ""
    content, embedded = strip_thinking(content)
    if embedded and not reasoning:
        reasoning = embedded
    finish = result["choices"][0].get("finish_reason", "?")
    return content, reasoning, finish


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("basename")
    ap.add_argument(
        "--out",
        help="output path (default: references-meta/<basename>-summary-qwen-test.txt)",
    )
    ap.add_argument("--endpoint", default="http://127.0.0.1:8099/v1")
    ap.add_argument("--model", help="model id (default: auto-discovered)")
    ap.add_argument("--temperature", type=float, default=0.3)
    ap.add_argument(
        "--no-vet",
        action="store_true",
        help="skip the second-pass self-verification step",
    )
    args = ap.parse_args()

    model = args.model or discover_model(args.endpoint)
    chunks = read_chunks(args.basename)

    system = SYSTEM_PROMPT
    user = (
        f"Paper basename: {args.basename}\n\n"
        "Below is the paper text in 3-page chunks. Produce the Draft 2 summary as specified.\n\n"
        f"{chunks}"
    )

    print(f"model:    {model}", file=sys.stderr)
    print(f"endpoint: {args.endpoint}/chat/completions", file=sys.stderr)
    print(f"system:   {len(system)} chars", file=sys.stderr)
    print(f"user:     {len(user)} chars ({len(chunks)} chars chunks)", file=sys.stderr)

    messages = [
        {"role": "system", "content": system},
        {"role": "user", "content": user},
    ]
    base_payload = {"model": model, "temperature": args.temperature, "stream": False}

    # ── pass 1: draft ─────────────────────────────────────────────────────
    print("pass 1: drafting...", file=sys.stderr)
    result1, dt1 = chat_complete(args.endpoint, {**base_payload, "messages": messages})
    draft, reasoning1, finish1 = extract_message(result1)
    usage1 = result1.get("usage", {})
    print(
        f"  prompt={usage1.get('prompt_tokens', '?')} completion={usage1.get('completion_tokens', '?')} "
        f"finish={finish1} time={dt1:.1f}s draft={len(draft)} chars",
        file=sys.stderr,
    )

    summary = draft
    reasoning2 = ""
    if not args.no_vet and draft:
        # ── pass 2: vet ──────────────────────────────────────────────────
        print("pass 2: vetting...", file=sys.stderr)
        messages.append({"role": "assistant", "content": draft})
        messages.append({"role": "user", "content": VERIFY_PROMPT})
        result2, dt2 = chat_complete(args.endpoint, {**base_payload, "messages": messages})
        verified, reasoning2, finish2 = extract_message(result2)
        usage2 = result2.get("usage", {})
        print(
            f"  prompt={usage2.get('prompt_tokens', '?')} completion={usage2.get('completion_tokens', '?')} "
            f"finish={finish2} time={dt2:.1f}s verified={len(verified)} chars",
            file=sys.stderr,
        )
        if verified:
            summary = verified

    out_path = (
        Path(args.out)
        if args.out
        else (ROOT / "references-meta" / f"{args.basename}-summary-qwen-test.txt")
    )
    out_path.write_text(summary)

    # If the verifier changed the draft, keep both for inspection.
    if not args.no_vet and draft and summary != draft:
        draft_path = out_path.with_name(out_path.stem + ".draft.txt")
        draft_path.write_text(draft)
        print(f"pre-vet draft saved to {draft_path}", file=sys.stderr)

    if reasoning1 or reasoning2:
        parts = []
        if reasoning1:
            parts.append(f"=== pass 1 (draft) ===\n{reasoning1}")
        if reasoning2:
            parts.append(f"=== pass 2 (vet) ===\n{reasoning2}")
        reasoning_path = out_path.with_name(out_path.stem + ".reasoning.txt")
        reasoning_path.write_text("\n\n".join(parts))
        print(
            f"reasoning saved to {reasoning_path} "
            f"(pass1={len(reasoning1)} pass2={len(reasoning2)})",
            file=sys.stderr,
        )

    print(f"wrote {out_path} ({len(summary)} chars)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

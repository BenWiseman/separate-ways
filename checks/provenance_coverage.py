#!/usr/bin/env python3
"""Which numbers in the manuscript have no row in CLAIMS.tsv.

OVER-REPORTS BY DESIGN; READ, DO NOT OBEY.  It cannot tell a computed result from an arXiv
identifier that happens to look like a decimal, a cited literature value, or a parameter
setting, so a clean run is not the goal and a long list is not an alarm.  What it is for is
the reading pass: go down the list and ask of each entry whether a script in the release
produces it.

Written 2026-09-24 after the availability note was checked against the table it describes.
The note says every quantitative claim is produced by a script and that CLAIMS.tsv maps each
claim to its file.  Of 131 numbers carrying three or more decimals, 65 had a row.  Thirty of
the remainder were genuinely computed.  Running every script and grepping its output
attributed 22 of those; the other eight were the interesting ones, and among them were a
table with no script anywhere in the release, a threshold whose last digit was wrong, and a
pair of numbers nothing produced.

    python3 checks/provenance_coverage.py [paper.md]

Add --near to list pairs of numbers that agree to about three significant figures and
differ beyond them.  Most pairs are unrelated quantities and a few are the same quantity
quoted at two precisions, both of which are fine; what the mode is for is the third case,
the same quantity quoted twice with different digits.  It found the decoherence threshold
written as 1.95374 in Section 2.1 and 1.95376 two paragraphs later, against a root of
1.9537646713, and the same split in the proper time.

Add --attribute to run every script in checks/calc and far_side and report which of them
prints each uncovered number.  That is slow, a few minutes, and it is what turns a list of
numbers into a list of rows to add.
"""
import glob
import os
import re
import subprocess
import sys

PAPER = "paper/COMPANION_v1.md"
CLAIMS = "checks/CLAIMS.tsv"


def uncovered(paper=PAPER, claims=CLAIMS):
    text = open(paper, encoding="utf-8").read()
    cut = text.find("\n## References")
    body = text[:cut] if cut > 0 else text
    table = open(claims, encoding="utf-8").read()
    flat = " ".join(body.split())
    out = {}
    for m in re.finditer(r"(?<![\w.^{])(\d+\.\d{3,})(?![\d}])", flat):
        n = m.group(1)
        if n in table:
            continue
        ctx = flat[max(0, m.start() - 95):m.start() + 55]
        if re.search(r"arXiv:|DOI|10\.1\d{3}", ctx):      # identifiers, not results
            continue
        out.setdefault(n, ctx)
    return out


def attribute(numbers):
    scripts = sorted(glob.glob("checks/calc/*.R") + glob.glob("checks/calc/*.py")
                     + glob.glob("checks/*.R"))
    found = {n: [] for n in numbers}
    for s in scripts:
        cmd = ["Rscript", s] if s.endswith(".R") else [".venv/bin/python", s]
        try:
            res = subprocess.run(cmd, capture_output=True, text=True, timeout=240).stdout
        except Exception:
            continue
        printed = [float(x) for x in re.findall(r"-?\d+\.\d+", res)]
        for n in numbers:
            # compare numerically at the precision the paper quotes, which is what the
            # gate does. Substring matching found "2.241" inside "32.241896" and
            # attributed a value to a script that never computed it.
            tol = 0.5 * 10 ** -len(n.split(".")[1])
            if any(abs(abs(p) - float(n)) <= tol for p in printed):
                found[n].append(os.path.basename(s))
    return found


def near_duplicates(paper=PAPER, rel_max=2e-3):
    """Numbers close enough to be the same quantity written twice."""
    text = open(paper, encoding="utf-8").read()
    cut = text.find("\n## References")
    flat = " ".join((text[:cut] if cut > 0 else text).split())
    seen = {}
    for m in re.finditer(r"(?<![\w.^{])(\d+\.\d{2,})(?![\d}])", flat):
        if re.search(r"arXiv:|DOI|10\.1\d{3}", flat[max(0, m.start() - 60):m.start()]):
            continue
        seen.setdefault(m.group(1), m.start())
    nums = sorted(seen, key=float)
    out = []
    for a, b in zip(nums, nums[1:]):
        fa, fb = float(a), float(b)
        if fa == 0 or a == b:
            continue
        rel = abs(fb - fa) / max(abs(fa), abs(fb))
        if 0 < rel < rel_max:
            out.append((a, b, rel, flat[max(0, seen[a] - 62):seen[a] + 18],
                        flat[max(0, seen[b] - 62):seen[b] + 18]))
    return out


if __name__ == "__main__":
    paper = next((a for a in sys.argv[1:] if not a.startswith("--")), PAPER)
    if "--near" in sys.argv:
        pairs = near_duplicates(paper)
        print(f"=== {paper}: {len(pairs)} near-duplicate pair(s)")
        print("  over-reports by design; look for the same quantity with different digits\n")
        for a, b, rel, ca, cb in pairs:
            print(f"  {a} vs {b}  (rel {rel:.1e})\n     A ...{ca}\n     B ...{cb}\n")
        sys.exit(0)
    miss = uncovered(paper)
    print(f"=== {paper} against {CLAIMS} ===")
    print(f"  {len(miss)} number(s) of three or more decimals with no row")
    print("  over-reports by design; read, do not obey\n")
    if "--attribute" in sys.argv:
        found = attribute(sorted(miss))
        for n in sorted(miss, key=float):
            who = ", ".join(found[n]) if found[n] else "NO SCRIPT PRINTS IT"
            print(f"  {n:<14} {who}")
    else:
        for n, ctx in sorted(miss.items(), key=lambda kv: float(kv[0])):
            print(f"  {n:<14} ...{ctx[-82:]}")

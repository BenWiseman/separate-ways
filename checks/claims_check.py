#!/usr/bin/env python3
"""Verify every claim in CLAIMS.tsv, now that the manuscript carries no script paths.

Ben, 2026-09-22 and again on 2026-09-24: "you can't have script paths in a manuscript at all.
Nobody does that." He is right, and the provenance those paths carried had to survive the strip
rather than go with it. `checks/CLAIMS.tsv` is that provenance, one row per claim:

    script <TAB> anchor, the first words of the passage it backs <TAB> numbers claimed there

This checks both halves of each row. The anchor must still be findable in the manuscript, so a
rewrite that moves a passage out from under its script is caught. And every number listed must
appear in that script's output at the precision the paper states it, which is what the inline
citations used to buy.

    python3 checks/claims_check.py paper/COMPANION_v1.md
    python3 checks/claims_check.py paper/COMPANION_v1.md --reanchor

An anchor is a phrase, so it breaks whenever the passage it names is reworded. That is the point:
a passage moving out from under its script has to be noticed. But it must not become a wall, so
--reanchor shortens every broken anchor from the end, word by word, until it lands again, and
says which rows it touched. Read what it prints. If an anchor shortens to almost nothing, the
passage did not move, it went, and the row needs a person rather than a trim.
"""
import re, io, os, sys, subprocess, tempfile

# A sibling number_provenance.py wins over the one in the authoring toolkit, so the
# export repository runs on the copy it actually ships rather than on a private path
# that happens to exist on the author's machine.
sys.path.insert(0, os.path.expanduser("~/benlm/tools"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import number_provenance as np

# Paths are overridable so the same checker runs in the export repository, where the
# release sits under paper/ and checks/ instead of paper/ and checks/. The defaults
# are the working-tree layout, so nothing here changes for a plain run. Without this the
# export carried a README promising a pass that could not start.
TSV   = os.environ.get("CLAIMS_TSV", "checks/CLAIMS.tsv")
PAPERS = os.environ.get("CLAIMS_PAPERS_DIR", "pub/paper2")
PAPER = sys.argv[1] if len(sys.argv) > 1 else os.path.join(PAPERS, "COMPANION_v1.md")

# CLAIMS.tsv names each script by its working-tree path. CLAIMS_SCRIPT_MAP rewrites those
# for another layout, as "from=to,from=to", longest prefix first.
_MAP = [tuple(x.split("=", 1)) for x in
        os.environ.get("CLAIMS_SCRIPT_MAP", "").split(",") if "=" in x]
_MAP.sort(key=lambda kv: -len(kv[0]))

def remap(path):
    for a, b in _MAP:
        if path.startswith(a):
            return b + path[len(a):]
    return path

REANCHOR = "--reanchor" in sys.argv
rows = [l.rstrip("\n").split("\t") for l in io.open(TSV, encoding="utf-8")
        if l.strip() and not l.startswith("#")]
assert rows, "CLAIMS.tsv is empty: the check would pass by matching nothing"

# An anchor may live in EITHER manuscript. CLAIMS.tsv was written when only the companion was
# gated, so a row anchored in the cosmology paper reported "anchor no longer in the paper" and
# looked like a regression; that happened three times in one session before this was changed.
# Both files are searched, which keeps the moved-anchor plant working because a planted anchor is
# in neither.
SIBLING = os.path.join(PAPERS, "PAPER2_v4_draft.md" if "COMPANION" in PAPER
                              else "COMPANION_v1.md")
paper = " ".join(io.open(PAPER, encoding="utf-8").read().split())
if os.path.exists(SIBLING):
    paper = paper + "  " + " ".join(io.open(SIBLING, encoding="utf-8").read().split())
cache = tempfile.mkdtemp(prefix="claims_")
scripts = sorted({remap(r[0]) for r in rows})
print(f"  running {len(scripts)} scripts behind {len(rows)} claims...")
np.run_all(scripts, cache)

if REANCHOR:
    paper_flat = paper
    lines = io.open(TSV, encoding="utf-8").read().rstrip("\n").split("\n")
    out, fixed = [], []
    for l in lines:
        if l.startswith("#") or not l.strip():
            out.append(l); continue
        sc, an, nu = l.split("\t")
        if an in paper_flat:
            out.append(l); continue
        w = an.split()
        while w and " ".join(w) not in paper_flat:
            w = w[:-1]
        if w:
            fixed.append((sc, len(an.split()), len(w), " ".join(w)))
            out.append(f"{sc}\t{' '.join(w)}\t{nu}")
        else:
            out.append(l)
    io.open(TSV, "w", encoding="utf-8").write("\n".join(out) + "\n")
    print(f"  re-anchored {len(fixed)} row(s):")
    for sc, n0, n1, a in fixed:
        flag = "   <-- shortened to almost nothing; check the passage still exists" if n1 < 4 else ""
        print(f"    {sc}: {n0} words -> {n1}{flag}\n      \"{a[:74]}\"")
    if not fixed:
        print("    none needed")

lost_anchor, not_found, checked = [], [], 0
for script0, anchor, nums in rows:
    script = remap(script0)
    if " ".join(anchor.split()) not in paper:
        lost_anchor.append((script, anchor)); continue
    out = os.path.join(cache, script.replace("/", "_") + ".out")
    if not os.path.exists(out):
        not_found.append((script, "(script produced no output)")); continue
    outs = np.nums_from_out(io.open(out, encoding="utf-8", errors="replace").read())
    for lit in [x for x in nums.split(";") if x]:
        got = np.nums_from_text(lit)
        if not got: continue
        val, dec, _, sci = got[0]
        checked += 1
        if not np.matches(val, dec, sci, outs):
            not_found.append((script, lit))

if "--selftest" in sys.argv:
    # Gate 4 is the gate that matters most and nothing validated it, which is how the
    # sibling cross-document check sat unreachable for a day. Three things have to hold:
    # a corrupted number must be caught, a moved anchor must be caught, and a figure the
    # script merely quotes from the paper must not count as evidence.
    ok = True
    bad_val = np.matches(9.999, 3, False, np.nums_from_out("computed 1.234 (paper: 9.999)"))
    good_val = np.matches(1.234, 3, False, np.nums_from_out("computed 1.234 (paper: 9.999)"))
    print("  claims_check: a quoted figure is not evidence: %s" % ("yes" if not bad_val else "NO"))
    print("  claims_check: the computed one on the same line still is: %s"
          % ("yes" if good_val else "NO"))
    ok &= (not bad_val) and good_val
    # a corrupted claim must fail against its own script
    sc, an, nu = next((r for r in rows if "." in r[2].split(";")[0]), rows[0])
    lit = nu.split(";")[0]
    val, dec, _, sci = np.nums_from_text(lit)[0]
    outs = np.nums_from_out(io.open(os.path.join(cache, remap(sc).replace("/", "_") + ".out"),
                                    encoding="utf-8", errors="replace").read())
    caught = not np.matches(val * 1.7 + 0.31, dec, sci, outs)
    print("  claims_check: a corrupted number is caught: %s" % ("yes" if caught else "NO"))
    ok &= caught
    moved = " ".join("a planted anchor that is not in the manuscript".split()) not in paper
    print("  claims_check: a missing anchor is caught: %s" % ("yes" if moved else "NO"))
    ok &= moved
    print("  claims_check: the real table passes: %s"
          % ("yes" if not (lost_anchor or not_found) else "NO"))
    ok &= not (lost_anchor or not_found)
    sys.exit(0 if ok else 1)

print(f"  {len(rows)} claims, {checked} numbers checked against the script named for them")
for s, a in lost_anchor:
    print(f"      anchor no longer in the paper: {s}  <-- ISSUE\n        \"{a[:78]}\"")
for s, l in not_found:
    print(f"      {l:>20}  not in the output of {s}   <-- ISSUE")
if not lost_anchor and not not_found:
    print("  every anchor still lands and every number reproduces")
sys.exit(1 if (lost_anchor or not_found) else 0)

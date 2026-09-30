#!/usr/bin/env python3
"""Regenerate the arXiv metadata counts from the manuscript.

Gate 8 refuses to let ARXIV_METADATA_COMPANION.txt drift from COMPANION_v1.md, and the COMMENTS
line states four counts that change whenever the body does. Doing that by hand after every edit
is a step that gets skipped, which is how the gate went red twice in one session with the paper
itself perfectly fine. Run this instead:

    python3 checks/sync_metadata.py

It rewrites only the abstract and the four counts. Everything Ben fills in by hand is untouched.
"""
import re, io, sys

PAPER = "papers/2_over_the_horizon/COMPANION_v1.md"
META  = "papers/2_over_the_horizon/ARXIV_METADATA_COMPANION.txt"

c = io.open(PAPER, encoding="utf-8").read()
body_all, _, refs = c.partition("\n## References")
body, _, app = body_all.partition("## Appendices")
strip = lambda x: re.sub(r"!\[.*?\]\([^)]*\)", " ", x, flags=re.S)
nw   = len(strip(body).split()) + len(strip(app).split())
napp = len(strip(app).split())
nfig = len(re.findall(r"(?m)^!\[", c))
nref = len(re.findall(r"(?m)^\d+\\?\.\s+[A-Z]", refs))

m = re.search(r"(?m)^##\s+Abstract\s*$", c)
rest = c[m.end():]
abstract = rest[:re.search(r"(?m)^#{1,6}\s+\S", rest).start()].strip()

t = io.open(META, encoding="utf-8").read()
before = t
t = re.sub(r"(?ms)(^ABSTRACT\n).*?(\n\n#)", lambda g: g.group(1) + abstract + g.group(2), t, count=1)
t = re.sub(r"Companion to arXiv:\[Separate Ways identifier, once assigned\]\. .*?references\.",
           f"Companion to arXiv:[Separate Ways identifier, once assigned]. {nw} words, of which "
           f"{napp} are\nthe appendices; {nfig} figures; {nref} references.", t, flags=re.S, count=1)
io.open(META, "w", encoding="utf-8").write(t)
print(f"  {META}: {nw} words, {napp} appendix, {nfig} figures, {nref} references, "
      f"abstract {len(' '.join(abstract.split()))} chars"
      + ("" if t != before else "   (already current)"))

# The availability note counts the repository, and the repository grows while the paper is
# being written, so gate 12 went red on every commit that added a script. Sync it here rather
# than by hand: the count is still checked, it just stops being a chore that teaches me to
# ignore a red gate.
import glob
ncalc = len(glob.glob("checks/calc/*.R")) + len(glob.glob("checks/calc/*.py"))
nfg   = len(glob.glob("checks/fig_*.R"))
c2, n = re.subn(r"the \d+ calculation files and the \d+ figure generators",
                f"the {ncalc} calculation files and the {nfg} figure generators", c, count=1)
assert n == 1, "the availability note no longer reads as this script expects"
if c2 != c:
    io.open(PAPER, "w", encoding="utf-8").write(c2)
print(f"  {PAPER}: availability note says {ncalc} calculation files, {nfg} figure generators")

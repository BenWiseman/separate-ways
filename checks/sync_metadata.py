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
# "about N", not N. Gate 8 requires the approximate form, because a word count moves on every
# prose edit and a gate that fails on every edit is a gate somebody deletes. This script was
# writing the exact form and sync_release_counts.py was writing the approximate one into the
# same sentence; this one ran second, so the gate failed with both tools reporting success.
# One owner now, and it is this one.
t = re.sub(r"Companion to arXiv:\[Separate Ways identifier, once assigned\]\. .*?references\.",
           f"Companion to arXiv:[Separate Ways identifier, once assigned]. about {round(nw, -2)} "
           f"words, of which about {round(napp, -2)} are\nthe appendices; {nfig} figures; "
           f"{nref} references.", t, flags=re.S, count=1)
io.open(META, "w", encoding="utf-8").write(t)
print(f"  {META}: {nw} words, {napp} appendix, {nfig} figures, {nref} references, "
      f"abstract {len(' '.join(abstract.split()))} chars"
      + ("" if t != before else "   (already current)"))

# The availability note counts the repository, and the repository grows while the paper is
# being written, so gate 12 went red on every commit that added a script. Sync it here rather
# than by hand: the count is still checked, it just stops being a chore that teaches me to
# ignore a red gate.
import glob, os
calcs = sorted(glob.glob("checks/calc/*.R") + glob.glob("checks/calc/*.py"))
ncalc = len(calcs)
npy   = len([f for f in calcs if f.endswith(".py")])
nR    = ncalc - npy

# Only four of these five counts were ever synced, and the other four were typed. All four were
# wrong: the note claimed 21 generators where the tree has 22, said all of them were R where
# three are Python, said three of the Python calculations used the standard library alone where
# six do, and said 15 generators draw the figures here where 17 do. Counted now, so they cannot
# drift again. A fig_*.R that never opens a device is a helper, not a generator: fig_label.R
# masks labels for five others and draws nothing itself.
def _draws(f):
    return "dev.off()" in io.open(f, encoding="utf-8", errors="replace").read() \
        or "savefig(" in io.open(f, encoding="utf-8", errors="replace").read()
gens   = [f for f in sorted(glob.glob("checks/fig_*.R") + glob.glob("checks/fig_*.py"))
          if _draws(f)]
ngen   = len(gens)
ngen_R = len([f for f in gens if f.endswith(".R")])
# R that loads no package: base R is then enough to run it.
ngen_R_base = len([f for f in gens if f.endswith(".R")
                   and not re.search(r"(?m)^\s*(library|require)\s*\(",
                                     io.open(f, encoding="utf-8", errors="replace").read())])
# Python calculations importing nothing outside the standard library.
_SCI = re.compile(r"(?m)^\s*(?:import|from)\s+(numpy|scipy|sympy|mpmath|matplotlib|pandas)\b")
nstd = len([f for f in calcs if f.endswith(".py")
            and not _SCI.search(io.open(f, encoding="utf-8", errors="replace").read())])
# Generators that draw a figure this paper actually embeds.
embedded = set(re.findall(r"\]\((fig_[a-z0-9_]+\.pdf)\)", c))
# Match on the stem, not the full filename: three generators build the name as
# sprintf("papers/2_over_the_horizon/fig_companion_two_ends.%s", dev), so ".pdf" never appears in them. A
# loose substring fallback overcounted by one, because "contact" sits in two generator names.
stems = {pdf[:-4] for pdf in embedded}
ndraw = len({f for f in gens
             if any(stem in io.open(f, encoding="utf-8", errors="replace").read()
                    for stem in stems)})

c2 = c
subs = [
  (r"the \d+ calculation files and the \d+ figure generators",
   f"the {ncalc} calculation files and the {ngen} figure generators"),
  (r"Of those, \d+ calculations and (?:all \d+|\d+ of the) generators are R and load no package",
   f"Of those, {nR} calculations and {ngen_R_base} of the generators are R and load no package"),
  (r"the other \d+ calculations are Python, (?:three|four|five|six|seven|eight|nine|ten|\d+) using the standard library alone",
   f"the other {npy} calculations are Python, {nstd} using the standard library alone"),
  (r"\d+ of those generators draw the figures here",
   f"{ndraw} of those generators draw the figures here"),
]
WS = re.compile(r"[ ]")
for pat, rep in subs:
    # hard-wrapped prose: a space in the pattern may be a newline in the file
    c2, n = re.subn(WS.sub(r"\\s+", pat), rep, c2, count=1)
    assert n == 1, "the availability note no longer reads as this script expects: " + pat
if c2 != c:
    io.open(PAPER, "w", encoding="utf-8").write(c2)
print(f"  {PAPER}: availability note says {ncalc} calculations ({nR} R, {npy} Python, "
      f"{nstd} stdlib-only), {ngen} generators ({ngen_R} R, {ngen_R_base} needing base R only), "
      f"{ndraw} drawing figures here")

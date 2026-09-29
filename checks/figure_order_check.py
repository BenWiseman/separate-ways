#!/usr/bin/env python3
"""Every numbered figure must be cited in the text, and cited before it appears.

    python3 checks/figure_order_check.py [files...]
    python3 checks/figure_order_check.py --selftest

WHY. Checked on 2026-09-29 at Ben's instruction and three of the cosmology paper's five figures
failed: Figures 1 and 2 were never referred to in the text at all, and Figure 3 only after it
had already appeared. A figure nobody points at is a figure the reader meets without knowing
why, and journals ask for the reference to come first. Nothing measured this before.

A citation is any "Figure N" outside that figure's own caption. The caption itself does not
count, which is the whole trap: a naive search finds the caption and reports the figure cited.
"""

import io
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DEFAULT = [os.path.join(ROOT, "pub", "paper2", "PAPER2_v4_draft.md"),
           os.path.join(ROOT, "pub", "paper2", "COMPANION_v1.md")]


def audit(text):
    """Return [(label, n_citations, n_before_caption)] for every numbered figure."""
    flat = " ".join(text.split())
    out = []
    for cap in re.finditer(r"\*\*(Figure\s+\d+)\.\*\*", flat):
        label = " ".join(cap.group(1).split())
        n = label.split()[-1]
        span = (cap.start(), cap.end())
        cites = [m.start() for m in re.finditer(r"Figure\s+" + n + r"(?![.\d])", flat)
                 if not (span[0] <= m.start() < span[1])]
        out.append((label, len(cites), len([c for c in cites if c < span[0]])))
    return out


def report(paths):
    bad = 0
    for p in paths:
        if not os.path.exists(p):
            continue
        text = io.open(p, encoding="utf-8").read()
        rows = audit(text)
        name = os.path.basename(p)
        includes = len(re.findall(r"!\[", text))
        broken = [r for r in rows if r[2] == 0]
        print("   %-24s %d figures, %d numbered, %d cited before they appear"
              % (name, includes, len(rows), len(rows) - len(broken)))
        # A document whose figures carry no numbers cannot cite them, and this gate would
        # otherwise pass it by matching nothing. The companion is in that state: 17 figures,
        # all with markdown alt-text captions and no "Figure N" label. Journals want numbers.
        if includes and not rows:
            print("      %d figures and none numbered, so none can be cited  <-- NOT GATED"
                  % includes)
        for label, cites, before in broken:
            bad += 1
            why = "never cited in the text" if cites == 0 else "cited only after it appears"
            print("      %s: %s  <-- ISSUE" % (label, why))
    return bad


def _selftest():
    ok = True
    good = "As Figure 1 shows, the modes cross. ![](f.png) **Figure 1.** The crossing."
    r = audit(good)
    a = r == [("Figure 1", 1, 1)]
    print("   plant: a figure cited before it appears is spared: %s" % ("yes" if a else "NO"))
    ok = ok and a
    never = "The modes cross. ![](f.png) **Figure 1.** The crossing."
    b = audit(never) == [("Figure 1", 0, 0)]
    print("   plant: a figure never cited is caught: %s" % ("yes" if b else "NO"))
    ok = ok and b
    after = "![](f.png) **Figure 1.** The crossing. Figure 1 showed the crossing."
    r3 = audit(after)
    c = r3 == [("Figure 1", 1, 0)]
    print("   plant: a figure cited only afterwards is caught: %s" % ("yes" if c else "NO"))
    ok = ok and c
    # the trap: the caption must not count as its own citation
    d = audit("![](f.png) **Figure 1.** The crossing.") == [("Figure 1", 0, 0)]
    print("   plant: a caption does not count as a citation: %s" % ("yes" if d else "NO"))
    return ok and d


def main():
    if "--selftest" in sys.argv:
        sys.exit(0 if _selftest() else 1)
    paths = [a for a in sys.argv[1:] if not a.startswith("-")] or DEFAULT
    sys.exit(1 if report(paths) else 0)


if __name__ == "__main__":
    main()

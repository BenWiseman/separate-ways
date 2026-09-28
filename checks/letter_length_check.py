#!/usr/bin/env python3
"""The Letter has to fit PRL's length limit, and it will grow if nothing watches it.

PRL allows 3750 words in the core of the paper, between the abstract and the
acknowledgements, plus up to two pages of End Matter that do not count against it
(journals.aps.org/prl/authors, read 2026-09-28). Displayed equations and figures count
toward the 3750, so counting prose alone would report a comfortable margin that does
not exist.

The estimate here is ours and not PRL's. Their submission tool does the real count and
weighs a figure by its aspect ratio; we charge a flat 150 words for a single-column
figure and 16 for a displayed equation, which is the usual rule of thumb. Treat a pass
as "there is room", not as "PRL will agree". The margin is deliberately reported so a
draft creeping toward the cap is visible before it arrives.
"""

import io
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
LETTER = os.path.join(ROOT, "pub", "paper2", "LETTER_PRL_v1.md")

CAP = 3750
PER_EQUATION = 16
PER_FIGURE = 150


def measure(text):
    """Return (prose, equations, figures, charged total) for the Letter's core."""
    try:
        core = text.split("## The Letter", 1)[1]
    except IndexError:
        raise SystemExit("LETTER: no '## The Letter' heading, cannot find the core")
    core = core.split("## References", 1)[0]
    neq = core.count("$$") // 2
    # A figure is charged ONCE. It normally appears twice, as an image include and as a
    # numbered caption, and adding the two charged this Letter's single figure 300 words.
    nfig = max(len(re.findall(r"\*\*FIG\.\s*\d", core)),
               len(re.findall(r"^!\[", core, re.M)))
    prose = re.sub(r"\$\$.*?\$\$", " ", core, flags=re.S)
    prose = re.sub(r"\$[^$]*\$", " x ", prose)
    prose = re.sub(r"\*\*FIG\.\s*\d.*", " ", prose)
    words = len(prose.split())
    return words, neq, nfig, words + neq * PER_EQUATION + nfig * PER_FIGURE


def _selftest(text):
    """Make it fail on purpose: a cap that never fires is not a cap."""
    ok = True
    base = measure(text)[3]
    padded = text.replace("## References", ("filler " * 4000) + "\n\n## References", 1)
    over = measure(padded)[3]
    caught = over > CAP >= base
    print("   plant: 4000 words of filler is caught: %s" % ("yes" if caught else "NO"))
    ok = ok and caught
    # An equation must be charged, or the estimate is prose-only by another name.
    with_eq = text.replace("## References", "$$a=b$$\n\n## References", 1)
    charged = measure(with_eq)[3] - base == PER_EQUATION
    print("   plant: a displayed equation is charged %d words: %s"
          % (PER_EQUATION, "yes" if charged else "NO"))
    ok = ok and charged
    # And so must a figure.
    with_fig = text.replace("## References", "**FIG. 9.** A plant.\n\n## References", 1)
    charged_f = measure(with_fig)[3] - base == PER_FIGURE
    print("   plant: a figure is charged %d words: %s"
          % (PER_FIGURE, "yes" if charged_f else "NO"))
    ok = ok and charged_f
    # The figure already in the Letter has both an include and a caption, and charging
    # both put 150 phantom words on the count. One figure, one charge.
    once = measure(text)[2] == len(re.findall(r"^!\[", text.split("## The Letter", 1)[1],
                                              re.M))
    print("   plant: one figure is charged once, not twice: %s" % ("yes" if once else "NO"))
    return ok and once


def main():
    text = io.open(LETTER, encoding="utf-8").read()
    words, neq, nfig, total = measure(text)
    print("   core %d prose + %d eq x%d + %d fig x%d = %d charged, cap %d, margin %d"
          % (words, neq, PER_EQUATION, nfig, PER_FIGURE, total, CAP, CAP - total))
    dashes = text.count("—")
    print("   em dashes: %d" % dashes)
    ok = _selftest(text)
    if total > CAP or dashes or not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()

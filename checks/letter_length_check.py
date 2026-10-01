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

# The body heading is not a fixed string. It was "## The Letter" until 2026-09-29, when the
# manuscript went to Nature Physics as an Article and the Physical Review furniture was renamed:
# "## The Letter" became "## Main" and "## End Matter" became "## Appendices". Both checkers
# split on the old literal, so both raised SystemExit, and the length and number checks they
# exist to run did not run at all while the gate still printed. A checker that cannot find its
# subject must say so loudly rather than exit quietly mid-suite, and it must not be pinned to a
# journal's house style when the manuscript may be sent to any journal.
BODY_HEADINGS = ("## Main", "## The Letter", "## Letter")


def _find_heading(text, headings):
    """Split on a heading only where it starts a line.

    A plain `h in text` test is wrong and was: "## Appendix" is a substring of the subsection
    heading "### Appendix A.", so renaming the section heading still matched, the split landed
    mid-document and the region came back partial instead of raising. Caught by planting a
    renamed heading and finding the guard did not fire."""
    for h in headings:
        m = re.search(r"^" + re.escape(h) + r"\s*$", text, re.M)
        if m:
            return text[m.end():]
    return None


def split_body(text, what):
    """Return everything after whichever body heading this manuscript uses."""
    got = _find_heading(text, BODY_HEADINGS)
    if got is not None:
        return got
    raise SystemExit(
        "LETTER: none of %s found, cannot find the %s. If the heading was renamed again, add it "
        "to BODY_HEADINGS rather than letting this check pass over an empty string."
        % (", ".join(repr(h) for h in BODY_HEADINGS), what))



HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
LETTER = os.path.join(ROOT, "pub", "paper2", "LETTER_PRL_v1.md")

CAP = 3750
PER_EQUATION = 16
PER_FIGURE = 150


def measure(text):
    """Return (prose, equations, figures, charged total) for the Letter's core."""
    try:
        core = split_body(text, "core")
    except IndexError:
        raise SystemExit("LETTER: could not split the core")
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
    # And so must a figure. The old plant put "**FIG. 9.** A plant." immediately before
    # "## References", which measure() cuts away before it counts anything, so the delta was
    # always zero and this line printed NO for as long as it existed: a plant that cannot
    # fire in the checker that certifies the Letter fits PRL's cap. Two things were wrong.
    # It has to land inside the core, and it has to use the include form this Letter uses,
    # where the caption sits inside the image include and not in a separate **FIG.** line.
    # Third time. The first plant sat immediately before "## References", which measure()
    # cuts away. The second anchored on "\n### ", and this Letter's only "### " headings are
    # its appendices, which sit AFTER "## References" and are cut away too, so the plant went
    # on printing NO. Anchor on the core itself, which split_body already locates, and there
    # is no heading left to rename out from under it.
    _core = split_body(text, "core")
    _at = text.index(_core)
    with_fig = text[:_at] + "\n![FIG. 9.](plant.pdf)\n\n" + text[_at:]
    if with_fig == text:
        print("   plant: a figure is charged %d words: NO (the plant did not change the text)"
              % PER_FIGURE)
        return False
    b_words, _, b_nfig, b_tot = measure(text)
    p_words, _, p_nfig, p_tot = measure(with_fig)
    # Subtract the plant's own prose so the figure charge is measured on its own.
    charged_f = (p_nfig == b_nfig + 1
                 and p_tot - b_tot - (p_words - b_words) == PER_FIGURE)
    print("   plant: a figure is charged %d words: %s"
          % (PER_FIGURE, "yes" if charged_f else "NO"))
    ok = ok and charged_f
    # The figure already in the Letter has both an include and a caption, and charging
    # both put 150 phantom words on the count. One figure, one charge.
    once = measure(text)[2] == len(re.findall(r"^!\[", split_body(text, "core"),
                                              re.M))
    print("   plant: one figure is charged once, not twice: %s" % ("yes" if once else "NO"))
    return ok and once


def main():
    # Measure the file we are given. Defaulting to LETTER silently while accepting an
    # argument meant every run reported on LETTER_PRL_v1.md whatever was asked for.
    path = sys.argv[1] if len(sys.argv) > 1 else LETTER
    text = io.open(path, encoding="utf-8").read()
    words, neq, nfig, total = measure(text)
    print("   %s" % os.path.basename(path))
    print("   core %d prose + %d eq x%d + %d fig x%d = %d charged, cap %d, margin %d"
          % (words, neq, PER_EQUATION, nfig, PER_FIGURE, total, CAP, CAP - total))
    dashes = text.count("—")
    print("   em dashes: %d" % dashes)
    ok = _selftest(text)
    if total > CAP or dashes or not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()

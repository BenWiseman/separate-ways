#!/usr/bin/env python3
"""Every number in the Letter must already appear in the manuscript.

The Letter is carved out of PAPER2_v4_draft.md, whose numbers are checked against the
scripts that produce them by claims_check. Rather than duplicate that machinery, this
gate says the Letter quotes nothing the manuscript does not, so the existing check
covers it. A number in the Letter that is absent from the manuscript is either a typo
or a claim with no script behind it, and both have to be caught.

The trap this went through first: tokenising on bare digit runs splits
``$3\\times10^{-14}$`` into "3", "10" and "14", each of which appears somewhere in a
49-page paper, so a corrupted exponent passes. Numeric expressions are matched whole,
LaTeX and all, and the plants at the bottom include a corrupted exponent for that
reason.

What this gate cannot do, said plainly so nobody reads more into a pass than it carries.
It is a substring test against a 49-page paper, so a short number is effectively
unchecked: "52" corrupted to "53", or "2.1" to "2.7", both still match somewhere and
both pass. What it does catch is the distinctive kind, exponents, error bars and values
carried to three or four figures, which is where a transcription slip between the two
documents actually changes a claim. It is a guard against drift, not a proof of
agreement.
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
# The appendix heading moved with the body heading, and this one failed SILENTLY. The old code
# read `rest.split("## End Matter", 1)[1] if "## End Matter" in rest else ""`, so after the
# rename the conditional took the empty branch and the whole appendix stopped being checked
# while the gate still printed its numbers. Only the self-test caught it, by reporting that its
# planted value in that section was not found. An absent heading now raises.
APPENDIX_HEADINGS = ("## Appendices", "## End Matter", "## Appendix")


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
PAPER = os.path.join(ROOT, "pub", "paper2", "PAPER2_v4_draft.md")

# A numeric expression: a decimal, optionally times a power of ten, optionally with a
# +- error, optionally with a trailing unit word. The \times10^{..} part is glued on so
# an exponent cannot be checked apart from its mantissa.
NUM = re.compile(
    r"""
    (?<![\w.])
    \d[\d,]*(?:\.\d+)?
    (?: \s*\\pm\s*\d[\d,]*(?:\.\d+)? )?
    (?: \s*\\times\s*10\^\{?-?\d+\}? )?
    (?!\d)(?!\.\d)
    """,
    re.VERBOSE,
)

# Numbers that are not physics: reference volumes and pages, years, section numbers.
def _skip(tok, context):
    if re.search(r"\[\d+\]\s*[A-Z]", context):      # a reference-list line
        return True
    if re.match(r"^(19|20)\d\d$", tok.strip()) and re.search(r"\b" + re.escape(tok) + r"s?\b", context):
        return True      # a year or a decade, e.g. "since the 1960s"
    return False


def _norm(s):
    return re.sub(r"\s+", "", s).replace(",", "")


def _checked_region(letter_text):
    """Everything after the title that carries a claim: the body, the figure caption and
    End Matter. Only the reference list is exempt, and it is cut by its own headings
    rather than by a test on the remainder, which silently took End Matter with it."""
    try:
        body = split_body(letter_text, "body")
    except IndexError:
        raise SystemExit("LETTER: could not split the body")
    if "## References" in body:
        before, rest = body.split("## References", 1)
        after = _find_heading(rest, APPENDIX_HEADINGS)
        if after is None:
            raise SystemExit(
                "LETTER: references found but none of %s after them, so the appendix would go "
                "unchecked. Add the new heading to APPENDIX_HEADINGS rather than letting this "
                "check pass over an empty string."
                % ", ".join(repr(h) for h in APPENDIX_HEADINGS))
        return before + "\n" + after
    return body


def unmatched(letter_text, paper_text):
    """Numeric expressions in the checked region that the manuscript does not carry."""
    body = _checked_region(letter_text)
    hay = _norm(paper_text)
    out = []
    for m in NUM.finditer(body):
        tok = m.group(0)
        ctx = body[max(0, m.start() - 90):m.end() + 30].replace("\n", " ")
        if _skip(tok, ctx):
            continue
        if _norm(tok) in hay:
            continue
        out.append((tok, ctx))
    return out


def _counted(letter_text):
    """How many expressions the gate actually looked at, so the total cannot overstate it."""
    body = _checked_region(letter_text)
    return sum(1 for m in NUM.finditer(body)
               if not _skip(m.group(0), body[max(0, m.start() - 90):m.end() + 30]))


def _plants(letter_text, paper_text):
    """The gate has to fail on purpose before its pass means anything."""
    base = len(unmatched(letter_text, paper_text))
    cases = [
        # (description, from, to) -- each must raise the unmatched count
        ("a corrupted exponent", r"3\times10^{-14}", r"3\times10^{-11}"),
        ("a corrupted mantissa", "491.6", "491.7"),
        ("a corrupted plain number", "58.8 meV", "58.9 meV"),
        ("a corrupted error bar", r"245.8\pm1.0", r"245.8\pm1.4"),
        ("a corrupted small integer", "factor of 21.6", "factor of 21.7"),
        ("a sentence-final number", "moves to 5.04", "moves to 5.07"),
        # The appendix section is all small integers, which a substring test cannot discriminate,
        # so this plant inserts a distinctive value there instead of corrupting one.
        ("a number in the appendices", "Two of the four are what writing a metric theory means",
                                   "Two of the four, to $7.3194\\times10^{-8}$, are what writing a metric theory means"),
    ]
    ok = True
    for what, a, b in cases:
        if a not in letter_text:
            print("   plant NOT APPLICABLE, %s: %r absent from the Letter" % (what, a))
            ok = False
            continue
        n = len(unmatched(letter_text.replace(a, b, 1), paper_text))
        caught = n > base
        print("   plant: %s is caught: %s" % (what, "yes" if caught else "NO"))
        ok = ok and caught
    return ok


def main():
    letter = io.open(LETTER, encoding="utf-8").read()
    paper = io.open(PAPER, encoding="utf-8").read()
    bad = unmatched(letter, paper)
    total = _counted(letter)
    print("   %d numeric expressions in the Letter, %d not in the manuscript"
          % (total, len(bad)))
    for tok, ctx in bad:
        print("     %-22s ...%s..." % (tok, ctx))
    ok = _plants(letter, paper)
    if bad or not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()

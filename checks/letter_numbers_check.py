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
            # A Letter is allowed to have no appendix at all: the PRL submission uses no End
            # Matter, and this guard refused that outright. It exists to catch an appendix
            # hiding under a heading nobody added to the list, so look for one. If what
            # follows the references is reference entries and nothing else, there is no
            # appendix to miss.
            leftover = [l for l in rest.split("\n")
                        if l.strip() and not re.match(r"^\d+\\?\.\s", l.strip())
                        and not l.startswith(("  ", "\t"))]
            if any(l.lstrip().startswith("#") for l in leftover):
                raise SystemExit(
                    "LETTER: references are followed by a heading that is none of %s, so an "
                    "appendix would go unchecked. Add it to APPENDIX_HEADINGS rather than "
                    "letting this check pass over an empty string. Found: %s"
                    % (", ".join(repr(h) for h in APPENDIX_HEADINGS),
                       next(l.strip()[:48] for l in leftover if l.lstrip().startswith("#"))))
            return before
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
    # Every one of the seven above names a literal from LETTER_PRL_v1.md. Run against any
    # other Letter they all print NOT APPLICABLE, and the gate then validates nothing at all
    # while still looking like a gate. So build plants out of the document in hand: take one
    # numeric expression of each shape the checker distinguishes and corrupt its last digit.
    if not any(a in letter_text for _, a, _ in cases):
        body = _checked_region(letter_text)
        seen, cases = set(), []
        for m in NUM.finditer(body):
            tok = m.group(0)
            if _skip(tok, body[max(0, m.start() - 90):m.end() + 30]):
                continue
            # Small integers are not plantable and the original table said so: they are
            # everywhere in the corpus, so corrupting one lands on another real number and
            # the substring test cannot tell. "1,2" was picked as an integer and is a list.
            if "." not in tok and "times10" not in tok and "pm" not in tok:
                continue
            shape = ("exponent" if "times10" in tok else
                     "error bar" if "pm" in tok else "decimal")
            if shape in seen or "," in tok:
                continue
            bumped = None
            for i in range(len(tok) - 1, -1, -1):
                if tok[i].isdigit():
                    bumped = tok[:i] + str((int(tok[i]) + 3) % 10) + tok[i + 1:]
                    break
            if bumped and bumped != tok and letter_text.count(tok) == 1:
                seen.add(shape)
                cases.append(("a corrupted %s from this Letter (%s)" % (shape, tok), tok, bumped))
        # And one insertion, so the gate is shown to catch a number that is simply not in
        # any source, which is the failure it exists for.
        head = _checked_region(letter_text)[:400]
        anchor = next((l for l in head.split(". ") if len(l) > 40), None)
        if anchor and letter_text.count(anchor) == 1:
            cases.append(("a value present in no source", anchor,
                          anchor + ", to $7.3194\\times10^{-8}$"))
        if not cases:
            print("   plant: NO usable numeric expression to corrupt; the gate is unvalidated")
            return False

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
    # Take the file named on the command line. This read LETTER unconditionally, so every
    # run against LETTER2_CONTACT_v2.md silently checked LETTER_PRL_v1.md instead and
    # reported a clean pass on a document nobody had asked about. Its sibling
    # letter_length_check.py had the same defect. Refuse a path that is not there, rather
    # than falling back to the default and passing over the wrong paper again.
    path = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("-") else LETTER
    if not os.path.exists(path):
        raise SystemExit("LETTER: no such file: %s" % path)
    letter = io.open(path, encoding="utf-8").read()
    print("   reading %s" % os.path.basename(path))
    # The corpus a Letter's numbers may come from is not one manuscript. "On a road to
    # nowhere" was carved out of the companion, not out of PAPER2_v4_draft.md, and checking
    # it against the cosmology paper alone reported three numbers untraced that were all
    # sound: 0.89M and 6e-4 sit in COMPANION_v1.md, and 0.469 is printed by the Letter's own
    # Penrose generator and belongs in no manuscript at all. Read the companion and the
    # figure generators too, so a real orphan stands out instead of being lost in three
    # false ones.
    corpus = [PAPER, os.path.join(ROOT, "pub", "paper2", "COMPANION_v1.md")]
    corpus += [os.path.join(HERE, f) for f in sorted(os.listdir(HERE))
               if f.startswith("fig_letter_") and f.endswith((".py", ".R"))]
    paper = "\n".join(io.open(f, encoding="utf-8").read()
                      for f in corpus if os.path.exists(f))
    print("   against %d source(s): %s"
          % (len([f for f in corpus if os.path.exists(f)]),
             ", ".join(os.path.basename(f) for f in corpus if os.path.exists(f))))
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

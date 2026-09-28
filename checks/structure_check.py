#!/usr/bin/env python3
"""Structural faults that only surface when somebody tries to build the paper.

Three of them, all cheap, none covered by the other gates.

Unbalanced maths. Thirty edits in a day, every one splicing LaTeX into hard-wrapped prose, and an
odd number of dollar signs in one paragraph turns the rest of the document into mathematics.
Nothing else here reads the file as a renderer would.

A figure reference with no file. The arXiv build is the first thing that would notice, which on
the wrong day is the submission itself.

A figure file nobody references. Harmless to a build and a sign that something moved and left its
picture behind, which is how the dimension figure came to sit two thousand lines from its claim.

    python3 checks/structure_check.py paper/COMPANION_v1.md
    python3 checks/structure_check.py paper/COMPANION_v1.md --selftest
"""
import io, os, re, sys

PAPER = sys.argv[1] if len(sys.argv) > 1 else "paper/COMPANION_v1.md"


def maths_faults(text):
    out = []
    for n, para in enumerate(text.split("\n\n")):
        flat = " ".join(para.split())
        if flat.count("$$") % 2:
            out.append("unclosed $$ in paragraph %d: %s" % (n, flat[:70]))
        if len(re.findall(r"(?<!\$)\$(?!\$)", flat)) % 2:
            out.append("unclosed $ in paragraph %d: %s" % (n, flat[:70]))
    return out


def figure_faults(text, paper_path):
    here = os.path.dirname(paper_path) or "."
    refs = sorted(set(re.findall(r"\]\(([^)]+\.pdf)\)", text)))
    out = ["figure referenced but absent: %s" % r
           for r in refs if not os.path.exists(os.path.join(here, r))]
    present = sorted(f for f in os.listdir(here)
                     if f.startswith("fig_companion") and f.endswith(".pdf"))
    out += ["figure present but unreferenced: %s" % f for f in present if f not in refs]
    return refs, out


def layout_faults(text):
    """Two figure blocks with no prose between them.

    Invisible in markdown and unmissable in a rendered PDF, where the two captions run
    together and the reader loses which picture is being described. It happened in the
    introduction, where the silence figure and the Kruskal figure sat back to back while
    the argument the second one illustrates was still two paragraphs away.
    """
    paras = [q for q in text.split("\n\n") if q.strip()]
    isfig = lambda q: q.lstrip().startswith("![")
    return ["two figures with no prose between them: %s..." % " ".join(paras[i].split())[:54]
            for i in range(len(paras) - 1) if isfig(paras[i]) and isfig(paras[i + 1])]


def quote_faults(text):
    """A blockquote marker stranded in the middle of a line.

    Found on 2026-09-28 by reading the built PDF: page 46 of the companion carried "The >
    derivative acting on the", and the line after it, having lost its own marker, fell out of
    the quote block and set flush left. Invisible in markdown, where the stray character reads
    as punctuation, and unmissable once typeset. A ">" inside a quoted line, followed by a
    space and a lowercase word, is the signature; a mathematical one is written with a digit or
    a symbol after it and is spared.
    """
    out = []
    for i, line in enumerate(text.split("\n"), 1):
        if not line.lstrip().startswith(">"):
            continue
        body = line.lstrip()[1:]
        if re.search(r"\s>\s+[a-z]", body):
            out.append("line %d: a blockquote marker stranded mid-line, \"%s\""
                       % (i, " ".join(line.split())[:70]))
    return out


if __name__ == "__main__":
    text = io.open(PAPER, encoding="utf-8").read()

    if "--selftest" in sys.argv:
        ok = True
        planted = text + "\n\nA planted paragraph with $one unclosed delimiter in it.\n"
        got = maths_faults(planted)
        print("  structure_check: catches an unclosed $: %s" % ("yes" if got else "NO"))
        ok &= bool(got)
        planted2 = text + "\n\n![planted](fig_companion_does_not_exist.pdf)\n"
        _, got2 = figure_faults(planted2, PAPER)
        hit = any("absent" in g for g in got2)
        print("  structure_check: catches a missing figure: %s" % ("yes" if hit else "NO"))
        ok &= hit
        stacked = text + "\n\n![a](fig_companion_silence.pdf)\n\n![b](fig_companion_neff.pdf)\n"
        got3 = layout_faults(stacked)
        print("  structure_check: catches two stacked figures: %s" % ("yes" if got3 else "NO"))
        ok &= bool(got3)
        quoted = text + "\n\n> a quoted line with a stranded > marker in the middle of it\n"
        got4 = quote_faults(quoted)
        print("  structure_check: catches a stranded blockquote marker: %s" % ("yes" if got4 else "NO"))
        ok &= bool(got4)
        spared = text + "\n\n> a quoted line saying that $A$ > 1 is required\n"
        got5 = quote_faults(spared)
        print("  structure_check: spares a mathematical one: %s" % ("yes" if not got5 else "NO"))
        ok &= not got5
        clean = (maths_faults(text) + figure_faults(text, PAPER)[1] + layout_faults(text)
                 + quote_faults(text))
        print("  structure_check: the real file is clean: %s" % ("yes" if not clean else "NO"))
        ok &= not clean
        sys.exit(0 if ok else 1)

    faults = maths_faults(text)
    refs, figfaults = figure_faults(text, PAPER)
    faults += figfaults
    faults += layout_faults(text)
    faults += quote_faults(text)
    npara = len(text.split("\n\n"))
    if not faults:
        print("   %d paragraphs balanced, %d figures all present and all used" % (npara, len(refs)))
        sys.exit(0)
    print("   %d structural fault(s)  <-- ISSUE" % len(faults))
    for f in faults:
        print("      %s" % f)
    sys.exit(1)

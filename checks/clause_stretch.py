#!/usr/bin/env python3
"""No body sentence may make the reader hold more than N words without a hard break.

    python3 checks/clause_stretch.py FILE.md [FILE.md ...]
    python3 checks/clause_stretch.py --selftest

WHY. Ben, 2026-09-27: the paper has to read exceptionally well and get through peer review.
The word cap in sentence_length.py is the coarse bar and a long sentence broken by semicolons
and colons is not the problem; the reader gets a place to put the first half down. What costs a
reader is a long run of clauses with no hard break in it, which is what forces a re-read. This
measures that run and caps it. Appendices are exempt for the same reason as the word cap.

The caps are set at what each manuscript reaches today, with a little room, so they ratchet:
prose can get easier and cannot get harder without someone deciding to move the number.
"""
import io, re, sys

CAPS = {"PAPER2_v4_draft.md": 60, "COMPANION_v1.md": 70}
DEFAULT = 60


def stretches(path):
    """(words, text) for every run of prose between hard breaks, body only."""
    c = io.open(path, encoding="utf-8").read()
    body = c.split("## Appendi")[0]
    if "## 1. Introduction" in body:
        body = body[body.index("## 1. Introduction"):]
    b = re.sub(r"\$\$.*?\$\$", " EQ ", body, flags=re.S)
    b = re.sub(r"\$[^$]*\$", " x ", b)
    b = re.sub(r"!\[.*?\]\([^)]*\)", " ", b, flags=re.S)
    b = re.sub(r"(?m)^\|.*$", " ", b)
    b = re.sub(r"(?m)^#{1,6} .*$", " ", b)                # nor a heading, which has no full stop
    b = " ".join(b.split())
    return [(len(t.split()), t) for t in re.split(r"[.;:!?]\s+", b)]


if "--selftest" in sys.argv:
    import tempfile, os
    d = tempfile.mkdtemp()
    hard = " ".join(["word"] * 75) + "."
    soft = " ".join(["word"] * 40) + "; " + " ".join(["word"] * 40) + "."
    for name, s, want in (("bad.md", hard, True), ("good.md", soft, False)):
        f = os.path.join(d, name)
        io.open(f, "w", encoding="utf-8").write("## 1. Introduction\n\n" + s + "\n")
        over = [w for w, _ in stretches(f) if w > DEFAULT]
        assert bool(over) == want, "plant: %s gave %r, wanted over=%r" % (name, over, want)
    f = os.path.join(d, "app.md")
    io.open(f, "w", encoding="utf-8").write("## 1. Introduction\n\ntext.\n\n## Appendices\n\n" + hard + "\n")
    assert not [w for w, _ in stretches(f) if w > DEFAULT], "plant: the appendix exemption is not applied"
    # a semicolon must genuinely count as a break, or the check is just the word cap again
    assert max(w for w, _ in stretches(os.path.join(d, "good.md"))) == 40, \
        "plant: the semicolon did not split, so this measures nothing new"
    print("   plant: an unbroken 75-word run is caught, the same words split by a semicolon")
    print("          are spared, and the appendix exemption is real: yes")
    sys.exit(0)

bad = 0
for path in [a for a in sys.argv[1:] if a.endswith(".md")]:
    name = path.split("/")[-1]
    cap = CAPS.get(name, DEFAULT)
    rows = stretches(path)
    over = sorted(((w, t) for w, t in rows if w > cap), reverse=True)
    if over:
        bad = 1
        for w, t in over:
            print('   %s: %d words with no hard break, cap %d  <-- ISSUE' % (name, w, cap))
            print('      "%s..."' % t[:110])
    else:
        print("   %-22s longest unbroken run %d words, cap %d"
              % (name, max(w for w, _ in rows), cap))
sys.exit(bad)

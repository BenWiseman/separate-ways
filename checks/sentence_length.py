#!/usr/bin/env python3
"""No body sentence may run past a word ceiling.

    python3 checks/sentence_length.py FILE.md [FILE.md ...] [--max N]
    python3 checks/sentence_length.py --selftest

WHY. Ben, 2026-09-27: readability and accessibility are the priority. Nothing measured sentence
length, and the cosmology paper's body carried one of 113 words, a four-item list run together
with semicolons and a long parenthetical inside the first item. Four of the worst were split;
this stops the next one. Appendices are exempt: a derivation carrying its own qualifiers is a
different kind of sentence from one in the body, and holding them to the same ceiling would cost
precision for no reader's benefit.
"""
import io, re, sys

CAP = 80


def sentences(path):
    c = io.open(path, encoding="utf-8").read()
    body = c.split("## Appendi")[0]
    if "## 1. Introduction" in body:
        body = body[body.index("## 1. Introduction"):]
    b = re.sub(r"\$\$.*?\$\$", " EQ ", body, flags=re.S)   # display maths is not prose
    b = re.sub(r"\$[^$]*\$", " x ", b)                     # nor is inline maths, but it is a word
    b = re.sub(r"!\[.*?\]\([^)]*\)", " ", b, flags=re.S)   # nor a figure caption's image line
    b = re.sub(r"(?m)^\|.*$", " ", b)                      # nor a table row
    b = re.sub(r"(?m)^#{1,6} .*$", " ", b)                # nor a heading, which has no full stop
    return [s for s in re.split(r"(?<=[.!?]) +", " ".join(b.split())) if len(s.split()) > 2]


if "--selftest" in sys.argv:
    import tempfile, os
    d = tempfile.mkdtemp()
    long_one = " ".join(["word"] * (CAP + 5)) + "."
    ok_one = " ".join(["word"] * (CAP - 5)) + "."
    for name, s, want in (("bad.md", long_one, 1), ("good.md", ok_one, 0)):
        f = os.path.join(d, name)
        io.open(f, "w", encoding="utf-8").write("## 1. Introduction\n\n" + s + "\n")
        over = [x for x in sentences(f) if len(x.split()) > CAP]
        assert len(over) == want, "plant: %s gave %d over the cap, wanted %d" % (name, len(over), want)
    # and the appendix exemption has to be real
    f = os.path.join(d, "app.md")
    io.open(f, "w", encoding="utf-8").write("## 1. Introduction\n\ntext.\n\n## Appendices\n\n" + long_one + "\n")
    assert not [x for x in sentences(f) if len(x.split()) > CAP], "plant: the appendix exemption is not applied"
    print("   plant: a sentence over the cap is caught, one under it is spared, and the")
    print("          appendix exemption is real: yes")
    sys.exit(0)

cap = CAP
if "--max" in sys.argv:
    cap = int(sys.argv[sys.argv.index("--max") + 1])
bad = 0
for path in [a for a in sys.argv[1:] if a.endswith(".md")]:
    ss = sentences(path)
    over = sorted(((len(s.split()), s) for s in ss if len(s.split()) > cap), reverse=True)
    name = path.split("/")[-1]
    if over:
        bad = 1
        for n, s in over:
            print('   %s: a %d-word sentence, cap %d  <-- ISSUE' % (name, n, cap))
            print('      "%s..."' % s[:110])
    else:
        ls = sorted(len(s.split()) for s in ss)
        print("   %-22s %4d body sentences, longest %d, mean %.1f, cap %d"
              % (name, len(ss), ls[-1], sum(ls) / len(ls), cap))
sys.exit(bad)

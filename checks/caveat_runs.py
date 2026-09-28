#!/usr/bin/env python3
"""No paragraph of either body may run four or more caveats in a row.

    python3 checks/caveat_runs.py FILE.md [FILE.md ...]
    python3 checks/caveat_runs.py --selftest

WHY. REORIENT names the drift to watch for as premature reasonableness: qualifying a real result
into nothing. Ben's form of the rule is that a result with a caveat is a result and a result under
four caveats reads as an apology. Nothing measured it. Overall caveat DENSITY turns out to say
nothing useful, since both bodies sit near one sentence in five and there is no baseline to judge
that against, but consecutive runs are exactly the failure mode and they are countable.

Measured on 2026-09-28: neither body has a run of four, so the bar is zero and any run is new.
The two that exist are in appendices, where a numbered series of distinctions is the content
rather than a hedge, and appendices are exempt for the same reason gate 21 exempts them.
"""
import io, re, sys

CAVEAT = re.compile(r"\b(?:is not|are not|does not|do not|cannot|is no|nothing (?:here|in this)|"
                    r"neither|not a claim|we do not|this paper does not|makes no|claims no)\b",
                    re.I)
CAP = 4


def runs(path):
    c = io.open(path, encoding="utf-8").read()
    body = c.split("## Appendi")[0]
    if "## 1. Introduction" in body:
        body = body[body.index("## 1. Introduction"):]
    b = re.sub(r"\$\$.*?\$\$", " ", body, flags=re.S)
    b = re.sub(r"\$[^$]*\$", " x ", b)
    b = re.sub(r"!\[.*?\]\([^)]*\)", " ", b, flags=re.S)
    out = []
    for para in b.split("\n\n"):
        ss = [s for s in re.split(r"(?<=[.!?]) +", " ".join(para.split())) if len(s.split()) > 4]
        run, first = 0, ""
        for s in ss:
            if CAVEAT.search(s):
                if run == 0:
                    first = s
                run += 1
                if run >= CAP:
                    out.append((run, first))
            else:
                run = 0
    # keep only the longest report per starting sentence
    seen, keep = set(), []
    for n, f in sorted(out, reverse=True):
        if f not in seen:
            seen.add(f); keep.append((n, f))
    return keep


if "--selftest" in sys.argv:
    import os, tempfile
    d = tempfile.mkdtemp()
    caveat = "It is not the case that anything follows from this at all."
    plain = "The computed value comes out at one half of the horizon radius."
    for name, body, want in (("bad.md", " ".join([caveat]*4), 1),
                             ("ok.md", " ".join([caveat]*3 + [plain] + [caveat]*3), 0),
                             ("split.md", " ".join([caveat]*2 + [plain] + [caveat]*2), 0)):
        f = os.path.join(d, name)
        io.open(f, "w", encoding="utf-8").write("## 1. Introduction\n\n" + body + "\n")
        got = len(runs(f))
        assert got == want, "plant: %s gave %d runs, wanted %d" % (name, got, want)
    f = os.path.join(d, "app.md")
    io.open(f, "w", encoding="utf-8").write("## 1. Introduction\n\ntext here is fine and plain.\n\n"
                                            "## Appendices\n\n" + " ".join([caveat]*6) + "\n")
    assert not runs(f), "plant: the appendix exemption is not applied"
    print("   plant: four caveats in a row are caught, three are spared, a break resets the run,")
    print("          and the appendix exemption is real: yes")
    sys.exit(0)

bad = 0
for path in [a for a in sys.argv[1:] if a.endswith(".md")]:
    r = runs(path)
    name = path.split("/")[-1]
    if r:
        bad = 1
        for n, first in r:
            print("   %s: %d caveats in a row  <-- ISSUE" % (name, n))
            print('      starts "%s..."' % first[:110])
    else:
        print("   %-22s no paragraph runs %d caveats in a row" % (name, CAP))
sys.exit(bad)

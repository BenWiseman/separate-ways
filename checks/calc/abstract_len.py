#!/usr/bin/env python3
"""abstract_len.py -- arXiv caps the abstract at 1920 characters. Measure it by finding the
heading, not by a fixed line window: an earlier ad-hoc check used lines 10:34 and silently
stopped counting when the abstract grew past line 34, reporting 1911 for a 2096-character
abstract. A checker that cannot fail is worse than no checker, so --validate plants an
over-long abstract and a truncation trap and requires both to be caught.
"""
import io, re, sys

CAP = 1920

def abstract_of(path):
    s = io.open(path, encoding="utf-8").read()
    m = re.search(r"(?m)^##\s+Abstract\s*$", s)
    if not m:
        raise SystemExit(f"{path}: no '## Abstract' heading")
    rest = s[m.end():]
    n = re.search(r"(?m)^#{1,6}\s+\S", rest)
    body = rest[:n.start()] if n else rest
    return " ".join(body.split())

def check(path, quiet=False):
    a = abstract_of(path)
    ok = len(a) <= CAP
    if not quiet:
        print(f"  {path}: abstract {len(a)} chars, cap {CAP} -> {'OK' if ok else 'OVER by %d' % (len(a)-CAP)}")
        print(f"    starts: {a[:70]}...")
        print(f"    ends:   ...{a[-70:]}")
    return ok

def validate():
    import tempfile, os
    print("=== VALIDATION ===")
    good = "## Abstract\n\n" + ("word " * 100) + "\n\n## 1. Introduction\n\nbody\n"
    over = "## Abstract\n\n" + ("word " * 500) + "\n\n## 1. Introduction\n\nbody\n"
    # the truncation trap: an abstract whose tail sits far below the heading
    trap = "## Abstract\n\n" + ("word word word word word word word word\n" * 60) + "\n\n## 1. Introduction\n\nbody\n"
    fails = []
    for name, text, want in (("short", good, True), ("over-long", over, False), ("multi-line over-long", trap, False)):
        p = os.path.join(tempfile.gettempdir(), "_abs_probe.md")
        io.open(p, "w", encoding="utf-8").write(text)
        got = check(p, quiet=True)
        status = "as expected" if got == want else "WRONG"
        print(f"  {name:22s} -> {'under' if got else 'over':5s} cap, {status}")
        if got != want:
            fails.append(name)
    print("  " + ("every case handled" if not fails else "BLIND on: " + ", ".join(fails)))
    return 1 if fails else 0

if __name__ == "__main__":
    if "--validate" in sys.argv:
        sys.exit(validate())
    paths = [a for a in sys.argv[1:] if not a.startswith("--")] or ["papers/2_over_the_horizon/COMPANION_v1.md"]
    sys.exit(0 if all(check(p) for p in paths) else 1)

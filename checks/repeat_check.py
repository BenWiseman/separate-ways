#!/usr/bin/env python3
"""No sentence should appear twice in the manuscript.

`closer_audit.py` is a screen: it lists paragraph-closing shapes so a repeated flourish becomes
visible to a person. Most of what it lists is content and needs judgement. One case needs none.
"Nothing was tuned to arrange that" stood in this paper four times, twice inside one paragraph of
the introduction, and a reader meeting it the third time stops reading the physics and starts
noticing the writer. That is mechanical and so it is gated here.

Sentences of six words or more, normalised for whitespace and case, must be unique. Short ones are
exempt because "That is the point." is allowed to recur and "The proof is one substitution." is a
statement of fact rather than a flourish.

    python3 checks/repeat_check.py papers/2_over_the_horizon/COMPANION_v1.md
    python3 checks/repeat_check.py papers/2_over_the_horizon/COMPANION_v1.md --selftest
"""
import io, re, sys
from collections import defaultdict

MINWORDS = 6


def sentences(text):
    out = []
    for block in text.split("\n\n"):
        b = block.strip()
        if not b or b.startswith(("#", "$$", "|", ">")):
            continue
        flat = re.sub(r"\s+", " ", b)
        flat = re.sub(r"\$[^$]*\$", " MATH ", flat)          # a shared formula is not a shared sentence
        for sent in re.split(r"(?<=[.!?]) ", flat):
            sent = sent.strip(" *_")
            if len(sent.split()) >= MINWORDS:
                out.append(sent)
    return out


def duplicates(text):
    seen = defaultdict(int)
    for s in sentences(text):
        seen[s.lower()] += 1
    return {k: v for k, v in seen.items() if v > 1}


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "papers/2_over_the_horizon/COMPANION_v1.md"
    text = io.open(path, encoding="utf-8").read()

    if "--selftest" in sys.argv:
        plant = text + "\n\nThis planted sentence exists to make the checker fail on purpose.\n\n" \
                       "This planted sentence exists to make the checker fail on purpose.\n"
        got = duplicates(plant)
        ok = any("planted sentence exists" in k for k in got)
        print("  repeat_check: catches a planted duplicate: %s" % ("yes" if ok else "NO"))
        clean = duplicates(text)
        print("  repeat_check: the real file is clean: %s" % ("yes" if not clean else "NO"))
        sys.exit(0 if ok else 1)

    dup = duplicates(text)
    n = len(sentences(text))
    if not dup:
        print("   %d sentences of %d+ words, none repeated" % (n, MINWORDS))
        sys.exit(0)
    print("   %d repeated sentence(s) out of %d  <-- ISSUE" % (len(dup), n))
    for k, v in sorted(dup.items(), key=lambda kv: -kv[1]):
        print("      %dx  %s" % (v, k[:96]))
    sys.exit(1)

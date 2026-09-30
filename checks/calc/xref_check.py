#!/usr/bin/env python3
"""Every internal cross-reference in the companion must resolve to a heading that exists.

Four renumberings in one session broke references twice, and both times a human read
caught it rather than a check. This is the check. It knows three kinds of reference:

  Section N / Sections N and M / SS N   -> a "## N." heading
  §N or §N.M                            -> a "## N." or "### N.M" heading
  A.N                                   -> an "**A.N " appendix heading

References that name the other paper ("the cosmology paper's §2.1") are skipped, because
that paper's headings are not in this file.

--validate plants a reference to a section that does not exist and checks the tool reports it.
"""
import re, sys

def headings(s):
    sec = set(re.findall(r"(?m)^##\s+(\d+)\.\s", s))
    sub = set(re.findall(r"(?m)^###\s+(\d+\.\d+)\s", s))
    app = set(re.findall(r"(?m)^\*\*(A\.\d+)\s", s))
    return sec, sub, app

def check(path, quiet=False):
    s = open(path, encoding="utf-8").read()
    sec, sub, app = headings(s)
    bad = []
    # strip references that name the other paper, so they are not tested against this one
    # line breaks fall between "paper's" and the reference, so the stripper must span them
    # Blank each match to spaces of the SAME length, so the offsets used for the context
    # line below still point at the right place in the original text.
    blank = lambda m: " " * len(m.group(0))
    s_own = re.sub(r"(?i)the\s+cosmology\s+paper's\s+(§\s*[\d.]+|[AC]\.\d+)", blank, s)
    # and the appendix heading legitimately names the other paper's A.1 to A.5
    s_own = re.sub(r"numbered from A\.6 because A\.1 to A\.5 belong to the cosmology paper\.",
                   blank, s_own)
    # A reference qualified as somebody else's is not a reference to this document. Both papers
    # discuss other people's sections by number, and on the cosmology paper that produced five
    # dangling reports where nothing was wrong, which is the shape that teaches a reader to skip
    # the gate. The qualifier has to be explicit: "of their paper", "of that paper", or a
    # parenthetical naming the companion.
    s_own = re.sub(r"(?i)(?:§\s*[\d.]+|Sections?\s+[\d.]+)\s+of\s+(?:their|that|the\s+other)"
                   r"\s+paper", blank, s_own)
    s_own = re.sub(r"(?i)companion,\s*(?:§\s*)?[A-Z]\.\d+", blank, s_own)
    for m in re.finditer(r"§\s*(\d+(?:\.\d+)?)", s_own):
        r = m.group(1)
        if ("." in r and r not in sub) or ("." not in r and r not in sec):
            bad.append(("§"+r, m.start()))
    for m in re.finditer(r"\bSections?\s+(\d+(?:\.\d+)?)", s_own):
        r = m.group(1)
        if ("." in r and r not in sub) or ("." not in r and r not in sec):
            bad.append(("Section "+r, m.start()))
    for m in re.finditer(r"\b(A\.\d+)\b", s_own):
        if m.group(1) not in app:
            bad.append((m.group(1), m.start()))
    if not quiet:
        print(f"=== {path} ===")
        print(f"  headings: {len(sec)} sections, {len(sub)} subsections, {len(app)} appendices")
        if bad:
            for r, pos in bad:
                ctx = re.sub(r"\s+", " ", s[max(0, pos-70):pos+40])
                print(f"  DANGLING {r:<12} ...{ctx}...")
        print(f"  {len(bad)} dangling reference(s)")
    return bad

def check_external(companion="papers/2_over_the_horizon/COMPANION_v1.md", cosmology="papers/1_separate_ways/PAPER2_v4_draft.md",
                   quiet=False):
    """Each paper's references to the other must name a heading the other paper has.

    Added 2026-09-23 after the companion was found attributing a claim to the cosmology
    paper's 4.2 that the de-noising rewrite had removed. The internal checker cannot see
    this, because 4.2 exists in the companion too.

    It then sat below the __main__ block's sys.exit for a day, so it could not run even if
    something had called it, and nothing did. Found 2026-09-24 while adding the reverse
    direction. Both directions run now and the gate calls them.
    """
    c = open(companion, encoding="utf-8").read()
    m = open(cosmology, encoding="utf-8").read()
    secs = set(re.findall(r"(?m)^#{2,3}\s+(\d+(?:\.\d+)?)", m))
    apps = set(re.findall(r"\*\*([A-Z]\.\d+)\s", m))
    bad = []
    for mm in re.finditer(r"cosmology paper'?s\s+(?:Appendix\s+)?(?:§\s*)?([\d.]+|[A-Z]\.\d+)", c):
        r = mm.group(1).rstrip(".")
        if r not in secs and r not in apps:
            bad.append((r, mm.start()))
    print(f"=== external cross-references, {companion} -> {cosmology} ===")
    print(f"  cosmology paper has {len(secs)} sections and {len(apps)} appendices")
    for r, pos in bad:
        print(f"  MISSING {r:8s} ...{re.sub(r'[ \n]+',' ', c[max(0,pos-70):pos+50])}...")
    if not quiet:
        print(f"  {len(bad)} reference(s) to a heading the cosmology paper does not have")

    # the other direction, which nothing checked: "(companion, A.10)" and "the companion's 5.2"
    csecs = set(re.findall(r"(?m)^#{2,3}\s+(\d+(?:\.\d+)?)", c))
    capps = set(re.findall(r"(?m)^\*\*([A-Z]\.\d+)\s", c))
    back = []
    for mm in re.finditer(r"(?i)companion(?:'s)?,?\s+(?:Appendix\s+)?(?:§\s*)?"
                          r"([A-Z]\.\d+|\d+(?:\.\d+)?)\b", m):
        r = mm.group(1).rstrip(".")
        if r not in csecs and r not in capps:
            back.append((r, mm.start()))
    if not quiet:
        print(f"=== external cross-references, {cosmology} -> {companion} ===")
        print(f"  companion has {len(csecs)} sections and {len(capps)} appendices")
        for r, pos in back:
            print(f"  MISSING {r:8s} ...{re.sub(r'[ \n]+',' ', m[max(0,pos-70):pos+50])}...")
        print(f"  {len(back)} reference(s) to a heading the companion does not have")
    return bad + back


def validate_external():
    """Plant a bad reference in each direction and confirm both are caught."""
    import tempfile, os
    c = open("papers/2_over_the_horizon/COMPANION_v1.md", encoding="utf-8").read()
    m = open("papers/1_separate_ways/PAPER2_v4_draft.md", encoding="utf-8").read()
    d = tempfile.mkdtemp()
    cpath, mpath = os.path.join(d, "c.md"), os.path.join(d, "m.md")
    open(cpath, "w", encoding="utf-8").write(c + "\n\nSee the cosmology paper's A.97 for this.\n")
    open(mpath, "w", encoding="utf-8").write(m + "\n\nSee the companion, A.96, for this.\n")
    planted = check_external(cpath, mpath, quiet=True)
    names = {r for r, _ in planted}
    clean = check_external(quiet=True)
    print("=== VALIDATION, external ===")
    print(f"  planted A.97 one way and A.96 the other; caught {sorted(names & {'A.97', 'A.96'})}")
    print(f"  the real pair reports {len(clean)}")
    ok = {"A.97", "A.96"} <= names and not clean
    print("  both directions catch their plant and the real pair is clean" if ok
          else "  A PLANTED CASE WAS MISSED OR THE REAL PAIR IS DIRTY  <-- BLIND")
    return 0 if ok else 1

if __name__ == "__main__":
    if "--external" in sys.argv:
        sys.exit(1 if check_external() else 0)
    if "--validate-external" in sys.argv:
        sys.exit(validate_external())
    if "--validate" in sys.argv:
        src = open("papers/2_over_the_horizon/COMPANION_v1.md", encoding="utf-8").read()
        probe = "/tmp/_xref_probe.md"
        open(probe, "w", encoding="utf-8").write(src + "\n\nA planted reference to Section 99 and to A.99 and to §42.7.\n")
        bad = check(probe, quiet=True)
        names = {b[0] for b in bad}
        want = {"Section 99", "A.99", "§42.7"}
        missing = want - names
        print("=== VALIDATION ===")
        print("  planted:", sorted(want))
        print("  caught: ", sorted(names & want))
        if missing:
            print("  MISSED:", sorted(missing), " <-- BLIND")
            sys.exit(1)
        print("  every planted dangling reference was caught")
        sys.exit(0)
    sys.exit(1 if check(sys.argv[1] if len(sys.argv) > 1 else "papers/2_over_the_horizon/COMPANION_v1.md") else 0)

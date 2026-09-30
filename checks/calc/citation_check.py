#!/usr/bin/env python3
"""Does each citation point at the paper the sentence says it does?

Written 2026-09-24 after finding that it did not. The companion's citation numbers had been
written against an OLDER version of the cosmology paper's bibliography, which was later
renumbered, and a few later citations against the new one. So "as in the work of Halliwell,
Hawking and Kiefer [17-19,45]" resolved, in the current list, to four Turok and Boyle papers.
Nothing caught it because nothing was checking the names.

The check: for each citation, take the prose in front of it, pull out the capitalised surnames,
and require at least one of them to appear in the entry cited. Citations with no surname in
front of them (a bare number after a formula, a collaboration named only by acronym) are
reported separately rather than counted as failures, because the check has no grip there.

    python3 checks/calc/citation_check.py [paper.md]
    python3 checks/calc/citation_check.py --validate
"""
import io, re, sys, os

STOP = {"The","This","That","These","Those","Section","Sections","Appendix","Figure","Write","With",
        "At","In","On","So","And","But","For","From","Two","Three","Both","Its","It","We","A","An",
        "GeV","PeV","TeV","CPT","KMS","WKB","Big","Bang","Kruskal","Hadamard","Hilbert","Killing",
        "Schwarzschild","Kerr","Sitter","Bunch","Davies","Gibbons-Hawking","Planck","Euclidean",
        "Legendre","Nariai","Smarr","Weyl-","Formation","Positive","Electroweak","Note","Read",
        "Every","One","Under","Where","When","Whether","What","Nothing","Here","There","Now"}

def refs_of(text):
    i = text.find("\n## References")
    if i < 0: return {}, text
    body, rest = text[:i], text[i:]
    ent = {int(m.group(1)): m.group(2) for m in re.finditer(r"(?m)^(\d+)\\\.\s+(.*)$", rest)}
    return ent, body

def is_surname(word, ent, skip=None):
    """Does this word sit in a surname position in some OTHER entry?

    A surname in these lists follows an initial or a comma: "J. M. Weisberg, Y. Huang".
    Requiring that keeps "Bogoliubov coefficients", "Pulsar timing" and "Lorentzian" from
    counting as citations to people, which is what buried the real mismatches on the first pass.
    """
    pat = re.compile(r"(?:^|[.,]\s+|\s)%s\b" % re.escape(word), re.I)
    sur = re.compile(r"(?:[A-Z]\.\s*|,\s*)%s\b" % re.escape(word))
    for k, v in ent.items():
        if k == skip: continue
        if sur.search(v): return True
    return False

def expand(g):
    out = []
    for p in g.split(","):
        p = p.strip()
        if "-" in p and all(x.strip().isdigit() for x in p.split("-", 1)):
            a, b = p.split("-", 1); out += list(range(int(a), int(b) + 1))
        elif p.isdigit(): out.append(int(p))
    return out

def check(path, quiet=False):
    text = io.open(path, encoding="utf-8").read()
    ent, body = refs_of(text)
    if not ent:
        print("  %s: no reference list" % path); return False
    # strip maths first: $w\in[0,1]$ is an interval, not a citation to references 0 and 1
    body = re.sub(r"\$\$.*?\$\$", " ", body, flags=re.S)
    body = re.sub(r"(?<!\$)\$(?:[^$\n]|\n(?!\s*\n))+?\$(?!\$)", " ", body)
    flat = re.sub(r"\s+", " ", body)
    bad, blind, ok = [], 0, 0
    for m in re.finditer(r"\[([\d,\s\-]+)\]", flat):
        ctx = flat[max(0, m.start() - 130):m.start()]
        # stop at the previous citation: the names before it belong to that one, not this
        prev = ctx.rfind("]")
        if prev >= 0: ctx = ctx[prev + 1:]
        ctx = re.sub(r"\$[^$]*\$", " ", ctx)
        names = [w.strip(".,'’s") for w in re.findall(r"\b[A-Z][a-zA-Zàéèíóúü'’-]{2,}\b", ctx)]
        names = [n for n in names if n not in ("Collaboration","Collaborations")]
        names = [n for n in names if n not in STOP and len(n) > 2]
        if not names: blind += 1; continue
        for n in expand(m.group(1)):
            if n not in ent: bad.append((n, "no such entry", ctx[-60:])); continue
            if any(s.lower() in ent[n].lower() for s in names): ok += 1
            else:
                # A prose word is only evidence if it is a surname SOMEWHERE in the
                # bibliography. "Bogoliubov coefficients" and "Pulsar timing" are not
                # citations to Bogoliubov or to Pulsar, and counting them as mismatches
                # buries the real ones.
                real = [s for s in names if len(s) > 3 and is_surname(s, ent, skip=n)]
                if real: bad.append((n, ent[n][:56], " ".join(real[-4:])))
                else: blind += 1
    if not quiet:
        print("  %s" % path)
        print("    %d citations matched a surname in the entry, %d with no surname in front of them"
              % (ok, blind))
        if bad:
            print("    %d MISMATCHED:" % len(bad))
            for n, e, ctx in bad[:12]:
                print("      [%d] prose says %-38s entry is %s" % (n, ctx, e))
        else:
            print("    0 mismatched")
    return not bad

def validate():
    import tempfile
    src = io.open("papers/2_over_the_horizon/COMPANION_v1.md", encoding="utf-8").read()
    i = src.find("\n## References")
    body, rest = src[:i], src[i:]
    # swap two entries so the numbers point at the wrong papers
    swapped = re.sub(r"(?m)^1\\\. (.*)$", r"1\\. J. Doe, \"An unrelated paper,\" (1999).", rest)
    p = os.path.join(tempfile.gettempdir(), "_cite_probe.md")
    io.open(p, "w", encoding="utf-8").write(body + swapped)
    print("=== VALIDATION: replacing entry 1 with an unrelated paper ===")
    caught = not check(p, quiet=True)
    print("  detector fired: %s" % ("yes" if caught else "NO, THE CHECK IS BLIND"))
    print("=== and the real file, which must pass ===")
    clean = check("papers/2_over_the_horizon/COMPANION_v1.md", quiet=True)
    print("  real file passes: %s" % ("yes" if clean else "no"))
    return 0 if (caught and clean) else 1

if __name__ == "__main__":
    if "--validate" in sys.argv: sys.exit(validate())
    args = [a for a in sys.argv[1:] if not a.startswith("--")] or ["papers/2_over_the_horizon/COMPANION_v1.md"]
    sys.exit(0 if all(check(a) for a in args) else 1)

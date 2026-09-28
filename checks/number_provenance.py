#!/usr/bin/env python3
"""Does every number the paper states beside a script citation appear in that script's output?

The paper's central rigour claim is that every quantitative claim is backed by a script cited
at the point it is used. Nothing checked the other half of that: whether the script actually
PRODUCES the number the sentence quotes. This class of drift has bitten repeatedly, and each
time it was found by accident rather than by a check:

  - the inverted-ordering neutrino sum, 98.9 where the stated inputs give 100.9
  - the DESI margin, quoted as 64.0 in one sentence and 64.2 in the next
  - the same PBH fraction written 0.0101 and 0.01014 two paragraphs apart
  - a duty cycle held fixed at its mean where the growth law varies it per object

Running it the first time found two real provenance gaps, both since closed: endpoint_power.R
printed its k=3 false-alarm rate as 0.0000 while 3.2 quoted 3.2e-5 and three sample-wide rates
derived from it, and separable_constrained_proof.R printed raw integrals where 3.1 quotes them
divided by pi^2.

MATCHING IS AT THE TEXT'S OWN PRECISION. If the paper writes 0.98, a script printing 0.978
matches; if the paper writes 0.978, a script printing 0.98 does not. That asymmetry is the
point: the paper is claiming the digits it prints.

Usage:
    python3 tools/number_provenance.py            # runs the cited scripts (a few minutes)
    python3 tools/number_provenance.py --cache D  # reuse outputs already in D
    python3 tools/number_provenance.py --validate # plant a wrong number, check it is caught

A hit is a thing to LOOK at. Legitimate reasons a number is absent: it is an input rather than
a result, it is a literature value, or the script prints a quantity the paper converts. Each
one that survives review belongs in EXEMPT below with its reason.
"""
import sys, re, io, os, sys, math, subprocess, tempfile, shutil
from concurrent.futures import ThreadPoolExecutor

DOCS   = tuple(a for a in sys.argv[1:] if not a.startswith("--")) or \
         ("paper/PAPER2_v3.md", "paper/COMPANION_v1.md")
WINDOW = 320          # the sentence the citation sits in, roughly

# script -> [numbers that are legitimately absent, with the reason]
EXEMPT = {}
_UNITS = {w: i for i, w in enumerate(
    "zero one two three four five six seven eight nine ten eleven twelve thirteen fourteen "
    "fifteen sixteen seventeen eighteen nineteen".split())}
_TENS = {w: (i + 2) * 10 for i, w in enumerate(
    "twenty thirty forty fifty sixty seventy eighty ninety".split())}


def _words_to_number(phrase):
    """'twenty-four' -> 24, 'ninety-five' -> 95, 'a hundred' -> 100. None if not a number."""
    p = phrase.lower().replace("-", " ").split()
    if p and p[0] in ("a", "an"):
        p = p[1:]
    total = 0
    for w in p:
        if w in _UNITS:
            total += _UNITS[w]
        elif w in _TENS:
            total += _TENS[w]
        elif w == "hundred":
            total = (total or 1) * 100
        elif w == "thousand":
            total = (total or 1) * 1000
        else:
            return None
    return total or None


def nums_from_text(s):
    """Numbers a reader would take as claims. Returns (value, decimals, literal).

    Spelled-out numbers count too. The checker read digits only, so "twenty-four orders of
    magnitude", "nineteen times the propagated width" and "twelve thousand events" were
    invisible to it: roughly a hundred quantitative claims in this paper sat outside the audit.
    Only unambiguous quantity forms are taken, N per cent / orders / times / widths, because a
    bare "one" is usually the article.
    """
    out, sci_spans = [], []
    _WORD = (r"\b((?:a |an )?(?:one|two|three|four|five|six|seven|eight|nine|ten|eleven|twelve|"
             r"thirteen|fourteen|fifteen|sixteen|seventeen|eighteen|nineteen|twenty|thirty|"
             r"forty|fifty|sixty|seventy|eighty|ninety|hundred|thousand)"
             r"(?:[- ](?:one|two|three|four|five|six|seven|eight|nine|hundred|thousand))*)"
             r"\s+(per cent|orders|times|widths)\b")
    for m in re.finditer(_WORD, s, re.I):
        v = _words_to_number(m.group(1))
        if v is not None:
            out.append((float(v), 0, m.group(0), False))
    # LaTeX scientific: 2.1\times10^{-53}  or  3\times10^5
    for m in re.finditer(r'(\d+(?:\.\d+)?)\s*\\times\s*10\^\{?(-?\d+)\}?', s):
        mant, exp = m.group(1), int(m.group(2))
        dec = len(mant.split('.')[1]) if '.' in mant else 0
        out.append((float(mant) * 10**exp, dec, m.group(0), True))
        sci_spans.append(m.span())   # the mantissa is NOT a separate claim
    # plain decimals, but not inside a citation [12] or a section ref
    for m in re.finditer(r'(?<![\w^{\u00a7.])(\d+\.\d+)(?![\d}])', s):
        lit = m.group(1)
        if any(a <= m.start() < b for a, b in sci_spans):
            continue
        pre = s[max(0, m.start()-12):m.start()].lower()
        if any(k in pre for k in ('section', 'appendix', 'table', 'figure', 'eq.')):
            continue
        out.append((float(lit), len(lit.split('.')[1]), lit, False))
    return out

# A script that prints the manuscript's own figure makes every number self-certifying.
# wdw_weights.R printed "aH = 1.95376 (paper: 1.95374) diff 2.5e-05" for days: the paper's
# 1.95374 was wrong, the script said so on every run, and the check passed because the wrong
# number appeared in the output, inside the parenthetical that existed to flag it. Lines that
# quote the paper are therefore not evidence and are dropped before the numbers are read.
_ECHO = re.compile(r"paper\s*:|paper says|paper quotes|\bquotes\b|quoted:|manuscript says"
                   r"|A\.\d+\s+quotes|\u00a7\s*[\d.]+\s+quotes", re.I)

def nums_from_out(s):
    vals = []
    for line in s.splitlines():
        hit = _ECHO.search(line)
        if hit:
            # keep what the script computed and drop the echo that follows it, rather
            # than the whole line: most of these print both on one line, as in
            # "C = 0.9581   the paper says 0.958".
            line = line[:hit.start()]
        for m in re.finditer(r'-?\d+\.?\d*(?:[eE][-+]?\d+)?', line):
            try: vals.append(float(m.group(0)))
            except ValueError: pass
    return vals

def matches(claim, dec, sci, outs):
    """Does some output value agree with the claim at the claim's own precision?

    Scripts print fractions where the paper writes percentages (0.6985 against 69.9 per
    cent), so a factor of 100 either way counts as found.
    """
    for v0 in outs:
      for v in (v0, v0 * 100.0, v0 / 100.0):
        if v == 0 and claim == 0: return True
        if claim == 0: continue
        if sci:
            if v == 0: continue
            if math.copysign(1, v) != math.copysign(1, claim): continue
            e = math.floor(math.log10(abs(claim)))
            try:
                if abs(v/10**e - claim/10**e) <= 0.5*10**(-dec) + 1e-9: return True
            except (ValueError, OverflowError): pass
        else:
            if abs(v - claim) <= 0.5*10**(-dec) + 1e-9*max(1, abs(claim)): return True
    return False


def cited_scripts():
    cite = set()
    for doc in DOCS:
        if not os.path.exists(doc): continue
        t = io.open(doc, encoding="utf-8").read()
        for m in re.finditer(r'`([^`\n]*?\.R)`', t):
            cite.add(m.group(1))
    return sorted(c for c in cite if os.path.exists(c))


def interpreter(script):
    """How to run one script. R was the only language when this was written; a repository
    that also holds Python has to be able to cite it, and a hard-coded Rscript silently
    produced an empty output file for every .py row instead of an error, so every number in
    those rows failed to match and looked like a manuscript defect."""
    if not script.endswith(".py"):
        return ["Rscript", script]
    for venv in (".venv/bin/python", "venv/bin/python"):
        if os.path.exists(venv):
            return [os.path.abspath(venv), script]
    return [sys.executable, script]


def run_all(scripts, outdir):
    os.makedirs(outdir, exist_ok=True)
    def one(s):
        dst = os.path.join(outdir, s.replace('/', '_') + ".out")
        with open(dst, "w") as fh:
            subprocess.run(interpreter(s), stdout=fh, stderr=subprocess.STDOUT, timeout=1800)
    with ThreadPoolExecutor(max_workers=min(30, (os.cpu_count() or 4))) as ex:
        list(ex.map(one, scripts))


def uncited(paper_text, min_dec=3):
    """High-precision numbers with no .R citation in reach.

    audit() walks citations and looks backward, so a number with no citation near it is never
    examined and passes silently. A blind audit of the appendices found four such numbers,
    including "four fundamental domains of one involution give -0.63, -2.64, -4.17 and +3.18"
    with no formula, no units and no script. Those gaps are invisible to a checker that only
    verifies what is already cited, so they are listed here instead of being passed over.

    Only numbers carrying min_dec decimals or an exponent are reported: at one decimal the
    tolerance is +-0.05 and the check has no power anyway, so flagging those would be noise.
    """
    flat = re.sub(r'\s+', ' ', paper_text)
    cits = [m.start() for m in re.finditer(r'`[^`\n]*?\.R`', flat)]
    # scan with spans. A first version looked each literal up with flat.find(), which returns
    # the FIRST occurrence every time, so a number appearing twice was reported twice and both
    # reports tested the citation distance at the same wrong place.
    spans = []
    for m in re.finditer(r'(\d+(?:\.\d+)?)\s*\\times\s*10\^\{?(-?\d+)\}?', flat):
        mant = m.group(1)
        dec = len(mant.split('.')[1]) if '.' in mant else 0
        spans.append((m.start(), m.group(0), dec, True))
    sci_spans = [(a, a + len(l)) for a, l, _, _ in spans]
    for m in re.finditer(r'(?<![\w^{\u00a7.])(\d+\.\d+)(?![\d}])', flat):
        if any(a <= m.start() < b for a, b in sci_spans): continue
        pre = flat[max(0, m.start()-12):m.start()].lower()
        if any(k in pre for k in ('section', 'appendix', 'table', 'figure', 'eq.')): continue
        spans.append((m.start(), m.group(1), len(m.group(1).split('.')[1]), False))
    out = []
    for at, lit, dec, sci in sorted(spans):
        if dec < min_dec and not sci: continue
        if any(abs(c - at) <= WINDOW for c in cits): continue
        out.append((lit, flat[max(0, at - 60):at + len(lit) + 30]))
    return out


def audit(outdir, paper_text):
    # 2026-09-24: the window used to be cut at the previous SENTENCE boundary, so a citation
    # covered only the sentence it closed. The papers do not cite that way: a citation closes
    # the passage it belongs to, and the companion's own availability note says every claim
    # "names the script that produces it, at the point it is used". Splitting one long sentence
    # into three therefore dropped its numbers out of the checked set without changing a single
    # citation, and the count fell from 98 to 85 on a day when nothing lost its provenance. The
    # cut is now the PARAGRAPH boundary. A number paired this way that the script does not
    # produce still reports as not found, which is the outcome worth having.
    paper_text = paper_text.replace("\n\n", " \x01 ")
    flat = re.sub(r'\s+', ' ', paper_text)
    cits = [(m.start(), m.end(), m.group(1)) for m in re.finditer(r'`([^`\n]*?\.R)`', flat)]
    groups, cur = [], []
    for c in cits:
        if cur and c[0] - cur[-1][1] <= 120: cur.append(c)
        else:
            if cur: groups.append(cur)
            cur = [c]
    if cur: groups.append(cur)

    missing, checked, seen = [], 0, []
    for grp in groups:
        outs = []
        for g in grp:
            path = os.path.join(outdir, g[2].replace('/', '_') + ".out")
            if os.path.exists(path):
                outs += nums_from_out(io.open(path, encoding="utf-8", errors="replace").read())
        if not outs: continue
        window = flat[max(0, grp[0][0] - WINDOW):grp[0][0]]
        cut = window.rfind('\x01')
        if cut >= 0: window = window[cut+1:]
        for val, dec, lit, sci in nums_from_text(window):
            names = [g[2] for g in grp]
            if any(lit in EXEMPT.get(n, []) for n in names): continue
            checked += 1
            seen.append(lit)
            if not matches(val, dec, sci, outs):
                missing.append((", ".join(names), lit))
    return checked, missing, seen


def main():
    cache = None
    if "--cache" in sys.argv: cache = sys.argv[sys.argv.index("--cache") + 1]
    validate = "--validate" in sys.argv

    scripts = cited_scripts()
    tmp = None
    if cache is None:
        tmp = tempfile.mkdtemp(prefix="numprov_")
        print(f"  running {len(scripts)} cited scripts (use --cache to reuse)...")
        run_all(scripts, tmp); cache = tmp

    # BOTH documents. The companion carries 48 script citations of its own and went
    # unaudited in the first version of this tool, which is half the work unchecked.
    missing, seen, checked, paper = [], [], 0, None
    for doc in DOCS:
        if not os.path.exists(doc): continue
        txt = io.open(doc, encoding="utf-8").read()
        if paper is None: paper = txt
        c, m, sn = audit(cache, txt)
        checked += c; seen += sn
        missing += [(f"{os.path.basename(doc)}: {sc}", lit) for sc, lit in m]
        print(f"  {os.path.basename(doc):22} {c:4d} claims checked, {len(m)} not found")

    strong = [l for l in seen if re.fullmatch(r'\d+\.\d{3,}', l) or '\\times' in l]
    print(f"\n  {checked} number-claims checked against the output of the script cited beside them")
    print(f"  {len(strong)} of them carry three or more decimals or an exponent, which is where")
    print(f"  this check has power: at one decimal the tolerance is +-0.05 and almost anything matches.")
    print(f"  {len(missing)} not found\n")
    for sc, lit in missing:
        print(f"    {lit:>18}   {sc}")

    for doc in DOCS:
        u = uncited(io.open(doc, encoding="utf-8").read())
        print(f"\n  {os.path.basename(doc):22} {len(u):4d} high-precision numbers with no script in reach")
        for lit, ctx in u[:12]:
            print(f"    {lit:>18}   ...{ctx.strip()[-64:]}")
        if len(u) > 12:
            print(f"    ... and {len(u)-12} more")

    rc = 1 if missing else 0
    if validate:
        # The echo rule, first. A script that prints the manuscript's own figure makes every
        # number self-certifying; wdw_weights.R printed "aH = 1.95376 (paper: 1.95374)" for
        # days while the paper's 1.95374 was wrong, and the check passed on the parenthetical
        # that existed to flag it. These five cases pin both halves: the echo is not evidence,
        # and the computed value on the same line still is.
        echo_cases = [
            ("computed 1.234   (paper: 9.999)  diff 8.8", 9.999, False),
            ("computed 1.234   (paper: 9.999)  diff 8.8", 1.234, True),
            ("  C = 0.9581   the paper says 0.958", 0.958, True),
            ("  A.7 quotes exact -261.627 against -261.799", 261.627, False),
            ("  magnitudes: 261.627 and 261.799, gap 0.066", 261.627, True),
        ]
        bad = [c for c in echo_cases if matches(c[1], 3, False, nums_from_out(c[0])) != c[2]]
        if bad:
            print(f"  ECHO RULE BROKEN on {len(bad)} case(s) <-- BLIND")
            rc = 1
        else:
            print("  echo rule: a quoted figure is not evidence, a computed one still is")

        # A checker that matches nothing reports success. Plant a wrong digit and require a catch.
        # Plant on a number the checker DEMONSTRABLY looks at. A first version corrupted a
        # number with no citation in its window, so the checker never examined it and the
        # plant "failed" while the checker was fine. Targeting `seen` makes that impossible.
        # a plain decimal only: corrupting a LaTeX scientific literal mangles its exponent
        # and produces something that no longer parses as a claim, so the plant fails rather
        # than the checker. That happened twice while this was being written.
        # and a HIGH-PRECISION one. At one decimal place the tolerance is +-0.05, and a
        # script's output holds enough numbers that something usually lands in that window,
        # so a 1-dp plant is missed by a checker that is working correctly. The check has
        # real power at three or more decimals; that is stated in the header and the plant
        # is chosen to match where the power is.
        target = next((l for l in seen
                       if re.fullmatch(r'\d+\.\d{3,}', l)), None)
        if target is None:
            print("\n  VALIDATE: nothing checked, so nothing to plant   <-- FIX THE CHECK")
            rc = 2
        else:
            head, _, tail = target.partition('.')
            bogus = head + '.' + ''.join('9' if c != '9' else '1' for c in tail)
            planted = paper.replace(target, bogus, 1)
            c2, m2, _ = audit(cache, planted)
            ok = any(lit == bogus for _, lit in m2)
            print(f"\n  VALIDATE: corrupted a checked claim, {target} -> {bogus}")
            print(f"            {'CAUGHT' if ok else 'MISSED   <-- THE CHECK IS BLIND'}")
            rc = 0 if (ok and not missing) else 1
    if tmp: shutil.rmtree(tmp, ignore_errors=True)
    return rc


if __name__ == "__main__":
    sys.exit(main())

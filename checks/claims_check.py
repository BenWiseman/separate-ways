#!/usr/bin/env python3
"""Verify every claim in CLAIMS.tsv, now that the manuscript carries no script paths.

Ben, 2026-09-22 and again on 2026-09-24: "you can't have script paths in a manuscript at all.
Nobody does that." He is right, and the provenance those paths carried had to survive the strip
rather than go with it. `checks/CLAIMS.tsv` is that provenance, one row per claim:

    script <TAB> anchor, the first words of the passage it backs <TAB> numbers claimed there

This checks both halves of each row. The anchor must still be findable in the manuscript, so a
rewrite that moves a passage out from under its script is caught. And every number listed must
appear in that script's output at the precision the paper states it, which is what the inline
citations used to buy.

    python3 checks/claims_check.py papers/2_over_the_horizon/COMPANION_v1.md
    python3 checks/claims_check.py papers/2_over_the_horizon/COMPANION_v1.md --reanchor

An anchor is a phrase, so it breaks whenever the passage it names is reworded. That is the point:
a passage moving out from under its script has to be noticed. But it must not become a wall, so
--reanchor shortens every broken anchor from the end, word by word, until it lands again, and
says which rows it touched. Read what it prints. If an anchor shortens to almost nothing, the
passage did not move, it went, and the row needs a person rather than a trim.
"""
import re, io, os, sys, subprocess, tempfile

# A sibling number_provenance.py wins over the one in the authoring toolkit, so the
# export repository runs on the copy it actually ships rather than on a private path
# that happens to exist on the author's machine.
sys.path.insert(0, os.path.expanduser("~/benlm/tools"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import number_provenance as np

# Paths are overridable so the same checker runs in the export repository, where the
# release sits under paper/ and checks/ instead of pub/paper2/ and checks/. The defaults
# are the working-tree layout, so nothing here changes for a plain run. Without this the
# export carried a README promising a pass that could not start.
TSV   = os.environ.get("CLAIMS_TSV", "checks/CLAIMS.tsv")
PAPERS = os.environ.get("CLAIMS_PAPERS_DIR", "pub/paper2")
PAPER = sys.argv[1] if len(sys.argv) > 1 else os.path.join(PAPERS, "COMPANION_v1.md")

# CLAIMS.tsv names each script by its working-tree path. CLAIMS_SCRIPT_MAP rewrites those
# for another layout, as "from=to,from=to", longest prefix first.
_MAP = [tuple(x.split("=", 1)) for x in
        os.environ.get("CLAIMS_SCRIPT_MAP", "").split(",") if "=" in x]
_MAP.sort(key=lambda kv: -len(kv[0]))

def remap(path):
    for a, b in _MAP:
        if path.startswith(a):
            return b + path[len(a):]
    return path

REANCHOR = "--reanchor" in sys.argv
rows = [l.rstrip("\n").split("\t") for l in io.open(TSV, encoding="utf-8")
        if l.strip() and not l.startswith("#")]
assert rows, "CLAIMS.tsv is empty: the check would pass by matching nothing"

# An anchor may live in EITHER manuscript. CLAIMS.tsv was written when only the companion was
# gated, so a row anchored in the cosmology paper reported "anchor no longer in the paper" and
# looked like a regression; that happened three times in one session before this was changed.
# Both files are searched, which keeps the moved-anchor plant working because a planted anchor is
# in neither.
def _find(name):
    """Locate a manuscript by name. The two no longer share a directory: the export splits
    them into papers/1_separate_ways/ and papers/2_over_the_horizon/. Joining a single
    PAPERS directory silently produced a path that was not there, and every anchor in the
    sibling then reported as missing from the paper, eleven of them at once."""
    here = os.path.dirname(os.path.abspath(PAPER))
    roots = [here, os.path.dirname(here), PAPERS]
    for r in roots:
        c = os.path.join(r, name)
        if os.path.exists(c):
            return c
    for r in roots:                              # one level down, papers/<paper>/<name>
        if not os.path.isdir(r):
            continue
        for d in sorted(os.listdir(r)):
            c = os.path.join(r, d, name)
            if os.path.exists(c):
                return c
    return os.path.join(PAPERS, name)            # report the path we looked for

SIBLING = _find("PAPER2_v4_draft.md" if "COMPANION" in PAPER else "COMPANION_v1.md")
paper = " ".join(io.open(PAPER, encoding="utf-8").read().split())
if os.path.exists(SIBLING):
    paper = paper + "  " + " ".join(io.open(SIBLING, encoding="utf-8").read().split())
cache = tempfile.mkdtemp(prefix="claims_")
scripts = sorted({remap(r[0]) for r in rows})
print(f"  running {len(scripts)} scripts behind {len(rows)} claims...")
np.run_all(scripts, cache)

if REANCHOR:
    paper_flat = paper
    lines = io.open(TSV, encoding="utf-8").read().rstrip("\n").split("\n")
    out, fixed = [], []
    for l in lines:
        if l.startswith("#") or not l.strip():
            out.append(l); continue
        sc, an, nu = l.split("\t")
        if an in paper_flat:
            out.append(l); continue
        w = an.split()
        while w and " ".join(w) not in paper_flat:
            w = w[:-1]
        if w:
            fixed.append((sc, len(an.split()), len(w), " ".join(w)))
            out.append(f"{sc}\t{' '.join(w)}\t{nu}")
        else:
            out.append(l)
    io.open(TSV, "w", encoding="utf-8").write("\n".join(out) + "\n")
    print(f"  re-anchored {len(fixed)} row(s):")
    for sc, n0, n1, a in fixed:
        flag = "   <-- shortened to almost nothing; check the passage still exists" if n1 < 4 else ""
        print(f"    {sc}: {n0} words -> {n1}{flag}\n      \"{a[:74]}\"")
    if not fixed:
        print("    none needed")

lost_anchor, not_found, checked = [], [], 0
for script0, anchor, nums in rows:
    script = remap(script0)
    if " ".join(anchor.split()) not in paper:
        lost_anchor.append((script, anchor)); continue
    out = os.path.join(cache, script.replace("/", "_") + ".out")
    if not os.path.exists(out):
        not_found.append((script, "(script produced no output)")); continue
    outs = np.nums_from_out(io.open(out, encoding="utf-8", errors="replace").read())
    abs_outs = [abs(v) for v in outs]
    for lit in [x for x in nums.split(";") if x]:
        # nums_from_claim, not nums_from_text: the field holds numbers and nothing else,
        # and the prose-safe parser silently dropped every bare integer in it.
        got = np.nums_from_claim(lit)
        # A field the parser cannot read was skipped here without a word, so the row looked
        # checked and was not: "1e-15" parsed to nothing for as long as this line said
        # `continue`. Say so instead, because an unread claim is worse than a failing one.
        if not got:
            not_found.append((script, lit + "   (no number the parser could read)")); continue
        # And check every number in the field, not just the first. One row listed two and the
        # second went unexamined.
        for val, dec, lit_tok, sci in got:
            checked += 1
            # A bare count is not a fraction written as a percentage, so the hundredfold
            # allowance has no business applying to it: without this a claim of 4 is satisfied
            # by any 400 in the output, and dec = 0 would make integer rows nearly vacuous.
            pct = not re.match(r"^-?\d+$", str(lit_tok).strip())
            # A row written |1.76| claims the magnitude, so compare against magnitudes.
            cand = abs_outs if str(lit_tok).startswith("|") else outs
            if not np.matches(val, dec, sci, cand, pct=pct):
                not_found.append((script, lit))
                break

if "--selftest" in sys.argv:
    # Gate 4 is the gate that matters most and nothing validated it, which is how the
    # sibling cross-document check sat unreachable for a day. Three things have to hold:
    # a corrupted number must be caught, a moved anchor must be caught, and a figure the
    # script merely quotes from the paper must not count as evidence.
    ok = True
    bad_val = np.matches(9.999, 3, False, np.nums_from_out("computed 1.234 (paper: 9.999)"))
    good_val = np.matches(1.234, 3, False, np.nums_from_out("computed 1.234 (paper: 9.999)"))
    print("  claims_check: a quoted figure is not evidence: %s" % ("yes" if not bad_val else "NO"))
    print("  claims_check: the computed one on the same line still is: %s"
          % ("yes" if good_val else "NO"))
    ok &= (not bad_val) and good_val
    # a corrupted claim must fail against its own script
    sc, an, nu = next((r for r in rows if "." in r[2].split(";")[0]), rows[0])
    lit = nu.split(";")[0]
    val, dec, _, sci = np.nums_from_text(lit)[0]
    outs = np.nums_from_out(io.open(os.path.join(cache, remap(sc).replace("/", "_") + ".out"),
                                    encoding="utf-8", errors="replace").read())
    caught = not np.matches(val * 1.7 + 0.31, dec, sci, outs)
    print("  claims_check: a corrupted number is caught: %s" % ("yes" if caught else "NO"))
    ok &= caught
    moved = " ".join("a planted anchor that is not in the manuscript".split()) not in paper
    print("  claims_check: a missing anchor is caught: %s" % ("yes" if moved else "NO"))
    ok &= moved
    # Programmer notation. Two rows were written as 1.58e31 and 1e-15; the first was checked
    # as the plain number 1.58 and failed, the second parsed to nothing and was skipped in
    # silence. Both had to be made to fail on purpose before the fix counted.
    e31 = np.nums_from_text("1.58e31")
    e15 = np.nums_from_text("1e-15")
    reads = (len(e31) == 1 and abs(e31[0][0] - 1.58e31) < 1e20 and e31[0][3]
             and len(e15) == 1 and abs(e15[0][0] - 1e-15) < 1e-25)
    print("  claims_check: 1.58e31 and 1e-15 read as numbers, mantissa not double-counted: %s"
          % ("yes" if reads else "NO"))
    ok &= reads
    # A count quoted to the unit must match a computed 36.90, and must NOT be satisfied by a
    # hundredfold coincidence. Both halves planted, because the first alone would pass with the
    # allowance still on and the row would look checked when it was only lucky.
    c37 = np.nums_from_claim("37")[0]
    rounds = np.matches(c37[0], c37[1], c37[3], np.nums_from_out("ratio 36.9036"), pct=False)
    nolucky = not np.matches(c37[0], c37[1], c37[3], np.nums_from_out("scale 3700.0"), pct=False)
    print("  claims_check: 37 matches a computed 36.90 and not a stray 3700: %s"
          % ("yes" if (rounds and nolucky) else "NO"))
    ok &= rounds and nolucky
    fake = np.nums_from_out("tau = 1.579e+31 s")
    hit = np.matches(e31[0][0], e31[0][1], e31[0][3], fake)
    miss = not np.matches(2.58e31, 2, True, fake)
    print("  claims_check: e-notation matches its script and a corrupted one does not: %s"
          % ("yes" if (hit and miss) else "NO"))
    ok &= hit and miss
    print("  claims_check: the real table passes: %s"
          % ("yes" if not (lost_anchor or not_found) else "NO"))
    ok &= not (lost_anchor or not_found)
    sys.exit(0 if ok else 1)

print(f"  {len(rows)} claims, {checked} numbers checked against the script named for them")
for s, a in lost_anchor:
    print(f"      anchor no longer in the paper: {s}  <-- ISSUE\n        \"{a[:78]}\"")
def _died(script):
    """True when the script exited non-zero, so its output is a traceback, not numbers."""
    rc = os.path.join(cache, remap(script).replace("/", "_") + ".out.rc")
    try:
        return io.open(rc).read().strip() not in ("0", "")
    except OSError:
        return False

dead = sorted({s for s, _ in not_found if _died(s)})
for s in dead:
    tail = ""
    try:
        lines = [l.rstrip() for l in io.open(
            os.path.join(cache, remap(s).replace("/", "_") + ".out"),
            encoding="utf-8", errors="replace") if l.strip()]
        tail = lines[-1][:88] if lines else ""
    except OSError:
        pass
    n_here = len([1 for t, _ in not_found if t == s])
    print(f"      DID NOT RUN: {s}  ({n_here} claim(s) unverifiable)   <-- ISSUE")
    if tail:
        print(f"        {tail}")
for s, l in not_found:
    if s in dead:
        continue
    print(f"      {l:>20}  not in the output of {s}   <-- ISSUE")
if not lost_anchor and not not_found:
    print("  every anchor still lands and every number reproduces")
sys.exit(1 if (lost_anchor or not_found) else 0)
